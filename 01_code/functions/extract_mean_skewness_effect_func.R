

ExtractMeanSkewEffect <- function(
    file,
    model_inputs,
    effect = c("temperature", "time", "mpa"),
    mpa_year = NULL,
    mpa_ref = 0) {
  
  effect <- match.arg(effect)
  
  # Read fit
  fit <- readRDS(file)
  
  # Species name from parent folder
  species <- basename(dirname(file))
  
  # Species-specific scaling values
  sst_sd  <- model_inputs[[species]]$scales$zS_sd
  time_sd <- model_inputs[[species]]$scales$dT_sd
  time_mean <- model_inputs[[species]]$scales$dT_mean
  
  # --------------------------------------------------
  # Parameters needed
  # --------------------------------------------------
  
  if (effect == "temperature") {
    
    vars <- c(
      "beta_mu_zS",
      "beta_alpha_zS"
    )
    
  } else if (effect == "time") {
    
    vars <- c(
      "beta_mu_dT",
      "beta_alpha_dT"
    )
    
    # If looking at time effects inside MPAs,
    # interaction terms are also needed
    if (mpa_ref == 1) {
      vars <- c(
        vars,
        "beta_mu_MPA_dT",
        "beta_alpha_MPA_dT"
      )
    }
    
  } else if (effect == "mpa") {
    
    vars <- c(
      "beta_mu_MPA",
      "beta_alpha_MPA",
      "beta_mu_MPA_dT",
      "beta_alpha_MPA_dT"
    )
  }
  
  # Extract posterior draws
  dr <- fit$draws(
    variables = vars,
    format = "df"
  )
  
  # --------------------------------------------------
  # Construct effective beta coefficients
  # --------------------------------------------------
  
  if (effect == "temperature") {
    
    beta_mu <-dr$beta_mu_zS / sst_sd
    
    beta_alpha <-dr$beta_alpha_zS / sst_sd
    
  } else if (effect == "time") {
    
    # Outside MPA
    if (mpa_ref == 0) {
      
      beta_mu <- dr$beta_mu_dT / time_sd
      
      beta_alpha <- dr$beta_alpha_dT / time_sd
    }
    
    # Inside MPA
    if (mpa_ref == 1) {
      
      beta_mu <- (dr$beta_mu_dT +dr$beta_mu_MPA_dT) / time_sd
      
      beta_alpha <- (dr$beta_alpha_dT +dr$beta_alpha_MPA_dT) / time_sd
    }
    
  } else if (effect == "mpa") {
    
    # If no year supplied, evaluate at centred year
    if (is.null(mpa_year)) {
      
      dT_ref <- 0
      
    } else {
      
      # Convert actual calendar year to scaled dT
      dT_ref <-
        (mpa_year - time_mean) / time_sd
    }
    
    beta_mu <-
      dr$beta_mu_MPA +
      dT_ref * dr$beta_mu_MPA_dT
    
    beta_alpha <-
      dr$beta_alpha_MPA +
      dT_ref * dr$beta_alpha_MPA_dT
  }
  
  # --------------------------------------------------
  # Convert to % change in mean and skewness
  # --------------------------------------------------
  
  mean_effect <-
    100 * (exp(beta_mu) - 1)
  
  skew_effect <-
    100 * (exp(-0.5 * beta_alpha) - 1)
  
  # --------------------------------------------------
  # Posterior median and 80% credible intervals
  # --------------------------------------------------
  
  tibble(
    species = species,
    effect = effect,
    
    mean = median(mean_effect),
    mean_lower = quantile(mean_effect, 0.1),
    mean_upper = quantile(mean_effect, 0.9),
    
    skewness = median(skew_effect),
    skewness_lower = quantile(skew_effect, 0.1),
    skewness_upper = quantile(skew_effect, 0.9), 
    
    sst_sd = sst_sd,
    time_sd = time_sd,
    time_mean = time_mean
  )
}


