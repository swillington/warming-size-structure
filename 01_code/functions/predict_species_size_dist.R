#get size bins

master_mid <- sort(unique(
  unlist(
    lapply(
      model_inputs,
      function(x) x$stan_data$mid))))

MakeMasterBins <- function(mid) {
  
  mid <- sort(mid)
  
  # Boundaries between adjacent bin midpoints
  internal_boundaries <-
    (mid[-length(mid)] + mid[-1]) / 2
  
  first_boundary <- 0
  
  final_boundary <-
    mid[length(mid)] +
    (mid[length(mid)] - mid[length(mid) - 1]) / 2
  
  boundaries <- c(
    first_boundary,
    internal_boundaries,
    final_boundary )
  
  tibble(
    mid = mid,
    lower = boundaries[-length(boundaries)],
    upper = boundaries[-1],
    width = upper - lower)
}


master_bins <- MakeMasterBins(master_mid)

master_bins


GetSpeciesFitFile <- function(species, fit_files) {
  
  candidate <-
    fit_files[
      basename(dirname(fit_files)) == species]
  
  if (length(candidate) != 1) {
    
    stop(
      "Expected exactly one fit file for ",
      species,
      ", but found ",
      length(candidate))
  }
  candidate
}

MakeCurveScenarios <- function(
    species,
    model_inputs,
    time_years = c(1990, 2000, 2010, 2020),
    temp_change = c(1, 5, 10)) {
  
  inp <- model_inputs[[species]]
  
  zS_mean <- inp$scales$zS_mean
  zS_sd   <- inp$scales$zS_sd
  
  dT_mean <- inp$scales$dT_mean
  dT_sd   <- inp$scales$dT_sd
  
  
  # -----------------------------------------
  # Time effect
  # exact years, inside and outside MPAs
  # SST held at species mean: zS = 0
  # -----------------------------------------
  
  time_scenarios <- tidyr::crossing(
    year = time_years,
    MPA = c(0, 1)) |>
    mutate(
      category = "Time",
      zS = 0,
      dT = (year - dT_mean) / dT_sd,
      MPA_label = if_else(MPA == 0, "Non-MPA", "MPA"),
      scenario_order = match(year, time_years),
      actual_SST = zS_mean,
      actual_time = year,
      temp_label = NA_character_)
  
  
  # -----------------------------------------
  # Temperature effect
  # +1, +5, +10 C
  # time held at species mean: dT = 0
  # show common SST effect only
  # -----------------------------------------
  
  temp_scenarios <- tibble(
    delta_C = temp_change ) |>
    mutate(
      category = "Temperature",
      zS = delta_C / zS_sd,
      dT = 0,
      MPA = NA_real_,
      MPA_label = NA_character_,
      year = NA_real_,
      scenario_order = seq_along(temp_change),
      actual_SST = zS_mean + delta_C,
      actual_time = dT_mean,
      temp_label = case_when(
        delta_C == 0 ~ paste0(
          "Reference (",
          round(zS_mean, 1),
          "°C)" ),
        TRUE ~ paste0("+", delta_C, "°C")))
  
  
  bind_rows(
    time_scenarios,
    temp_scenarios)
}

PredictOneCurveScenario <- function(
    draws,
    scenario,
    x_grid,
    ci_width = 0.80) {
  
  zS <- scenario$zS
  dT <- scenario$dT
  MPA <- scenario$MPA
  category <- scenario$category
  
  alpha_lower <- (1 - ci_width) / 2
  alpha_upper <- 1 - alpha_lower
  
  
  # -----------------------------------------
  # Time panel
  # full time + MPA + MPA:time structure
  # SST fixed at zS = 0
  # -----------------------------------------
  
  if (category == "Time") {
    
    eta_mu <-
      draws$mu_0 +
      draws$beta_mu_dT * dT +
      draws$beta_mu_MPA * MPA +
      draws$beta_mu_MPA_dT * MPA * dT
    
    eta_alpha <-
      draws$alpha_0 +
      draws$beta_alpha_dT * dT +
      draws$beta_alpha_MPA * MPA +
      draws$beta_alpha_MPA_dT * MPA * dT
  }
  
  
  # -----------------------------------------
  # Temperature panel
  # common SST effect only
  # time fixed at dT = 0
  # -----------------------------------------
  
  if (category == "Temperature") {
    
    eta_mu <-
      draws$mu_0 +
      draws$beta_mu_zS * zS
    
    eta_alpha <-
      draws$alpha_0 +
      draws$beta_alpha_zS * zS
  }
  
  
  mu <- exp(eta_mu)
  alpha <- exp(eta_alpha)
  gamma_scale <- mu / alpha
  
  
  density_matrix <- vapply(
    x_grid,
    function(x) {
      dgamma(
        x,
        shape = alpha,
        scale = gamma_scale)
    },
    numeric(length(alpha)))
  
  
  # Central fitted Gamma distribution
  mu_median <- median(mu)
  alpha_median <- median(alpha)
  
  scale_median <- mu_median / alpha_median
  
  central_density <- dgamma(
    x_grid,
    shape = alpha_median,
    scale = scale_median)
  
  curve_df <- tibble(
    length = x_grid,
   # Smooth Gamma curve based on posterior median parameters
    density = central_density,
    # 80% credible interval based on all posterior draws
    density_lower = apply( density_matrix,2, quantile,probs = alpha_lower ),
   density_upper = apply(density_matrix,2,quantile,probs = alpha_upper) ) |> 
    arrange(length)
  
  scenario_info <- scenario |>
    select(
      category,
      year,
      MPA,
      MPA_label,
      temp_label,
      delta_C,
      actual_SST,
      actual_time,
      scenario_order)
  
  bind_cols(
    curve_df,
    scenario_info[rep(1, nrow(curve_df)), ])
}

PredictSpeciesCurves <- function(
    species,
    fit_files,
    model_inputs,
    master_bins,
    time_years = c(1990, 2000, 2010, 2020),
    temp_change = c(0, 1, 5, 10),
    n_draws = 1000,
    grid_step = 0.25,
    seed = 123,
    ci_width = 0.80) {
  
  fit_file <- GetSpeciesFitFile(
    species,
    fit_files)
  
  fit <- readRDS(fit_file)
  
  inp <- model_inputs[[species]]
  
  max_species_mid <- max(inp$stan_data$mid)
  
  species_bins <- master_bins |>
    filter(mid <= max_species_mid)
  
  x_grid <- seq(
    0,
    max(species_bins$upper),
    by = grid_step )
  
  
  pars <- c(
    "mu_0",
    "beta_mu_zS",
    "beta_mu_dT",
    "beta_mu_MPA",
    "beta_mu_MPA_dT",
    "alpha_0",
    "beta_alpha_zS",
    "beta_alpha_dT",
    "beta_alpha_MPA",
    "beta_alpha_MPA_dT" )
  
  draws <- fit$draws(
    variables = pars,
    format = "df" ) |>
    as_tibble()
  
  set.seed(seed)
  
  if (nrow(draws) > n_draws) {
    draws <- draws |>
      slice_sample(n = n_draws)
  }
  
  
  scenarios <- MakeCurveScenarios(
    species = species,
    model_inputs = model_inputs,
    time_years = time_years,
    temp_change = temp_change )
  
  
  curve_df <- map_dfr(
    seq_len(nrow(scenarios)),
    function(i) {
      PredictOneCurveScenario(
        draws = draws,
        scenario = scenarios[i, ],
        x_grid = x_grid,
        ci_width = ci_width)
    } )
  
  
  list(
    species = species,
    curve = curve_df,
    x_max = max(species_bins$upper) )
}