
# === Prepare data for Stan fitting ============================================
# inputs: df = data frame containing the following columns:
#                location, year, month, sst_mean_annual, size_class, total
# outputs: a 4-element list
#            stan_data: list of processed data to pass to stan for fitting
#            meta:      data frame of raw sample info
#            scales:    means and sd of temperature and time when processing
#            init_fun:  parameter values to pass to stan as starting guesses
make_stan_data <- function(df) {
  # 1. Create year-month index
  df <- df |>
    dplyr::select(location, year, month, sst_mean_annual, size_class, total, MPA) |>
    arrange(location, year, month, MPA, size_class) |>
    mutate(
      ym = paste(sprintf("%04d", year), sprintf("%02d", month), sep = "-"),
      location = factor(location)
    )
  
  # 2. Aggregate
  df_agg <- df |>
    group_by(location, ym, MPA, size_class) |>
    summarise(
      total = sum(total, na.rm = TRUE),
      sst = mean(sst_mean_annual, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 3. Sample index
  df_agg <- df_agg |>
    mutate(sample_id = interaction(location, ym, MPA, drop = TRUE))
  df_agg$sample_id <- factor(  df_agg$sample_id)
  
  # 4. Wide format
  df_wide <- df_agg |>
    pivot_wider(
      id_cols = sample_id,
      names_from = size_class,
      values_from = total,
      values_fill = 0,
      names_sort = TRUE
    ) |>
    arrange(sample_id)
  
  y <- df_wide |>
    select(-sample_id) |>
    as.matrix()
  
  mid <- sort(as.numeric(colnames(y)))
  
  # 5. Remove empty bins
  keep <- colSums(y) > 0
  y <- y[, keep]
  mid <- mid[keep]
  
  # 6. Metadata
  df_meta <- df_agg |>
    group_by(sample_id) |>
    summarise(
      location = first(location),
      ym = first(ym), 
      MPA = first(MPA),
      zS_raw = mean(sst),
      .groups = "drop"
    ) |>
    arrange(sample_id)
  
  df_meta <- df_meta |>
    mutate(r = as.integer(factor(location)))
  
  # 7. Time
  df_meta <- df_meta |>
    separate(ym, into = c("year", "month"), convert = TRUE) |>
    mutate(time = year + (month - 0.5) / 12)
  
  # 8. Standardise predictors
  zS_scaled <- scale(df_meta$zS_raw)
  dT_scaled <- scale(df_meta$time)
  
  zS <- zS_scaled[,1]
  dT <- dT_scaled[,1]
  
  zS_mean <- attr(zS_scaled, "scaled:center")
  zS_sd   <- attr(zS_scaled, "scaled:scale")
  
  dT_mean <- attr(dT_scaled, "scaled:center")
  dT_sd   <- attr(dT_scaled, "scaled:scale")
  
  # 9. Bin widths
  J <- length(mid)
  bin_width <- numeric(J)
  
  for (j in 2:(J-1)) {
    bin_width[j] <- (mid[j+1] - mid[j-1]) / 2
  }
  
  bin_width[1] <- 0.5*(mid[1] + mid[2]) # range: [0, mid-point of first two mids]
  bin_width[J] <- mid[J] - mid[J-1]
  
  log_bin_width <- log(bin_width)
  
  # 10. Totals
  yT <- rowSums(y)
  
  # 11. Initial guess
  # total counts
  total_counts <- sum(y)
  # weighted mean
  mu_hat       <- sum(y * rep(mid, each = nrow(y))) / total_counts
  # expand mid to match matrix shape
  mid_mat      <- matrix(rep(mid, each = nrow(y)), nrow = nrow(y))
  var_hat      <- sum(y * (mid_mat - mu_hat)^2) / total_counts
  alpha_hat    <- mu_hat^2 / var_hat
  # initial guesses for parameters
  mu_0_init    <- log(mu_hat)
  alpha_0_init <- log(alpha_hat)
  
  # recover parameters
  mu_init    <- exp(mu_0_init)
  alpha_init <- exp(alpha_0_init)
  beta_init  <- alpha_init / mu_init
  
  # compute lambda_ij
  N <- nrow(y)
  log_lambda <- matrix(NA, N, J)
  
  for (i in 1:N) {
    for (j in 1:J) {
      log_lambda[i, j] <-
        log(yT[i]) +
        dgamma(mid[j], shape = alpha_init, rate = beta_init, log = TRUE) +
        log_bin_width[j]
    }
  }
  
  lambda <- exp(log_lambda)
  
  # method-of-moments estimate for phi
  num <- sum(lambda^2)
  den <- sum((y - lambda)^2 - lambda)
  
  phi_init <- num / den
  phi_init <- max(phi_init, 1e-2)     # avoid near-zero
  phi_init <- min(phi_init, 100)      # avoid extreme values
  
  init_fun <- function() {
    list(
      mu_0          = mu_0_init,
      alpha_0       = alpha_0_init, 
      beta_mu_zS    = 0.0,
      beta_mu_dT    = 0.0,
      beta_alpha_zS = 0.0,
      beta_alpha_dT = 0.0,
      phi_nb        = phi_init
    )
  }
  
  # Output
  list(
    stan_data = list(
      sample_id = df_wide$sample_id,
      N      = N,
      R      = length(unique(df_meta$r)),
      J      = length(mid),
      mid    = mid,
      log_bw = log_bin_width,
      r      = df_meta$r,
      zS     = zS,
      dT     = dT,
      MPA    = df_meta$MPA,
      n      = y,
      nT     = yT
    ),
    
    meta = df_meta,
    
    scales = list(
      zS_mean = zS_mean,
      zS_sd   = zS_sd,
      dT_mean = dT_mean,
      dT_sd   = dT_sd
    ),
    
    init_fun <- init_fun    
  )
}