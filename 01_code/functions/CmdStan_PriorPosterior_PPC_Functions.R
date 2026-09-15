# =============================================================================
# CmdStan_PriorPosterior_PPC_Functions.R
#
# Prior-posterior plots and posterior predictive checks for CmdStanR.
#
# This file was adapted from PriorPosteriorFunctions.R and tailored to the
# fish_length_NegBinomial_fixed(2).stan model. The prior parser is still
# reasonably general for scalar real, vector, and matrix parameters.
#
# Required package for plots: ggplot2
# The fitted model must be a CmdStanR CmdStanMCMC object.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Small distribution helpers
# -----------------------------------------------------------------------------

dlaplace <- function(x, mu = 0, b = 1) {
  1 / (2 * b) * exp(-abs(x - mu) / b)
}

plaplace <- function(x, mu = 0, b = 1) {
  ifelse(
    x < mu,
    0.5 * exp((x - mu) / b),
    1 - 0.5 * exp(-(x - mu) / b)
  )
}

qlaplace <- function(p, mu = 0, b = 1) {
  ifelse(
    p < 0.5,
    mu + b * log(2 * p),
    mu - b * log(2 * (1 - p))
  )
}

# -----------------------------------------------------------------------------
# 2. Helpers for reading priors from Stan code
# -----------------------------------------------------------------------------

.strip_comments <- function(code) {
  code <- gsub("//.*", "", code, perl = TRUE)
  code <- gsub("(?s)/\\*.*?\\*/", "", code, perl = TRUE)
  code
}

.extract_block <- function(code, block_name) {
  pattern <- paste0("\\b", block_name, "\\b\\s*\\{")
  m <- regexpr(pattern, code, perl = TRUE)
  if (m[1] == -1) return(NA_character_)

  brace_start <- m[1] + attr(m, "match.length") - 1
  chars <- strsplit(code, "")[[1]]

  depth <- 0
  end <- NA_integer_
  for (i in brace_start:length(chars)) {
    if (chars[i] == "{") depth <- depth + 1
    if (chars[i] == "}") depth <- depth - 1
    if (depth == 0) {
      end <- i
      break
    }
  }

  if (is.na(end)) {
    stop("Unbalanced braces while extracting the '", block_name, "' block.")
  }

  substr(code, brace_start + 1, end - 1)
}

.resolve_numeric <- function(x, data_list = NULL, par_name = NULL) {
  x <- trimws(x)
  val <- suppressWarnings(as.numeric(x))
  if (!is.na(val)) return(val)

  if (!is.null(data_list) && x %in% names(data_list)) {
    candidate <- data_list[[x]]
    if (length(candidate) == 1) {
      val <- suppressWarnings(as.numeric(candidate))
      if (!is.na(val)) return(val)
    }
  }

  warning(
    "Could not resolve value '", x, "'",
    if (!is.null(par_name)) paste0(" for parameter '", par_name, "'") else "",
    ". Supply it through data_list if it is a Stan data constant."
  )
  NA_real_
}

.parse_bounds <- function(bounds_str, data_list = NULL) {
  v_lwr <- NA_real_
  v_upr <- NA_real_

  if (!is.na(bounds_str) && nzchar(bounds_str)) {
    lwr_m <- regmatches(
      bounds_str,
      regexpr("lower\\s*=\\s*[^,>]+", bounds_str, perl = TRUE)
    )
    upr_m <- regmatches(
      bounds_str,
      regexpr("upper\\s*=\\s*[^,>]+", bounds_str, perl = TRUE)
    )

    if (length(lwr_m) > 0) {
      v_lwr <- .resolve_numeric(
        trimws(sub("lower\\s*=\\s*", "", lwr_m)), data_list
      )
    }
    if (length(upr_m) > 0) {
      v_upr <- .resolve_numeric(
        trimws(sub("upper\\s*=\\s*", "", upr_m)), data_list
      )
    }
  }

  list(v_lwr = v_lwr, v_upr = v_upr)
}

.parse_parameters_block <- function(block, data_list = NULL) {
  statements <- strsplit(block, ";")[[1]]
  statements <- trimws(gsub("\\s+", " ", statements))
  statements <- statements[statements != ""]

  results <- list()

  add_result <- function(pname, base, v_lwr, v_upr) {
    results[[length(results) + 1]] <<- data.frame(
      par = pname,
      base = base,
      v_lwr = v_lwr,
      v_upr = v_upr,
      stringsAsFactors = FALSE
    )
  }

  for (stmt in statements) {

    # vector<...>[N] beta
    m_vec <- regmatches(
      stmt,
      regexec(
        "^vector(<([^>]*)>)?\\s*\\[\\s*([^\\]]+)\\s*\\]\\s+(.+)$",
        stmt,
        perl = TRUE
      )
    )[[1]]

    if (length(m_vec) > 0) {
      bounds <- .parse_bounds(m_vec[3], data_list)
      size_expr <- trimws(m_vec[4])
      names_str <- sub("=.*$", "", m_vec[5])
      n_elem <- .resolve_numeric(size_expr, data_list)

      if (is.na(n_elem)) {
        warning("Skipping vector declaration because its dimension could not be resolved: ", stmt)
        next
      }

      n_elem <- as.integer(round(n_elem))
      for (pname in trimws(strsplit(names_str, ",")[[1]])) {
        if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", pname)) next
        for (i in seq_len(n_elem)) {
          add_result(
            paste0(pname, "[", i, "]"),
            pname,
            bounds$v_lwr,
            bounds$v_upr
          )
        }
      }
      next
    }

    # matrix<...>[R,C] beta
    m_mat <- regmatches(
      stmt,
      regexec(
        "^matrix(<([^>]*)>)?\\s*\\[\\s*([^,\\]]+)\\s*,\\s*([^\\]]+)\\s*\\]\\s+(.+)$",
        stmt,
        perl = TRUE
      )
    )[[1]]

    if (length(m_mat) > 0) {
      bounds <- .parse_bounds(m_mat[3], data_list)
      nrow_expr <- trimws(m_mat[4])
      ncol_expr <- trimws(m_mat[5])
      names_str <- sub("=.*$", "", m_mat[6])

      n_row <- .resolve_numeric(nrow_expr, data_list)
      n_col <- .resolve_numeric(ncol_expr, data_list)

      if (is.na(n_row) || is.na(n_col)) {
        warning("Skipping matrix declaration because its dimensions could not be resolved: ", stmt)
        next
      }

      n_row <- as.integer(round(n_row))
      n_col <- as.integer(round(n_col))

      # Stan element names use column-major order.
      for (pname in trimws(strsplit(names_str, ",")[[1]])) {
        if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", pname)) next
        for (j in seq_len(n_col)) {
          for (i in seq_len(n_row)) {
            add_result(
              paste0(pname, "[", i, ",", j, "]"),
              pname,
              bounds$v_lwr,
              bounds$v_upr
            )
          }
        }
      }
      next
    }

    # scalar real<...> a
    m_real <- regmatches(
      stmt,
      regexec("^real(<([^>]*)>)?\\s+(.+)$", stmt, perl = TRUE)
    )[[1]]

    if (length(m_real) > 0) {
      bounds <- .parse_bounds(m_real[3], data_list)
      names_str <- sub("=.*$", "", m_real[4])

      for (pname in trimws(strsplit(names_str, ",")[[1]])) {
        if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", pname)) next
        add_result(pname, pname, bounds$v_lwr, bounds$v_upr)
      }
      next
    }

    # Cholesky factors and other constrained types are intentionally skipped.
    if (grepl(
      "^(array|int|simplex|ordered|positive_ordered|row_vector|cholesky_factor|corr_matrix|cov_matrix|unit_vector)",
      stmt
    )) {
      warning("Unsupported parameter type; skipping: ", stmt)
    }
  }

  if (length(results) == 0) {
    stop("No supported real/vector/matrix parameters were found in the parameters block.")
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}

.parse_priors_block <- function(block, par_names) {
  pattern <- paste0(
    "(?:to_vector\\(\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*\\)|",
    "([A-Za-z_][A-Za-z0-9_]*))",
    "\\s*~\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*\\(([^)]*)\\)"
  )

  m <- gregexpr(pattern, block, perl = TRUE)
  matches <- regmatches(block, m)[[1]]
  priors <- list()

  if (length(matches) == 1 && identical(matches, character(0))) return(priors)

  for (mt in matches) {
    parts <- regmatches(mt, regexec(pattern, mt, perl = TRUE))[[1]]
    pname <- if (nzchar(parts[2])) parts[2] else parts[3]
    dist_stan <- parts[4]
    args_str <- parts[5]

    if (!(pname %in% par_names)) next
    if (pname %in% names(priors)) next

    args <- trimws(strsplit(args_str, ",")[[1]])
    priors[[pname]] <- list(dist_stan = dist_stan, args = args)
  }

  priors
}

.dist_map <- list(
  normal             = list(dist = "normal",      nargs = 2),
  lognormal          = list(dist = "log_normal",  nargs = 2),
  exponential        = list(dist = "exponential", nargs = 1),
  gamma              = list(dist = "gamma",       nargs = 2),
  uniform            = list(dist = "uniform",     nargs = 2),
  double_exponential = list(dist = "laplace",     nargs = 2),
  beta               = list(dist = "beta",        nargs = 2)
)

.prior_density <- function(x, dist, arg1, arg2) {
  switch(
    dist,
    normal      = dnorm(x, mean = arg1, sd = arg2),
    exponential = dexp(x, rate = arg1),
    gamma       = dgamma(x, shape = arg1, rate = arg2),
    log_normal  = dlnorm(x, meanlog = arg1, sdlog = arg2),
    laplace     = dlaplace(x, mu = arg1, b = arg2),
    uniform     = dunif(x, min = arg1, max = arg2),
    beta        = dbeta(x, shape1 = arg1, shape2 = arg2),
    rep(NA_real_, length(x))
  )
}

.prior_cdf <- function(x, dist, arg1, arg2) {
  switch(
    dist,
    normal      = pnorm(x, mean = arg1, sd = arg2),
    exponential = pexp(x, rate = arg1),
    gamma       = pgamma(x, shape = arg1, rate = arg2),
    log_normal  = plnorm(x, meanlog = arg1, sdlog = arg2),
    laplace     = plaplace(x, mu = arg1, b = arg2),
    uniform     = punif(x, min = arg1, max = arg2),
    beta        = pbeta(x, shape1 = arg1, shape2 = arg2),
    NA_real_
  )
}

.prior_quantile <- function(p, dist, arg1, arg2) {
  switch(
    dist,
    normal      = qnorm(p, mean = arg1, sd = arg2),
    exponential = qexp(p, rate = arg1),
    gamma       = qgamma(p, shape = arg1, rate = arg2),
    log_normal  = qlnorm(p, meanlog = arg1, sdlog = arg2),
    laplace     = qlaplace(p, mu = arg1, b = arg2),
    uniform     = qunif(p, min = arg1, max = arg2),
    beta        = qbeta(p, shape1 = arg1, shape2 = arg2),
    NA_real_
  )
}

# Build the prior table required by PriorPosteriorPlot().
# Supports priors written as `par ~ distribution(...)` and
# `to_vector(par) ~ distribution(...)`.
Create_df_priors_cmdstan <- function(stan_code, data_list = NULL) {

  if (is.character(stan_code) && length(stan_code) == 1 && file.exists(stan_code)) {
    code <- paste(readLines(stan_code, warn = FALSE), collapse = "\n")
  } else if (is.character(stan_code)) {
    code <- paste(stan_code, collapse = "\n")
  } else {
    stop("stan_code must be a path to a .stan file or a character string containing Stan code.")
  }

  code <- .strip_comments(code)
  parameters_block <- .extract_block(code, "parameters")
  model_block <- .extract_block(code, "model")

  if (is.na(parameters_block)) stop("No parameters block found in Stan code.")
  if (is.na(model_block)) stop("No model block found in Stan code.")

  df_params <- .parse_parameters_block(parameters_block, data_list)
  priors <- .parse_priors_block(model_block, unique(df_params$base))

  rows <- lapply(seq_len(nrow(df_params)), function(i) {
    pname <- df_params$par[i]
    base <- df_params$base[i]
    v_lwr <- df_params$v_lwr[i]
    v_upr <- df_params$v_upr[i]

    dist <- NA_character_
    arg1 <- NA_real_
    arg2 <- NA_real_

    if (!is.null(priors[[base]])) {
      dist_stan <- priors[[base]]$dist_stan
      args <- priors[[base]]$args

      if (dist_stan %in% names(.dist_map)) {
        mapping <- .dist_map[[dist_stan]]
        dist <- mapping$dist

        if (length(args) >= 1 && nzchar(args[1])) {
          arg1 <- .resolve_numeric(args[1], data_list, base)
        }
        if (mapping$nargs == 2 && length(args) >= 2 && nzchar(args[2])) {
          arg2 <- .resolve_numeric(args[2], data_list, base)
        }
      } else {
        warning(
          "Unsupported prior distribution '", dist_stan,
          "' for parameter '", base, "'."
        )
      }
    }

    # Proper implicit uniform prior when both hard bounds are finite.
    if (is.na(dist) && !is.na(v_lwr) && !is.na(v_upr)) {
      dist <- "uniform"
      arg1 <- v_lwr
      arg2 <- v_upr
    }

    # Plot range: use hard bounds when present and otherwise the central
    # 95% of the *constrained* prior distribution.
    v_min <- v_lwr
    v_max <- v_upr

    if (!is.na(dist) && !is.na(arg1)) {
      lower_cdf <- if (is.na(v_lwr)) 0 else .prior_cdf(v_lwr, dist, arg1, arg2)
      upper_cdf <- if (is.na(v_upr)) 1 else .prior_cdf(v_upr, dist, arg1, arg2)

      if (is.na(v_min)) {
        v_min <- .prior_quantile(
          lower_cdf + 0.025 * (upper_cdf - lower_cdf),
          dist, arg1, arg2
        )
      }
      if (is.na(v_max)) {
        v_max <- .prior_quantile(
          lower_cdf + 0.975 * (upper_cdf - lower_cdf),
          dist, arg1, arg2
        )
      }
    }

    data.frame(
      par = pname,
      v_min = v_min,
      v_max = v_max,
      v_lwr = v_lwr,
      v_upr = v_upr,
      dist = dist,
      arg1 = arg1,
      arg2 = arg2,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Backward-compatible alias using the name from the template file.
Create_df_priors <- Create_df_priors_cmdstan

# -----------------------------------------------------------------------------
# 3. Prior-posterior plotting for CmdStanR
# -----------------------------------------------------------------------------

.check_cmdstan_fit <- function(stan_fit) {
  if (!inherits(stan_fit, "CmdStanMCMC")) {
    stop("stan_fit must be a CmdStanR CmdStanMCMC object.")
  }
}

.expand_parameter_names <- function(pars, all_par_names) {
  pars <- unique(pars)

  out <- unlist(lapply(pars, function(p) {
    if (p %in% all_par_names) return(p)

    matches <- all_par_names[
      grepl(paste0("^", p, "\\["), all_par_names, perl = TRUE)
    ]

    if (length(matches) == 0) {
      warning("Parameter '", p, "' was not found in the parsed prior table; dropping it.")
      return(character(0))
    }

    matches
  }))

  unique(out)
}

.make_prior_curves <- function(df_priors) {
  curves <- lapply(seq_len(nrow(df_priors)), function(i) {
    row <- df_priors[i, ]

    if (is.na(row$dist) || is.na(row$v_min) || is.na(row$v_max)) return(NULL)

    x <- seq(row$v_min, row$v_max, length.out = 250)
    y <- .prior_density(x, row$dist, row$arg1, row$arg2)

    # Normalise for any hard parameter bounds (e.g. half-normal priors for
    # sigma_loc and sigma_loctime). Do NOT renormalise merely because the
    # display range shows only the central 95%.
    lower_prob <- if (is.na(row$v_lwr)) 0 else .prior_cdf(row$v_lwr, row$dist, row$arg1, row$arg2)
    upper_prob <- if (is.na(row$v_upr)) 1 else .prior_cdf(row$v_upr, row$dist, row$arg1, row$arg2)
    normalising_mass <- upper_prob - lower_prob

    if (!is.finite(normalising_mass) || normalising_mass <= 0) {
      stop("Invalid prior normalisation for parameter ", row$par)
    }

    y <- y / normalising_mass

    data.frame(par = as.character(row$par), x = x, y = y)
  })

  do.call(rbind, curves)
}

# Generic prior-posterior plot for a CmdStanR fit.
PriorPosteriorPlot <- function(
    stan_fit,
    stan_code,
    species,
    pars,
    data_list = NULL,
    ncol = NA,
    nbins = 30) {

  .check_cmdstan_fit(stan_fit)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  df_priors <- Create_df_priors_cmdstan(stan_code, data_list = data_list)
  pars_use <- .expand_parameter_names(pars, df_priors$par)

  if (length(pars_use) == 0) stop("No valid parameters remained after matching pars.")

  df_priors_use <- df_priors[df_priors$par %in% pars_use, , drop = FALSE]
  df_priors_use <- df_priors_use[match(pars_use, df_priors_use$par), , drop = FALSE]

  missing_prior <- is.na(df_priors_use$dist)
  if (any(missing_prior)) {
    warning(
      "No supported prior was detected for: ",
      paste(df_priors_use$par[missing_prior], collapse = ", "),
      ". These parameters will not have a prior curve."
    )
  }

  base_pars <- unique(sub("\\[.*$", "", pars_use))
  draws_df <- as.data.frame(stan_fit$draws(variables = base_pars, format = "df"))

  missing_draw_cols <- setdiff(pars_use, names(draws_df))
  if (length(missing_draw_cols) > 0) {
    stop(
      "The following requested posterior columns were not found in the CmdStanR draws: ",
      paste(missing_draw_cols, collapse = ", ")
    )
  }

  posterior_long <- do.call(
    rbind,
    lapply(pars_use, function(p) {
      data.frame(par = p, val = draws_df[[p]], stringsAsFactors = FALSE)
    })
  )

  posterior_long$par <- factor(posterior_long$par, levels = pars_use)

  prior_curves <- .make_prior_curves(df_priors_use)
  if (!is.null(prior_curves) && nrow(prior_curves) > 0) {
    prior_curves$par <- factor(prior_curves$par, levels = pars_use)
  }

  if (is.na(ncol)) ncol <- ceiling(sqrt(length(pars_use)))

  p <- ggplot2::ggplot()

  if (!is.null(prior_curves) && nrow(prior_curves) > 0) {
    p <- p + ggplot2::geom_ribbon(
      data = prior_curves,
      ggplot2::aes(x = x, ymin = 0, ymax = y),
      fill = "steelblue",
      alpha = 0.65
    )
  }
  species_title <- gsub("_", " ", species)
  p +
    ggplot2::geom_histogram(
      data = posterior_long,
      ggplot2::aes(x = val, y = ggplot2::after_stat(density)),
      bins = nbins,
      fill = "grey85",
      colour = "black",
      linewidth = 0.15
    ) +
    ggplot2::facet_wrap(~par, scales = "free", ncol = ncol) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", colour = "black")
    ) +
    ggplot2::labs(title = species_title,
      x = "Parameter value",
      y = "Density",
      subtitle = "Blue = prior; grey histogram = posterior"
    )
}

# Convenience wrapper for the main interpretable parameters in
# fish_length_NegBinomial_fixed(2).stan.
FishLengthPriorPosteriorPlot <- function(
    stan_fit,
    stan_code,
    species,
    data_list = NULL,
    pars = c(
      "mu_0",
      "alpha_0",
      "beta_mu_zS",
      "beta_mu_dT",
      "beta_mu_MPA", 
      "beta_mu_MPA_dT",
      "beta_alpha_zS",
      "beta_alpha_dT",
      "beta_alpha_MPA", 
      "beta_alpha_MPA_dT",
      "sigma_loc",
      "sigma_loctime",
      "phi_nb"
    ),
    ncol = 3,
    nbins = 30) {

  PriorPosteriorPlot(
    stan_fit = stan_fit,
    stan_code = stan_code,
    pars = pars,
    data_list = data_list,
    ncol = ncol,
    nbins = nbins, 
    species = species
  )
}

# -----------------------------------------------------------------------------
# 4. Posterior predictive simulation for fish_length_NegBinomial_fixed(2).stan
# -----------------------------------------------------------------------------

.log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

.validate_fish_length_data <- function(stan_data) {
  required <- c("N", "R", "J", "mid", "log_bw", "r", "zS", "dT", "n", "nT")
  missing <- setdiff(required, names(stan_data))
  if (length(missing) > 0) {
    stop("stan_data is missing: ", paste(missing, collapse = ", "))
  }

  N <- as.integer(stan_data$N)
  J <- as.integer(stan_data$J)

  if (length(stan_data$mid) != J) stop("length(mid) must equal J.")
  if (length(stan_data$log_bw) != J) stop("length(log_bw) must equal J.")
  if (length(stan_data$r) != N) stop("length(r) must equal N.")
  if (length(stan_data$zS) != N) stop("length(zS) must equal N.")
  if (length(stan_data$dT) != N) stop("length(dT) must equal N.")
  if (length(stan_data$nT) != N) stop("length(nT) must equal N.")

  n_obs <- stan_data$n
  if (length(dim(n_obs)) != 2 || !all(dim(n_obs) == c(N, J))) {
    stop("stan_data$n must be an N x J matrix/array.")
  }

  if (any(rowSums(n_obs) != as.numeric(stan_data$nT))) {
    warning(
      "For at least one sample, rowSums(n) does not equal nT. ",
      "The Stan likelihood uses nT to set the expected total count, so check the data preparation."
    )
  }

  invisible(TRUE)
}

# Generate replicated datasets from posterior draws using exactly the same
# likelihood as the Stan model.
#
# This is a CONDITIONAL posterior predictive check for the observations in the
# fitted dataset: it uses posterior draws of loc_re and loctime_re. It therefore
# asks whether the fitted model can reproduce data like the observations at the
# fitted locations and location x time units.
SimulatePosteriorPredictiveCmdStan <- function(
    stan_fit,
    stan_data,
    ndraws = 200,
    seed = 123) {

  .check_cmdstan_fit(stan_fit)
  .validate_fish_length_data(stan_data)

  if (!is.numeric(ndraws) || length(ndraws) != 1 || ndraws < 1) {
    stop("ndraws must be a positive integer.")
  }
  ndraws <- as.integer(ndraws)

  needed <- c(
    "mu_0", "alpha_0",
    "beta_mu_zS", "beta_mu_dT",
    "beta_alpha_zS", "beta_alpha_dT",
    "loc_re", "loctime_re", "phi_nb"
  )

  draws_df <- as.data.frame(stan_fit$draws(variables = needed, format = "df"))
  total_draws <- nrow(draws_df)

  if (ndraws > total_draws) {
    warning("ndraws exceeds available posterior draws; using all draws.")
    ndraws <- total_draws
  }

  set.seed(seed)
  draw_rows <- sort(sample.int(total_draws, ndraws, replace = FALSE))
  draws_use <- draws_df[draw_rows, , drop = FALSE]

  N <- as.integer(stan_data$N)
  J <- as.integer(stan_data$J)
  mid <- as.numeric(stan_data$mid)
  log_bw <- as.numeric(stan_data$log_bw)
  r <- as.integer(stan_data$r)
  zS <- as.numeric(stan_data$zS)
  dT <- as.numeric(stan_data$dT)
  nT <- as.numeric(stan_data$nT)

  # Column names needed for the transformed random effects.
  loc_mu_cols <- paste0("loc_re[1,", r, "]")
  loc_alpha_cols <- paste0("loc_re[2,", r, "]")
  loctime_mu_cols <- paste0("loctime_re[1,", seq_len(N), "]")
  loctime_alpha_cols <- paste0("loctime_re[2,", seq_len(N), "]")

  required_draw_cols <- unique(c(
    "mu_0", "alpha_0",
    "beta_mu_zS", "beta_mu_dT",
    "beta_alpha_zS", "beta_alpha_dT",
    "phi_nb",
    loc_mu_cols, loc_alpha_cols,
    loctime_mu_cols, loctime_alpha_cols
  ))

  missing_cols <- setdiff(required_draw_cols, names(draws_use))
  if (length(missing_cols) > 0) {
    stop(
      "Required posterior variables are missing from the fit: ",
      paste(head(missing_cols, 20), collapse = ", "),
      if (length(missing_cols) > 20) " ..." else ""
    )
  }

  yrep <- array(
    0L,
    dim = c(ndraws, N, J),
    dimnames = list(
      draw = seq_len(ndraws),
      sample = seq_len(N),
      bin = seq_len(J)
    )
  )

  for (d in seq_len(ndraws)) {
    loc_mu <- as.numeric(unlist(draws_use[d, loc_mu_cols, drop = FALSE], use.names = FALSE))
    loc_alpha <- as.numeric(unlist(draws_use[d, loc_alpha_cols, drop = FALSE], use.names = FALSE))
    loctime_mu <- as.numeric(unlist(draws_use[d, loctime_mu_cols, drop = FALSE], use.names = FALSE))
    loctime_alpha <- as.numeric(unlist(draws_use[d, loctime_alpha_cols, drop = FALSE], use.names = FALSE))

    log_mu <-
      draws_use$mu_0[d] +
      draws_use$beta_mu_zS[d] * zS +
      draws_use$beta_mu_dT[d] * dT +
      loc_mu +
      loctime_mu

    log_alpha <-
      draws_use$alpha_0[d] +
      draws_use$beta_alpha_zS[d] * zS +
      draws_use$beta_alpha_dT[d] * dT +
      loc_alpha +
      loctime_alpha

    mu <- exp(log_mu)
    alpha <- exp(log_alpha)
    beta <- alpha / mu
    phi <- draws_use$phi_nb[d]

    for (i in seq_len(N)) {
      log_bin_weight <-
        dgamma(mid, shape = alpha[i], rate = beta[i], log = TRUE) +
        log_bw

      log_weight_sum <- .log_sum_exp(log_bin_weight)
      bin_prob <- exp(log_bin_weight - log_weight_sum)
      lambda <- nT[i] * bin_prob

      # Stan neg_binomial_2(mean=lambda, precision=phi) corresponds to
      # R's rnbinom(size=phi, mu=lambda).
      yrep[d, i, ] <- stats::rnbinom(J, size = phi, mu = lambda)
    }
  }

  out <- list(
    yrep = yrep,
    observed = unname(as.matrix(stan_data$n)),
    nT = nT,
    mid = mid,
    log_bw = log_bw,
    draw_rows = draw_rows,
    seed = seed
  )
  class(out) <- "fish_length_ppc"
  out
}

# -----------------------------------------------------------------------------
# 5. Posterior predictive plots
# -----------------------------------------------------------------------------

.check_ppc_object <- function(ppc) {
  if (!inherits(ppc, "fish_length_ppc")) {
    stop("ppc must be the output from SimulatePosteriorPredictiveCmdStan().")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
}

# Overall length-frequency shape across bins.
# scale = "proportion" is generally most useful because it focuses on the
# allocation of fish among length bins rather than sample size.
PPCPlotBins <- function(
    ppc,
    scale = c("proportion", "mean_count", "sum_count"),
    interval = 0.90) {

  .check_ppc_object(ppc)
  scale <- match.arg(scale)

  if (interval <= 0 || interval >= 1) stop("interval must be between 0 and 1.")
  alpha_q <- (1 - interval) / 2

  yrep <- ppc$yrep
  obs <- ppc$observed
  D <- dim(yrep)[1]
  N <- dim(yrep)[2]
  J <- dim(yrep)[3]

  if (scale == "proportion") {
    obs_den <- rowSums(obs)
    obs_prop <- sweep(obs, 1, obs_den, "/")
    obs_stat <- colMeans(obs_prop, na.rm = TRUE)

    rep_stat <- matrix(NA_real_, nrow = D, ncol = J)
    for (d in seq_len(D)) {
      yy <- yrep[d, , ]
      rep_den <- rowSums(yy)
      keep <- rep_den > 0
      if (any(keep)) {
        rep_stat[d, ] <- colMeans(
          sweep(yy[keep, , drop = FALSE], 1, rep_den[keep], "/"),
          na.rm = TRUE
        )
      }
    }
    ylab <- "Mean proportion per sample"
  }

  if (scale == "mean_count") {
    obs_stat <- colMeans(obs)
    rep_stat <- t(vapply(
      seq_len(D),
      function(d) colMeans(yrep[d, , ]),
      numeric(J)
    ))
    ylab <- "Mean count per sample"
  }

  if (scale == "sum_count") {
    obs_stat <- colSums(obs)
    rep_stat <- t(vapply(
      seq_len(D),
      function(d) colSums(yrep[d, , ]),
      numeric(J)
    ))
    ylab <- "Total count"
  }

  qs <- apply(
    rep_stat,
    2,
    stats::quantile,
    probs = c(alpha_q, 0.5, 1 - alpha_q),
    na.rm = TRUE
  )

  df <- data.frame(
    bin = seq_len(J),
    mid = ppc$mid,
    observed = obs_stat,
    lower = qs[1, ],
    median = qs[2, ],
    upper = qs[3, ]
  )

  ggplot2::ggplot(df, ggplot2::aes(x = mid)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.25
    ) +
    ggplot2::geom_line(ggplot2::aes(y = median), linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(y = observed), size = 1.7) +
    ggplot2::geom_line(ggplot2::aes(y = observed), linetype = "dashed", linewidth = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Length-bin midpoint",
      y = ylab,
      title = "Posterior predictive check across length bins",
      subtitle = paste0(
        round(interval * 100),
        "% posterior predictive interval; points/dashed line = observed"
      )
    )
}

# Check a single sample's observed length-frequency distribution.
PPCPlotSample <- function(ppc, sample_id = 1, interval = 0.90) {
  .check_ppc_object(ppc)

  N <- dim(ppc$yrep)[2]
  if (sample_id < 1 || sample_id > N || sample_id %% 1 != 0) {
    stop("sample_id must be an integer from 1 to N.")
  }
  if (interval <= 0 || interval >= 1) stop("interval must be between 0 and 1.")

  alpha_q <- (1 - interval) / 2
  rep_sample <- ppc$yrep[, sample_id, , drop = FALSE]
  rep_sample <- matrix(rep_sample, nrow = dim(ppc$yrep)[1])

  qs <- apply(
    rep_sample,
    2,
    stats::quantile,
    probs = c(alpha_q, 0.5, 1 - alpha_q),
    na.rm = TRUE
  )

  df <- data.frame(
    mid = ppc$mid,
    observed = ppc$observed[sample_id, ],
    lower = qs[1, ],
    median = qs[2, ],
    upper = qs[3, ]
  )

  ggplot2::ggplot(df, ggplot2::aes(x = mid)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.25
    ) +
    ggplot2::geom_line(ggplot2::aes(y = median), linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(y = observed), size = 1.8) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Length-bin midpoint",
      y = "Count",
      title = paste("Posterior predictive check: sample", sample_id),
      subtitle = paste0(
        round(interval * 100),
        "% posterior predictive interval; points = observed"
      )
    )
}

# Compare observed sample totals with totals in replicated datasets.
# Because nT is explicitly used by the likelihood to set expected counts,
# this is a secondary diagnostic; the bin-shape checks above are more
# informative for this model.
PPCPlotSampleTotals <- function(ppc, interval = 0.90) {
  .check_ppc_object(ppc)

  if (interval <= 0 || interval >= 1) stop("interval must be between 0 and 1.")
  alpha_q <- (1 - interval) / 2

  D <- dim(ppc$yrep)[1]
  N <- dim(ppc$yrep)[2]

  rep_totals <- matrix(NA_real_, nrow = D, ncol = N)
  for (d in seq_len(D)) {
    rep_totals[d, ] <- rowSums(ppc$yrep[d, , ])
  }

  qs <- apply(
    rep_totals,
    2,
    stats::quantile,
    probs = c(alpha_q, 0.5, 1 - alpha_q),
    na.rm = TRUE
  )

  df <- data.frame(
    sample = seq_len(N),
    observed = rowSums(ppc$observed),
    lower = qs[1, ],
    median = qs[2, ],
    upper = qs[3, ]
  )

  ggplot2::ggplot(df, ggplot2::aes(x = observed, y = median)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      width = 0,
      alpha = 0.35
    ) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Observed total count",
      y = "Posterior predictive total count",
      title = "Posterior predictive check of sample totals",
      subtitle = paste0("Points = predictive medians; bars = ", round(interval * 100), "% intervals")
    )
}

# Simple numerical discrepancy checks. Posterior predictive probabilities near
# 0 or 1 indicate that the observed statistic is unusual under the fitted model.
PPCStatistics <- function(ppc) {
  if (!inherits(ppc, "fish_length_ppc")) {
    stop("ppc must be the output from SimulatePosteriorPredictiveCmdStan().")
  }

  obs <- ppc$observed
  yrep <- ppc$yrep
  D <- dim(yrep)[1]

  obs_stats <- c(
    zero_fraction = mean(obs == 0),
    max_cell_count = max(obs),
    variance_cell_count = stats::var(as.numeric(obs))
  )

  rep_stats <- t(vapply(
    seq_len(D),
    function(d) {
      yy <- yrep[d, , ]
      c(
        zero_fraction = mean(yy == 0),
        max_cell_count = max(yy),
        variance_cell_count = stats::var(as.numeric(yy))
      )
    },
    numeric(3)
  ))

  data.frame(
    statistic = names(obs_stats),
    observed = as.numeric(obs_stats),
    predictive_median = apply(rep_stats, 2, stats::median, na.rm = TRUE),
    predictive_q05 = apply(rep_stats, 2, stats::quantile, probs = 0.05, na.rm = TRUE),
    predictive_q95 = apply(rep_stats, 2, stats::quantile, probs = 0.95, na.rm = TRUE),
    posterior_predictive_probability = vapply(
      seq_along(obs_stats),
      function(k) mean(rep_stats[, k] >= obs_stats[k], na.rm = TRUE),
      numeric(1)
    ),
    row.names = NULL
  )
}

# One-call convenience function. It returns the simulation, three standard
# plots, and the numerical discrepancy table. Add sample_id to also get an
# individual-sample PPC plot.
PosteriorPredictiveCheck <- function(
    stan_fit,
    stan_data,
    ndraws = 200,
    seed = 123,
    sample_id = NULL,
    interval = 0.90) {

  ppc <- SimulatePosteriorPredictiveCmdStan(
    stan_fit = stan_fit,
    stan_data = stan_data,
    ndraws = ndraws,
    seed = seed
  )

  out <- list(
    simulation = ppc,
    bins = PPCPlotBins(ppc, scale = "proportion", interval = interval),
    bins_counts = PPCPlotBins(ppc, scale = "mean_count", interval = interval),
    totals = PPCPlotSampleTotals(ppc, interval = interval),
    statistics = PPCStatistics(ppc)
  )

  if (!is.null(sample_id)) {
    out$sample <- PPCPlotSample(ppc, sample_id = sample_id, interval = interval)
  }

  out
}

