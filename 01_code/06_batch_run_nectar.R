# ======================================================================
# RUN STAN MODEL SEPARATELY FOR ALL FISH SPECIES
#
# For each species this script:
#
#   1. Reads the species CSV
#   2. Runs make_stan_data()
#   3. Extracts $stan_data and $init_fun from the returned list
#   4. Fits the model using fit_model()
#   5. Saves the complete CmdStanR fit as an .rds
#   6. Saves the data-preparation output
#   7. Saves model summary
#   8. Marks the species as COMPLETE
#
# The Stan model is compiled ONCE before the species loop.
#
# If the script is stopped and restarted, species that have already
# completed will automatically be skipped.
# ======================================================================


# ----------------------------------------------------------------------
# 1. LOAD PACKAGES -----
# ----------------------------------------------------------------------

library(cmdstanr)
library(tidyverse)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# ----------------------------------------------------------------------
# 2. SOURCE DATA-PREPARATION FUNCTION -----
# ----------------------------------------------------------------------

#source("functions/make_stan_data_func.R")

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

# ----------------------------------------------------------------------
# 3. SET PATHS -----
# ----------------------------------------------------------------------

# Stan model file
stan_file <- "fish_length_NegBinomial_parallelised_MPA.stan"

# Folder containing one CSV file per species
species_data_dir <- "species_data"

# Main folder where all fitted models will be stored
output_root <- "stan_fits"


# Create output directory if it doesn't already exist
dir.create(
  output_root,
  showWarnings = FALSE,
  recursive = TRUE
)


# ----------------------------------------------------------------------
# 4. DEFINE FIT_MODEL() ----
# ----------------------------------------------------------------------

fit_model <- function(stan_data, init_fun, mod, output_dir) {
  
  fit <- mod$sample(
    
    # Stan data for this particular species
    data = stan_data,
    
    # MCMC chains
    chains = 4,
    
    # Run all chains simultaneously
    parallel_chains = 4,
    
    # Number of threads available to reduce_sum() within each chain
    threads_per_chain = 2,
    
    # Warmup iterations per chain
    iter_warmup = 1500,
    
    # Posterior sampling iterations per chain
    iter_sampling = 3000,
    
    # Species-specific initialisation function created by make_stan_data()
    init = init_fun,
    
    output_dir = output_dir,
    
    # Sampling controls
    adapt_delta = 0.99,
    
    max_treedepth = NULL
  )
  
  # Return the CmdStanMCMC fit object
  fit
}


# ----------------------------------------------------------------------
# 5. COMPILE THE STAN MODEL ONCE -----
# ----------------------------------------------------------------------

cat("\nPreparing Stan model...\n")


mod <- cmdstan_model(
  stan_file,
  cpp_options = list(
    stan_threads = TRUE
  ),
  force_recompile = FALSE
)


cat("Stan model ready.\n\n")


# ----------------------------------------------------------------------
# 6. SHOULD EXISTING FITS BE OVERWRITTEN? ----
# ----------------------------------------------------------------------

# FALSE is safest.
#
# If both the fit .rds and COMPLETE.txt exist for a species,
# that species will be skipped.
#
# This means the script can safely be restarted after an interruption.

overwrite_existing <- FALSE


# ----------------------------------------------------------------------
# 7. FIND ALL SPECIES CSV FILES ----
# ----------------------------------------------------------------------

csv_files <- list.files(
  path = species_data_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)


# Stop if no files were found
if (length(csv_files) == 0) {
  
  stop(
    "No CSV files found in: ",
    species_data_dir
  )
}


# ----------------------------------------------------------------------
# 8. GET SPECIES NAMES FROM CSV FILENAMES ---- 
# ----------------------------------------------------------------------

# Example:
#
# species_data/Notolabrus_tetricus.csv
#
# becomes:
#
# Notolabrus_tetricus

species_names <- tools::file_path_sans_ext(
  basename(csv_files)
)


cat(
  "\nFound",
  length(csv_files),
  "species CSV files.\n\n"
)


# ----------------------------------------------------------------------
# 9. SET UP RUN LOG ----
# ----------------------------------------------------------------------

# This keeps track of which species completed, failed, or were skipped.

run_log <- tibble(
  species = character(),
  status = character(),
  start_time = character(),
  end_time = character(),
  message = character()
)


# ----------------------------------------------------------------------
# 10. LOOP THROUGH ALL SPECIES -----
# ----------------------------------------------------------------------

for (i in seq_along(csv_files)) {
  
  
  # Current CSV file
  csv_file <- csv_files[i]
  
  
  # Current species name
  species <- species_names[i]
  
  
  cat("\n")
  cat("============================================================\n")
  
  cat(
    "Species",
    i,
    "of",
    length(csv_files),
    ":",
    species,
    "\n"
  )
  
  cat("============================================================\n")
  
  
  # --------------------------------------------------------------------
  # 11. CREATE A SEPARATE OUTPUT FOLDER FOR THIS SPECIES ----
  # --------------------------------------------------------------------
  
  species_dir <- file.path(
    output_root,
    species
  )
  
  
  dir.create(
    species_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  
  # --------------------------------------------------------------------
  # 12. DEFINE OUTPUT FILES ----
  # --------------------------------------------------------------------
  
  # Full CmdStanR fit object
  fit_file <- file.path(
    species_dir,
    paste0(
      "fit_",
      species,
      ".rds"
    )
  )
  
  
  # Complete output returned by make_stan_data()
  #
  # This may contain:
  #
  # $stan_data
  # $init_fun
  # $scales
  #
  # and any other objects returned by make_stan_data().
  
  model_inputs_file <- file.path(
    species_dir,
    paste0(
      "model_inputs_",
      species,
      ".rds"
    )
  )
  
  
  # Posterior summary
  summary_file <- file.path(
    species_dir,
    paste0(
      "summary_",
      species,
      ".csv"
    )
  )
  
  
  # File indicating successful completion
  complete_file <- file.path(
    species_dir,
    "COMPLETE.txt"
  )
  
  
  # File containing error information if fitting fails
  error_file <- file.path(
    species_dir,
    "ERROR.txt"
  )
  
  
  # --------------------------------------------------------------------
  # 13. CHECK WHETHER THIS SPECIES HAS ALREADY BEEN COMPLETED ----
  # --------------------------------------------------------------------
  
  if (
    !overwrite_existing &&
    file.exists(fit_file) &&
    file.exists(complete_file)
  ) {
    
    
    cat(
      "Already completed - skipping this species.\n"
    )
    
    
    run_log <- bind_rows(
      run_log,
      tibble(
        species = species,
        status = "skipped - already complete",
        start_time = "",
        end_time = "",
        message = ""
      )
    )
    
    
    # Save run log immediately
    write_csv(
      run_log,
      file.path(
        output_root,
        "run_log.csv"
      )
    )
    
    
    next
  }
  
  
  # --------------------------------------------------------------------
  # 14. RECORD START TIME ----
  # --------------------------------------------------------------------
  
  start_time <- Sys.time()
  
  # Temporary CmdStan output directory for this species
  cmdstan_temp_dir <- file.path(
    tempdir(),
    paste0("cmdstan_", species)
  )
  
  # Remove any stale temporary directory from an earlier interrupted run
  if (dir.exists(cmdstan_temp_dir)) {
    unlink(
      cmdstan_temp_dir,
      recursive = TRUE,
      force = TRUE
    )
  }
  
  # Create a fresh temporary directory for this species
  dir.create(
    cmdstan_temp_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  cat(
    "Starting:",
    as.character(start_time),
    "\n"
  )
  
  
  # --------------------------------------------------------------------
  # 15. TRY TO FIT THIS SPECIES ----
  #
  # If one species fails, tryCatch() allows the loop to continue
  # with the next species rather than terminating the whole batch.
  # --------------------------------------------------------------------
  
  tryCatch({
    
    
    # ------------------------------------------------------------------
    # 16. READ THE CSV FILE ----
    # ------------------------------------------------------------------
    
    cat("Reading data...\n")
    
    
    df <- read_csv(
      csv_file,
      show_col_types = FALSE
    )
    
    
    cat(
      "Rows:",
      nrow(df),
      "\n"
    )
    
    
    # ------------------------------------------------------------------
    # 17. PREPARE THE DATA ----
    # ------------------------------------------------------------------
    
    cat("Creating Stan data...\n")
    
    
    # make_stan_data() returns a list.
    #
    # We save the complete returned list as "model_inputs".
    
    model_inputs <- make_stan_data(df)
    
    
    # For example:
    #
    # model_inputs$stan_data
    # model_inputs$init_fun
    # model_inputs$scales
    
    
    # ------------------------------------------------------------------
    # 18. SAVE COMPLETE DATA-PREPARATION OUTPUT ----
    # ------------------------------------------------------------------
    
    # Saving the entire list means we preserve not only the Stan data,
    # but also scaling information and anything else generated by
    # make_stan_data().
    
    saveRDS(
      model_inputs,
      model_inputs_file
    )
    
    
    # ------------------------------------------------------------------
    # 19. FIT THE STAN MODEL ----
    # ------------------------------------------------------------------
    
    cat("Running Stan model...\n")
    
    
    # Pass the already-compiled model object "mod" to fit_model().
    #
    # No call to cmdstan_model() happens here.
    
    fit <- fit_model(
      
      stan_data = model_inputs$stan_data,
      
      init_fun = model_inputs$init_fun,
      
      mod = mod, 
      
      output_dir = cmdstan_temp_dir
    )
    
    
    cat(
      "Stan sampling finished.\n"
    )
    
    
    # ------------------------------------------------------------------
    # 20. CREATE POSTERIOR SUMMARY ----
    # ------------------------------------------------------------------
    
    # Calculate the summary while the original CmdStan output files
    # are definitely still available.
    
    fit_summary <- fit$summary()
    
    
    # ------------------------------------------------------------------
    # 21. SAVE COMPLETE FIT ----
    # ------------------------------------------------------------------
    
    cat(
      "Saving fit object...\n"
    )
    
    # save_object() creates a self-contained RDS containing the
    # posterior draws and diagnostics. The CmdStan CSV files are
    # therefore no longer needed after this succeeds.
    
    fit$save_object(
      file = fit_file
    )
    
    
    # ------------------------------------------------------------------
    # 22. VERIFY FIT WAS SAVED ----
    # ------------------------------------------------------------------
    
    if (!file.exists(fit_file)) {
      
      stop(
        "Fit RDS was not created. CmdStan CSV files will NOT be deleted."
      )
    }
    
    
    # ------------------------------------------------------------------
    # 23. SAVE POSTERIOR SUMMARY ----
    # ------------------------------------------------------------------
    
    write_csv(
      fit_summary,
      summary_file
    )
    
    # ------------------------------------------------------------------
    # DELETE TEMPORARY CMDSTAN FILES ----
    # ------------------------------------------------------------------
    
    cat("Deleting temporary CmdStan files...\n")
    
    unlink(
      cmdstan_temp_dir,
      recursive = TRUE,
      force = TRUE
    )
    
    if (dir.exists(cmdstan_temp_dir)) {
      
      warning(
        "Temporary CmdStan directory could not be deleted: ",
        cmdstan_temp_dir
      )
      
    } else {
      
      cat("Temporary CmdStan files deleted successfully.\n")
    }
    
    # ------------------------------------------------------------------
    # 24. MARK SPECIES AS COMPLETE ----
    # ------------------------------------------------------------------
    
    end_time <- Sys.time()
    
    
    writeLines(
      c(
        
        paste(
          "Species:",
          species
        ),
        
        paste(
          "Input file:",
          csv_file
        ),
        
        paste(
          "Started:",
          start_time
        ),
        
        paste(
          "Completed:",
          end_time
        )
      ),
      
      complete_file
    )
    
    
    # Remove an ERROR.txt left over from an earlier failed attempt
    if (file.exists(error_file)) {
      
      file.remove(
        error_file
      )
    }
    
    
    # ------------------------------------------------------------------
    # 25. UPDATE RUN LOG ----
    # ------------------------------------------------------------------
    
    run_log <- bind_rows(
      run_log,
      
      tibble(
        species = species,
        status = "complete",
        start_time = as.character(start_time),
        end_time = as.character(end_time),
        message = ""
      )
    )
    
    
    write_csv(
      run_log,
      file.path(
        output_root,
        "run_log.csv"
      )
    )
    
    
    # Calculate elapsed time
    elapsed_minutes <- as.numeric(
      difftime(
        end_time,
        start_time,
        units = "mins"
      )
    )
    
    
    cat("\n")
    
    cat(
      "Completed:",
      species,
      "\n"
    )
    
    
    cat(
      "Elapsed time:",
      round(
        elapsed_minutes,
        1
      ),
      "minutes\n"
    )
    
    
    # ------------------------------------------------------------------
    # 26. FREE MEMORY BEFORE THE NEXT SPECIES ----
    # ------------------------------------------------------------------
    
    # IMPORTANT:
    #
    # Do NOT remove "mod".
    #
    # We want the compiled Stan model to remain in memory so that the
    # next species can use it.
    
    rm(
      fit,
      fit_summary,
      model_inputs,
      df
    )
    
    gc()
    
    
    # --------------------------------------------------------------------
    # 27. IF THIS SPECIES FAILS ----
    # --------------------------------------------------------------------
    
  }, error = function(e) {
    
    
    end_time <- Sys.time()
    
    
    error_message <- conditionMessage(e)
    
    # Remove any partial CmdStan output left by the failed species
    if (dir.exists(cmdstan_temp_dir)) {
      
      cat(
        "Removing temporary CmdStan files from failed run...\n"
      )
      
      unlink(
        cmdstan_temp_dir,
        recursive = TRUE,
        force = TRUE
      )
    }
    
    cat("\n")
    
    cat(
      "ERROR fitting:",
      species,
      "\n"
    )
    
    
    cat(
      error_message,
      "\n"
    )
    
    
    # Save error information to the species directory
    writeLines(
      c(
        
        paste(
          "Species:",
          species
        ),
        
        paste(
          "Time:",
          end_time
        ),
        
        paste(
          "Error:",
          error_message
        )
      ),
      
      error_file
    )
    
    
    # Add error to run log
    run_log <<- bind_rows(
      run_log,
      
      tibble(
        species = species,
        status = "ERROR",
        start_time = as.character(start_time),
        end_time = as.character(end_time),
        message = error_message
      )
    )
    
    
    write_csv(
      run_log,
      file.path(
        output_root,
        "run_log.csv"
      )
    )
    
    
    # Clear unused memory before moving to next species
    gc()
    
  })
  
}


# ----------------------------------------------------------------------
# 28. FINISHED ----
# ----------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("ALL SPECIES HAVE BEEN PROCESSED\n")
cat("============================================================\n\n")


print(run_log)