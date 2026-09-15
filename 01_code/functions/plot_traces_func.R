
plot_traces <- function(fit,
                        params = c("beta_mu_zS",    "beta_mu_dT",    
                                   "beta_alpha_zS", "beta_alpha_dT")) {
  
  # extract draws
  draws <- fit$draws()
  
  # convert to array (needed for bayesplot)
  draws_array <- as_draws_array(draws)
  
  mcmc_trace(
    draws_array,
    pars = params,facet_args = list(ncol = 2)) 
    # ggplot2::labs(
    #   title = "MCMC Trace Plots",
    #   subtitle = "Check for mixing and convergence across chains" )
}