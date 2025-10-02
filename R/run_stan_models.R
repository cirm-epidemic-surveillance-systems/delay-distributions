run_stan_models <- function(sim_df,
                            max_gi){
  
  C0<- sim_df |>
    pivot_wider(id_cols = infector_id,
                names_from = tau,
                names_prefix = "tau_",
                values_from = n_contacts_tau)|> as.matrix() 
  C <- C0[, 2:(max_gi+1)]
  
  NI0 <-sim_df |>
    pivot_wider(id_cols = infector_id,
                names_from = tau,
                names_prefix = "tau_",
                values_from = n_pos_contacts_tau)|> as.matrix() 
  NI <- NI0[, 2:(max_gi+1)]
  n_infectors <- length(unique(sim_df$infector_id))
  # Binomial fit--------------------------------------------
  stan_data <- list(
    N_infectors = n_infectors,
    max_gi = max_gi,
    C = C,
    N = NI,
    tau_vec = 1:max_gi
  )
  
  
  # Fit the model using contacts
  fit_binomial<- model$sample(
    data = stan_data
  )
  
  all_draws <- fit_binomial$draws()
  logmean_estimate_binomial<- all_draws|>
    spread_draws(logmean_gi) |>
    mutate(draw = `.draw`,
           name = "logmean_gi",
    ) |>
    rename(value = logmean_gi) 
  
  
  logsd_estimate_binomial <- all_draws |>
    spread_draws(logsd_gi) |>
    mutate(draw = `.draw`,
           name = "logsd_gi"
    ) |>
    rename(value = logsd_gi)
  gi_binomial <- all_draws|>
    spread_draws(gi[tau]) |>
    sample_draws(ndraws = 100) |>
    mutate(draw = `.draw`) |>
    mutate(
      name = "gi",
    ) |>
    rename(value = gi) 
  # Fit to only positive -----------------------------------
  transmission_times <- c()
  sim_df_filtered <- sim_df[sim_df$n_pos_contacts_tau>0, ]
  for( j in 1:nrow(sim_df_filtered)){
    transmission_times <- c(transmission_times, 
                            rep(sim_df_filtered$tau[j],
                                sim_df_filtered$n_pos_contacts_tau[j]))
  }
  
  stan_data_po <- list(
    N = length(transmission_times),
    transmission_times = transmission_times,
    max_gi = max_gi,
    tau_vec = 1:max_gi
  )
  # Fit the model using only positives 
  fit_po<- po_model$sample(
    data = stan_data_po
  )
  
  all_draws <- fit_po$draws()
  logmean_estimate_po<- all_draws|>
    spread_draws(logmean_gi) |>
    mutate(draw = `.draw`,
           name = "logmean_gi",
    ) |>
    rename(value = logmean_gi) 
  
  logsd_estimate_po <- all_draws |>
    spread_draws(logsd_gi) |>
    mutate(draw = `.draw`,
           name = "logsd_gi"
    ) |>
    rename(value = logsd_gi)
  
  gi_po<- all_draws|>
    spread_draws(gi[tau]) |>
    sample_draws(ndraws = 100) |>
    mutate(draw = `.draw`) |>
    mutate(
      name = "gi",
    ) |>
    rename(value = gi) 
  return(list(sim_df = sim_df, 
              # Estimates from all contacts (binomial obs model)
              logmean_estimate_binomial = logmean_estimate_binomial, 
              logsd_estimate_binomial = logsd_estimate_binomial, 
              gi_binomial = gi_binomial,
              logmean_median = fit_binomial$summary("logmean_gi")$median,
              logmean_mean = mean(logmean_estimate_binomial$value),
              logmean_lower_95 = as.numeric(quantile(logmean_estimate_binomial$value, 0.025)),
              logmean_upper_95 = as.numeric(quantile(logmean_estimate_binomial$value, 0.975)),
              logsd_median = fit_binomial$summary("logsd_gi")$median,
              # Estimates from positive contact only 
              logmean_estimate_po = logmean_estimate_po, 
              logsd_estimate_po = logsd_estimate_po, 
              gi_po = gi_po,
              logmean_median_po = fit_po$summary("logmean_gi")$median,
              logmean_mean_po = mean(logmean_estimate_po$value),
              logmean_lower_95_po = as.numeric(quantile(logmean_estimate_po$value, 0.025)),
              logmean_upper_95_po = as.numeric(quantile(logmean_estimate_po$value, 0.975)),
              logsd_median_po = fit_po$summary("logsd_gi")$median
  )
  )
}