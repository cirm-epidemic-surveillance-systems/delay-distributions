# This script will be used to run the stan model for various parameter settings
library(ggplot2)
library(glue)
library(dplyr)
library(fitdistrplus)
library(cmdstanr)
library(tidybayes)
library(tidyverse)
library(patchwork)

list.files(file.path("R"), full.names = TRUE) |>
  walk(source)

# Compile stan model
model_file_path <- file.path("inst", "stan", "binomial_obs_model.stan")
model <- cmdstan_model(model_file_path)
model$compile()

# Set parameter values
max_gi <- 28
sigma_gi <- 2
vary_ind_infectiousness <- FALSE
phi_alpha <- 0.1
vary_ind_contact_rate <- FALSE

gi <- dlnorm(x = 1:max_gi, 
             meanlog = convert_to_logmean(mean_gi, sigma_gi), 
             sdlog = convert_to_logsd(mean_gi, sigma_gi))
meanlog <- convert_to_logmean(mean_gi, sigma_gi)
sdlog <- convert_to_logsd(mean_gi, sigma_gi)
#R0 = \sum_{\tau = 1}^{max_gi}C_mean*GI(\tau)*alpha_mean = C_mean*alpha_mean

# Eventually we will want to do a sensitivity analyses where we vary 
# some of these variables and look at the empirical estimates that 
# result and the inferred estimates we get from applying observation process 
# degradation
run_stan_model <- function(n_infectors, 
                           alpha_mean,
                           mean_gi, 
                           R0, 
                           phi_ind_C,
                           phi_C,
                           max_gi){

    sim_df <- simulate_infector_data(max_gi = max_gi,
                                   mean_gi = mean_gi,
                                   alpha_mean = alpha_mean,
                                   vary_ind_infectiousness = vary_ind_infectiousness,
                                   phi_alpha = phi_alpha,
                                   R0 = R0,
                                   phi_ind_C = phi_ind_C,
                                   vary_ind_contact_rate = vary_ind_contact_rate,
                                   phi_C = phi_C,
                                   n_infectors = n_infectors)
  
  
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
  
  # A very simple example of fitting a binomial regression in stan
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
  
  return(list(sim_df = sim_df, 
              logmean_estimate_binomial = logmean_estimate_binomial, 
              logsd_estimate_binomial = logsd_estimate_binomial, 
              logmean_median = fit_binomial$summary("logmean_gi")$median,
              logsd_median = fit_binomial$summary("logsd_gi")$median))
}


n_infectors <- c(100, 50, 25, 5)
alpha_mean <- c(0.1, 0.2, 0.4) # Mean number of contacts per person/day (20, 10, 5)
mean_gi <- 5 #c(2, 4, 14, 21) 
max_gi <- 28 #c(14, 20, 36, 60) 
R0 <- 2  #c(2, 4, 8, 16) #2 (flu), 4 (ancestral covid), 8 (omicron), 16 (measles)
phi_ind_C <- 2 #c(0.1, 2, 20) # Overdispersion in average daily contacts of individuals
phi_C <- 0.1 #c(0.1, 2, 20) # dispersion in daily contact


res <- data.frame(n_infectors = NA, alpha_mean = NA, 
                  mean_gi = NA, R0 = NA,
                  phi_ind_C = NA, phi_C = NA,
                  logmean_median = NA, logsd_median = NA)
n <- 0 
for(i in 1:length(n_infectors)){
  for(j in 1:length(alpha_mean)){
    for(k in 1:length(mean_gi)){
      for(l in 1:length(R0)){
        for(m in 1:length(phi_ind_C)){
          for(o in 1:length(phi_C)){
            n <- n + 1
            mod <- run_stan_model(n_infectors = n_infectors[i], 
                                alpha_mean = alpha_mean[j], 
                                mean_gi = mean_gi[k],
                                max_gi = max_gi[k],
                                R0 = R0[l],
                                phi_ind_C = phi_ind_C[m],
                                phi_C = phi_C[o])
            # Save results
            res[n, ] <- c(n_infectors[i], alpha_mean[j], mean_gi[k], 
                          R0 = R0[l], phi_ind_C[m], phi_C[o],
                          mod$logmean_median, mod$logsd_median)
          }
        }
      }
    }
  }
}


results <- res |> mutate(true_logmean = meanlog, true_sdlog = sdlog) 


write.csv(results, file.path("results", "results_all_contacts.csv"))



