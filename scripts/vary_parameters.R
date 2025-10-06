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

# Compile stan models
model_file_path <- file.path("inst", "stan", "binomial_obs_model.stan")
model <- cmdstan_model(model_file_path)
model$compile()

po_model_file_path <- file.path("inst", "stan", "transmission_times_only.stan")
po_model <- cmdstan_model(po_model_file_path)
po_model$compile()

#R0 = \sum_{\tau = 1}^{max_gi}C_mean*GI(\tau)*alpha_mean = C_mean*alpha_mean

# Eventually we will want to do a sensitivity analyses where we vary 
# some of these variables and look at the empirical estimates that 
# result and the inferred estimates we get from applying observation process 
# degradation


# Base scenario + varying within-indiviudual contacts over course of days since infection-----------------------
n_infectors <- c(100, 50, 25, 5)
n_replicates <- 10
alpha_mean <- 0.1 #c(0.1, 0.2, 0.4) #P(inf|contact) means that --> Mean number of contacts per person/day (20, 10, 5) = R0/alpha_mean
phi_alpha <- 0.1 # Variation in intrinsic infectiousness across individuals 
vary_ind_infectiousness <- FALSE # Whether individuals vary in their infectiousness
vary_ind_contact_rate <- FALSE # Whether individuals vary in their average number of contacts
mean_gi <- 4 # c(4, 14, 21) 
sigma_gi <- 2 #c(2, 5, 5)
max_gi <- 28 # c(28, 32, 42) #c(14, 20, 36, 60) 
R0 <- 2  #c(2, 4, 8, 16) #2 (flu), 4 (ancestral covid), 8 (omicron), 16 (measles)
phi_ind_C <- 2 #c(0.1, 2, 20) # Overdispersion in average daily contacts across individuals
phi_C <- c(0.1, 2, 20) # dispersion in daily contact within an individual across the days since infection 



n <- 0 
# Iterate through all combinations defined above
for(i in 1:length(n_infectors)){
  n_infectors_n <- n_infectors[i]
  for(j in 1:length(alpha_mean)){
    alpha_mean_n <- alpha_mean[j]
    phi_alpha_n <- phi_alpha
    for(k in 1:length(mean_gi)){
      mean_gi_n <- mean_gi[k]
      max_gi_n <- max_gi[k]
      sigma_gi_n <- sigma_gi[k]
      for(l in 1:length(R0)){
        R0_n <- R0[l]
        for(m in 1:length(phi_ind_C)){
          phi_ind_C_n <- phi_ind_C[m]
          for(o in 1:length(phi_C)){
            phi_C_n <- phi_C[o]
            n <- n + 1
            for(p in 1:n_replicates){
            
            # True GI for this iteration
             C_bar <- R0_n/alpha_mean_n # Mean number of contacts per day
             true_gi <- dlnorm(x = 1:max_gi_n, 
                          meanlog = convert_to_logmean(mean_gi_n, sigma_gi_n), 
                          sdlog = convert_to_logsd(mean_gi_n, sigma_gi_n))
             meanlog <- convert_to_logmean(mean_gi_n, sigma_gi_n)
             sdlog <- convert_to_logsd(mean_gi_n, sigma_gi_n)

            # Simulate an outbreak save infector infectee pairs in a dataframr
            sim_df <- simulate_infector_data(max_gi = max_gi_n,
                                             mean_gi = mean_gi_n,
                                             sigma_gi = sigma_gi_n,
                                             alpha_mean = alpha_mean_n,
                                             phi_alpha = phi_alpha_n,
                                             R0 = R0_n,
                                             phi_ind_C = phi_ind_C_n,
                                             phi_C = phi_C_n,
                                             n_infectors = n_infectors_n,
                                             vary_ind_infectiousness = vary_ind_infectiousness,
                                             vary_ind_contact_rate = vary_ind_contact_rate)
            #Must have some positive contacts
            if(nrow(sim_df[sim_df$n_pos_contacts_tau>0, ]) >0){
              
              # Fit both models using the same outbreak data -- one using all contacts
              # and another using positive only contacts
              mod <- run_stan_models(
                sim_df = sim_df,
                max_gi = max_gi_n)
              
              # Save results in a table alongside the metadata 
              res_n <- data.frame(
                replicate = p,
                n_infectors = n_infectors_n, 
                alpha_mean = alpha_mean_n, 
                phi_alpha = phi_alpha_n,
                true_mean_gi = mean_gi_n, 
                true_logmean_gi = meanlog,
                true_sdlog_gi = sdlog,
                R0 = R0_n,
                phi_ind_C = phi_ind_C_n, 
                phi_C = phi_C_n,
                C_bar = C_bar,
                vary_ind_infectiousness = vary_ind_infectiousness,
                vary_ind_contact_rate = vary_ind_contact_rate,
                logmean_median = mod$logmean_median, 
                logsd_median = mod$logsd_median,
                logmean_mean = mod$logmean_mean, 
                logmean_upper_95 = mod$logmean_upper_95,
                logmean_lower_95 = mod$logmean_lower_95,
                logmean_median_po = mod$logmean_median_po,
                logsd_median_po = mod$logsd_median_po,
                logmean_mean_po = mod$logmean_mean_po,
                logmean_upper_95_po = mod$logmean_upper_95_po,
                logmean_lower_95_po = mod$logmean_lower_95_po)
              
              # Create a plot that compares true GI to estimated GI from both methods
              # plot_gi_comparison(gi_true = gi_true,
              #                    gi_binomial = mod$gi_binomial,
              #                    gi_po = mod$gi_po,
              #                    mean_gi = mean_gi_n,
              #                    n_infectors = n_infectors_n,
              #                    mean_n_contacts_per_day = C_bar,
              #                    R0 = R0_n,
              #                    disp_contacts_day_to_day = phi_C_n,
              #                    disp_contacts_across_ind = phi_ind_C_n,
              #                    disp_ind_infectiousness = phi_alpha,
              #                    output_dir = file.path("results", "figs", "gi_comparisons"))
              if(n == 1){
                res <- res_n
              }else{
                res <- bind_rows(res, res_n)
              }
            }
          }
        }
      }
    }
  }
  }
}

# Base scenario + vary daily contacts between individuals across n_infectors----

n_infectors <- c(100, 50, 25, 5)
n_replicates <- 10
alpha_mean <- 0.1 #c(0.1, 0.2, 0.4) #P(inf|contact) means that --> Mean number of contacts per person/day (20, 10, 5) = R0/alpha_mean
phi_alpha <- 0.1 # Variation in intrinsic infectiousness across individuals 
vary_ind_infectiousness <- FALSE # Whether individuals vary in their infectiousness
vary_ind_contact_rate <- TRUE # Whether individuals vary in their average number of contacts
mean_gi <- 4 # c(4, 14, 21) 
sigma_gi <- 2 #c(2, 5, 5)
max_gi <- 28 # c(28, 32, 42) #c(14, 20, 36, 60) 
R0 <- 2  #c(2, 4, 8, 16) #2 (flu), 4 (ancestral covid), 8 (omicron), 16 (measles)
phi_ind_C <- c(0.1, 2, 20) # Overdispersion in average daily contacts across individuals
phi_C <- 2 # dispersion in daily contact within an individual across the days since infection 



n <- 0 
# Iterate through all combinations defined above
for(i in 1:length(n_infectors)){
  n_infectors_n <- n_infectors[i]
  for(j in 1:length(alpha_mean)){
    alpha_mean_n <- alpha_mean[j]
    phi_alpha_n <- phi_alpha
    for(k in 1:length(mean_gi)){
      mean_gi_n <- mean_gi[k]
      max_gi_n <- max_gi[k]
      sigma_gi_n <- sigma_gi[k]
      for(l in 1:length(R0)){
        R0_n <- R0[l]
        for(m in 1:length(phi_ind_C)){
          phi_ind_C_n <- phi_ind_C[m]
          for(o in 1:length(phi_C)){
            phi_C_n <- phi_C[o]
            n <- n + 1
            for(p in 1:n_replicates){
              
              # True GI for this iteration
              C_bar <- R0_n/alpha_mean_n # Mean number of contacts per day
              true_gi <- dlnorm(x = 1:max_gi_n, 
                                meanlog = convert_to_logmean(mean_gi_n, sigma_gi_n), 
                                sdlog = convert_to_logsd(mean_gi_n, sigma_gi_n))
              meanlog <- convert_to_logmean(mean_gi_n, sigma_gi_n)
              sdlog <- convert_to_logsd(mean_gi_n, sigma_gi_n)
              
              # Simulate an outbreak save infector infectee pairs in a dataframr
              sim_df <- simulate_infector_data(max_gi = max_gi_n,
                                               mean_gi = mean_gi_n,
                                               sigma_gi = sigma_gi_n,
                                               alpha_mean = alpha_mean_n,
                                               phi_alpha = phi_alpha_n,
                                               R0 = R0_n,
                                               phi_ind_C = phi_ind_C_n,
                                               phi_C = phi_C_n,
                                               n_infectors = n_infectors_n,
                                               vary_ind_infectiousness = vary_ind_infectiousness,
                                               vary_ind_contact_rate = vary_ind_contact_rate)
              #Must have some positive contacts
              if(nrow(sim_df[sim_df$n_pos_contacts_tau>0, ]) >0){
                
                # Fit both models using the same outbreak data -- one using all contacts
                # and another using positive only contacts
                mod <- run_stan_models(
                  sim_df = sim_df,
                  max_gi = max_gi_n)
                
                # Save results in a table alongside the metadata 
                res_n <- data.frame(
                  replicate = p,
                  n_infectors = n_infectors_n, 
                  alpha_mean = alpha_mean_n, 
                  phi_alpha = phi_alpha_n,
                  true_mean_gi = mean_gi_n, 
                  true_logmean_gi = meanlog,
                  true_sdlog_gi = sdlog,
                  R0 = R0_n,
                  phi_ind_C = phi_ind_C_n, 
                  phi_C = phi_C_n,
                  C_bar = C_bar,
                  vary_ind_infectiousness = vary_ind_infectiousness,
                  vary_ind_contact_rate = vary_ind_contact_rate,
                  logmean_median = mod$logmean_median, 
                  logsd_median = mod$logsd_median,
                  logmean_mean = mod$logmean_mean, 
                  logmean_upper_95 = mod$logmean_upper_95,
                  logmean_lower_95 = mod$logmean_lower_95,
                  logmean_median_po = mod$logmean_median_po,
                  logsd_median_po = mod$logsd_median_po,
                  logmean_mean_po = mod$logmean_mean_po,
                  logmean_upper_95_po = mod$logmean_upper_95_po,
                  logmean_lower_95_po = mod$logmean_lower_95_po)
                
                if(n == 1){
                  res <- res_n
                }else{
                  res <- bind_rows(res, res_n)
                }
              }
            }
          }
        }
      }
    }
  }
}




write.csv(res, file.path("results", "results_vary_daily_contacts.csv"))

# Base scenario + vary R0-------------------------------------------------------

n_infectors <- c(100, 50, 25, 5)
n_replicates <- 10
alpha_mean <- 0.1 #c(0.1, 0.2, 0.4) #P(inf|contact) means that --> Mean number of contacts per person/day (20, 10, 5) = R0/alpha_mean
phi_alpha <- 0.1 # Variation in intrinsic infectiousness across individuals 
vary_ind_infectiousness <- FALSE # Whether individuals vary in their infectiousness
vary_ind_contact_rate <- FALSE # Whether individuals vary in their average number of contacts
mean_gi <- 4 # c(4, 14, 21) 
sigma_gi <- 2 #c(2, 5, 5)
max_gi <- 28 # c(28, 32, 42) #c(14, 20, 36, 60) 
R0 <- c(2, 4, 8, 16)  #c(2, 4, 8, 16) #2 (flu), 4 (ancestral covid), 8 (omicron), 16 (measles)
phi_ind_C <- 2 # Overdispersion in average daily contacts across individuals
phi_C <- 2 # dispersion in daily contact within an individual across the days since infection 



n <- 0 
# Iterate through all combinations defined above
for(i in 1:length(n_infectors)){
  n_infectors_n <- n_infectors[i]
  for(j in 1:length(alpha_mean)){
    alpha_mean_n <- alpha_mean[j]
    phi_alpha_n <- phi_alpha
    for(k in 1:length(mean_gi)){
      mean_gi_n <- mean_gi[k]
      max_gi_n <- max_gi[k]
      sigma_gi_n <- sigma_gi[k]
      for(l in 1:length(R0)){
        R0_n <- R0[l]
        for(m in 1:length(phi_ind_C)){
          phi_ind_C_n <- phi_ind_C[m]
          for(o in 1:length(phi_C)){
            phi_C_n <- phi_C[o]
            n <- n + 1
            for(p in 1:n_replicates){
              
              # True GI for this iteration
              C_bar <- R0_n/alpha_mean_n # Mean number of contacts per day
              true_gi <- dlnorm(x = 1:max_gi_n, 
                                meanlog = convert_to_logmean(mean_gi_n, sigma_gi_n), 
                                sdlog = convert_to_logsd(mean_gi_n, sigma_gi_n))
              meanlog <- convert_to_logmean(mean_gi_n, sigma_gi_n)
              sdlog <- convert_to_logsd(mean_gi_n, sigma_gi_n)
              
              # Simulate an outbreak save infector infectee pairs in a dataframr
              sim_df <- simulate_infector_data(max_gi = max_gi_n,
                                               mean_gi = mean_gi_n,
                                               sigma_gi = sigma_gi_n,
                                               alpha_mean = alpha_mean_n,
                                               phi_alpha = phi_alpha_n,
                                               R0 = R0_n,
                                               phi_ind_C = phi_ind_C_n,
                                               phi_C = phi_C_n,
                                               n_infectors = n_infectors_n,
                                               vary_ind_infectiousness = vary_ind_infectiousness,
                                               vary_ind_contact_rate = vary_ind_contact_rate)
              #Must have some positive contacts
              if(nrow(sim_df[sim_df$n_pos_contacts_tau>0, ]) >0){
                
                # Fit both models using the same outbreak data -- one using all contacts
                # and another using positive only contacts
                mod <- run_stan_models(
                  sim_df = sim_df,
                  max_gi = max_gi_n)
                
                # Save results in a table alongside the metadata 
                res_n <- data.frame(
                  replicate = p,
                  n_infectors = n_infectors_n, 
                  alpha_mean = alpha_mean_n, 
                  phi_alpha = phi_alpha_n,
                  true_mean_gi = mean_gi_n, 
                  true_logmean_gi = meanlog,
                  true_sdlog_gi = sdlog,
                  R0 = R0_n,
                  phi_ind_C = phi_ind_C_n, 
                  phi_C = phi_C_n,
                  C_bar = C_bar,
                  vary_ind_infectiousness = vary_ind_infectiousness,
                  vary_ind_contact_rate = vary_ind_contact_rate,
                  logmean_median = mod$logmean_median, 
                  logsd_median = mod$logsd_median,
                  logmean_mean = mod$logmean_mean, 
                  logmean_upper_95 = mod$logmean_upper_95,
                  logmean_lower_95 = mod$logmean_lower_95,
                  logmean_median_po = mod$logmean_median_po,
                  logsd_median_po = mod$logsd_median_po,
                  logmean_mean_po = mod$logmean_mean_po,
                  logmean_upper_95_po = mod$logmean_upper_95_po,
                  logmean_lower_95_po = mod$logmean_lower_95_po)
                if(n == 1){
                  res <- res_n
                }else{
                  res <- bind_rows(res, res_n)
                }
              }
            }
          }
        }
      }
    }
  }
}




write.csv(res, file.path("results", "results_vary_R0.csv"))



