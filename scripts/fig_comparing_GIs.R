# Make a plot comparing the estimates from both methods for 3 GIs

library(ggplot2)
library(glue)
library(dplyr)
library(cmdstanr)
library(tidybayes)

list.files(file.path("R"), full.names = TRUE) |>
  walk(source)

max_gi <- 28
mean_gi <- c(2, 14, 21)
sigma_gi <- c(2, 5, 5)
max_gi <- c(28, 32, 42) #c(14, 20, 36, 60) 
alpha_mean <- 0.2 # scale factor on GI
vary_ind_infectiousness <- FALSE
phi_alpha <- 0.1
R0 <- 2
C_mean <- R0/alpha_mean # Mean number of contacts per day (starting at 10)
phi_ind_C <- 2
vary_ind_contact_rate <- FALSE
phi_C <- 0.1 # dispersion in daily contacts 
n_infectors <- 20
gi <- dlnorm(x = 1:max_gi, 
             meanlog = convert_to_logmean(mean_gi, sigma_gi), 
             sdlog = convert_to_logsd(mean_gi, sigma_gi))
meanlog <- convert_to_logmean(mean_gi, sigma_gi)
sdlog <- convert_to_logsd(mean_gi, sigma_gi)

# Compile stan models
model_file_path <- file.path("inst", "stan", "binomial_obs_model.stan")
model <- cmdstan_model(model_file_path)
model$compile()

po_model_file_path <- file.path("inst", "stan", "transmission_times_only.stan")
po_model <- cmdstan_model(po_model_file_path)
po_model$compile()

# Iterate through all combinations defined above
for(i in 1:length(mean_gi)){
      mean_gi_n <- mean_gi[i]
      max_gi_n <- max_gi[i]
      sigma_gi_n <- sigma_gi[i]
              
      true_gi <- dlnorm(x = 1:max_gi_n, 
                        meanlog = convert_to_logmean(mean_gi_n, sigma_gi_n), 
                        sdlog = convert_to_logsd(mean_gi_n, sigma_gi_n))
      meanlog <- convert_to_logmean(mean_gi_n, sigma_gi_n)
      sdlog <- convert_to_logsd(mean_gi_n, sigma_gi_n)
              
      # Simulate an outbreak save infector infectee pairs in a dataframr
      sim_df <- simulate_infector_data(max_gi = max_gi_n,
                                       mean_gi = mean_gi_n,
                                       sigma_gi = sigma_gi_n,
                                       alpha_mean = alpha_mean,
                                       phi_alpha = phi_alpha,
                                       R0 = R0,
                                       phi_ind_C = phi_ind_C,
                                       phi_C = phi_C,
                                       n_infectors = n_infectors,
                                       vary_ind_infectiousness = vary_ind_infectiousness,
                                       vary_ind_contact_rate = vary_ind_contact_rate)
      #Must have some positive contacts
      if(nrow(sim_df[sim_df$n_pos_contacts_tau>0, ]) >0){
    
    # Fit both models using the same outbreak data -- one using all contacts
    # and another using positive only contacts
    mod <- run_stan_models(
      sim_df = sim_df,
      max_gi = max_gi_n)
      median_binomial <- mod$gi_binomial |>
        group_by(tau) |>
        summarise(median_gi = median(value))
      median_po <- mod$gi_po |>
        group_by(tau) |>
        summarise(median_gi = median(value))
      
      binom_df <- mod$gi_binomial |>
        mutate(estimate_type = "all contacts") |>
        left_join(median_binomial, by = "tau")
      po_df <- mod$gi_po |>
        mutate(estimate_type = "positives only") |>
        left_join(median_po, by = "tau")
      
  
      gi_df <- bind_rows(binom_df, po_df) |>
        mutate(true_mean_gi = mean_gi_n) |>
        left_join(data.frame(tau = 1:max_gi_n,
                             true_gi = true_gi),
                  by = "tau") 
      if(i ==1){
        gi_res <- gi_df
      }else{
        gi_res <- bind_rows(gi_res, gi_df)
      }
    }
                

}

ggplot(gi_res) + 
  geom_line(aes(x = tau, y = value, group = draw, color = estimate_type),
            alpha = 0.1) +
  geom_line(aes(x = tau, y = median_gi, color = estimate_type), linewidth = 1)+
  geom_line(aes(x = tau, y = true_gi), color = "black") +
  theme_bw() +
  theme(legend.title = element_blank()) + 
  ylab("Generation interval PMF") +
  xlab("Time since infection") + 
  facet_wrap(~true_mean_gi, scale = "free", nrow = 3)

