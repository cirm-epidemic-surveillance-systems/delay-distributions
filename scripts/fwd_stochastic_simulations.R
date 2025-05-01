# This script will be used to define the parameters we would like to create the 
# number of contacts as a function of time since infection (\tau) and the number
# of positive contacts as a function of time since infection
library(ggplot2)
library(glue)
library(dplyr)

list.files(file.path("R"), full.names = TRUE) |>
  walk(source)

max_gi <- 28
mean_gi <- 5
sigma_gi <- 2
alpha_mean <- 0.2 # scale factor on GI
vary_ind_infectiousness <- FALSE
phi_alpha <- 0.1
R0 <- 2
C_mean <- R0/alpha_mean # Mean number of contacts per
phi_ind_C <- 2
vary_ind_contact_rate <- FALSE
phi_C <- 0.1 # dispersion in daily contacts 
n_infectors <- 100
#R0 = \sum_{\tau = 1}^{max_gi}C_mean*GI(\tau)*alpha_mean = C_mean*alpha_mean

# Eventually we will want to do a sensitivity analyses where we vary 
# some of these variables and look at the empirical estimates that 
# result and the inferred estimates we get from applying observation process 
# degradation
sim_df <- simulate_infector_data(max_gi = max_gi,
                                 mean_gi = mean_gi,
                                 alpha_mean = alpha_mean,
                                 vary_ind_infectiousness = vary_inf_infectiousness,
                                 phi_alpha = phi_alpha,
                                 R0 = R0,
                                 phi_ind_C = phi_ind_C,
                                 vary_ind_contact_rate = vary_ind_contact_rate,
                                 phi_C = phi_C,
                                 n_infectors = n_infectors)

# Make some plots from the simulated data 
ggplot(sim_df) +
  geom_line(aes(x = tau, y = n_pos_contacts_tau, color = as.factor(infector_id), 
                group = as.factor(infector_id)), 
            show.legend = FALSE) +
  xlab('Time since infection') + ylab('Number of infected contacts') +
  theme_bw() + coord_cartesian(xlim = c(0, 15))

# Empirical estimate of the proportion infected among total contacts per day
# since infection
summary_sim_df <- sim_df |>
  group_by(tau) |>
  summarise(n_total_pos_contacts_tau = sum(n_pos_contacts_tau),
             n_total_contacts_tau = sum(n_contacts_tau)) |>
  mutate(gi_from_contacts = to_simplex( n_total_pos_contacts_tau/n_total_contacts_tau),
         gi_ground_truth = gi,
         gi_transmission_pairs = to_simplex(n_total_pos_contacts_tau)) |>
  pivot_longer(cols = starts_with("gi_"),
               names_to = "gi_type",
               values_to = "gi_value")

# Empirical estimates of GI from both methods
ggplot(summary_sim_df)+
  geom_line(aes(x = tau, y = gi_value, color = gi_type))+
  theme_bw() +
  coord_cartesian(xlim = c(0, 15)) +
  xlab('Time since infection') +
  ylab('Infectiousness profile') +
  ggtitle(glue("Empirical GI from {n_infectors} infectors"))

# Distribution of number of contacts across individual-days
ggplot(sim_df) + 
  geom_histogram(aes(x = n_contacts_tau)) +
  theme_bw() +
  xlab("Number of contacts per day") +
  ylab("Frequency") 




