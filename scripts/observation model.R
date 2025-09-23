## The observation model ##
# The script runs the simulation and adds the observations

library(tidyverse)
library(sigmoid)
library(glue)

list.files(file.path("R"), full.names = TRUE) |>
  walk(source)

# Initiate the simulation
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


# Adding observed data to the simulation 
obs_data <- sim_df %>% group_by(infector_id) %>% 
    mutate(d = round(rgamma(1,15,3))) %>% # days from infection to detection
           group_by(infector_id, tau) %>% 
    mutate(p_contact = case_when(tau < d ~ 1-sigmoid(-5+d-tau), 
                                tau >= d ~ 0), # probability of being recalled contact traced
           p_pos_contact = p_contact,
           rep_n_contacts_tau = case_when(tau < d ~ rbinom(1, n_contacts_tau, p_contact),
                                          tau >= d ~ 0),
           rep_n_pos_contacts_tau = case_when(tau < d ~ rbinom(1, n_pos_contacts_tau, p_contact),
                                          tau >= d ~ 0)) # reported positive number contacts

# Plot of days since infection to detection
obs_data %>% ungroup %>% select(infector_id, d) %>% distinct() %>% 
  ggplot() + geom_histogram(aes(d)) + 
  xlab("Number days since infection") +
  ylab("Frequency") +
  theme_bw()
ggsave("figures/days_to_detection.png")

# Probability of being a recalled contact
data.frame(p = 1-sigmoid(-5+1:28), days = 1:28) %>% 
  ggplot() + 
  geom_point(aes(x = days, y = p)) + 
  xlab("Number days before detection") +
  ylab("Probability of being a recalled contact")+
  theme_bw()
ggsave("figures/p_recalled_contact.png")


# Make some plots from the observation data 
ggplot(obs_data) +
  geom_line(aes(x = tau, y = rep_n_pos_contacts_tau, color = as.factor(infector_id), 
                group = as.factor(infector_id)), 
            show.legend = FALSE) +
  xlab('Time since infection') + ylab('Reported number of infected contacts') +
  theme_bw() + coord_cartesian(xlim = c(0, 15))
ggsave("figures/time_rep_infection.png")


# Distribution of number reported recalled contacts across individual-days
ggplot(obs_data) + 
  geom_histogram(aes(x = rep_n_contacts_tau)) +
  theme_bw() +
  xlab("Number of recalled reported contacts per day") +
  ylab("Number") 
