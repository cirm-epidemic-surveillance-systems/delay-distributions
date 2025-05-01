# This script will be used to define the parameters we would like to create the 
# number of contacts as a function of time since infection (\tau) and the number
# of positive contacts as a function of time since infection
library(ggplot2)
library(glue)
library(dplyr)
library(fitdistrplus)
library(cmdstanr)
library(tidybayes)

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

# Estimate the GI in two ways using a MLE estimate of a parametric fit to a 
# lognormal 

generation_intervals <- c()
sim_df_filtered <- sim_df[sim_df$n_pos_contacts_tau>0, ]
for( j in 1:nrow(sim_df_filtered)){
  generation_intervals <- c(generation_intervals, 
                            rep(sim_df_filtered$tau[j],
                                sim_df_filtered$n_pos_contacts_tau[j]))
}

# Eventually replace with the primary censored dist implementation, but for now
# leave as is
fit_lnorm <- fitdist(generation_intervals, "lnorm")

# Extract parameter estimates from fit to generation times
meanlog_estimate <- fit_lnorm$estimate["meanlog"]
sdlog_estimate <- fit_lnorm$estimate["sdlog"]

df_gi_estimates <-data.frame(
  tau = 1:max_gi,
  gi_ground_truth = gi,
  gi_trad = dlnorm(1:max_gi,
                   meanlog = meanlog_estimate,
                   sdlog = sdlog_estimate)
)|> pivot_longer(cols = starts_with("gi"),
                 values_to = "gi_value",
                 names_to = "estimation_type")

print(glue("True meanlog = {meanlog}, estimated meanlog = {meanlog_estimate}"))
print(glue("True sdlog = {sdlog}, estimated meanlog = {sdlog_estimate}"))

ggplot(df_gi_estimates) +
  geom_line(aes(x = tau, y = gi_value, color = estimation_type)) +
  theme_bw() +
  xlab("Time since infection") +
  ylab("discrete GI")

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

model_file_path <- file.path("inst", "stan", "binomial_obs_model.stan")
model <- cmdstan_model(model_file_path)
model$compile()

# Fit the model
fit_binomial<- model$sample(
  data = stan_data
)

all_draws <- fit_binomial$draws()
logmean_estimate_binomial<- all_draws|>
  spread_draws(logmean_gi) |>
  mutate(draw = `.draw`) |>
  mutate(
    name = "logmean_gi",
  ) |>
  rename(value = logmean_gi) 
ggplot(logmean_estimate_binomial) +
  aes(x = value) +
  stat_halfeye() +
  xlab("Estimated logmean from positive contacts")

gi_binomial<- all_draws|>
  spread_draws(gi[tau]) |>
  sample_draws(ndraws = 100) |>
  mutate(draw = `.draw`) |>
  mutate(
    name = "gi",
  ) |>
  rename(value = gi) 
ggplot(gi_binomial)+
  geom_line(aes(x = tau, y = value, group = draw), alpha = 0.1) +
  xlab("Estimated GI from positive contacts")







