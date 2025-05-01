# This script will be used to define the parameters we would like to create the 
# number of contacts as a function of time since infection (\tau) and the number
# of positive contacts as a function of time since infection
library(ggplot2)

list.files(file.path("R"), full.names = TRUE) |>
  walk(source)

max_gi <- 28
mean_gi <- 5
sigma_gi <- 2
alpha_mean <- 0.2 # scale factor on GI
R0 <- 2
C_mean <- R0/alpha_mean # Mean number of contacts per
phi_C <- 10 # dispersion in daily contacts 
n_infectors <- 100
#R0 = \sum_{\tau = 1}^{max_gi}C_mean*GI(\tau)*alpha_mean = C_mean*alpha_mean


gi <- dlnorm(x = 1:max_gi, 
             meanlog = convert_to_logmean(mean_gi, sigma_gi), 
             sdlog = convert_to_logsd(mean_gi, sigma_gi))

for(i in 1:n_infectors){
  # Get the individuals number of contacts per time since infection
  C_i <- C_mean # Assume each individual has same mean daily contacts
  C_i_tau <- rnbinom(n = max_gi, mu = C_i, size = phi_C)
  
  # Intrinsic transmissibility of this infectee (AUC of viral load)
  alpha_i <- alpha_mean
  
  # Get the number of secondary infections at each time since infection
  N_i_tau <- rbinom(n = max_gi, size = C_i_tau, prob = gi*alpha_i)
  
  df_i <- data.frame(
    infector = i,
    n_contacts_tau = C_i_tau,
    tau = 1:max_gi,
    n_pos_contacts_tau = N_i_tau,
    # keep the individuals metadata 
    ind_mean_num_contacts_ind = C_i,
    pop_mean_num_contacts = C_mean,
    ind_transmissibility = alpha_i,
    pop_transmissibility = alpha_mean
  )
  
  if(i ==1){
    sim_df <- df_i
  }else{
    sim_df <- rbind(sim_df, df_i)
  }
  
}


# Make some plots from the simulated data 

ggplot


