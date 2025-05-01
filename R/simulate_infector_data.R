simulate_infector_data <- function(alpha_mean = 0.2, # scale factor on GI
                                   phi_alpha = 0.1,
                                   R0 = 2,
                                   n_infectors = 100,
                                   phi_ind_C =  2,
                                   phi_C =  0.1, # dispersion in daily contacts 
                                   vary_ind_infectiousness = FALSE,
                                   vary_ind_contact_rate = FALSE,
                                   max_gi = 28,
                                   mean_gi = 5,
                                   sigma_gi = 2){
  C_mean <- R0/alpha_mean # Mean number of contacts per
  gi <- dlnorm(x = 1:max_gi, 
               meanlog = convert_to_logmean(mean_gi, sigma_gi), 
               sdlog = convert_to_logsd(mean_gi, sigma_gi))
  
  for(i in 1:n_infectors){
    # Get the individuals number of contacts per time since infection
    if(vary_ind_contact_rate){
      C_i <- rnbinom(n = 1, mu = C_mean, size = phi_ind_C)
    }else{
      C_i <- C_mean # Assume each individual has same mean daily contacts
    }
    C_i_tau <- rnbinom(n = max_gi, mu = C_i, size = phi_C)
    
    # Intrinsic transmissibility of this infectee (AUC of viral load)
    if(vary_ind_infectiousness){
      alpha_i <- rnbinom(n = 1, mu = alpha_mean, size = phi_alpha)
    }else{
      alpha_i <- alpha_mean
    }
    
    
    # Get the number of secondary infections at each time since infection
    N_i_tau <- rbinom(n = max_gi, size = C_i_tau, prob = gi*alpha_i)
    
    df_i <- data.frame(
      infector_id = i,
      n_contacts_tau = C_i_tau,
      tau = 1:max_gi,
      n_pos_contacts_tau = N_i_tau,
      gi = gi,
      # keep the individuals and population metadata 
      ind_mean_num_contacts_ind = C_i,
      pop_mean_num_contacts = C_mean,
      ind_transmissibility = alpha_i,
      pop_transmissibility = alpha_mean,
      # Metadata on this simulation
      vary_ind_infectiousness = vary_ind_infectiousness,
      vary_ind_contact_rate = vary_ind_contact_rate,
      n_infectors = n_infectors,
      R0 = R0
    )
    
    if(i ==1){
      sim_df <- df_i
    }else{
      sim_df <- rbind(sim_df, df_i)
    }
    
  }
  
  return(sim_df)
  
}