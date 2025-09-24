functions{
/**
  * Calculate lognormal density at a vector of indices
  *
  * @param x Vector of values at which to calculate the density
  * @param mu Mean parameter of the underlying normal distribution
  * @param sigma Standard deviation parameter of the underlying normal distribution
  * @return Vector of lognormal density values
  */
  vector lognormal_density_vector(vector x, real mu, real sigma) {
    int n = num_elements(x);
    vector[n] densities;
    
    for (i in 1:n) {
      // Use the built-in lognormal density function
      densities[i] = lognormal_lpdf(x[i] | mu, sigma);
      
      // Convert from log density to regular density
      densities[i] = exp(densities[i]);
    }
    
    return densities;
  }
  
}
 
data {
  int<lower=0> N;                    // Number of transmission pairs
  vector<lower=0>[N] transmission_times; // Vector of transmission times
  int<lower=0> max_gi;
  vector<lower=0> [max_gi] tau_vec;
}

parameters {
  real logmean_gi;                          // Mean of underlying normal distribution
  real<lower=0> logsd_gi;             // Standard deviation of underlying normal distribution
}

model {
  // Priors
  logmean_gi ~ normal(1, 0.5);             // Prior for log-mean parameter
  logsd_gi ~ normal(0, 0.5);          // Prior for log-sd parameter (truncated at 0)
  
  // Likelihood
  transmission_times ~ lognormal(logmean_gi, logsd_gi);
}

generated quantities {
  vector<lower=0>[max_gi] gi;
  gi = lognormal_density_vector(tau_vec, logmean_gi, logsd_gi);
}
