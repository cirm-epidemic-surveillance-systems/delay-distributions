data {
  int<lower=0> N;                    // Number of transmission pairs
  vector<lower=0>[N] transmission_times; // Vector of transmission times
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
