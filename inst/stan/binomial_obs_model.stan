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

data{
  int<lower=0> N_infectors;
  int<lower=0> max_gi;
  array [N_infectors, max_gi] int C;
  array [N_infectors, max_gi] int N;
  vector<lower=0> [max_gi] tau_vec;
 }
 
 parameters{
  real logmean_gi;
  real <lower=0>logsd_gi;
  real <lower=0> alpha;

 }
 
 transformed parameters{
   vector<lower=0>[max_gi] p;
   vector<lower=0>[max_gi] gi;
   
   gi = lognormal_density_vector(tau_vec, logmean_gi, logsd_gi);
   
   for(k in 1:max_gi){
        p[k] = gi[k]*alpha;
   }
   
 
 }
 
 model{
   // Priors
   logmean_gi ~ normal(1, 0.5);
   logsd_gi ~ normal(0, 0.5);
   alpha ~ normal(0.1, 0.1); // This should maybe be something else(possibly beta?) but leaving for now
   
   // Likelihood
   for(i in 1:N_infectors){
     for(k in 1:max_gi){
        N[i,k] ~ binomial(C[i,k], p[k]);
     }
   }
 }
 
 