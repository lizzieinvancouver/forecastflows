
data {
  int<lower=1> N;
  int<lower=1> Nsp;
  // int<lower=1> Npop;
  vector[N] y;
  vector[N] x;
  // array[Npop] int<lower=1,upper=Nsp> spid_perpop;
  // array[N] int<lower=1,upper=Npop> popid;
  
  array[N] int<lower=1,upper=Nsp> spid;
}

/*transformed data{
  vector[N] logy = log(y);
}*/

parameters {
  real mu_alpha1; // overall species intercept
  real<lower=0> sigma_alpha1; 
  
  // real mu_alpha2;
  // real<lower=0> sigma_alpha2;
  
  array[Nsp] real mu_alpha_sp; // species-level intercepts
  // array[Nsp] real<lower=0> sig_species_alpha;
  
  // array[Npop] real alpha_tilde_pop_sp; // non-centered population intercepts
  
  real mu_beta1; // overall species trend
  real<lower=0> sigma_beta1;
  
  // real mu_beta2;
  // real<lower=0> sigma_beta2;
  
  array[Nsp] real mu_beta_sp; // species-level slopes
  // array[Nsp] real<lower=0> sig_species_beta;
  
  // array[Npop] real beta_tilde_pop_sp; // non-centered population slopes
  
  real<lower=0> sigma;
   
}

// transformed parameters{
  
  // array[Npop] real alpha_pop_sp;
  // array[Npop] real beta_pop_sp;
  
  // for (p in 1:Npop) {
    
      // alpha_pop_sp[p] = mu_alpha_sp[spid_perpop[p]] + sig_species_alpha[spid_perpop[p]] *  alpha_tilde_pop_sp[p];
      // beta_pop_sp[p] = mu_beta_sp[spid_perpop[p]] + sig_species_beta[spid_perpop[p]] *  beta_tilde_pop_sp[p];

  // }
  
// }

model {
  
  vector[N] mu;
  
  mu_alpha_sp ~ normal(mu_alpha1, sigma_alpha1);
  // sig_species_alpha ~ normal(mu_alpha2, sigma_alpha2);
  
  mu_beta_sp ~ normal(mu_beta1, sigma_beta1);
  // sig_species_beta ~ normal(mu_beta2, sigma_beta2);
  
  // for (p in 1:Npop) {
    
    // alpha_tilde_pop_sp[p] ~ normal(0,1);
    // beta_tilde_pop_sp[p] ~ normal(0,1);
    
  // }
  
  for(i in 1:N){
    
    // mu[i] = alpha_pop_sp[popid[i]] + beta_pop_sp[popid[i]] * (year[i] - 1980);
    mu[i] = mu_alpha_sp[spid[i]] + mu_beta_sp[spid[i]] * (x[i]);
    y[i] ~ normal(mu[i], sigma);
    
  }
  
  mu_alpha1 ~ normal(0, 5 / 2.57); // -log(500) < mu_alpha1 < log(500)
  sigma_alpha1 ~ normal(0, 2 / 2.57); // 0 < sigma_alpha1 < log(100)
  
  mu_beta1 ~ normal(0, 5 / 2.57); // -2 < mu_beta1 < 2
  sigma_beta1 ~ normal(0, 2 / 2.57); // 0 < sigma_beta1 < 2
  
  // mu_beta2 ~ normal(0, 2 / 2.57); //-2 < mu_beta1 < 2
  // sigma_beta2 ~ normal(0, 2 / 2.57); // 0 < sigma_beta1 < 2
  
  sigma ~ normal(0, 5 / 2.57);
  
}

generated quantities {
  
  vector[N] y_pred;
  
  for(i in 1:N){
    
    // real mu = alpha_pop_sp[popid[i]] + beta_pop_sp[popid[i]] * (year[i] - 1980);
    real mu = mu_alpha_sp[spid[i]] + mu_beta_sp[spid[i]] * (x[i]);
    y_pred[i] = normal_rng(mu, sigma);
  }
  
  
  
}

