
data {
  int<lower=1> N;
  int<lower=1> Nsp;
  int<lower=1> Npop;
  vector[N] y;
  vector[N] year;
  array[Npop] int<lower=1,upper=Nsp> spid_perpop;
  array[N] int<lower=1,upper=Npop> popid;
}

transformed data{
  vector[N] logy = log(y);
}

parameters {
  real mu_alpha1;
  real<lower=0> sigma_alpha1;
  
  real mu_alpha2;
  real<lower=0> sigma_alpha2;
  
  array[Nsp] real mu_tilde_alpha_sp;
  array[Nsp] real<lower=0> sig_tilde_species_alpha;
  
  array[Npop] real alpha_pop_sp;
  
  
  real mu_beta1;
  real<lower=0> sigma_beta1;
  
  real mu_beta2;
  real<lower=0> sigma_beta2;
  
  array[Nsp] real mu_tilde_beta_sp;
  array[Nsp] real<lower=0> sig_tilde_species_beta;
  
  array[Npop] real beta_pop_sp;
  
  real<lower=0> sigma;
   
}

transformed parameters{
  
  array[Nsp] real mu_alpha_sp;
  array[Nsp] real sig_species_alpha;
  
  array[Nsp] real mu_beta_sp;
  array[Nsp] real sig_species_beta;
  
  for (s in 1:Nsp) {
      mu_alpha_sp[s] = mu_alpha1 + sigma_alpha1 * mu_tilde_alpha_sp[s];
      sig_species_alpha[s] = mu_alpha2 + sigma_alpha2 * sig_tilde_species_alpha[s];
      
      mu_beta_sp[s] = mu_beta1 + sigma_beta1 * mu_tilde_beta_sp[s];
      sig_species_beta[s] = mu_beta2 + sigma_beta2 * sig_tilde_species_beta[s];
  }
  
}

model {
  
  vector[N] mu;
  
  for (s in 1:Nsp) {
      mu_tilde_alpha_sp[s] ~ normal(0,1);
      sig_tilde_species_alpha[s] ~ normal(0,1);
      
      mu_tilde_beta_sp[s] ~ normal(0,1);
      sig_tilde_species_beta[s] ~ normal(0,1);
  }
  
  alpha_pop_sp ~ normal(mu_alpha_sp[spid_perpop], sig_species_alpha[spid_perpop]);
  beta_pop_sp ~ normal(mu_beta_sp[spid_perpop], sig_species_beta[spid_perpop]);
  
  for(i in 1:N){
    
    mu[i] = alpha_pop_sp[popid[i]] + beta_pop_sp[popid[i]] * (year[i] - 1980);
    logy[i] ~ normal(mu[i], sigma);
    
  }
  
  mu_alpha1 ~ normal(0, log(500) / 2.57); // -log(500) < mu_alpha1 < log(500)
  sigma_alpha1 ~ normal(0, log(50) / 2.57); // 0 < sigma_alpha1 < log(50)
  
  mu_alpha2 ~ normal(0, log(500) / 2.57); // -log(500) < mu_alpha1 < log(500)
  sigma_alpha2 ~ normal(0, log(50) / 2.57); // 0 < sigma_alpha1 < log(50)
  
  mu_beta1 ~ normal(0, 2 / 2.57); // -2 < mu_beta1 < 2
  sigma_beta1 ~ normal(0, 2 / 2.57); // 0 < sigma_beta1 < 2
  
  mu_beta2 ~ normal(0, 2 / 2.57); //-2 < mu_beta1 < 2
  sigma_beta2 ~ normal(0, 2 / 2.57); // 0 < sigma_beta1 < 2
  
  sigma ~ normal(0, log(100) / 2.57);
  
}

