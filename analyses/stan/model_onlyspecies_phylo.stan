
data {
  int<lower=1> N;
  int<lower=1> Nsp;
  vector[N] y;
  vector[N] x;
  
  array[N] int<lower=1,upper=Nsp> spid;
  
  corr_matrix[Nsp] Cphy; // phylogenetic relationship matrix (fixed)
}

parameters {
  real mu_alpha; // overall species intercept
  real<lower=0, upper=1> lambda_alpha; // phylogenetic structure      
  real<lower=0> sigma_alpha; // overall rate of change (brownian motion?)
  vector[Nsp] mu_alpha_sp; // species-level intercepts
  
  real mu_beta; // overall species trend
  real<lower=0, upper=1> lambda_beta; // phylogenetic structure      
  real<lower=0> sigma_beta; // overall rate of change (brownian motion?)
  vector[Nsp] mu_beta_sp; // species-level trends
  
  real<lower=0> sigma;
   
}

model {
  
  vector[N] mu;
  
  // move all corr/cov matrix to model block, to save RAM
  matrix[Nsp,Nsp] C_alpha = lambda_alpha * Cphy; // previously defined as corr_matrix, but not working in model block?
  C_alpha = C_alpha - diag_matrix(diagonal(C_alpha)) + diag_matrix(diagonal(Cphy));
  matrix[Nsp,Nsp]  C_beta = lambda_beta * Cphy;
  C_beta = C_beta - diag_matrix(diagonal(C_beta)) + diag_matrix(diagonal(Cphy));
  
  // more numerically stable and more efficient to use pre-factored covariance matrices (i.e. multi_normal_cholesky in the following
  matrix[Nsp,Nsp] L_alpha = cholesky_decompose(sigma_alpha^2*C_alpha);
  matrix[Nsp,Nsp] L_beta =  cholesky_decompose(sigma_beta^2*C_beta); 
  
  mu_alpha_sp ~ multi_normal_cholesky(rep_vector(mu_alpha,Nsp), L_alpha); 
  mu_beta_sp ~ multi_normal_cholesky(rep_vector(mu_beta,Nsp), L_beta); 
  
  for(i in 1:N){
    
    mu[i] = mu_alpha_sp[spid[i]] + mu_beta_sp[spid[i]] * x[i];
    y[i] ~ normal(mu[i], sigma);
    
  }
  
  mu_alpha ~ normal(0, 5 / 2.57); // -log(500) < mu_alpha1 < log(500)
  sigma_alpha ~ normal(0, 2 / 2.57); // 0 < sigma_alpha1 < log(100)
  
  mu_beta ~ normal(0, 5 / 2.57); // -2 < mu_beta1 < 2
  sigma_beta ~ normal(0, 2 / 2.57); // 0 < sigma_beta1 < 2

  sigma ~ normal(0, 5 / 2.57);
  
  lambda_alpha ~ beta(1.5, 1.5);
  lambda_beta ~ beta(1.5, 1.5);
}

generated quantities {
  
  vector[N] y_pred;
  
  for(i in 1:N){
    
    real mu = mu_alpha_sp[spid[i]] + mu_beta_sp[spid[i]] * x[i];
    y_pred[i] = normal_rng(mu, sigma);
  }
  
  
  
}

