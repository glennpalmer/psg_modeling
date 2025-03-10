//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> k;
  // Sleep stage data
  int Y[4*N,3];
  matrix[4*N,4] X;
  matrix[4*N,4*N] Z;
  
  // Event count data
  matrix[N,2] Y_seconds;
  int<lower=0> Y_eventcounts[N,2];
}


parameters {
  //////// Markov model /////////
  matrix[4,2] mu_tau;
  matrix[4*N,2] gamma_alpha;
  //<lower=0> Sigma[4,2];
  
  //////// Poisson model ////////
  real<lower=0> lambda_R;
  real<lower=0> lambda_N;
  real phi_R[N];
  real phi_N[N];
  //real<lower=0> psi2_R;
  //real<lower=0> psi2_N;
  
  // Factor model parameters
  matrix[10,k] Lambda;
  matrix[k,N] Eta;
  vector<lower=0>[10] Sigma;
}


transformed parameters {
  matrix[4,3] mu_tau_full = append_col(rep_vector(0,4), mu_tau);
  matrix[4*N,3] gamma_alpha_full = append_col(rep_vector(0,4*N), gamma_alpha);
}


model {
  ////////// Markov model //////////
  
  // priors
  
  // fixed effects
  for (i in 1:4) {
    for (j in 1:2) {
      mu_tau[i,j] ~ normal(0,sqrt(10)); // might need to make these more informative
    }
  }
  
 // random effects
 // for (n in 1:N) {
//    for (i in 1:4) {
//      for (j in 1:2) {
//        gamma_alpha[(n-1)*4 + i, j] ~ normal(0, sqrt(Sigma[i,j]));
//      }
//    }
//  }
  
  // random effects
  for (n in 1:N) {
    append_row(append_row(to_vector(gamma_alpha[((n-1)*4 + 1):(n*4), 1:2]), phi_R[n]), phi_N[n]) ~ multi_normal(Lambda * to_vector(Eta[,n]), diag_matrix(Sigma));
  }
  
  // random effect variances
  //for (i in 1:4) {
  //  for (j in 1:2) {
  //    Sigma[i,j] ~ inv_gamma(1,1);
  //  }
  //}
  
  // random effect variances
  for (i in 1:10) {
    Sigma[i] ~ inv_gamma(1,1);
  }
  for (i in 1:10) {
    for (j in 1:k) {
      Lambda[i,j] ~ normal(0,sqrt(10));
    }
  }
  for (i in 1:N) {
    for (j in 1:k) {
      Eta[j,i] ~ normal(0,1);
    }
  }
  
  // model
  for (i in 1:(4*N)) {
    Y[i] ~ multinomial(softmax(to_vector(X[i]*mu_tau_full + Z[i]*gamma_alpha_full)));
  }
  
  //////////// Poisson model //////////
  // priors
  lambda_R ~ gamma(1,10);
  lambda_N ~ gamma(1,10);
  
  //for (i in 1:N) {
  //  phi_R[i] ~ normal(0, sqrt(psi2_R));
  //  phi_N[i] ~ normal(0, sqrt(psi2_N));
  //}
  
  //psi2_R ~ inv_gamma(1,1);
  //psi2_N ~ inv_gamma(1,1);
  
  // model
  for (i in 1:N) {
    Y_eventcounts[i,1] ~ poisson(Y_seconds[i,1] * lambda_R * exp(phi_R[i]));
    Y_eventcounts[i,2] ~ poisson(Y_seconds[i,2] * lambda_N * exp(phi_N[i]));
  }
}


