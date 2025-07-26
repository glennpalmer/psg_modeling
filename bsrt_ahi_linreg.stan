// Linear regression model -- for 2nd stage analysis of AHI

data {
  int<lower=0> N;
  real<lower=0> coef_prior_sd;
  vector[N] bsrt;
  vector[N] age;
  vector[N] sexMale;
  vector[N] AHI;
}

parameters {
  real intercept;
  real<lower=0> sigma2;
  real beta_age;
  real beta_sexMale;
  real beta_AHI;
}

model {
  // linear regression
  bsrt ~ normal(intercept + age*beta_age + sexMale*beta_sexMale + AHI*beta_AHI, sqrt(sigma2));
  
  // priors
  sigma2 ~ inv_gamma(1,1);
  intercept ~ normal(0,coef_prior_sd);
  beta_age ~ normal(0,coef_prior_sd);
  beta_sexMale ~ normal(0,coef_prior_sd);
  beta_AHI ~ normal(0,coef_prior_sd);
}

