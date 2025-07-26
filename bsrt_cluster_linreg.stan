// Linear regression model -- for 2nd stage analysis of clusters

data {
  int<lower=0> N;
  real<lower=0> coef_prior_sd;
  vector[N] bsrt;
  vector[N] age;
  vector[N] sexMale;
  vector[N] cluster1;
  vector[N] cluster2;
  vector[N] cluster3;
}

parameters {
  real intercept;
  real<lower=0> sigma2;
  real beta_age;
  real beta_sexMale;
  real beta_cluster1;
  real beta_cluster2;
  real beta_cluster3;
}

model {
  // linear regression
  bsrt ~ normal(intercept + age*beta_age + sexMale*beta_sexMale + cluster1*beta_cluster1 + cluster2*beta_cluster2 + cluster3*beta_cluster3, sqrt(sigma2));
  
  // priors
  sigma2 ~ inv_gamma(1,1);
  intercept ~ normal(0,coef_prior_sd);
  beta_age ~ normal(0,coef_prior_sd);
  beta_sexMale ~ normal(0,coef_prior_sd);
  beta_cluster1 ~ normal(0,coef_prior_sd);
  beta_cluster2 ~ normal(0,coef_prior_sd);
  beta_cluster3 ~ normal(0,coef_prior_sd);
}

