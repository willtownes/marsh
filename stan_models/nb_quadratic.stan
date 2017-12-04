data {
  int<lower=0> N; //number of total observations
  int<lower=1> D; //dimensionality of covariates, including intercept
  int<lower=1> K; //number of clusters
  matrix[N,D] X; //covariate design matrix
  int<lower=0> y[N]; //nonnegative count outcome
  int<lower=1,upper=K> id[N]; //cluster id
  int<lower=1,upper=2> phi_group[N]; //indicator of salamander species
}
parameters {
  //real intercept;
  vector[D] beta; //fixed effects
  real<lower=0> sigma; //stdev of random intercepts
  vector[K] rand_ints;
  //real logphi[2];
  real<lower=0> phi[2];
}
model {
  vector[N] eta;
  vector[N] rand_ints_vec;
  vector[N] phi_vec;
  //omit beta[1] since it is assumed to be intercept
  for(d in 2:D) beta[d]~cauchy(0,10); 
  phi~gamma(2,.1);
  // prior for random effect stdev
  sigma~cauchy(0,1);
  //prior for random effects
  for(k in 1:K){
    rand_ints[k]~normal(0,sigma);
  }
  for(n in 1:N){
    rand_ints_vec[n]=rand_ints[id[n]];
    phi_vec[n]=phi[phi_group[n]];
  }
  eta= X*beta + rand_ints_vec;
  y ~ neg_binomial_2_log(eta, phi_vec);
}
