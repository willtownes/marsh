data {
  int<lower=1> K; //number of clusters (dates)
  int<lower=0> N1; //number of endemic observations
  int<lower=0> N2; //number of red-backed observations
  int<lower=1,upper=K> id1[N1]; //cluster id for endemic
  int<lower=1,upper=K> id2[N2]; //cluster id for red backed
  int<lower=0> y1[N1]; //nonnegative count outcome for endemic
  int<lower=0> y2[N2]; //nonnegative count outcome for red backed
  vector[N1] soilWC_e;
  vector[N1] temp_e;
  vector[N1] temp2_e;
  vector[N2] soilWC_r;
  vector[N2] temp_r;
  vector[N2] temp2_r;
  //_e means endemic (group1) and _r means red-backed (group2)
  //covariates: intercept,soilWC,temp,temp2. Coefficients: c,w,b,a
}
parameters {
  //real intercept;
  vector[2] c; //intercepts
  real w; //coefficient for soilWC
  vector[2] b; //coefficients for temp
  real<lower=0> neg_a[2]; //negative of coefficients for temp2
  real<lower=0> sigma; //stdev of random intercepts
  vector[K] rand_ints;
  //real logphi[2];
  real<lower=0> phi[2]; //dispersion parameters for each species
}
transformed parameters {
  real<upper=0> a[2]; //coefficients for temp2
  for(j in 1:2){
    a[j] = -neg_a[j];
  }
}
model {
  vector[N1] eta1;
  vector[N2] eta2;
  vector[N1] rand_ints_vec1;
  vector[N2] rand_ints_vec2;
  //no prior for intercepts
  //weak priors for linear terms
  w ~ cauchy(0,2);
  b ~ cauchy(0,2);
  neg_a ~ gamma(2,1); //force negative concavity of quadratic term
  //beta ~ cauchy(0,5);
  //logphi ~ cauchy(0,1);
  //phi~lognormal(0,1);
  phi~gamma(2,.5); //dispersion params weak prior
  sigma~cauchy(0,1); // prior for random effect stdev
  //prior for random effects
  rand_ints ~ normal(0,sigma);
  //for(k in 1:K){
  //  rand_ints[k]~normal(0,sigma);
  //}
  for(n in 1:N1){
    rand_ints_vec1[n]=rand_ints[id1[n]];
  }
  for(n in 1:N2){
    rand_ints_vec2[n]=rand_ints[id2[n]];
  }
  eta1= a[1]*temp2_e + b[1]*temp_e + c[1] + w*soilWC_e + rand_ints_vec1; 
  eta2= a[2]*temp2_r + b[2]*temp_r + c[2] + w*soilWC_r + rand_ints_vec2;
  y1 ~ neg_binomial_2_log(eta1, phi[1]);
  y2 ~ neg_binomial_2_log(eta2, phi[2]);
}
