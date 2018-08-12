# Claas Heuer, 2015

library(rstan)

hcode = 
"
// the data, see 'dat'

data { 

  int<lower=0> Px;
  int<lower=0> Pz;
  int<lower=0> N;
  matrix[N,Px] X;
  matrix[N,Pz] Z;
  real y[N];
  real df_BayesA;

}

parameters {

  vector[Px] beta; // flat prior 

// vector of marker effects
  vector[Pz] u; 
// this gives the lower bound for the variance components
  real<lower=0> sigma;

// this is our vector of marker-specific variances
  real<lower=0> tau[Pz];

// the commom scale parameter. we want to estimate it from data (see: Matthew Stephens et al. 
// 'Polygenic Modeling with Bayesian Sparse Linear Mixed Models', 2013)
  real<lower=0> scale_BayesA;

}

model {

// the likelihood (vector expression)
    y ~ normal(X * beta + Z * u, sigma);

// priors

// stan allows vectorized expressions. 'u' and '0' and 'tau' are vectors
    u ~ normal(0,tau);

// weakly informative cauchy prior for residual variance (see stan manual)
    sigma ~ cauchy(0,5);

// the marker specific variances have scaled inv-chi-square hyperpriors (see Meuwissen et al 2001)
// df_BayesA is fixed here (4 or 5) and scale_BayesA gets estimated from data (hopefully)
    tau ~ scaled_inv_chi_square(df_BayesA, scale_BayesA);

// stan doesnt care about conjugacy, so we put a cauchy prior on scale_BayesA
    scale_BayesA ~ cauchy(0,5);

}

"

# some random data
n = 400
p = 20
y <- rnorm(n)
X <- matrix(1,n,1)
Z <- matrix(rnorm(n*p),n,p)


dat = list(N=length(y),Px = ncol(X), Pz = ncol(Z), Z=Z,y=y,X=X, df_BayesA = 5)

# run the model
mod = stan(model_name="BayesA", model_code = hcode, data=dat , iter = 20000, warmup = 5000, thin = 1, verbose = FALSE, chains=1)



