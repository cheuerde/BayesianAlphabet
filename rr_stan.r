# Claas Heuer, 2015

library(rstan)

hcode = 
"

data { 

  int<lower=0> Px;
  int<lower=0> Pz;
  int<lower=0> N;
  matrix[N,Px] X;
  matrix[N,Pz] Z;
  real y[N];

}

parameters {

  vector[Px] beta; // flat prior 
  vector[Pz] u; 
// this gives the lower bound for the variance component
  real<lower=0> sigma;
  real<lower=0> tau;

}

model {

// vectorized
    u ~ normal(0,tau);

    y ~ normal(X * beta + Z * u, sigma);

// weakly informative cauchy priors for variance components (see stan manual)
    sigma ~ cauchy(0,5);
    tau ~ cauchy(0,5);

}

"

# some random data
n = 100
p = 20
y <- rnorm(n)
X <- matrix(1,n,1)
Z <- matrix(rnorm(n*p),n,p)


dat = list(N=length(y),Px = ncol(X), Pz = ncol(Z), Z=Z,y=y,X=X)

# run the model
mod = stan(model_name="Ridge Regression", model_code = hcode, data=dat , iter = 10000, warmup = 5000, thin = 1, verbose = FALSE, chains=1)



