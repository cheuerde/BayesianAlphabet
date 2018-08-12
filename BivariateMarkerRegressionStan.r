# Claas Heuer, July 2016
#
# Stan model for bivariate mixed model with two
# random effects that can be specified by
# their according design matrices

library(pacman)
pacman::p_load(rstan, BGLR, shinystan)

mvntest <- "

data { 

	int<lower=0> Px;
	int<lower=0> PzA;
	int<lower=0> PzB;
	int<lower=0> N;
	matrix[N,Px] X;
	matrix[N,PzA] ZA;
	matrix[N,PzB] ZB;
	matrix[N,2] Y;

	# the prior scale for the cauchy prior
	real priorscaletau;

}

parameters {

	matrix[Px,2] beta;
	matrix[PzA,2] uA;
	matrix[PzB,2] uB;

	// this gives the lower bound for the variance component
	real<lower=0> tauE1;
	real<lower=0> tauE2;

	real<lower=0> tauA1;
	real<lower=0> tauA2;

	real<lower=0> tauB1;
	real<lower=0> tauB2;

	// here I attempt to put a uniform prior on the correlation coefficients
	// between -1 and 1
	real<lower=-1,upper=1> rhoE;
	real<lower=-1,upper=1> rhoA;
	real<lower=-1,upper=1> rhoB;

}

transformed parameters {

	// correlations
	real sE1;
	real sE2;

	real sA1;
	real sA2;

	real sB1;
	real sB2;

	sE1 = sqrt(tauE1);
	sE2 = sqrt(tauE2);

	sA1 = sqrt(tauA1);
	sA2 = sqrt(tauA2);

	sB1 = sqrt(tauB1);
	sB2 = sqrt(tauB2);

}

model {


	// create a covariance matrix out of the taus and correlations
	matrix[2,2] covE;
	matrix[2,2] covA;
	matrix[2,2] covB;

	vector[2] musUA;
	vector[2] musUB;

	musUA[1] = 0;
	musUA[2] = 0;

	musUB[1] = 0;
	musUB[2] = 0;

	covE[1,1] = tauE1;
	covE[2,2] = tauE2;
	covE[1,2] = rhoE * (sE1 * sE2);
	covE[2,1] = covE[1,2];

	covA[1,1] = tauA1;
	covA[2,2] = tauA2;
	covA[1,2] = rhoA * (sA1 * sA2);
	covA[2,1] = covA[1,2];

	covB[1,1] = tauB1;
	covB[2,2] = tauB2;
	covB[1,2] = rhoB * (sB1 * sB2);
	covB[2,1] = covB[1,2];

	// sample random effects from mvn 
	for(i in 1:PzA) uA[i,] ~ multi_normal(musUA, covA);

	// sample random effects from mvn 
	for(i in 1:PzB) uB[i,] ~ multi_normal(musUB, covB);

	// MVN likelihood
	for(i in 1:N) Y[i,] ~ multi_normal(X[i,] * beta + ZA[i,] * uA + ZB[i,] * uB, covE);


	// weakly informative cauchy priors for variance components (see stan manual)
	tauE1 ~ cauchy(0,priorscaletau);
	tauE2 ~ cauchy(0,priorscaletau);

	# the random effects
	tauA1 ~ cauchy(0,priorscaletau);
	tauA2 ~ cauchy(0,priorscaletau);

	tauB1 ~ cauchy(0,priorscaletau);
	tauB2 ~ cauchy(0,priorscaletau);

}

"

# the model 
modobj = stan_model(model_name="mvntest", model_code = mvntest)

# the BGLR data
data(wheat)
Y <- wheat.Y[,1:2]

# two different subsets of the markers as random effects
ZA <- wheat.X[,1:100]
ZB <- wheat.X[,101:200]

N = nrow(Y)

# add some variance to the scaled phenotypes to see whether our estimates make sense
Y[,1] <- Y[,1] + 5 + rnorm(N,0, 2)
Y[,2] <- Y[,2] + 15 + rnorm(N,0, 3)

# fixed effects. only intercept
X <- array(1, dim = c(N,1))

Px = ncol(X)
PzA = ncol(ZA)
PzB = ncol(ZB)

priorscaletau <- 5

dat = list(N = N, Px = Px, PzA = PzA, PzB = PzB, Y = Y, X = X, ZA = ZA, ZB = ZB, priorscaletau = priorscaletau) 

# run the model
niter = 1000
burnin = 500
mod <- sampling(object = modobj, data=dat , iter = niter, warmup = burnin, thin = 1, verbose = FALSE, chains=1)
