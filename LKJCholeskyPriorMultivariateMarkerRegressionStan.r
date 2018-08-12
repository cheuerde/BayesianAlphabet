# Claas Heuer, July 2016
#
# Stan model for multivariate mixed model with two
# random effects that can be specified by
# their according design matrices
# Help for the LKJ prior from here: http://stla.github.io/stlapblog/posts/StanLKJprior.html

library(pacman)
pacman::p_load(rstan, BGLR, shinystan)

mvntest <- "

data { 

	// number of columns for fixed effects design matrix
	int<lower=0> Px;

	// number of columns for first random effect
	int<lower=0> PzA;

	// number of columns for second random effect
	int<lower=0> PzB;

	// number of observations (rows of Y)
	int<lower=0> N;

	// number of traits
	int<lower=0> Ntraits;

	// the design matrices
	matrix[N,Px] X;
	matrix[N,PzA] ZA;
	matrix[N,PzB] ZB;

	// the matrix of phenotypes
	matrix[N,Ntraits] Y;

	// the prior scale for the cauchy prior
	real priorscaletau;

	// vector[Ntraits] musUA;
	// vector[Ntraits] musUB;

}

transformed data {

	// expectaion vectors for random effects
	vector[Ntraits] musUA;
	vector[Ntraits] musUB;

	for(i in 1:Ntraits) {

		musUA[i] = 0;
		musUB[i] = 0;

	}

}	

parameters {

	// solutuion matrices for fixed and random effects
	matrix[Px,Ntraits] beta;
	matrix[PzA,Ntraits] uA;
	matrix[PzB,Ntraits] uB;

	// the choleskies of the correaltion matrices
	cholesky_factor_corr[Ntraits] OmegaCholE;
	cholesky_factor_corr[Ntraits] OmegaCholA;
	cholesky_factor_corr[Ntraits] OmegaCholB;

	// the standard deviations
	vector<lower=0>[Ntraits] sigmaE;
	vector<lower=0>[Ntraits] sigmaA;
	vector<lower=0>[Ntraits] sigmaB;

}

model {

	// sample random effects from mvn 
	for(i in 1:PzA) uA[i,] ~ multi_normal_cholesky(musUA, diag_pre_multiply(sigmaA, OmegaCholA));

	// sample random effects from mvn 
	for(i in 1:PzB) uB[i,] ~ multi_normal_cholesky(musUB, diag_pre_multiply(sigmaB, OmegaCholB));

	// MVN likelihood
	for(i in 1:N) Y[i,] ~ multi_normal_cholesky(X[i,] * beta + ZA[i,] * uA + ZB[i,] * uB, diag_pre_multiply(sigmaE, OmegaCholE));

	// cauchy prior on the standard deviations
	sigmaE ~ cauchy(0, priorscaletau); 
	sigmaA ~ cauchy(0, priorscaletau); 
	sigmaB ~ cauchy(0, priorscaletau); 

	// LKJ prior on the correlation matrices (1 is uniform in the interval -1 to 1)
	OmegaCholE ~ lkj_corr_cholesky(1); 
	OmegaCholA ~ lkj_corr_cholesky(1); 
	OmegaCholB ~ lkj_corr_cholesky(1); 

}

generated quantities {

	// the correlation matrices
	matrix[Ntraits,Ntraits] OmegaE;
	matrix[Ntraits,Ntraits] OmegaA;
	matrix[Ntraits,Ntraits] OmegaB;

	// the covariance matrices
	matrix[Ntraits,Ntraits] covE;
	matrix[Ntraits,Ntraits] covA;
	matrix[Ntraits,Ntraits] covB;

	// Cor = LL'
	OmegaE = multiply_lower_tri_self_transpose(OmegaCholE);
	OmegaA = multiply_lower_tri_self_transpose(OmegaCholA);
	OmegaB = multiply_lower_tri_self_transpose(OmegaCholB);

	// from correlation- to covariance matrices
	covE = quad_form_diag(OmegaE, sigmaE);
	covA = quad_form_diag(OmegaA, sigmaA);
	covB = quad_form_diag(OmegaB, sigmaB);

}

"

# the model 
modobj = stan_model(model_name="mvntest", model_code = mvntest)

# the BGLR data
data(wheat)
Y <- wheat.Y[,1:4, drop = FALSE]
Ntraits <- ncol(Y)

ZA <- wheat.X[,1:10]
ZB <- wheat.X[,101:110]

N = nrow(Y)

X <- array(1, dim = c(N,1))

Px = ncol(X)
PzA = ncol(ZA)
PzB = ncol(ZB)

priorscaletau <- 5

dat = list(N = N, Px = Px, PzA = PzA, PzB = PzB, Y = Y, X = X, ZA = ZA, ZB = ZB, 
	   priorscaletau = priorscaletau, Ntraits = ncol(Y)) 

# run the model

# needed for running chains in parallel
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# options(mc.cores = 4)

chains = 1
niter = 1000
burnin = 500
mod <- sampling(object = modobj, data=dat , iter = niter, warmup = burnin, thin = 1, verbose = FALSE, chains=chains)

# diagnose
launch_shinystan(mod)
