#' @name StanBR
#' @title Bayesian Regression models from the "Bayesian Alphabet" implemented in Stan
#' @author Claas Heuer
#' @details April 2016
#' @usage 
#' JagsBR(y, X = NULL, Z, model = "RidgeRegression", niter = 1000, burnin = 500, verbose = FALSE,
#'		   priorscaleE = 10, priorscaleU = 5, chains = 1, ...)
#' @description
#' @param y 
#' \code{numeric vector} of responses
#' @param X 
#' Optional design matrix for fixed effects. If omitted, a column vector of
#' ones will be assigned.
#' @param Z 
#' Design Matrix for the random effects, usually a matrix of marker covariates.
#' @param model 
#' Character string specifying the model to be run. Must be in: \code{RidgeRegression}, 
#' \code{BayesA}, \code{HorseShow}
#' @param niter 
#' Number of iterations
#' @param burnin 
#' Number of samples to discard before making inferences on the posterior
#' @param verbose 
#' Print progress to screen
#' @param priorscaleE 
#' Scale parameter for the cauchy prior on residual variance
#' @param priorscaleU
#' Scale parameter for the cauchy prior on marker effects variance
#' @param chains
#' Number of chains to run
#' @return \code{list} with two elements: vectors of weights and combined and weighted similarity matrix
#' @references 
#' Meuwissen, T., B. J. Hayes, and M. E. Goddard. 
#' “Prediction of Total Genetic Value Using Genome-Wide Dense Marker Maps.” 
#' Genetics 157, no. 4 (2001): 1819–29.
#' 
#' de los Campos, G., H. Naya, D. Gianola, J. Crossa, A. Legarra, E. Manfredi, K. Weigel, and J. M. Cotes. 
#' “Predicting Quantitative Traits With Regression Models for Dense Molecular Markers and Pedigree.” 
#' Genetics 182, no. 1 (May 1, 2009): 375–85. doi:10.1534/genetics.109.101501.
#' 
#' Gianola, D., G. de los Campos, W. G. Hill, E. Manfredi, and R. Fernando. 
#' “Additive Genetic Variability and the Bayesian Alphabet.” 
#' Genetics 183, no. 1 (September 1, 2009): 347–63. doi:10.1534/genetics.109.103952.
#'
#' Habier, David, Rohan L. Fernando, Kadir Kizilkaya, and Dorian J. Garrick. 
#' “Extension of the Bayesian Alphabet for Genomic Selection.” 
#' BMC Bioinformatics 12, no. 1 (2011): 186.
#'
#' Gianola, D. 
#' “Priors in Whole-Genome Regression: The Bayesian Alphabet Returns.” 
#' Genetics 194, no. 3 (July 1, 2013): 573–96. doi:10.1534/genetics.113.151753.
#'
#' Moser, Gerhard, Sang Hong Lee, Ben J. Hayes, Michael E. Goddard, Naomi R. Wray, and Peter M. Visscher. 
#' “Simultaneous Discovery, Estimation and Prediction Analysis of Complex Traits Using a Bayesian Mixture Model.” 
#' PLOS Genet 11, no. 4 (July 4, 2015): e1004969. doi:10.1371/journal.pgen.1004969.
#' 
#' Park, Trevor, and George Casella. 
#' “The Bayesian Lasso.” 
#' Journal of the American Statistical Association 103, no. 482 (June 2008): 681–86. doi:10.1198/016214508000000337.
#'
#' Zhou, Xiang, Peter Carbonetto, and Matthew Stephens. 
#' “Polygenic Modeling with Bayesian Sparse Linear Mixed Models.” 
#' PLoS Genet 9, no. 2 (February 7, 2013): e1003264. doi:10.1371/journal.pgen.1003264.
#' 
#' Gelman, Andrew, Daniel Lee, and Jiqiang Guo. 
#' “Stan A Probabilistic Programming Language for Bayesian Inference and Optimization.” 
#' Journal of Educational and Behavioral Statistics, 2015, 1076998615606113.
#' 
#' Guo, Jiqiang, Daniel Lee, Krzysztof Sakrejda, Jonah Gabry, Ben Goodrich, Joel de Guzman (Boost), 
#' Eric Niebler (Boost), Thomas Heller (Boost), and John Fletcher (Boost). 
#' Rstan: R Interface to Stan (version 2.9.0-3), 2016. https://cran.r-project.org/web/packages/rstan/index.html.
#' 
#' http://andrewgelman.com/2015/02/17/bayesian-survival-analysis-horseshoe-priors/
#' https://groups.google.com/forum/#!topic/stan-users/mCcCg7cpW30
#' https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
#' https://groups.google.com/forum/#!msg/stan-users/4gv3fNCqSNk/VonPUkdcZuUJ
#' https://ariddell.org/horseshoe-prior-with-stan.html
#' @examples
#'
#' library(BGLR)
#' data(wheat)
#' 
#' mod <- StanBR(y = wheat.Y[,1], Z = wheat.X)

library(rstan)
#' @export
StanBR <- function(y, X = NULL, Z, model = "RidgeRegression", niter = 1000, burnin = 500, verbose = FALSE,
		   priorscaleE = 10, priorscaleU = 5, chains = 1, ...) {


	###############################
	### First the model strings ###
	###############################

	modelstringBRR = 
	"

	data { 

		int<lower=0> Px;
		int<lower=0> Pz;
		int<lower=0> N;
		real priorscaleE;
		real priorscaleU;

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
		sigma ~ cauchy(0,priorscaleE);
		tau ~ cauchy(0,priorscaleU);

	}

	"

	# horshoe prior on u
	# https://ariddell.org/horseshoe-prior-with-stan.html
	modelstringHorseShoe = 
	"

	data { 

		int<lower=0> Px;
		int<lower=0> Pz;
		int<lower=0> N;
		real priorscaleE;
		real priorscaleU;

		matrix[N,Px] X;
		matrix[N,Pz] Z;
		real y[N];

	}

	parameters {

		vector[Px] beta; // flat prior 
		vector[Pz] u; 
		vector<lower=0>[Pz] lambda;
		// this gives the lower bound for the variance component
		real<lower=0> sigma;
		real<lower=0> tau;

	}

	model {

		# horseshoe prior
		lambda ~ cauchy(0, 1);
		tau ~ cauchy(0, priorscaleU);

		// vectorized
		u ~ normal(0,lambda * tau);
		y ~ normal(X * beta + Z * u, sigma);

		// weakly informative cauchy priors for variance components (see stan manual)
		sigma ~ cauchy(0,priorscaleE);

	}

	"

	# BayesA
	modelstringBayesA = 
	"

	data { 

		int<lower=0> Px;
		int<lower=0> Pz;
		int<lower=0> N;
		real priorscaleE;
		real priorscaleU;

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
		sigma ~ cauchy(0,priorscaleE);

		// the marker specific variances have scaled inv-chi-square hyperpriors (see Meuwissen et al 2001)
		// df_BayesA is fixed here (4 or 5) and scale_BayesA gets estimated from data (hopefully)
		tau ~ scaled_inv_chi_square(df_BayesA, scale_BayesA);

		// stan doesnt care about conjugacy, so we put a cauchy prior on scale_BayesA
		scale_BayesA ~ cauchy(0, priorscaleU);

	}

	"

	# some checks
	if(!(is.numeric(y) | is.numeric(y) | is.vector(y))) stop("y must be a numeric/integer vector")
	# if(threshold & length(unique(y)) > 2) stop("y must be 0/1 for the threshold model")

	if(!is.matrix(Z)) stop("Z must be a numeric matrix")
	if(nrow(Z) != length(y)) stop("Z must have as many rows as elements in y")

	if(is.null(X)) X <- array(1, dim = c(nrow(Z), 1))
	if(!is.matrix(X)) stop("X must be a numeric matrix")
	if(nrow(X) != length(y)) stop("X must have as many rows as elements in y")

	# check model
	allowedModels <- c("RidgeRegression", "HorseShoe", "BayesA")
	if(!model %in% allowedModels) stop(paste("model must be one of ", paste(allowedModels, collapse = ","), sep = ""))

	if(model == "RidgeRegression") {

		modelstring = modelstringBRR
		model_name = "Ridge Regression"

	}

	if(model == "HorseShoe") {

		modelstring = modelstringHorseShoe
		model_name = "Horse Shoe"

	}

	if(model == "BayesA") {

		modelstring = modelstringBayesA
		model_name = "BayesA"

	}



	# make the data.frame
	N = length(y)
	Pz <- ncol(Z)
	Px <- ncol(X)

	dat = list(N = N, Px = Px, Pz = Pz, Z = Z, y = y, X = X, 
		   priorscaleE = priorscaleE, priorscaleU = priorscaleU, df_BayesA = 5)

	# run the model
	mod = stan(model_name = model_name, model_code = modelstring, data=dat, 
		   iter = niter, warmup = burnin, thin = 1, verbose = verbose, chains = chains, ...)

	return(mod)

}

