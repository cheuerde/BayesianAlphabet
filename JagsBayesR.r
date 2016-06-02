#' @name JagsBayesR
#' @title Bayesian Regression for genomic prediction using a 4 component mixture of normals
#' on the marker effects
#' @author Claas Heuer
#' @details June 2016
#' @usage 
#' JagsBayesR(y, X = NULL, Z = NULL, niter = 1000, burnin = 500, thin = 1)
#' @description
#' @param y 
#' \code{numeric vector} of responses
#' @param X 
#' Optional design matrix for fixed effects. If omitted, a column vector of
#' ones will be assigned.
#' @param Z 
#' Design Matrix for the random effects, usually a matrix of marker covariates.
#' @param niter 
#' Number of iterations
#' @param burnin 
#' Number of samples to discard before making inferences on the posterior
#' @param thin 
#' Keep every \code{thin} samples after burnin
#' @return \code{list} with two elements: vectors of weights and combined and weighted similarity matrix
#' @references 
#' Moser, Gerhard, Sang Hong Lee, Ben J. Hayes, Michael E. Goddard, Naomi R. Wray, and Peter M. Visscher. 
#' “Simultaneous Discovery, Estimation and Prediction Analysis of Complex Traits Using a Bayesian Mixture Model.” 
#' PLOS Genet 11, no. 4 (July 4, 2015): e1004969. doi:10.1371/journal.pgen.1004969.
#'
#' http://doingbayesiandataanalysis.blogspot.com/2012/06/mixture-of-normal-distributions.html
#' Inspired from: https://darrenjw.wordpress.com/2012/11/20/getting-started-with-bayesian-variable-selection-using-jags-and-rjags/
#' @examples
#'
#' library(BGLR)
#' data(wheat)
#' 
#' mod <- JagsBayesR(y = wheat.Y[,1], Z = wheat.X)

library(rjags)

JagsBayesR <- function(y, X = NULL, Z = NULL, niter = 1000, burnin = 500, thin = 1) {

	# some checks
	if(!(is.numeric(y) | is.numeric(y) | is.vector(y))) stop("y must be a numeric/integer vector")

	if(!is.matrix(Z)) stop("Z must be a numeric matrix")
	if(nrow(Z) != length(y)) stop("Z must have as many rows as elements in y")

	if(is.null(X)) X <- array(1, dim = c(nrow(Z), 1))
	if(!is.matrix(X)) stop("X must be a numeric matrix")
	if(nrow(X) != length(y)) stop("X must have as many rows as elements in y")

	# the data
	N = length(y)
	Nb <- ncol(X)
	Nu <- ncol(Z)

	Nclust <- 4
	clust <- rep(as.numeric(NA), Nu)

	# here the factors with which varA is to be multiplied depending on cluster.
	# NOTE: 0 is an invalid precision parameter for dnorm
	varFac <- c(0.00000000001, 0.0001, 0.001, 0.01)

	dataList = list(
			y = y,
			N = N,
			Nclust = Nclust,
			clust = clust,
			onesRepNclust = rep(1,Nclust),
			Nu = Nu,
			Nb = Nb,
			X = X,
			Z = Z,
			varFac = varFac
			)

	# the Jags model
	modelstring <- "


	model {

		# Likelihood:
		for(i in 1 : N) {

			y[i] ~ dnorm(inprod(X[i,], b) + inprod(Z[i,], u), tauE)

		}

		# prior on marker effects: 4 component mixture of normals
		for (i in 1:Nu) {

			clust[i] ~ dcat(pClust[1:Nclust])
			u[i] ~ dnorm(0, tauTransform[clust[i]])

		}

		# relatively flat prior on fixed effects
		for (i in 1:Nb) {

			b[i] ~ dnorm(0, 0.0001)

		}


		# dirichelt prior on cluster probabilities (uniform)
		pClust[1:Nclust] ~ ddirch(onesRepNclust)

		# variance components
		tauE ~ dgamma(0.01, 0.01)
		varE <- 1/tauE

		tauA ~ dgamma(0.01, 0.01)
		varA <- 1/tauA

		# variance to precision
		for(i in 1:Nclust) {

			tauTransform[i] <- 1 / (varFac[i] * varA)

		}

	}
	"


	# run the model
	load.module("glm")
	model=jags.model(textConnection(modelstring), data=dataList)

	update(model, n.iter = burnin)

	output = coda.samples(model = model,
			      variable.names = c("pClust", "varA", "varE", "b", "u", "clust"),
			      n.iter = niter - burnin, thin = 1)

	return(output)

}
