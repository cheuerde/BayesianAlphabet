# Claas Heuer, 2015
# 
# Inspired from: https://darrenjw.wordpress.com/2012/11/20/getting-started-with-bayesian-variable-selection-using-jags-and-rjags/

if (!require("pacman")) install.packages("pacman")
pacman::p_load(rjags)

BayesCJags <- function(y, Z, threshold = FALSE, niter = 5000, burnin = 2500, betaPriorAlpha = 2, betaPriorBeta = 5) {


	###############################
	### First the model strings ###
	###############################

	# probit threshold model
	modelstringThreshold = "

	model {

		for (i in 1:n) {

			y[i] ~ dbern(prob[i])
			probit(prob[i]) <- alpha + inprod(X[i,],beta)

		}

		for (j in 1:p) {

			ind[j] ~ dbern(pind)
			betaT[j] ~ dnorm(0,tau)
			beta[j] <- ind[j] * betaT[j]

		}

		tau ~ dgamma(1,0.001)
		alpha ~ dnorm(0,0.0001)
		# pind ~ dunif(0,1)
		pind ~ dbeta(betaPriorAlpha, betaPriorBeta)

	}

	"

	# and the gaussian model
	modelstringGaussian = "

	model {

		for (i in 1:n) {

			y[i] ~ dnorm(alpha + inprod(X[i,],beta), tauE)

		}

		for (j in 1:p) {

			ind[j] ~ dbern(pind)
			betaT[j] ~ dnorm(0,tauB)
			beta[j] <- ind[j] * betaT[j]

		}

		tauB ~ dgamma(1,0.001)
		alpha ~ dnorm(0,0.0001)
		# pind ~ dunif(0,1)
		pind ~ dbeta(betaPriorAlpha, betaPriorBeta)
		tauE ~ dgamma(1, 0.001)

	}

	"

	# the default (gaussian)
	modelstring <- modelstringGaussian

	if(threshold) {

		modelstring = modelstringThreshold

	}


	# some checks
	if(!(is.numeric(y) | is.numeric(y) | is.vector(y))) stop("y must be a numeric/integer vector")
	if(threshold & length(unique(y)) > 2) stop("y must be 0/1 for the threshold model")

	if(!is.matrix(Z)) stop("Z must be a numeric matrix")
	if(nrow(Z) != length(y)) stop("Z must have as many rows as elements in y")

	# make the data.frame
	n = length(y)
	p <- ncol(Z)

	data = list(y = y, X=Z ,n = n,p = p, betaPriorAlpha = betaPriorAlpha, betaPriorBeta = betaPriorBeta)
	init = list(alpha = 0, betaT = rep(0,p), pind = 0, ind = rep(0,p), tauB = 1)
	if(!threshold) init$tauE = 1

	load.module("glm")
	model = jags.model(textConnection(modelstring),
			 data = data, inits = init)

	# burnin
	update(model, n.iter = burnin)

	# sampling
	output = coda.samples(model = model,
			    variable.names = names(init),
			    n.iter = (niter - burnin), thin = 1)

	return(output)

}

