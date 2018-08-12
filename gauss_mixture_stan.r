# Claas Heuer, 2015
#
# adopted from: https://gist.github.com/breakbee/1456db5dad0132c3f116

########################
### Gaussian Mixture ###
########################

# some random data
set.seed(100)
N <- 1000
theta <- 0.3
mu1 <- -10
mu2 <- 10
sd1 <- 1
sd2 <- 20
y <- ifelse(runif(N) < theta, rnorm(N, mu1, sd1), rnorm(N, mu2, sd2))
 
dat = list(N=length(y),y=y, k=2)

library(rstan)

hcode = 
"
data {
    int<lower=1> N;
    int<lower=1> k;
    real y[N];
}
parameters {
    simplex[k] theta;
    real mu[k];
    real<lower=0> tau[k];
}
model {
    real ps[k];
    for (i in 1:k){
        mu[i] ~ normal(0, 1.0e+2);
    }
    for(i in 1:N){
        for(j in 1:k){
            ps[j] <- log(theta[j]) + normal_log(y[i], mu[j], tau[j]);
        }
        increment_log_prob(log_sum_exp(ps));
    }

    tau ~ cauchy(0,5);

}
"


hlm = stan(model_name="Gauss mixture", model_code = hcode, data=dat , iter = 10000, burnin = 5000, verbose = FALSE, chains=1)


