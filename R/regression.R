#' Simulate from posterior distribution in linear regression
#'
#' Simulate from the posterior distribution of the Gaussian linear regression model
#' \deqn{y = \mathbf{x}^\top\boldsymbol{\beta} + \varepsilon, \varepsilon \sim N(0,\sigma^2)}{y = x'*beta + eps, eps ~ N(0, sigma2)}
#' with conjugate prior
#' \deqn{\beta | \sigma^22 ~ N(\mu_0, \sigma^2Omega_0^{-1})}{beta | sigma2 ~ N(mu_0, sigma2*inv(Omega_0))}
#' \deqn{\sigma^2 \sim \mathrm{Inv-}\chi^2(v_0,\sigma^2_0)}{sigma2 ~ Inv-Chi2(v_0,sigma2_0)}
#'
#' @param y n-by-1 vector with response data observations
#' @param X n-by-nCovs matrix with covariates, first column should be ones if you want an intercept.
#' @param mu_0 prior mean for beta
#' @param Omega_0 prior precision matrix for beta
#' @param v_0 degrees of freedom in the prior for sigma2
#' @param sigma2_0 location ("best guess") in the prior for sigma2
#' @param nIter Number of draws from the posterior
#' @return list with a vector of sigma2 draws (sigma2Sample) and a matrix with beta draws (betaSample)
#' @export
#' @examples
#' library(bayeslearn)
#' X = cbind(1,cars$speed)
#' y = cars$speed
#' mu_0 = c(0,0); Omega_0 = 0.1*diag(2);v_0 = 5; sigma2_0 = 10;
#' postRes = BayesLinReg(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter = 1000)
#' hist(postRes$sigma2Sample, breaks = 30, main = expression(paste("Posterior ", sigma^2)), xlab = "")
#' hist(postRes$betaSample[,1], breaks = 30, main = expression(paste("Posterior ", beta[0])), xlab = "")
#' hist(postRes$betaSample[,2], breaks = 30, main = expression(paste("Posterior ", beta[1])), xlab = "")



BayesLinReg <- function(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter){

  # Compute posterior hyperparameters
  n = length(y) # Number of observations
  nCovs = dim(X)[2] # Number of covariates
  XX = t(X)%*%X
  betaHat <- solve(XX,t(X)%*%y)
  Omega_n = XX + Omega_0
  mu_n = solve(Omega_n,XX%*%betaHat+Omega_0%*%mu_0)
  v_n = v_0 + n
  sigma2_n = as.numeric((v_0*sigma2_0 + ( t(y)%*%y + t(mu_0)%*%Omega_0%*%mu_0 - t(mu_n)%*%Omega_n%*%mu_n))/v_n)
  invOmega_n = solve(Omega_n)

  # The actual sampling
  sigma2Sample = rep(NA, nIter)
  betaSample = matrix(NA, nIter, nCovs)
  for (i in 1:nIter){

    # Simulate from p(sigma2 | y, X)
    sigma2 = rScaledInvChi2(n=1, v = v_n, tau2 = sigma2_n)
    sigma2Sample[i] = sigma2

    # Simulate from p(beta | sigma2, y, X)
    beta_ = mvtnorm::rmvnorm(n=1, mean = mu_n, sigma = sigma2*invOmega_n)
    betaSample[i,] = beta_

  }
  return(results = list(sigma2Sample = sigma2Sample, betaSample=betaSample))
}
