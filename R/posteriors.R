#' Log posterior for logistic regression with normal prior.
#'
#' \deqn{\mathrm{Pr}(y=1|\mathbf{x}) = \frac{1}{1 +\exp(-\mathbf{x}^\top\boldsymbol{\beta})}}{Pr(y = 1 | x) = 1/(1 + exp(-x*beta))}
#' with multivariate normal prior
#' \deqn{\beta ~ N(\mu, \Sigma}{beta ~ N(mu,Sigma)}
#' @param beta_ vector of regression coefficients
#' @param y n-vector with binary responses
#' @param X n-by-p matrix with covariates, including a column of ones for intercept
#' @param mu prior mean vector for the beta coefficients
#' @param Sigma prior covariance matrix for the beta coefficients
#' @return log posterior density at beta.
#' @export
#' @examples
#' library(bayeslearn)
#' y = SUdatasets::titanic$survived
#' X = cbind(1, SUdatasets::titanic$age)
#' logPostLogisticReg(c(0,1), y, X, mu = c(0,0), Sigma = 10^2*diag(2))
logPostLogisticReg <- function(beta_, y, X, mu, Sigma){
  nPara = length(beta_)
  linPred = X%*%beta_
  logLik = sum( linPred*y -log(1 + exp(linPred)))
  logPrior = mvtnorm::dmvnorm(beta_, mu, Sigma, log=TRUE)
  return(logLik + logPrior)
}

#' Log posterior for probit regression with normal prior.
#'
#' \deqn{\mathrm{Pr}(y=1|\mathbf{x}) = \Phi(\mathbf{x}^\top\boldsymbol{\beta})}{Pr(y = 1 | x) = Phi(x*beta)}
#' with multivariate normal prior
#' \deqn{\beta ~ N(\mu, \Sigma}{beta ~ N(mu,Sigma)}
#' @param beta_ vector of regression coefficients
#' @param y n-vector with binary responses
#' @param X n-by-p matrix with covariates, including a column of ones for intercept
#' @param mu prior mean vector for the beta coefficients
#' @param Sigma prior covariance matrix for the beta coefficients
#' @return log posterior density at beta.
#' @export
#' @examples
#' library(bayeslearn)
#' y = SUdatasets::titanic$survived
#' X = cbind(1, SUdatasets::titanic$age)
#' logPostProbitReg(c(0,1), y, X, mu = c(0,0), Sigma = 10^2*diag(2))
logPostProbitReg <- function(beta_, y, X, mu, Sigma){
  nPara = length(beta_)
  linPred = X%*%beta_
  logLik = sum(y*pnorm(linPred, log.p = TRUE) + (1-y)*pnorm(linPred, log.p = TRUE, lower.tail = FALSE))
  logPrior = mvtnorm::dmvnorm(beta_, mu, Sigma, log=TRUE)
  return(logLik + logPrior)
}
