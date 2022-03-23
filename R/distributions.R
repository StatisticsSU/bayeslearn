#' Simulates from the Scaled Inv-Chi2 distribution.
#'
#' @param n the number of draws
#' @param v degrees of freedom parameter in the Scaled Inv-Chi2 distribution
#' @param tau2 scale parameter in the Scaled Inv-Chi2 distribution
#' @return vector with n draws
#' @export
#' @examples
#' library(bayeslearn)
#' draws = rScaledInvChi2(n = 10000, v = 5, tau2 = 10)
#' mean(draws) # Should be close to the theoretical (v/(v-2))*\tau2 = 16.667
rScaledInvChi2 <- function(n, v, tau2){
  return((v*tau2)/rchisq(n, df=v))
}

#' Simulates from the Dirichlet distribution.
#'
#' @param n number of draws
#' @param alpha vector with hyperparameters
#' @return n-by-p matrix where each row is a draw.
#' @export
#' @examples
#' library(bayeslearn)
#' draws = rDirichlet(n = 1000, alpha = c(3,5,10))
#' colMeans(draws) # Close to theoretical c(3,5,10)/18 = 0.167 0.278 0.556
rDirichlet <- function(n, alpha){
  p <- length(alpha) # dimension of the generated vectors
  thetaDraws <- as.data.frame(matrix(NA, n, p)) # Storage.
  for (j in 1:p){
    thetaDraws[,j] <- rgamma(n,alpha[j],1)
  }
  for (i in 1:n){
    thetaDraws[i,] = thetaDraws[i,]/sum(thetaDraws[i,])
  }
  return(thetaDraws)
}
