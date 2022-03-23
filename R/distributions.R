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
