#' Consistency probability calculation for one MRCT
#'
#' conProb returns the consistency probability given the fraction calculation
#'
#' @param alpha The type I error.
#' @param power Power.
#' @param pi The threshold ratio in Japan's criterion I.
#' @param regionFraction The regional fraction.
#' @returns consistencyProbability The consistency probability.
#' @examples
#' conProb(0.025,0.8,0.5,0.23)
#' @importFrom stats dnorm pnorm qnorm rbinom rnorm sd uniroot
#' @export
conProb <- function(alpha, power, pi, regionFraction){
  integerand <- function(x){
    return(pnorm((1-pi)*(x+qnorm(1-alpha)+qnorm(power))/sqrt(1/regionFraction-1))*dnorm(x))
  }
  cP <- stats::integrate(integerand,lower = -qnorm(power),upper = Inf)
  cP <- cP$value/power
  cP <- list(consistencyProbability = cP)
  return(cP)
}
