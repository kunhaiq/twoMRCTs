#' Regional fraction calculation for two homogeneous MRCTs when the regional fractions are equal
#'
#'
#' @param alpha The type I error.
#' @param power Power.
#' @param pi The threshold ratio in Japan's criterion I.
#' @param regionFraction1 The regional fraction for MRCT 1.
#' @param regionFraction2 The regional fraction for MRCT 2.
#'
#' @returns consistencyProbability The consistency probability.
#'
#' @examples
#' conProb2MRCTsH(0.025,0.8,0.5,0.127,0.127)
#' @export
conProb2MRCTsH <- function(alpha, power, pi, regionFraction1, regionFraction2){
  integerand <- function(x){
    return(pnorm((1-pi)*(sqrt(2)*x+2*qnorm(1-alpha)+2*qnorm(power))/sqrt(1/regionFraction1+1/regionFraction2-2))*(2*pnorm(x+sqrt(2)*qnorm(power))-1)*dnorm(x))
  }
  cP <- stats::integrate(integerand,lower = -sqrt(2)*qnorm(power),upper = Inf)
  cP <- cP$value/power^2
  cP <- list(consistencyProbability = cP)
  return(cP)
}



