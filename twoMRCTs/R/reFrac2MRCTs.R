#' Regional fraction calculation for two homogeneous MRCTs when the regional fractions are equal
#'
#'
#' @param alpha The type I error.
#' @param power Power.
#' @param pi The threshold ratio in Japan's criterion I.
#' @param consistencyProbability The consistency probability.
#'
#' @returns regionFraction The regional fraction.
#'
#' @examples
#' reFrac2MRCTs(0.025,0.8,0.5,0.8)
#' @export
reFrac2MRCTs <- function(alpha, power, pi, consistencyProbability){
  f <- function(x){
    return(twoMRCTs::conProb2MRCTsH(alpha, power, pi, regionFraction1 = x, regionFraction2 = x)$consistencyProbability-consistencyProbability)
  }
  rF <- uniroot(f,interval = c(0,1))$root
  rF <- list(regionFraction = rF)
  return(rF)
}



