#' Regional fraction calculation for one MRCT
#'
#' reFrac returns the fraction calculation given the consistency probability
#'
#' @param alpha The type I error.
#' @param power Power.
#' @param pi The threshold ratio in Japan's criterion I.
#' @param consistencyProbability The consistency probability.
#' @returns regionFraction The regional fraction.
#' @examples
#' reFrac(0.025,0.8,0.5,0.8)
#' @export
reFrac <- function(alpha, power, pi, consistencyProbability){
  f <- function(x){
    return(twoMRCTs::conProb(alpha, power, pi, regionFraction = x)$consistencyProbability-consistencyProbability)
  }
  rF <- uniroot(f,interval = c(0,1))$root
  rF <- list(regionFraction = rF)
  return(rF)
}
