#' Consistency probability calculation for two MRCTs
#'
#' @param alpha The type I error.
#' @param power Power.
#' @param pi The threshold ratio in Japan's criterion I.
#' @param regionFraction1 The regional fraction for MRCT 1.
#' @param regionFraction2 The regional fraction for MRCT 2.
#' @param d1 The true mean of difference of response for MRCT 1
#' @param d2 The true mean of difference of response for MRCT 2
#' @param sigmaTrt1 The standard deviation of response in the treatment group for MRCT1
#' @param sigmaCtrl1 The standard deviation of response in the control group for MRCT1
#' @param sigmaTrt2 The standard deviation of response in the treatment group for MRCT2
#' @param sigmaCtrl2 The standard deviation of response in the control group for MRCT2
#' @param ratio The ratio between the treatment group and control group
#'
#' @returns consistencyProbability The consistency probability.
#'
#' @examples
#' conProb2MRCTs(0.025,0.8,0.5,0.127,0.127,1,1,4,4,4,4,1)
#' @export
conProb2MRCTs <- function(alpha, power, pi, regionFraction1, regionFraction2, d1, d2, sigmaTrt1, sigmaCtrl1, sigmaTrt2, sigmaCtrl2, ratio){
  N1 <- (ratio+1)*(sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d1^2)
  N2 <- (ratio+1)*(sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d2^2)
  sigmad1 <- sqrt((ratio+1)*(sigmaTrt1^2/ratio+sigmaCtrl1^2)/N1)
  sigmad2 <- sqrt((ratio+1)*(sigmaTrt2^2/ratio+sigmaCtrl2^2)/N2)
  w1 <- N1/(N1+N2)
  w2 <- N2/(N1+N2)
  integerand <- function(x){
    return(pnorm((1-pi)*(w1*sigmad1*x[1]+w2*sigmad2*x[2]+w1*d1+w2*d2)/sqrt((1/regionFraction1-1)*w1^2*sigmad1^2+(1/regionFraction2-1)*w2^2*sigmad2^2))*dnorm(x[1])*dnorm(x[2]))
  }
  cP <- cubature::adaptIntegrate(integerand,lower = c(-qnorm(power),-qnorm(power)),upper = c(Inf,Inf))
  cP <- cP$integral/(power)^2
  cP <- list(consistencyProbability = cP)
  return(cP)
}

