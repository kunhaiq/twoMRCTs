#' Regional fraction calculation for one MRCT
#'
#' table2MRCTs returns the tables shown in the reference
#'
#' @param numSimu The number of simulation times.
#' @param choice Table Number.
#' @param seed Random seed
#' @returns Tables
#' @examples
#' table2MRCTs(100000,1,1)
#' @export
table2MRCTs <- function(numSimu,choice,seed) {
  #Table1
  table1 <- function(numSimu){
    ratio <- 1
    pi <- 0.5
    consistencyProbability <- 0.8
    alpha <- 0.025
    consistencyProb <- data.frame()
    for (power in c(0.8,0.9)) {
      for (d in c(0.1,0.15)) {
        for (pCtrl in c(0.5,0.6,0.7,0.8)) {
          sigmaTrt <- sqrt((pCtrl+d)*(1-(pCtrl+d)))
          sigmaCtrl <- sqrt(pCtrl*(1-pCtrl))
          ctrlN <- (sigmaTrt^2/ratio+sigmaCtrl^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
          ctrlN <- ceiling(ctrlN)
          trtN <- ceiling(ratio*ctrlN)
          N <- trtN + ctrlN
          regionFraction <- twoMRCTs::reFrac(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
          count <- matrix(0,nrow = numSimu,ncol = 2)
          regionCtrlN <- ceiling(regionFraction*ctrlN)
          regionTrtN <- ceiling(regionFraction*trtN)
          actualRegionFraction <- regionCtrlN/ctrlN
          trtOutcome <- matrix(rbinom(trtN*numSimu,1,pCtrl+d),nrow = numSimu,ncol = trtN)
          ctrlOutcome <- matrix(rbinom(ctrlN*numSimu,1,pCtrl),nrow = numSimu,ncol = trtN)
          regionTrtOutcome <- trtOutcome[,1:regionTrtN]
          regionCtrlOutcome <- ctrlOutcome[,1:regionCtrlN]
          D <- apply(trtOutcome,1,mean) - apply(ctrlOutcome,1,mean)
          varD <- apply(trtOutcome,1,sd)^2/trtN + apply(ctrlOutcome,1,sd)^2/ctrlN
          z <- D/sqrt(varD)
          count[,1] <- ifelse(z > qnorm(1-alpha),1,0)
          D1 <- apply(regionTrtOutcome,1,mean) - apply(regionCtrlOutcome,1,mean)
          count[,2] <- ifelse(D1 > pi*D & z > qnorm(1-alpha),1,0)
          row <- c(power, d, pCtrl, N, ceiling(regionFraction*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3), ceiling(actualRegionFraction*1000)/1000)
          consistencyProb <- rbind(consistencyProb, row)
        }
      }
      for (d in c(0.2)) {
        for (pCtrl in c(0.5,0.6,0.7)) {
          sigmaTrt <- sqrt((pCtrl+d)*(1-(pCtrl+d)))
          sigmaCtrl <- sqrt(pCtrl*(1-pCtrl))
          ctrlN <- (sigmaTrt^2/ratio+sigmaCtrl^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
          ctrlN <- ceiling(ctrlN)
          trtN <- ceiling(ratio*ctrlN)
          N <- trtN + ctrlN
          regionFraction <- twoMRCTs::reFrac(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
          count <- matrix(0,nrow = numSimu,ncol = 2)
          regionCtrlN <- ceiling(regionFraction*ctrlN)
          regionTrtN <- ceiling(regionFraction*trtN)
          actualRegionFraction <- regionCtrlN/ctrlN
          trtOutcome <- matrix(rbinom(trtN*numSimu,1,pCtrl+d),nrow = numSimu,ncol = trtN)
          ctrlOutcome <- matrix(rbinom(ctrlN*numSimu,1,pCtrl),nrow = numSimu,ncol = trtN)
          regionTrtOutcome <- trtOutcome[,1:regionTrtN]
          regionCtrlOutcome <- ctrlOutcome[,1:regionCtrlN]
          D <- apply(trtOutcome,1,mean) - apply(ctrlOutcome,1,mean)
          varD <- apply(trtOutcome,1,sd)^2/trtN + apply(ctrlOutcome,1,sd)^2/ctrlN
          z <- D/sqrt(varD)
          count[,1] <- ifelse(z > qnorm(1-alpha),1,0)
          D1 <- apply(regionTrtOutcome,1,mean) - apply(regionCtrlOutcome,1,mean)
          count[,2] <- ifelse(D1 > pi*D & z > qnorm(1-alpha),1,0)
          row <- c(power, d, pCtrl, N, ceiling(regionFraction*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3), ceiling(actualRegionFraction*1000)/1000)
          consistencyProb <- rbind(consistencyProb, row)
        }
      }
    }
    colnames(consistencyProb) <- c("power","d", "pCtrl", "N","f_k","empirical CP","actual f_k")
    return(consistencyProb)
  }
  #Table2
  table2 <- function(numSimu){
    sigmaTrt <- 4
    sigmaCtrl <- sigmaTrt
    ratio <- 1
    pi <- 0.5
    consistencyProbability <- 0.8
    alpha <- 0.025
    consistencyProb <- data.frame()
    for (power in c(0.8,0.9)) {
      for (d in c(1,1.25,1.5,2)) {
        ctrlN <- (sigmaTrt^2/ratio+sigmaCtrl^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
        ctrlN <- ceiling(ctrlN)
        trtN <- ceiling(ratio*ctrlN)
        N <- trtN + ctrlN
        regionFraction <- twoMRCTs::reFrac(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
        count <- matrix(0,nrow = numSimu,ncol = 2)
        regionCtrlN <- ceiling(regionFraction*ctrlN)
        regionTrtN <- ceiling(regionFraction*trtN)
        actualRegionFraction <- regionCtrlN/ctrlN
        trtOutcome <- matrix(rnorm(trtN*numSimu,mean = d,sd = sigmaTrt),nrow = numSimu,ncol = trtN)
        ctrlOutcome <- matrix(rnorm(ctrlN*numSimu,mean = 0,sd = sigmaCtrl),nrow = numSimu,ncol = ctrlN)
        regionTrtOutcome <- trtOutcome[,1:regionTrtN]
        regionCtrlOutcome <- ctrlOutcome[,1:regionCtrlN]
        D <- apply(trtOutcome,1,mean) - apply(ctrlOutcome,1,mean)
        varD <- apply(trtOutcome,1,sd)^2/trtN + apply(ctrlOutcome,1,sd)^2/ctrlN
        z <- D/sqrt(varD)
        count[,1] <- ifelse(z > qnorm(1-alpha),1,0)
        D1 <- apply(regionTrtOutcome,1,mean) - apply(regionCtrlOutcome,1,mean)
        count[,2] <- ifelse(D1 > pi*D & z > qnorm(1-alpha),1,0)
        row <- c(power, d, N, ceiling(regionFraction*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3), ceiling(actualRegionFraction*1000)/1000)
        consistencyProb <- rbind(consistencyProb, row)
      }
    }
    colnames(consistencyProb) <- c("power","d", "N","f_k","empirical CP","actual f_k")
    return(consistencyProb)
  }
  #Table3
  table3 <- function(numSimu){
    ratio <- 1
    pi <- 0.5
    consistencyProbability <- 0.8
    alpha <- 0.025
    consistencyProb <- data.frame()
    for (power in c(0.8,0.9)) {
      for (d in c(0.1,0.15)) {
        for (pCtrl1 in c(0.5)) {
          for (pCtrl2 in c(0.5,0.8)) {
            sigmaTrt1 <- sqrt((pCtrl1+d)*(1-(pCtrl1+d)))
            sigmaCtrl1 <- sqrt(pCtrl1*(1-pCtrl1))
            sigmaTrt2 <- sqrt((pCtrl2+d)*(1-(pCtrl2+d)))
            sigmaCtrl2 <- sqrt(pCtrl2*(1-pCtrl2))
            ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN1 <- ceiling(ctrlN1)
            trtN1 <- ceiling(ratio*ctrlN1)
            N1 <- trtN1 + ctrlN1
            ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN2 <- ceiling(ctrlN2)
            trtN2 <- ceiling(ratio*ctrlN2)
            N2 <- trtN2 + ctrlN2
            reFracTable3 <- function(alpha, power, pi, consistencyProbability){
              f <- function(x){
                return(twoMRCTs::conProb2MRCTs(alpha, power, pi, regionFraction1 = x, regionFraction2 = x, d1 = d, d2 = d, sigmaTrt1, sigmaCtrl1, sigmaTrt2, sigmaCtrl2, ratio)$consistencyProbability-consistencyProbability)
              }
              rF <- uniroot(f,interval = c(0,1))$root
              rF <- list(regionFraction = rF)
              return(rF)
            }
            regionFraction <- reFracTable3(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
            count <- matrix(0,nrow = numSimu,ncol = 2)
            regionCtrlN1 <- ceiling(regionFraction*ctrlN1)
            regionTrtN1 <- ceiling(regionFraction*trtN1)
            regionCtrlN2 <- ceiling(regionFraction*ctrlN2)
            regionTrtN2 <- ceiling(regionFraction*trtN2)
            w1 <- N1/(N1+N2)
            w2 <- 1 - w1
            actualRegionFraction1 <- regionCtrlN1/ctrlN1
            actualRegionFraction2 <- regionCtrlN2/ctrlN2
            trtOutcome1 <- matrix(rbinom(trtN1*numSimu,1,pCtrl1+d),nrow = numSimu,ncol = trtN1)
            ctrlOutcome1 <- matrix(rbinom(ctrlN1*numSimu,1,pCtrl1),nrow = numSimu,ncol = ctrlN1)
            trtOutcome2 <- matrix(rbinom(trtN2*numSimu,1,pCtrl2+d),nrow = numSimu,ncol = trtN2)
            ctrlOutcome2 <- matrix(rbinom(ctrlN2*numSimu,1,pCtrl2),nrow = numSimu,ncol = ctrlN2)
            regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
            regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
            regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
            regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
            D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
            D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
            varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
            varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
            z1 <- D1/sqrt(varD1)
            z2 <- D2/sqrt(varD2)
            count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
            D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
            D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
            Dpool <- w1*D1 + w2*D2
            Dpool1 <- w1*D11 + w2*D21
            count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
            row <- c(power,d,pCtrl1, pCtrl2,N1,N2,ceiling(regionFraction*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
            # consistencyProb equals the sum of the 1st col divided by the 2nd col, the last entry above
            consistencyProb <- rbind(consistencyProb, row)
          }
        }
        for (pCtrl1 in c(0.8)) {
          for (pCtrl2 in c(0.8)) {
            sigmaTrt1 <- sqrt((pCtrl1+d)*(1-(pCtrl1+d)))
            sigmaCtrl1 <- sqrt(pCtrl1*(1-pCtrl1))
            sigmaTrt2 <- sqrt((pCtrl2+d)*(1-(pCtrl2+d)))
            sigmaCtrl2 <- sqrt(pCtrl2*(1-pCtrl2))
            ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN1 <- ceiling(ctrlN1)
            trtN1 <- ceiling(ratio*ctrlN1)
            N1 <- trtN1 + ctrlN1
            ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN2 <- ceiling(ctrlN2)
            trtN2 <- ceiling(ratio*ctrlN2)
            N2 <- trtN2 + ctrlN2
            reFracTable3 <- function(alpha, power, pi, consistencyProbability){
              f <- function(x){
                return(twoMRCTs::conProb2MRCTs(alpha, power, pi, regionFraction1 = x, regionFraction2 = x, d1 = d, d2 = d, sigmaTrt1, sigmaCtrl1, sigmaTrt2, sigmaCtrl2, ratio)$consistencyProbability-consistencyProbability)
              }
              rF <- uniroot(f,interval = c(0,1))$root
              rF <- list(regionFraction = rF)
              return(rF)
            }
            regionFraction <- reFracTable3(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
            count <- matrix(0,nrow = numSimu,ncol = 2)
            regionCtrlN1 <- ceiling(regionFraction*ctrlN1)
            regionTrtN1 <- ceiling(regionFraction*trtN1)
            regionCtrlN2 <- ceiling(regionFraction*ctrlN2)
            regionTrtN2 <- ceiling(regionFraction*trtN2)
            w1 <- N1/(N1+N2)
            w2 <- 1 - w1
            actualRegionFraction1 <- regionCtrlN1/ctrlN1
            actualRegionFraction2 <- regionCtrlN2/ctrlN2
            trtOutcome1 <- matrix(rbinom(trtN1*numSimu,1,pCtrl1+d),nrow = numSimu,ncol = trtN1)
            ctrlOutcome1 <- matrix(rbinom(ctrlN1*numSimu,1,pCtrl1),nrow = numSimu,ncol = ctrlN1)
            trtOutcome2 <- matrix(rbinom(trtN2*numSimu,1,pCtrl2+d),nrow = numSimu,ncol = trtN2)
            ctrlOutcome2 <- matrix(rbinom(ctrlN2*numSimu,1,pCtrl2),nrow = numSimu,ncol = ctrlN2)
            regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
            regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
            regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
            regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
            D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
            D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
            varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
            varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
            z1 <- D1/sqrt(varD1)
            z2 <- D2/sqrt(varD2)
            count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
            D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
            D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
            Dpool <- w1*D1 + w2*D2
            Dpool1 <- w1*D11 + w2*D21
            count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
            row <- c(power,d,pCtrl1, pCtrl2,N1,N2,ceiling(regionFraction*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
            # consistencyProb equals the sum of the 1st col divided by the 2nd col, the last entry above
            consistencyProb <- rbind(consistencyProb, row)
          }
        }
      }
    }
    colnames(consistencyProb) <- c("power","d", "pCtrl1", "pCtrl2", "N1", "N2","f_k","empirical CP","actual f_k^1","actual f_k^2")
    return(consistencyProb)
  }
  #Table4
  table4 <- function(numSimu){
    sigmaTrt1 <- 4
    sigmaCtrl1 <- sigmaTrt1
    sigmaTrt2 <- 4
    sigmaCtrl2 <- sigmaTrt2
    ratio <- 1
    pi <- 0.5
    consistencyProbability <- 0.8
    alpha <- 0.025
    consistencyProb <- data.frame()
    for (power in c(0.8,0.9)) {
      for (d1 in c(1)) {
        for (d2 in c(1,2)) {
          ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d1^2)
          ctrlN1 <- ceiling(ctrlN1)
          trtN1 <- ceiling(ratio*ctrlN1)
          N1 <- trtN1 + ctrlN1
          ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d2^2)
          ctrlN2 <- ceiling(ctrlN2)
          trtN2 <- ceiling(ratio*ctrlN2)
          N2 <- trtN2 + ctrlN2
          reFracTable3 <- function(alpha, power, pi, consistencyProbability){
            f <- function(x){
              return(twoMRCTs::conProb2MRCTs(alpha, power, pi, regionFraction1 = x, regionFraction2 = x, d1, d2, sigmaTrt1, sigmaCtrl1, sigmaTrt2, sigmaCtrl2, ratio)$consistencyProbability-consistencyProbability)
            }
            rF <- uniroot(f,interval = c(0,1))$root
            rF <- list(regionFraction = rF)
            return(rF)
          }
          regionFraction <- reFracTable3(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
          count <- matrix(0,nrow = numSimu,ncol = 2)
          regionCtrlN1 <- ceiling(regionFraction*ctrlN1)
          regionTrtN1 <- ceiling(regionFraction*trtN1)
          regionCtrlN2 <- ceiling(regionFraction*ctrlN2)
          regionTrtN2 <- ceiling(regionFraction*trtN2)
          w1 <- N1/(N1+N2)
          w2 <- 1 - w1
          actualRegionFraction1 <- regionCtrlN1/ctrlN1
          actualRegionFraction2 <- regionCtrlN2/ctrlN2
          trtOutcome1 <- matrix(rnorm(trtN1*numSimu,mean = d1,sd = sigmaTrt1),nrow = numSimu,ncol = trtN1)
          ctrlOutcome1 <- matrix(rnorm(ctrlN1*numSimu,mean = 0,sd = sigmaCtrl1),nrow = numSimu,ncol = ctrlN1)
          trtOutcome2 <- matrix(rnorm(trtN2*numSimu,mean = d2,sd = sigmaTrt2),nrow = numSimu,ncol = trtN2)
          ctrlOutcome2 <- matrix(rnorm(ctrlN2*numSimu,mean = 0,sd = sigmaCtrl2),nrow = numSimu,ncol = ctrlN2)
          regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
          regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
          regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
          regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
          D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
          D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
          varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
          varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
          z1 <- D1/sqrt(varD1)
          z2 <- D2/sqrt(varD2)
          count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
          D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
          D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
          Dpool <- w1*D1 + w2*D2
          Dpool1 <- w1*D11 + w2*D21
          count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
          row <- c(power,d1,d2,N1,N2,ceiling(regionFraction*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
          # consistencyProb equals the sum of the 1st col divided by the 2nd col, the last entry above
          consistencyProb <- rbind(consistencyProb, row)
        }
      }
      for (d1 in c(1.5)) {
        for (d2 in c(1.5)) {
          ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d1^2)
          ctrlN1 <- ceiling(ctrlN1)
          trtN1 <- ceiling(ratio*ctrlN1)
          N1 <- trtN1 + ctrlN1
          ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d2^2)
          ctrlN2 <- ceiling(ctrlN2)
          trtN2 <- ceiling(ratio*ctrlN2)
          N2 <- trtN2 + ctrlN2
          reFracTable3 <- function(alpha, power, pi, consistencyProbability){
            f <- function(x){
              return(twoMRCTs::conProb2MRCTs(alpha, power, pi, regionFraction1 = x, regionFraction2 = x, d1, d2, sigmaTrt1, sigmaCtrl1, sigmaTrt2, sigmaCtrl2, ratio)$consistencyProbability-consistencyProbability)
            }
            rF <- uniroot(f,interval = c(0,1))$root
            rF <- list(regionFraction = rF)
            return(rF)
          }
          regionFraction <- reFracTable3(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
          count <- matrix(0,nrow = numSimu,ncol = 2)
          regionCtrlN1 <- ceiling(regionFraction*ctrlN1)
          regionTrtN1 <- ceiling(regionFraction*trtN1)
          regionCtrlN2 <- ceiling(regionFraction*ctrlN2)
          regionTrtN2 <- ceiling(regionFraction*trtN2)
          w1 <- N1/(N1+N2)
          w2 <- 1 - w1
          actualRegionFraction1 <- regionCtrlN1/ctrlN1
          actualRegionFraction2 <- regionCtrlN2/ctrlN2
          trtOutcome1 <- matrix(rnorm(trtN1*numSimu,mean = d1,sd = sigmaTrt1),nrow = numSimu,ncol = trtN1)
          ctrlOutcome1 <- matrix(rnorm(ctrlN1*numSimu,mean = 0,sd = sigmaCtrl1),nrow = numSimu,ncol = ctrlN1)
          trtOutcome2 <- matrix(rnorm(trtN2*numSimu,mean = d2,sd = sigmaTrt2),nrow = numSimu,ncol = trtN2)
          ctrlOutcome2 <- matrix(rnorm(ctrlN2*numSimu,mean = 0,sd = sigmaCtrl2),nrow = numSimu,ncol = ctrlN2)
          regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
          regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
          regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
          regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
          D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
          D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
          varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
          varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
          z1 <- D1/sqrt(varD1)
          z2 <- D2/sqrt(varD2)
          count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
          D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
          D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
          Dpool <- w1*D1 + w2*D2
          Dpool1 <- w1*D11 + w2*D21
          count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
          row <- c(power,d1,d2,N1,N2,ceiling(regionFraction*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
          # consistencyProb equals the sum of the 1st col divided by the 2nd col, the last entry above
          consistencyProb <- rbind(consistencyProb, row)
        }
      }
      for (d1 in c(2)) {
        for (d2 in c(2)) {
          ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d1^2)
          ctrlN1 <- ceiling(ctrlN1)
          trtN1 <- ceiling(ratio*ctrlN1)
          N1 <- trtN1 + ctrlN1
          ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d2^2)
          ctrlN2 <- ceiling(ctrlN2)
          trtN2 <- ceiling(ratio*ctrlN2)
          N2 <- trtN2 + ctrlN2
          reFracTable3 <- function(alpha, power, pi, consistencyProbability){
            f <- function(x){
              return(twoMRCTs::conProb2MRCTs(alpha, power, pi, regionFraction1 = x, regionFraction2 = x, d1, d2, sigmaTrt1, sigmaCtrl1, sigmaTrt2, sigmaCtrl2, ratio)$consistencyProbability-consistencyProbability)
            }
            rF <- uniroot(f,interval = c(0,1))$root
            rF <- list(regionFraction = rF)
            return(rF)
          }
          regionFraction <- reFracTable3(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
          count <- matrix(0,nrow = numSimu,ncol = 2)
          regionCtrlN1 <- ceiling(regionFraction*ctrlN1)
          regionTrtN1 <- ceiling(regionFraction*trtN1)
          regionCtrlN2 <- ceiling(regionFraction*ctrlN2)
          regionTrtN2 <- ceiling(regionFraction*trtN2)
          w1 <- N1/(N1+N2)
          w2 <- 1 - w1
          actualRegionFraction1 <- regionCtrlN1/ctrlN1
          actualRegionFraction2 <- regionCtrlN2/ctrlN2
          trtOutcome1 <- matrix(rnorm(trtN1*numSimu,mean = d1,sd = sigmaTrt1),nrow = numSimu,ncol = trtN1)
          ctrlOutcome1 <- matrix(rnorm(ctrlN1*numSimu,mean = 0,sd = sigmaCtrl1),nrow = numSimu,ncol = ctrlN1)
          trtOutcome2 <- matrix(rnorm(trtN2*numSimu,mean = d2,sd = sigmaTrt2),nrow = numSimu,ncol = trtN2)
          ctrlOutcome2 <- matrix(rnorm(ctrlN2*numSimu,mean = 0,sd = sigmaCtrl2),nrow = numSimu,ncol = ctrlN2)
          regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
          regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
          regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
          regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
          D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
          D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
          varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
          varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
          z1 <- D1/sqrt(varD1)
          z2 <- D2/sqrt(varD2)
          count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
          D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
          D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
          Dpool <- w1*D1 + w2*D2
          Dpool1 <- w1*D11 + w2*D21
          count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
          row <- c(power,d1,d2,N1,N2,ceiling(regionFraction*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
          # consistencyProb equals the sum of the 1st col divided by the 2nd col, the last entry above
          consistencyProb <- rbind(consistencyProb, row)
        }
      }
    }
    colnames(consistencyProb) <- c("power","d1", "d2", "N1", "N2","f_k","empirical CP","actual f_k^1","actual f_k^2")
    return(consistencyProb)
  }
  #Table5
  table5 <- function(numSimu){
    ratio <- 1
    pi <- 0.5
    consistencyProbability <- 0.8
    alpha <- 0.025
    consistencyProb <- data.frame()
    for (power in c(0.8)) {
      for (d in c(0.1,0.15)) {
        for (pCtrl1 in c(0.5)) {
          for (pCtrl2 in c(0.5)) {
            sigmaTrt1 <- sqrt((pCtrl1+d)*(1-(pCtrl1+d)))
            sigmaCtrl1 <- sqrt(pCtrl1*(1-pCtrl1))
            sigmaTrt2 <- sqrt((pCtrl2+d)*(1-(pCtrl2+d)))
            sigmaCtrl2 <- sqrt(pCtrl2*(1-pCtrl2))
            ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN1 <- ceiling(ctrlN1)
            trtN1 <- ceiling(ratio*ctrlN1)
            N1 <- trtN1 + ctrlN1
            ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN2 <- ceiling(ctrlN2)
            trtN2 <- ceiling(ratio*ctrlN2)
            N2 <- trtN2 + ctrlN2
            regionFraction <- twoMRCTs::reFrac2MRCTs(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
            for (regionFraction1 in c(0.1,0.08)) {
              regionFraction2 <- 1/(2/regionFraction-1/regionFraction1)
              count <- matrix(0,nrow = numSimu,ncol = 2)
              regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
              regionTrtN1 <- ceiling(regionFraction1*trtN1)
              regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
              regionTrtN2 <- ceiling(regionFraction2*trtN2)
              w1 <- N1/(N1+N2)
              w2 <- 1 - w1
              actualRegionFraction1 <- regionCtrlN1/ctrlN1
              actualRegionFraction2 <- regionCtrlN2/ctrlN2
              trtOutcome1 <- matrix(rbinom(trtN1*numSimu,1,pCtrl1+d),nrow = numSimu,ncol = trtN1)
              ctrlOutcome1 <- matrix(rbinom(ctrlN1*numSimu,1,pCtrl1),nrow = numSimu,ncol = ctrlN1)
              trtOutcome2 <- matrix(rbinom(trtN2*numSimu,1,pCtrl2+d),nrow = numSimu,ncol = trtN2)
              ctrlOutcome2 <- matrix(rbinom(ctrlN2*numSimu,1,pCtrl2),nrow = numSimu,ncol = ctrlN2)
              regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
              regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
              regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
              regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
              D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
              D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
              varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
              varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
              z1 <- D1/sqrt(varD1)
              z2 <- D2/sqrt(varD2)
              count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
              D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
              Dpool <- w1*D1 + w2*D2
              Dpool1 <- w1*D11 + w2*D21
              count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              row <- c(power,d,pCtrl1, pCtrl2,N1,N2,ceiling(regionFraction1*1000)/1000,ceiling(regionFraction2*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
              consistencyProb <- rbind(consistencyProb, row)
            }
          }
        }
        for (pCtrl1 in c(0.8)) {
          for (pCtrl2 in c(0.8)) {
            sigmaTrt1 <- sqrt((pCtrl1+d)*(1-(pCtrl1+d)))
            sigmaCtrl1 <- sqrt(pCtrl1*(1-pCtrl1))
            sigmaTrt2 <- sqrt((pCtrl2+d)*(1-(pCtrl2+d)))
            sigmaCtrl2 <- sqrt(pCtrl2*(1-pCtrl2))
            ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN1 <- ceiling(ctrlN1)
            trtN1 <- ceiling(ratio*ctrlN1)
            N1 <- trtN1 + ctrlN1
            ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN2 <- ceiling(ctrlN2)
            trtN2 <- ceiling(ratio*ctrlN2)
            N2 <- trtN2 + ctrlN2
            regionFraction <- twoMRCTs::reFrac2MRCTs(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
            for (regionFraction1 in c(0.1,0.08)) {
              regionFraction2 <- 1/(2/regionFraction-1/regionFraction1)
              count <- matrix(0,nrow = numSimu,ncol = 2)
              regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
              regionTrtN1 <- ceiling(regionFraction1*trtN1)
              regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
              regionTrtN2 <- ceiling(regionFraction2*trtN2)
              w1 <- N1/(N1+N2)
              w2 <- 1 - w1
              actualRegionFraction1 <- regionCtrlN1/ctrlN1
              actualRegionFraction2 <- regionCtrlN2/ctrlN2
              trtOutcome1 <- matrix(rbinom(trtN1*numSimu,1,pCtrl1+d),nrow = numSimu,ncol = trtN1)
              ctrlOutcome1 <- matrix(rbinom(ctrlN1*numSimu,1,pCtrl1),nrow = numSimu,ncol = ctrlN1)
              trtOutcome2 <- matrix(rbinom(trtN2*numSimu,1,pCtrl2+d),nrow = numSimu,ncol = trtN2)
              ctrlOutcome2 <- matrix(rbinom(ctrlN2*numSimu,1,pCtrl2),nrow = numSimu,ncol = ctrlN2)
              regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
              regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
              regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
              regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
              D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
              D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
              varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
              varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
              z1 <- D1/sqrt(varD1)
              z2 <- D2/sqrt(varD2)
              count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
              D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
              Dpool <- w1*D1 + w2*D2
              Dpool1 <- w1*D11 + w2*D21
              count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              row <- c(power,d,pCtrl1, pCtrl2,N1,N2,ceiling(regionFraction1*1000)/1000, ceiling(regionFraction2*1000)/1000,format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
              consistencyProb <- rbind(consistencyProb, row)
            }
          }
        }
      }
      for (d in c(0.2)) {
        for (pCtrl1 in c(0.5)) {
          for (pCtrl2 in c(0.5)) {
            sigmaTrt1 <- sqrt((pCtrl1+d)*(1-(pCtrl1+d)))
            sigmaCtrl1 <- sqrt(pCtrl1*(1-pCtrl1))
            sigmaTrt2 <- sqrt((pCtrl2+d)*(1-(pCtrl2+d)))
            sigmaCtrl2 <- sqrt(pCtrl2*(1-pCtrl2))
            ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN1 <- ceiling(ctrlN1)
            trtN1 <- ceiling(ratio*ctrlN1)
            N1 <- trtN1 + ctrlN1
            ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN2 <- ceiling(ctrlN2)
            trtN2 <- ceiling(ratio*ctrlN2)
            N2 <- trtN2 + ctrlN2
            regionFraction <- twoMRCTs::reFrac2MRCTs(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
            for (regionFraction1 in c(0.1,0.08)) {
              regionFraction2 <- 1/(2/regionFraction-1/regionFraction1)
              count <- matrix(0,nrow = numSimu,ncol = 2)
              regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
              regionTrtN1 <- ceiling(regionFraction1*trtN1)
              regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
              regionTrtN2 <- ceiling(regionFraction2*trtN2)
              w1 <- N1/(N1+N2)
              w2 <- 1 - w1
              actualRegionFraction1 <- regionCtrlN1/ctrlN1
              actualRegionFraction2 <- regionCtrlN2/ctrlN2
              trtOutcome1 <- matrix(rbinom(trtN1*numSimu,1,pCtrl1+d),nrow = numSimu,ncol = trtN1)
              ctrlOutcome1 <- matrix(rbinom(ctrlN1*numSimu,1,pCtrl1),nrow = numSimu,ncol = ctrlN1)
              trtOutcome2 <- matrix(rbinom(trtN2*numSimu,1,pCtrl2+d),nrow = numSimu,ncol = trtN2)
              ctrlOutcome2 <- matrix(rbinom(ctrlN2*numSimu,1,pCtrl2),nrow = numSimu,ncol = ctrlN2)
              regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
              regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
              regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
              regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
              D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
              D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
              varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
              varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
              z1 <- D1/sqrt(varD1)
              z2 <- D2/sqrt(varD2)
              count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
              D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
              Dpool <- w1*D1 + w2*D2
              Dpool1 <- w1*D11 + w2*D21
              count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              row <- c(power,d,pCtrl1, pCtrl2,N1,N2,ceiling(regionFraction1*1000)/1000,ceiling(regionFraction2*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
              consistencyProb <- rbind(consistencyProb, row)
            }
          }
        }
        for (pCtrl1 in c(0.7)) {
          for (pCtrl2 in c(0.7)) {
            sigmaTrt1 <- sqrt((pCtrl1+d)*(1-(pCtrl1+d)))
            sigmaCtrl1 <- sqrt(pCtrl1*(1-pCtrl1))
            sigmaTrt2 <- sqrt((pCtrl2+d)*(1-(pCtrl2+d)))
            sigmaCtrl2 <- sqrt(pCtrl2*(1-pCtrl2))
            ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN1 <- ceiling(ctrlN1)
            trtN1 <- ceiling(ratio*ctrlN1)
            N1 <- trtN1 + ctrlN1
            ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN2 <- ceiling(ctrlN2)
            trtN2 <- ceiling(ratio*ctrlN2)
            N2 <- trtN2 + ctrlN2
            regionFraction <- twoMRCTs::reFrac2MRCTs(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
            for (regionFraction1 in c(0.1,0.08)) {
              regionFraction2 <- 1/(2/regionFraction-1/regionFraction1)
              count <- matrix(0,nrow = numSimu,ncol = 2)
              regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
              regionTrtN1 <- ceiling(regionFraction1*trtN1)
              regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
              regionTrtN2 <- ceiling(regionFraction2*trtN2)
              w1 <- N1/(N1+N2)
              w2 <- 1 - w1
              actualRegionFraction1 <- regionCtrlN1/ctrlN1
              actualRegionFraction2 <- regionCtrlN2/ctrlN2
              trtOutcome1 <- matrix(rbinom(trtN1*numSimu,1,pCtrl1+d),nrow = numSimu,ncol = trtN1)
              ctrlOutcome1 <- matrix(rbinom(ctrlN1*numSimu,1,pCtrl1),nrow = numSimu,ncol = ctrlN1)
              trtOutcome2 <- matrix(rbinom(trtN2*numSimu,1,pCtrl2+d),nrow = numSimu,ncol = trtN2)
              ctrlOutcome2 <- matrix(rbinom(ctrlN2*numSimu,1,pCtrl2),nrow = numSimu,ncol = ctrlN2)
              regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
              regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
              regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
              regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
              D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
              D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
              varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
              varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
              z1 <- D1/sqrt(varD1)
              z2 <- D2/sqrt(varD2)
              count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
              D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
              Dpool <- w1*D1 + w2*D2
              Dpool1 <- w1*D11 + w2*D21
              count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              row <- c(power,d,pCtrl1, pCtrl2,N1,N2,ceiling(regionFraction1*1000)/1000,ceiling(regionFraction2*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
              consistencyProb <- rbind(consistencyProb, row)
            }
          }
        }
      }
    }
    for (power in c(0.9)) {
      for (d in c(0.1,0.15)) {
        for (pCtrl1 in c(0.5)) {
          for (pCtrl2 in c(0.5)) {
            sigmaTrt1 <- sqrt((pCtrl1+d)*(1-(pCtrl1+d)))
            sigmaCtrl1 <- sqrt(pCtrl1*(1-pCtrl1))
            sigmaTrt2 <- sqrt((pCtrl2+d)*(1-(pCtrl2+d)))
            sigmaCtrl2 <- sqrt(pCtrl2*(1-pCtrl2))
            ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN1 <- ceiling(ctrlN1)
            trtN1 <- ceiling(ratio*ctrlN1)
            N1 <- trtN1 + ctrlN1
            ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN2 <- ceiling(ctrlN2)
            trtN2 <- ceiling(ratio*ctrlN2)
            N2 <- trtN2 + ctrlN2
            regionFraction <- twoMRCTs::reFrac2MRCTs(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
            for (regionFraction1 in c(0.09,0.08)) {
              regionFraction2 <- 1/(2/regionFraction-1/regionFraction1)
              count <- matrix(0,nrow = numSimu,ncol = 2)
              regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
              regionTrtN1 <- ceiling(regionFraction1*trtN1)
              regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
              regionTrtN2 <- ceiling(regionFraction2*trtN2)
              w1 <- N1/(N1+N2)
              w2 <- 1 - w1
              actualRegionFraction1 <- regionCtrlN1/ctrlN1
              actualRegionFraction2 <- regionCtrlN2/ctrlN2
              trtOutcome1 <- matrix(rbinom(trtN1*numSimu,1,pCtrl1+d),nrow = numSimu,ncol = trtN1)
              ctrlOutcome1 <- matrix(rbinom(ctrlN1*numSimu,1,pCtrl1),nrow = numSimu,ncol = ctrlN1)
              trtOutcome2 <- matrix(rbinom(trtN2*numSimu,1,pCtrl2+d),nrow = numSimu,ncol = trtN2)
              ctrlOutcome2 <- matrix(rbinom(ctrlN2*numSimu,1,pCtrl2),nrow = numSimu,ncol = ctrlN2)
              regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
              regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
              regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
              regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
              D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
              D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
              varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
              varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
              z1 <- D1/sqrt(varD1)
              z2 <- D2/sqrt(varD2)
              count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
              D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
              Dpool <- w1*D1 + w2*D2
              Dpool1 <- w1*D11 + w2*D21
              count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              row <- c(power,d,pCtrl1, pCtrl2,N1,N2,ceiling(regionFraction1*1000)/1000,ceiling(regionFraction2*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
              consistencyProb <- rbind(consistencyProb, row)
            }
          }
        }
        for (pCtrl1 in c(0.8)) {
          for (pCtrl2 in c(0.8)) {
            sigmaTrt1 <- sqrt((pCtrl1+d)*(1-(pCtrl1+d)))
            sigmaCtrl1 <- sqrt(pCtrl1*(1-pCtrl1))
            sigmaTrt2 <- sqrt((pCtrl2+d)*(1-(pCtrl2+d)))
            sigmaCtrl2 <- sqrt(pCtrl2*(1-pCtrl2))
            ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN1 <- ceiling(ctrlN1)
            trtN1 <- ceiling(ratio*ctrlN1)
            N1 <- trtN1 + ctrlN1
            ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN2 <- ceiling(ctrlN2)
            trtN2 <- ceiling(ratio*ctrlN2)
            N2 <- trtN2 + ctrlN2
            regionFraction <- twoMRCTs::reFrac2MRCTs(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
            for (regionFraction1 in c(0.09,0.08)) {
              regionFraction2 <- 1/(2/regionFraction-1/regionFraction1)
              count <- matrix(0,nrow = numSimu,ncol = 2)
              regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
              regionTrtN1 <- ceiling(regionFraction1*trtN1)
              regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
              regionTrtN2 <- ceiling(regionFraction2*trtN2)
              w1 <- N1/(N1+N2)
              w2 <- 1 - w1
              actualRegionFraction1 <- regionCtrlN1/ctrlN1
              actualRegionFraction2 <- regionCtrlN2/ctrlN2
              trtOutcome1 <- matrix(rbinom(trtN1*numSimu,1,pCtrl1+d),nrow = numSimu,ncol = trtN1)
              ctrlOutcome1 <- matrix(rbinom(ctrlN1*numSimu,1,pCtrl1),nrow = numSimu,ncol = ctrlN1)
              trtOutcome2 <- matrix(rbinom(trtN2*numSimu,1,pCtrl2+d),nrow = numSimu,ncol = trtN2)
              ctrlOutcome2 <- matrix(rbinom(ctrlN2*numSimu,1,pCtrl2),nrow = numSimu,ncol = ctrlN2)
              regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
              regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
              regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
              regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
              D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
              D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
              varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
              varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
              z1 <- D1/sqrt(varD1)
              z2 <- D2/sqrt(varD2)
              count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
              D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
              Dpool <- w1*D1 + w2*D2
              Dpool1 <- w1*D11 + w2*D21
              count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              row <- c(power,d,pCtrl1, pCtrl2,N1,N2,ceiling(regionFraction1*1000)/1000, ceiling(regionFraction2*1000)/1000,format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
              consistencyProb <- rbind(consistencyProb, row)
            }
          }
        }
      }
      for (d in c(0.2)) {
        for (pCtrl1 in c(0.5)) {
          for (pCtrl2 in c(0.5)) {
            sigmaTrt1 <- sqrt((pCtrl1+d)*(1-(pCtrl1+d)))
            sigmaCtrl1 <- sqrt(pCtrl1*(1-pCtrl1))
            sigmaTrt2 <- sqrt((pCtrl2+d)*(1-(pCtrl2+d)))
            sigmaCtrl2 <- sqrt(pCtrl2*(1-pCtrl2))
            ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN1 <- ceiling(ctrlN1)
            trtN1 <- ceiling(ratio*ctrlN1)
            N1 <- trtN1 + ctrlN1
            ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN2 <- ceiling(ctrlN2)
            trtN2 <- ceiling(ratio*ctrlN2)
            N2 <- trtN2 + ctrlN2
            regionFraction <- twoMRCTs::reFrac2MRCTs(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
            for (regionFraction1 in c(0.09,0.08)) {
              regionFraction2 <- 1/(2/regionFraction-1/regionFraction1)
              count <- matrix(0,nrow = numSimu,ncol = 2)
              regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
              regionTrtN1 <- ceiling(regionFraction1*trtN1)
              regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
              regionTrtN2 <- ceiling(regionFraction2*trtN2)
              w1 <- N1/(N1+N2)
              w2 <- 1 - w1
              actualRegionFraction1 <- regionCtrlN1/ctrlN1
              actualRegionFraction2 <- regionCtrlN2/ctrlN2
              trtOutcome1 <- matrix(rbinom(trtN1*numSimu,1,pCtrl1+d),nrow = numSimu,ncol = trtN1)
              ctrlOutcome1 <- matrix(rbinom(ctrlN1*numSimu,1,pCtrl1),nrow = numSimu,ncol = ctrlN1)
              trtOutcome2 <- matrix(rbinom(trtN2*numSimu,1,pCtrl2+d),nrow = numSimu,ncol = trtN2)
              ctrlOutcome2 <- matrix(rbinom(ctrlN2*numSimu,1,pCtrl2),nrow = numSimu,ncol = ctrlN2)
              regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
              regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
              regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
              regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
              D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
              D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
              varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
              varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
              z1 <- D1/sqrt(varD1)
              z2 <- D2/sqrt(varD2)
              count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
              D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
              Dpool <- w1*D1 + w2*D2
              Dpool1 <- w1*D11 + w2*D21
              count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              row <- c(power,d,pCtrl1, pCtrl2,N1,N2,ceiling(regionFraction1*1000)/1000,ceiling(regionFraction2*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
              consistencyProb <- rbind(consistencyProb, row)
            }
          }
        }
        for (pCtrl1 in c(0.7)) {
          for (pCtrl2 in c(0.7)) {
            sigmaTrt1 <- sqrt((pCtrl1+d)*(1-(pCtrl1+d)))
            sigmaCtrl1 <- sqrt(pCtrl1*(1-pCtrl1))
            sigmaTrt2 <- sqrt((pCtrl2+d)*(1-(pCtrl2+d)))
            sigmaCtrl2 <- sqrt(pCtrl2*(1-pCtrl2))
            ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN1 <- ceiling(ctrlN1)
            trtN1 <- ceiling(ratio*ctrlN1)
            N1 <- trtN1 + ctrlN1
            ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
            ctrlN2 <- ceiling(ctrlN2)
            trtN2 <- ceiling(ratio*ctrlN2)
            N2 <- trtN2 + ctrlN2
            regionFraction <- twoMRCTs::reFrac2MRCTs(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
            for (regionFraction1 in c(0.09,0.08)) {
              regionFraction2 <- 1/(2/regionFraction-1/regionFraction1)
              count <- matrix(0,nrow = numSimu,ncol = 2)
              regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
              regionTrtN1 <- ceiling(regionFraction1*trtN1)
              regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
              regionTrtN2 <- ceiling(regionFraction2*trtN2)
              w1 <- N1/(N1+N2)
              w2 <- 1 - w1
              actualRegionFraction1 <- regionCtrlN1/ctrlN1
              actualRegionFraction2 <- regionCtrlN2/ctrlN2
              trtOutcome1 <- matrix(rbinom(trtN1*numSimu,1,pCtrl1+d),nrow = numSimu,ncol = trtN1)
              ctrlOutcome1 <- matrix(rbinom(ctrlN1*numSimu,1,pCtrl1),nrow = numSimu,ncol = ctrlN1)
              trtOutcome2 <- matrix(rbinom(trtN2*numSimu,1,pCtrl2+d),nrow = numSimu,ncol = trtN2)
              ctrlOutcome2 <- matrix(rbinom(ctrlN2*numSimu,1,pCtrl2),nrow = numSimu,ncol = ctrlN2)
              regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
              regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
              regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
              regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
              D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
              D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
              varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
              varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
              z1 <- D1/sqrt(varD1)
              z2 <- D2/sqrt(varD2)
              count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
              D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
              Dpool <- w1*D1 + w2*D2
              Dpool1 <- w1*D11 + w2*D21
              count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
              row <- c(power,d,pCtrl1, pCtrl2,N1,N2,ceiling(regionFraction1*1000)/1000,ceiling(regionFraction2*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
              consistencyProb <- rbind(consistencyProb, row)
            }
          }
        }
      }
    }
    colnames(consistencyProb) <- c("power","d", "pCtrl1", "pCtrl2", "N1", "N2","f_k^1","f_k^2","empirical CP","actual f_k^1","actual f_k^2")
    return(consistencyProb)
  }
  #Table6
  table6 <- function(numSimu){
    sigmaTrt1 <- 4
    sigmaCtrl1 <- sigmaTrt1
    sigmaTrt2 <- 4
    sigmaCtrl2 <- sigmaTrt2
    ratio <- 1
    pi <- 0.5
    consistencyProbability <- 0.8
    alpha <- 0.025
    consistencyProb <- data.frame()
    for (power in c(0.8)) {
      for (d in c(1,1.5,2)) {
        d1 <- d
        d2 <- d
        ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d1^2)
        ctrlN1 <- ceiling(ctrlN1)
        trtN1 <- ceiling(ratio*ctrlN1)
        N1 <- trtN1 + ctrlN1
        ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d2^2)
        ctrlN2 <- ceiling(ctrlN2)
        trtN2 <- ceiling(ratio*ctrlN2)
        N2 <- trtN2 + ctrlN2
        reFracTable3 <- function(alpha, power, pi, consistencyProbability){
          f <- function(x){
            return(twoMRCTs::conProb2MRCTs(alpha, power, pi, regionFraction1 = x, regionFraction2 = x, d1, d2, sigmaTrt1, sigmaCtrl1, sigmaTrt2, sigmaCtrl2, ratio)$consistencyProbability-consistencyProbability)
          }
          rF <- uniroot(f,interval = c(0,1))$root
          rF <- list(regionFraction = rF)
          return(rF)
        }
        regionFraction <- reFracTable3(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
        for (regionFraction1 in c(0.1,0.08)) {
          regionFraction2 <- 1/(2/regionFraction-1/regionFraction1)
          count <- matrix(0,nrow = numSimu,ncol = 2)
          regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
          regionTrtN1 <- ceiling(regionFraction1*trtN1)
          regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
          regionTrtN2 <- ceiling(regionFraction2*trtN2)
          w1 <- N1/(N1+N2)
          w2 <- 1 - w1
          actualRegionFraction1 <- regionCtrlN1/ctrlN1
          actualRegionFraction2 <- regionCtrlN2/ctrlN2
          trtOutcome1 <- matrix(rnorm(trtN1*numSimu,mean = d1,sd = sigmaTrt1),nrow = numSimu,ncol = trtN1)
          ctrlOutcome1 <- matrix(rnorm(ctrlN1*numSimu,mean = 0,sd = sigmaCtrl1),nrow = numSimu,ncol = ctrlN1)
          trtOutcome2 <- matrix(rnorm(trtN2*numSimu,mean = d2,sd = sigmaTrt2),nrow = numSimu,ncol = trtN2)
          ctrlOutcome2 <- matrix(rnorm(ctrlN2*numSimu,mean = 0,sd = sigmaCtrl2),nrow = numSimu,ncol = ctrlN2)
          regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
          regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
          regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
          regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
          D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
          D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
          varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
          varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
          z1 <- D1/sqrt(varD1)
          z2 <- D2/sqrt(varD2)
          count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
          D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
          D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
          Dpool <- w1*D1 + w2*D2
          Dpool1 <- w1*D11 + w2*D21
          count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
          row <- c(power,d1,d2,N1,N2,ceiling(regionFraction1*1000)/1000,ceiling(regionFraction2*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
          consistencyProb <- rbind(consistencyProb, row)
        }
      }
    }
    for (power in c(0.9)) {
      for (d in c(1,1.5,2)) {
        d1 <- d
        d2 <- d
        ctrlN1 <- (sigmaTrt1^2/ratio+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power))^2/(d1^2)
        ctrlN1 <- ceiling(ctrlN1)
        trtN1 <- ceiling(ratio*ctrlN1)
        N1 <- trtN1 + ctrlN1
        ctrlN2 <- (sigmaTrt2^2/ratio+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power))^2/(d2^2)
        ctrlN2 <- ceiling(ctrlN2)
        trtN2 <- ceiling(ratio*ctrlN2)
        N2 <- trtN2 + ctrlN2
        reFracTable3 <- function(alpha, power, pi, consistencyProbability){
          f <- function(x){
            return(twoMRCTs::conProb2MRCTs(alpha, power, pi, regionFraction1 = x, regionFraction2 = x, d1, d2, sigmaTrt1, sigmaCtrl1, sigmaTrt2, sigmaCtrl2, ratio)$consistencyProbability-consistencyProbability)
          }
          rF <- uniroot(f,interval = c(0,1))$root
          rF <- list(regionFraction = rF)
          return(rF)
        }
        regionFraction <- reFracTable3(alpha = alpha, power = power, pi = pi, consistencyProbability = consistencyProbability)$regionFraction
        for (regionFraction1 in c(0.09,0.08)) {
          regionFraction2 <- 1/(2/regionFraction-1/regionFraction1)
          count <- matrix(0,nrow = numSimu,ncol = 2)
          regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
          regionTrtN1 <- ceiling(regionFraction1*trtN1)
          regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
          regionTrtN2 <- ceiling(regionFraction2*trtN2)
          w1 <- N1/(N1+N2)
          w2 <- 1 - w1
          actualRegionFraction1 <- regionCtrlN1/ctrlN1
          actualRegionFraction2 <- regionCtrlN2/ctrlN2
          trtOutcome1 <- matrix(rnorm(trtN1*numSimu,mean = d1,sd = sigmaTrt1),nrow = numSimu,ncol = trtN1)
          ctrlOutcome1 <- matrix(rnorm(ctrlN1*numSimu,mean = 0,sd = sigmaCtrl1),nrow = numSimu,ncol = ctrlN1)
          trtOutcome2 <- matrix(rnorm(trtN2*numSimu,mean = d2,sd = sigmaTrt2),nrow = numSimu,ncol = trtN2)
          ctrlOutcome2 <- matrix(rnorm(ctrlN2*numSimu,mean = 0,sd = sigmaCtrl2),nrow = numSimu,ncol = ctrlN2)
          regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
          regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
          regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
          regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
          D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
          D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
          varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
          varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
          z1 <- D1/sqrt(varD1)
          z2 <- D2/sqrt(varD2)
          count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
          D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
          D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
          Dpool <- w1*D1 + w2*D2
          Dpool1 <- w1*D11 + w2*D21
          count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
          row <- c(power,d1,d2,N1,N2,ceiling(regionFraction1*1000)/1000,ceiling(regionFraction2*1000)/1000, format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
          consistencyProb <- rbind(consistencyProb, row)
        }
      }
    }
    colnames(consistencyProb) <- c("power","d1", "d2", "N1", "N2","f_k^1","f_k^2","empirical CP","actual f_k^1","actual f_k^2")
    return(consistencyProb)
  }
  ####
  functionList <- list(table1,table2,table3,table4,table5,table6)
  if (choice >= 1 && choice <= length(functionList)) {
    set.seed(seed)
    functionList[[choice]](numSimu)
  } else {
    stop("Invalid choice!")
  }
}
