#R code for simulations for paper 'Reference based multiple imputation - what is the right variance and how to estimate it'
#Jonathan Bartlett

library(dejaVu)
library(bootImpute)

#perform MI under the copy reference assumption or J2R
impM <- function(obsData,M,proper=TRUE,impMethod="CR") {
  library(dejaVu)
  #first import obsData. Make ids all unique
  obsData$Id <- 1:dim(obsData)[1]
  obsData <- ImportSim(MakeDejaData(subset(obsData,select=c("Id", "arm")), arm="arm", Id="Id"),
                       event.times=expandEventCount(count=obsData$observed.events, time=obsData$censored.time),
                       status="dropout",
                       study.time=365, censored.time=obsData$censored.time, allow.beyond.study=FALSE)
  #fit model to observed data
  fit <- Simfit(obsData)
  if (impMethod=="CR") {
    #impute using copy reference
    imps <- Impute(fit, copy_reference(proper=proper), M)
  } else {
    #impute using J2R
    imps <- Impute(fit, weighted_j2r(trt.weight=0,proper=proper), M)
  }
  imps_list <- vector(mode = "list", length = M)
  for (i in 1:M) {
    imps_list[[i]] <- GetImputedDataSet(imps,index=i)
  }
  imps_list
}

analyseImp <- function(singleImpData) {
  library(dejaVu)
  fitted <- Simfit(singleImpData,family="negbin")
  crMIEst <- log(summary(fitted)$treatment.effect)
  crMIEst
}
  
#define function to run simulations
runSim <- function(nSim=100,nBoot=500,dropoutRate=0.0025,evRates=c(0.01,0.005)) {

  rubinCR <- array(0, dim=c(nSim,4))
  rubinJ2R <- array(0, dim=c(nSim,4))
  vonHippelCR <- array(0, dim=c(nSim,4))
  vonHippelJ2R <- array(0, dim=c(nSim,4))
  
  for (i in 1:nSim) {
    print(i)
    
    complete <- SimulateComplete(study.time=365,
                                 number.subjects=500,
                                 event.rates=evRates,
                                 dispersions=0.25)
    #impose some dropout
    with.MCAR.dropout <- SimulateDropout(complete,
                                         drop.mechanism=ConstantRateDrop(rate=dropoutRate,
                                                                         var=0))
    
    fit <- Simfit(with.MCAR.dropout)
    imps <- Impute(fit, copy_reference(), 10)
    impfitted <- Simfit(imps,
                     family="negbin")
    result <- summary(impfitted)
    rubinCR[i,] <- c(log(result$treatment.effect),
                   result$se,
                   log(result$treatment.effect)-qt(0.975,df=result$adjusted.df)*result$se,
                   log(result$treatment.effect)+qt(0.975,df=result$adjusted.df)*result$se)
    
    imps <- Impute(fit, weighted_j2r(trt.weight=0, proper=TRUE), 10)
    impfitted <- Simfit(imps,
                        family="negbin")
    result <- summary(impfitted)
    rubinJ2R[i,] <- c(log(result$treatment.effect),
                     result$se,
                     log(result$treatment.effect)-qt(0.975,df=result$adjusted.df)*result$se,
                     log(result$treatment.effect)+qt(0.975,df=result$adjusted.df)*result$se)
  
    #von Hippel and Bartlett
    
    #copy reference
    
      #save seed because with multiple cores bootImpute requires us to set the seed
      oldseed <- .Random.seed
      runBootImpute <- bootImpute(with.MCAR.dropout$data, impM, nBoot=nBoot, nImp=2, M=2, proper=FALSE,
                                  nCores=10, seed=4126, impMethod="CR")
      #restore seed to where it was
      .Random.seed <- oldseed
    
      runBootImputeAnalyse <- bootImputeAnalyse(runBootImpute, analyseImp,nCores=1, quiet=TRUE)
      vonHippelCR[i,] <- c(runBootImputeAnalyse$ests,
                         runBootImputeAnalyse$var^0.5,
                         runBootImputeAnalyse$ci)
      
    #jump to reference
      #save seed because with multiple cores bootImpute requires us to set the seed
      oldseed <- .Random.seed
      runBootImpute <- bootImpute(with.MCAR.dropout$data, impM, nBoot=nBoot, nImp=2, M=2, proper=FALSE,
                                  nCores=10, seed=4126, impMethod="J2R")
      #restore seed to where it was
      .Random.seed <- oldseed
      
      runBootImputeAnalyse <- bootImputeAnalyse(runBootImpute, analyseImp,nCores=1, quiet=TRUE)
      vonHippelJ2R[i,] <- c(runBootImputeAnalyse$ests,
                           runBootImputeAnalyse$var^0.5,
                           runBootImputeAnalyse$ci)
  }
  
  list(rubinCR=rubinCR,rubinJ2R=rubinJ2R, vonHippelCR=vonHippelCR, vonHippelJ2R=vonHippelJ2R)
}

#probability of not dropping out before end of fup is
exp(-365*0.00025)
exp(-365*0.0025)

set.seed(72347218)
nSim <- 5000
#alternative hypothesis
alt1 <- runSim(nSim=nSim, nBoot=200, dropoutRate=0.00025)
alt2 <- runSim(nSim=nSim, nBoot=200, dropoutRate=0.0025)

#null hypothesis
null1 <- runSim(nSim=nSim, nBoot=200, dropoutRate=0.00025, evRates=c(0.01,0.01))
null2 <- runSim(nSim=nSim, nBoot=200, dropoutRate=0.0025, evRates=c(0.01,0.01))

#construct results table for Beamer presentation
createSmallTab <- function(resSet) {
  resTab <- array(0, dim=c(4,4))
  #row.names(resTab) <- rep(c("Rubin J2R", "Bootstrap J2R","Rubin CR", "Bootstrap CR"))
  colnames(resTab) <- c("Est. log", "Emp.  SE", "Est. SE", "SE ratio")
  resTab[1,] <- c(mean(resSet$rubinJ2R[,1]),
                       sd(resSet$rubinJ2R[,1]),
                       mean(resSet$rubinJ2R[,2]),
                       mean(resSet$rubinJ2R[,2])/sd(resSet$rubinJ2R[,1]))
  resTab[2,] <- c(mean(resSet$vonHippelJ2R[,1]),
                       sd(resSet$vonHippelJ2R[,1]),
                       mean(resSet$vonHippelJ2R[,2]),
                       mean(resSet$vonHippelJ2R[,2])/sd(resSet$vonHippelJ2R[,1]))
  resTab[3,] <- c(mean(resSet$rubinCR[,1]),
                       sd(resSet$rubinCR[,1]),
                       mean(resSet$rubinCR[,2]),
                       mean(resSet$rubinCR[,2])/sd(resSet$rubinCR[,1]))
  resTab[4,] <- c(mean(resSet$vonHippelCR[,1]),
                       sd(resSet$vonHippelCR[,1]),
                       mean(resSet$vonHippelCR[,2]),
                       mean(resSet$vonHippelCR[,2])/sd(resSet$vonHippelCR[,1]))
  resTab <- cbind(c("Rubin J2R", "Bootstrap J2R",
          "Rubin CR", "Bootstrap CR"), format(round(resTab,3), nsmall=3))
  colnames(resTab)[1] <- "Method"
  resTab
}

rbind(createSmallTab(null1), createSmallTab(null2))
rbind(createSmallTab(alt1), createSmallTab(alt2))

library(xtable)
print(xtable(rbind(createSmallTab(null1), createSmallTab(null2))), include.rownames=FALSE)
print(xtable(rbind(createSmallTab(alt1), createSmallTab(alt2))),  include.rownames=FALSE)