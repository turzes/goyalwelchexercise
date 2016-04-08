#1-Step-Ahead (Expanding Window) Forecast
expWindowOSAforecast <- function(dataEstim, dataPred, rolling, estimFunc) {
 
  cat("\n")
  #Perform one step ahead expanding window forecasts
  predictions = NULL #data.table(date=as.Date("100001", format="%Y%m"), forecasts=0)
  for(i in 1:dim(dataPred)[1]) {
    cat(".")
    #reg <- lm(y ~ . , data=dataEstim, na.action=na.omit)      #Estimate Model
    reg <- estimFunc(dataEstim)
    #forecast <- forecast(reg, newdata=dataPred[i, colnames(dataPred[]) != "y"], h=1)$mean       #Predict t using t-1 (as all regressors are lagged)
    #forecast <- predict(reg, newdata=dataPred[i, colnames(dataPred) != "y"])       #Predict t using t-1 (as all regressors are lagged)
    forecast <- predict(reg, newdata=dataPred[i, ])  
    
    #Update Forecast Vector
    pred <- data.table(date=as.Date(index(dataPred[i]), format="%Ym"), forecasts=forecast)
    if(is.null(predictions)) {
      predictions <- pred
    } else {
     predictions <- rbind(predictions, pred)
    }
    
    #Update Data
    dataEstim = rbind(dataEstim, dataPred[i])                      
    
    #If rolling window, delete first observation
    if(rolling) {
      dataEstim <- dataEstim[2:nrow(dataEstim), ]
    }
  }
  #predictions <- predictions[2:nrow(predictions), ]
  
  return (predictions)
  
}

calcForecasts <- function(predictions) {
  #Calculate statistics
  predictions[,trueValues:=dataPred$y]
  pred_errors <- predictions$forecasts - dataPred$y
  predictions[,forecastErrors:=pred_errors]
  return (predictions)
}


#Plots Residual, ACF, PACF and QQNorm Plots
diagnosticPlots <- function(reg, title) {
  par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(reg$residuals, ylab="Residuals", main="Residuals") #Plot the Residuals
  abline(0, 0)
  qqnorm(reg$residuals)
  qqline(reg$residuals)
  Acf(reg$residuals, main="ACF")
  Pacf(reg$residuals, main="PACF")
  title(title, outer=TRUE)
}

#Calculate OOS Statistics
calculateOOSStatistics <-function(trueValues, forecasts, mseN_OS) {
  mse_OS <- mse(trueValues, forecasts)
  return ( calculateOOSStatisticsMSE(mse_OS, mseN_OS) )
}
calculateOOSStatisticsMSE <-function(mse_OS, mseN_OS) {
  DRMSE_OOS = round(sqrt(mseN_OS) - sqrt(mse_OS), 9)
  R2_OOS = (1 - mse_OS / mseN_OS)
  PR2_OOS = round(R2_OOS*100, 5)
  list = list(DRMSE=DRMSE_OOS, DRMSS=determineDRMSESignificance(DRMSE_OOS), PR2=PR2_OOS, RS=determineAdjR2Significance(R2_OOS) )
  return (list)
}

#Run a Univariate Regression for each Predicor
runUnivariateRegressions <- function(tsData) {

  regressors = colnames(tsData)
  regressors = regressors[regressors[] != "yyyy" & regressors[] != "y"]
  univResults <- data.table(X="NA", b="NA",  se_b="NA", T_Stat="NA", P_Val="NA", PVS="NA", PAR2="NA", ARS="NA", MSE="NA")

  for(i in 1:length(regressors)) {
    #Run Univariate regression for predictor i
    reg <- lm(y ~ eval(parse(text=regressors[i])), data=tsData, na.action=na.omit)    

    #Get relevant statistics
    results <- data.table(X=regressors[i], 
                          b=round(coefficients(reg)[2], 4), 
                          se_b=round(coef(summary(reg))[2, "Std. Error"],4), 
                          T_Stat=round(coef(summary(reg))[2, "t value"],4), 
                          P_Val=round(coef(summary(reg))[2, "Pr(>|t|)"],4),  
                          PVS=determineCoefficientSignificance(coef(summary(reg))[2, "Pr(>|t|)"]),
                          PAR2=round(summary(reg)$adj.r.squared,4)*100, 
                          ARS=determineAdjR2Significance(summary(reg)$adj.r.squared),
                          MSE=round(sum(summary(reg)$residuals^2)/nrow(model.frame(reg)),4))
    
    univResults <- rbind(univResults, results)
    
  }
    
  univResults <- univResults[2:nrow(univResults),]
  return (univResults)
}

#Run the Kitchen Sink Regression
runKitchenSink <- function(tsData) {
  regressors = colnames(tsData)
  regressors = regressors[regressors[] != "yyyy" & regressors[] != "y"]
  
  #Run a regression including all (lagged) predictors
  x = zoo(x=tsData[, regressors], order.by=index(tsData))
  tData = cbind(tsData[, "y"], x)
  colnames(tData)[1] = "y"
  reg <- lm(y ~ . , data=tData, na.action=na.omit)

  return (reg)
}

runPrincipalComponentAnalysis <- function(dataEstimTemp, dataPred, rolling) {
  pred <- NULL
  
  #========================================
  #Cross-Validate
  #========================================
  cvm <- vector()
  
  #Loop over alphas
  cat("\n... cross validating Principal Components Analysis")
  for (p in 1:15) {
    cat("\n # of PC = ", p)
    #Split estimation data into training (70%) and validation (30%) sample
    split = as.integer(dim(dataEstimTemp)[1]*0.7)
    dataTrain <- dataEstimTemp[1:split]
    dataVal <- dataEstimTemp[(split + 1):nrow(dataEstimTemp)]
    pred <- NULL
    
    #Loop over validation sample
    for(i in 1:nrow(dataVal)) {
      cat(".")
      
      #Estimate
      pc <- prcomp(dataTrain[, colnames(dataTrain) != "y"])
      reg <- lm(dataTrain[, "y"]~pc$x[,1:p])
      
      #Predict
      pc <- predict(pc, newdata=dataVal[i, colnames(dataVal) != "y"])
      y = reg$coef[1] + pc[1:p]%*%reg$coef[2:(p+1)]
      
      if(is.null(pred)) {
        pred <- y
      } else {
        pred <- rbind(pred, y)
      }
      
      dataTrain <- rbind(dataTrain, dataPred[i, ])
      if(rolling) {
        dataTrain <- dataTrain[2:nrow(dataTrain), ]
      }
    }
    
    predictions <- zoo(x=pred, order.by=(index(dataPred)))
    
    #Calculate results for partiucular # of PC
    cvm[p] <- mse(dataPred[, "y"], predictions)
  }
  
  #Select number of principal components with min MSE
  p <- which(cvm == min(cvm), arr.ind=TRUE)
  cat("\n # of Principal Components yielding Min MSE: ", p)
  
  cat("\n... forecasting using Principal Components Analysis\n")
  
  #========================================
  #Predict
  #========================================
  pred <- NULL
  for(i in 1:nrow(dataPred)) {
    cat(".")
    
    #Estimate
    pc <- prcomp(dataEstimTemp[, colnames(dataEstimTemp) != "y"])
    reg <- lm(dataEstimTemp[, "y"]~pc$x[,1:p])
    
    #Predict
    pc <- predict(pc, newdata=dataPred[i, colnames(dataPred) != "y"])
    y = reg$coef[1] + pc[1:p]%*%reg$coef[2:(p+1)]
    
    if(is.null(pred)) {
      pred <- y
    } else {
      pred <- rbind(pred, y)
    }
    
    dataEstimTemp <- rbind(dataEstimTemp, dataPred[i, ])
    if(rolling) {
      dataEstimTemp <- dataEstimTemp[2:nrow(dataEstimTemp), ]
    }
  }
  
  predictions <- zoo(x=pred, order.by=(index(dataPred)))
  
  return (predictions)
}


runQuantilePrincipalComponentAnalysis <- function(dataEstimTemp, dataPred, rolling) {
  pred <- NULL
  
  #========================================
  #Cross-Validate
  #========================================
  cvm <- vector()
  
  #Loop over alphas
  cat("\n... cross validating Principal Components Analysis")
  for (p in 1:15) {
    cat("\n # of PC = ", p)
    #Split estimation data into training (70%) and validation (30%) sample
    split = as.integer(dim(dataEstimTemp)[1]*0.7)
    dataTrain <- as.matrix(dataEstimTemp[1:split])
    dataVal <- dataEstimTemp[(split + 1):nrow(dataEstimTemp)]
    pred <- NULL
    
    #Loop over validation sample
    for(i in 1:nrow(dataVal)) {
      cat(".")
      
      #Estimate
      pc <- prcomp(dataTrain[, colnames(dataTrain) != "y"])
      reg <- rq(dataTrain[, "y"]~pc$x[,1:p], tau=c(.1, .25, .5, .75, .9))

      #weight the quantile forecasts for each predictor to get a prediction of the central tendency
      comb <- (0.05 * reg$coeff[,1] + 
                 0.25 * reg$coeff[,2]+ 
                 0.40 * reg$coeff[,3] +
                 0.25 * reg$coeff[,4] +
                 0.05 * reg$coeff[,5])
      
      #Predict
      pc <- predict(pc, newdata=dataVal[i, colnames(dataVal) != "y"])
      y = comb[1] + pc[1:p]%*%comb[2:(p+1)]
      
      if(is.null(pred)) {
        pred <- y
      } else {
        pred <- rbind(pred, y)
      }
      
      dataTrain <- rbind(dataTrain, dataPred[i, ])
      if(rolling) {
        dataTrain <- dataTrain[2:nrow(dataTrain), ]
      }
    }
    
    predictions <- zoo(x=pred, order.by=(index(dataPred)))
    
    #Calculate results for partiucular # of PC
    cvm[p] <- mse(dataPred[, "y"], predictions)
  }
  
  #Select number of principal components with min MSE
  p <- which(cvm == min(cvm), arr.ind=TRUE)
  cat("\n # of Principal Components yielding Min MSE: ", p)
  
  cat("\n... forecasting using Principal Components Analysis\n")
  
  #========================================
  #Predict
  #========================================
  pred <- NULL
  for(i in 1:nrow(dataPred)) {
    cat(".")
    
    #Estimate
    pc <- prcomp(dataEstimTemp[, colnames(dataEstimTemp) != "y"])
    reg <- rq(dataEstimTemp[, "y"]~pc$x[,1:p], tau=c(.1, .25, .5, .75, .9))
    
    #weight the quantile forecasts for each predictor to get a prediction of the central tendency
    comb <- (0.05 * reg$coeff[,1] + 
               0.25 * reg$coeff[,2]+ 
               0.40 * reg$coeff[,3] +
               0.25 * reg$coeff[,4] +
               0.05 * reg$coeff[,5])
    
    #Predict
    pc <- predict(pc, newdata=dataPred[i, colnames(dataPred) != "y"])
    y = comb[1] + pc[1:p]%*%comb[2:(p+1)]
    
    if(is.null(pred)) {
      pred <- y
    } else {
      pred <- rbind(pred, y)
    }
    
    dataEstimTemp <- rbind(dataEstimTemp, dataPred[i, ])
    if(rolling) {
      dataEstimTemp <- dataEstimTemp[2:nrow(dataEstimTemp), ]
    }
  }
  
  predictions <- zoo(x=pred, order.by=(index(dataPred)))
  
  return (predictions)
}

#Estimates, Cross-Validates and Preficsts using Adaptive Elastic Net
runAdaptiveElasticNet <- function(dataEstim, dataPred, regressors, mseN_OS, alphaSteps, rolling) {
  
  #Drop rows with NA
  dataEstimComp <- dataEstim[complete.cases(dataEstim), ]
  conRegressors = regressors[regressors[] != "ik" & regressors[] != "csp"]
  
  #------------------------------------------------
  #Cross Validation
  #------------------------------------------------
  #alpha =1 for lasso only and can blend with ridge penalty down to alpha=0 ridge only
  alphalist<-seq(0,1,by=alphaSteps)
  #lambdalist <- 10^seq(10,-2,length=100)
  lambdalist <- NULL
  alpha = alphalist

  lambda = alphalist
  mse = alphalist
  cvm = rbind(alpha, lambda)
  cvm = rbind(cvm, mse)
  cvm[2:3,] = 0
  
  #Loop over alphas
  cat("\n... cross validating Adaptive Elastic Net")
  for (a in 1:length(alphalist)) {
    cat("\n alpha = ", alphalist[a])
    #Split estimation data into training (70%) and validation (30%) sample
    split = as.integer(dim(dataEstimComp)[1]*0.7)
    dataTrain <- dataEstimComp[1:split]
    dataVal <- dataEstimComp[(split + 1):nrow(dataEstimComp)]
    
    cv <- NULL
    #Loop over validation sample
    for(i in 1:nrow(dataVal)) {
      
      #Estimate on training sample
      if(is.null(lambdalist)) {
        glmmod<-glmnet(y=dataTrain[, "y"], x=dataTrain[, conRegressors],alpha=alphalist[a],family='gaussian')
        lambdalist <- glmmod$lambda  
      } else {
        glmmod<-glmnet(y=dataTrain[, "y"], x=dataTrain[, conRegressors],alpha=alphalist[a],family='gaussian', lambda=lambdalist)
      }
      
      #Predict next observation (validation sample) and add to cv vector
      pred <- predict(glmmod, newx=as.matrix(dataVal[i, conRegressors]), s=lambdalist, method="ls")
      if (is.null(cv)) { cv <- pred
      } else { cv <- rbind(cv, pred) }
    
      #Add observation to estimation data
      dataTrain <- rbind(dataTrain, dataVal[i, ])
      
      #If rolling window, delete first observation
      if(rolling) {
        dataTrain <- dataTrain[2:nrow(dataTrain), ]
      }
    }
    
    #Calculate results for partiucular a
    mse <- array(0, dim=c(1,ncol(cv)))
    for(i in 1:ncol(cv)) {
      mse[i] <- mse(dataVal[, "y"],  cv[, i])
    }
    cvm[2, a] <- colnames(cv)[which(mse == min(mse), arr.ind=TRUE)[2]] #save lambda with min cv
    cvm[3, a] <- min(mse) #save corresponding mse
    
  }
  
  #Use alpha/lambda combination that produces min crossval mse
  minMse = min(cvm[3,])
  lambdaMin <- lambdalist[as.integer(cvm[2, which(cvm[3,] == minMse, arr.ind=TRUE)])] #save lambda with min cv
  minAlpha <- as.numeric(cvm[1, which(cvm[3,] == minMse, arr.ind=TRUE)])
  cat("\n")
  print(paste("Alpha Min", minAlpha, "Lambda Min", lambdaMin, "Min MSE", minMse))
  
  #------------------------------------------------
  #Re-Estimate using whole estimation sample
  #------------------------------------------------
  cat("... forecasting using Adaptive Elastic Net")
  glmmod<-glmnet(y=dataEstimComp[, "y"], x=as.matrix(dataEstimComp[, conRegressors],alpha=minAlpha,family='gaussian'))
  eKSResults <- coefficients(glmmod, s=lambdaMin)
  
  #------------------------------------------------
  #Expanding Window Forecast
  #------------------------------------------------
  compDataPred <- dataPred[complete.cases(dataPred), ] #TODO: this removes too many cases...
  mCompDataPred <- as.matrix(compDataPred)
  eksPred <- zoo(x=NA, order.by=index(compDataPred[1:nrow(mCompDataPred),]))
  dataEstimComp <- as.matrix(dataEstimComp)
  cat("\n")
  
  #collect coefficient vector
  coeff <- NULL
  for(i in 1:nrow(mCompDataPred)) {  
    cat(".")
    pred <- predict(glmmod, newx=as.matrix(rbind(mCompDataPred[i, conRegressors]), 1), s=lambdaMin)
    eksPred[index(compDataPred[i])] = pred[1]
    dataEstimComp <- rbind(dataEstimComp, mCompDataPred[i, ])
    #If rolling window, delete first observation
    if(rolling) {
      dataEstimComp <- dataEstimComp[2:nrow(dataEstimComp), ]
    }
    glmmod<-glmnet(y=dataEstimComp[, "y"], x=dataEstimComp[, conRegressors], alpha=minAlpha,family='gaussian')
    
    #Add coefficients to coefficient vector
    if(is.null(coeff)) {
      coeff <- coefficients(glmmod, s=lambdaMin)
    } else {
      coeff <- cbind(coeff, coefficients(glmmod, s=lambdaMin))
    }
    colnames(coeff)[i] <- as.character(index(compDataPred)[i])
  }
  
  
  #------------------------------------------------
  #Calculate Test Statistics
  #------------------------------------------------
  oosStats <- calculateOOSStatistics(dataPred$y, eksPred, mseN_OS)
  oosStatsDT = data.table(Model="Adaptive Elastic Net", DRMSE=oosStats$DRMSE, DRMSS=oosStats$DRMSS, PR2=oosStats$PR2, R2S=oosStats$RS)
  
  returnVals <- list(eKSResults, lambdaMin, minAlpha, eksPred, oosStatsDT, t(coeff))
  return (returnVals)
}


#Estimates, Cross-Validates and Preficsts using Adaptive Elastic Net 
#If lambda and alpha values are passed cross validation will not be performed
runQuantileAdaptiveElasticNet <- function(dataEstim, dataPred, regressors, mseN_OS, alphaSteps, rolling, pLambda, pAlpha) {
  
  #Drop rows with NA
  dataEstimComp <- dataEstim[complete.cases(dataEstim), ]
  conRegressors = regressors[regressors[] != "ik" & regressors[] != "csp"]
  
  #------------------------------------------------
  #Cross Validation
  #------------------------------------------------
  if(is.null(pLambda) | is.null(pAlpha)) {
    #alpha =1 for lasso only and can blend with ridge penalty down to alpha=0 ridge only
    alphalist<-seq(0,1,by=alphaSteps)
    #lambdalist <- 10^seq(10,-2,length=100)
    lambdalist <- NULL
    alpha = alphalist
    
    lambda = alphalist
    mse = alphalist
    cvm = rbind(alpha, lambda)
    cvm = rbind(cvm, mse)
    cvm[2:3,] = 0
    
    #Loop over alphas
    cat("\n... cross validating Quantile Adaptive Elastic Net")
    for (a in 1:length(alphalist)) {
      cat("\n alpha = ", alphalist[a], "\n")
      #Split estimation data into training (70%) and validation (30%) sample
      split = as.integer(dim(dataEstimComp)[1]*0.7)
      dataTrain <- dataEstimComp[1:split]
      dataVal <- dataEstimComp[(split + 1):nrow(dataEstimComp)]
      
      cv <- NULL
      #Loop over validation sample
      for(i in 1:nrow(dataVal)) {
        cat(".")
        #Estimate on training sample
        if(is.null(lambdalist)) {
          reg10 <- hqreg(X=as.matrix(dataTrain[, conRegressors]), y=dataTrain[, "y"], 
                        method = "quantile", tau=0.1, alpha=alphalist[a])
          lambdalist <- reg10$lambda
        } else {
          reg10 <- hqreg(X=as.matrix(dataTrain[, conRegressors]), y=dataTrain[, "y"], 
                        method = "quantile", tau=0.1, alpha=alphalist[a], lambda=lambdalist)
        }
        reg25 <- hqreg(X=as.matrix(dataTrain[, conRegressors]), y=dataTrain[, "y"], 
              method = "quantile", tau=0.25, alpha=alphalist[a], lambda=lambdalist)
        reg50 <-hqreg(X=as.matrix(dataTrain[, conRegressors]), y=dataTrain[, "y"], 
              method = "quantile", tau=0.5, alpha=alphalist[a], lambda=lambdalist)
        reg75 <-hqreg(X=as.matrix(dataTrain[, conRegressors]), y=dataTrain[, "y"], 
              method = "quantile", tau=0.75, alpha=alphalist[a], lambda=lambdalist)
        reg90 <-hqreg(X=as.matrix(dataTrain[, conRegressors]), y=dataTrain[, "y"], 
              method = "quantile", tau=0.9, alpha=alphalist[a], lambda=lambdalist)
          
        #Predict next observation (validation sample) and add to cv vector
        pred10 <- predict(reg10 , X=as.matrix(dataVal[i, conRegressors]), s=lambdalist)
        pred25 <- predict(reg25, X=as.matrix(dataVal[i, conRegressors]), s=lambdalist)
        pred50 <- predict(reg50, X=as.matrix(dataVal[i, conRegressors]), s=lambdalist)
        pred75 <- predict(reg75, X=as.matrix(dataVal[i, conRegressors]), s=lambdalist)
        pred90 <- predict(reg90, X=as.matrix(dataVal[i, conRegressors]), s=lambdalist)
        
        #weight the quantile forecasts for each predictor to get a prediction of the central tendency
        pred  <- ( 0.05 * pred10  + 
                   0.25 * pred25 + 
                   0.40 * pred50 +
                   0.25 * pred75 +
                   0.05 * pred90 )
        
        if (is.null(cv)) { cv <- pred
        } else { cv <- rbind(cv, pred) }
        
        #Add observation to estimation data
        dataTrain <- rbind(dataTrain, dataVal[i, ])
        
        #If rolling window, delete first observation
        if(rolling) {
          dataTrain <- dataTrain[2:nrow(dataTrain), ]
        }
      }
      
      #Calculate results for partiucular a
      mse <- array(0, dim=c(1,ncol(cv)))
      for(i in 1:ncol(cv)) {
        mse[i] <- mse(dataVal[, "y"],  cv[, i])
      }
      cvm[2, a] <- colnames(cv)[which(mse == min(mse), arr.ind=TRUE)[2]] #save lambda with min cv
      cvm[3, a] <- min(mse) #save corresponding mse
      
    }
    
    #Use alpha/lambda combination that produces min crossval mse
    minMse    <- min(cvm[3,])
    lambdaMin <- cvm[2, which(cvm[3,] == minMse, arr.ind=TRUE)] #save lambda with min cv
    minAlpha  <- as.numeric(cvm[1, which(cvm[3,] == minMse, arr.ind=TRUE)])
    cat("\n")
    print(paste("Alpha Min", minAlpha, "Lambda Min", lambdaMin, "Min MSE", minMse))
    
  } else {
    lambdaMin <- pLambda
    minAlpha  <- pAlpha
  }
  
  #------------------------------------------------
  #Expanding Window Forecast
  #------------------------------------------------
  cat("... forecasting using Quantile Adaptive Elastic Net\n")
  eksPred <- zoo(x=NA, order.by=index(dataPred[1:nrow(dataPred),]))
  for(i in 1:nrow(dataPred)) {  
    cat(".")
    
    #Estimate
    reg10 <- hqreg(X=as.matrix(dataEstimComp[, conRegressors]), y=dataEstimComp[, "y"], 
                   method = "quantile", tau=0.1, alpha=alphalist[a], lambda=lambdalist)
    reg25 <- hqreg(X=as.matrix(dataEstimComp[, conRegressors]), y=dataEstimComp[, "y"], 
                   method = "quantile", tau=0.25, alpha=alphalist[a], lambda=lambdalist)
    reg50 <- hqreg(X=as.matrix(dataEstimComp[, conRegressors]), y=dataEstimComp[, "y"], 
                   method = "quantile", tau=0.5, alpha=alphalist[a], lambda=lambdalist)
    reg75 <- hqreg(X=as.matrix(dataEstimComp[, conRegressors]), y=dataEstimComp[, "y"], 
                   method = "quantile", tau=0.75, alpha=alphalist[a], lambda=lambdalist)
    reg90 <- hqreg(X=as.matrix(dataEstimComp[, conRegressors]), y=dataEstimComp[, "y"], 
                   method = "quantile", tau=0.9, alpha=alphalist[a], lambda=lambdalist)
    
    #Predict next observation (validation sample) and add to cv vector
    pred10 <- predict(reg10 , X=as.matrix(dataPred[i, conRegressors]), s=lambdalist)
    pred25 <- predict(reg25, X=as.matrix(dataPred[i, conRegressors]), s=lambdalist)
    pred50 <- predict(reg50, X=as.matrix(dataPred[i, conRegressors]), s=lambdalist)
    pred75 <- predict(reg75, X=as.matrix(dataPred[i, conRegressors]), s=lambdalist)
    pred90 <- predict(reg90, X=as.matrix(dataPred[i, conRegressors]), s=lambdalist)
    
    #weight the quantile forecasts for each predictor to get a prediction of the central tendency
    pred  <-  ( 0.05 * pred10  + 
                0.25 * pred25 + 
                0.40 * pred50 +
                0.25 * pred75 +
                0.05 * pred90 )
    
    #Add predictions to prediction vector
    eksPred[index(dataPred[i])] <- pred[, lambdaMin]
    
    #Add observation to estimation data
    dataEstimComp <- rbind(dataEstimComp, dataPred[i, ])
    
    #If rolling window, delete first observation
    if(rolling) {
      dataEstimComp <- dataEstimComp[2:nrow(dataEstimComp), ]
    }
    
  }
  
  #------------------------------------------------
  #Calculate Test Statistics
  #------------------------------------------------
  oosStats <- calculateOOSStatistics(dataPred$y, eksPred, mseN_OS)
  oosStatsDT = data.table(Model="QR-aen", DRMSE=oosStats$DRMSE, DRMSS=oosStats$DRMSS, PR2=oosStats$PR2, R2S=oosStats$RS)
  
  returnVals <- list(lambdaMin, minAlpha, eksPred, oosStatsDT)
  return (returnVals)
}



determineCoefficientSignificance <- function(pval) {
  significance = ""
  if(pval < 0.10) { significance = "*" }
  if(pval < 0.05) { significance = "**" }
  if(pval < 0.01) { significance = "***" }
  return (significance)
}
determineAdjR2Significance <- function(adj_r2) {
  significance = ""
  if(!is.null(adj_r2) & !is.na(adj_r2)) {
    if(adj_r2 > 0) { significance = "*" }
    if(adj_r2 > 2) { significance = "**" }
    if(adj_r2 > 5) { significance = "***" }
  }
  return (significance)
}
determineDRMSESignificance <- function(drmse) {
  significance <- determineAdjR2Significance(drmse)
  return (significance)
}

#Adds squares and lags to the regressors
getExpandedRegressorSet <- function(data, conRegressors) {
  squares = data[, conRegressors]^2
  colnames(squares) = paste(conRegressors, "^2", sep="")
  cubes = data[, conRegressors]^3
  colnames(cubes) = paste(conRegressors, "^3", sep="")
  lags = apply(data[, conRegressors], 2, Lag)
  colnames(lags) = paste("L", conRegressors, sep="")
  lags2 = apply(lags, 2, Lag)
  colnames(lags2) = paste("LL", conRegressors, sep="")
  #TOOD: create interactions
  #interactions = dataEstim[, conRegressors]^2
  retData = cbind(data, squares, cubes, lags, lags2)
  return (retData)
}




#Calculates all statistics: from http://christophj.github.io/replicating/r/replicating-goyal-welch-2008/
getStatistics <- function(ts_df, indep, dep, h=1, start=1872, end=2005, est_periods_OOS = 20) {
  #### IS ANALYSIS
  #1. Historical mean model
  avg   <- mean(window(ts_df, start, end)[, dep], na.rm=TRUE)
  IS_error_N <- (window(ts_df, start, end)[, dep] - avg)
  #2. OLS model
  reg <- dyn$lm(eval(parse(text=dep)) ~ eval(parse(text=indep)), data=window(ts_df, start, end))
  IS_error_A <- reg$residuals
  ### 
  
  ####OOS ANALYSIS
  OOS_error_N <- numeric(end - start - est_periods_OOS)
  OOS_error_A <- numeric(end - start - est_periods_OOS)
  #Only use information that is available up to the time at which the forecast is made
  j <- 0
  for (i in (start + est_periods_OOS):(end-1)) {
    j <- j + 1
    
    #Get the actual ERP that you want to predict
    actual_ERP <- as.numeric(window(ts_df, i+1, i+1)[, dep])
    
    #1. Historical mean model
    OOS_error_N[j] <- actual_ERP - mean(window(ts_df, start, i)[, dep], na.rm=TRUE)
    
    #2. OLS model
    reg_OOS <- dyn$lm(eval(parse(text=dep)) ~ eval(parse(text=indep)), 
                      data=window(ts_df, start, i))
    
    #Compute_error
    df <- data.frame(x=as.numeric(window(ts_df, i, i)[, indep]))
    names(df) <- indep
    pred_ERP   <- predict.lm(reg_OOS, newdata=df)
    OOS_error_A[j] <-  pred_ERP - actual_ERP
  }
  
  #Compute statistics 
  MSE_N <- mean(OOS_error_N^2)
  MSE_A <- mean(OOS_error_A^2)
  T <- length(!is.na(ts_df[, dep]))
  OOS_R2  <- 1 - MSE_A/MSE_N
  
  #Is the -1 enough (maybe -2 needed because of lag)?
  OOS_oR2 <- OOS_R2 - (1-OOS_R2)*(reg$df.residual)/(T - 1) 
  dRMSE <- sqrt(MSE_N) - sqrt(MSE_A)
  ##
  
  #### CREATE PLOT
  IS  <- cumsum(IS_error_N[2:length(IS_error_N)]^2)-cumsum(IS_error_A^2)
  OOS <- cumsum(OOS_error_N^2)-cumsum(OOS_error_A^2)
  df  <- data.frame(x=seq.int(from=start + 1 + est_periods_OOS, to=end), 
                    IS=IS[(1 + est_periods_OOS):length(IS)], 
                    OOS=OOS) #Because you lose one observation due to the lag
  
  #Shift IS errors vertically, so that the IS line begins 
  # at zero on the date of first OOS prediction. (see Goyal/Welch (2008, p. 1465))
  df$IS <- df$IS - df$IS[1] 
  df  <- melt(df, id.var="x") 
  plotGG <- ggplot(df) + 
    geom_line(aes(x=x, y=value,color=variable)) + 
    geom_rect(data=data.frame(),#Needed by ggplot2, otherwise not transparent
              aes(xmin=1973, xmax=1975,ymin=-0.2,ymax=0.2), 
              fill='red',
              alpha=0.1) + 
    scale_y_continuous('Cumulative SSE Difference', limits=c(-0.2, 0.2)) + 
    scale_x_continuous('Year')
  ##
  
  return(list(IS_error_N = IS_error_N,
              IS_error_A = reg$residuals,
              OOS_error_N = OOS_error_N,
              OOS_error_A = OOS_error_A,
              IS_R2 = summary(reg)$r.squared, 
              IS_aR2 = summary(reg)$adj.r.squared, 
              OOS_R2  = OOS_R2,
              OOS_oR2 = OOS_oR2,
              dRMSE = dRMSE,
              plotGG = plotGG))
}

#Calculates statistics
getDPEPStatistics <- function(tsData) {
  #Use DP STAT
  #Dividend-yield
  dp_stat <- getStatistics(tsData, "dp", "excRet", start=1872)
  dp_stat$plotGG
  
  #Earnings-price ratio
  ep_stat <- getStatistics(tsData, "ep", "excRet", start=1872)
  ep_stat$plotGG
}