#================================================
#   Setup Environment
#================================================
# setwd("~/Users/markhsalmon/Documents/working directory/daniel")

#Load Libraries
library(zoo)
library(quantmod)
library(forecast)
library(Metrics)
library(data.table)
library(ggplot2)
library(lubridate)
library(dyn)
library(reshape2)
library(glmnet)
library(e1071)
library(quantreg)
library(hqreg)
library(matrixStats)

#Remove variables and plots and shift console up
# rm(list = ls())
# if(!is.null(dev.list()["RStudioGD"])) { dev.off(dev.list()["RStudioGD"]) }
cat("\n\n\n\n\n\n\n\n\n\n\n\n\n")

#Compile functions
source("Regression.r")
source("DataUtils.r")
source("Settings.r")

#================================================
#   Load Data
#================================================
cat("\nLoading data...")
tsData <- loadAnnualData(data, logRetunrs)
#checkEquPremCalculation(tsData)

#Plot Equity premium
plot(tsData[, "y"], ylab="Equity premium", xlab="Year")
abline(0,0)

#Lag Regressors
y <- zoo(tsData[, "y"], order.by=index(tsData))
x <- zoo(x=apply(tsData[, colnames(tsData) != "y" & colnames(tsData) != "yyyy"], 2, Lag), order.by=index(tsData))
tsData <- cbind(y, x)

#Drop na cases
tsData <- tsData[, colnames(tsData) != "csp"]

# Split into Estimation / Prediction data
dataEstim <- tsData[index(tsData) >= estimStart & index(tsData) < predStart]
dataPred <- tsData[index(tsData) >= predStart & index(tsData) <= predEnd]

#================================================
#   IN-Sample Analysis
#================================================
cat("\nPerfmorning In Sampe Analysis...")
#1. Historical Mean 
cat("\n...Historical Mean Model")
REG_N <- lm(y ~ 1, data = dataEstim)    
IS_MSE_N <- mean(residuals(REG_N)^2)

#2. Univariate Models
cat("\n...Univariate Models")
#TODO: use bootstrapped t statistics instead
univResults <- runUnivariateRegressions(dataEstim)

#3. Kitchen Sink
cat("\n...Kitchen Sink")
ksReg <- runKitchenSink(dataEstim)
ksResults = summary(ksReg)
diagnosticPlots(ksReg, "Kitchen Sink")


#================================================
#   Out-Of-Sample Analysis
#================================================
cat("\nPerforming Out of Sample Analysis...")
#TODO: Instead of using > 0, test whether DRMSE is significantly diff. to 0

#1. Historical Mean
#------------------------------------------------
cat("\n...Historical Mean Model")
REG_N_OS <- lm(y ~ 1, data = dataEstim)    
pred_errors <- dataPred$y - REG_N_OS$coeff[1]
mseN_OS <- sum(pred_errors^2)/length(pred_errors)

#Setup prediction vector collecting forecasts
predVector <- dataPred[, "y"]
predVector <- cbind(predVector, REG_N_OS$coeff[1])
colnames(predVector)[1] = "y"
colnames(predVector)[2] = "avg"
predResults = data.table(Model="NA", DRMSE=0, DRMSS="", PR2=0, R2S="")

#2. Univariate Models
#------------------------------------------------
cat("\n...Univariate Models")
regressors = colnames(tsData)[colnames(tsData) != "y"]

estimFunc <- function(dataEstim) { return ( lm(y ~ . , data=dataEstim, na.action=na.omit) ) }

#Calculate Out of Sample forecasts for each univariate model
for(i in 1:length(regressors)) {
  cat("\n...", regressors[i])
  predictions <- expWindowOSAforecast(dataEstim[, c("y", regressors[i])], dataPred[, c("y", regressors[i])], rolling, estimFunc)    #Forecast
  oosStats <- calculateOOSStatistics(dataPred$y, predictions$forecasts, mseN_OS)  #Calculate OOS Statistics
  predResults <- addOosStats(oosStats, predResults,  regressors[i])
  predVector <- addPredictions(predictions$forecast, predVector, regressors[i])
}

predResults <- predResults[2:nrow(predResults), ] #Drop first row

#3. Kitchen Sink
#------------------------------------------------
cat("\n...Kitchen Sink")
#tms and de cause singularity and csp is not available after 2003 - drop
conRegressors = regressors[regressors[] != "tms" & regressors[] != "csp" & regressors[] != "de"]

#Predict
#predictions <- expWindowOSAforecast(dataEstim[, c("y", conRegressors)], dataPred[, c("y", conRegressors)], rolling)
predictions <- expWindowOSAforecast(dataEstim[, c("y", conRegressors)], 
                                    dataPred[, c("y", conRegressors)], 
                                    rolling, estimFunc)
oosStats <- calculateOOSStatistics(dataPred$y, predictions$forecasts, mseN_OS)  #Calculate OOS Statistics
predResults <- addOosStats(oosStats, predResults, "Kitchen Sink")
predVector <- addPredictions(predictions$forecast, predVector, "ks")


#4. Parameter Restrictions (Campbell & Thompson)
#------------------------------------------------
cat("\n...Parameter Restrictions")

if(runSimpleParamRestrictions) {
  
  #beta coeffiecient > 0
  estimFunc <- function(dataEstim) { 
    reg <- lm(y ~ . , data=dataEstim, na.action=na.omit)
    reg$coefficients[reg$coefficients < 0] <- 0
    return ( reg ) 
  }
  
  #Predict Univariate
  for(i in 1:length(regressors)) {
    cat("\n...", regressors[i])
    predictions <- expWindowOSAforecast(dataEstim[, c("y", regressors[i])], dataPred[, c("y", regressors[i])], rolling, estimFunc)    #Forecast
    oosStats <- calculateOOSStatistics(dataPred$y, predictions$forecasts, mseN_OS)  #Calculate OOS Statistics
    predResults <- addOosStats(oosStats, predResults,  paste(regressors[i], "B>0", sep=""))
    predVector <- addPredictions(predictions$forecast, predVector, paste(regressors[i], "B>0", sep=""))
  }
  
  #Predict Kitchen Sink
  cat("\n...", "Kitchen Sink")
  predictions <- expWindowOSAforecast(dataEstim[, c("y", conRegressors)], 
                                      dataPred[, c("y", conRegressors)], 
                                      rolling, estimFunc)
  oosStats <- calculateOOSStatistics(dataPred$y, predictions$forecasts, mseN_OS)  #Calculate OOS Statistics
  predResults <- addOosStats(oosStats, predResults, "ksB>0")
  predVector <- addPredictions(predictions$forecast, predVector, "ksB>0")
  
  #forecast > 0
  comb <- as.matrix(predVector[, c(regressors, "ks")])
  comb <- ifelse(comb < 0, 0, comb) #replace negative forecasts by 0
  for(i in 1:ncol(comb)) {
    oosStats <- calculateOOSStatistics(dataPred$y, comb[, i], mseN_OS)  #Calculate OOS Statistics
    predResults <- addOosStats(oosStats, predResults, paste(colnames(comb)[i], "F>0", sep=""))
    predVector <- addPredictions(comb[, i], predVector, paste(colnames(comb)[i], "F>0", sep=""))
  }
  
  #forecast + beta > 0
  #TODO: this is where the bug comes in - all forecasts are the same
  comb <- as.matrix(predVector[, c(paste(regressors, "B>0", sep=""), "ksB>0")])
  comb <- ifelse(comb < 0, 0, comb) #replace negative forecasts by 0
  for(i in 1:ncol(comb)) {
    oosStats <- calculateOOSStatistics(dataPred$y, comb[, i], mseN_OS)  #Calculate OOS Statistics
    predResults <- addOosStats(oosStats, predResults, paste(colnames(comb)[i], "F>0", sep=""))
    predVector <- addPredictions(comb[, i], predVector, paste(colnames(comb)[i], "F>0", sep=""))
  }

}
#5. Forecast Combinations
#------------------------------------------------
cat("\n...Forecast Combinations")

#Mean Combination
cat("\n... Mean")
comb <- predVector[, regressors]
meanComb <- zoo(x=rowSums(comb)/ncol(comb), order.by=index(predVector))
oosStats <- calculateOOSStatistics(dataPred$y, meanComb, mseN_OS)
predResults <- addOosStats(oosStats, predResults, "Mean")

#Median Combination
cat("\n... Median")
mat <- as.matrix(comb)
medianComb <- zoo(x=rowMedians(mat), order.by=index(predVector))
oosStats <- calculateOOSStatistics(dataPred$y, medianComb, mseN_OS)
predResults <- addOosStats(oosStats, predResults, "Median")

#Trimmed Mean Combination
cat("\n... Trimmed Mean")
mat <- as.matrix(comb)
mat <- apply(mat, 1, sort)
mat <- mat[2:(nrow(mat)-1),]
trimmedComb <- colSums(mat)/nrow(mat)
oosStats <- calculateOOSStatistics(dataPred$y, trimmedComb, mseN_OS)
predResults <- addOosStats(oosStats, predResults, "Trimmed Mean")

#Add to predVector
predVector <- addPredictions(predictions$forecast, predVector, "mean")
predVector <- addPredictions(predictions$forecast, predVector, "median")
predVector <- addPredictions(predictions$forecast, predVector, "trim-mean")


#5. Principal Component Analysis
#------------------------------------------------
cat("\n...Principal Component Analysis")
dataEstimTemp <- dataEstim[2:nrow(dataEstim), ] #ik is missing in 1947, drop that row
predictions <- runPrincipalComponentAnalysis(dataEstimTemp, dataPred, rolling)
oosStats <- calculateOOSStatistics(dataPred$y, predictions, mseN_OS)  #Calculate OOS Statistics
predResults <- addOosStats(oosStats, predResults, "Principal Component Analysis")
predVector <- addPredictions(predictions, predVector, "pca")


#5. Efficient Kitchen Sink (Wang et. al.) - (=Adaptive Elastic Net)
#------------------------------------------------
cat("\n...Adaptive Elastic Net")
eKSResults <- runAdaptiveElasticNet(dataEstim[, c("y", conRegressors)], dataPred[, c("y", conRegressors)], conRegressors, mseN_OS, 0.1, rolling) 
predResults <- rbind(predResults, eKSResults[[5]]) #Update statistics vector with model results
predVector <- cbind(predVector, eKSResults[[4]])
colnames(predVector)[ncol(predVector)] = "aen"

#Create plot of how betas vary over time
plot(zoo(x=as.matrix(eKSResults[[6]]), order.by=index(dataPred)), main="AEN Coefficients")

#AEN with Squares and Lags
#-------------------------
# dataEstimTemp <- getExpandedRegressorSet(dataEstim, conRegressors)
# dataPredTemp <- getExpandedRegressorSet(dataPred, conRegressors)
# conRegressors2 = colnames(dataPredTemp[, colnames(dataPredTemp) != "y"])
# eKSResults2 <- runAdaptiveElasticNet(dataEstimTemp, dataPredTemp, conRegressors2, mseN_OS, 0.1, rolling) 


#5. Support Vector Regression
#------------------------------------------------
cat("\n...Support Vector Regression")

#Linear Kernel - Essentially a linear support vector regression
cat("\n... Linear Kernel")
estimFunc <- function(dataEstim) { return (svm(y ~ . , data=dataEstim, na.action=na.omit, kernel="linear")) }
predictions <- expWindowOSAforecast(dataEstim[, c("y", conRegressors)], 
                                    dataPred[, c("y", conRegressors)], 
                                    rolling, estimFunc)

oosStats <- calculateOOSStatistics(dataPred$y, predictions$forecasts, mseN_OS)  #Calculate OOS Statistics
predResults <- addOosStats(oosStats, predResults, "Support Vector Regression (Linaer)")
predVector <- addPredictions(predictions$forecast, predVector, "svrl")

#Polynomial Kernel - Extends linear SVR by polynomials as in OLS
cat("\n... Polynomial Kernel")
estimFunc <- function(dataEstim) { return (svm(y ~ . , data=dataEstim, na.action=na.omit, kernel="polynomial")) }
predictions <- expWindowOSAforecast(dataEstim[, c("y", conRegressors)], 
                                    dataPred[, c("y", conRegressors)], 
                                    rolling, estimFunc)

oosStats <- calculateOOSStatistics(dataPred$y, predictions$forecasts, mseN_OS)  #Calculate OOS Statistics
predResults <- addOosStats(oosStats, predResults, "Support Vector Regression (Polynomial)")
predVector <- addPredictions(predictions$forecast, predVector, "svrp")

#Radial Kernel - Good to pick up non-linear relationship
cat("\n... Radial Kernel")
estimFunc <- function(dataEstim) { return (svm(y ~ . , data=dataEstim, na.action=na.omit, kernel="radial")) }
predictions <- expWindowOSAforecast(dataEstim[, c("y", conRegressors)], 
                                    dataPred[, c("y", conRegressors)], 
                                    rolling, estimFunc)

oosStats <- calculateOOSStatistics(dataPred$y, predictions$forecasts, mseN_OS)  #Calculate OOS Statistics
predResults <- addOosStats(oosStats, predResults, "Support Vector Regression (Radial)")
predVector <- addPredictions(predictions$forecast, predVector, "svrr")


#X. Quantile Regression
#------------------------------------------------
cat("\n...Quantile Regression")
regressors = colnames(tsData)
regressors = colnames(tsData)[colnames(tsData) != "y"]
conRegressors = regressors[regressors[] != "y" & regressors[] != "tms" & regressors[] != "csp"]

estimFunc <- function(dataEstim) { return (rq(y ~ ., tau=c(.1, .25, .5, .75, .9), data = dataEstim)) }

#Quantile Univariate
predQR <- list()
for (i in 1:length(conRegressors)) {
 
  cat("\n...", conRegressors[i])

  #generate quantile forecasts for each predictor
  predictions <- expWindowOSAforecast(dataEstim[, c("y", conRegressors[i])], 
                                      dataPred[, c("y", conRegressors[i])], 
                                      rolling, estimFunc)
  #save quantile predictions
  predQR[[i]] <- predictions
  names(predQR)[i] <-conRegressors[i]
}

#Combine qunatiles
predQR2 <- NULL
for(i in 1:length(predQR)) {
  
  #weight the quantile forecasts for each predictor to get a prediction of the central tendency
  pred <- (0.05 * predQR[[i]]$"forecasts.tau= 0.1" + 
             0.25 * predQR[[i]]$"forecasts.tau= 0.25"+ 
             0.40 * predQR[[i]]$"forecasts.tau= 0.5" +
             0.25 * predQR[[i]]$"forecasts.tau= 0.75" +
             0.05 * predQR[[i]]$"forecasts.tau= 0.9"
  )
  
  #collect forecasts for each predictor
  if(is.null(predQR2)) {
    predQR2 <- as.matrix(pred)
  } else {
    predQR2 <- cbind(predQR2, pred)
  }
  colnames(predQR2)[ncol(predQR2)] <- conRegressors[i]
}

#Univariate Quantile Regressions
for(i in 1:ncol(predQR2)) {
  #calculate test statistics
  oosStats <- calculateOOSStatistics(dataPred$y, predQR2[i], mseN_OS)
  predResults <- addOosStats(oosStats, predResults, paste("QR-", colnames(predQR2)[i], sep=""))
  predVector <- addPredictions(predQR2[i], predVector, paste("qr-", colnames(predQR2)[i], sep=""))
}

#Quantile MEAN combination (RFC FW3) - combine the weighted forecasts for each regressor
cat("\n... Forecast Combinations")
predictions <- rowSums(predQR2)/ncol(predQR2)
predictions <- zoo(x=predictions, order.by=index(dataPred))
oosStats <- calculateOOSStatistics(dataPred$y, predictions, mseN_OS)
predResults <- addOosStats(oosStats, predResults, "QR-Mean")
predVector <- addPredictions(predictions, predVector, "qr-mean")

#Median Combination
mat <- as.matrix(comb)
medianComb <- zoo(x=rowMedians(predQR2), order.by=index(predVector))
oosStats <- calculateOOSStatistics(dataPred$y, medianComb, mseN_OS)
predResults <- addOosStats(oosStats, predResults, "QR-Median")
predVector <- addPredictions(medianComb, predVector, "qr-median")

#Trimmed Mean Combination
mat <- as.matrix(predQR2)
mat <- apply(predQR2, 1, sort)
mat <- mat[2:(nrow(mat)-1),]
trimmedComb <- colSums(mat)/nrow(mat)
oosStats <- calculateOOSStatistics(dataPred$y, trimmedComb, mseN_OS)
predResults <- addOosStats(oosStats, predResults, "QR-TrimMean")
predVector <- addPredictions(trimmedComb, predVector, "qr-TrimMean")

#Quantile Kitchen Sink Regression
cat("\n... Kitchen Sink")
conRegressors2 <- conRegressors[conRegressors != "de"]
predictions <- expWindowOSAforecast(dataEstim[, c("y", conRegressors2)], 
                                      dataPred[, c("y", conRegressors2)], 
                                      rolling, estimFunc)

pred <- (0.05 * predictions$"forecasts.tau= 0.1" + 
           0.25 * predictions$"forecasts.tau= 0.25"+ 
           0.40 * predictions$"forecasts.tau= 0.5" +
           0.25 * predictions$"forecasts.tau= 0.75" +
           0.05 * predictions$"forecasts.tau= 0.9"
      )

oosStats <- calculateOOSStatistics(dataPred$y, pred, mseN_OS)
predResults <- addOosStats(oosStats, predResults, "QR-KS")
predVector <- addPredictions(pred, predVector, "QR-ks")


#Adaptive Elastic Net + Quantile Regression
if(runAENQR) {
  cat("\n... Quantile Regression using Adaptive Elastic Net")
  eKSResults <- runQuantileAdaptiveElasticNet(dataEstim[, c("y", conRegressors)], dataPred[, c("y", conRegressors)], conRegressors, mseN_OS, 0.1, rolling, NULL, NULL) 
  predResults <- rbind(predResults, eKSResults[[4]]) #Update statistics vector with model results
  predVector <- cbind(predVector, eKSResults[[3]])
  colnames(predVector)[ncol(predVector)] = "QR-aen"
}

#------------------------------------------------
#TODO: Quantile Regression + PCA
#------------------------------------------------
cat("\n... Principal Component Analysis")
dataEstimTemp <- dataEstim[2:nrow(dataEstim), ] #ik is missing in 1947, drop that row
predictions <- runQuantilePrincipalComponentAnalysis(dataEstimTemp, dataPred, rolling)
oosStats <- calculateOOSStatistics(dataPred$y, predictions, mseN_OS)  #Calculate OOS Statistics
predResults <- addOosStats(oosStats, predResults, "QR-pca")
predVector <- addPredictions(predictions, predVector, "QR-pca")



#X. Forecast Restrictions on Better models...
#------------------------------------------------
models <- colnames(predVector[, !(colnames(predVector) %in% paste(regressors, "B>0", sep="") 
                                  | colnames(predVector) %in% paste(regressors, "B>0F>0", sep="") 
                                  | colnames(predVector) %in% paste(regressors, "F>0", sep="") 
                                  | colnames(predVector) %in% regressors
                                  | colnames(predVector) %in% c("y", "avg", "tavg", "ks", "ksB>0"))])
comb <- as.matrix(predVector[, models])
comb <- ifelse(comb < 0, 0, comb) #replace negative forecasts by 0
for(i in 1:length(models)) {
  oosStats <- calculateOOSStatistics(dataPred$y, comb[, i], mseN_OS)  #Calculate OOS Statistics
  predResults <- addOosStats(oosStats, predResults, paste(colnames(comb)[i], "_F>0", sep=""))
  predVector <- addPredictions(comb[, i], predVector, paste( colnames(comb)[i], "_F>0", sep=""))
}


#X. Time Varying Historical Mean
#------------------------------------------------
#because the AEN shrinks all betas to 0
#this suggests using a time varying historical mean as a bnechmark instead
cat("\n...Time Varying Historical Mean Model")
estimFunc <- function(dataEstim) { return ( lm(y ~ 1 , data=dataEstim, na.action=na.omit) ) }
predictions <- expWindowOSAforecast(cbind(dataEstim[, "y"], 1), cbind(dataPred[, "y"], 1), rolling, estimFunc)    #Forecast
oosStats <- calculateOOSStatistics(dataPred$y, predictions$forecasts, mseN_OS)  #Calculate OOS Statistics
predResults <- addOosStats(oosStats, predResults,  "Time Varying Mean")
predVector <- addPredictions(predictions$forecast, predVector, "tavg")

#Calculate OOS Statistics for each model relative to time varying mean
stats = data.table(Model="NA", DRMSE=0, DRMSS="", PR2=0, R2S="")
cPredVector <- predVector[, colnames(predVector) != "tavg" & colnames(predVector) != "y"]
mse_tavg <- mse(predVector[,"y"], predVector[,"tavg"])
for(i in 1:ncol(cPredVector)) {
    mse_model <- mse(predVector[,"y"], cPredVector[,i])
    oosStats <- calculateOOSStatisticsMSE(mse_model, mse_tavg)  
    stats <- addOosStats(oosStats, stats, colnames(cPredVector)[i])
}
stats <- stats[2:nrow(stats)]
stats$Model[1] <- "Historical Mean"
stats$Model[2:nrow(stats)] <- predResults$Model[1:(length(predResults$Model)-1)]







#================================================
#   Print Summary Statistics
#================================================
#Note: In sample R2 are absolute, OOS R2 are in %

#Print Results
cat("\n\n\n\n")
cat("============================\n")
cat("   In-Sample Results:\n")
cat("============================\n")
cat("\n   Univariate Results:\n")
cat("--------------------------------\n")
print(univResults)
cat("\n\n   Kitchen Sink Results:\n")
cat("--------------------------------\n")
print(ksResults)

cat("\n")
cat("============================\n")
cat("   Out-of-Sample Results:\n")
cat("============================\n")

cat("\n   Unrestricted Models:\n")
cat("--------------------------------\n")
print(predResults[predResults$Model %in% regressors])

cat("\n Beta > 0 Restriction:\n")
cat("--------------------------------\n")
print(predResults[predResults$Model %in% paste(regressors, "B>0", sep="") | predResults$Model == "ksB>0"])

cat("\n Forecast > 0 Restriction:\n")
cat("--------------------------------\n")
print(predResults[predResults$Model %in% paste(regressors, "F>0", sep="") | predResults$Model == "ksF>0"])

cat("\n Beta & Forecast > 0 Restriction:\n")
cat("--------------------------------\n")
print(predResults[predResults$Model %in% paste(regressors, "B>0F>0", sep="") | predResults$Model == "ksB>0F>0"])

cat("\n Forecast Combinations:\n")
cat("--------------------------------\n")
print(predResults[predResults$Model == "Mean" | predResults$Model == "Median" | predResults$Model == "Trimmed Mean"])

cat("\n Sophisticated Statistical Techniques:\n")
cat("--------------------------------\n")
print(predResults[predResults$Model == "Support Vector Regression (Linaer)" 
                  | predResults$Model == "Support Vector Regression (Polynomial)" 
                  | predResults$Model == "Support Vector Regression (Radial)"
                  | predResults$Model == "Adaptive Elastic Net" 
                  | predResults$Model == "Adaptive Elastic Net" 
                  | predResults$Model == "Principal Component Analysis" ])

cat("\n Quantile Regression:\n")
cat("--------------------------------\n")
print(predResults[grepl("QR-", predResults$Model) & !grepl("F>0", predResults$Model)])

cat("\n Forecast > 0 Restriction:\n")
cat("--------------------------------\n")
print(predResults[grepl("_F>0", predResults$Model)])

cat("\n Time Varying Mean:\n")
cat("--------------------------------\n")
print(predResults[predResults$Model == "Time Varying Mean"])

#cat("\n All:\n")
#cat("--------------------------------\n")
#print(predResults, nrows=200)

cat("\n\n")
cat("======================================================================\n")
cat("   Out-of-Sample Results relative to time varying historical mean:\n")
cat("======================================================================\n")
print(stats, nrows=200)

#Produce some more plots
plot(predVector)  #Plot Actual and Predicted Values
#pairs(cor(predVector[, colnames(predVector) != "y" & colnames(predVector) != "avg"]))  #Plot Correlations between forecasts


# A = sqrt(mse(predictions$trueValues, predictions$forecasts))
# N = sqrt(sum((0.05380816-dataPred[, "y"])^2)/dim(dataPred)[1])
# plot(cumsum((0.05380816-dataPred[, "y"])^2)-cumsum(predictions$forecastErrors^2))
