#Correct calculation from http://christophj.github.io/replicating/r/replicating-goyal-welch-2008/
loadAnnualData <- function(dataFreq, logRetunrs) {
 
  data  <- read.csv2(paste("data/data_", dataFreq, ".csv", sep=""), 
                       na.strings="NaN", stringsAsFactors=FALSE)
  data  <- as.data.table(data)
  
  #Return
  #TOOD: they use continuously compounded returns!!
  if(logReturns) {
    returns <- log(data$Index) - log(Lag(data$Index)) #rt = log(Pt) - log(Pt-1)
  } else {
    returns <- (data$Index - Lag(data$Index)) / Lag(data$Index) #Rt = (Pt - Pt-1)/Pt-1
  }
  excReturns <- returns - data$Rfree
  
  #Rgressors 
  if(logReturns) {
    df <- data.frame(y = excReturns,
                     dp = log(data$D12) - log(data$Index), #Dividend Price Ratio is the difference between the log of dividends and the log of prices.
                     dy = log(data$D12) - log(Lag(data$Index)), #Dividend Yield (d/y) is the difference between the log of dividends and the log of lagged prices.
                     ep <- log(data$E12) - log(data$Index), #Earnings Price Ratio (e/p) is the difference between log of earnings and log of prices. 
                     de <- log(data$D12) - log(data$E12), #Dividend Payout Ratio (d/e) is the difference between log of dividends and log of earnings.
                     svar <- data$svar, #Stock Variance
                     csp <- data$csp, #Cross-Sectional yium
                     bm <- data$b.m, #Book to Market Ratio (b/m)
                     ntis <- data$ntis, #Net Equity Expansion
                     tbl <- data$tbl, #Treasury Bills
                     lty <- data$lty, #Long Term Yield
                     ltr <- data$ltr, #Long Term Rate of Return 
                     tms <- data$lty - data$tbl,  #Term Spread is the difference between the long term yield on government bonds and the T-bill.
                     dfy <- data$BAA - data$AAA, #Default Yield Spread is the difference between BAA- and AAA- rated corporate bond yields.
                     dfr <- data$corpr - data$ltr, #The Default Return Spread is the difference between the return on long-term corporate bonds and returns on the long-term government bonds.
                     infl <- data$infl, #Inflation
                     ik <- data$ik #Investment to Capital Ratio
                    )
    } else {
      df <- data.frame(y = excReturns,
                       dp = data$D12 / data$Index, #Dividend Price Ratio is the difference between the log of dividends and the log of prices.
                       dy = data$D12 / Lag(data$Index), #Dividend Yield (d/y) is the difference between the log of dividends and the log of lagged prices.
                       ep <- data$E12 / data$Index, #Earnings Price Ratio (e/p) is the difference between log of earnings and log of prices. 
                       de <- data$D12 / data$E12, #Dividend Payout Ratio (d/e) is the difference between log of dividends and log of earnings.
                       svar <- data$svar, #Stock Variance
                       csp <- data$csp, #Cross-Sectional yium
                       bm <- data$b.m, #Book to Market Ratio (b/m)
                       ntis <- data$ntis, #Net Equity Expansion
                       tbl <- data$tbl, #Treasury Bills
                       lty <- data$lty, #Long Term Yield
                       ltr <- data$ltr, #Long Term Rate of Return 
                       tms <- data$lty - data$tbl,  #Term Spread is the difference between the long term yield on government bonds and the T-bill.
                       dfy <- data$BAA - data$AAA, #Default Yield Spread is the difference between BAA- and AAA- rated corporate bond yields.
                       dfr <- data$corpr - data$ltr, #The Default Return Spread is the difference between the return on long-term corporate bonds and returns on the long-term government bonds.
                       infl <- data$infl, #Inflation
                       ik <- data$ik #Investment to Capital Ratio
      )
    }
  
  #Update column names
  colnames(df) = c("y", "dp", "dy", "ep", "de", "svar", "csp", "bm", "ntis", "tbl", "lty", "ltr", "tms", "dfy", "dfr", "infl", "ik")
 
  #ntis is only available for annual data
  if(dataFreq == "annually") { 
    ntis <- data$eqis #Percent Equity Issuing - only available for data data
    cbind(df, ntis)
    colnames(df)[ncol(df)] = "ntis"
  }
  
  #Convert to time series object and return
  if(dataFreq == "quarterly") { tsdata <- zoo(x=df, order.by=as.yearqtr(as.character(data$yyyyq), format="%Y%q"))
  } else if(dataFreq == "monthly") { tsdata <- zoo(x=df, order.by=as.yearmon(as.character(data$yyyyq), format="%Y%m"))
  } else if(dataFreq == "annually") { tsdata <- zoo(x=df, order.by=as.date(as.character(data$yyyyq), format="%Y"))
  } else { tsdata <- null }
  return (tsdata)
  
}

#Check excess return calculation is correct (comparison values from paper p. 1457 on the right)
checkEquPremCalculation <- function(tsData) {
  cat("============================================================\n")
  cat(" Comparing calculated excess return data to Goyal and Welch\n")
  cat("============================================================\n")
  cat(mean(tsData[, "y"], na.rm=TRUE) * 100, "% == 4.85 %\n")           
  cat(sd(tsData[, "y"], na.rm=TRUE) * 100, "% == 17.79 %\n")                    
  cat(mean(tsData[index(tsData) >= 1927, "y"], na.rm=TRUE) * 100, "% == 6.04 %\n")    
  cat(sd(tsData[index(tsData) >= 1927, "y"], na.rm=TRUE) * 100, "% == 19.17 %\n")    
  cat(mean(tsData[index(tsData) >= 1965, "y"], na.rm=TRUE) * 100, "% == 4.03 %\n")    
  cat(sd(tsData[index(tsData) >= 1965, "y"], na.rm=TRUE) * 100, "% == 15.70 %\n")    
}


addOosStats <- function(oosStats, predResults, model) {
  oosStatsDT = data.table(Model=model, DRMSE=oosStats$DRMSE, DRMSS=oosStats$DRMSS, PR2=oosStats$PR2, R2S=oosStats$RS)
  predResults <- rbind(predResults, oosStatsDT) #Update statistics vector with model results
  return (predResults)
}

addPredictions <- function(predictions, predVector, model) {
  predVector <- cbind(predVector, predictions)
  colnames(predVector)[ncol(predVector)] = model
  return (predVector)
}