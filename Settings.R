#================================================
#   Set Global Parameters
#================================================
rolling     <- TRUE                # If TRUE use rolling window, FALSE use expanding window
data        <- "quarterly"         # Which data to load: monthly, quarterly or annually
logReturns  <- TRUE                # If TRUE use log returns, FALSE use normal returns
runAENQR    <- TRUE                # If TRUE runs the adaptive elastic net quantile regression - set to FALSE to shorten execution time
runSimpleParamRestrictions <- TRUE # Does not run parameter restrictions on univariate models and KS - clutters results
estimStart <- 1947                 # Start date of estimation
predStart  <- 1965                 # Start date of predicitons
predEnd    <- 2014                 # End date of predictions


#Goyal and Welch
# estimStart <- 1927
# predStart  <- 1965 
# predEnd    <- 2005

#Greek
# estimStart <- 1947
# predStart  <- 1965 
# predEnd    <- 2010

#CT (Parameter restrictions)
# estimStart <- 1872
# predStart  <- 1927 
# predEnd    <- 2010