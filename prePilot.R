# script for trials and errors

# set the sequence of lambda
lambdaRatio = 1e-04
lambdaMax <- 1
nLambda <- 100
lambdaMin <- lambdaMax * lambdaRatio
loghi <- log(lambdaMax)
loglo <- log(lambdaMin)
logrange <- loghi - loglo
interval <- -logrange/(nLambda - 1)
seq.lambda <- exp(seq.int(from = loghi, to = loglo, by = interval))
seq.lambda

# set up the parameters for generating the time series
Sigma <- matrix(c(2, 1,
                  1, 2), ncol = 2)
matP <- matrix(c(0.4, 0.0,
                 0.0, 0.3), ncol = 2, byrow = TRUE)
vInit <- c(100, 120)

# run the experiment
start.time <- Sys.time()
iter <- 500
origin <- 5
h <- 12

# collect all results in this list
collect.all <- vector("list", iter)
for (j in 1:iter) {
  
  # a matrix to collect the result from cv
  collect.cv <- matrix(NA, ncol = 2, nrow = nLambda)
  
  # generated time series
  y <- sim.ves(model = "ANN", obs = 36, nsim = 1, nvariate = 2,
               persistence = matP, initial = vInit, bounds = "admissible",
               randomizer = "rnorm", mean = rep(0, 2), sd = Sigma)$data
  
  # grid search
  for (i in 1:nLambda) {
    
    # collect rmse for each k in cv
    collect.rmse.cv <- matrix(NA, ncol = origin, nrow = 1)
    for (k in 1:origin) {
      yTrain <- window(y, end = 20+k)
      yTest <- window(y, start = 20+k+1, end = 20+k+1+1)
      
      Etype <<- "A"
      Ttype <<- NULL
      Stype <<- NULL
      damped <<- FALSE
      penalty <<- "offdiag.only"
      lambda <<- seq.lambda[i]
      hyperparam <<- 1
      
      fit1 <- ves(yTrain, model = "ANN", loss = loss.shrVES, bounds = "admissible",
                  persistence = "dependent", h = 2, holdout = TRUE)
      
      collect.rmse.cv[k] <- mean(sqrt(sum((yTest[1,] - fit1$forecast[1,])^2)))
    }
    
    # calculate the mean and the standard dev of the cv rmse
    collect.cv[i,] <- c(mean(collect.rmse.cv), ((sd(collect.rmse.cv))^2)/sqrt(origin))
  
  } 
  
  # find the optimal lambda
  minLambda <- seq.lambda[which.min(collect.cv[,1])]
  
  # find the second best lambda
  rank2 <- order(collect.cv[,1])[2]
  oseLambda <- seq.lambda[rank2]
  
  ## applying the ves-shrinkage with the lambdas
  yTrain <- window(y, end = 20)
  yTest <- window(y, start = 21)
  
  lambda <<- minLambda
  fit.minL <- ves(yTrain, model = "ANN", loss = loss.shrVES, bounds = "admissible",
                  persistence = "dependent", h = h, holdout = TRUE)
  lambda <<- oseLambda
  fit.1seL <- ves(yTrain, model = "ANN", loss = loss.shrVES, bounds = "admissible",
                  persistence = "dependent", h = h, holdout = TRUE)
  # benchmark
  fit.bcm <- ves(yTrain, model = "ANN", loss = "likelihood", bounds = "admissible",
                  persistence = "dependent", h = h, holdout = TRUE)
  
  # calculate the error
  error.minL <- yTest[1:h,] - fit.minL$forecast
  error.1seL <- yTest[1:h,] - fit.1seL$forecast
  error.bcm <- yTest[1:h,] - fit.bcm$forecast
  
  # calculate the rmse
  collect.accuracy <- cbind(rowMeans(apply(error.minL^2, 2, cumsum)/(1:h)),
                            rowMeans(apply(error.1seL^2, 2, cumsum)/(1:h)),
                            rowMeans(apply(error.bcm^2, 2, cumsum)/(1:h)))
  # name column and row
  colnames(collect.accuracy) <- c("lambda.min", "lambda.1se", "benchmark")
  rownames(collect.accuracy) <- paste0("t+1-", 1:h)
  
  collect.all[[j]] <- list(lambda.min = minLambda,
                           lambda.1se = oseLambda,
                           result = collect.accuracy,
                           data = y,
                           lossValue.minL = fit.minL$lossValue,
                           lossValue.1se = fit.1seL$lossValue,
                           lossValue.bcm = fit.bcm$lossValue) 
}
end.time <- Sys.time()
end.time - start.time
save(collect.all, file = "collect_all.RData")

# junks
freqTableMinLambda <- NULL 
for (i in 1:nLambda) {
  freqTableMinLambda <- c(freqTableMinLambda, 
                          sum(sapply(collect.all, function(x) x$lambda.min) == seq.lambda[i]))
}

barplot(freqTableMinLambda, names.arg = seq.lambda)
plot(cbind(exp(seq.lambda), freqTableMinLambda), type = "l")

sapply(collect.all, function(x) x$lambda.1se)
median(sapply(collect.all, function(x) x$result[1,1]))
median(sapply(collect.all, function(x) x$result[1,2]))
median(sapply(collect.all, function(x) x$result[1,3]))
