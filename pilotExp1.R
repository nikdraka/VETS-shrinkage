## Running the experiments

# require if you havent installed the latest version of packages below
devtools::install_github("config-i1/legion")
devtools::install_github("config-i1/smooth")
devtools::install_github("config-i1/greybox")

## Packages to run
library(nloptr)
library(parallel)
library(snow)
library(clusterGeneration)
library(rlecuyer)
library(smooth)
library(legion)

## Function to run

genmultiseries <- function(model = model, nObs = nObs, nVariate = nVariate,
                           persistence.type = "dependent") {
  
  Sigma <- matrix(c(2, 1,
                    1, 2), ncol = 2)
  
  vInit <- c(100, 120)
  
  if (persistence.type == "dependent") {
    
    matP <- matrix(c(0.4, 0.01,
                     0.01, 0.3), ncol = 2, byrow = TRUE)
    
  } else if (persistence.type == "independent") {
    
    matP <- matrix(c(0.4, 0.0,
                     0.0, 0.3), ncol = 2, byrow = TRUE)
    
  }
  
  y <- sim.ves(model = model, obs = nObs, nsim = 1, nvariate = nVariate,
               persistence = matP, initial = vInit, bounds = "admissible",
               randomizer = "rnorm", mean = rep(0, nVariate), sd = Sigma)$data
  
  return(y)
  
}

scalerVES <- function(distribution="dnorm", Etype, obsInSample, other=NULL,
                      errors, yFitted=NULL, normalizer=1, loss="likelihood"){
  if(loss=="likelihood"){
    scaleValue <- (errors / normalizer) %*% t(errors / normalizer) / obsInSample;
    return(scaleValue*normalizer^2);
  }
  else{
    scaleValue <- diag(rowSums(errors^2) / obsInSample);
    return(scaleValue);
  }
}

dmvnormInternal <- function(q, mean=0, Sigma=1, log=FALSE){
  # The function returns PDF of multivariate normal distribution
  # q should contain obs in columns and series in rows
  if(!is.null(ncol(q))){
    obs <- ncol(q);
    nSeries <- nrow(q);
  }
  else{
    return(dnorm(x=q, mean=mean, sd=sqrt(Sigma), log=log));
  }
  # If dims of mean differ from the q, create the matrix
  # if(!all(dim(q)==dim(mean))){
  #     mean <- matrix(mean,nSeries,obs,byrow=TRUE);
  # }
  # Take invert. If it doesn't work, return NAs
  SigmaInv <- try(chol2inv(chol(Sigma)), silent=TRUE);
  if(inherits(SigmaInv,"try-error")){
    SigmaInv <- try(solve(Sigma, diag(nSeries), tol=1e-20), silent=TRUE);
    if(inherits(SigmaInv,"try-error")){
      return(rep(NA,obs));
    }
  }
  # Defin X and X transposed
  xt <- q - mean;
  x <- t(xt);
  # Calculate PDF
  mvnormReturn <- vector("numeric",obs);
  for(i in 1:obs){
    mvnormReturn[i] <- x[i,,drop=FALSE] %*% SigmaInv %*% xt[,i,drop=FALSE];
  }
  mvnormReturn[] <- exp(-0.5 * mvnormReturn) * (2*pi)^{-nSeries/2} *det(Sigma)^{-0.5};
  if(log){
    mvnormReturn[] <- log(mvnormReturn);
  }
  return(mvnormReturn);
}

norm.offdiag <- function(x) {
  diag(x) <- NA
  return(sqrt(sum(c(x[!is.na(x)])^2)))
}

norm.diag <- function(x) {sqrt(sum(c(diag(x))^2))}

norm.vec <- function(x) {sqrt(sum(c(x)^2))}

loss.shrVES <- function(actual, fitted, B) {
  
  nObs <- nrow(actual)
  nSeries <- ncol(actual)
  nState <- is.character(Etype) + is.character(Ttype) + is.character(Stype)
  
  error <- actual - fitted
  normalizer <- sum(colMeans(abs(diff(t(actual))),na.rm=TRUE))
  
  if (loss.type == "likelihood") {
    scaleValue <- scalerVES("dnorm", Etype, nObs, NULL,
                            error, NULL, normalizer=normalizer,
                            loss="likelihood");
    
    cfRes <- -sum(switch(Etype,
                         "A"=dmvnormInternal(error, 0, scaleValue, log=TRUE),
                         "M"=dmvnormInternal(error, 0, scaleValue, log=TRUE)-
                           colSums(log(actual))));
  } else if (loss.type == "trace") {
    cfRes <- sum(colSums(error^2)) / nObs;
  }
  
  mPersistance <- matrix(B[1:(nState*(nSeries^2))], ncol = nSeries, byrow = TRUE)
  
  if (nState == 1) {
    
    mPersist.Level <- mPersistance[seq(1, nState*nSeries, nState),]
    listPersistance <- list(level = mPersist.Level)
    
  } else if (nState == 2) {
    
    mPersist.Level <- mPersistance[seq(1, nState*nSeries, nState),]
    mPersist.Trend <- mPersistance[seq(2, nState*nSeries, nState),]
    listPersistance <- list(level = mPersist.Level,
                            trend = mPersist.Trend)
    
  } else if (nState == 3) {
    
    mPersist.Level <- mPersistance[seq(1, nState*nSeries, nState),]
    mPersist.Trend <- mPersistance[seq(2, nState*nSeries, nState),]
    mPersist.Season <- mPersistance[seq(3, nState*nSeries, nState),]
    listPersistance <- list(level = mPersist.Level,
                            trend = mPersist.Trend,
                            season = mPersist.Season)
    
  }
  
  if (penalty == "double.penalty") {
    
    norm.penalty <- sapply(listPersistance, norm.vec) + sapply(listPersistance, norm.offdiag)
    
  } else if (penalty == "offdiag.only") {
    
    norm.penalty <- sapply(listPersistance, norm.offdiag)
    
  } else if (penalty == "diag.only") {
    
    norm.penalty <- sapply(listPersistance, norm.diag)
    
  } else if (penalty == "all") {
    
    norm.penalty <- sapply(listPersistance, norm.vec)
    
  }
  
  reg.cfRes <- (1-lambda) * cfRes - lambda * (t(hyperparam) %*% norm.penalty)
  
  return(reg.cfRes)
  
}

cvLambda <- function(y, model = model, end.date = NULL, h = 2, loss.type = "likelihood",
                     Etype, Ttype, Stype, damped, penalty, hyperparam,
                     seq.lambda = seq.lambda) {
  
  model <- model
  end.date <<- end.date
  h <- h
  loss.type <<- loss.type
  seq.lambda <- seq.lambda
  origin <- 5
  
  # global environment
  Etype <<- Etype
  Type <<- Ttype
  Stype <<- Stype
  damped <<- damped
  penalty <<- penalty
  hyperparam <<- hyperparam
  loss.type <<- loss.type
  
  collect.cv <- matrix(NA, ncol = 2, nrow = nLambda)
  
  for (i in 1:nLambda) {
    
    collect.rmse.cv <- matrix(NA, ncol = origin, nrow = 1)
    
    for (k in 1:origin) {
      yTrain <- window(y, end = end.date+k)
      yTest <- window(y, start = end.date+k+1, end = end.date+k+h) # h = 1
      
      lambda <<- seq.lambda[i]
      fit1 <- ves(yTrain, model = model, loss = loss.shrVES, bounds = "admissible",
                  persistence = "dependent", h = h, holdout = TRUE)
      
      collect.rmse.cv[k] <- mean(sqrt(sum((yTest[1,] - fit1$forecast[1,])^2)))
    }
    
    collect.cv[i,] <- c(mean(collect.rmse.cv), ((sd(collect.rmse.cv))^2)/sqrt(origin))
    
  }
  
  return(collect.cv)
  
}

analyseVES <- function(y, model = model, end.date = end.date, h = 2, 
                       loss.type = "likelihood",
                       Etype, Ttype, Stype, damped, penalty, hyperparam, nLambda) {
  
  y <- y
  model <- model
  end.date <<- end.date
  h <- h
  loss.type <<- loss.type
  
  Etype <<- Etype
  Ttype <<- Ttype
  Stype <<- Stype
  damped <<- damped
  penalty <<- penalty
  hyperparam <<- hyperparam
  nLambda <<- nLambda
  
  # sequence of lambdas
  lambdaRatio = 1e-04
  lambdaMax <- 0.9999
  nLambda <- nLambda
  lambdaMin <- lambdaMax * lambdaRatio
  loghi <- log(lambdaMax)
  loglo <- log(lambdaMin)
  logrange <- loghi - loglo
  interval <- -logrange/(nLambda - 1)
  seq.lambda <<- exp(seq.int(from = loghi, to = loglo, by = interval))
  
  collect.lambda <- cvLambda(y, model = model, end.date = end.date, h = 2, 
                             loss.type = loss.type,
                             Etype, Ttype, Stype, damped, penalty, hyperparam,
                             seq.lambda = seq.lambda)
  
  minLambda <- seq.lambda[which.min(collect.lambda[,1])]
  rank2 <- order(collect.lambda[,1])[2]
  oseLambda <- seq.lambda[rank2]
  
  yTrain <- window(y, end = end.date)
  yTest <- window(y, start = end.date+1)
  
  lambda <<- minLambda
  fit.minL <- ves(yTrain, model = model, loss = loss.shrVES, bounds = "admissible",
                  persistence = "dependent", h = h, holdout = TRUE)
  eigen.minL <- eigen(fit.minL$transition - fit.minL$persistence %*% fit.minL$measurement,
                      only.values=TRUE, symmetric=TRUE)$values
  
  lambda <<- oseLambda
  fit.1seL <- ves(yTrain, model = model, loss = loss.shrVES, bounds = "admissible",
                  persistence = "dependent", h = h, holdout = TRUE)
  eigen.1seL <- eigen(fit.1seL$transition - fit.1seL$persistence %*% fit.1seL$measurement,
                      only.values=TRUE, symmetric=TRUE)$values
  
  # benchmark
  fit.bcm1 <- ves(yTrain, model = model, loss = "trace", bounds = "admissible",
                  persistence = "dependent", h = h, holdout = TRUE)
  eigen.bcm1 <- eigen(fit.bcm1$transition - fit.bcm1$persistence %*% fit.bcm1$measurement,
                      only.values=TRUE, symmetric=TRUE)$values
  
  # benchmark
  fit.bcm2 <- ves(yTrain, model = model, loss = "likelihood", bounds = "admissible",
                  persistence = "dependent", h = h, holdout = TRUE)
  eigen.bcm2 <- eigen(fit.bcm2$transition - fit.bcm2$persistence %*% fit.bcm2$measurement,
                      only.values=TRUE, symmetric=TRUE)$values
  
  # benchmark
  fit.bcm3 <- apply(yTrain, 2, 
                    function(x) adam(x, h = h, bounds = "admissible", 
                                     loss = "likelihood"))
  
  eigenValues <- cbind(eigen.minL, eigen.1seL, eigen.bcm1, eigen.bcm2)
  
  error.minL <- yTest[1:h,] - fit.minL$forecast
  error.1seL <- yTest[1:h,] - fit.1seL$forecast
  error.bcm1 <- yTest[1:h,] - fit.bcm1$forecast
  error.bcm2 <- yTest[1:h,] - fit.bcm2$forecast
  error.bcm3 <- yTest[1:h,] - cbind(fit.bcm3$Series1$forecast, fit.bcm3$Series2$forecast)
  
  collect.rmse <- cbind(rowMeans(apply(error.minL^2, 2, cumsum)/(1:h)),
                        rowMeans(apply(error.1seL^2, 2, cumsum)/(1:h)),
                        rowMeans(apply(error.bcm1^2, 2, cumsum)/(1:h)),
                        rowMeans(apply(error.bcm2^2, 2, cumsum)/(1:h)),
                        rowMeans(apply(error.bcm3^2, 2, cumsum)/(1:h)))
  
  collect.me <- cbind(rowMeans(apply(error.minL, 2, cumsum)/(1:h)),
                      rowMeans(apply(error.1seL, 2, cumsum)/(1:h)),
                      rowMeans(apply(error.bcm1, 2, cumsum)/(1:h)),
                      rowMeans(apply(error.bcm2, 2, cumsum)/(1:h)),
                      rowMeans(apply(error.bcm3, 2, cumsum)/(1:h)))
  
  colnames(collect.rmse) <- c("lambda.min", "lambda.1se", "benchmark1", "benchmark2", "benchmark3")
  colnames(collect.me) <- c("lambda.min", "lambda.1se", "benchmark1", "benchmark2", "benchmark3")
  rownames(collect.rmse) <- paste0("t+1-", 1:h)
  rownames(collect.me) <- paste0("t+1-", 1:h)
  
  return(list(acc.rmse = collect.rmse,
              acc.me = collect.me,
              eigenValues = eigenValues,
              data = y,
              seq.lambda = seq.lambda,
              loss.type = loss.type))
} 

compile.all <- function(model = model, nObs = 100, nVariate = 2, matP = matP, vInit = vInit, Sigma = Sigma, persistence.type = "dependent",
                        end.date = end.date, h = 2, loss.type = "likelihood", Etype, Ttype, Stype, damped, penalty, hyperparam, nLambda) {
  
  model <- model
  nObs <- nObs
  nVariate <- nVariate
  persistence.type <- persistence.type
  
  end.date <<- end.date
  h <- h
  loss.type <<- loss.type
  
  Etype <<- Etype
  Ttype <<- Ttype
  Stype <<- Stype
  damped <<- damped
  penalty <<- penalty
  hyperparam <<- hyperparam
  nLambda <<- nLambda
  
  y <- genmultiseries(model = "ANN", nObs = 100, nVariate = 2, persistence.type = persistence.type)
  
  result <- analyseVES(y, model = model, end.date = 20, h = 12, loss.type = loss.type,
                       Etype = Etype, Ttype = Ttype, Stype = Stype, damped = damped,
                       penalty = penalty, hyperparam = hyperparam, nLambda = nLambda)
  
  return(result)
  
}


# Running the parallel ----------------------------------------------------

crs <- detectCores()
cl <- makeCluster(getOption("cl.cores", crs))
writeLines(paste("Running with", crs, 'cores'))
# Load packages to cluster
invisible(clusterCall(cl, function(pkgs) {
  library(clusterGeneration)
  library(smooth)
  library(legion)
  library(nloptr)
  library(snow)
}))

invisible(clusterExport(cl, "genmultiseries"))
invisible(clusterExport(cl, "scalerVES"))
invisible(clusterExport(cl, "dmvnormInternal"))
invisible(clusterExport(cl, "norm.offdiag"))
invisible(clusterExport(cl, "norm.diag"))
invisible(clusterExport(cl, "norm.vec"))
invisible(clusterExport(cl, "loss.shrVES"))
invisible(clusterExport(cl, "cvLambda"))
invisible(clusterExport(cl, "analyseVES"))
invisible(clusterExport(cl, "compile.all"))

clusterSetupRNG(cl, seed = 020193)

runs <- 5

# with Trace, dependent persistence
system.time({dgpANN_modelANN_S20_OffDiag_TR_PDP <- clusterApplyLB(cl, 1:runs, function(x) compile.all(model = "ANN", nObs = 36, nVariate = 2, matP = matP, vInit = vInit, Sigma = Sigma, persistence.type = "dependent",
                                                                                                      end.date = 20, h = 12, loss.type = "trace", Etype = "A", Ttype = NULL, Stype = NULL,
                                                                                                      damped = FALSE, penalty = "offdiag.only", hyperparam = 1, nLambda = 100))})

system.time({dgpANN_modelANN_S240_OffDiag_TR_PDP <- clusterApplyLB(cl, 1:runs, function(x) compile.all(model = "ANN", nObs = 240, nVariate = 2, matP = matP, vInit = vInit, Sigma = Sigma, persistence.type = "dependent",
                                                                                                       end.date = 20, h = 12, loss.type = "trace", Etype = "A", Ttype = NULL, Stype = NULL,
                                                                                                       damped = FALSE, penalty = "offdiag.only", hyperparam = 1, nLambda = 100))})

# with Likelihood, dependent persistence
system.time({dgpANN_modelANN_S20_OffDiag_LH_PDP <- clusterApplyLB(cl, 1:runs, function(x) compile.all(model = "ANN", nObs = 36, nVariate = 2, matP = matP, vInit = vInit, Sigma = Sigma, persistence.type = "dependent",
                                                                                                      end.date = 20, h = 12, loss.type = "likelihood", Etype = "A", Ttype = NULL, Stype = NULL,
                                                                                                      damped = FALSE, penalty = "offdiag.only", hyperparam = 1, nLambda = 100))})

system.time({dgpANN_modelANN_S240_OffDiag_LH_PDP <- clusterApplyLB(cl, 1:runs, function(x) compile.all(model = "ANN", nObs = 240, nVariate = 2, matP = matP, vInit = vInit, Sigma = Sigma, persistence.type = "dependent",
                                                                                                       end.date = 20, h = 12, loss.type = "likelihood", Etype = "A", Ttype = NULL, Stype = NULL,
                                                                                                       damped = FALSE, penalty = "offdiag.only", hyperparam = 1, nLambda = 100))})


# with Trace, independent persistence
system.time({dgpANN_modelANN_S20_OffDiag_TR_PID <- clusterApplyLB(cl, 1:runs, function(x) compile.all(model = "ANN", nObs = 36, nVariate = 2, matP = matP, vInit = vInit, Sigma = Sigma, persistence.type = "independent",
                                                                                                      end.date = 20, h = 12, loss.type = "trace", Etype = "A", Ttype = NULL, Stype = NULL,
                                                                                                      damped = FALSE, penalty = "offdiag.only", hyperparam = 1, nLambda = 100))})

system.time({dgpANN_modelANN_S240_OffDiag_TR_PID <- clusterApplyLB(cl, 1:runs, function(x) compile.all(model = "ANN", nObs = 240, nVariate = 2, matP = matP, vInit = vInit, Sigma = Sigma, persistence.type = "independent",
                                                                                                       end.date = 20, h = 12, loss.type = "trace", Etype = "A", Ttype = NULL, Stype = NULL,
                                                                                                       damped = FALSE, penalty = "offdiag.only", hyperparam = 1, nLambda = 100))})

# with Likelihood, independent persistence
system.time({dgpANN_modelANN_S20_OffDiag_LH_PID <- clusterApplyLB(cl, 1:runs, function(x) compile.all(model = "ANN", nObs = 36, nVariate = 2, matP = matP, vInit = vInit, Sigma = Sigma, persistence.type = "independent",
                                                                                                      end.date = 20, h = 12, loss.type = "likelihood", Etype = "A", Ttype = NULL, Stype = NULL,
                                                                                                      damped = FALSE, penalty = "offdiag.only", hyperparam = 1, nLambda = 100))})

system.time({dgpANN_modelANN_S240_OffDiag_LH_PID <- clusterApplyLB(cl, 1:runs, function(x) compile.all(model = "ANN", nObs = 240, nVariate = 2, matP = matP, vInit = vInit, Sigma = Sigma, persistence.type = "independent",
                                                                                                       end.date = 20, h = 12, loss.type = "likelihood", Etype = "A", Ttype = NULL, Stype = NULL,
                                                                                                       damped = FALSE, penalty = "offdiag.only", hyperparam = 1, nLambda = 100))})

stopCluster(cl)
