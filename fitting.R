
Sigma <- matrix(c(2, 1,
                  1, 2), ncol = 2)
matP <- matrix(c(0.4, 0.0,
                 0.0, 0.3), ncol = 2, byrow = TRUE)
vInit <- c(100, 120)

Etype <<- "A"
Ttype <<- NULL
Stype <<- NULL
damped <<- FALSE
penalty <<- "offdiag.only"
lambda <<- 0.01
hyperparam <<- 1

y <- sim.ves(model = "ANN", obs = 24, nsim = 1, nvariate = 2,
             persistence = matP, initial = vInit, bounds = "admissible",
             randomizer = "rnorm", mean = rep(0, 2), sd = Sigma)$data

fit1 <- ves(y, model = "ANN", loss = loss.shrVES, bounds = "admissible",
            persistence = "dependent", h = 4, holdout = TRUE)
fit2 <- ves(y, model = "ANN", loss = "likelihood", bounds = "admissible",
            persistence = "dependent", h = 4, holdout = TRUE)

c(fit1$lossValue, fit2$lossValue)
fit1$persistence
fit2$persistence

eigen(fit1$transition - fit1$persistence %*% fit1$measurement, only.values=TRUE, symmetric=TRUE)$value
eigen(fit2$transition - fit2$persistence %*% fit2$measurement, only.values=TRUE, symmetric=TRUE)$value

mean(sqrt(colMeans((fit1$holdout - fit1$forecast)^2)))/mean(sqrt(colMeans((fit2$holdout - fit2$forecast)^2)))


# real data
fit1 <- ves(yTrain[,1:2], model = "ANN", loss = loss.shrVES, bounds = "admissible",
            persistence = "dependent", h = 4, holdout = TRUE)
fit2 <- ves(yTrain[,1:2], model = "ANN", loss = "likelihood", bounds = "admissible",
            persistence = "dependent", h = 4, holdout = TRUE)

c(fit1$lossValue, fit2$lossValue)
fit1$persistence
fit2$persistence

eigen(fit1$transition - fit1$persistence %*% fit1$measurement, only.values=TRUE, symmetric=TRUE)$value
eigen(fit2$transition - fit2$persistence %*% fit2$measurement, only.values=TRUE, symmetric=TRUE)$value

mean(sqrt(colMeans((fit1$holdout - fit1$forecast)^2)))/mean(sqrt(colMeans((fit2$holdout - fit2$forecast)^2)))

