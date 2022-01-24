load("scandUnemp.Rdata")

library(smooth)
library(greybox)
library(tsutils)
library(corrplot)
library(legion)

library("devtools")
install_github("config-i1/legion")

# library(legion)

ts.plot(Y.bottom)
corrplot(cor(Y.bottom))

yTrain <- window(Y.bottom, end = 2013.917, frequency = 12)
yTest <- window(Y.bottom, start = 2014.00, frequency = 12)

collect.univariate.model <- NULL
for (i in 1:16) {
  collect.univariate.model <- c(collect.univariate.model, adam(yTrain[,i])$model)
}


fit.P.IDP <- ves(yTrain, model = "ANN", persistence = "independent", 
                 transition = "independent", loss = "likelihood")
fit.P.ALL <- ves(yTrain, model = "ANN", persistence = "dependent", 
                 transition = "independent", loss = "likelihood")

fit.P.ALL$B[1:ncol(yTrain)^2]

fit.P.ALL$B[(ncol(yTrain)^2 + 1):length(fit.P.ALL$B)]

error.P.IDP <- yTest - forecast(fit.P.IDP, h = 12)$mean
error.P.ALL <- yTest - forecast(fit.P.ALL, h = 12)$mean

fit.P.IDP.S <- ves(yTrain, model = "ANA", persistence = "independent", 
                   transition = "independent", loss = "likelihood", initialSeason = "common")
fit.P.ALL.S <- ves(yTrain, model = "ANA", persistence = "dependent", 
                   transition = "independent", loss = "likelihood", initialSeason = "common")

error.P.IDP.S <- yTest - forecast(fit.P.IDP.S, h = 12)$mean
error.P.ALL.S <- yTest - forecast(fit.P.ALL.S, h = 12)$mean

mean(apply(error.P.IDP, 2, function(x) sqrt(mean(x^2)))/ apply(yTrain, 2, function(x) mean(abs(x))))
mean(apply(error.P.ALL, 2, function(x) sqrt(mean(x^2)))/ apply(yTrain, 2, function(x) mean(abs(x))))

mean(apply(error.P.IDP.S, 2, function(x) sqrt(mean(x^2)))/ apply(yTrain, 2, function(x) mean(abs(x))))
mean(apply(error.P.ALL.S, 2, function(x) sqrt(mean(x^2)))/ apply(yTrain, 2, function(x) mean(abs(x))))

mean(apply(error.P.IDP, 2, function(x) mean(x))/ apply(yTrain, 2, function(x) mean(abs(x))))
mean(apply(error.P.ALL, 2, function(x) mean(x))/ apply(yTrain, 2, function(x) mean(abs(x))))

mean(apply(error.P.IDP.S, 2, function(x) mean(x))/ apply(yTrain, 2, function(x) mean(abs(x))))
mean(apply(error.P.ALL.S, 2, function(x) mean(x))/ apply(yTrain, 2, function(x) mean(abs(x))))

fit.P.IDP.S$initialSeason
fit.P.ALL.S$initialSeason

fit.P.ALL.S$persistence

