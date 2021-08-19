#Enjoy the joy of Copulas
install.packages("scatterplot3d")
library(copula)
library(mvtnorm)
library(scatterplot3d)
set.seed(1)
#classes
myCop.norm <- ellipCopula(family = "normal", dim = 3, dispstr = "ex",  param = 0.4)
myCop.t <- ellipCopula(family = "t", dim = 3, dispstr = "toep",  param = c(0.8, 0.5), df = 8)
myCop.clayton <- archmCopula(family = "clayton", dim = 3, param = 2)
#the mvdc class
myMvd <- mvdc(copula = myCop.clayton, margins = c("norm", "norm",
                                                  "norm"), paramMargins = list(list(mean = 0, sd = 2), list(mean = 0, sd = 1), list(mean = 0, sd = 2)))
u <- rCopula(4,myCop.t)
u
cbind(dCopula(u, myCop.t), pCopula(u, myCop.t))

x <- rMvdc(4,myMvd)
x
cbind(dMvdc(x, myMvd), pMvdc(x, myMvd))

##graphics
par(mfrow = c(1, 2), mar = c(2, 2, 1, 1), oma = c(1, 1, 0, 0),
    mgp = c(2, 1, 0))
u <- rCopula(200,myCop.norm)
scatterplot3d(u)
v <- rCopula(200,myCop.norm)
scatterplot3d(v)

myMvd1 <- mvdc(copula = archmCopula(family = "clayton", param = 2),
                  margins = c("norm", "norm"), paramMargins = list(list(mean = 0,
                                                                          sd = 1), list(mean = 0, sd = 1)))
myMvd2 <- mvdc(copula = archmCopula(family = "frank", param = 5.736),
                   margins = c("norm", "norm"), paramMargins = list(list(mean = 0,
                                                                           sd = 1), list(mean = 0, sd = 1)))
myMvd3 <- mvdc(copula = archmCopula(family = "gumbel", param = 2),
                  margins = c("norm", "norm"), paramMargins = list(list(mean = 0,
                                                                           sd = 1), list(mean = 0, sd = 1)))
par(mfrow = c(1, 3), mar = c(2, 2, 1, 1), oma = c(1, 1, 0, 0),
        mgp = c(2, 1, 0))
contour(myMvd1, dMvdc, xlim = c(-3, 3), ylim = c(-3, 3))
contour(myMvd2, dMvdc, xlim = c(-3, 3), ylim = c(-3, 3))
contour(myMvd3, dMvdc, xlim = c(-3, 3), ylim = c(-3, 3))

#fit a model
 myMvd <- mvdc(copula = ellipCopula(family = "normal", param = 0.5),
                 margins = c("gamma", "gamma"), paramMargins = list(list(shape = 2,
                  scale = 1), list(shape = 3, scale = 2)))
n <- 200
dat <- rMvdc(n,myMvd)

loglikMvdc(c(2, 1, 3, 2, 0.5), dat, myMvd)

 mm <- apply(dat, 2, mean)
 vv <- apply(dat, 2, var)
 b1.0 <- c(mm[1]^2/vv[1], vv[1]/mm[1])
 b2.0 <- c(mm[2]^2/vv[2], vv[2]/mm[2])
 a.0 <- sin(cor(dat[, 1], dat[, 2], method = "kendall") * pi/2)
 start <- c(b1.0, b2.0, a.0)
 fit <- fitMvdc(dat, myMvd, start = start,
 optim.control = list(trace = TRUE, maxit = 2000))
fit

loglik.marg <- function(b, x) sum(dgamma(x, shape = b[1], scale = b[2], log = TRUE))
ctrl <- list(fnscale = -1)
b1hat <- optim(b1.0, fn = loglik.marg, x = dat[, 1], control = ctrl)$par
b2hat <- optim(b2.0, fn = loglik.marg, x = dat[, 2], control = ctrl)$par
udat <- cbind(pgamma(dat[, 1], shape = b1hat[1], scale = b1hat[2]),
                  pgamma(dat[, 2], shape = b2hat[1], scale = b2hat[2]))
fit.ifl <- fitCopula(myMvd@copula,udat, start = a.0) 

c(b1hat, b2hat, fit.ifl@estimate)
fit.ifl

eu <- cbind((rank(dat[, 1]) - 0.5)/n, (rank(dat[, 2]) - 0.5)/n)
fit.cml <- fitCopula(myMvd@copula, eu,start = a.0)
fit.cml
