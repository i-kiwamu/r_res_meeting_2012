### colors
k.red <- "#FF0000"
k.yellow <- "darkgoldenrod1"
k.blue <- "#0080FF"
k.orange <- "#EB9507"
k.green <- "darkgreen"
k.purple <- "#FF00FF"
k.cyan <- "#00FFFF"
k.orange.back <- "#ffe5cc"
k.green.back <- "#ccffcc"
k.cyan.back <- "#ccffff"
k.blue.back <- "#cce6ff"
k.purble.back <- "#ffccff"
k.red.back <- "#ffcccc"
k.yellow.back <- "#ffffcc"


### RStan
library(rstan)
schools_dat <- list(J = 8,
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
fit1 <- stan(file='8schools.stan', data=schools_dat,
             iter=10000, chains=4)

par(cex=1.5)
plot(fit1)
par(cex=1.5)
rstan::traceplot(fit1)

fit1.coef <- NULL
fit1.coef$eta1 <- extract(fit1, "eta[1]", inc_warmup=TRUE)
fit1.coef$theta1 <- extract(fit1, "theta[1]", inc_warmup=TRUE)

### JAGS
library(R2jags)
schools_init <- function(){
    list(mu = runif(1),
         tau = runif(1),
         eta = rnorm(8))
}
schools_params <- c("mu", "tau", "eta", "theta")
fit2 <- jags(schools_dat, schools_init, schools_params,
             model.file="8schools.jags", n.iter=10000, n.chain=4,
             n.burnin=0, n.thin=1)
fit2.coef <- NULL
attach.jags(fit2)
fit2.coef$eta1 <- eta[,1]
fit2.coef$theta1 <- theta[,1]
detach.jags()

### Nile data
N <- length(Nile)

## dlm package
library(dlm)
system.time({
    nile.dlm <- dlmGibbsDIG(Nile, mod=dlmModPoly(1),
                            shape.y=0.1, rate.y=0.1,
                            shape.theta=0.1, rate.theta=0.1,
                            n.sample=8000)
})
#   user  system elapsed
#115.857   0.145 115.993
burn <- 1:4000
sigma.dlm <- sqrt(matrix(unlist(nile.dlm[c("dV", "dW")]),
                         ncol=2)[-burn,])
colnames(sigma.dlm) <- c("sigma y", "sigma theta")
sigma.mcmc.dlm <- mcmc.list(mcmc(sigma.dlm))
plot(sigma.mcmc.dlm)


## JAGS
nile.1 <- Nile[1]
nile.data.jags <- list("N", "Nile", "nile.1")
nile.inits.jags <- function() {
    list(theta=rnorm(N), sigma.theta=runif(1), sigma.y=runif(1),
         F=rnorm(1), G=rnorm(1))
}
nile.params.jags <- c("theta", "sigma.y", "sigma.theta", "F", "G")
system.time({
    nile.jags <- jags(nile.data.jags, nile.inits.jags,
                      nile.params.jags, model.file="nile.jags",
                      n.iter=8000, n.burnin=4000, n.chains=4)
})
#   user  system elapsed
#  8.595   0.263  11.647
sigma.mcmc.jags <-
    as.mcmc(lapply(1:4, function(i) as.mcmc(nile.jags)[[i]][,c(5,4)]))
plot(sigma.mcmc.jags)


## RStan
nile.data.stan <- list(N=N, Nile=Nile)
system.time({
    nile.stan <- stan(file="nile.stan", data=nile.data.stan,
                      iter=8000)
})
#   user  system elapsed
# 24.355   0.449  24.944
system.time({
    nile.stan <- stan(fit=nile.stan, data=nile.data.stan,
                      iter=8000)
})
nile.stan.extract <- extract(nile.stan, inc_warmup=FALSE)
sigma.mcmc.stan <- mcmc.list(mcmc(nile.stan.extract[,1,103:104]),
                             mcmc(nile.stan.extract[,2,103:104]),
                             mcmc(nile.stan.extract[,3,103:104]),
                             mcmc(nile.stan.extract[,4,103:104]))
plot(sigma.mcmc.stan)


## how converged?
nile.stan.extract.all <- extract(nile.stan, inc_warmup=TRUE)
plot(1, xlim=c(50, 2500), ylim=c(0, 400), type="n",
     xlab="sigma y", ylab="sigma theta")
abline(v=119.5, lty=2)
abline(h=50.7, lty=2)
xx <- seq(-50, 3000, length=1000)
yy <- seq(-50, 450, length=1000)
polygon(c(rep(92.8, 1000), rep(144.3, 1000)), c(yy, rev(yy)),
        border=NA, col=rgb(1.0, 0.8, 0.8, alpha=0.75))
polygon(c(xx, rev(xx)), c(rep(26.4, 1000), rep(92.6, 1000)),
        border=NA, col=rgb(0.8, 0.9, 1.0, alpha=0.75))

points(nile.stan.extract.all[1:50,1,103],
       nile.stan.extract.all[1:50,1,104],
       type="o", col=rgb(0, 0.5, 1, alpha=0.75))
points(nile.stan.extract.all[1:50,2,103],
       nile.stan.extract.all[1:50,2,104],
       type="o", col=rgb(0.18, 0.56, 0.15, alpha=0.75))
points(nile.stan.extract.all[1:50,3,103],
       nile.stan.extract.all[1:50,3,104],
       type="o", col=rgb(0.98, 0.71, 0.15, alpha=0.75))
points(nile.stan.extract.all[1:50,4,103],
       nile.stan.extract.all[1:50,4,104],
       type="o", col=rgb(1, 0, 0, alpha=0.75))

points(nile.stan.extract.all[51:500,1,103],
       nile.stan.extract.all[51:500,1,104],
       type="o", col=rgb(0, 0.5, 1, alpha=0.75))
points(nile.stan.extract.all[51:500,2,103],
       nile.stan.extract.all[51:500,2,104],
       type="o", col=rgb(0.18, 0.56, 0.15, alpha=0.75))
points(nile.stan.extract.all[51:500,3,103],
       nile.stan.extract.all[51:500,3,104],
       type="o", col=rgb(0.98, 0.71, 0.15, alpha=0.75))
points(nile.stan.extract.all[51:500,4,103],
       nile.stan.extract.all[51:500,4,104],
       type="o", col=rgb(1, 0, 0, alpha=0.75))

points(nile.stan.extract.all[501:5000,1,103],
       nile.stan.extract.all[501:5000,1,104],
       type="o", col=rgb(0, 0.5, 1, alpha=0.75))
points(nile.stan.extract.all[501:5000,2,103],
       nile.stan.extract.all[501:5000,2,104],
       type="o", col=rgb(0.18, 0.56, 0.15, alpha=0.75))
points(nile.stan.extract.all[501:5000,3,103],
       nile.stan.extract.all[501:5000,3,104],
       type="o", col=rgb(0.98, 0.71, 0.15, alpha=0.75))
points(nile.stan.extract.all[501:5000,4,103],
       nile.stan.extract.all[501:5000,4,104],
       type="o", col=rgb(1, 0, 0, alpha=0.75))

## how fit?
nile.stan.coef <- NULL
nile.stan.coef$theta <-
    apply(matrix(extract(nile.stan, "theta", inc_warmup=FALSE),
                 ncol=100), 2, median)
nile.stan.coef$F <- median(extract(nile.stan, "F", inc_warmup=FALSE))
nile.stan.coef$nile.est <-
    mcmc(matrix(extract(nile.stan, "F", inc_warmup=FALSE)[,,1], ncol=1)[,1] * matrix(extract(nile.stan, "theta", inc_warmup=FALSE), ncol=100))
nile.stan.hpd <- HPDinterval(nile.stan.coef$nile.est)
plot(Nile, type="n", cex=1.5, cex.axis=1.5, cex.lab=1.5)
polygon(c(1871:1970, 1970:1871),
        c(nile.stan.hpd[,1], rev(nile.stan.hpd[,2])),
        border=NA, col=rgb(1.0, 0.8, 0.8, alpha=0.75))
lines(1871:1970, nile.stan.coef$F * nile.stan.coef$theta,
      col=k.blue, lwd=1.5)
points(1871:1970, Nile, type="o", cex=1.5)


## how wide?
nile.inits.jags2 <- function() {
    list(theta=Nile, sigma.theta=50, sigma.y=120,
         F=rnorm(1, mean=1, sd=0.1), G=rnorm(1, mean=1, sd=0.1))
}
nile.jags <- jags(nile.data.jags, nile.inits.jags2,
                  nile.params.jags, model.file="nile.jags",
                  n.iter=8000, n.burnin=0, n.chains=4, n.thin=1)
nile.jags.coef <- NULL
attach.jags(nile.jags)
nile.jags.coef$theta2 <- theta[,2]
nile.jags.coef$theta3 <- theta[,3]
detach.jags()

nile.inits.stan <- function(){
    list(theta=Nile, sigma_theta=50, sigma_y=120,
         F=rnorm(1, mean=1, sd=0.1), G=rnorm(1, mean=1, sd=0.1))
}
nile.stan <- stan(fit=nile.stan, data=nile.data.stan,
                  inits=nile.inits.stan, iter=8000)
nile.stan.coef <- NULL
nile.stan.coef$theta2 <- extract(nile.stan, "theta[2]", inc_warmup=TRUE)
nile.stan.coef$theta3 <- extract(nile.stan, "theta[3]", inc_warmup=TRUE)

par(mai=c(1.02, 1.02, 0.82, 0.42))
plot(nile.stan.coef$theta2, nile.stan.coef$theta3,
     xlab=expression(theta[2]), ylab=expression(theta[3]),
     cex=.1, cex.lab=1.5, cex.axis=1.5)
abline(v = Nile[2], h=Nile[3], lty=2)
points(nile.jags.coef$theta2, nile.jags.coef$theta3, cex=.1, col="red")
par(mai=c(1.02, 0.82, 0.82, 0.42))



### Inverse Gaussian
y <- rinvGauss(100, nu=5, lambda=50)
inv.gauss.data <- list(N = 100, y = y)
inv.gauss.stan <- stan(file="inv_gauss.stan", data=inv.gauss.data)
