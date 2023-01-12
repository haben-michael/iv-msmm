## 40 check formulas (p.26 reverse)

## 40a dnormmin
require(mvtnorm)
source('misc.R')
n <- 1e5
mean <- runif(2)
sigma <- with(list(root=matrix(4*runif(4),nrow=2)),root%*%t(root))
x <- rmvnorm(n,mean=mean,sigma=sigma)
op <- par(mfrow=c(1,2))
hist(pmin(x[,1],x[,2]),prob=TRUE)
curve(dminnorm(x,mean=mean,sd=sqrt(diag(sigma)),cov=sigma[1,2]),add=TRUE)
plot(ecdf(pmin(x[,1],x[,2])))
curve(pminnorm(x,mean=mean,sd=sqrt(diag(sigma)),cov=sigma[1,2]),add=TRUE,col='red',lty=2)
par(op)

## 40b 
alpha <- runif(1); beta <- runif(1); gamma <- runif(1)
## integrate(function(u)dnorm(alpha*u)*dnorm(beta*u)*dnorm(gamma*u),-Inf,Inf)$val
## 1/(2*pi*sqrt(alpha^2+beta^2+gamma^2))
## integrate(function(u)u*dnorm(alpha*u)*dnorm(gamma*u)*pnorm(beta*u),-Inf,Inf)$val
## 1/2/pi*beta/(alpha^2+gamma^2)/sqrt(alpha^2+beta^2+gamma^2)
## integrate(function(alpha)1/2/pi*beta/(alpha^2+gamma^2)/sqrt(alpha^2+beta^2+gamma^2),5,10)$val
## with(list(alpha=c(5,10)), diff(1/gamma/2/pi*atan(alpha*beta/gamma/sqrt(alpha^2+beta^2+gamma^2))
## ))
integrate(function(u)pnorm(alpha*u)*pnorm(beta*u)*dnorm(gamma*u),-Inf,Inf)$val
1/gamma/2/pi*atan(alpha*beta/gamma/sqrt(alpha^2+beta^2+gamma^2)) + 1/4/gamma
sum(replicate(1e2,{alpha <- runif(1); beta <- runif(1); gamma <- runif(1)
;integrate(function(u)pnorm(alpha*u)*pnorm(beta*u)*dnorm(gamma*u),-Inf,Inf)$val-(1/gamma/2/pi*atan(alpha*beta/gamma/sqrt(alpha^2+beta^2+gamma^2)) + 1/4/gamma)}))


samps <- replicate(1e2,{alpha <- runif(1); beta <- runif(1); gamma <- runif(1)
;c(integrate(function(u)pnorm(alpha*u)*pnorm(beta*u)*dnorm(gamma*u),-Inf,Inf)$val,
   1/gamma/2/pi*atan(alpha*beta/gamma/sqrt(alpha^2+beta^2+gamma^2)))})
plot(samps[1,],samps[2,]);abline(a=0,b=1);abline(lm(samps[2,] ~ samps[1,]))
lm(samps[2,] ~ samps[1,])

alpha <- runif(1); beta <- runif(1); gamma <- runif(1)
deriv <- function(alpha,beta,gamma)1/2/pi*beta/(alpha^2+gamma^2)/sqrt(alpha^2+beta^2+gamma^2)
## deriv(alpha,beta,gamma)
## integrate(function(u)u*dnorm(alpha*u)*dnorm(gamma*u)*pnorm(beta*u),-Inf,Inf)$val
start <- runif(1); end <- runif(1)
integrate(function(alpha)deriv(alpha,beta,gamma),start,end)$val
integrate(function(u)pnorm(end*u)*pnorm(beta*u)*dnorm(gamma*u),-Inf,Inf)$val-integrate(function(u)pnorm(start*u)*pnorm(beta*u)*dnorm(gamma*u),-Inf,Inf)$val
with(list(alpha=c(start,end)), diff(1/gamma/2/pi*atan(alpha*beta/gamma/sqrt(alpha^2+beta^2+gamma^2))))

beta <- runif(1); gamma <- runif(1);alpha <- runif(1);
integrate(function(u)pnorm(alpha*u)*pnorm(beta*u)*dnorm(gamma*u),-Inf,Inf)$val-1/gamma/2/pi*atan(alpha*beta/gamma/sqrt(alpha^2+beta^2+gamma^2))-1/4/gamma
c(integrate(function(u)pnorm(alpha*u)*pnorm(beta*u)*dnorm(gamma*u),-Inf,Inf)$val,1/gamma/2/pi*atan(alpha*beta/gamma/sqrt(alpha^2+beta^2+gamma^2))+1/4/gamma)



## 40c
source('misc.R')
## supp.S <- c(1,5)
## dS <- function(s)dunif(s,supp.S[1],supp.S[2])
## pS <- function(q)punif(q,supp.S[1],supp.S[2])
## rS <- function(n)runif(n,supp.S[1],supp.S[2])
## E.S2 <- diff(supp.S)^2/12 + mean(supp.S)^2
attach(S.env(supp.S=c(1,5)))
n <- 5e3
s <- rS(n)
w1.var <- 2*E.S2^2/(s[1]-s[2])^2-2*E.S2/n+sum(s[1:4]^2)/n^2
w2.var <- 2*E.S2^2/(s[3]-s[4])^2-2*E.S2/n+sum(s[1:4]^2)/n^2
w.var <- c(w1.var,w2.var)
cov.w <- -2*E.S2/n+sum(s[1:4]^2)/n^2
rho.w <- cov.w/sqrt(w1.var*w2.var)
## sqrt(w.var)/2/pi*atan(sqrt(w.var)*sqrt(n)/sqrt(1/n*sum(s[5:n]^2))
## 1/2*pi/gamma*atan(alpha*beta/gamma/sqrt(alpha^2+beta^2+gamma^2))+1/4/gamma
E.1234.cond.S <- integrate(function(u)pnorm(sqrt(n)/sqrt(1/n*sum(s[5:n]^2))*u)*dminnorm(u,sd=sqrt(c(w1.var,w2.var)),cov=cov.w),-Inf,Inf)$val
n*(E.1234.cond.S-1/4)
E.1234.cond.S
alpha <- sqrt(n)/sqrt(1/n*sum(s[5:n]^2)); gamma <- 1/sqrt(w.var); beta <- 1/sqrt(1-rho.w^2)*(rho.w/sqrt(w.var)-1/sqrt(rev(w.var)))
1/sqrt(w1.var)*integrate(function(u)pnorm(alpha*u)*pnorm(beta[1]*u)*dnorm(gamma[1]*u),-Inf,Inf)$val + 1/sqrt(w2.var)*integrate(function(u)pnorm(alpha*u)*pnorm(beta[2]*u)*dnorm(gamma[2]*u),-Inf,Inf)$val
n*(sum((1/2/pi*atan(alpha*beta*sqrt(w.var)/sqrt(alpha^2+beta^2+gamma^2))+1/4))-1/4)
n*(sum((1/2/pi*atan(beta*sqrt(w.var)/sqrt(1+(beta^2+1/w.var)/alpha^2))+1/4))-1/4)
n*(1/2/pi*sum(atan((rho.w-sqrt(w.var/rev(w.var)))/sqrt(1-rho.w^2)/sqrt(1+  (beta^2+gamma^2)/alpha^2  )  ))+1/2-1/4)
n*(1/2/pi*sum(atan((rho.w-sqrt(w.var/rev(w.var)))/sqrt(1-rho.w^2)/sqrt(1+       (beta^2+1/w.var)/alpha^2            )  ))+1/2-1/4)
n*(1/2/pi*sum(atan((rho.w-sqrt(w.var/rev(w.var)))/sqrt(1-rho.w^2)/sqrt(1+       E.S2/n*(sum(1/w.var)-2*rho.w/sqrt(prod(w.var)))            )  ))+1/2-1/4)
n*(1/2/pi*sum(atan((rho.w-sqrt(w.var/rev(w.var)))/sqrt(1-rho.w^2)/sqrt(1+       E.S2/n*(sum(1/w.var))            )  ))+1/2-1/4)
u <- 1/n; k <- E.S2; delta <- c((s[1]-s[2]),(s[3]-s[4]))
cov.1234 <- 1/2/pi*sum((1+rev(delta^2)/delta^2)^(-1)*( -(1/2*prod(abs(delta))/k+1/2/k*abs(rev(delta^3))/abs(delta)) + 1/4/k*abs(rev(delta))/abs(delta)*sum(delta^2)))
cov.1234

n*(1/2/pi*sum(atan(-(2*E.S2/n/sqrt(prod(w.var))+sqrt(w.var/rev(w.var)))/sqrt(      (1+E.S2/n*(sum(1/w.var)))*(1-4*E.S2^2/prod(w.var)/n^2)            )  ))+1/2-1/4)
n*(1/2/pi*sum(atan(-     (2*E.S2/n/sqrt(prod(w.var))+sqrt(w.var/rev(w.var)))  /    sqrt(1+E.S2/n*(sum(1/w.var)))       ))+1/2-1/4)
1/2/pi*sum(1/(1+w.var/rev(w.var)) * ( -2*E.S2/sqrt(prod(w.var)) + sqrt(w.var/rev(w.var))*(-1/2)*E.S2*(sum(1/w.var))))

1/2/pi*sum(atan(-     (2*E.S2/n/sqrt(prod(w.var))+sqrt(w.var/rev(w.var)))  /    sqrt(1+E.S2/n*(sum(1/w.var)))       ))
-1/4 + 1/2/pi*1/n*sum( 1/(1+w.var/rev(w.var)) * (-2*E.S2/sqrt(prod(w.var)) + sqrt(w.var/rev(w.var))*(-1/2)*E.S2*sum(1/w.var)))

u <- 1/n; k <- E.S2; delta <- c((s[1]-s[2]),(s[3]-s[4]))
sum(atan(-(2*E.S2/n/sqrt(prod(w.var))+sqrt(w.var/rev(w.var)))/sqrt(      (1+E.S2/n*(sum(1/w.var)))*(1-4*E.S2^2/prod(w.var)/n^2)            )  ))+pi/2
n*(sum(atan(- (    u/sqrt(k^2/prod(delta^2) - u*k*sum(1/delta^2)) + sqrt((k/delta^2-u)/(k/rev(delta^2)-u))   ) / sqrt(      (1+k*u*(sum(1/w.var)))*(1-4*k^2/prod(w.var)*u^2)            )  ))+pi/2)
n*(sum(atan(- (    u/sqrt(k^2/prod(delta^2) - u*k*sum(1/delta^2)) + sqrt((k/delta^2-u)/(k/rev(delta^2)-u))   ) / sqrt(      (1+k*u*(1/(2*k)*sum(1/(k/delta^2-u))
))        )  ))+pi/2)
n*(sum(atan(- (    u/sqrt(k^2/prod(delta^2) - u*k*sum(1/delta^2)) + sqrt((k/delta^2-u)/(k/rev(delta^2)-u))   ) / sqrt(      (1+u/2*sum(1/(k/delta^2-u))
)        )  ))+pi/2)

f <- Vectorize(function(u)(sum(1/2/pi*atan(- (    u/sqrt(k^2/prod(delta^2) - u*k*sum(1/delta^2)) + sqrt((k/delta^2-u)/(k/rev(delta^2)-u))   ) / sqrt(      (1+u/2*sum(1/(k/delta^2-u)))        )  ))+1/4))
curve(f,-1/2,1/2)
slope <- 1/2/pi*sum((1+rev(delta^2)/delta^2)^(-1)*( -(1/2*prod(abs(delta))/k+1/2/k*abs(rev(delta^3))/abs(delta)) + 1/4/k*abs(rev(delta))/abs(delta)*sum(delta^2)))
abline(a=f(0),b=slope,col='red')

n*f(1/n)
slope

n*(1/2/pi*sum(atan((rho.w-sqrt(w.var/rev(w.var)))/sqrt(1-rho.w^2)/sqrt(1+  1e-10 )  ))+1/2-1/4)
n*(1/2/pi*sum(atan((rho.w-sqrt(w.var/rev(w.var)))/sqrt(1-rho.w^2)/sqrt(1+  1e-7  )  ))+1/2-1/4)
n*(1/2/pi*sum(atan((rho.w-sqrt(w.var/rev(w.var)))/sqrt(1-rho.w^2)))+1/2-1/4)





n*(1/2/pi*sum(atan((rho.w-sqrt(w.var/rev(w.var)))/sqrt(1-rho.w^2)/sqrt(1+  1e-10 )  ))+1/2-1/4)
n*(1/2/pi*sum(atan((rho.w-sqrt(w.var/rev(w.var)))/sqrt(1-rho.w^2)/sqrt(1+  1e-7  )  ))+1/2-1/4)

## rm(list=names(S.env))
## n*cov.1234
n*(1/2/pi*sum(atan( (rho.w-sqrt(w.var/rev(w.var))) / sqrt(1-rho.w^2)  ))+1/4)




source('misc.R')
n <- 5e1
s <- S.env(supp.S=c(1,5))$rS(n)
## V <- runif(2,1,4)
## cov.W <- runif(1,0,.001)
## stopifnot(cov.W<sqrt(prod(V)))
V <- sum(s^2)^2/sum(s[5:n]^2) * (2/(c(s[1]-s[2],s[3]-s[4]))^2 - 2/sum(s^2) + sum(s[1:4]^2)/sum(s^2)^2)
cov.W <- -sum(s^2)/sum(s[5:n]^2)*(2- sum(s[1:4]^2)/sum(s^2))
rho.W <- cov.W/sqrt(prod(V))
integrate(function(u)pnorm(u)*dminnorm(u,sd=sqrt(V),cov=cov.W),-Inf,Inf)$val
sum(sapply(1:2, function(j)1/sqrt(V[j])*integrate(function(u)pnorm(u)*dnorm(u/sqrt(V[j]))*pnorm(u/sqrt(1-rho.W^2)*(cov.W/V[j]/sqrt(V[3-j]) - 1/sqrt(V[3-j]))),-Inf,Inf)$val))
beta <- 1/sqrt(1-rho.W^2) * (rho.W/sqrt(V)-1/sqrt(rev(V))); gamma <- 1/sqrt(V)
sum(sapply(1:2, function(j)1/sqrt(V[j])*integrate(function(u)pnorm(u)*dnorm(u*gamma[j])*pnorm(u*beta[j]),-Inf,Inf)$val))
## sum(sapply(1:2, function(j)gamma[j]*(1/2/pi/gamma[j]*atan(beta[j]/gamma[j]/sqrt(1+beta[j]^2+gamma[j]^2))+1/4/gamma[j])  ))
1/2 + sum(1/2/pi*atan(beta/gamma/sqrt(1+beta^2+gamma^2)))

atan(beta/gamma/sqrt(1+beta^2+gamma^2))
atan((rho.W-sqrt(V/rev(V)))/sqrt(1-rho.W^2) / sqrt(1 + 1/(1-rho.W^2)*(rho.W^2/V-2*rho.W/sqrt(prod(V)) + 1/rev(V)) + 1/V) )

delta <- c(s[1]-s[2],s[3]-s[4])
beta/gamma/sqrt(1+beta^2+gamma^2)
-2/delta^2 * sum(s^2)^2/sum(s[5:n]^2)/ sqrt(prod(V)-cov.W^2 + 2*sum(1/delta^2)*sum(s^2)^2/sum(s[5:n]^2))
-2/delta^2 / sqrt(4/prod(delta^2)-4*sum(1/delta^2)*(1/sum(s^2)-sum(s[1:4]^2)/2/sum(s^2)^2 -1/2*sum(s[5:n]^2)/sum(s^2)^2))
-1/delta^2 / sqrt(1/prod(delta^2)-1*sum(1/delta^2)*1/2*(1/sum(s^2)))
    
beta/gamma
-2/delta^2/sqrt(prod(V)-cov.W^2) * sum(s^2)^2/sum(s[5:n]^2)
1+beta^2+gamma^2
1+ 2*sum(1/delta^2)/(prod(V)-cov.W)* sum(s^2)^2/sum(s[5:n]^2)

-2/delta^2 * sum(s^2)^2/sum(s[5:n]^2)/ sqrt(prod(V)-cov.W^2 + 2*sum(1/delta^2)*sum(s^2)^2/sum(s[5:n]^2))

dd

## 40d marginalize E(Cov.1234|s) over s
source('misc.R')
attach(S.env(supp.S=c(0,5)))
n <- 1e3
s <- matrix(rS(4*n),ncol=4)
w1.var <- 2*E.S2^2/(s[,1]-s[,2])^2-2*E.S2/n+rowSums(s^2)/n^2#sum(s[1:4]^2)/n^2
w2.var <- 2*E.S2^2/(s[,3]-s[,4])^2-2*E.S2/n+rowSums(s^2)/n^2#sum(s[1:4]^2)/n^2
w.var <- cbind(w1.var,w2.var)
cov.w <- -2*E.S2/n+rowSums(s^2)/n^2#sum(s[1:4]^2)/n^2
rho.w <- cov.w/sqrt(w1.var*w2.var)
alpha <- sqrt(n)/sqrt(1/n*sum(rS(n-4)^2)); gamma <- 1/sqrt(w.var); beta <- 1/sqrt(1-rho.w^2)*(rho.w/sqrt(w.var)-1/sqrt(w.var[,2:1]))
E.1234 <- mean(rowSums((1/2/pi*atan(alpha*beta*sqrt(w.var)/sqrt(alpha^2+beta^2+gamma^2))+1/4)))
n*(E.1234-1/4)
## E.1234 <- mean(1/2/pi*rowSums(atan((rho.w-sqrt(w.var/w.var[,2:1]))/sqrt(1-rho.w^2)))+1/2)
## n*(E.1234-1/4)
cov.1234 <- -mean.S.pair^2/4/pi/E.S2
cov.1234
k <- E.S2; delta <- cbind(abs(s[,1]-s[,2]),abs(s[,3]-s[,4]))
cov.1234 <- mean( 1/4/k/pi*rowSums((1+delta[,2:1]^2/delta^2)^(-1)*( -(delta[,1]*delta[,2]+delta[,2:1]^3/delta) + 1/2*delta[,2:1]/delta*rowSums(delta^2))) )
cov.1234
cov.1234 <- mean( 1/4/k/pi*rowSums((1+delta[,2:1]^2/delta^2)^(-1)*( -delta[,1]*delta[,2]-delta[,2:1]^3/delta + 1/2*delta[,2:1]/delta*rowSums(delta^2))) )
cov.1234


## 40e taylor approx [error]
source('misc.R')
attach(S.env(supp.S=c(0,5)))

n <- 50
pairs <- replicate(1e3, {
    s <- rS(n)
    w1.var <- 2*E.S2^2/(s[1]-s[2])^2-2*E.S2/n+sum(s[1:4]^2)/n^2
    w2.var <- 2*E.S2^2/(s[3]-s[4])^2-2*E.S2/n+sum(s[1:4]^2)/n^2
    w.var <- c(w1.var,w2.var)
    cov.w <- -2*E.S2/n+sum(s[1:4]^2)/n^2
    rho.w <- cov.w/sqrt(w1.var*w2.var)
    c(true=1/2/pi*sum(atan((rho.w-sqrt(w.var/rev(w.var)))/sqrt(1-rho.w^2)))+1/4,
      ## approx=rho.w/2/pi * sum(1/(1+(w.var/rev(w.var)))))
      taylor=cov.w/2/pi/sum(w.var)*(sqrt(w1.var/w2.var)+sqrt(w2.var/w1.var)),
      approx=-abs(s[1]-s[2])*abs(s[3]-s[4])/2/pi/n/E.S2)
})
op <- par(mfrow=c(1,2))
plot(pairs[1,],pairs[2,],main='taylor'); abline(a=0,b=1)
cor(t(pairs))[1,2]
plot(pairs['true',],pairs['approx',],main='approx'); abline(a=0,b=1)
cor(t(pairs))
par(op)

## 41 formula for correction
source('misc.R')
require(parallel)
supp.S <- c(1,2)
dS <- function(s)dunif(s,supp.S[1],supp.S[2])
pS <- function(q)punif(q,supp.S[1],supp.S[2])
rS <- function(n)runif(n,supp.S[1],supp.S[2])
E.S2 <- diff(supp.S)^2/12 + mean(supp.S)^2
mean.S.pair <- with(list(b=supp.S[2],a=supp.S[1]), 1/(b-a)*((a^2+b^2)/3-2/3*a*b))
n <- 3e2
ns <- round(seq(1e1,5e2,len=10))
## by.n <- mclapply(ns, mc.cores=detectCores()-2, FUN=function(n) {
by.n <- sapply(ns, FUN=function(n) {
    tau.stats <- mclapply(1:5e4, mc.cores=detectCores()-2, FUN=function(jj) {
        s <- rS(n)
        z <- rnorm(n)
        y <- z/s
        ## theta.fe <- sum(y*s^2)/sum(s^2)
        tau(z,s)
    })
    tau.stats <- simplify2array(tau.stats)
    var(tau.stats)
})
## by.n <- simplify2array(by.n)
## save.image('210130.RData')
ns*by.n - 4/9
cov.1234 <- -mean.S.pair^2/4/pi/E.S2/ns
4*(ns-2)*(ns-3)/(ns-1)*cov.1234
## png('210130.png')
plot(ns,ns*by.n-4/9,type='l',main='bias')
lines(ns,4*(ns-2)*(ns-3)/(ns-1)*cov.1234,lty=2)
legend('topright',lty=1:2,legend=c('observed','formula'))
## dev.off()
## plot(ns,ns*by.n-4/9 - 4*(ns-2)*(ns-3)/(ns-1)*cov.1234)

## 41a other S distributions
source('misc.R')
require(parallel)
rS <- rexp
rho <- 1/2
rSs <- c(uniform=function(n)runif(n,0,1), exponential=rexp, gamma=function(n)rgamma(n,shape=.54), beta=function(n)rbeta(n,.15,.39), pareto=function(n)rpareto(n,location=1,shape=2.52687))
rhos <- c(uniform=1/3,exponential=1/2,gamma=.56,beta=.66,pareto=.1390052)
rS <- rSs[['uniform']]; rho <- rhos['uniform']
n <- 3e2
ns <- round(seq(1e1,5e2,len=10))
## by.n <- mclapply(ns, mc.cores=detectCores()-2, FUN=function(n) {
by.n <- lapply(ns, FUN=function(n) {
    tau.stats <- mclapply(1:5e4, mc.cores=detectCores()-2, FUN=function(jj) {
        s <- rS(n)
        z <- rnorm(n)
        y <- z/s
        ## theta.fe <- sum(y*s^2)/sum(s^2)
        tau(z,s)
    })
    tau.stats <- simplify2array(tau.stats)
    var(tau.stats)
})
by.n <- simplify2array(by.n)
## save.image('210130.RData')
ns*by.n - 4/9
cov.1234 <- -rho/4/pi/ns
4*(ns-2)*(ns-3)/(ns-1)*cov.1234
## png('210130.png')
plot(ns,ns*by.n-4/9,type='l',main='bias')
lines(ns,4*(ns-2)*(ns-3)/(ns-1)*cov.1234,lty=2)
legend('topright',lty=1:2,legend=c('observed','formula'))
## dev.off()
## plot(ns,ns*by.n-4/9 - 4*(ns-2)*(ns-3)/(ns-1)*cov.1234)


## 41b type 1 error rate
source('misc.R')
require(parallel)
B <- 3e3
q <- qnorm(1-.05/2)
rS <- rexp
rho <- 1/2
rSs <- c(uniform=function(n)runif(n,0,1), exponential=rexp, gamma=function(n)rgamma(n,shape=.54), beta=function(n)rbeta(n,.15,.39), pareto=function(n)rpareto(n,location=1,shape=2.52687))
rhos <- c(uniform=1/3,exponential=1/2,gamma=.56,beta=.66,pareto=.1390052)
biases <- rhos/pi
distributions <- structure(names(rSs),names=names(rSs))
## distributions <- distributions[1:2]
ns <- round(seq(1e1,1.5e2,len=20))
by.distr <- lapply(distributions, function(distr) {
    rS <- rSs[[distr]]; bias <- biases[distr]
    ## rS <- runif; bias <- 1/3/pi
    by.n <- sapply(ns, FUN=function(n) {
        tau.stats <- mclapply(1:B, mc.cores=detectCores()-2, FUN=function(jj) {
            s <- rS(n)
            z <- rnorm(n)
            s.md.hat <- mean(abs(outer(s,s,`-`)))
            s2.hat <- mean(s^2)
            bias.hat <- s.md.hat^2/s2.hat/pi
            c(tau=tau(z,s),bias.hat=bias.hat)
        })
        tau.stats <- simplify2array(tau.stats)
        tau.hat <- tau.stats['tau',]; bias.hat <- tau.stats['bias.hat',]
        ## var(tau.stats)
        c(power=mean(sqrt(9*n/4)*abs(tau.hat)>q),power.debiased=mean(sqrt(n/(4/9-bias))*abs(tau.hat)>q),power.debiased.hat=mean(sqrt(n/(4/9-bias.hat))*abs(tau.hat)>q))
    })
})
## plot(ns,by.n['power',])
## lines(ns,by.n['power.debiased',])
## abline(h=0.05)

## load('210313b.RData')
power <- sapply(by.distr,function(m)m['power',])
power.debiased <- sapply(by.distr,function(m)m['power.debiased',])
power.debiased.hat <- sapply(by.distr,function(m)m['power.debiased.hat',])
matplot(power,x=ns,type='l',lty=1:length(distributions),col=1,xlab='number of studies',ylab='FPR',ylim=range(unlist(by.distr)))
matplot(ns,power.debiased,lty=1:length(distributions),type='l',col=2,add=TRUE)
matplot(ns,power.debiased.hat,lty=1:length(distributions),type='l',col=3,add=TRUE)
abline(h=.05)
legend('topright',lty=1:length(distributions),legend=distributions)
## save.image('210313b.RData')

require(xtable)
ns.idx <- which(ns %in% c(25,76,150))
out <- matrix(paste0(t(round(power,2)[ns.idx,]),', ',t(round(power.debiased,2)[ns.idx,]),', ',t(round(power.debiased.hat,2)[ns.idx,])),nrow=length(distributions))
attr(out,'dimnames') <- list(distribution=distributions,n=ns[ns.idx])
addtorow <- list(pos=list(0, 0),
                 command=c("& \\multicolumn{3}{c}{meta-analysis size} \\\\\n",
                      paste0('precision distribution &',paste0(ns[ns.idx],collapse='&'), '\\\\\n')))
sink('ms/210313_table.tex')
print.xtable(xtable(out),add.to.row = addtorow,floating=FALSE,latex.environment=NULL,include.colnames = FALSE)
## print.xtable(xtable(out),add.to.row = addtorow,include.colnames = FALSE)
sink()

dd

## ## power simulation [error]
## require(parallel)
## source('misc.R')
## alpha <- .01
## q <- qnorm(1-alpha/2)
## n <- 5e2
## a <- .0
## b <- 1
## as <- seq(0,3-.1,length.out=4)
## ## bs <- seq(.1,5,length.out=20)
## ## as <- seq(0,1-.1,length.out=20)
## ## bs <- seq(.5,5,length.out=3)
## by.endpoint <- sapply(as, function(a) {
##     ## attach(S.env(c(b-b[1],b)))
##     attach(S.env(c(a,a+.1)))
##     rejected <- mclapply(1:1e4, mc.cores=detectCores()-2, FUN=function(jj) {
##         s <- rS(n)
##         z <- rnorm(n,1)
##         y <- z/s
##         tau.stat <- tau(z,s)
##         sqrt(9/4*n)*abs(tau.stat) > q
##     })
##     rejected <- simplify2array(rejected)
##     detach(S.env(c(a,a+.1)))
##     mean(rejected)    
## })
## plot(as,by.endpoint)


## 42
n <- 1e4
## s <- matrix(runif(n*2),ncol=2)
s <- matrix(rbinom(n*2,1,1/10),ncol=2)
hist(abs(diff(s)),prob=TRUE)
abline(v=mean(abs(diff(s))),col='red')
abline(v=sqrt(mean(abs(diff(s)^2))),col='blue')
mean(abs(diff(s))) / sqrt(mean(abs(diff(s)^2)))

ps <- seq(0.01,.99,len=20)
by.p <- sapply(ps, function(p){
    ## s <- matrix(rbinom(n*2,1,p),ncol=2)
    ## s <- matrix(rnorm(n*2,0,p+1),ncol=2)
    ## s <- matrix(rexp(n*2,1/p),ncol=2)
    ## s <- matrix(rbeta(n*2,5-p*3,2+p*3),ncol=2)    
    ## s <- matrix(rbeta(n*2,7-p*6,1+p*6),ncol=2)    
    ## s <- matrix(rbeta(n*2,5*p,5*p),ncol=2)    
    s <- matrix(rbeta(n*2,20*p,20),ncol=2)    
    c(pair=mean(abs(diff(s))), sd=sqrt(mean(abs(diff(s)^2))))
    ## c(pair=mean(abs(diff(s)))^2, sd=mean(abs(diff(s)^2)))
})
plot(ps,by.p['pair',],ylim=range(c(by.p,by.p['pair',]/by.p['sd',])))
points(ps,by.p['sd',],col='red')
lines(ps,by.p['pair',]/by.p['sd',],col='blue')
## pair=var for bernoulli:
## curve(2*x*(1-x),add=TRUE)
## curve(sqrt(2*x*(1-x)),add=TRUE)

plot(0,type='n',xlim=c(0,1),ylim=c(0,8))
## for(p in ps) curve(dbeta(x,7-p*6,1+p*6),add=TRUE)
for(p in ps) curve(dbeta(x,5*p,5*p),add=TRUE)


## 43 beta maximization

ps <- seq(0.01,.99,len=20)
beta1 <- function(p)3*p; beta2 <- function(p)3
## beta1 <- function(p)3*p; beta2 <- function(p)3
## beta1 <- function(p)5-p*3; beta2 <- function(p)2+p*3
## beta1 <- function(p)7-p*6; beta2 <- function(p)1+p*6
## beta1 <- beta2 <- function(p)5*p
by.p <- sapply(ps, function(p){
    s <- matrix(rbeta(n*2,beta1(p),beta2(p)),ncol=2)    
    c(pair=mean(abs(diff(s)))^2, s2=mean(s^2))
})
op <- par(mfrow=c(1,2))
plot(ps,by.p['pair',],ylim=range(c(by.p,by.p['pair',]/by.p['s2',])))
points(ps,by.p['s2',],col='red')
lines(ps,by.p['pair',]/by.p['s2',],col='blue')
plot(0,type='n',xlim=c(0,1),ylim=c(0,8))
for(p in ps) curve(dbeta(x,beta1(p),beta2(p)),add=TRUE)
p.star <- ps[which.max(by.p['pair',]/by.p['s2',])]
curve(dbeta(x,beta1(p.star),beta2(p.star)),add=TRUE,col='blue')
curve(dbeta(x,beta1(ps[1]),beta2(ps[1])),add=TRUE,col='green',lty=2)
curve(dbeta(x,beta1(rev(ps)[1]),beta2(rev(ps)[1])),add=TRUE,col='red',lty=2)
par(op)

max <- 0
argmax <- NULL
beta.md <- function(betas)4/sum(betas)*beta(sum(betas),sum(betas))/beta(betas[1],betas[1])/beta(betas[2],betas[2])
## with(list(betas=abs(runif(2))), mean(abs(rbeta(1e5,betas[1],betas[2])-rbeta(1e5,betas[1],betas[2])) - beta.md(betas)))
beta.s2 <- function(betas)(betas[1]/sum(betas))^2 + prod(betas)/sum(betas)^2/(sum(betas)+1)
beta.s2 <- function(betas)betas[1]/sum(betas)^2*(betas[1]^2+prod(betas)+sum(betas))/(sum(betas)+1)
## with(list(betas=abs(runif(2))), mean(rbeta(1e5,betas[1],betas[2])^2) - beta.s2(betas))
ratio <- function(betas)prod(betas>0)*beta.md(betas)^2/beta.s2(betas)
invisible(replicate(1e6,{
    betas <- runif(2)*1
    try <- ratio(betas)#beta.md(betas)^2/beta.s2(betas)
    if(try>max){max <<- try;argmax <<- betas}
}))
print(max);print(argmax)
## png('./ms/fig1.png')
curve(dbeta(x,argmax[1],argmax[2]),xlab='',ylab='')
## dev.off()
## beta(argmax[1],argmax[2])


argmax <- optim(c(.5,.5),ratio,control=list(fnscale=-1))$par



alpha <- argmax[1]
f <- Vectorize(function(beta)beta.md(c(alpha,beta))^2/beta.s2(c(alpha,beta)))

f <- function(beta) 4*(alpha+beta+1)*(alpha*(alpha+beta+1)+beta)*(digamma(alpha+beta)-digamma(2*(alpha+beta))-digamma(beta)+digamma(2*beta)) - (alpha+1)
uniroot(f,c(0.1,.9),extendInt='y')

beta <- argmax[2]
f <- function(alpha) 4*(alpha+beta+1)*(alpha^3+alpha^2*beta+alpha^2+alpha*beta)*(digamma(alpha+beta)-digamma(2*(alpha+beta))-digamma(alpha)+digamma(2*alpha)) - (alpha+beta+1)*(2*alpha^2+2*alpha*beta+2*alpha+beta) + alpha*beta
uniroot(f,c(01,.5),extendInt='y')


f <- function(beta) digamma(beta)-digamma(2*beta)-digamma(alpha)+digamma(2*alpha) - (2*alpha+1)/4/alpha/(alpha+1)
uniroot(f,c(.1,.5),extendInt='n')

series <- Vectorize(function(x,n=1e3)sum(x/(0:n+x)/(0:n+2*x)))
series(alpha,n=1e6)-series(beta,n=1e6) -  1/4*(1/alpha + 1/(alpha+1))
series(beta^2,n=1e6)-series(beta,n=1e6) -  1/4*(1/beta^2 + 1/(beta^2+1))
curve(series(x)-series(beta),.04,1)
curve(series(x)-1/4*(1/alpha+1/(alpha+1)),add=TRUE,lty=2)
abline(h=series(beta))

curve(digamma(x)-digamma(2*x)-(2*x+1)/(4*x*(x+1)),.01,40)

require(rgl)
ratio2 <- ratio <- function(beta1,beta2)(beta1>0)*(beta2>0)*beta.md(c(beta1,beta2))^2/beta.s2(c(beta1,beta2))
ratio2 <- Vectorize(ratio2)
alpha <- beta <- seq(0,1,len=30)
z <- outer(alpha,beta,function(x,y)ratio2(x,y))
persp3d(x,y,z)



## 44 gamma maximization
k <- 6
(gamma(k+1)/gamma(k))^2*integrate(function(s)1/s*(1-pgamma(s,shape=k+1))*dgamma(s,shape=k+1),0,Inf)$val
integrate(Vectorize(function(s2)integrate(function(s1)s1*dgamma(s1,shape=k),s2,Inf)$val*dgamma(s2,shape=k)),0,Inf)$val - integrate(Vectorize(function(s2)integrate(function(s1)s1*dgamma(s1,shape=k),0,s2)$val*dgamma(s2,shape=k)),0,Inf)$val
mean(abs(apply(matrix(rgamma(2*1e4,shape=k),ncol=2),1,diff)))
gamma.md <- function(k)2*k*(1 - 2*integrate(function(s2)(pgamma(s2,shape=k+1))*dgamma(s2,shape=k),0,Inf)$val)
2*(k - 2*k^2*integrate(function(s2)(1/s2*pgamma(s2,shape=k+1))*dgamma(s2,shape=k+1),0,Inf)$val)
gamma.s2 <- function(k) k + k^2
ratio <- Vectorize(function(k) gamma.md(k)^2/gamma.s2(k))
curve(ratio,0.1,10)
## curve(ratio,0.558,.5583)

curve(dgamma(1,shape=x)*(log(1)-digamma(x)),.01,10)

ff <- Vectorize(function(k)gamma.md(k)^2)
curve(ff(x)/x/(x+1),0,10)
curve(ff(x),0,.5)
ff <- Vectorize(function(k)(1-2*integrate(function(s)pgamma(s,shape=k+1)*dgamma(s,shape=k),0,Inf)$val)^2)
curve(ff(x),0,20)
curve(x/(x+1),0,10)

argmax <- optimize(ratio,lower=.01,upper=1,maximum=TRUE)$obj
max.bias <- optimize(ratio,lower=.01,upper=1,maximum=TRUE)$max / pi
curve(dgamma(x,shape=argmax))

k*(4*pbeta(1/2,k+1,k)-2)

curve(log(pbeta(1/2,x,x+1)-1/2),0,2)
abline(a=log(pbeta(1/2,k,k+1)-1/2)-slope*k,b=slope,col='red')


delta <- .01;grid <- seq(0,20,by=delta)
plot(grid[2:length(grid)],1/delta*diff(log(pbeta(1/2,grid,grid+1)-1/2)),type='l',ylim=c(-3,0))
curve(-1/2/x/(x+1),add=TRUE)

curve(-1/2*log(x/(x+1)),0,4,ylim=c(-2,2))
abline(v=argmax)
curve(1.7+log(pbeta(1/2,x,x+1)-1/2),add=T,lty=2)

op <- par(mfrow=c(1,3))
curve(log(sqrt(x/(x+1))),0,10)
curve(log(pbeta(1/2,x,x+1)-1/2),0,10)
curve(log(sqrt(x/(x+1))) + log(pbeta(1/2,x,x+1)-1/2),0,10)
abline(v=argmax)
par(op)

CF.approx <- function(cf)Reduce(function(x,y)y/(x+1),c(rev(cf),1))
n <- 30
n <- 2*n+1
cf <- numeric(n)
k <- .2
j <- 0:((n-1)/2)
cf[2*j+1] <- -(k+j)*(2*k+1+j)/2 / ((k+2*j)*(k+2*j+1))
j <- 1:(n/2)
cf[2*j] <- j*(k+1-j)/2 / ((k+2*j-1)*(k+2*j))
pbeta(1/2,k,k+1)
1/k * 1/2^(2*k+1) * CF.approx(cf) / beta(k,k+1)


curve(sqrt(x/(x+1))*(pbeta(1/2,x,x+1)-1/2),ylim=c(-1,1)); abline(v=argmax)
curve((pbeta(1/2,x,x+1)-1/2),add=TRUE,col='red')
curve(sqrt(x/(x+1)),add=TRUE,col='blue')

ibeta.cf <- function(k,acc=5) {
    n <- 2*acc+1
    cf <- numeric(n)
    j <- 0:((n-1)/2)
    cf[2*j+1] <- -(k+j)*(2*k+1+j)/2 / ((k+2*j)*(k+2*j+1))
    j <- 1:(n/2)
    cf[2*j] <- j*(k+1-j)/2 / ((k+2*j-1)*(k+2*j))
    CF.approx(cf) 
}
k <- 6
pbeta(1/2,k,k+1)
1/k * 1/2^(2*k+1) * ibeta.cf(k,acc=3) / beta(k,k+1)

k <- pi
1/beta(k,k+1)
2^(2*k-1)*2*k/pi*integrate(function(theta)sin(theta)^(2*k),0,pi)$val

poch <- Vectorize(function(x,n) (n>0)*prod(x:(x+n-1))+(n==0))

k <- .4
pbeta(1/2,k,k+1)*beta(k,k+1)
ns <- 0:30
-1/2^k*sum(choose(ns-k-1,ns)/(k+ns)/2^ns)
1/2^k*sum(poch(-k,ns)/(factorial(ns))/(k+ns)/2^ns)


k <- .4
ns <- 0:5
pbeta(1/2,k,k+1)
k*2^k/pi*integrate(function(theta)sin(theta)^(2*k),0,pi)$val * sum(choose(ns-k-1,ns)/(k+ns)/2^ns)
p1 <- Vectorize(function(k)k*2^k/pi*integrate(function(theta)sin(theta)^(2*k),0,pi)$val)
p2 <- Vectorize(function(k)  sum(choose(ns-k-1,ns)/(k+ns)/2^ns))
p3 <- function(k)sqrt(k/(k+1))
p1(k)*p2(k)
curve(p1,0,5,ylim=c(0,2))
curve(p2,add=TRUE)
curve(p3,add=TRUE)
curve(p1(x)*p2(x)-1/2,0,3)

p1 <- Vectorize(function(k)1/2/pi*integrate(function(theta)sin(theta)^(2*k),0,pi)$val)
p2 <- Vectorize(ibeta.cf)
p3 <- function(k)sqrt(k/(k+1))
curve(p3(x)*(p1(x)*p2(x)-1/2))


curve(p2(x,acc=20))
curve(p2(x,acc=1),add=TRUE)
curve((x^2+7/2*x+2)/(x+1),add=TRUE,col='red')
curve((7*x^2+19*x+10)/(x^2+6*x+5),add=TRUE,col='blue')


## 44a pareto maximization

require(EnvStats)
n <- 1e5
shape <- runif(1,2,10)
mean(abs(rpareto(n,location=1,shape=shape) - rpareto(n,location=1,shape=shape)))
2*shape/(shape-1)/(2*shape-1)
mean(rpareto(n,location=1,shape=shape)^2)
shape/(shape-1)^2/(shape-2) + (shape/(shape-1))^2

curve( (2*shape/(shape-1)/(2*shape-1))^2 / (shape/(shape-1)^2/(shape-2) + (shape/(shape-1))^2), xname='shape',2,4)

k <- shape
(2*shape/(shape-1)/(2*shape-1))^2 / (shape/(shape-1)^2/(shape-2) + (shape/(shape-1))^2)
4*k^2/(2*k-1)^2*(k-2)/(k^3-2*k^2+k)

curve( (2*k/(k-1)/(2*k-1))^2 / (k/(k-1)^2/(k-2) + (k/(k-1))^2), xname='k',2,4)

argmax <- uniroot(function(k)-2*k^3+6*k^2-2*k-1,c(2,3))$root
curve(-2*k^3+6*k^2-2*k-1,xname='k',2,3)
1+2*sqrt(2/3)*cos(1/3*atan(sqrt(101/27)))
rho.max <- 4*argmax*(argmax-2)/((2*argmax-1)*(argmax-1))^2


n <- 1e4; mean(n/sum(rchisq(n,df=1)^2))

n <- 1e4; mean(n/sum(runif(n)^2))

n <- 1e5; mean(n/sum(rcauchy(n)^2))

source('misc.R')
a <- runif(1,0,10); b <- a+runif(1,0,3)
n <- 1e4
with(S.env(c(a,b)), {
    print(c(mean(rS(n)^2),E.S2))
    print(c(mean.S.pair, mean(abs(rS(n)-rS(n)))))
    ## print(c( mean(abs(rS(n)-rS(n)))^2 / mean(rS(n)^2), 1/9/(a^2+a+1/3)  ))
    print(c( mean(abs(rS(n)-rS(n)))^2 / mean(rS(n)^2), (b-a)^2/(3*(b^2+a^2+a*b))  ))
    })

curve((b-a)^2/3/(a^2+b^2+a*b),xname='a',-100,100)
abline(v=-b)



maxfound <- a.max <- b.max <- 0
for(jj in 1:1e5) {
    a <- runif(1,0,20); b <- a+runif(1,0,20)
    try <- (b-a)^2/(3*(b^2+a^2+a*b))
    if(try>maxfound){a.max <- a;b.max <- b;maxfound <- try}
}
print(c(a.max,b.max))
print(maxfound)



## 45 cov term sign
n <- 10
## cors <- replicate(1e3,{
pairs <- replicate(1e5, {
    ## z <- runif(n)
    z <- rnorm(n)
    s <- runif(n,1,2)
    ## c((z[1]-z[2])/(s[1]-s[2]) - mean(z[1:4]*s[1:4]),(z[3]-z[4])/(s[3]-s[4]) - mean(z[1:4]*s[1:4]))
    ## c((z[1])/(s[1]-s[2]),(z[1]*s[1]))
    c(1/(s[1]-s[2]),s[1])
})
## plot(t(pairs))
cov(t(pairs))[1,2]
## })
## mean(cors)+c(-1,1)*3*sd(cors)/sqrt(length(cors))
## hist(cors)
## abline(v=mean(cors),col='red')
s1 <- runif(1e6,1,4);s2 <- runif(1e6,1,4)
## -4/n*mean(s1/(s1-s2)) + 4/n^2*mean(s1^2)
cov(1/(s1-s2),s1)



require(mvtnorm)
rho <- -.5
plot(rmvnorm(1e5,sigma=matrix(c(1,rho,rho,1),nrow=2)),col=rgb(0,0,1,.7),pch='.')
abline(a=0,b=1);abline(a=0,b=-1);abline(h=0);abline(v=0)

n <- 1e5
rho <- .10
rdist <- runif
rdist <- rnorm
rdist <- function(n)rbeta(n,.1,.1)
W <- rdist(n)
Z <- rdist(n)
X <- rho*Z-W
Y <- rho*W-Z
qs <- seq(-2,2,len=40)
by.q <- sapply(qs, function(q) {
    c(ind=mean(c(X,Y)>q),
    dep=mean(X>q & Y>q)/mean(Y>q))
})
plot(qs,by.q['ind',],type='l',ylim=c(0,1))
lines(qs,by.q['dep',],lty=2)
abline(v=0)
## plot(qs,by.q['dep',]-by.q['ind',],type='l')


W <- rdist(n)
Z <- rdist(n)
X <- rho*Z-W
Y <- rho*W-Z
cov(X>.05,Y>.05)



require(mvtnorm)
## rho <- -1
## xy <- rmvnorm(1e5,sigma=matrix(c(1,rho,rho,1),nrow=2))
n <- 1e5
m <- 100
supp.S <- c(0,3)
E.S2 <- diff(supp.S)^2/12 + mean(supp.S)^2
z <- matrix(rnorm(n*4),ncol=4)
s <- matrix(runif(n*4,supp.S[1],supp.S[2]),ncol=4)
xy <- cbind(  (z[,1]-z[,2])/(s[,1]-s[,2]),(z[,3]-z[,4])/(s[,3]-s[,4]) ) - rowSums(z*s)
th <- (rowSums(z*s)+rnorm(m))/(rowSums(s^2)+runif(m,supp.S[1],supp.S[2]))#+rnorm(nrow(xy),sd=.10)
xy <- cbind((z[,1]-z[,2]-th*(s[,1]-s[,2]))*(s[,1]-s[,2])>0,(z[,3]-z[,4]-th*(s[,3]-s[,4]))*(s[,3]-s[,4])>0)
cov(xy)

dd

sd <- 1/100
a <- mean(pnorm(apply(xy,1,min),sd=sd))
b <- mean(1-pnorm(apply(xy,1,max),sd=sd))
c <- mean((pnorm(xy[,1],sd=sd)-pnorm(xy[,2],sd=sd))*(xy[,2]>xy[,1]))
d <- mean((pnorm(xy[,2],sd=sd)-pnorm(xy[,1],sd=sd))*(xy[,2]>xy[,1]))
c(a,b,c,d)
a-(a+c)*(a+d)

mean(xy>th)^2
mean(apply(xy,1,min)>th)


n <- 10
pairs <- replicate(1e5, {
    z <- rnorm(n)
    s <- runif(n,0,5)
    th <- sum(z*s)/sum(s^2)
    th <- sum((z*s)[1:4])/sum(s[1:4]^2)
    c((z[1]-z[2]-th*(s[1]-s[2]))*(s[1]-s[2])>0, (z[3]-z[4]-th*(s[3]-s[4]))*(s[3]-s[4])>0)
    ## c((z[1]-z[2])/(s[1]-s[2])>th, (z[3]-z[4])/(s[3]-s[4])>th)
})
n*cov(t(pairs))[1,2]



## xy <- rmvnorm(1e5,sigma=matrix(c(1,rho,rho,1),nrow=2))
n <- 1e6
m <- 100
rS <- function(n)runif(n,0,.5)
## rS <- rexp
rZ <- function(n)rnorm(n)
## rZ <- function(n)runif(n,-1,1)
## rZ <- function(n)rexp(n)-1
## rZ <- function(n)sample(c(-10,10/9),size=n,replace=TRUE,prob=c(1/10,9/10))
z <- matrix(rZ(n*4),ncol=4)
s <- matrix(rS(n*4),ncol=4)
th <- with(list(s.th=matrix(rS(n*m),ncol=m),z.th=matrix(rZ(n*m),ncol=m)),
           (1*rowSums(z*s)+rowSums(z.th*s.th))/(rowSums(s^2)+rowSums(s.th^2)))
## th <- 0
## xy <- cbind(  (z[,1]-z[,2])/(s[,1]-s[,2]) > th  ,(z[,3]-z[,4])/(s[,3]-s[,4]) < th)
x <- (z[,1]-z[,2])/(s[,1]-s[,2]) 
y <- (z[,3]-z[,4])/(s[,3]-s[,4]) 
2*mean((x>th)*(y<th))
## m*cov(x>th,y>th)
## th <- rnorm(n)
mean(x*y-(x+y)<0)
mean(x*y-(x+y)+th^2<0)
mean(x*y-th*(x+y)+th^2<0)

a <- mean((x>th)*(y<th))
b <- mean((x<th)*(y>th))
c <- mean((x>th)*(y>th))
d <- mean((x<th)*(y<th))
c - (a+c)*(b+c)
c - (a+c)^2

n <- 1e5
rZ <- function(n)runif(n,-1,1)
## rZ <- rnorm
x <- rZ(n); y <- rZ(n)
## target <- x*y-(x-(2*rbinom(n,1,.5)-1)*y)
target <- x*y+(x-y)*(2*rbinom(n,1,.5)-1)
summary(target)
hist(target)
mean(target>0)


mu <- runif(2,1,10)
n <- 1e4
z <- matrix(rnorm(2*n),ncol=2)
a <- (z[,1]+sqrt(sum(mu^2)))^2 + z[,2]^2
b <- a + 2*rowSums(t(mu*t(z))) - 2*z[,1]*sqrt(sum(mu^2))
qqplot(a,b);abline(a=0,b=1)

require(mvtnorm)
n <- 5e3
rho <- .9
qZ <- qnorm
## qZ <- qcauchy
qZ <- function(p)qunif(p,-.01,.01)
pZ <- pnorm
## pZ <- function(q)punif(q,-.01,.01)
xy <- rmvnorm(n,sigma=matrix(c(1,rho,rho,1),nrow=2))
u1 <- pnorm(xy[,1]); u2 <- pnorm(xy[,2])
x <- qZ(u1); y <- qZ(u2)
z <- rnorm(n)
mean(x*y - z*(x+y) > 0)
## cor(xy[,1]-z > 0, xy[,2]-z > 0)
## c(mean((z < x*y/(x+y))*(x+y>0)*(x*y>0)),mean((z < x*y/(x+y))*(x+y>0)*(x*y<0)),
##   mean((z > x*y/(x+y))*(x+y<0)*(x*y>0)),mean((z > x*y/(x+y))*(x+y<0)*(x*y<0)))
## c(mean(pnorm(x*y/(x+y))*(x+y>0)*(x*y>0)),mean(pnorm(x*y/(x+y))*(x+y>0)*(x*y<0)),
##   mean((1-pnorm(x*y/(x+y)))*(x+y<0)*(x*y>0)),mean((1-pnorm(x*y/(x+y)))*(x+y<0)*(x*y<0)))
mean(pZ(x*y/(x+y))[(x+y>0)&(x*y<0)])*mean(x*y<0) + mean(pZ(x*y/(x+y))[(x+y>0)&(x*y>0)])*mean(x*y>0)
c(mean(pZ(x*y/(x+y))[(x+y>0)&(x*y<0)]),mean(pZ(x*y/(x+y))[(x+y>0)&(x*y>0)]))-1/2
c(mean(x*y<0),mean(x*y>0))
## mean((1-pnorm(x*y/(x+y)))[(x+y<0)&(x*y<0)])
## summary(((x*y)/(x+y))[x+y>0])

## exponential
n <- 1e6
lambda <- runif(1,0,10)
lambda
mean(abs(rexp(n,1/lambda)-rexp(n,1/lambda))) 

n <- 1e5
p <- 1
pZ <- pcauchy
pZ <- function(q)punif(q,-1,1)
## pZ <- function(q)pnorm(q,mean=0)
s <- matrix(runif(4*n,0,5),ncol=4)
z <- matrix(rnorm(4*n,sd=2),ncol=4)
x <- (z[,1]-z[,2])/(s[,1]-s[,2]) - rowSums(z*s)/rowSums(s^2)/p
y <- (z[,3]-z[,4])/(s[,3]-s[,4]) - rowSums(z*s)/rowSums(s^2)/p
c(mean(pZ((x*y/(x+y))[x+y>0])),mean(1-pZ((x*y/(x+y))[x+y<0])))
mean(c(mean(pZ((x*y/(x+y))[x+y>0])),mean(1-pZ((x*y/(x+y))[x+y<0]))))
c(mean(pZ((x*y/(x+y))[(x+y>0) & (x*y<0)])),mean(pZ((x*y/(x+y))[(x+y>0) & (x*y>0)])))
sum(c(mean(pZ((x*y/(x+y))[(x+y>0) & (x*y<0)])),mean(pZ((x*y/(x+y))[(x+y>0) & (x*y>0)]))))
mean(pZ(x*y/(x+y))*(x+y>0))
mean(pZ(x*y/(x+y))*(x+y<0))

mean(x*y/(x+y)>0)
c(mean((x+y>0) & (x*y<0)),mean((x+y>0) & (x*y>0)))

c(mean(pZ((x*y/(x+y))[(x+y<0) & (x*y<0)])),mean(pZ((x*y/(x+y))[(x+y<0) & (x*y>0)])))
## 2*sum(c(mean((x+y>0) & (x*y<0)),mean((x+y>0) & (x*y>0))) *
##  c(mean(pZ((x*y/(x+y))[x+y>0 & x*y<0])),mean(pZ((x*y/(x+y))[x+y>0 & x*y>0]))))
## mean(c(mean(pZ((x*y/(x+y))[x+y>0 & x*y<0])),mean(pZ((x*y/(x+y))[x+y>0 & x*y>0]))))

mean(x+y>0)
sum((x+y>0) - (((1/x+1/y<0) & (x*y<0)) | ((1/x+1/y>0) & (-y<0) & (x>0))))
mean(pZ((x*y/(x+y))[x+y>0]))
mean(pZ((x*y/(x+y))[x+y>0 & x*y>0]))*mean(x*y>0 & x+y>0)/mean(x+y>0)+mean(pZ((x*y/(x+y))[x+y>0 & x*y<0]))*mean(x*y<0 & x+y>0)/mean(x+y>0)

require(mvtnorm)
n <- 5e3
x <- rnorm(n);y <- rnorm(n);z <- rnorm(n)
mean(z/100<pmin(x,y))
p <- 2
rho <- -.5
xy <- rmvnorm(n,sigma=matrix(c(1,rho,rho,1),nrow=2))
x <- xy[,1];y <- xy[,2]
## mean(z/100<pmin(x,y))
## mean(pnorm(x*y/(x+y)*50)*(x+y>0))
## mean(pnorm(x*y/(x+y)*p)*(x+y>0)*(x*y>0))
## mean(pnorm(x*y/(x+y)*p)*(x+y<0)*(x*y>0))
mean(pnorm(pmin(x,y),sd=.1))

dd





## cov term sign
n <- 1e5; p <- 10
x <- rnorm(n);y <- rnorm(n);z <- 0*rnorm(n)/sqrt(p)+(x+y)/p
mean((z<x)&(z<y)) + mean((z>x)&(z>y))
c(mean((z<x)&(z<y)),mean((z>x)&(z>y)),mean((z<x)&(z>y)), mean((z>x)&(z<y)))
cov((x>z),(y>z))


n <- 1e5
z <- matrix(rnorm(3*n),ncol=3)
s <- matrix(runif(3*n),ncol=3)
mean(pmin((z[,1]-z[,2])/(s[,1]-s[,2]),(z[,1]-z[,3])/(s[,1]-s[,3]))>0)

a <- rnorm(2);b <- rnorm(2)

n <- 1e6
p <- 5e2
z <- matrix(rnorm(4*n),ncol=4)
z <- matrix(runif(4*n,-5,5),ncol=4)
s <- matrix(runif(4*n,0,5),ncol=4)
## mean(pmin((z[,1]-z[,2])/(s[,1]-s[,2]),(z[,3]-z[,4])/(s[,3]-s[,4])) - (rowSums(z*s)+rnorm(n)*(p-4))/(rowSums(s^2)+runif(n)*(p-4))>0)
## mean(pmin((z[,1]-z[,2])/(s[,1]-s[,2]),(z[,3]-z[,4])/(s[,3]-s[,4]))>1/p)
## mean(pmin((z[,1]-z[,2])/(s[,1]-s[,2]),(z[,3]-z[,4])/(s[,3]-s[,4]))>0)
## mean((z[,1]-z[,2])/(s[,1]-s[,2])*(z[,3]-z[,4])/(s[,3]-s[,4])>0)
## mean(((s[,1]-1/(s[,1]-s[,2]))*s[,1]) + ((s[,2]+1/(s[,1]-s[,2]))*s[,2]))
## summary((s[,1]-1/(s[,1]-s[,2]))*s[,1])
## summary((s[,2]+1/(s[,1]-s[,2]))*s[,2])
## mean((s[,1]-1/(s[,1]-s[,2]))*s[,1] + (s[,2]+1/(s[,1]-s[,2]))*s[,2] >0)
## mean(z[,1]^2*(s[,1]-1/(s[,1]-s[,2]))*s[,1] + z[,2]^2*(s[,2]+1/(s[,1]-s[,2]))*s[,2] >0)
## mean((z[,1]-z[,2])/(s[,1]-s[,2])*(z[,3]-z[,4])/(s[,3]-s[,4]) > (rowSums(z*s))+runif(n)/p) - 1/2
mean((z[,1]-z[,2] - 1/p*(a[1]*z[,3]+a[2]*z[,4]) > 0) & (z[,3]-z[,4] - 1/p*(b[1]*z[,1]+b[2]*z[,2])>0))

mean(pmin((z[,1]-z[,2])/(s[,1]-s[,2]),(z[,3]-z[,4])/(s[,3]-s[,4]))>rowSums(z*s))
mean(rowSums(cbind(s[,1]-1/(s[,1]-s[,2]),s[,2]+1/(s[,1]-s[,2]),s[,3],s[,4]) * z)*
rowSums(cbind(s[,1],s[,2],s[,3]-1/(s[,3]-s[,4]),s[,4]+1/(s[,3]-s[,4])) * z)>0)/2


## counterex to P(z1-z2-(a1z3+a2z4)/p>0,z3-z4-(b1z1+b2z2)/p>0)<1/4
a <- c( 0.8616181, 0.6051115)
b <- c(-0.6209407 , 1.4389778)
ps <- seq(10,200,len=10)
by.p <- sapply(ps, function(p) {
    z <- matrix(runif(4*n,-5,5),ncol=4)
    s <- matrix(runif(4*n,0,5),ncol=4)
    mean((z[,1]-z[,2] - 1/p*(a[1]*z[,3]+a[2]*z[,4]) > 0) & (z[,3]-z[,4] - 1/p*(b[1]*z[,1]+b[2]*z[,2])>0))
})
plot(ps,by.p)

s <- runif(4,1,5)
v1 <- matrix(c(s[1]-1/(s[1]-s[2]),s[2]+1/(s[1]-s[2]),s[3],s[4]),ncol=1)
v2 <- matrix(c(s[1],s[2],s[3]-1/(s[3]-s[4]),s[4]+1/(s[3]-s[4])),ncol=1)
sum(s^2)-2
w <- v1%*%t(v2)
ee <- eigen(w)
ee$vectors%*%diag(ee$values)%*%t(ee$vectors)

(v1%*%t(v2) - v1%*%t(v1)/sum(v1^2)*as.numeric(t(v1)%*%v2))/w

n <- 1e2
B <- 1e3
eigens <- replicate(B,{
    s <- runif(n)
    z <- runif(n,-1,1)
    z <- rcauchy(n)
    vx <- -(s[1]-s[2])^2*s/sum(s^2); vx[1] <- vx[1]+s[1]-s[2]; vx[2] <- vx[2]-s[1]+s[2]
    vy <- -(s[3]-s[4])^2*s/sum(s^2); vx[3] <- vx[3]+s[3]-s[4]; vx[4] <- vx[3]-s[4]+s[2]
    m <- matrix(vx,ncol=1)%*%matrix(vy,nrow=1)
    m <- m+t(m)
    sort(eigen(m)$val)[c(1,n)]
})
summary(eigens[1,]);summary(eigens[2,])



n <- 1e5
z <- matrix(runif(2*n,-1,1),ncol=2)
hist(z[,1]-z[,2]-z[,1]*z[,2])
hist(z[,1]-z[,2]+z[,1]*z[,2],add=TRUE,col='red')


n <- 1
B <- 1e6
z <- matrix(runif(4*B,-1,1),ncol=4)
## z <- matrix(rnorm(4*B),ncol=4)
z <- matrix(rbeta(4*B,.5,.5),ncol=4)-1/2
s <- matrix(runif(4*B,0,1),ncol=4)
## mean(pnorm(sqrt(n)*pmin(z[,1]-z[,2],z[,3]-z[,4]) - rowSums(z)/sqrt(n)))
## mean(pnorm(sqrt(n)*pmin((z[,1]-z[,2])/(s[,1]-s[,2]),(z[,3]-z[,4])/(s[,3]-s[,4])) - rowSums(z*s)/sqrt(n)))
summary( z[,1]*z[,2] + 1/n * (z[,1]+z[,2]))
summary( z[,1]*z[,2] + z[,3]/n * (z[,1]+z[,2]))

plot(x=rbeta(1e5,.5,.5),y=rbeta(1e5,.5,.5),col=rgb(0,0,1,.03))

B <- 1e5
lim <- 4
c <- -1/5
x <- runif(B,-lim,lim)*.2; y <- runif(B,-lim,lim)
x <- runif(B,-lim,lim); y <- runif(B,-lim,lim)
idx <- (x*y+c*(x+y))< 0
plot(x[idx],y[idx]);abline(h=0,v=0,col='red')
## curve(-c*x/(x+c),add=TRUE,col='green')

op <- par(mfrow=c(1,2))
idx <- (x*y+c*(x+y))< -.2
plot(x[idx],y[idx],xlim=c(-lim,lim),asp=1);abline(h=0,v=0,col='red')
idx <- (x*y+c*(x+y))< -.253
plot(x[idx],y[idx],xlim=c(-lim,lim),asp=1);abline(h=0,v=0,col='red')
par(op)

B <- 1e5
x <- rexp(B)-log(2)
y <- rexp(B)-log(2)
median(x)
median(x+y)
plot(x,y,col=rgb(0,0,1,.2))
abline(a=0,b=-1,col='red')


B <- 1e6
n <- 10
z <- matrix(rnorm(B*n),ncol=n)
## z <- matrix(runif(B*n,-.5,.5),ncol=n)
s <- matrix(runif(B*n,0,5),ncol=n)
## s <- matrix(rexp(B*n,rate=.01),ncol=n)
theta.4 <- rowSums((z*s)[,1:4])/rowSums(s^2)
theta.rest <- rowSums((z*s)[,5:n])/rowSums(s^2)
theta <- theta.4+theta.rest
x <- (z[,1]-z[,2])/(s[,1]-s[,2])
y <- (z[,3]-z[,4])/(s[,3]-s[,4])
median(x*y - theta.rest*(x+y))
median(-theta.4*(x+y)+theta^2)
median(x*y - theta.rest*(x+y) -theta.4*(x+y)+theta^2)

idx <- x*y - theta.rest*(x+y) < 0
mean(idx)
plot(x[idx],y[idx])

require(parallel)
B <- 1e6
ns <- unique(round(seq(6,1e2,len=30)))
medians <- mclapply(ns, mc.cores=detectCores()-2,FUN=function(n) {          
  ## medians <- sapply(ns, FUN=function(n) {          
    z <- matrix(rnorm(B*n),ncol=n)
    ## z <- matrix(runif(B*n,-.5,.5),ncol=n)
    s <- matrix(runif(B*n,0,5),ncol=n)
    ## s <- matrix(rexp(B*n,rate=.01),ncol=n)
    theta.4 <- rowSums((z*s)[,1:4])/rowSums(s^2)
    theta.rest <- rowSums((z*s)[,5:n])/rowSums(s^2)
    theta.rest.ind <- rowSums((z*s)[,5:n])/rowSums(s[,5:n]^2)
    theta <- theta.4+theta.rest
    x <- (z[,1]-z[,2])/(s[,1]-s[,2])
    y <- (z[,3]-z[,4])/(s[,3]-s[,4])
    c(main=median(x*y - theta.rest*(x+y)),
    second=median(-theta.4*(x+y)+theta^2),
    full=median(x*y - theta.rest*(x+y) -theta.4*(x+y)+theta^2),
    prob=mean(x*y - theta.rest*(x+y) < 0))
})
medians <- simplify2array(medians)
plot(ns,medians['main',])
points(ns,medians['second',],col='red')
points(ns,medians['full',],col='blue')

plot(ns,medians['full',]*ns^2)

## using a loop, sapply taking too much memory
## decay of prob(A) with n
B <- 1e6
ns <- round(seq(6,300,len=40))
prob.A <- numeric()
#medians <- mclapply(ns, mc.cores=detectCores()-2,FUN=function(n) {          
## medians <- sapply(ns, FUN=function(n) {
for(n in ns) {
    print(n)
    z <- matrix(rnorm(B*n),ncol=n)
    ## z <- matrix(runif(B*n,-.5,.5),ncol=n)
    ## s <- matrix(runif(B*n,0,5),ncol=n)
    s <- matrix(rexp(B*n,rate=1),ncol=n)
    theta.4 <- rowSums((z*s)[,1:4])/rowSums(s^2)
    theta.rest <- rowSums((z*s)[,5:n])/rowSums(s^2)
    theta.rest.ind <- rowSums((z*s)[,5:n])/rowSums(s[,5:n]^2)
    theta <- theta.4+theta.rest
    x <- (z[,1]-z[,2])/(s[,1]-s[,2])
    y <- (z[,3]-z[,4])/(s[,3]-s[,4])
    ## x <- runif(B,-1,1);y <- runif(B,-1,1)
    prob.A.full <- mean(x*y - theta.rest*(x+y) < 0)
    prob.A.1 <- mean( (x>abs(theta.rest)) & (y>abs(theta.rest)) & (y< x*abs(theta.rest)/(x-abs(theta.rest))))
    prob.A <- cbind(prob.A,c(A.full=prob.A.full,A.1=prob.A.1))
}

op <- par(mfrow=c(1,2))
A.full <- prob.A['A.full',]-1/2
A.1 <- prob.A['A.1',]
for(g in c(function(n)n^(-1/4),function(n)log(n))) {
plot(ns,A.full,ylim=range(c(A.full,A.1)))
fit <- lm(A.full ~ I(g(ns)))
lines(ns,coef(fit)[1]+coef(fit)[2]*g(ns))
points(ns,A.1,col=2)
fit <- lm(A.1 ~ g(ns))
lines(ns,coef(fit)[1]+coef(fit)[2]*g(ns),col=2)
}
par(op)


## decay of prob(B) with n
B <- 1e6
ns <- round(seq(6,100,len=20))
prob.B <- numeric()
#medians <- mclapply(ns, mc.cores=detectCores()-2,FUN=function(n) {          
## medians <- sapply(ns, FUN=function(n) {
for(n in ns) {
    print(n)
    z <- matrix(rnorm(B*n),ncol=n)
    ## z <- matrix(runif(B*n,-.5,.5),ncol=n)
    ## s <- matrix(runif(B*n,0,5),ncol=n)
    s <- matrix(rexp(B*n,rate=1),ncol=n)
    theta.4 <- rowSums((z*s)[,1:4])/rowSums(s^2)
    theta.rest <- rowSums((z*s)[,5:n])/rowSums(s^2)
    theta.rest.ind <- rowSums((z*s)[,5:n])/rowSums(s[,5:n]^2)
    theta.hat <- theta.4+theta.rest
    x <- (z[,1]-z[,2])/(s[,1]-s[,2])
    y <- (z[,3]-z[,4])/(s[,3]-s[,4])
    term.A <- x*y - theta.rest.ind*(x+y)
    term.B <- theta.4*(theta.4-1)*(x+y)+theta.hat^2
    term.A <- x*y; term.B <- theta.hat^2 #!!!!!!
    prob.B <- c(prob.B, mean( abs(term.A)<abs(term.B)  ))
}

plot(ns,prob.B)
g <- function(n)n^(-1)
fit <- lm(prob.B ~ I(g(ns)))
lines(ns,coef(fit)[1]+coef(fit)[2]*g(ns))


## checking fitting procedure as a measure of convergence rate
B <- 1e2
ns <- seq(1,200,len=40)
## means <- sapply(ns, function(n) mean(replicate(B,abs(mean(rnorm(n))))))
## means <- sapply(ns, function(n) mean(replicate(B,mean(rnorm(n))^2)))
means <- sapply(ns, function(n) mean(replicate(B,mean(rnorm(n))^2 > .1)))
plot(ns,means)
g <- function(n)n^(-1)
fit <- lm(means ~ I(g(ns)))
lines(ns,coef(fit)[1]+coef(fit)[2]*g(ns))

dd






## prob(A) for uniform x,y--exact formulas
B <- 1e6
cs <- -seq(.001,.05,len=20)
prob.A <- numeric()
#medians <- mclapply(ns, mc.cores=detectCores()-2,FUN=function(n) {          
## medians <- sapply(ns, FUN=function(n) {
for(c in cs) {
    print(c)
    x <- runif(B,-1,1);y <- runif(B,-1,1)
    prob.A.full <- mean(x*y - c*(x+y) < 0)
    prob.A.1 <- mean( (x>abs(c)) & (y>abs(c)) & (y< x*abs(c)/(x-abs(c))))
    ## prob.A.1 <- mean( (x>abs(c)) & (y>abs(c)) & (y< c^2/(x-abs(c))+abs(c)))
    prob.A.2 <- mean( (x>0) & (y>c*x/(c-x)) & (y< pmin(abs(c),x))  )
    prob.A <- cbind(prob.A,c(A.full=prob.A.full,A.1=prob.A.1,A.2=prob.A.2))
}

A.full <- prob.A['A.full',]-1/2; A.1 <- prob.A['A.1',]; A.2 <- prob.A['A.2',]#; A.3 <- prob.A['A.3',]
plot(cs,A.full,ylim=range(c(A.full,A.1)))
points(cs,A.1,col=2)
points(cs,A.2,col=3)
legend('topright',col=1:3,legend=c('full','2','3'),pch=1)
## lines(cs,A.1+2*A.2) # A.full = A.1+2*A.2
abline(h=0)
g  <- function(n)abs(n)^(-1/4)
g <- function(n)log(abs(n))
fit <- lm(rev(A.full) ~ I(g(cs)))
## lines(rev(cs),coef(fit)[1]+coef(fit)[2]*g(cs))
lines(cs,1/4*cs^2*(1+2*abs(log(abs(cs)))))
lines(cs,1/4*(-cs^2/2+cs^2*log((1+abs(c))/abs(c))))

idx <-  (x>abs(c)) & (y>abs(c)) & (y< x*abs(c)/(x-abs(c)))
plot(x[idx],y[idx],pch='.')

dd

c <- 3
delta <- c/10
curve( (x+delta)*c / (x-(c-delta)) - delta,from=c-delta,to=c-delta+1)
curve(x*(c-delta)/(x-(c-delta)),add=TRUE)

require(mvtnorm)
B <- 1e6
## x <- runif(B,-1,1)
## y <- runif(B,-1,1)
rho <- -.2
xy <- rmvnorm(B,sigma=matrix(c(1,rho,rho,1),nrow=2))
x  <- xy[,1]; y <- xy[,2]
mean(punif(pmin(x,y),-1,1))-1/2

## uniform distr bound too fast, using A_n density

B <- 1e5
lim <- 4
c <- -.5
x <- runif(B,-lim,lim); y <- runif(B,-lim,lim)
A <- (x*y+c*(x+y))
plot(x[A<0],y[A<0]);abline(h=0,v=0,col='red')
points(x[-c^2<A & A<c^2],y[-c^2<A & A<c^2],col=rgb(1,0,0,.1))
curve(-c*x/(x+c),add=TRUE,col='green'); abline(h=abs(c),col='green')


B <- 1e5
lim <- 1
c <- -.5
cs <- seq(.005,.15,len=50)
areas <- sapply(cs, function(c) {
x <- runif(B,-lim,lim); y <- runif(B,-lim,lim)
## x <- rnorm(B); y <- rnorm(B)
A <- (x*y+c*(x+y))
c(mean(A<0) - 1/2, mean(-c^2<A & A<c^2))
})
plot(cs,areas[1,])
points(cs,areas[2,],col='red')

curve(1/4*((2*c^2+c*(c-1)^2)/(1+c) - 2*c^2*log(2*c^2)),xname='c',add=TRUE)





B <- 1e5
n <- 10
ns <- round(seq(6,300,len=40))
prob.A <- numeric()
##medians <- mclapply(ns, mc.cores=detectCores()-2,FUN=function(n) {     
## medians <- sapply(ns, FUN=function(n) {
## for(n in ns) {
## print(n)
z <- matrix(rnorm(B*n),ncol=n)
z <- matrix(runif(B*n,-.5,.5),ncol=n)
s <- matrix(runif(B*n,0,5),ncol=n)
s <- matrix(rexp(B*n,rate=1),ncol=n)
## theta.4 <- rowSums((z*s)[,1:4])/rowSums(s^2)
## theta.rest <- rowSums((z*s)[,5:n])/rowSums(s^2)
## theta.rest.ind <- rowSums((z*s)[,5:n])/rowSums(s[,5:n]^2)
## theta <- theta.4+theta.rest
theta.hat <- rowSums(z*s)/rowSums(s^2)
theta.del  <- rowSums((z*s)[,5:n])/rowSums(s[,1:n]^2)
theta.ind <- rowSums((z*s)[,5:n])/rowSums(s[,5:n]^2)
x <- (z[,1]-z[,2])/(s[,1]-s[,2])
y <- (z[,3]-z[,4])/(s[,3]-s[,4])
## x <- runif(B,-1,1);y <- runif(B,-1,1)
## prob.A.full <- mean(x*y - theta.rest*(x+y) < 0)
## prob.A.1 <- mean( (x>abs(theta.rest)) & (y>abs(theta.rest)) & (y< x*abs(theta.rest)/(x-abs(theta.rest))))
cc <- -.5
## median(x*y - theta.hat*(x+y))
## median(x*y - theta.del*(x+y))
median(x*y - theta.ind*(x+y))
median(theta.hat^2)
median(-(theta.hat - theta.ind)*(x+y))
median(-theta.hat*(x+y))
median(-(theta.hat - theta.ind)*(x+y) + theta.hat^2)
median(-theta.hat*(x+y) + theta.hat^2)
median(x*y - theta.hat*(x+y) + theta.hat^2)
mean(x*y - theta.ind*(x+y) < 0)-1/2
mean(-(theta.hat - theta.ind)*(x+y) + theta.hat^2 < 0) - 1/2
mean(x*y - theta.hat*(x+y) +theta.hat^2 < 0)-1/2
mean(x*y - theta.ind*(x+y) +theta.hat^2 < 0)-1/2



B <- 3e5
n <- 10
ns <- round(seq(6,150,len=20))
probs <- numeric()
## probs <- matrix(nrow=2,ncol=0)
##medians <- mclapply(ns, mc.cores=detectCores()-2,FUN=function(n) {     
## medians <- sapply(ns, FUN=function(n) {
for(n in ns) {
    print(n)
    z <- matrix(rnorm(B*n),ncol=n)
    ## z <- matrix(runif(B*n,-.5,.5),ncol=n)
    ## s <- matrix(runif(B*n,0,5),ncol=n)
    s <- t(matrix(rbeta(n,runif(1),runif(1)),nrow=n,ncol=B)) # condit on s
    ## s <- matrix(rexp(B*n,rate=1),ncol=n)
    theta.hat <- rowSums(z*s)/rowSums(s^2)
    ## theta.del  <- rowSums((z*s)[,5:n])/rowSums(s[,1:n]^2)
    theta.ind <- rowSums((z*s)[,5:n])/rowSums(s[,5:n]^2)
    x <- (z[,1]-z[,2])/(s[,1]-s[,2])
    y <- (z[,3]-z[,4])/(s[,3]-s[,4])
    term1 <- x*y + theta.ind*(x+y)
    term2 <- (theta.ind-theta.hat)*(x+y)+theta.hat^2
    ## probs <- c(probs,mean((theta.ind-theta.hat)*(x+y)+theta.hat^2 > 0))
    ## probs <- c(probs, median((theta.ind-theta.hat)*(x+y)+theta.hat^2))
    ## probs <- cbind(probs,c(mean(term1<0),mean((term1<0) & (term2>0) & (abs(term1)<abs(term2)))))
    ## probs <- cbind(probs,mean((term1<0) & (abs(term1)>abs(term2))))
    probs <- cbind(probs,mean(term1 + term2 < 0))
    ## probs <- cbind(probs,median(x*y - theta.hat*(x+y) + theta.hat^2))
    ## probs <- cbind(probs,)
}
plot(ns,probs)

plot(ns,probs[1,]-1/2)
points(ns,probs[2,],col='red')






## quad form
n <- 3
s <- runif(4)
sqrt((2*n^2+(s[1]-s[2])^2*(sum(s^2)-2*n))*(2*n^2+(s[3]-s[4])^2*(sum(s^2)-2*n)))
op <- par(mfrow=c(1,2))
curve((-2*n + sum(s^2)+  sqrt((2*n^2+(s[1]-s[2])^2*(sum(s^2)-2*n))*(2*n^2+(s[3]-s[4])^2*(sum(s^2)-2*n)))/((s[1]-s[2])*(s[3]-s[4])))  /n^2,xname='n',1,20)
curve((-2*n + sum(s^2)-  sqrt((2*n^2+(s[1]-s[2])^2*(sum(s^2)-2*n))*(2*n^2+(s[3]-s[4])^2*(sum(s^2)-2*n)))/((s[1]-s[2])*(s[3]-s[4])))  /n^2,xname='n',1,20)
par(op)


S1 <- c(1/(s[1] - s[2]) - s[1]/n, -1/(s[1] - s[2]) - s[2]/n, -s[3]/n, -s[4]/n)
S2 <- c(-s[1]/n, -s[2]/n, 1/(s[3] - s[4]) - s[3]/n, -1/(s[3] - s[4]) - s[4]/n)
SS <- (S1%*%t(S2) + t(S1%*%t(S2)))/2
sort(eigen(SS)$val)
(-2*n + sum(s^2)-  sqrt((2*n^2+(s[1]-s[2])^2*(sum(s^2)-2*n))*(2*n^2+(s[3]-s[4])^2*(sum(s^2)-2*n)))/((s[1]-s[2])*(s[3]-s[4])))  /(2*n^2)
(-2*n + sum(s^2)+  sqrt((2*n^2+(s[1]-s[2])^2*(sum(s^2)-2*n))*(2*n^2+(s[3]-s[4])^2*(sum(s^2)-2*n)))/((s[1]-s[2])*(s[3]-s[4])))  /(2*n^2)

eig <- eigen(SS)
vv <- c(-(n+s[1]*s[2]-s[2]^2+s[3]*s[4]-s[4]^2)/((s[1]-s[2])*(s[3]-s[4])),-(n-s[1]^2+s[1]*s[2]+s[3]*s[4]-s[4]^2),0,1)
vv/sqrt(sum(vv^2))

B <- 1e6
Z <- matrix(runif(4*B,-3,3),nrow=4)
## Z <- matrix(rnorm(4*B),nrow=4)
summary(as.numeric(eig$vec[,1]%*%Z))
summary(as.numeric(eig$vec[,4]%*%Z))
qqplot(as.numeric(eig$vec[,1]%*%matrix(rnorm(4*B),nrow=4)),as.numeric(eig$vec[,4]%*%matrix(rnorm(4*B),nrow=4)));abline(a=0,b=1,col='red')

n <- 6
B <- 1e3
data <- split(matrix(runif(B*n),nrow=B),1:B)
data <- split(matrix(rexp(B*n),nrow=B),1:B)
evs <- sapply(data,function(s) {
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    (( ev1 <- eig$vec[,1])); ((ev2 <- eig$vec[,n] ))
    ## sum(sign(ev1)*sign(ev2)>0)
    ## qqplot(ev1%*%Z,ev2%*%Z); abline(a=0,b=1,col=2)
    cbind(ev1=ev1,ev2=ev2)
},simplify=FALSE)
evs <- simplify2array(evs)
signs <- abs(colSums(sign(evs)[,1,]*sign(evs)[,2,]))
table(signs)

idx <- which(signs!=n-4)[1]
ev1 <- evs[,1,idx]; ev2 <- evs[,2,idx]
B <- 1e6
Z <- matrix(runif(n*B,-1,1),nrow=n)
Z <- matrix(rexp(n*B),nrow=n)-1
proj1 <- as.numeric(ev1%*%Z); proj2 <- as.numeric(ev2%*%Z)
summary(abs(proj1)-abs(proj2))

plot(proj1,proj2,asp=1)
points(proj2,proj1,col=rgb(1,0,0,.2))
abline(h=0,v=0); abline(a=0,b=1)

## (Z^Tev1,Z^Tev2) not exchangeable when signs of ev1 and ev2 dont
## line up

n <- 4
B <- 1e5
Z <- matrix(runif(n*B,-1,1),nrow=n)
## Z <- matrix(rcauchy(n*B),nrow=n)
## Z <- matrix(rt(n*B,df=2),nrow=n)
## Z <- matrix(runif(n*B,0,1),nrow=n)
## Z <- matrix(rexp(n*B,2),nrow=n)-1/2
## Z <- matrix(rgamma(n*B,1,30),nrow=n)-1/30
## s <- runif(n)
s <- rexp(n)
Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
eig <- eigen(SSmat)
eig$vec[,1]; eig$vec[,n]

stmp <- s[3]; s[3] <- s[4]; s[4] <- stmp
stmp <- s[1]; s[1] <- s[2]; s[2] <- stmp
(( ev1 <- eig$vec[,1])); ((ev2 <- eig$vec[,n] ))
sign(ev1)*sign(ev2)
## qqplot(ev1%*%Z,ev2%*%Z); abline(a=0,b=1,col=2)

proj1 <- as.numeric(ev1%*%Z); proj2 <- as.numeric(ev2%*%Z)
## mean((proj1+proj2>2)[proj1>0 & proj2>0])
## ev1+ev2
## mean(((ev1+ev2)%*%Z > 2)[ev1%*%Z>0 & ev2%*%Z>0])
plot(proj1,proj2,asp=1)
points(proj2,proj1,col=rgb(1,0,0,.3))
## abline(a=0,b=-1);abline(a=0,b=1)

hist(as.numeric(eig$vec[,1]%*%Z))
hist(-as.numeric(eig$vec[,1]%*%Z),add=TRUE,col=rgb(1,0,0,.1))


## checking if exchangeable [seems no]
s <- runif(n)
Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
(( eig <- eigen(SSmat) ))
a <- 3;b <- 4; stmp <- s[a]; s[a] <- s[b]; s[b] <- stmp
## stmp <- s[1]; s[1] <- s[2]; s[2] <- stmp
Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
(( eig <- eigen(SSmat) ))

pairs <- replicate(1e3,{
    s <- runif(n,1,2)
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    ## c(eig$val[1],eig$vec[1,1],eig$vec[1,4])
    c(eig$val[1],sum(eig$vec[,1]),sum(eig$vec[,4]))
})
plot(pairs[1,],pairs[2,])
points(pairs[1,],pairs[3,],col='red')

qqplot(pairs[2,],pairs[3,]);abline(a=0,b=1)


## magnitude of proj by n
ns <- round(seq(4,40,len=20))
B <- 1e6
by.n <- sapply(ns, simplify=FALSE, FUN=function(n) {
    projections <- replicate(1e2, {
    s <- runif(n,1,2)
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    ev1 <- eig$vect[,1]; ev2 <- eig$vect[,n]
    Z <- runif(n,-1,1)
    ## Z <- rexp(n)-1
    proj1 <- as.numeric(ev1%*%Z); proj2 <- as.numeric(ev2%*%Z)
    c(proj1,proj2)
    })
})

diffs <- sapply(by.n,function(m)m[1,]^2-m[2,]^2)
matplot(t(diffs),col=1,pch='.',cex=2)
plot(apply(diffs,2,sd),type='l')

lens <- sapply(by.n,function(m)abs(m[1,]))
matplot(t(lens),col=1,pch='.',cex=2)


## median of diff of squared projections
ns <- round(seq(6,200,len=30))
by.n <- sapply(ns, function(n) {
    cat('.')
    diffs <- replicate(4e2, {
        ## s <- runif(n,0,3)
        ## s <- rexp(n)
        s <- rgamma(n,shape=.54)
        ## z <- rnorm(n)
        z <- runif(n,-4,4)
        ## z <- rexp(n)-1
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        (z%*%eig$vec[,1])^2 - (z%*%eig$vec[,n])^2
        ## (z%*%eig$vec[,1])^2*eig$val[1]+ (z%*%eig$vec[,n])^2*eig$val[n]
        ## (z%*%Svec1)*(z%*%Svec2)
    })
    mean(diffs<0)-1/2
})
plot(ns,by.n)

## median of diff of squared projections -- by B [both A and B terms > 1/2? for cauchy/t <1/2 though.]
require(parallel)
n <- 20
Bs <- round(seq(50,3e3,len=50))
by.B <- mclapply(Bs, mc.cores=detectCores()-2,FUN=function(B) {
    cat('.')
    diffs <- replicate(B, {
        s <- runif(n,0,3)
        ## s <- rexp(n)
        ## s <- rgamma(n,shape=.54)
        ## z <- rnorm(n)
        ## z <- runif(n,-4,4)
        ## z <- rt(n,df=4)
        z <- rexp(n)-1/2
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        N <- sum(s^2)
        A <- -1/(2*N)
        B <- sqrt(N^2*(-2*N+(s[1]-s[2])^2)*(-2*N+(s[3]-s[4])^2))/(s[1]-s[2])/(s[3]-s[4])/2/N^2
        projsqr1 <- (z%*%eig$vec[,1])^2; projsqr2 <- (z%*%eig$vec[,2])^2
        c(first=A*(projsqr1+projsqr2),
          ## second=abs(B)*(projsqr1-projsqr2),
          second=(projsqr1-projsqr2),
          full=A*(projsqr1+projsqr2)+abs(B)*(projsqr1-projsqr2))
        ## abs(B)*((z%*%eig$vec[,1])^2 - (z%*%eig$vec[,n])^2)
        ## (z%*%eig$vec[,1])^2*eig$val[1]+ (z%*%eig$vec[,n])^2*eig$val[n]
        ## (z%*%Svec1)*(z%*%Svec2)
    })
    rowMeans(diffs<0)-1/2
})
by.B <- simplify2array(by.B)
plot(Bs,by.B['second',])
points(Bs,by.B['full',],col=2)
abline(lm(by.B['second',]~Bs))
abline(h=0,lty=2)

mm <- sapply(Bs, function(B)mean(rnorm(B)))
plot(Bs,mm)
abline(lm(mm~Bs))
abline(h=0,lty=2)





## checking fla for eigenvalues
n <- 4
for(jj in 1:1e4) {
    s <- rgamma(n,shape=.54)
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    N <- sum(s^2)
    A <- -1/(2*N)
    B <- sqrt(N^2*(-2*N+(s[1]-s[2])^2)*(-2*N+(s[3]-s[4])^2))/(s[1]-s[2])/(s[3]-s[4])/2/N^2
    z <- rnorm(4)
    projsqr1 <- (z%*%eig$vec[,1])^2; projsqr2 <- (z%*%eig$vec[,n])^2
    stopifnot(abs(A*(projsqr1+projsqr2) +   abs(B)*(projsqr1-projsqr2) - 
                  t(z)%*%SSmat%*%z) < 1e-8)
    theta.hat <- sum(z*s)/sum(s^2)
    stopifnot(abs(((z[1]-z[2])/(s[1]-s[2])-theta.hat)*((z[3]-z[4])/(s[3]-s[4])-theta.hat) - t(z)%*%SSmat%*%z)<1e-8)
    ## stopifnot(abs(A+B - eig$val[1]) < 1e-8)
    ## stopifnot(abs(A-B - eig$val[n]) < 1e-8)
    ## jcf1 <- c((-s[1]*s[2]-s[2]^2-s[3]^2-s[4]^2)/((s[1]-s[2])*s[4]),
    ## (s[1]*s[2]+s[1]^2+s[3]^2+s[4]^2)/((s[1]-s[2])*s[4]),
    ## s[3]/s[4],1)
    ## jcf2 <- c(-s[1]*s[4]/sum(s^2),-s[2]*s[4]/sum(s^2),-1+s[3]*(1/(s[3]-s[4])-s[4]/sum(s^2)), s[3]/(s[4]-s[3])+(s[1]^2+s[2]^2+s[3]^2)/sum(s^2))
    ## jcf.val <- -1/sum(s^2)
    ## (z%*%jcf1)*(z%*%jcf2)*jcf.val
}

## jcf [offers no advantage]
require(parallel)
n <- 4
projections <- replicate(1e5, {
    s <- runif(n,0,3)
    ## s <- rexp(n)
    ## s <- rgamma(n,shape=.54)
    ## z <- rnorm(n)
    z <- runif(n,-4,4)
    ## z <- rexp(n)-1/2
    jcf1 <- c((-s[1]*s[2]-s[2]^2-s[3]^2-s[4]^2)/((s[1]-s[2])*s[4]),
    (s[1]*s[2]+s[1]^2+s[3]^2+s[4]^2)/((s[1]-s[2])*s[4]),
    s[3]/s[4],1)
    jcf2 <- c(-s[1]*s[4]/sum(s^2),-s[2]*s[4]/sum(s^2),-1+s[3]*(1/(s[3]-s[4])-s[4]/sum(s^2)), s[3]/(s[4]-s[3])+(s[1]^2+s[2]^2+s[3]^2)/sum(s^2))
    c(proj1=jcf1%*%z,proj2=jcf2%*%z)
})
mean(projections[1,]*projections[2,]>0)

dd


## cant detect bias, trying by B
require(parallel)
## B <- 1e3
n <- 5
Bs <- unique(round(seq(1e2,4e3,len=1e2)))
biases <- mclapply(Bs, mc.cores=detectCores()-2,FUN=function(B) {
    ## medians <- sapply(ns, FUN=function(n) {          
    z <- matrix(rnorm(B*n),ncol=n)
    ## z <- matrix(runif(B*n,-10,10),ncol=n)
    ## a <- .1;z <- matrix(rbeta(B*n,a,a),ncol=n)-1/2
    ## s <- matrix(runif(B*n,0,1),ncol=n)
    s <- matrix(rgamma(B*n,shape=.54),ncol=n)
    theta.hat <- rowSums(z*s)/rowSums(s^2)
    ## median(((z[,1]-z[,2])-theta.hat*(s[,1]-s[,2]))*((z[,3]-z[,4])-theta.hat*(s[,3]-s[,4]))*(s[,1]-s[,2])*(s[,3]-s[,4]))
    mean(((z[,1]-z[,2])-theta.hat*(s[,1]-s[,2]))*((z[,3]-z[,4])-theta.hat*(s[,3]-s[,4]))*(s[,1]-s[,2])*(s[,3]-s[,4])<0)-1/2
    ## x <- ((z[,1]-z[,2])-theta.hat*(s[,1]-s[,2]))*((z[,3]-z[,4])-theta.hat*(s[,3]-s[,4]))*(s[,1]-s[,2])*(s[,3]-s[,4])
    ## x <- split(x,rep(1:1e3,each=B/1e3))
    ## mean(sapply(x,median))
})
biases <- simplify2array(biases)
## plot(Bs,biases['main',])
## points(Bs,biases['second',],col='red')
plot(Bs,biases,col='blue')


require(parallel)
B <- 1e3
ns <- unique(round(seq(6,3e2,len=1e3)))
medians <- mclapply(ns, mc.cores=detectCores()-2,FUN=function(n) {          
    ## medians <- sapply(ns, FUN=function(n) {          
    ## z <- matrix(rnorm(B*n),ncol=n)
    z <- matrix(runif(B*n,-10,10),ncol=n)
    ## a <- .1;z <- matrix(rbeta(B*n,a,a),ncol=n)-1/2
    s <- matrix(runif(B*n,0,1),ncol=n)
    ## s <- matrix(rgamma(B*n,shape=.54),ncol=n)
    theta.hat <- rowSums(z*s)/rowSums(s^2)
    ## median(((z[,1]-z[,2])-theta.hat*(s[,1]-s[,2]))*((z[,3]-z[,4])-theta.hat*(s[,3]-s[,4]))*(s[,1]-s[,2])*(s[,3]-s[,4]))
    mean(((z[,1]-z[,2])-theta.hat*(s[,1]-s[,2]))*((z[,3]-z[,4])-theta.hat*(s[,3]-s[,4]))*(s[,1]-s[,2])*(s[,3]-s[,4])<0)-1/2
    ## x <- ((z[,1]-z[,2])-theta.hat*(s[,1]-s[,2]))*((z[,3]-z[,4])-theta.hat*(s[,3]-s[,4]))*(s[,1]-s[,2])*(s[,3]-s[,4])
    ## x <- split(x,rep(1:1e3,each=B/1e3))
    ## mean(sapply(x,median))
})
medians <- simplify2array(medians)
## plot(ns,medians['main',])
## points(ns,medians['second',],col='red')
plot(ns,medians*ns^(0),col='blue')



## B term properties
require(parallel)
n <- 5
Bs <- round(seq(50,3e3,len=50))
by.B <- mclapply(Bs, mc.cores=detectCores()-2,FUN=function(B) {
    cat('.')
    diffs <- replicate(B, {
        s <- runif(n,0,3)
        ## s <- rgamma(n,shape=.54)
        ## z <- rnorm(n)
        z <- runif(n,-4,4)
        ## z <- rcauchy(n)
        ## z <- rexp(n)-1/2
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        ## eig$vec[1:4,1]; eig$vec[1:4,n]        
        idx <- !((1:n)%in%(c(2,3)))
        idx <- !((1:n)%in%(c(1,4)))
        idx <- (1:n)>4
        ## eig$vec[idx,1] - eig$vec[idx,n] # [>0]
        idx <- (1:n)<=4
        ## idx <- 1:n
        proj1 <- z[idx]%*%eig$vec[idx,1]; proj2 <- z[idx]%*%eig$vec[idx,n]
        abs(proj1) -  abs(proj2)
    })
    mean(diffs<0)
})
by.B <- simplify2array(by.B)
plot(Bs,by.B); abline(h=1/2,lty=2)

points(Bs,by.B,col='red')


require(parallel)
n <- 6
## Bs <- round(seq(50,3e3,len=50))
## by.B <- mclapply(Bs, mc.cores=detectCores()-2,FUN=function(B) {
ns <- round(seq(6,100,len=10))
by.n <- mclapply(ns, mc.cores=detectCores()-2,FUN=function(n) {
    cat('.')
    eigvecs <- replicate(1e2, {
        s <- runif(n,0,3)
        ## s <- rgamma(n,shape=.54)
        ## z <- rnorm(n)
        ## z <- runif(n,-4,4)
        ## z <- rexp(n)-1/2
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        cbind(eig$vec[,1], eig$vec[,n]        )
    })
    apply(abs(eigvecs),1:2,mean)
})

op <- par(mfrow=c(1,2))
firstfour <- sapply(by.n,function(m)m[1:4,2])
matplot(t(firstfour),col=1,pch=1)
overfour <- sapply(1:length(by.n),function(j)cbind(ns[j],by.n[[j]][5:nrow(by.n[[j]]),1]))
overfour <- do.call(rbind,overfour)
plot(overfour[,1],overfour[,2])
par(op)

overfour <- sapply(1:length(by.n),function(j)cbind(ns[j],by.n[[j]][5:nrow(by.n[[j]]),2]))
overfour <- do.call(rbind,overfour)
plot(overfour[,1],overfour[,2]*overfour[,1]^(1)) ## [1/n decay]

overfoursum <- sapply(by.n,function(m)sum(m[5:nrow(m),1]))
plot(ns,overfoursum)

## difference in magnitude of eigenvector components
require(parallel)
n <- 6
Bs <- round(seq(50,3e3,len=50))
by.B <- mclapply(Bs, mc.cores=detectCores()-2,FUN=function(B) {
    cat('.')
    diffs <- replicate(B, {
        s <- runif(n,0,3)
        ## s <- rgamma(n,shape=.54)
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        ## eig$vec[1:4,1]; eig$vec[1:4,n]        
        ## proj1 <- z[idx]%*%eig$vec[idx,1]; proj2 <- z[idx]%*%eig$vec[idx,2]
        abs(eig$vec[4,1])-abs(eig$vec[1,1])
    })
    mean(diffs<0)
})
by.B <- simplify2array(by.B)
plot(Bs,by.B);abline(h=1/2)

## boxplots first 4 eigenvec components
require(parallel)
n <- 6
B <- 3e3
ns <- round(seq(50,3e2,len=3))
## by.n <- mclapply(Bs, mc.cores=detectCores()-2,FUN=function(B) {
##     cat('.')
    firstfour <- replicate(B, {
        s <- runif(n,0,1/3)
        ## s <- rgamma(n,shape=.54)
        ## s <- rexp(n)
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        ## eig$vec[1:4,1]; eig$vec[1:4,n]        
        ## proj1 <- z[idx]%*%eig$vec[idx,1]; proj2 <- z[idx]%*%eig$vec[idx,2]
        cbind(eig$vec[1:4,1],eig$vec[1:4,n])
    })
op <- par(mfrow=c(2,1))
boxplot(t(abs(firstfour[,1,])));abline(h=1/2)
boxplot(t(abs(firstfour[,2,])));abline(h=1/2)
par(op)
op <- par(mfrow=c(2,4))
for(j in 1:2) for(k in 1:4)hist(firstfour[k,j,],prob=TRUE)
par(op)

plot(firstfour[1,1,],firstfour[1,2,])

## boxplots mod of first 4 and over 4 parts
require(parallel)
n <- 30
B <- 3e3
ns <- round(seq(50,3e2,len=3))
## by.n <- mclapply(Bs, mc.cores=detectCores()-2,FUN=function(B) {
##     cat('.')
    twoparts <- replicate(B, {
        s <- runif(n,0,1)
        ## s <- rgamma(n,shape=.54)
        s <- rexp(n)
        z <- runif(n,-1,1)
        ## z <- rcauchy(n)
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        ## eig$vec[1:4,1]; eig$vec[1:4,n]        
        sapply(list(lower=1:4,upper=5:n), function(idx) c(abs(z[idx]%*%eig$vec[idx,1]),abs(z[idx]%*%eig$vec[idx,n])))        
    })
twoparts <- cbind(t(twoparts[,'lower',]),t(twoparts[,'upper',]))
boxplot(twoparts)

## hist of projections
require(parallel)
n <- 6
B <- 5e3
ns <- round(seq(50,3e2,len=3))
## by.n <- mclapply(Bs, mc.cores=detectCores()-2,FUN=function(B) {
##     cat('.')
projections <- replicate(B, {
    s <- runif(n,0,1)
    ## s <- rgamma(n,shape=.54)
    s <- rexp(n)
    ## z <- runif(n,-1,1)
    z <- rbeta(n,1,1)-1/2
    ## z <- rt(n,df=40)
    ## z <- rnorm(n)
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    c(pos=eig$vec[,1]%*%z,neg=eig$vec[,n]%*%z)
})
op <- par(mfrow=c(1,2))
hist(projections['pos',],prob=TRUE)
hist(projections['neg',],add=TRUE,prob=TRUE,col=rgb(1,0,0,.1))
hist(projections['pos',]^2,prob=TRUE)
hist(projections['neg',]^2,add=TRUE,prob=TRUE,col=rgb(1,0,0,.1))
par(op)
plot(ecdf(projections[1,]^2))
lines(ecdf(projections[2,]^2),col='red')
## op <- par(mfrow=c(2,1))
## plot(projections['pos',],projections['neg',])
## plot(projections['pos',]^2,projections['neg',]^2)
## par(op)
plot(projections['pos',]^2,projections['neg',]^2,asp=1,col=rgb(0,0,0,.5));abline(0,1)
points(projections['neg',]^2,projections['pos',]^2,col=rgb(1,0,0,.5))

## upper and lower parts of projections
require(parallel)
n <- 50
B <- 3e4
ns <- round(seq(10,1e2,len=100))
by.n <- mclapply(ns, mc.cores=detectCores()-2,FUN=function(n) {
    cat('.')
    proj <- replicate(B, {
        s <- runif(n,0,1)
        ## s <- rgamma(n,shape=.54)
        s <- rexp(n)
        ## z <- runif(n,-1,1)
        z <- rbeta(n,.5,.5)-1/2
        ## z <- rt(n,df=40)
        ## z <- rnorm(n)
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        c(pos4=eig$vec[1:4,1]%*%z[1:4],neg4=eig$vec[1:4,n]%*%z[1:4],posn=eig$vec[5:n,1]%*%z[5:n],negn=eig$vec[5:n,n]%*%z[5:n])
    })
    ## boxplot(t(abs(proj)))
    possq <- cbind(pos4sq=proj['pos4',]^2,poscross=proj['pos4',]*proj['posn',],posnsq=proj['posn',]^2)
    negsq <- cbind(neg4sq=proj['neg4',]^2,negcross=proj['neg4',]*proj['negn',],negnsq=proj['negn',]^2)
    c(colMeans(possq),colMeans(negsq))
})
by.n <- simplify2array(by.n)
plot(ns,by.n['neg4sq',]-by.n['pos4sq',],type='l')
lines(ns,by.n['negcross',]-by.n['poscross',],col=2)
lines(ns,by.n['negnsq',]-by.n['posnsq',],col=3)
lines(ns,by.n['neg4sq',]+by.n['negcross',]+by.n['negnsq',] - (by.n['pos4sq',]+by.n['poscross',]+by.n['posnsq',]),col=4,lty=2)
abline(h=0,lty=2)
mean(by.n['pos4sq',]<by.n['neg4sq',])
mean(by.n['neg4sq',]+by.n['negcross',]+by.n['negnsq',] > by.n['pos4sq',]+by.n['poscross',]+by.n['posnsq',])

matplot(t(by.n),pch=1)

## upper and lower parts of projections--by B
require(parallel)
n <- 40
B <- 1e2
Bs <- round(seq(5e2,1e4,len=100))
by.B <- mclapply(Bs, mc.cores=detectCores()-2,FUN=function(B) {
    cat('.')
    proj <- replicate(B, {
        s <- runif(n,0,1)
        ## s <- rgamma(n,shape=.54)
        s <- rexp(n)
        ## z <- runif(n,-1,1)
        z <- rbeta(n,.5,.5)-1/2
        ## z <- rt(n,df=40)
        ## z <- rnorm(n)
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        c(pos4=eig$vec[1:4,1]%*%z[1:4],neg4=eig$vec[1:4,n]%*%z[1:4],posn=eig$vec[5:n,1]%*%z[5:n],negn=eig$vec[5:n,n]%*%z[5:n])
    })
    ## boxplot(t(abs(proj)))
    possq <- cbind(pos4sq=proj['pos4',]^2,poscross=proj['pos4',]*proj['posn',],posnsq=proj['posn',]^2)
    negsq <- cbind(neg4sq=proj['neg4',]^2,negcross=proj['neg4',]*proj['negn',],negnsq=proj['negn',]^2)
    c(colMeans(possq),colMeans(negsq))
})
by.B <- simplify2array(by.B)
plot(Bs,by.B['neg4sq',]-by.B['pos4sq',],type='l')
lines(Bs,by.B['negcross',]-by.B['poscross',],col=2)
lines(Bs,by.B['negnsq',]-by.B['posnsq',],col=3)
lines(Bs,by.B['neg4sq',]+by.B['negcross',]+by.B['negnsq',] - (by.B['pos4sq',]+by.B['poscross',]+by.B['posnsq',]),col=4,lty=2)
abline(h=0,lty=2)
mean(by.B['pos4sq',]<by.B['neg4sq',])
mean(by.B['neg4sq',]+by.B['negcross',]+by.B['negnsq',] > by.B['pos4sq',]+by.B['poscross',]+by.B['posnsq',])

dd

## decay of lower eigenvector components [ans: 1/n] 
n <- 40
B <- 1e5
ns <- round(seq(5,200,len=20))
by.n <- sapply(ns, function(n) {
    lower.components <- replicate(1e2, {
        ## s <- runif(n)
        s <- rexp(n)
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        ## eig$vec[,1]; eig$vec[,n]
        abs(eig$vec[5:n,1])
    })
    mean(lower.components)
})
plot(ns,by.n*ns)

## asymptotic behavior of first four components [like 1/2-1/n]
n <- 40
B <- 1e5
ns <- round(seq(5,100,len=20))
by.n <- sapply(ns, function(n) {
    firstfour <- replicate(1e2, {
        ## s <- runif(n)
        s <- rexp(n)
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        ## eig$vec[,1]; eig$vec[,n]
        eig$vec[1:4,1]
    })
    mean(abs(firstfour))-1/2
})
plot(ns,by.n*ns^(1))

## over four domination
B <- 1e4
ns <- round(seq(5,100,len=20))
overfour <- replicate(B, {
    n <- round(runif(1,5,200))

    n <- 20
    s <- runif(n,0,1)
    s[1] <- s[2]+.2
     s[3] <- s[4]+.05
    ## s <- rexp(n)
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    eig$vec[n,1]/eig$vec[n,n]
    eig$vec[,c(1,n)]
    
    ratios <- eig$vec[5:n,1]/eig$vec[5:n,n]
    if(abs(max(ratios)-min(ratios))>1e-8)cat(abs(max(ratios)-min(ratios)),n,'\n')
    ## stopifnot(abs(max(ratios)-min(ratios))<1e-7)
    if(abs(min(abs(ratios))-1)<1e-4){browser()}
    ratios[1]
    ## c(sum(eig$vec[5:n,1] - eig$vec[5:n,n] > 0),        sum(eig$vec[5:n,1] - eig$vec[5:n,n] < 0))
})

summary(abs(overfour)) ## ratios >= 1

## sum(apply(overfour,2,min)!=0)
## mean(overfour[1,]==0) + c(-1,1)/2/sqrt(B)

## over four domination--ratio by diffs
## [ abs((s3-s4)/(s1-s2)) = factor
## for abs(factor) > 1:
## abs(ratio) = 1+2/(factor-1) ]
## abs(ev[>4,1]/ev[>4,n]) = [(s3-s4)+(s1-s2)]/[(s3-s4)-(s1-s2)]
## [ for 0 < factor < 1:
## abs(ratio) = 1+2*factor/(1-factor) ]
## abs(ev[>4,1]/ev[>4,n]) = -[(s3-s4)+(s1-s2)]/[(s3-s4)-(s1-s2)]
## ]
n <- 8
factors <- seq(2,10,len=100)
ratios <- sapply(factors, function(factor) {
    s <- runif(n,0,1)
    s <- rexp(n)
    s[1] <- s[2]+.3
    s[3] <- s[4]+.3*factor
    ## s[4] <- s[3]+.3*factor
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    eig$vec[n,1]/eig$vec[n,n]
    ## ((s[3]-s[4])+(s[1]-s[2]))/((s[3]-s[4])-(s[1]-s[2]))
    ## (abs(s[3]-s[4])+abs(s[1]-s[2]))/(abs(s[3]-s[4])-abs(s[1]-s[2]))
})
plot(factors,abs(ratios))
curve(1+2/(x-1),add=TRUE)
## plot(factors,(abs(ratios)))
## curve(1+2*x/(1-x),add=T)

## abs(ratio)=(abs(s3-s4)+abs(s1-s2))/abs(abs(s3-s4)-abs(s1-s2)) > 1
n <- 20
ratios <- replicate(1e3, {
    factor <- runif(1,-10,10)
    s <- runif(n,0,1)
    s <- rexp(n)
    s[1] <- s[2]+.3
    s[3] <- s[4]+.3*factor
    ## s[4] <- s[3]+.3*factor
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    c(obs=eig$vec[n,1]/eig$vec[n,n],
      formula=(abs(s[3]-s[4])+abs(s[1]-s[2]))/(abs(s[3]-s[4])-abs(s[1]-s[2]))
      )
    ## c(obs=eig$vec[n,1]/eig$vec[n,n],
    ##   formula=(abs(s[3]-s[4])+abs(s[1]-s[2]))/(max((s[3]-s[4]),(s[1]-s[2]))-min((s[3]-s[4]),(s[1]-s[2])))
    ##    )
})
plot(ratios['obs',],ratios['formula',],asp=1);abline(0,1);abline(0,-1)

## sign seems random
s5s <- seq(.001,20,len=200)
signs <- sapply(s5s, function(s5) {
    s <- c(0,1,0,2,s5)
    ## s[4] <- s[3]+.3*factor
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    (eig$vec[n,1]/eig$vec[n,n])
})
plot(s5s,signs)


## sign of ratio
B <- 3e3
n <- 20
evs <- replicate(B, {
    factor <- runif(1,0,10)
    s <- runif(n,0,1)
    s <- rexp(n)
    ## s[4] <- s[3]+.3*factor
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    (eig$vec[n,1]/eig$vec[n,n])   
    d1 <- s[1]-s[2]; d2 <- s[3]-s[4]
    ## c(obs=eig$vec[n,1]/eig$vec[n,n],
    ##   formula=-(abs(d1)+abs(d2))/(abs(d1)-abs(d2)) 
    ##   )
    ## c(sign(eig$vec[n,1]/eig$vec[n,n]),paste(rank(s),collapse=''))
    ## c(sign(eig$vec[n,1]/eig$vec[n,n]),paste(sign(outer(s,s,`-`)),collapse=''})
    ## cbind(eig$vec[,1],eig$vec[,n])    
    ## eig$vec[,n]  <- eig$vec[,n]* sign(eig$vec[n,n]/eig$vec[n,1])
    cbind(eig$vec[,n] , eig$vec[,1])
})
## split(ratios[1,],ratios[2,])
## plot(ratios['obs',],ratios['formula',],asp=1);abline(0,1);abline(0,-1)
evs <- lapply(1:B, function(j){
    ## ratio.sign <- unique(sign(evs[5:n,1,j]*evs[5:n,2,j]))
    ## stopifnot(length(ratio.sign)==1)
    cbind(evs[,1,j]*sign(evs[5,1,j]),evs[,2,j]*sign(evs[5,2,j]))
})
## evs <- simplify2array(evs)    
signs <- simplify2array(lapply(evs,function(m)sign(m[1:4,1]*m[1:4,2])))
stopifnot(sum(signs[1,]!=signs[2,])==0)
stopifnot(sum(signs[3,]!=signs[4,])==0)
stopifnot(sum(signs[1,]==signs[3,])==0)
signs <- signs[1,]
mean(signs==1)





dd

## sign of ratio [somehow related to components > 4]
require(parallel)
n <- 5
B <- 1e3
factor <- 3
smat <- matrix(runif(n*B,0,1),nrow=n)
smat[2,] <- smat[1,]+.4
smat[4,] <- smat[3,]+.4*factor
signs <- mclapply(1:B, mc.cores=detectCores()-3,FUN=function(j) {
    s <- smat[,j]
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    sign(eig$vec[n,1]/eig$vec[n,n])
})
signs <- simplify2array(signs)
mean(signs>0)
boxplot(smat[5,] ~ signs)
d1 <- (smat[2,]-smat[1,]); d2 <- (smat[4,]-smat[3,])
boxplot(smat[5,]* (d1*d2)~ signs)

## abs(s[5]/ev[5]) is a parabola in s5, 2/(d1+d2)*s5^2 + ...
d1 <- 3
d2 <- 2
s <- c(1,1+d1,1,1+d2)
## 5th component
## horiz asymptote of abs(ev[5])*s[5] at [|s1-s2| + |s3-s4|]/2
require(parallel)
n <- 5
B <- 1e3
## s <- runif(5)
## s <- rexp(5)
## s <- (1:4)
## s <- c(1,3,4,5)
s5 <- seq(-2,2,len=1e3)
fifth <- mclapply(s5, mc.cores=detectCores()-3,FUN=function(x) {
    s[5] <- x
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    eig$vec[5,1]
})
fifth <- simplify2array(fifth)
## plot(s5,abs(fifth*s5^2),type='l',asp=1)
plot(s5,abs(s5/fifth),type='l',asp=1)
## (abs(s[3]-s[4])+abs(s[1]-s[2]))/(abs(s[3]-s[4])-abs(s[1]-s[2]))
(( lm0 <- lm(I(abs(s5/fifth)) ~ I(s5^2)) ))
curve(coef(lm0)[1]+x^2*coef(lm0)[2],add=TRUE,col=2)
curve(coef(lm0)[1]+x^2*2/(abs(s[2]-s[1])+abs(s[4]-s[3])),add=TRUE,col=3)
2/(d1+d2)

summary(abs(s5/fifth) - s5^2*2/(abs(s[2]-s[1])+abs(s[4]-s[3])))

dd

## looking for leading coeff [ans: 2/(d1+d2)]
d1 <- runif(1,0,10)
d2s <- seq(1,20,len=100)
d2s <- 1:50
by.d2 <- sapply(d2s, function(d2) {
    s <- c(1,1+d1,2,2+d2)
    s <- c(2,2+d1,3,3+d2)
    s <- runif(4)
    ## s <- c(1,NA,runif(1),NA)
    s[2] <- s[1]+d1; s[4] <- s[3]+d2
    require(parallel)
    n <- 5
    B <- 20
    s5 <- seq(-2,2,len=1e1)
    fifth <- mclapply(s5, mc.cores=detectCores()-3,FUN=function(x) {
        s[5] <- x
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        eig$vec[5,1]
    })
    fifth <- simplify2array(fifth)
    coef(lm(I(abs(1/fifth*s5)) ~ I(s5^2)))[2]
})
plot(d2s,by.d2,asp=1,xlim=c(1,2),type='l')
## curve(1/x^(2/3)*1/d1^(2/3),add=T,col=2)
## curve(2*(d1+1)^(-2/3)*(x+1)^(-2/3),add=T,col=2)
curve(2/(d1+x),add=T,col=2)

dd

## abs(s[5]/ev[5]) is a parabola in s5, look for constant term
d1 <- 1#runif(1,0,10)
d2s <- seq(.01,100,len=50)
## d2s <- 1:70
s <- c(1,1+d1,1,1)
## s <- runif(4)
## s <- rexp(4)
s[2] <- s[1]+d1
by.d2 <- sapply(d2s, function(d2) {
    ## s <- c(1,30,1,30)

    s[4] <- s[3]+d2
    require(parallel)
    n <- 5
    B <- 1e1
    s5 <- seq(-2,2,len=B)
    ## s5 <- 1:10
    fifth <- mclapply(s5, mc.cores=detectCores()-3,FUN=function(x) {
        s[5] <- x
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        eig$vec[5,1]
    })
    fifth <- simplify2array(fifth)
    unname(coef(lm(I(abs(s5/fifth)) ~ I(s5^2)))[1])
    abs(s5/fifth) - 2*s5^2/(d1+d2)
    2*sum(s[1:4]^2)/(d1+d2)

})

plot(d2s,by.d2,type='l')
curve(2*(sum(s[1:3]^2)+(s[3]+x)^2)/(d1+x),add=TRUE,col=2)
## plot(d2s,(d1+d2s)/by.d2,type='l')

dd

## ## constant term really is constant wrt s[5]?
## ## s <- c(1,30,1,30)
## s <- c(1,1+runif(1),1,1+runif(1))
## s <- s/sqrt(sum(s^2))
## require(parallel)
## B <- 1e2
## remainder <- replicate(B, {
## s[5] <- runif(1)
## Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
## Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
## SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
## eig <- eigen(SSmat)
## ## eig$vec[5,1]
## abs(s[5]/eig$vec[5,1]) - s[5]^2*2/(abs(s[2]-s[1])+abs(s[3]-s[4]))
## })
## summary(remainder)
## 2*sum(s[1:4]^2)/(abs(s[2]-s[1])+abs(s[3]-s[4]))

## d1 <- .1;d2 <- .3;s1 <- .13
## d1 <- runif(1);d2 <- runif(1);
## uniroot(function(s2)s2^2+(s2+d2)^2+s1^2+(s1+d1)^2-1,c(.1,5))$root
## s2 

## d1 <- .33967; d2 <- .09426708; s1 <-  
## s <- c(0.03407867,0.03407867+d1,0.6066143,0.6066143+d2)
## s <- c(0.4704088,0.4704088+d1,0.1958315,0.1958315+d2)


## s <- c(1,1+runif(1),1,1+runif(1))
## s <- c(0,1/sqrt(2),0,1/sqrt(2))
## Q <- qr.Q(qr(matrix(rnorm(16),nrow=4)))
## s <- s%*%Q
require(parallel)
n <- 5
B <- 1e3
## s.norm <- runif(1)
rem <- replicate(B, {
    s <- rexp(4)
    ## s <- c(1,2,1,2)
    s <- s/sqrt(sum(s^2))
    s[5] <- 1#runif(1)
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    ## eig$vec[5,1]
    c(d1=abs(s[2]-s[1]),d2=abs(s[3]-s[4]),normsqr=sum(s[1:4]^2),const=abs(s[5]/eig$vec[5,1]) - 2*s[5]^2/(abs(s[2]-s[1])+abs(s[3]-s[4])))
})
formula = rem['normsqr',]*2/(rem['d1',]+rem['d2',])
formula = rem['normsqr',]*2/(rem['d1',]+rem['d2',])/(1+rem['d1',]^2)
## plot(rem['d1',],rem['const',])
## require(rgl)
## plot3d(rem['d1',],rem['d2',],rem['const',])
plot(rem['const',],rem['normsqr',]*2/(rem['d1',]+rem['d2',]))
abline(0,1)
max(abs(rem['const',] - formula))

plot(rem['normsqr',],abs(rem['const',] - formula))
plot(rem['d2',],formula/abs(rem['const',])-1,asp=1)

rem0 <- rem[,rem['d2',]>.1]
formula0 <- formula[rem['d2',]>.1]
plot(rem0['d2',],(formula0/abs(rem0['const',])-1)/rem0['d2',]^(1/2),xlim=c(.2,1))


require(rgl)
grid = seq(0,1,len=100)
xy = expand.grid(grid,grid)
plot3d(rem['d1',],rem['d2',],abs(rem['const',]-formula))
points3d(xy[,1],xy[,2],xy[,1]*xy[,2],col='red')

plot3d(rem['d1',],rem['d2',],formula/abs(rem['const',])-1)
points3d(xy[,1],xy[,2],(xy[,2]^2+xy[,1]^2)*.33,col=2)

plot3d(rem['d1',],rem['d2',],formula/abs(rem['const',])-1)
points3d(xy[,1],xy[,2],((1+xy[,1]^2)*(1+xy[,2]^2)-1)*.33,col=rgb(1,0,0,.3))

## discrepancy
require(parallel)
n <- 5
B <- 1e3
## s.norm <- runif(1)
rem <- replicate(B, {
    s <- rexp(4)
    sizes <- 1:100
    deltas <- sapply(sizes, function(size) {
        s <- c(1.4,2,1,8)
        s <- s/sqrt(sum(s^2))*size
    s[5] <- runif(1)
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    ## eig$vec[5,1]
    ## c(d1=abs(s[2]-s[1]),d2=abs(s[3]-s[4]),normsqr=sum(s[1:4]^2),const=abs(s[5]/eig$vec[5,1]) - 2*s[5]^2/(abs(s[2]-s[1])+abs(s[3]-s[4])))
    const=abs(s[5]/eig$vec[5,1]) - 2*s[5]^2/(abs(s[2]-s[1])+abs(s[3]-s[4]))
    formula <- 2*sum(s[1:4]^2)/(abs(s[2]-s[1])+abs(s[3]-s[4]))
    delta <- abs(const-formula)
        c(true=const,formula=formula,delta=delta)
    })    
    plot(sizes,deltas['delta',]/deltas['delta',1]);abline(0,1)
    ## lm(deltas['delta',] ~ sizes)

})
## plot(rem['d1',],rem['const',])
## require(rgl)
## plot3d(rem['d1',],rem['d2',],rem['const',])
plot(rem['const',],rem['normsqr',]*2/(rem['d1',]+rem['d2',]))
abline(0,1)
max(abs(rem['const',]-rem['normsqr',]*2/(rem['d1',]+rem['d2',])))



## abs(s5) \approx (|s1-s2|+|s3-s4|)/(sum_{j=1}^5 s_j^2) * s_5/2
require(parallel)
n <- 5
B <- 1e3
fifth <- replicate(B, {
    s <- rexp(n)
    s <- runif(n,0,.1)
    ## s <- runif(n,3,5) ## formula tight here
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    ## c(obs=eig$vec[5,1],formula=(abs(s[2]-s[1])+abs(s[3]-s[4]))/(s[5]+2*sum(s[1:4]^2)/s[5]))
    ## d1 <- abs(s[2]-s[1]);d2 <- abs(s[3]-s[4])
    c(obs=eig$vec[5,1],formula=(abs(s[2]-s[1])+abs(s[3]-s[4]))/sum(s^2) * s[5]/2)
})
plot(abs(fifth['obs',]),fifth['formula',]); abline(0,1)

dd

## n>5 case
require(parallel)
n <- 50
B <- 1e3
ns <- 5:50
## by.n <- sapply(ns, function(n) {
fifth <- replicate(B, {
    s <- rexp(n)
    ## s <- runif(n,3,5) ## formula tight here
    ## s <- runif(n,0,1)
    Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
    SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
    eig <- eigen(SSmat)
    ## c(obs=eig$vec[5,1],formula=(abs(s[2]-s[1])+abs(s[3]-s[4]))/(s[5]+2*sum(s[1:4]^2)/s[5]))
    ## d1 <- abs(s[2]-s[1]);d2 <- abs(s[3]-s[4])
    c(obs=eig$vec[6,1],formula=(abs(s[2]-s[1])+abs(s[3]-s[4]))/sum(s^2) * s[6]/2)
})
## coef(lm(fifth['formula',] ~ abs(fifth['obs',])))
## })
## plot(ns,by.n[2,])
plot(abs(fifth['obs',]),fifth['formula',]); abline(0,1)







dd



## 46 checking formulas p.31, asymptotic variance taking into account O(1/n) factor in the variances
n <- 10
s <- runif(n)
sigma <- 1/s
sigma.bar <- sqrt(sigma^2-1/sum(1/sigma^2))
cov <- sum(s^2)/sum(s[5:n]^2)*(-2+sum(s[1:4]^2)/sum(s^2))
v <- c(sum(s^2)^2/sum(s[5:n]^2)*(sigma[1]^2*sigma.bar[2]^2+sigma.bar[1]^2*sigma[2]^2)/(sigma.bar[2]-sigma.bar[1])^2 + cov, sum(s^2)^2/sum(s[5:n]^2)*(sigma[3]^2*sigma.bar[4]^2+sigma.bar[3]^2*sigma[4]^2)/(sigma.bar[4]-sigma.bar[3])^2 + cov)
rho <- cov/sqrt(prod(v))
beta <- 1/sqrt(1-rho^2) * (rho/sqrt(v) - 1/sqrt(rev(v))); gamma <- 1/sqrt(v)
## beta/gamma
## sqrt(prod(v)/(prod(v)-cov^2))*(cov/sqrt(prod(v)) - sqrt(v/rev(v)))
## cov-v
## sum(s^2)^2/sum(s[5:n]^2)*(sigma[1]^2*sigma.bar[2]^2+sigma.bar[1]^2*sigma[2]^2)/(sigma.bar[2]-sigma.bar[1])^2
beta/gamma/sqrt(1+beta^2+gamma^2)
(cov-v)/sqrt(prod(v)-cov^2-2*cov+sum(v))
## u <- 1/sum(s[5:n]^2)
u.bar <- function(u)1/(1/u+sum(s[1:4]^2))
## f <- function(u)c(-(s[1]^2+s[2]^2)/sum(s^2),-(s[3]^2+s[4]^2)/sum(s^2))
## g <- c(-2*s[1]^2*s[2]^2/sum(s^2) + 2*s[1]*s[2]*(1-sqrt(1-(s[1]^2+s[2]^2)/sum(s^2)+s[1]^2*s[2]^2/sum(s^2)^2)),   -2*s[3]^2*s[4]^2/sum(s^2) + 2*s[3]*s[4]*(1-sqrt(1-(s[3]^2+s[4]^2)/sum(s^2)+s[3]^2*s[4]^2/sum(s^2)^2)) )
f <- function(u)c(-(s[1]^2+s[2]^2),-(s[3]^2+s[4]^2))*u.bar(u)
g <- function(u)c(-2*s[1]^2*s[2]^2*u.bar(u) + 2*s[1]*s[2]*(1-sqrt(1-(s[1]^2+s[2]^2)*u.bar(u)+s[1]^2*s[2]^2*u.bar(u)^2)),   -2*s[3]^2*s[4]^2*u.bar(u) + 2*s[3]*s[4]*(1-sqrt(1-(s[3]^2+s[4]^2)*u.bar(u)+s[3]^2*s[4]^2*u.bar(u)^2)) )
delta <- c(s[1]-s[2],s[3]-s[4])
## (2+f(u)) / (delta^2+g(u)) / sqrt((2+f(u))/(delta^2+g(u))*(2+rev(f(u)))/(rev(delta)^2+rev(g(u))) - 1/sum(s^2)*((2+f(u))/(delta^2+g(u)) + (2+rev(f(u)))/(rev(delta^2)+rev(g(u))))  )
## fn <- Vectorize(function(u)((2+f(u))/(delta^2+g(u)))[1])
## curve(fn,-1,3)
## abline(a=fn(0),b=-(s[1]+s[2])^2/(s[1]-s[2])^2)
## fn <- Vectorize(function(u)(-(s[1]^2+s[2]^2)*u.bar(u))[1])
## curve(fn,-.2,3)
## abline(a=fn(0),b=-(s[1]^2+s[2]^2))
## u <- 1/sum(s[5:n]^2)
u.bar <- function(u)1/(1/u+sum(s[1:4]^2))
ratio <- function(u)(2+f(u))/(delta^2+g(u))
cov <- function(u)1/u.bar(u)*u*(-2+sum(s[1:4]^2)*u.bar(u))
v <- function(u) 1/u.bar(u)^2*u*ratio(u) + cov(u)
rho <- function(u)cov(u)/sqrt(prod(v(u)))
beta <- function(u)1/sqrt(1-rho(u)^2) * (rho(u)/sqrt(v(u)) - 1/sqrt(rev(v(u)))); gamma <- function(u)1/sqrt(v(u))
ratio.prime <- -c((s[1]+s[2])^2,(s[3]+s[4])^2)/delta^2
fn <- Vectorize(function(u) atan(beta(u)/gamma(u)/sqrt(1+beta(u)^2+gamma(u)^2))[1] )
curve(fn)
fn0 <- Vectorize(function(u) atan(-ratio(u)/sqrt(prod(ratio(u)) - u.bar(u)*sum(ratio(u))  ))[1])
curve(fn0,col='red',add=TRUE,lty=3)
slope <- -(1+ratio(0)[1]/ratio(0)[2])^(-1)*(ratio.prime[1]/sqrt(prod(ratio(0))) - ratio(0)[1]/2*prod(ratio(0))^(-3/2)*(sum(ratio.prime*rev(ratio(0))) - sum(ratio(0))))
abline(a=fn0(0),b=slope)
slope
-(1+delta[2]^2/delta[1]^2)^(-1)*(-(s[1]+s[2])^2/delta[1]^2*prod(abs(delta))/2 - 1/2*(2/delta[1]^2)^(-1/2)*(2/delta[2]^2)^(-3/2)*((-2*(s[3]+s[4])^2-2*(s[1]+s[2])^2)/prod(delta^2)-2*sum(delta^2)/prod(delta^2)))
-delta[1]^2/sum(delta^2)/2*delta[2]/delta[1]*(s[3]^2+s[4]^2-2*s[1]*s[2])
fn <- Vectorize(function(u) sum(atan(beta(u)/gamma(u)/sqrt(1+beta(u)^2+gamma(u)^2)) ))
curve(fn)
slope <- prod(delta)/2
abline(a=fn(.001),b=slope)





## 47 stochastic ordering for means as function of n
n <- 1
B <- 1e5
th1 <- runif(1); th2 <- runif(1)
th1 <- 0.3871228; th2 <- 0.2064894
rS <- function(B)rbeta(B,th1,th2)
rS <- function(B)rbinom(B,1,.1)-1
xn <- rowMeans(matrix(rS(n*B),ncol=n))
xn1 <- rowMeans(matrix(rS((n+1)*B),ncol=n+1))
plot(ecdf(xn))
lines(ecdf(xn1),col='red')
abline(v=mean(xn1))
## counterexample to P(\bar{X}_{n}<y)>P(\bar{X}_{n+1}<y) for y<E(X)
## > th1
## [1] 0.3871228
## > th2
## [1] 0.2064894
## > 



## concentration for S^2 density with infinite slope at 0
n <- 1e5
hist(runif(n)^(2/3),prob=TRUE)
curve(3/2*x^(1/2),add=TRUE,col='red')
rS <- function(n)runif(n)^(2/3)
hist(replicate(n,mean(rS(10))))


## 48 convergence of bias to asy bias

## 48a checking rho by simulations
n <- 1e4
mean(abs(rexp(n)-rexp(n)))^2 / mean(rexp(n)^2)


## 48b 
source('misc.R')
require(parallel)
require(EnvStats)
ns <- round(seq(1e1,1e3,len=10))
reps <- 5e2
rSs <- c(uniform=function(n)runif(n,0,1), exponential=rexp, gamma=function(n)rgamma(n,shape=.54), beta=function(n)rbeta(n,.15,.39), pareto=function(n)rpareto(n,location=1,shape=2.52687))
rhos <- c(uniform=1/3,exponential=1/2,gamma=.56,beta=.66,pareto=.1390052)
distributions <- names(rSs)
distributions <- c('uniform','exponential')
by.distr <- lapply(distributions, function(distribution) {
    rS <- rSs[distribution]; rho <- rhos[distribution]
    by.n <- sapply(ns, FUN=function(n) {
        tau.stats <- mclapply(1:reps, mc.cores=detectCores()-2, FUN=function(jj) {
            s <- rS(n)
            z <- rnorm(n)
            ## y <- z/s
            ## theta.fe <- sum(y*s^2)/sum(s^2)
            tau(z,s)
        })
        tau.stats <- simplify2array(tau.stats)
        var(tau.stats)
    })
    cov.1234 <- -rho/pi/4/ns
    observed <- ns*by.n-4/9
    asymptotic <- 4*(ns-2)*(ns-3)/(ns-1)*cov.1234
    ## 4*(ns-2)*(ns-3)/(ns-1)*cov.1234
    ## plot(ns,observed,type='l',main='bias')
    ## lines(ns,asymptotic,lty=2)
    ## legend('topright',lty=1:2,legend=c('observed','formula'))
    rbind(observed=observed,asymptotic=asymptotic)
})
diff <- sapply(by.distr, function(m)abs(m['observed',]-m['asymptotic',]))
matplot(ns,diff,type='l',lty=1:length(distributions),col=1,xlab=n,ylab='asy -observed variance times n')
legend('topright',lty=1:length(distributions),legend=distributions)

dd


## 48a rejection rate varying true mu, verifying no change
require(parallel)
source('misc.R')
rSs <- c(uniform=function(n)runif(n,0,1), exponential=rexp, gamma=function(n)rgamma(n,shape=.54), beta=function(n)rbeta(n,.15,.39), pareto=function(n)rpareto(n,location=1,shape=2.52687))
rSs <- c(uniform=function(n)runif(n,0,1), exponential=rexp, gamma=function(n)rgamma(n,shape=.54), beta=function(n)rbeta(n,.15,.39))
ns <- c(30,1e2,200)
alpha <- .05
q <- qnorm(1-alpha/2)
reps <- 1e2
mus <- seq(0,3,len=10)
by.n <- lapply(ns, function(n) {
    by.distr <- lapply(rSs, function(rS) {
        by.mu <- sapply(mus, function(mu) {
            rejected <- mclapply(1:reps, mc.cores=detectCores()-2, FUN=function(jj) {
                s <- rS(n)
                y <- rnorm(n,mu,sd=1/s)
                z <- (y-mu)*s
                tau.stat <- tau(z,s)
                sqrt(9/4*n)*abs(tau.stat) > q
            })
            rejected <- simplify2array(rejected)
            mean(rejected)    
        })
        ## plot(mus,by.mu)
    })
})
op <- par(mfrow=c(1,3))
lapply(1:length(by.n), function(jj) {
    by.distr <- simplify2array(by.n[[jj]])
    matplot(mus,by.distr,type='l',lty=1:length(rSs),col=1,xlab=expression(mu),ylab='rejection rate',main=paste0('n=',ns[jj]))
    ## abline(h=.05)
    legend('bottomright',lty=1:length(rSs),legend=names(rSs))
})
par(op)







## 49 begg '94 simulation replication
require(parallel)
require(reshape2)
tau <- function(y,v) {
    ## y <- z/s
    theta.fe <- sum(y/v)/sum(1/v)
    ## list(tau=cor((y-theta.fe)/sqrt(v-1/sum(1/v)),v,method='kendall'),theta.fe=theta.fe)
    cor((y-theta.fe)/sqrt(v-1/sum(1/v)),v,method='kendall')
}
w <- function(t,v,strong=FALSE){
    a <- strong*1.5+(!strong)*3
    b <- 4
    p <- pnorm(-t/sqrt(v))
    exp(-b*p^a)
}
## curve(w(-x,1),-2,2)
## curve(w(-x,1,strong=TRUE),add=TRUE,lty=2)
B <- 3e3
q <- qnorm(1-.05/2)
deltas <- seq(0,3,by=.5)
param.configs <- expand.grid(strong.select=c(TRUE,FALSE),vs=list(small=c(.5,1,2),large=c(.1,1,10)),delta=deltas,ns=list(small=c(8,9,8),large=c(25,25,25)))
reject <- mclapply(1:B, mc.cores=detectCores()-4, FUN= function(jj){
## reject <- lapply(1:B, FUN= function(jj){
    cat('.')
    ## reject <- lapply(1:B, FUN= function(jj) {
    apply(param.configs,1,function(params) {
        ## with(params, {        
        delta <- params$delta; vs <- unlist(params$vs); strong.select=params$strong.select; ns <- unlist(params$ns)
        total.sampled <- 0
        data <- lapply(1:3, function(j) {
            v <- vs[j];  n <- ns[j]
            t <- numeric()
            v.jittered <- numeric()
            while(length(t)<n) {
                total.sampled <<- total.sampled+1
                v.try <- jitter(v)
                t.try <- rnorm(1,mean=delta,sd=sqrt(v.try))
                if(rbinom(1,1,prob=w(t.try,v.try,strong=strong.select))) {
                    t <- c(t,t.try)
                    v.jittered <- c(v.jittered,v.try)
                }
            }
            rbind(t=t,v=v.jittered)
        })
        data <- do.call(cbind,data)
        t <- data['t',]; v <- data['v',]
         ## theta.fe <- sum(t/v)/sum(1/v);plot(v,(t-theta.fe)/sqrt(v))
        c(reject=sqrt(9*length(t)/4)*abs(tau(t,v)) > q,total.sampled=total.sampled)
        ## })
    })
})
total.sampled <- rowMeans(sapply(reject,function(m)m['total.sampled',]))
sample.rate <- sapply(param.configs$ns,sum) / total.sampled
power <- rowMeans(sapply(reject,function(m)m['reject',]))
tables <- lapply(list(sample.rate=sample.rate,power=power), function(data){
    out <- cbind(param.configs,value=data)
    out$vs <- names(out$vs)
    out$ns <- names(out$ns)
    ## out <- melt(out,id.vars=c('ID','power'))
    out.smallstudy <- out[out$ns=='small',]
    out.smallstudy <- acast(out.smallstudy, strong.select ~ vs ~ delta)
})
power <- tables[['power']]; sample.rate <- tables[['sample.rate']]
attr(sample.rate,'dimnames') <- attr(power,'dimnames')  <-  structure(attr(power,'dimnames'),names=c('strong.select','v.range','delta'))
power
## save.image('210312.RData')
## load('210312.RData')



## 50 power analysis
require(parallel)
require(reshape2)
require(EnvStats)
tau <- function(y,v) {
    ## y <- z/s
    theta.fe <- sum(y/v)/sum(1/v)
    ## list(tau=cor((y-theta.fe)/sqrt(v-1/sum(1/v)),v,method='kendall'),theta.fe=theta.fe)
    cor((y-theta.fe)/sqrt(v-1/sum(1/v)),v,method='kendall')
}
w <- function(t,v,a=1,b=4){
    p <- pnorm(-t/sqrt(v))
    exp(-b*p^a)
}
## curve(w(-x,1),-2,2)
## curve(w(-x,1,strong=TRUE),add=TRUE,lty=2)
B <- 1e3
q <- qnorm(1-.05/2)
ns <- c(30,1e2,200)
rSs <- c(uniform=function(n)runif(n,0,1), exponential=rexp, gamma=function(n)rgamma(n,shape=.54), beta=function(n)rbeta(n,.15,.39), pareto=function(n)rpareto(n,location=1,shape=2.52687))
select.strengths <- seq(0,5,len=10)
by.n <- mclapply(ns, mc.cores=detectCores()-2, FUN=function(n) {
    by.distr <- lapply(rSs, function(rS) {
        by.select <- sapply(select.strengths, function(select.strength) {
            rowMeans(replicate(B, {
                total.sampled <- 0
                t <- v <- numeric()
                while(length(t)<n) {
                    total.sampled <- total.sampled+1
                    v.try <- 1/rS(1)^2
                    t.try <- rnorm(1,sd=sqrt(v.try))
                    if(rbinom(1,1,prob=w(t.try,v.try,b=select.strength))) {
                        t <- c(t,t.try)
                        v <- c(v,v.try)
                    }
                }
                ## theta.fe <- sum(t/v)/sum(1/v);plot(v,(t-theta.fe)/sqrt(v))
                c(reject=sqrt(9*length(t)/4)*abs(tau(t,v)) > q,total.sampled=total.sampled)
            }))
        })
    })
})
op <- par(mfrow=c(1,3))
for(jj in 1:length(by.n)) {
    by.distr <- by.n[[jj]]#simplify2array(by.n[[jj]])
    power.curves <- sapply(by.distr,function(m)m['reject',])
    matplot(select.strengths,power.curves,type='l',lty=1:length(rSs),col=1,xlab='selection strength',ylab='rejection rate')
    ## abline(h=.05)
    legend('bottomright',lty=1:length(rSs),legend=colnames(power.curves))
}
par(op)
## save.image('210313a.RData')



## 50a power analysis, with bias corrected terms
require(parallel)
require(reshape2)
require(EnvStats)
tau <- function(y,v) {
    ## y <- z/s
    theta.fe <- sum(y/v)/sum(1/v)
    ## list(tau=cor((y-theta.fe)/sqrt(v-1/sum(1/v)),v,method='kendall'),theta.fe=theta.fe)
    cor((y-theta.fe)/sqrt(v-1/sum(1/v)),v,method='kendall')
}
w <- function(t,v,a=1,b=4){
    p <- pnorm(-t/sqrt(v))
    exp(-b*p^a)
}
## curve(w(-x,1),-2,2)
## curve(w(-x,1,strong=TRUE),add=TRUE,lty=2)
B <- 1e3
q <- qnorm(1-.05/2)
ns <- c(25,75,150)
## ns <- round(seq(10,1.5e2,len=10))
rSs <- c(uniform=function(n)runif(n,0,1), exponential=rexp, gamma=function(n)rgamma(n,shape=.54), beta=function(n)rbeta(n,.15,.39), pareto=function(n)rpareto(n,location=1,shape=2.52687))
biases <- c(uniform=.11, beta=.21, exponential=.16, gamma=.18, pareto=.04)
distributions <- names(rSs)
## distributions <- c('beta','uniform')
names(distributions) <- distributions
select.strengths <- seq(0,5,len=10)
by.n <- mclapply(ns, mc.cores=detectCores()-2, FUN=function(n) {
    by.distr <- lapply(distributions, function(distr) {
        rS <- rSs[[distr]]; bias <- biases[distr]
        by.select <- sapply(select.strengths, function(select.strength) {
            rowMeans(replicate(B, {
                total.sampled <- 0
                t <- v <- numeric()
                while(length(t)<n) {
                    total.sampled <- total.sampled+1
                    v.try <- 1/rS(1)^2
                    t.try <- rnorm(1,sd=sqrt(v.try))
                    if(rbinom(1,1,prob=w(t.try,v.try,b=select.strength))) {
                        t <- c(t,t.try)
                        v <- c(v,v.try)
                    }
                }
                s <- sqrt(1/v)
                s.md.hat <- mean(abs(outer(s,s,`-`)))
                s2.hat <- mean(s^2)
                bias.hat <- s.md.hat^2/s2.hat/pi
                ## theta.fe <- sum(t/v)/sum(1/v);plot(v,(t-theta.fe)/sqrt(v))
                c(reject=unname(sqrt(n/(4/9))*abs(tau(t,v)) > q), reject.debiased=unname(sqrt(n/(4/9-bias))*abs(tau(t,v)) > q),reject.debiased.hat=unname(sqrt(n/(4/9-bias.hat))*abs(tau(t,v)) > q),total.sampled=unname(total.sampled))
            }))
        })
    })
})
op <- par(mfrow=c(1,3))
for(jj in 1:length(by.n)) {
    by.distr <- by.n[[jj]]#simplify2array(by.n[[jj]])
    power.curves <- sapply(by.distr,function(m)m['reject',])
    matplot(select.strengths,power.curves,type='l',lty=1:length(rSs),col=1,xlab='selection strength',ylab='rejection rate')
    power.curves.debiased <- sapply(by.distr,function(m)m['reject.debiased',])
    matplot(select.strengths,power.curves.debiased,type='l',lty=1:length(rSs),col=2,xlab='selection strength',ylab='rejection rate',add=TRUE)
    power.curves.debiased.hat <- sapply(by.distr,function(m)m['reject.debiased.hat',])
    matplot(select.strengths,power.curves.debiased.hat,type='l',lty=1:length(rSs),col=3,xlab='selection strength',ylab='rejection rate',add=TRUE)
    abline(h=.05)
    legend('bottomright',lty=1:length(rSs),legend=distributions)
}
par(op)
## save.image('210313c.RData')
## load('210313c.RData')


op <- par(mfrow=c(length(distributions),length(ns)))
for(j in 1:length(distributions)) {
    power <- lapply(by.n,function(lst)lst[[distributions[j]]])
    for(k in 1:length(ns)) {
        ## plot(select.strengths,power[[k]]['reject',],ty)
        percent.filtered <- 1-ns[k]/power[[k]][4,]
        matplot(percent.filtered,t(power[[k]][-4,]),type='l',col=1,lty=1:3,ylab='power',xlab='proportion of studies not selected')
        }
}
par(op)

require(ggplot2)
by.distr <- lapply(distributions,function(distribution)simplify2array(lapply(by.n,function(lst)lst[[distribution]])))
by.distr <- lapply(by.distr,function(lst){
    dimnames(lst) <- list(estimator=c('standard','debiased','debiased.hat','total.sampled'),select.strength=round(select.strengths,2),n=ns)
    ## for(j in 1:dim(lst)[3]) lst[,,j]
    lst
})
by.distr <- simplify2array(by.distr)
names(dimnames(by.distr))[4] <- 'distribution'
by.distr <- melt(by.distr)

dat <- merge(by.distr[by.distr$estimator!='total.sampled',],by.distr[by.distr$estimator=='total.sampled',colnames(by.distr)!='estimator'],by=c('select.strength','n','distribution'))
colnames(dat)[colnames(dat)=='value.x'] <- 'power'
colnames(dat)[colnames(dat)=='value.y'] <- 'total.sampled'
dat$select.strength <- round(dat$n / dat$total.sampled,2)
plt <- ggplot(dat[dat$distribution!='pareto',],aes(x=select.strength,y=power,group=estimator,linetype=estimator))+geom_line()+facet_wrap(. ~ distribution + n,nrow=4) + scale_x_reverse() + theme_classic() + labs(x='% studies selected',y='power')
ggsave('ms/210314a.png')

dd



## 51. convergence of estimated bias term
require(parallel)
require(EnvStats)
n <- 1e1
B <- 1e3
rSs <- c(uniform=function(n)runif(n,0,1), exponential=rexp, gamma=function(n)rgamma(n,shape=.54), beta=function(n)rbeta(n,.15,.39), pareto=function(n)rpareto(n,location=1,shape=2.52687))
biases <- c(uniform=.11, beta=.21, exponential=.16, gamma=.18, pareto=.04)
s.mds <- c(uniform=1/3,exponential=1)
s2s <- c(uniform=1/3,exponential=2)
distributions <- names(rSs)
distributions <- c('uniform','exponential')
names(distributions) <- distributions
ns <- round(seq(10,2e2,len=10))
## distr <- 'uniform'
by.distr <- lapply(distributions, function(distr) {
    rS <- rSs[[distr]]; s.md <- unname(s.mds[distr]); s2 <- unname(s2s[distr])
    by.n <- lapply(ns, function(n) {
        bias.hats <- mclapply(1:B,mc.cores=detectCores()-2,FUN=function(jj) {
            s <- rS(n)
            s.md.hat <- mean(abs(outer(s,s,`-`)))
            s2.hat <- mean(s^2)
            bias.hat <- s.md.hat^2/s2.hat/pi
            bias <- s.md^2/s2/pi
            c(md.bias=s.md.hat-s.md,s2.bias=s2.hat-s2,bias.bias=bias.hat-bias)
        })
        bias.hats <- simplify2array(bias.hats)
        ## hist(bias.hats);abline(v=bias,col='red')
        ## c(mean=mean(bias.hats-bias),sd=sd(bias.hats))
        list(means=rowMeans(bias.hats),sds=apply(bias.hats,1,sd))
    })
})
bias.means <- lapply(by.distr,function(lst)sapply(lst,function(lst1)lst1$means))
op <- par(mfrow=c(1,3)) 
for (bias.type in c('md.bias','s2.bias','bias.bias')) {
matplot(ns,sapply(bias.means,function(m)m[bias.type,]),col=1,type='l',lty=1:length(distributions))
legend('bottomright',lty=1:length(distributions),legend=distributions)
abline(h=0,col='red')
}
par(op)
