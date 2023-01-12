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

beta <- runif(1); gamma <- runif(1)

alpha <- runif(1);
integrate(function(u)pnorm(alpha*u)*pnorm(beta*u)*dnorm(gamma*u),-Inf,Inf)$val-1/gamma/2/pi*atan(alpha*beta/gamma/sqrt(alpha^2+beta^2+gamma^2))



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
supp.S <- c(0,5)
dS <- function(s)dunif(s,supp.S[1],supp.S[2])
pS <- function(q)punif(q,supp.S[1],supp.S[2])
rS <- function(n)runif(n,supp.S[1],supp.S[2])
E.S2 <- diff(supp.S)^2/12 + mean(supp.S)^2
mean.S.pair <- with(list(b=supp.S[2],a=supp.S[1]), 1/(b-a)*((a^2+b^2)/3-2/3*a*b))
n <- 3e2
ns <- round(seq(1e2,1e3,len=20))
## by.n <- mclapply(ns, mc.cores=detectCores()-2, FUN=function(n) {
by.n <- sapply(ns, FUN=function(n) {
    tau.stats <- mclapply(1:1e3, mc.cores=detectCores()-2, FUN=function(jj) {
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
plot(ns,ns*by.n-4/9,type='l',main='correction')
lines(ns,4*(ns-2)*(ns-3)/(ns-1)*cov.1234,lty=2)
legend('topright',lty=1:2,legend=c('observed','formula'))
## dev.off()

## level calc
require(parallel)
source('misc.R')
alpha <- .01
q <- qnorm(1-alpha/2)
n <- 5e1
a <- .0
b <- 1
as <- seq(0,1-.1,length.out=20)
bs <- seq(.1,5,length.out=20)
as <- seq(0,1-.1,length.out=20)
bs <- seq(.5,5,length.out=20)
by.endpoint <- sapply(bs, function(b) {
    attach(S.env(c(b-b[1],b)))
    rejected <- mclapply(1:1e4, mc.cores=detectCores()-2, FUN=function(jj) {
        s <- rS(n)
        z <- rnorm(n,1)
        y <- z/s
        tau.stat <- tau(z,s)
        sqrt(9/4*n)*abs(tau.stat) > q
    })
    rejected <- simplify2array(rejected)
    detach(S.env(c(b-b[1],b)))
    mean(rejected)    
})
plot(bs,by.endpoint)
