## strange R behavior: seems assigning null doesn't mask the
## function. inside f a has value null, but still evaluates an
## argument. as if the dispatch goes up to the global environment when
## it encounters a null function.
a <- exp
f <- function(){
    a <- NULL
    a(3)
    }
f()


## checking order of auc.hat.beta.hat - auc.hat.gamma.hat at beta==gamma. does appear to
## be O(1/n).
source('misc.R')
set.seed(1)
p <- 3
n <- 1e3
ns <- round(seq(1e2,1e4,len=2e1))
## by.n <- sapply(ns, function(n) {
## terms <- replicate(2e1, {
## tryCatch({
beta <- runif(p+1)/4
beta[p+1] <- 0 
## gamma <- beta; gamma[p+1] <- 0
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta.sampler <- lda.sampler.init(beta,Sigma)        
## mu.diff <- c(0,mu.0-mu.1)
Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
by.n <- sapply(ns, function(n) {
    pairs <- replicate(1e2, {
        list2env(beta.sampler$sample(n),globalenv())
        ## data <- beta.sampler$sample(n)
        ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
        beta.hat <- lda.coefs(x.0,x.1)
        gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
        gamma.hat <- c(gamma.hat,0)
        ## scores.0 <- cbind(x.0%*%beta.hat,x.0%*%gamma.hat)
        ## scores.1 <- cbind(x.1%*%beta.hat,x.1%*%gamma.hat)
        ## var.diff <- delong.var(x=scores.0,y=scores.1)
        auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat.continuous(x.0%*%gamma.hat,x.1%*%gamma.hat)
        auc.hat.beta.hat-auc.hat.gamma.hat
    })
})
matplot(ns,ns*t(by.n),pch=1,col=1);abline(h=0)
## matplot(ns,sqrt(ns)*t(by.n),pch=1,col=1);abline(h=0)



## order of auc.hajek.beta.hat - auc.hajek.gamma.hat at
## beta==gamma. didnt appear at first to be O(1/n), but with larger n
## (>=2000 for set.seed(1)) it does. however the difference of the
## difference of the hatted estimates and the difference of the hajek
## estimates also appears to be O(1/n) (saved image for 40min sim)
start <- Sys.time()
source('misc.R')
set.seed(1)
p <- 3
## n <- 1e3
ns <- round(seq(1e2,3e4,len=2e1))
## by.n <- sapply(ns, function(n) {
## terms <- replicate(2e1, {
## tryCatch({
beta <- runif(p+1)/4
beta[p+1] <- 0 
## gamma <- beta; gamma[p+1] <- 0
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta.sampler <- lda.sampler.init(beta,Sigma)        
## mu.diff <- c(0,mu.0-mu.1)
Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
by.n <- lapply(ns, function(n) {
    pairs <- replicate(1e2, {
        list2env(beta.sampler$sample(n),globalenv())
        ## data <- beta.sampler$sample(n)
        ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
        beta.hat <- lda.coefs(x.0,x.1)
        gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
        gamma.hat <- c(gamma.hat,0)
        ## scores.0 <- cbind(x.0%*%beta.hat,x.0%*%gamma.hat)
        ## scores.1 <- cbind(x.1%*%beta.hat,x.1%*%gamma.hat)
        ## var.diff <- delong.var(x=scores.0,y=scores.1)
        auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat.continuous(x.0%*%gamma.hat,x.1%*%gamma.hat)
        hat.diff <- auc.hat.beta.hat-auc.hat.gamma.hat
        auc.hajek.beta.hat <- with(beta.sampler,auc.scores.hajek(x.0,x.1,c(1,mu.0),c(1,mu.1),beta.hat,Sigma.diff/2))
        auc.hajek.gamma.hat <- with(beta.sampler,auc.scores.hajek(x.0,x.1,c(1,mu.0),c(1,mu.1),gamma.hat,Sigma.diff/2))
        hajek.diff <- auc.hajek.beta.hat - auc.hajek.gamma.hat
        ## hat.diff - hajek.diff
        c(hat.diff=hat.diff,hajek.diff=hajek.diff)
    })
})
Sys.time() - start
## save.image('sessions/220930.RData')
hajek.diffs <- sapply(by.n,function(df)df['hajek.diff',])
hat.diffs <- sapply(by.n,function(df)df['hat.diff',])
matplot(ns,ns^1*t(hajek.diffs),pch=1,col=1);abline(h=0)
matplot(ns,ns^1*t(hat.diffs),pch=1,col=1);abline(h=0)
matplot(ns,ns^2*t(hat.diffs-hajek.diffs),pch=1,col=1);abline(h=0)
## matplot(ns,sqrt(ns)*t(by.n),pch=1,col=1);abline(h=0)
plot(ns, colMeans(hat.diffs-hajek.diffs))





## as above, but trying noise as beta and gamma perturbations. still
## seems O(1/n). so issue is with using diff of hajek projections.
start <- Sys.time()
source('misc.R')
set.seed(1)
p <- 3
## n <- 1e3
ns <- round(seq(1e2,1e4,len=4e1))
## by.n <- sapply(ns, function(n) {
## terms <- replicate(2e1, {
## tryCatch({
beta <- runif(p+1)/4
## gamma <- beta; gamma[p+1] <- 0
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta.sampler <- lda.sampler.init(beta,Sigma)        
## mu.diff <- c(0,mu.0-mu.1)
Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
by.n <- lapply(ns, function(n) {
    pairs <- replicate(1e2, {
        list2env(beta.sampler$sample(n),globalenv())
        ## data <- beta.sampler$sample(n)
        ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
        beta.hat <- beta+rnorm(length(beta))/sqrt(n)#lda.coefs(x.0,x.1)
        gamma.hat <- beta+rnorm(length(beta))/sqrt(n)
        ## scores.0 <- cbind(x.0%*%beta.hat,x.0%*%gamma.hat)
        ## scores.1 <- cbind(x.1%*%beta.hat,x.1%*%gamma.hat)
        ## var.diff <- delong.var(x=scores.0,y=scores.1)
        auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat.continuous(x.0%*%gamma.hat,x.1%*%gamma.hat)
        hat.diff <- auc.hat.beta.hat-auc.hat.gamma.hat
        auc.hajek.beta.hat <- with(beta.sampler,auc.scores.hajek(x.0,x.1,c(1,mu.0),c(1,mu.1),beta.hat,Sigma.diff/2))
        auc.hajek.gamma.hat <- with(beta.sampler,auc.scores.hajek(x.0,x.1,c(1,mu.0),c(1,mu.1),gamma.hat,Sigma.diff/2))
        hajek.diff <- auc.hajek.beta.hat - auc.hajek.gamma.hat
        ## hat.diff - hajek.diff
        c(hat.diff=hat.diff,hajek.diff=hajek.diff)
    })
})
Sys.time() - start
## save.image('220930a.RData')
hajek.diffs <- sapply(by.n,function(df)df['hajek.diff',])
hat.diffs <- sapply(by.n,function(df)df['hat.diff',])
matplot(ns,ns^1*t(hajek.diffs),pch=1,col=1);abline(h=0)
matplot(ns,ns^1*t(hat.diffs),pch=1,col=1);abline(h=0)
matplot(ns,ns^2*t(hat.diffs-hajek.diffs),pch=1,col=1);abline(h=0)
## matplot(ns,sqrt(ns)*t(by.n),pch=1,col=1);abline(h=0)
plot(ns, ns^1*colMeans(hat.diffs-hajek.diffs))



source('misc.R')
set.seed(1)
p <- 3
n <- 5e2
beta <- runif(p+1)
beta[p+1] <- 0
gamma <- beta
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta.sampler <- lda.sampler.init(beta,Sigma)
ns <- round(seq(1e2,1e4,len=20))
by.n <- sapply(ns, function(n) {
    pairs <- replicate(1e2, {
        ##     tryCatch({
        ##     beta <- runif(p+1)
        beta.hat <- beta+rnorm(p+1)/sqrt(n)
        ##     Sigma <- matrix(rnorm(p^2),nrow=p)
        ##     Sigma <- Sigma%*%t(Sigma)
        ## list2env(lda.sampler.init(beta,Sigma),globalenv())
        mu.diff <- with(beta.sampler,c(0,mu.1-mu.0))
        Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
        list2env(beta.sampler$sample(n),globalenv())
        beta.hat <- lda.coefs(x.0,x.1)
        ## gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
        ## gamma.hat <- c(gamma.hat,0)
        auc.beta <- auc.scores(beta,mu.diff,Sigma.diff)
        auc.beta.hat <- auc.scores(beta.hat,mu.diff,Sigma.diff)
        auc.hat.beta <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
        auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
        ## c(auc.hat.beta,auc.hat.beta.hat)
        ## auc.scores(coefs=beta,mu.diff,Sigma.diff)
        deriv.2 <- auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)
        ## observed <- auc.hat.beta.hat
        ## taylor <- auc.hat.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
        observed <- auc.hat.beta.hat - auc.hat.beta
        taylor <-  -t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
        ## observed <- auc.beta.hat - auc.beta
        ## taylor <-  t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
        ## observed <- auc.hat.beta.hat - auc.scores(beta,mu.diff,Sigma.diff)#auc.hat.beta
        ## observed <- auc.hat.beta.hat - auc.beta
        ## taylor <- auc.hat.beta-auc.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
        ##----
        ## observed <- auc.hat.beta.hat - auc.hat.beta
        ## observed <- auc.beta.hat - auc.beta
        ## taylor <- t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
        ## observed <- auc.hat.beta.hat
        ## taylor <- auc.beta.hat
        ## observed <- auc.hat.beta
        ## taylor <- auc.beta
        ## observed <-  auc.hat.beta.hat - auc.hat.beta 
        ## taylor <- auc.beta.hat - auc.beta
        ## observed <-  auc.hat.beta.hat - auc.beta.hat
        ## taylor <-  auc.hat.beta  - auc.beta
        c(observed=observed,taylor=taylor)
        ## },error=function(e)c(NA,NA))
    })
    pairs['observed',]-pairs['taylor',]
})
matplot(ns,ns^1*t(by.n),pch=1,col=1)
plot(ns,ns^1*apply(by.n,2,sd))




source('misc.R')
set.seed(1)
p <- 3
n <- 5e2
beta <- runif(p+1)
beta[p+1] <- 0
gamma <- beta
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta.sampler <- lda.sampler.init(beta,Sigma)
ns <- round(seq(1e2,1e4,len=20))
by.n <- sapply(ns, function(n) {
    pairs <- replicate(1e2, {
        ##     tryCatch({
        ##     beta <- runif(p+1)
        beta.hat <- beta+rnorm(p+1)/sqrt(n)
        gamma.hat <- beta+rnorm(p+1)/sqrt(n)
        ##     Sigma <- matrix(rnorm(p^2),nrow=p)
        ##     Sigma <- Sigma%*%t(Sigma)
        ## list2env(lda.sampler.init(beta,Sigma),globalenv())
        ## mu.diff <- with(beta.sampler,c(0,mu.1-mu.0))
        ## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
        list2env(beta.sampler$sample(n),globalenv())
        ## beta.hat <- lda.coefs(x.0,x.1)
        ## gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
        ## gamma.hat <- c(gamma.hat,0)
        ## auc.beta <- auc.scores(beta,mu.diff,Sigma.diff)
        ## auc.beta.hat <- auc.scores(beta.hat,mu.diff,Sigma.diff)
        auc.hat.beta <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
        auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat.continuous(x.0%*%gamma.hat,x.1%*%gamma.hat)
        ## c(auc.hat.beta,auc.hat.beta.hat)
        ## auc.scores(coefs=beta,mu.diff,Sigma.diff)
        deriv.2 <- auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)
        ## observed <- auc.hat.beta.hat
        ## taylor <- auc.hat.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
        observed <- auc.hat.beta.hat - auc.hat.beta - (auc.hat.gamma.hat - auc.hat.gama)
        taylor <-  -t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2 - (-t(gamma.hat-beta)%*%deriv.2%*%(gamma.hat-beta)/2)
        ## observed <- auc.beta.hat - auc.beta
        ## taylor <-  t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
        ## observed <- auc.hat.beta.hat - auc.scores(beta,mu.diff,Sigma.diff)#auc.hat.beta
        ## observed <- auc.hat.beta.hat - auc.beta
        ## taylor <- auc.hat.beta-auc.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
        ##----
        ## observed <- auc.hat.beta.hat - auc.hat.beta
        ## observed <- auc.beta.hat - auc.beta
        ## taylor <- t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
        ## observed <- auc.hat.beta.hat
        ## taylor <- auc.beta.hat
        ## observed <- auc.hat.beta
        ## taylor <- auc.beta
        ## observed <-  auc.hat.beta.hat - auc.hat.beta 
        ## taylor <- auc.beta.hat - auc.beta
        ## observed <-  auc.hat.beta.hat - auc.beta.hat
        ## taylor <-  auc.hat.beta  - auc.beta
        c(observed=observed,taylor=taylor)
        ## },error=function(e)c(NA,NA))
    })
    pairs['observed',]-pairs['taylor',]
})
matplot(ns,ns^1*t(by.n),pch=1,col=1)
plot(ns,ns^1*apply(by.n,2,sd))


dd


##  consistency of lda and logistic coef estimates under normal model from efron '75
source('misc.R')
require(mvtnorm)
set.seed(1)
auc.continuous <- function(x,y)mean(outer(x,y,'<'))
lda <- function(x.0,x.1) {
    n.0 <- nrow(x.0); n.1 <- nrow(x.1); n <- n.0+n.1
    pi.0.hat <- n.0/n; pi.1.hat <- n.1/n
    mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
    Sigma.hat <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n
    Sigma.hat.inv <- solve(Sigma.hat)
    beta.hat <- c(log(pi.1.hat/pi.0.hat) - 1/2*(t(mu.1.hat)%*%Sigma.hat.inv%*%mu.1.hat - t(mu.0.hat)%*%Sigma.hat.inv%*%mu.0.hat), t(mu.1.hat-mu.0.hat)%*%Sigma.hat.inv)
    return(beta.hat)
}
p <- 3
n <-5e2
beta <- runif(p+1)
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
Sigma.inv <- solve(Sigma)
mu.0 <- rep(0,p)
mu.1 <- mu.0+Sigma%*%beta[-1]
pi.1 <- plogis(beta[1]+1/2*(t(mu.1)%*%Sigma.inv%*%mu.1 - t(mu.0)%*%Sigma.inv%*%mu.0))
pi.0 <- 1-pi.1
ns <- round(seq(1e2,3e2,len=10))
by.n <- sapply(ns, function(n) {
    beta.hats <- replicate(1e2, {
        n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
        x.0 <- rmvnorm(n.0,mean=mu.0,sigma=Sigma)
        x.1 <- rmvnorm(n-n.0,mean=mu.1,sigma=Sigma)
        beta.hat.lda <- lda(x.0,x.1)
        g <- rep(0:1,c(n.0,n.1))
        beta.hat.logistic <- unname(coef(glm(g ~ rbind(x.0,x.1), family=binomial)))
        rbind(lda=beta.hat.lda,logistic=beta.hat.logistic)
    })
    c(lda=mean((beta.hats['lda',,]-beta)^2),logistic=mean((beta.hats['logistic',,]-beta)^2))
})
plot(ns,by.n['lda',],ylim=range(by.n),type='l')
lines(ns,by.n['logistic',],col=2)



##  consistency of lda and logistic coef estimates under normal model from efron '75
##  refactored
source('misc.R')
require(mvtnorm)
set.seed(1)
p <- 3
n <-1e3
beta <- runif(p+1)
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
## Sigma.inv <- solve(Sigma)
## mu.0 <- rep(0,p)
## mu.1 <- mu.0+Sigma%*%beta[-1]
## pi.1 <- plogis(beta[1]+1/2*(t(mu.1)%*%Sigma.inv%*%mu.1 - t(mu.0)%*%Sigma.inv%*%mu.0))
## pi.0 <- 1-pi.1
lda.sampler <- lda.sampler.init(beta,Sigma)
ns <- round(seq(1e2,3e3,len=10))
by.n <- sapply(ns, function(n) {
    beta.hats <- replicate(1e2, {
        ## n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
        ## x.0 <- rmvnorm(n.0,mean=mu.0,sigma=Sigma)
        ## x.1 <- rmvnorm(n-n.0,mean=mu.1,sigma=Sigma)
        list2env(lda.sampler$sample(n),globalenv())
        beta.hat.lda <- lda.coefs(x.0,x.1)
        g <- rep(0:1,c(n.0,n.1))
        beta.hat.logistic <- unname(coef(glm(g ~ rbind(x.0[,-1],x.1[,-1]), family=binomial)))
        rbind(lda=beta.hat.lda,logistic=beta.hat.logistic)
    })
    c(lda=mean((beta.hats['lda',,]-beta)^2),logistic=mean((beta.hats['logistic',,]-beta)^2))
})
plot(ns,by.n['lda',],ylim=range(by.n),type='l')
lines(ns,by.n['logistic',],col=2)




##  consistency of lda and logistic coef estimates when estimation
##  uses the reduced model and last coef is 0.
source('misc.R')
## set.seed(6)
require(mvtnorm)
## auc.continuous <- function(x,y)mean(outer(x,y,'<'))
## lda <- function(x.0,x.1) {
##     n.0 <- nrow(x.0); n.1 <- nrow(x.1); n <- n.0+n.1
##     pi.0.hat <- n.0/n; pi.1.hat <- n.1/n
##     mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
##     Sigma.hat <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n
##     Sigma.hat.inv <- solve(Sigma.hat)
##     beta.hat <- c(log(pi.1.hat/pi.0.hat) - 1/2*(t(mu.1.hat)%*%Sigma.hat.inv%*%mu.1.hat - t(mu.0.hat)%*%Sigma.hat.inv%*%mu.0.hat), t(mu.1.hat-mu.0.hat)%*%Sigma.hat.inv)
##     return(beta.hat)
## }
p <- 3
## n <-5e2
beta <- c(runif(p),0)
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
## Sigma <- diag(p)
Sigma.inv <- solve(Sigma)
mu.0 <- rep(0,p)
mu.1 <- mu.0+Sigma%*%beta[-1]
pi.1 <- plogis(beta[1]+1/2*(t(mu.1)%*%Sigma.inv%*%mu.1 - t(mu.0)%*%Sigma.inv%*%mu.0))
pi.0 <- 1-pi.1
ns <- round(seq(4e2,3e3,len=10))
by.n <- sapply(ns, function(n) {
    beta.hats <- replicate(1e2, {
        n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
        x.0 <- rmvnorm(n.0,mean=mu.0,sigma=Sigma)
        x.1 <- rmvnorm(n-n.0,mean=mu.1,sigma=Sigma)
        beta.hat.lda <- lda.coefs(x.0[,-p],x.1[,-p])
        g <- rep(0:1,c(n.0,n.1))
        beta.hat.logistic <- unname(coef(glm(g ~ rbind(x.0[,-p],x.1[,-p]), family=binomial)))
        rbind(lda=beta.hat.lda,logistic=beta.hat.logistic)
    })
    c(lda=mean((beta.hats['lda',,]-beta[-(p+1)])^2),logistic=mean((beta.hats['logistic',,]-beta[-(p+1)])^2))
})
plot(ns,by.n['lda',],ylim=range(by.n),type='l')
lines(ns,by.n['logistic',],col=2)

## collapsibility of lda
t(mu.1-mu.0)%*%solve(Sigma)
t(mu.1[-p]-mu.0[-p])%*%solve(Sigma[-p,-p])

beta <- c(runif(p-1),0)
mu <- Sigma%*%beta
solve(Sigma)%*%mu
solve(Sigma[-p,-p])%*%mu[-p]

require(expm)
p <- 5
v <- runif(p)
v[p] <- 0
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
(solve(Sigma))%*%v
(solve(Sigma[-p,-p]))%*%v[-p]


source('misc.R')
## set.seed(6)
require(mvtnorm)
p <- 3
n <-5e2
beta <- c(runif(p),0)
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
## Sigma <- diag(p)
Sigma.inv <- solve(Sigma)
mu.0 <- rep(0,p)
mu.1 <- mu.0+Sigma%*%beta[-1]
pi.1 <- plogis(beta[1]+1/2*(t(mu.1)%*%Sigma.inv%*%mu.1 - t(mu.0)%*%Sigma.inv%*%mu.0))
pi.0 <- 1-pi.1
ns <- round(seq(4e2,6e3,len=10))
by.n <- sapply(ns, function(n) {
    ## z.stats <- replicate(1e2, {
    ##     n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
    ##     x.0 <- rmvnorm(n.0,mean=mu.0,sigma=Sigma)
    ##     x.1 <- rmvnorm(n-n.0,mean=mu.1,sigma=Sigma)
    ##     x <- cbind(1,rbind(x.0,x.1))
    ##     beta.hat.full <- lda(x.0,x.1)
    ##     beta.hat.reduced <- lda(x.0[,-p],x.1[,-p])
    ##     g <- rep(0:1,c(n.0,n.1))
    ##     ## beta.hat.logistic <- unname(coef(glm(g ~ rbind(x.0[,-p],x.1[,-p]), family=binomial)))
    ##     ## rbind(lda=beta.hat.lda,logistic=beta.hat.logistic)
    ##     marker.full <- x%*%beta.hat.full
    ##     marker.reduced <- x[,-p]%*%beta.hat.reduced
    ##     delong.test(g=g,xy=cbind(marker.full,marker.reduced))
    ## })
    ## p.vals <- 1-pnorm(z.stats)
    ## plot(ecdf(p.vals))
    ## abline(0,1)
    beta.hats <- replicate(1e2, {
        n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
        x.0 <- rmvnorm(n.0,mean=mu.0,sigma=Sigma)
        x.1 <- rmvnorm(n-n.0,mean=mu.1,sigma=Sigma)
        x <- cbind(1,rbind(x.0,x.1))
        beta.hat.full <- lda.coefs(x.0,x.1)
        beta.hat.reduced <- lda.coefs(x.0[,-p],x.1[,-p])
        ## g <- rep(0:1,c(n.0,n.1))
        ## ## beta.hat.logistic <- unname(coef(glm(g ~ rbind(x.0[,-p],x.1[,-p]), family=binomial)))
        ## ## rbind(lda=beta.hat.lda,logistic=beta.hat.logistic)
        ## marker.full <- x%*%beta.hat.full
        ## marker.reduced <- x[,-(p+1)]%*%beta.hat.reduced
        ## delong.test(g=g,xy=cbind(marker.full,marker.reduced))
        rbind(full=beta.hat.full,reduced=c(beta.hat.reduced,0))
    })
    c(full=mean((beta.hats['full',,]-beta)^2),reduced=mean((beta.hats['reduced',,]-beta)^2))
})
plot(ns,by.n['full',],ylim=range(by.n),type='l')
lines(ns,by.n['reduced',],col=2)

## delong test using lda scorees as markers. outcome: markers for full
## and reduced models look the same by qqplot or ks test, but delong
## test pvalues are not uniform
source('misc.R')
set.seed(6)
require(mvtnorm)
p <- 3
n <-1e3
beta <- c(runif(p),0)
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
lda.sampler <- lda.sampler.init(beta,Sigma)
ns <- round(seq(3e2,6e3,len=10))
z.stats <- replicate(3e2, {
    ## n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
    ## x.0 <- rmvnorm(n.0,mean=mu.0,sigma=Sigma)
    ## x.1 <- rmvnorm(n-n.0,mean=mu.1,sigma=Sigma)
    list2env(lda.sampler(n),globalenv())
    x <- cbind(1,rbind(x.0,x.1))
    beta.hat.full <- lda.coefs(x.0,x.1)
    beta.hat.reduced <- lda.coefs(x.0[,-p],x.1[,-p])
    beta.hat.reduced <- c(beta.hat.reduced,0)
    ## beta.hat.full <- beta+rnorm(p+1)/n^(1.3)
    ## beta.hat.reduced <- beta+rnorm(p+1)/n^(1.3)
    g <- rep(0:1,c(n.0,n.1))
    ## ## beta.hat.logistic <- unname(coef(glm(g ~ rbind(x.0[,-p],x.1[,-p]), family=binomial)))
    ## ## rbind(lda=beta.hat.lda,logistic=beta.hat.logistic)
    marker.full <- x%*%beta.hat.full
    marker.reduced <- x%*%beta.hat.reduced
    delong.test(g=g,xy=cbind(marker.full,marker.reduced))
})
p.vals <- 1-pnorm(z.stats)
plot(ecdf(p.vals))
abline(0,1)


dd




require(MASS)
p <- 1
n <- 1e1
x <- matrix(rnorm(n*p),nrow=n)
g <- sample(0:1,n,replace=TRUE)
g <- rep(0:1,each=n/2)
## compute LDA
x.demeaned <- scale(x,scale=FALSE)#t(t(x) - colMeans(x))
x.grouped <- split.data.frame(x,g)
W.i <- lapply(x.grouped, function(x.i) {
    x.i.demeaned <- scale(x.i,scale=FALSE)#t(t(x.i)-colMeans(x.i))
    t(x.i.demeaned)%*%x.i.demeaned
})
group.means <- do.call(cbind,lapply(split.data.frame(x,g), colMeans))
d <- apply(group.means,1,diff)
W <- Reduce(`+`,W.i)
ns <- sapply(x.grouped,nrow)
## (group.means - colMeans(x))
a <- solve(W)%*%d
proj <- t(a) %*% (t(x) - rowMeans(group.means))
g.hat <- as.numeric(proj < 0)
lda(x,g,method='mle')
S.i <- lapply(x.grouped,cov)
S <- (S.i[['0']]*ns['0'] + S.i[['1']]*ns['1']) / (n-2)
## Sigma <- cov(x)
w <- solve(S)%*%d
c <- w%*%rowMeans(group.means)
w












## lda.dbg <-
##     function(x, grouping, prior = proportions, tol = 1.0e-4,
##              method = c("moment", "mle", "mve", "t"),
##              CV = FALSE, nu = 5, ...)
## {
##     browser()
##     if(is.null(dim(x))) stop("'x' is not a matrix")
##     x <- as.matrix(x)
##     if(any(!is.finite(x)))
##         stop("infinite, NA or NaN values in 'x'")
##     n <- nrow(x)
##     p <- ncol(x)
##     if(n != length(grouping))
##         stop("nrow(x) and length(grouping) are different")
##     g <- as.factor(grouping)
##     lev <- lev1 <- levels(g)
##     counts <- as.vector(table(g))
##     if(!missing(prior)) {
##         if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
##         if(length(prior) != nlevels(g)) stop("'prior' is of incorrect length")
##         prior <- prior[counts > 0L]
##     }
##     if(any(counts == 0L)) {
##         empty <- lev[counts == 0L]
##         warning(sprintf(ngettext(length(empty),
##                                  "group %s is empty",
##                                  "groups %s are empty"),
##                         paste(empty, collapse = " ")), domain = NA)
##         lev1 <- lev[counts > 0L]
##         g <- factor(g, levels = lev1)
##         counts <- as.vector(table(g))
##     }
##     proportions <- counts/n
##     ng <- length(proportions)
##     names(prior) <- names(counts) <- lev1
##     method <- match.arg(method)
##     if(CV && !(method == "moment" || method == "mle"))
##         stop(gettext("cannot use leave-one-out CV with method %s",
##                      sQuote(method)), domain = NA)
##     ## drop attributes to avoid e.g. matrix() methods
##     group.means <- tapply(c(x), list(rep(g, p), col(x)), mean)
##     f1 <- sqrt(diag(var(x - group.means[g,  ])))
##     if(any(f1 < tol)) {
##         const <- format((1L:p)[f1 < tol])
##         stop(sprintf(ngettext(length(const),
##                               "variable %s appears to be constant within groups",
##                               "variables %s appear to be constant within groups"),
##                      paste(const, collapse = " ")),
##              domain = NA)
##     }
##                                         # scale columns to unit variance before checking for collinearity
##     scaling <- diag(1/f1, , p)
##     if(method == "mve") {
##                                         # adjust to "unbiased" scaling of covariance matrix
##         cov <- n/(n - ng) * cov.rob((x - group.means[g,  ]) %*% scaling)$cov
##         sX <- svd(cov, nu = 0L)
##         rank <- sum(sX$d > tol^2)
##         if(rank == 0L) stop("rank = 0: variables are numerically constant")
##         if(rank < p) warning("variables are collinear")
##         scaling <- scaling %*% sX$v[, 1L:rank] %*%
##             diag(sqrt(1/sX$d[1L:rank]),,rank)
##     } else if(method == "t") {
##         if(nu <= 2) stop("'nu' must exceed 2")
##         w <- rep(1, n)
##         repeat {
##             w0 <- w
##             X <- x - group.means[g, ]
##             sX <- svd(sqrt((1 + p/nu)*w/n) * X, nu = 0L)
##             X <- X %*% sX$v %*% diag(1/sX$d,, p)
##             w <- 1/(1 + drop(X^2 %*% rep(1, p))/nu)
##             print(summary(w))
##             group.means <- tapply(w*x, list(rep(g, p), col(x)), sum)/
##                 rep.int(tapply(w, g, sum), p)
##             if(all(abs(w - w0) < 1e-2)) break
##         }
##         X <-  sqrt(nu/(nu-2)*(1 + p/nu)/n * w) * (x - group.means[g,  ]) %*% scaling
##         X.s <- svd(X, nu = 0L)
##         rank <- sum(X.s$d > tol)
##         if(rank == 0L) stop("rank = 0: variables are numerically constant")
##         if(rank < p) warning("variables are collinear")
##         scaling <- scaling %*% X.s$v[, 1L:rank] %*% diag(1/X.s$d[1L:rank],,rank)
##     } else {
##         fac <- if(method == "moment") 1/(n-ng) else 1/n
##         X <- sqrt(fac) * (x - group.means[g,  ]) %*% scaling
##         X.s <- svd(X, nu = 0L)
##         rank <- sum(X.s$d > tol)
##         if(rank == 0L) stop("rank = 0: variables are numerically constant")
##         if(rank < p) warning("variables are collinear")
##         scaling <- scaling %*% X.s$v[, 1L:rank] %*% diag(1/X.s$d[1L:rank],,rank)
##     }
##                                         # now have variables scaled so that W is the identity
##     if(CV) {
##         x <- x %*% scaling
##         dm <- group.means %*% scaling
##         K <- if(method == "moment") ng else 0L
##         dist <- matrix(0, n, ng)
##         for(i in 1L:ng) {
##             dev <- x - matrix(dm[i,  ], n, rank, byrow = TRUE)
##             dist[, i] <- rowSums(dev^2)
##         }
##         ind <- cbind(1L:n, g)
##         nc <- counts[g]
##         cc <- nc/((nc-1)*(n-K))
##         dist2 <- dist
##         for(i in 1L:ng) {
##             dev <- x - matrix(dm[i,  ], n, rank, byrow = TRUE)
##             dev2 <- x - dm[g, ]
##             tmp <- rowSums(dev*dev2)
##             dist[, i] <- (n-1L-K)/(n-K) * (dist2[, i] +  cc*tmp^2/(1 - cc*dist2[ind]))
##         }
##         dist[ind] <- dist2[ind] * (n-1L-K)/(n-K) * (nc/(nc-1))^2 /
##             (1 - cc*dist2[ind])
##         dist <- 0.5 * dist - matrix(log(prior), n, ng, byrow = TRUE)
##         dist <- exp(-(dist - min(dist, na.rm = TRUE)))
##         cl <- factor(lev1[max.col(dist)], levels = lev)
##         ##  convert to posterior probabilities
##         posterior <- dist/drop(dist %*% rep(1, length(prior)))
##         dimnames(posterior) <- list(rownames(x), lev1)
##         return(list(class = cl, posterior = posterior))
##     }
##     xbar <- colSums(prior %*% group.means)
##     fac <- if(method == "mle") 1/ng else 1/(ng - 1)
##     X <- sqrt((n * prior)*fac) * scale(group.means, center = xbar, scale = FALSE) %*% scaling
##     X.s <- svd(X, nu = 0L)
##     rank <- sum(X.s$d > tol * X.s$d[1L])
##     if(rank == 0L) stop("group means are numerically identical")
##     scaling <- scaling %*% X.s$v[, 1L:rank]
##     if(is.null(dimnames(x)))
##         dimnames(scaling) <- list(NULL, paste("LD", 1L:rank, sep = ""))
##     else {
##         dimnames(scaling) <- list(colnames(x), paste("LD", 1L:rank, sep = ""))
##         dimnames(group.means)[[2L]] <- colnames(x)
##     }
##     cl <- match.call()
##     cl[[1L]] <- as.name("lda")
##     structure(list(prior = prior, counts = counts, means = group.means,
##                    scaling = scaling, lev = lev, svd = X.s$d[1L:rank],
##                    N = n, call = cl),
##               class = "lda")
## }
## lda.dbg(x,g)






## replicate pROC's delong test
require(pROC)
set.seed(123)
df <- data.frame(disease.status = rbinom(n=100, size=1, prob=0.20),
                 test1 = rnorm(100, mean=15, sd=4),
                 test2 = rnorm(100, mean=30, sd=2),
                 test3 = rnorm(100, mean=50, sd=3))
##create roc object for test1, test2, test3roc.out <- test1.roc(df$disease.status, df$test1, plot=TRUE, smooth = FALSE)
roc.out.test1<-roc(df$disease.status, df$test1, plot=TRUE, smooth = FALSE)
roc.out.test2 <- roc(df$disease.status, df$test2, plot=TRUE, smooth = FALSE)
roc.out.test3 <- roc(df$disease.status, df$test3, plot=TRUE, smooth = FALSE)
                                        #compare the AUC of test1 and test 2
roc.test(roc.out.test1, roc.out.test2, reuse.auc=TRUE, method="delong", na.rm=TRUE)

x <- df[df$disease.status==0,c('test1','test2')]
y <- df[df$disease.status==1,c('test1','test2')]
m <- nrow(x); n <- nrow(y)
theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
V.10.1 <- sapply(x[,1], function(x.i)mean(x.i<y[,1]))
V.10.2 <- sapply(x[,2], function(x.i)mean(x.i<y[,2]))
V.01.1 <- sapply(y[,1], function(y.i)mean(x[,1]<y.i))
V.01.2 <- sapply(y[,2], function(y.i)mean(x[,2]<y.i))
S.10 <- cov(cbind(V.10.1,V.10.2))
S.01 <- cov(cbind(V.01.1,V.01.2))
S <- S.10/m + S.01/n
contrast <- matrix(c(1,-1),ncol=1)
z.stat <- diff(rev(theta.hats)) / sqrt(t(contrast)%*%S%*%contrast)

## check performance on null data--looks good
delong.test <- function(x,y) {
    m <- nrow(x); n <- nrow(y)
    theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
    V.10.1 <- sapply(x[,1], function(x.i)mean(x.i<y[,1]))
    V.10.2 <- sapply(x[,2], function(x.i)mean(x.i<y[,2]))
    V.01.1 <- sapply(y[,1], function(y.i)mean(x[,1]<y.i))
    V.01.2 <- sapply(y[,2], function(y.i)mean(x[,2]<y.i))
    S.10 <- cov(cbind(V.10.1,V.10.2))
    S.01 <- cov(cbind(V.01.1,V.01.2))
    S <- S.10/m + S.01/n
    contrast <- matrix(c(1,-1),ncol=1)
    z.stat <- diff(rev(theta.hats)) / sqrt(t(contrast)%*%S%*%contrast)
    return(z.stat)
}

m <- 15
n <- 20
z.stats <- replicate(1e2, {
    x <- matrix(rnorm(2*m),ncol=2)
    y <- matrix(rnorm(2*n),ncol=2)
    delong.test(x,y)
})
alpha <- .05
mean(abs(z.stats)>qnorm(1-alpha/2))

p.vals <- replicate(1e3, {
    x <- matrix(rnorm(2*m),ncol=2)
    y <- matrix(runif(2*n),ncol=2)
    z.stat <- delong.test(x,y)
    1-pnorm(z.stat)
})
plot(ecdf(p.vals))
abline(0,1)



## version accepting dataframe input
delong.test.old <- function(x,y) {
    m <- nrow(x); n <- nrow(y)
    theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
    V.10.1 <- sapply(x[,1], function(x.i)mean(x.i<y[,1]))
    V.10.2 <- sapply(x[,2], function(x.i)mean(x.i<y[,2]))
    V.01.1 <- sapply(y[,1], function(y.i)mean(x[,1]<y.i))
    V.01.2 <- sapply(y[,2], function(y.i)mean(x[,2]<y.i))
    S.10 <- cov(cbind(V.10.1,V.10.2))
    S.01 <- cov(cbind(V.01.1,V.01.2))
    S <- S.10/m + S.01/n
    contrast <- matrix(c(1,-1),ncol=1)
    z.stat <- diff(rev(theta.hats)) / sqrt(t(contrast)%*%S%*%contrast)
    return(z.stat)
}
delong.test <- function(xy,g) {
    g <- factor(g,labels=0:1)
    x <- xy[g==0,]; y <- xy[g==1,]
    m <- nrow(x); n <- nrow(y)
    theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
    V.10.1 <- sapply(x[,1], function(x.i)mean(x.i<y[,1]))
    V.10.2 <- sapply(x[,2], function(x.i)mean(x.i<y[,2]))
    V.01.1 <- sapply(y[,1], function(y.i)mean(x[,1]<y.i))
    V.01.2 <- sapply(y[,2], function(y.i)mean(x[,2]<y.i))
    S.10 <- cov(cbind(V.10.1,V.10.2))
    S.01 <- cov(cbind(V.01.1,V.01.2))
    S <- S.10/m + S.01/n
    contrast <- matrix(c(1,-1),ncol=1)
    z.stat <- diff(rev(theta.hats)) / sqrt(t(contrast)%*%S%*%contrast)
    return(z.stat)
}
p.vals <- replicate(1e3, {
    x <- matrix(rnorm(2*m),ncol=2)
    y <- matrix(runif(2*n),ncol=2)
    xy <- rbind(x,y)
    g <- rep(0:1,c(nrow(x),nrow(y)))
    z.stat <- delong.test(xy,g)
    stopifnot(z.stat==delong.test.old(x,y))
    1-pnorm(z.stat)
})
plot(ecdf(p.vals))
abline(0,1)


## now with x's and y's from scores
## require(MASS)
delong.test <- function(xy,g) {
    g <- factor(g,labels=0:1)
    x <- xy[g==0,]; y <- xy[g==1,]
    m <- nrow(x); n <- nrow(y)
    theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
    V.10.1 <- sapply(x[,1], function(x.i)mean(x.i<y[,1]))
    V.10.2 <- sapply(x[,2], function(x.i)mean(x.i<y[,2]))
    V.01.1 <- sapply(y[,1], function(y.i)mean(x[,1]<y.i))
    V.01.2 <- sapply(y[,2], function(y.i)mean(x[,2]<y.i))
    S.10 <- cov(cbind(V.10.1,V.10.2))
    S.01 <- cov(cbind(V.01.1,V.01.2))
    S <- S.10/m + S.01/n
    contrast <- matrix(c(1,-1),ncol=1)
    z.stat <- diff(rev(theta.hats)) / sqrt(t(contrast)%*%S%*%contrast)
    return(z.stat)
}
p <- 4
n <- 1e2
beta <- rep(1,p)
p.vals <- replicate(1e3, {
    c <- matrix(rnorm(n*p),ncol=p)
    true.probs <- plogis(c%*%beta)
    g <- rbinom(n,1,prob=true.probs)
    x.1 <- predict(MASS::lda(c[,1:2],g))$x
    x.2 <- predict(MASS::lda(c[,3:4],g))$x
    ## c.0 <- c[g==0,]; c.1 <- c[g==1,]
    ## beta.hat <- lda(c.0,c.1)
    ## gamma.hat <- lda(c.0[,1:3],c.1[,1:3])
    ## x.1 <- c%*%beta.hat[-1; x.2 <- c%*%gamma.hat
    xy <- cbind(x.1,x.2)
    z.stat <- delong.test(xy=xy,g=g)
    1-pnorm(z.stat)
})
plot(ecdf(p.vals))
abline(0,1)

require(mvtnorm)
n <- 1e3
p <- 2
beta <- c(rep(.1,p),0)
p.vals <- replicate(1e2, {
    Sigma <- matrix(rnorm((p+1)^2),nrow=p+1)
    Sigma <- Sigma%*%t(Sigma)
    ## Sigma.0 <- Sigma.1[1:p,1:p]
    ## mu.0 <- rep(0,p)
    ## mu.1 <- c(mu.0+1,0)
    ## x <- rmvnorm(n,mean=mu.0,sigma=Sigma.0)
    ## y <- rmvnorm(n,mean=mu.1,sigma=Sigma.1)
    x <- rmvnorm(n,sigma=Sigma)
    true.probs <- plogis(x%*%beta)
    g <- rbinom(n,1,prob=true.probs)
    Sigma.full <- Sigma
    x.bar.full <- sapply(split.data.frame(x,g),colMeans)
    d.full <- apply(x.bar.full,1,diff)
    a.full <- solve(Sigma.full)%*%d.full
    ## g.hat.full <- as.numeric(t(t(x) - rowMeans(x.bar.full))%*%a.full > 0)
    Sigma.small <- Sigma[1:p,1:p]
    x.bar.small <- x.bar.full[1:p,]
    d.small <- apply(x.bar.small,1,diff)
    a.small <- solve(Sigma.small)%*%d.small
    ## g.hat.small <- as.numeric(t(t(x[,1:p]) - rowMeans(x.bar.small))%*%a.small > 0)
    ## mean(g.hat.full==g.hat.small)
    x.hat.small <- x[,1:p]%*%a.small
    x.hat.full <- x%*%a.full
    ## plot(x.hat.small,x.hat.full,xlim=c(-1,1),ylim=c(-1,1)); abline(0,1)
    z.stat <- delong.test(xy=cbind(x.hat.small,x.hat.full),g=g)
    1-pnorm(z.stat)
})
hist(p.vals)


## separately estimated betahat
source('misc.R')
require(mvtnorm)
n <- 1e3
p <- 3
beta <- c(rep(.1,p),0)
p.vals <- replicate(1e2, {
    Sigma <- matrix(rnorm((p+1)^2),nrow=p+1)
    Sigma <- Sigma%*%t(Sigma)
    ## Sigma.0 <- Sigma.1[1:p,1:p]
    ## mu.0 <- rep(0,p)
    ## mu.1 <- c(mu.0+1,0)
    ## x <- rmvnorm(n,mean=mu.0,sigma=Sigma.0)
    ## y <- rmvnorm(n,mean=mu.1,sigma=Sigma.1)
    x <- rmvnorm(n,sigma=Sigma)
    true.probs <- plogis(x%*%beta)
    g <- rbinom(n,1,prob=true.probs)
    x.holdout <- rmvnorm(n,sigma=Sigma)
    true.probs <- plogis(x.holdout%*%beta)
    g.holdout <- rbinom(n,1,prob=true.probs)
    Sigma.full <- Sigma
    ## x.bar.full <- sapply(split.data.frame(x,g),colMeans)
    x.bar.full<- sapply(split.data.frame(x.holdout,g.holdout),colMeans)
    d.full <- apply(x.bar.full,1,diff)
    a.full <- solve(Sigma.full)%*%d.full
    ## g.hat.full <- as.numeric(t(t(x) - rowMeans(x.bar.full))%*%a.full > 0)
    Sigma.small <- Sigma[1:p,1:p]
    x.bar.small <- x.bar.full[1:p,]
    d.small <- apply(x.bar.small,1,diff)
    a.small <- solve(Sigma.small)%*%d.small
    ## g.hat.small <- as.numeric(t(t(x[,1:p]) - rowMeans(x.bar.small))%*%a.small > 0)
    ## mean(g.hat.full==g.hat.small)
    x.hat.small <- x[,1:p]%*%a.small
    x.hat.full <- x%*%a.full
    ## plot(x.hat.small,x.hat.full,xlim=c(-1,1),ylim=c(-1,1)); abline(0,1)
    z.stat <- delong.test(xy=cbind(x.hat.small,x.hat.full),g=g)
    1-pnorm(z.stat)
})
hist(p.vals)



## noise to generate betahat,gammahat, true beta and gamma are equal
p <- 3
n <-5e2
beta <- runif(p+1)
gamma <- beta
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
Sigma.inv <- solve(Sigma)
mu.0 <- rep(0,p)
mu.1 <- mu.0+Sigma%*%beta[-1]
pi.1 <- plogis(beta[1]+1/2*(t(mu.1)%*%Sigma.inv%*%mu.1 - t(mu.0)%*%Sigma.inv%*%mu.0))
pi.0 <- 1-pi.1
## ns <- round(seq(50,1e2,len=10))
## by.n <- sapply(ns, function(n) {
z.stats <- replicate(5e2, {
    n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
    x.0 <- rmvnorm(n.0,mean=mu.0,sigma=Sigma)
    x.1 <- rmvnorm(n-n.0,mean=mu.1,sigma=Sigma)
    x <- rbind(cbind(1,x.0),cbind(1,x.1))
    beta.hat <- beta+rnorm(p+1)/sqrt(n)
    gamma.hat <- gamma+rnorm(p+1)/sqrt(n)    
    ## beta.hat.lda <- lda(x.0,x.1)
    g <- rep(0:1,c(n.0,n.1))
    scores <- cbind(x%*%beta.hat,x%*%gamma.hat)
    delong.test(xy=scores,g=g)
    ## beta.hat.logistic <- unname(coef(glm(g ~ rbind(x.0,x.1), family=binomial)))
    ## rbind(lda=beta.hat.lda,logistic=beta.hat.logistic)
})
p.vals <- 1-pnorm(z.stats)
plot(ecdf(p.vals))
abline(0,1,col=2)



## source('misc.R')
## set.seed(6)
## require(mvtnorm)
## auc.continuous <- function(x,y)mean(outer(x,y,'<'))
## p <- 3
## n <-5e2
## beta <- runif(p)
## gamma <- runif(p)
## ## gamma <- beta;# gamma[p] <- 0
## beta <- gamma <- rep(1,p)#c(rep(1,p-1),0)#runif(p)#
## ## beta.gamma <- c(beta,gamma)
## Sigma <- matrix(rnorm(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
## mu.x <- rep(0,p)
## mu.y <- mu.x+Sigma%*%beta
## pi <- 1/2
## ## quad <- as.numeric(beta%*%Sigma%*%beta)
## ## D.prime.beta <- (Sigma%*%beta %*% (beta%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-beta%*%mu/sqrt(quad))
## ## D.prime.gamma <- (Sigma%*%beta %*% (beta%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-beta%*%mu/sqrt(quad))
## ## D.primes <- lapply(list(beta=beta,gamma=gamma), function(coefs) {
##     ## quad <- as.numeric(coefs%*%Sigma%*%coefs)
##     ## (Sigma%*%coefs %*% (coefs%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-coefs%*%mu/sqrt(quad))
## ## })
## ## D.prime <- c(D.primes$beta,D.primes$gamma)
## z.stats <- replicate(1e3, {
##     g <- rbinom(n,1,prob=pi)
    
##     beta.hat <- beta+rnorm(p)/sqrt(n)
##     gamma.hat <- gamma+rnorm(p)/sqrt(n)
##     ## beta.hat <- glm(g ~ covariates, family=binomial)
##     ## xy <- cbind(covariates%*%beta.hat, covariates[,]%*%gamma.hat[])
##     xy <- cbind(covariates%*%beta.hat, covariates%*%gamma.hat)
##     ## x <- rmvnorm(n,mu.x,Sigma.x)
##     ## y <- rmvnorm(n,mu.y,Sigma.y)
##     ## delong.test(x=cbind(x%*%beta.hat,x%*%gamma.hat),y=cbind(y%*%beta.hat,y%*%gamma.hat))
##     delong.test(g=g,xy=xy)
## })
## p.vals <- 1-pnorm(z.stats)
## plot(ecdf(p.vals))
## abline(0,1,col=2)


## ## glm to get beta.hat,gamma.hat
##   source('misc.R')
## ## set.seed(2)
## require(mvtnorm)
## auc.continuous <- function(x,y)mean(outer(x,y,'<'))
## p <- 3
## n <-1e2
## beta <- runif(p)
## gamma <- runif(p)
## ## gamma <- beta;# gamma[p] <- 0
## beta <- gamma <- rep(1,p)#runif(p)
## ## beta.gamma <- c(beta,gamma)
## mu.x <- rep(0,p)
## mu.y <- mu.x+0
## Sigma.x <- matrix(rnorm(p^2),nrow=p)
## Sigma.x <- Sigma.x%*%t(Sigma.x)
## Sigma.y <- matrix(rnorm(p^2),nrow=p)
## Sigma.y <- Sigma.y%*%t(Sigma.y)
## mu <- mu.x-mu.y
## Sigma <- Sigma.x+Sigma.y
## ## quad <- as.numeric(beta%*%Sigma%*%beta)
## ## D.prime.beta <- (Sigma%*%beta %*% (beta%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-beta%*%mu/sqrt(quad))
## ## D.prime.gamma <- (Sigma%*%beta %*% (beta%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-beta%*%mu/sqrt(quad))
## ## D.primes <- lapply(list(beta=beta,gamma=gamma), function(coefs) {
##     ## quad <- as.numeric(coefs%*%Sigma%*%coefs)
##     ## (Sigma%*%coefs %*% (coefs%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-coefs%*%mu/sqrt(quad))
## ## })
## ## D.prime <- c(D.primes$beta,D.primes$gamma)
## z.stats <- replicate(1e3, {
##     covariates.x <- rmvnorm(n,mu.x,Sigma)
##     covariates.y <- rmvnorm(n,mu.y,Sigma)
##     covariates <- rbind(covariates.x,covariates.y)
##     true.probs <- plogis(covariates%*%beta)
##     g <- rbinom(length(true.probs),1,prob=true.probs)
##     ## beta.hat <- coef(glm(g ~ covariates-1, family=binomial))
##     ## gamma.hat <- coef(glm(g ~ covariates[,-p]-1, family=binomial))
##     ## xy <- cbind(covariates%*%beta.hat, covariates[,-p]%*%gamma.hat)
##     beta.hat <- beta+rnorm(p)/sqrt(n)
##     gamma.hat <- gamma+rnorm(p)/sqrt(n)
##     xy <- cbind(covariates%*%beta.hat, covariates[,]%*%gamma.hat[])
##     ## x <- rmvnorm(n,mu.x,Sigma.x)
##     ## y <- rmvnorm(n,mu.y,Sigma.y)
##     ## delong.test(x=cbind(x%*%beta.hat,x%*%gamma.hat),y=cbind(y%*%beta.hat,y%*%gamma.hat))
##     delong.test(g=g,xy=xy)
## })
## p.vals <- 1-pnorm(z.stats)
## plot(ecdf(p.vals))
## abline(0,1,col=2)


## abline(0,1,col=2)


## checking root-n consistency of estimators
source('misc.R')
## set.seed(2)
require(mvtnorm)
auc.continuous <- function(x,y)mean(outer(x,y,'<'))
p <- 3
n <-1e2
beta <- runif(p)
gamma <- runif(p)
## gamma <- beta;# gamma[p] <- 0
beta <- gamma <- rep(0,p)#runif(p)
## beta.gamma <- c(beta,gamma)
mu.x <- rep(0,p)
mu.y <- mu.x+1
Sigma.x <- matrix(rnorm(p^2),nrow=p)
Sigma.x <- Sigma.x%*%t(Sigma.x)
Sigma.y <- matrix(rnorm(p^2),nrow=p)
Sigma.y <- Sigma.y%*%t(Sigma.y)
mu <- mu.x-mu.y
Sigma <- Sigma.x+Sigma.y
## quad <- as.numeric(beta%*%Sigma%*%beta)
## D.prime.beta <- (Sigma%*%beta %*% (beta%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-beta%*%mu/sqrt(quad))
## D.prime.gamma <- (Sigma%*%beta %*% (beta%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-beta%*%mu/sqrt(quad))
## D.primes <- lapply(list(beta=beta,gamma=gamma), function(coefs) {
    ## quad <- as.numeric(coefs%*%Sigma%*%coefs)
    ## (Sigma%*%coefs %*% (coefs%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-coefs%*%mu/sqrt(quad))
## })
## D.prime <- c(D.primes$beta,D.primes$gamma)
ns <- round(seq(10,1e2,len=10))
by.n <- sapply(ns, function(n)
    cat('.')
    z.stats <- replicate(1e2, {
        ## beta.hat <- beta+rnorm(p)/sqrt(n)
        ## gamma.hat <- gamma+rnorm(p)/sqrt(n)
        ## beta.gamma.hat <- c(beta.hat,gamma.hat)
        covariates.x <- rmvnorm(n,mu.x,Sigma)
        covariates.y <- rmvnorm(n,mu.y,Sigma)
        covariates <- rbind(covariates.x,covariates.y)
        true.probs <- plogis(covariates%*%beta)
        g <- rbinom(length(true.probs),1,prob=true.probs)
        beta.hat <- coef(glm(g ~ covariates-1, family=binomial))
        sqrt(n)*(beta.hat-beta)
        ## gamma.hat <- coef(glm(g ~ covariates[,-p]-1, family=binomial))
        ## xy <- cbind(covariates%*%beta.hat, covariates[,-p]%*%gamma.hat)
        ## x <- rmvnorm(n,mu.x,Sigma.x)
        ## y <- rmvnorm(n,mu.y,Sigma.y)
        ## delong.test(x=cbind(x%*%beta.hat,x%*%gamma.hat),y=cbind(y%*%beta.hat,y%*%gamma.hat))
        ## delong.test(g=g,xy=xy)
    })
)
matplot(t(by.n),pch='.',col=1,cex=3)

dd

## 4. taylor expansion, normal covariates

## 4a constant term
require(mvtnorm)
p <- 3
n <- 1e2
pairs <- replicate(1e2, {
    beta <- runif(p)
    mu.x <- rep(1,p)
    mu.y <- rep(0,p)
    Sigma.x <- matrix(rnorm(p^2),nrow=p)
    Sigma.x <- Sigma.x%*%t(Sigma.x)
    Sigma.y <- matrix(rnorm(p^2),nrow=p)
    Sigma.y <- Sigma.y%*%t(Sigma.y)
    mu <- mu.x-mu.y
    Sigma <- Sigma.x+Sigma.y
    x <- rmvnorm(n,mu.x,Sigma.x)
    y <- rmvnorm(n,mu.y,Sigma.y)
    obs <- mean((x-y)%*%beta < 0)
    fla <- pnorm(-(beta%*%mu)/sqrt(t(beta)%*%Sigma%*%beta))
    c(obs,fla)
})
plot(pairs[1,],pairs[2,])
abline(0,1)


auc.continuous <- function(x,y)mean(outer(x,y,'<'))
p <- 3
n <- 1e2
beta <- runif(p)
gamma <- runif(p)
err <- runif(p)/1e1
beta.hat <- beta+err
gamma.hat <- gamma+err
mu.x <- rep(1,p)
mu.y <- rep(0,p)
Sigma.x <- matrix(rnorm(p^2),nrow=p)
Sigma.x <- Sigma.x%*%t(Sigma.x)
Sigma.y <- matrix(rnorm(p^2),nrow=p)
Sigma.y <- Sigma.y%*%t(Sigma.y)
mu <- mu.x-mu.y
Sigma <- Sigma.x+Sigma.y
obs <- replicate(1e3, {
    x <- rmvnorm(n,mu.x,Sigma.x)
    y <- rmvnorm(n,mu.y,Sigma.y)
    ## obs <- mean((x-y)%*%beta < 0)
    ## fla <- pnorm(-(beta%*%mu)/sqrt(t(beta)%*%Sigma%*%beta))
    ## auc.beta.hat <- auc.continuous(x%*%beta.hat,y%*%beta.hat)
    aucs <- sapply(list(beta=beta,gamma=gamma,beta.hat=beta.hat,gamma.hat=gamma.hat), function(coefs)auc.continuous(x%*%coefs,y%*%coefs))
    aucs['beta']-aucs['gamma']
})
hist(obs)
abline(v=pnorm(-(beta%*%mu)/sqrt(t(beta)%*%Sigma%*%beta)) - pnorm(-(gamma%*%mu)/sqrt(t(gamma)%*%Sigma%*%gamma)))
abline(v=mean(obs),col=2)

## derivative formula [old]
## p <- 3
## beta <- runif(p)
## mu <- runif(p)
## Sigma <- matrix(rnorm(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
## E.D <- function(beta,mu,Sigma) pnorm(-(beta%*%mu)/sqrt(t(beta)%*%Sigma%*%beta))
## E.D.prime <-  function(coefs,mu,Sigma) {
##     quad <- as.numeric(coefs%*%Sigma%*%coefs)
##     (Sigma%*%coefs %*% (coefs%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-coefs%*%mu/sqrt(quad))
## }
## delta <- runif(p)
## ts <- seq(0.0001,1,len=20)
## newton <- (sapply(ts, function(t)E.D(beta+t*delta,mu,Sigma)) - c(E.D(beta,mu,Sigma))) / ts
## plot(newton)
## abline(h=delta%*%E.D.prime(beta,mu,Sigma))
## ## abline(a=E.D(beta,mu,Sigma),b=delta%*%E.D.prime(beta,mu,Sigma))

## ## -----moved to misc.R
## routine to compute mean auc between linear combn of markers
## source('misc.R')
## ## P(coefs%*%x.0 < coefs%*%x.1) where x.0-x.1 is multivariate normal
## ## with mean mu.diff and variance Sigma.diff
## auc.scores <- function(coefs,mu.diff,Sigma.diff) {
##     ## mu.diff <- c(0,mu.diff)
##     ## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
##     pnorm((coefs%*%mu.diff)/sqrt(t(coefs)%*%(Sigma.diff)%*%coefs))
## }
## p <- 3
## n <- 1e3
## beta <- runif(p+1)
## Sigma <- matrix(rnorm(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
## lda.sampler <- lda.sampler.init(beta,Sigma)
## aucs <- replicate(1e2, {
##     list2env(lda.sampler$sample(n),globalenv())
##     auc.continuous(x.0%*%beta,x.1%*%beta)
## })
## hist(aucs)
## mu.diff <- c(0,lda.sampler$mu.1-lda.sampler$mu.0)
## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
## abline(v=auc.scores(beta,mu.diff=mu.diff,Sigma.diff),col=2)
## abline(v=mean(aucs),col=3)

## ## 4b derivative of auc.scores
## source('misc.R')
## auc.scores.deriv <-  function(coefs,mu.diff,Sigma.diff) {
##     quad <- as.numeric(coefs%*%Sigma.diff%*%coefs)
##     -(Sigma.diff%*%coefs %*% (coefs%*%mu.diff/quad) - mu.diff)%*%(1/sqrt(quad)) %*% dnorm(coefs%*%mu.diff/sqrt(quad))
## }
## p <- 3
## n <- 1e3
## beta <- runif(p)
## Sigma <- matrix(rnorm(p^2),nrow=p)
## Sigma.diff <- Sigma%*%t(Sigma)
## mu.diff <- runif(p)
## delta <- runif(p)
## ts <- seq(0,.1,len=200)
## means <- sapply(ts, function(t)auc.scores(coefs=beta+t*delta,mu.diff=mu.diff,Sigma.diff=Sigma.diff))
## plot(ts,means,type='l')
## slope <- delta%*%auc.scores.deriv(coefs=beta,mu.diff=mu.diff,Sigma.diff=Sigma.diff)
## abline(a=means[1],b=slope,col=2)

## ## ---------moved to misc.R
## ## "standard situation" sampler--check lda.sammpler.init.std. written as a closure
## source('misc.R')
## auc.scores <- function(coefs,mu.diff,Sigma.diff) {
##     ## mu.diff <- c(0,mu.diff)
##     ## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
##     pnorm((coefs%*%mu.diff)/sqrt(t(coefs)%*%(Sigma.diff)%*%coefs))
## }
## p <- 3
## n <- 1e3
## Delta <- runif(1)
## lambda <- rlogis(1)
## list2env(lda.sampler.init.std(Delta,lambda,p),globalenv())
## aucs <- replicate(1e2, {
##     list2env(sample(n),globalenv())
##     auc.hat.continuous(x.0%*%beta,x.1%*%beta)
## })
## hist(aucs)
## ## mu.diff <- with(lda.sampler,c(0,mu.0-mu.1))
## ## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
## abline(v=auc.scores(beta,mu.diff=mu.diff,Sigma.diff),col=2)
## abline(v=mean(aucs),col=3)


## 4c derivative in lda case: slope is 0. when x.0 and x.1 are mv normal
## then P(beta^t x.0 < beta^t x.1) has a stationary point when beta
## are the lda coefficients.
source('misc.R')
p <- 3
n <- 1e3
beta <- runif(p+1)
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
lda.sampler <- lda.sampler.init(beta,Sigma)
delta <- runif(p+1)
ts <- seq(0,.01,len=200)
mu.diff <- c(0,lda.sampler$mu.1-lda.sampler$mu.0)
Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
means <- sapply(ts, function(t)auc.scores(coefs=beta+t*delta,mu.diff=mu.diff,Sigma.diff=Sigma.diff))
plot(ts,means,type='l')
deriv.1 <- delta%*%auc.scores.deriv(coefs=beta,mu.diff=mu.diff,Sigma.diff=Sigma.diff)
curve(means[1]+x*deriv.1,add=TRUE,col=2)      
## abline(a=means[1],b=slope,col=2)
deriv.2 <- t(delta)%*%auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)%*%delta
curve(means[1]+x*deriv.1+x^2/2*deriv.2,add=TRUE,col=2)      


source('misc.R')
p <- 5
n <- 1e2
beta <- runif(p+1)
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
lda.sampler <- lda.sampler.init(beta,Sigma)
delta <- runif(p+1)
mu.diff <- c(0,lda.sampler$mu.1-lda.sampler$mu.0)
Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
auc.scores.deriv(coefs=beta,mu.diff=mu.diff,Sigma.diff=Sigma.diff)

coefs <- beta
quad <- as.numeric(coefs%*%Sigma.diff%*%coefs)
(Sigma.diff%*%coefs %*% (coefs%*%mu.diff/quad) - mu.diff)%*%(1/sqrt(quad)) %*% dnorm(-coefs%*%mu.diff/sqrt(quad))
Sigma.diff%*%coefs %*% (coefs%*%mu.diff/quad) - mu.diff
(beta%*%mu.diff)/(beta[-1]%*%Sigma%*%beta[-1])



## what about when x.0 and x.1 aren't mv normal? not clear but seems possible from this example.
require(mvtnorm)
source('misc.R')
set.seed(1)
p <- 5
n <- 1e6
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
mu.0 <- runif(p)
mu.1 <- runif(p)
x.0 <- cbind(1,pnorm(rmvnorm(n,mu.0,Sigma)))
x.1 <- cbind(1,pnorm(rmvnorm(n,mu.1,Sigma)))
beta.lda <- lda.coefs(x.0,x.1)

n <- 1e3
delta <- runif(p+1)
ts <- seq(-3,3,len=1e2)
aucs <- sapply(ts, function(t) {
    beta <- beta.lda + delta*t
    x.0 <- cbind(1,pnorm(rmvnorm(n,mu.0,Sigma)))
    x.1 <- cbind(1,pnorm(rmvnorm(n,mu.1,Sigma)))
    auc.hat(x.0%*%beta, x.1%*%beta)
})
plot(ts,aucs,type='l')
abline(v=0,col=2)


## require(sdprisk)
source('misc.R')
p <- 5
n <- 3e3
mu.0 <- runif(p)
mu.1 <- runif(p)
mu.0 <- mu.0/sqrt(sum(mu.0^2))
mu.1 <- mu.1/sqrt(sum(mu.1^2))
beta.lda <- c(0,mu.1-mu.0)
## x.1 <- cbind(1,mu.1+matrix(rexp(n*p)-1,ncol=p))
## lda.coefs(x.0,x.1) - c(1,beta.lda)
## checking hypoexpo function
## beta <- runif(p)
## x <- matrix(rexp(n*p),ncol=p)%*%beta
## plot(ecdf(x))
## curve(phypoexp(x,rate=1/beta),add=TRUE,col=2)
delta <- c(0,runif(p))
ts <- seq(-3,3,len=1e2)
aucs <- sapply(ts, function(t) {
    beta <- (beta.lda+delta*t)[-1]
    ## browser()
    ## print(beta)
    ## print(sum(beta<0))
    ## P(beta*expo.0-beta*expo.1 < mu.1-mu.0)
    ## x.0 <- rhypoexp(n,rate=1/beta)
    ## x.1 <- rhypoexp(n,rate=1/beta)
    x.0 <- t(mu.0+t(matrix(rexp(n*p)-1,ncol=p)))%*%beta
    x.1 <- t(mu.1+t(matrix(rexp(n*p)-1,ncol=p)))%*%beta
    auc.hat(beta*mu.0+x.0, beta*mu.1+x.1)
})
plot(ts,aucs,type='l')
abline(v=0)

dd
## x.0 <- cbind(1,t(mu.0+t(matrix(rexp(n*p)-1,ncol=p))))
## x.1 <- cbind(1,t(mu.1+t(matrix(rexp(n*p)-1,ncol=p))))


## 4d formula for second derivative, general case and "standard situation"
auc.scores.deriv2.lda <- function(coefs,Sigma.diff) {
      quad <- as.numeric(coefs%*%Sigma.diff%*%coefs)
      -dnorm(sqrt(quad)/2)/sqrt(quad)*(Sigma.diff - (Sigma.diff%*%coefs)%*%t(Sigma.diff%*%coefs)/quad)/(-2)
      }
auc.scores.deriv2.lda.std <- function(Delta,p) {
      ## quad <- as.numeric(coefs%*%Sigma.diff%*%coefs)
    ## dnorm(sqrt(quad)/2)/sqrt(quad)*(Sigma.diff - (Sigma.diff%*%coefs)%*%t(Sigma.diff%*%coefs)/quad)/(-2)
    E.11 <- matrix(0,p,p); E.11[1,1] <- 1
    -rbind(0,cbind(0,dnorm(Delta/sqrt(2))*1/Delta/sqrt(2)*(diag(p)-E.11)))
      }
source('misc.R')
p <- 3
n <- 1e3
## beta <- runif(p+1)
## Sigma <- matrix(rnorm(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
## lda.sampler <- lda.sampler.init(beta,Sigma)
Delta <- runif(1)
lambda <- rlogis(1)
list2env(lda.sampler.init.std(Delta,lambda,p),globalenv())
auc.scores.deriv2.lda(beta,Sigma.diff)-auc.scores.deriv2.lda.std(Delta,p)
delta <- runif(p+1)
a <- runif(p+1)
ts <- seq(0,1,len=2000)
## mu.diff <- c(0,lda.sampler$mu.0-lda.sampler$mu.1)
## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
## means <- sapply(ts, function(t)auc.scores(coefs=beta+t*delta,mu.diff=mu.diff,Sigma.diff=Sigma.diff))
## plot(ts,means,type='l')
derivs <- sapply(ts, function(t)a%*%auc.scores.deriv(coefs=beta+t*delta,mu.diff=mu.diff,Sigma.diff=Sigma.diff))
plot(ts,derivs,type='l')
## curve(means[1]+x*deriv.1,add=TRUE,col=2)      
## abline(a=means[1],b=slope,col=2)
deriv.2 <- t(delta)%*%auc.scores.deriv2.lda.std(Delta,p)%*%a
curve(derivs[1]+x*deriv.2,add=TRUE,col=2)      





dd

## lda.sampler <- lda.sampler.init(beta,Sigma)
Delta <- runif(1)
lambda <- rlogis(1)
list2env(lda.sampler.init.std(Delta,lambda,p),globalenv())
auc.scores.deriv2.lda(beta,Sigma.diff)-auc.scores.deriv2.lda.std(Delta,p)

dd


## 4e taylor apprx. first derivative is 0 when expanding around the lda
## coefs, using second deriv.
source('misc.R')
p <- 3
n <- 3e2
pairs <- replicate(1e2, {
    tryCatch({
    beta <- runif(p+1)
    Sigma <- matrix(rnorm(p^2),nrow=p)
    Sigma <- Sigma%*%t(Sigma)
    list2env(lda.sampler.init(beta,Sigma),globalenv())
    mu.diff <- c(0,mu.1-mu.0)
    Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
    list2env(sample(n),globalenv())
    ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
    beta.hat <- lda.coefs(x.0,x.1)
    auc.beta <- auc.scores(beta,mu.diff,Sigma.diff)
    auc.hat.beta <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
    auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
    ## c(auc.hat.beta,auc.hat.beta.hat)
    ## auc.scores(coefs=beta,mu.diff,Sigma.diff)
    deriv.2 <- auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)
    observed <- auc.hat.beta.hat
    taylor <- auc.hat.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    ## observed <- auc.hat.beta.hat - auc.hat.beta
    ## taylor <- t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    ## observed <- auc.hat.beta.hat - auc.beta
    ## taylor <- auc.hat.beta-auc.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    c(observed=observed,taylor=taylor)
    },error=function(e)c(NA,NA))
})
plot(pairs['observed',],pairs['taylor',]); abline(0,1,col=2)

## "standard situation"
source('misc.R')
p <- 3
n <- 2e3
pairs <- replicate(1e2, {
    ## tryCatch({
    ## beta <- runif(p+1)
    ## Sigma <- matrix(rnorm(p^2),nrow=p)
    ## Sigma <- Sigma%*%t(Sigma)
    Delta <- runif(1)
    lambda <- rlogis(1)
    list2env(lda.sampler.init.std(Delta,lambda,p),globalenv())
    ## mu.diff <- c(0,mu.1-mu.0)
    ## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
    list2env(sample(n),globalenv())
    ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
    beta.hat <- lda.coefs(x.0,x.1)
    auc.beta <- auc.scores(beta,mu.diff,Sigma.diff)
    auc.hat.beta <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
    auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
    ## c(auc.hat.beta,auc.hat.beta.hat)
    ## auc.scores(coefs=beta,mu.diff,Sigma.diff)
    deriv.2 <- auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)
    observed <- auc.hat.beta.hat
    taylor <- auc.hat.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    ## observed <- auc.hat.beta.hat - auc.hat.beta
    ## taylor <- t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    ## observed <- auc.hat.beta.hat - auc.beta
    ## taylor <- auc.hat.beta-auc.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    c(observed=observed,taylor=taylor)
    ## },error=function(e)c(NA,NA))
})
plot(pairs['observed',],pairs['taylor',]); abline(0,1,col=2)



## 4f asy behavior. Since #1 n*(auc.hat(beta.hat)-auc.hat(beta)) is
## O_p(1), then #2 sqrt(n)*(auc.hat(beta.hat)-auc(beta)) is within
## o_p(1) of sqrt(n)*(auc.hat(beta)-auc(beta)). beta are the true lda
## coefs. 
source('misc.R')
p <- 3
n <- 3e2
ns <- round(seq(1e2,1e3,len=1e2))
by.n <- sapply(ns, function(n) {
    terms <- replicate(2e1, {
        tryCatch({
            beta <- runif(p+1)
            Sigma <- matrix(rnorm(p^2),nrow=p)
            Sigma <- Sigma%*%t(Sigma)
            list2env(lda.sampler.init(beta,Sigma),globalenv())
            mu.diff <- c(0,mu.1-mu.0)
            Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
            list2env(sample(n),globalenv())
            ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
            beta.hat <- lda.coefs(x.0,x.1)
            auc.beta <- auc.scores(beta,mu.diff,Sigma.diff)
            auc.beta.hat <- auc.scores(beta.hat,mu.diff,Sigma.diff)
            auc.hat.beta <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
            auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
            ## c(auc.hat.beta,auc.hat.beta.hat)
            ## auc.scores(coefs=beta,mu.diff,Sigma.diff)
            deriv.2 <- auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)
            ## observed <- auc.hat.beta.hat
            ## taylor <- auc.hat.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
            ## observed <- auc.hat.beta.hat - auc.hat.beta
            ## taylor <- t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
            ## observed <- auc.hat.beta.hat - auc.hat.beta
            ## taylor <- auc.hat.beta-auc.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
            ## c(auc.beta.hat-auc.beta,auc.hat.beta.hat - auc.hat.beta) #1
            c(auc.hat.beta.hat-auc.beta,auc.hat.beta - auc.beta) #2
        },error=function(e)NA)
    },simplify=FALSE)
    terms <- do.call(cbind,terms)
    rowMeans(terms,na.rm=TRUE)
})
n.power <- ns^(1/2)
plot(ns,n.power*by.n[1,],type='l')
lines(ns,n.power*by.n[2,],col=2)
abline(h=0)


## #1 is the general lda parameterization, #2 is the efron standard situation
source('misc.R')
set.seed(1)
p <- 3
n <- 3e3
## ns <- round(seq(1e2,1e3,len=1e2))
## by.n <- sapply(ns, function(n) {
## terms <- replicate(2e1, {
## tryCatch({
## beta <- runif(p+1)/4 #1
## Sigma <- matrix(rnorm(p^2),nrow=p) #1
## Sigma <- Sigma%*%t(Sigma) #1
## beta.sampler <- lda.sampler.init(beta,Sigma)    #1
Delta <- runif(1) #2
lambda <- rlogis(1) #2
beta.sampler <- lda.sampler.init.std(Delta=Delta,lambda,p) #2
beta <- beta.sampler$beta #2
gamma <- beta; gamma[p+1] <- 0
pairs <- replicate(1e2, {
    ## start <- Sys.time()
    list2env(beta.sampler$sample(n),globalenv())
    ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
    beta.hat <- lda.coefs(x.0,x.1)
    gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
    gamma.hat <- c(gamma.hat,0)
    ## scores.0 <- cbind(x.0%*%beta.hat,x.0%*%gamma.hat)
    ## scores.1 <- cbind(x.1%*%beta.hat,x.1%*%gamma.hat)
    ## var.diff <- delong.var(x=scores.0,y=scores.1)
    auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
    auc.hat.gamma.hat <- auc.hat.continuous(x.0%*%gamma.hat,x.1%*%gamma.hat)
    auc.hat.beta <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
    auc.hat.gamma <- auc.hat.continuous(x.0%*%gamma,x.1%*%gamma)
    ## c(auc.hat.beta.hat - auc.hat.gamma.hat,auc.hat.beta - auc.hat.gamma)
    c(auc.hat.beta.hat ,auc.hat.beta )
    ## c(auc.hat.gamma.hat ,auc.hat.gamma )
    ## Sys.time()-start
})
plot(pairs[1,],pairs[2,]); abline(0,1)



## 4ff ok to ignore parameter estimates when beta\neq gamma so maybe
## delong variance holds asymptotically [outcome: no, the first
## statement only true if beta is the lda coef, and if \beta\neq\gamma
## one of them isn't the lda coef, in which case the lienar term in
## the taylor exp needs to be accounted for
delong.var <- function(x=NULL,y=NULL,xy=NULL,g=NULL) {
    if(is.null(x)||is.null(y)) {
        if(is.null(xy)||is.null(g)) stop()
        g <- factor(g,labels=0:1)
        x <- xy[g==0,]; y <- xy[g==1,]
    }
    m <- nrow(x); n <- nrow(y)
    theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
    V.10.1 <- sapply(x[,1], function(x.i)mean(x.i<y[,1]))
    V.10.2 <- sapply(x[,2], function(x.i)mean(x.i<y[,2]))
    V.01.1 <- sapply(y[,1], function(y.i)mean(x[,1]<y.i))
    V.01.2 <- sapply(y[,2], function(y.i)mean(x[,2]<y.i))
    S.10 <- cov(cbind(V.10.1,V.10.2))
    S.01 <- cov(cbind(V.01.1,V.01.2))
    S <- S.10/m + S.01/n
    contrast <- matrix(c(1,-1),ncol=1)
    ## z.stat <- diff(rev(theta.hats)) / sqrt(t(contrast)%*%S%*%contrast)
    t(contrast)%*%S%*%contrast
    ## return(z.stat)
}
source('misc.R')
p <- 3
n <- 1e3
## ns <- round(seq(1e2,1e3,len=1e2))
## by.n <- sapply(ns, function(n) {
## terms <- replicate(2e1, {
## tryCatch({
beta <- runif(p+1)/4
## gamma <- beta; gamma[p+1] <- 0
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta.sampler <- lda.sampler.init(beta,Sigma)        
## mu.diff <- c(0,mu.0-mu.1)
## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
z.stats <- replicate(1e2, {
    list2env(beta.sampler$sample(n),globalenv())
    ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
    beta.hat <- lda.coefs(x.0,x.1)
    gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
    gamma.hat <- c(gamma.hat,0)
    scores.0 <- cbind(x.0%*%beta.hat,x.0%*%gamma.hat)
    scores.1 <- cbind(x.1%*%beta.hat,x.1%*%gamma.hat)
    var.diff <- delong.var(x=scores.0,y=scores.1)
    ## auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
    ## auc.hat.gamma.hat <- auc.hat.continuous(x.0%*%gamma.hat,x.1%*%gamma.hat)
    gamma <- beta; gamma[p+1] <- 0
    auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
    auc.hat.gamma.hat <- auc.hat.continuous(x.0%*%gamma,x.1%*%gamma)
    (auc.hat.beta.hat - auc.hat.gamma.hat) / sqrt(var.diff)
})
qqnorm(z.stats); abline(0,1)



## 4g diff of auc.hats, when beta==gamma. both constant term and linear
## term of taylor poly =0
source('misc.R')
## set.seed(1)
p <- 3
n <- 1e3
ns <- round(seq(1e2,1e3,len=1e1))
## by.n <- sapply(ns, function(n) {
## terms <- replicate(2e1, {
## tryCatch({
beta <- runif(p+1)/4
beta[p+1] <- 0 
## gamma <- beta; gamma[p+1] <- 0
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta.sampler <- lda.sampler.init(beta,Sigma)        
## mu.diff <- c(0,mu.0-mu.1)
## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
by.n <- sapply(ns, function(n) {
    pairs <- replicate(1e2, {
        list2env(beta.sampler$sample(n),globalenv())
        ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
        beta.hat <- lda.coefs(x.0,x.1)
        gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
        gamma.hat <- c(gamma.hat,0)
        ## scores.0 <- cbind(x.0%*%beta.hat,x.0%*%gamma.hat)
        ## scores.1 <- cbind(x.1%*%beta.hat,x.1%*%gamma.hat)
        ## var.diff <- delong.var(x=scores.0,y=scores.1)
        auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat.continuous(x.0%*%gamma.hat,x.1%*%gamma.hat)
        gamma <- beta; gamma[p+1] <- 0
        auc.hat.beta <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
        auc.hat.gamma <- auc.hat.continuous(x.0%*%gamma,x.1%*%gamma)
        ## c(auc.hat.beta.hat - auc.hat.gamma.hat,auc.hat.beta - auc.hat.gamma)
        ## c(auc.hat.beta.hat ,auc.hat.beta )
        ## c(auc.hat.gamma.hat ,auc.hat.gamma )
        ## c(auc.hat.beta.hat-auc.hat.gamma.hat,auc.hat.beta-auc.hat.gamma )
        auc.hat.beta.hat-auc.hat.gamma.hat
    })
})
matplot(t(by.n),pch=1,col=1);abline(h=0)

## diff of auc.hats, when beta==gamma.
source('misc.R')
## set.seed(1)
p <- 3
n <- 3e2
ns <- round(seq(1e2,1e3,len=1e1))
## by.n <- sapply(ns, function(n) {
## terms <- replicate(2e1, {
## tryCatch({
## parts <- replicate(3e2, {
beta <- runif(p+1)/4
beta[p+1] <- 0 
gamma <- beta
## gamma <- beta; gamma[p+1] <- 0
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta.sampler <- lda.sampler.init(beta,Sigma)        
mu.diff <- with(beta.sampler, c(0,mu.1-mu.0))
Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
## Delta <- runif(1) #2
## lambda <- rlogis(1) #2
## beta.sampler <- lda.sampler.init.std(Delta=Delta,lambda,p) #2
## beta <- beta.sampler$beta #2
## gamma <- beta; gamma[p+1] <- 0 #2
## by.n <- sapply(ns, function(n) {
parts <- replicate(1e2, {
    list2env(beta.sampler$sample(n),globalenv())
    ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
    ## gamma.hat <- gamma+rnorm(p+1)/sqrt(n)
    beta.hat <- lda.coefs(x.0,x.1)
    gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
    gamma.hat <- c(gamma.hat,0)
    auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
    auc.hat.gamma.hat <- auc.hat.continuous(x.0%*%gamma.hat,x.1%*%gamma.hat)
    auc.hat.beta <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
    auc.hat.gamma <- auc.hat.continuous(x.0%*%gamma,x.1%*%gamma) #==auc.hat.beta
    deriv.2 <- auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)
    ## observed <- auc.hat.beta.hat
    ## taylor <- auc.hat.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    ## observed <- auc.hat.beta.hat - auc.hat.beta
    ## taylor <-  t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    ## observed <- auc.hat.gamma.hat
    ## taylor <- auc.hat.gamma + t(gamma.hat-gamma)%*%deriv.2%*%(gamma.hat-gamma)/2
    ## observed <- auc.hat.beta.hat - auc.hat.gamma.hat
    ## taylor <-  t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2 - t(gamma.hat-gamma)%*%deriv.2%*%(gamma.hat-gamma)/2 
    ## c(observed,taylor)
    hat.beta <- auc.hat.beta.hat; hat.gamma <- auc.hat.gamma.hat
    const.beta <- auc.hat.beta; const.gamma <- auc.hat.gamma
    quad.beta <- t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    quad.gamma <- t(gamma.hat-gamma)%*%deriv.2%*%(gamma.hat-gamma)/2
    c(hat.beta=hat.beta,hat.gamma=hat.gamma,const.beta=const.beta,const.gamma=const.gamma,quad.beta=quad.beta,quad.gamma=quad.gamma)
})
## })
## matplot(t(by.n),pch=1,col=1);abline(h=0)
## plot(parts[1,],parts[2,]); abline(0,1)
## plot(pairs[1,]-pairs[2,]);abline(h=0)
parts <- as.data.frame(t(parts))
with(parts,plot(hat.beta,const.beta+quad.beta))
with(parts,plot(hat.gamma,const.beta+quad.gamma));abline(0,1)
with(parts,plot(hat.beta+hat.gamma,2*const.beta+quad.beta+quad.gamma));abline(0,1)
with(parts,plot(hat.beta-hat.gamma,quad.beta-quad.gamma));abline(0,1)
boxplot(parts)

dd

source('misc.R')
set.seed(1)
p <- 3
n <- 5e2
beta <- runif(p+1)
beta[p+1] <- 0
gamma <- beta
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta.sampler <- lda.sampler.init(beta,Sigma)        
pairs <- replicate(1e2, {
    ##     tryCatch({
    ##     beta <- runif(p+1)
    beta.hat <- beta+rnorm(p+1)/sqrt(n)
    ##     Sigma <- matrix(rnorm(p^2),nrow=p)
    ##     Sigma <- Sigma%*%t(Sigma)
    list2env(lda.sampler.init(beta,Sigma),globalenv())
    mu.diff <- c(0,mu.1-mu.0)
    Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
    list2env(sample(n),globalenv())
    beta.hat <- lda.coefs(x.0,x.1)
    ## gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
    ## gamma.hat <- c(gamma.hat,0)
    auc.beta <- auc.scores(beta,mu.diff,Sigma.diff)
    auc.beta.hat <- auc.scores(beta.hat,mu.diff,Sigma.diff)
    auc.hat.beta <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
    auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
    ## c(auc.hat.beta,auc.hat.beta.hat)
    ## auc.scores(coefs=beta,mu.diff,Sigma.diff)
    deriv.2 <- auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)
    observed <- auc.hat.beta.hat
    taylor <- auc.hat.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    observed <- auc.hat.beta.hat - auc.hat.beta
    taylor <-  -t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    ## observed <- auc.hat.beta.hat - auc.scores(beta,mu.diff,Sigma.diff)#auc.hat.beta
    ## observed <- auc.hat.beta.hat - auc.beta
    ## taylor <- auc.hat.beta-auc.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    ##----
    ## observed <- auc.hat.beta.hat - auc.hat.beta
    ## observed <- auc.beta.hat - auc.beta
    ## taylor <- t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    ## observed <- auc.hat.beta.hat
    ## taylor <- auc.beta.hat
    ## observed <- auc.hat.beta
    ## taylor <- auc.beta
    ## observed <-  auc.hat.beta.hat - auc.hat.beta 
    ## taylor <- auc.beta.hat - auc.beta
    ## observed <-  auc.hat.beta.hat - auc.beta.hat
    ## taylor <-  auc.hat.beta  - auc.beta
    c(observed=observed,taylor=taylor)
    ## },error=function(e)c(NA,NA))
})
plot(pairs['observed',],pairs['taylor',]); abline(0,1,col=2)





source('misc.R')
## set.seed(1)
p <- 3
n <- 1e3
beta <- runif(p+1)
beta[p+1] <- 0
gamma <- beta
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta.sampler <- list2env(lda.sampler.init(beta,Sigma),globalenv())
parts <- replicate(1e2, {
    ##     tryCatch({
    ##     beta <- runif(p+1)
    beta.hat <- beta+rnorm(p+1)/sqrt(n)
    ##     Sigma <- matrix(rnorm(p^2),nrow=p)
    ##     Sigma <- Sigma%*%t(Sigma)
    ## list2env(lda.sampler.init(beta,Sigma),globalenv())
    mu.diff <- c(0,mu.1-mu.0)
    Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
    list2env(sample(n),globalenv())
    beta.hat <- lda.coefs(x.0,x.1)
    gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
    gamma.hat <- c(gamma.hat,0)
    auc.beta <- auc.scores(beta,mu.diff,Sigma.diff)
    auc.beta.hat <- auc.scores(beta.hat,mu.diff,Sigma.diff)
    auc.hat.beta <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
    auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
    ## c(auc.hat.beta,auc.hat.beta.hat)
    ## auc.scores(coefs=beta,mu.diff,Sigma.diff)
    deriv.2 <- auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)
    ## observed <- auc.hat.beta.hat
    ## taylor <- auc.hat.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    ## observed <- auc.hat.beta.hat - auc.scores(beta,mu.diff,Sigma.diff)#auc.hat.beta
    ## observed <- auc.hat.beta.hat - auc.beta
    ## taylor <- auc.hat.beta-auc.beta + t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    ##----
    ## observed <- auc.hat.beta.hat - auc.hat.beta
    ## observed <- auc.beta.hat - auc.beta
    ## taylor <- t(beta.hat-beta)%*%deriv.2%*%(beta.hat-beta)/2
    observed <- auc.hat.beta.hat
    taylor <- auc.beta.hat
    observed <- auc.hat.beta
    taylor <- auc.beta
    observed <- auc.beta.hat - auc.beta
    taylor <- auc.hat.beta.hat - auc.hat.beta
    ## c(observed=observed,taylor=taylor)
    c(auc.hat.beta.hat=auc.hat.beta.hat,auc.beta.hat=auc.beta.hat,auc.hat.beta=auc.hat.beta,auc.beta=auc.beta)
    ## },error=function(e)c(NA,NA))
})
op <- par(mfrow=c(1,3))
plot(parts['auc.hat.beta.hat',],parts['auc.beta.hat',]); abline(0,1,col=2)
plot(parts['auc.hat.beta',],parts['auc.beta',]); abline(0,1,col=2)
plot(parts['auc.hat.beta.hat',] - parts['auc.hat.beta',],parts['auc.beta.hat',]-parts['auc.beta',]);abline(0,1)
par(op)

plot(parts['auc.hat.beta.hat',],parts['auc.hat.beta',]); abline(0,1,col=2)
plot(parts['auc.beta.hat',],parts['auc.beta',]); abline(0,1,col=2)



plot(parts['auc.beta.hat',]-parts['auc.beta',],parts['auc.hat.beta.hat',] - parts['auc.hat.beta',]);abline(0,1)

boxplot(t(parts))

op <- par(mfrow=c(1,3))
plot(parts['auc.hat.beta.hat',],parts['auc.beta.hat',],xlim=c(0,1),ylim=c(0,1)); abline(0,1,col=2)
plot(parts['auc.hat.beta',],parts['auc.beta',],xlim=c(0,1),ylim=c(0,1)); abline(0,1,col=2)
par(op)


op <- par(mfrow=c(1,2))
plot(parts['auc.hat.beta.hat',],parts['auc.hat.beta',]); abline(0,1,col=2)
plot(parts['auc.beta.hat',],parts['auc.beta',]); abline(0,1,col=2)
## plot(parts['auc.hat.beta.hat',] - parts['auc.hat.beta',],parts['auc.beta.hat',]-parts['auc.beta',]);abline(0,1)
par(op)

plot(parts['auc.hat.beta.hat',],parts['auc.beta.hat',]); abline(0,1,col=2)
points(parts['auc.hat.beta',],parts['auc.beta',],col=2)

op <- par(mfrow=c(1,2))
hist(parts['auc.hat.beta.hat',] - parts['auc.hat.beta',])
hist(parts['auc.beta.hat',]-parts['auc.beta',])
par(op)

cor(parts['auc.hat.beta.hat',],parts['auc.beta.hat',])
plot(parts['auc.hat.beta',],parts['auc.beta',])



## was debugging...taylor expansion seems to be negation of what is
## observed turns out auc.hat.beta.hat - auc.hat.beta approximates
## negation of auc.beta.hat - auc.beta. The approximation error in the
## auc.hat is stronger than the approximation in the beta.hat. trying
## hajek projection to be more careful.

## diff when gamma==beta
source('misc.R')
set.seed(1)
p <- 3
n <- 5e2
beta <- runif(p+1)
beta[p+1] <- 0
gamma <- beta
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta.sampler <- lda.sampler.init(beta,Sigma)        
pairs <- replicate(1e2, {
    ##     tryCatch({
    ##     beta <- runif(p+1)
    beta.hat <- beta+rnorm(p+1)/sqrt(n)
    ##     Sigma <- matrix(rnorm(p^2),nrow=p)
    ##     Sigma <- Sigma%*%t(Sigma)
    list2env(lda.sampler.init(beta,Sigma),globalenv())
    mu.diff <- c(0,mu.1-mu.0)
    Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
    list2env(sample(n),globalenv())
    beta.hat <- lda.coefs(x.0,x.1)
    gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
    gamma.hat <- c(gamma.hat,0)
    auc.beta <- auc.scores(beta,mu.diff,Sigma.diff)
    auc.beta.hat <- auc.scores(beta.hat,mu.diff,Sigma.diff)
    auc.gamma.hat <- auc.scores(gamma.hat,mu.diff,Sigma.diff)
    auc.hat.beta <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
    auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
    auc.hat.gamma.hat <- auc.hat.continuous(x.0%*%gamma.hat,x.1%*%gamma.hat)
    ## deriv.2 <- auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)
    ## observed <-  auc.hat.beta.hat - auc.beta.hat
    ## taylor <-  auc.hat.beta  - auc.beta
    ## observed <-  auc.hat.beta.hat  - auc.hat.beta
    ## taylor <-  auc.beta - auc.beta.hat
    ## observed <-  auc.hat.beta.hat - auc.beta
    ## taylor <-  auc.hat.beta - auc.beta.hat    
    ## observed <-  auc.hat.beta.hat - auc.beta.hat - auc.hat.beta #+ auc.beta.hat
    ## taylor <-  auc.hat.beta  - auc.beta - auc.hat.beta #+ auc.beta.hat
    observed <- auc.hat.beta.hat - auc.hat.gamma.hat
    taylor <- auc.beta.hat - auc.gamma.hat
    c(observed=observed,taylor=taylor)
    ## },error=function(e)c(NA,NA))
})
plot(pairs['observed',],pairs['taylor',]); abline(0,1,col=2)


## diff when gamma==beta. growth rate.
source('misc.R')
set.seed(1)
p <- 3
n <- 5e2
beta <- runif(p+1)
beta[p+1] <- 0
gamma <- beta
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
## beta.sampler <- lda.sampler.init(beta,Sigma)
list2env(lda.sampler.init(beta,Sigma),globalenv())
ns <- round(seq(1e2,1e3,len=10))
by.n <- sapply(ns, function (n) {
    pairs <- replicate(1e2, {
        ##     tryCatch({
        ##     beta <- runif(p+1)
        ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
        ## list2env(lda.sampler.init(beta,Sigma),globalenv())
        mu.diff <- c(0,mu.1-mu.0)
        Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
        list2env(sample(n),globalenv())
        beta.hat <- lda.coefs(x.0,x.1)
        gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
        gamma.hat <- c(gamma.hat,0)
        auc.beta <- auc.scores(beta,mu.diff,Sigma.diff)
        auc.beta.hat <- auc.scores(beta.hat,mu.diff,Sigma.diff)
        auc.gamma.hat <- auc.scores(gamma.hat,mu.diff,Sigma.diff)
        auc.hat.beta <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
        auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat.continuous(x.0%*%gamma.hat,x.1%*%gamma.hat)
        ## deriv.2 <- auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)
        ## observed <-  auc.hat.beta.hat - auc.beta.hat
        ## taylor <-  auc.hat.beta  - auc.beta
        ## observed <-  auc.hat.beta.hat  - auc.hat.beta
        ## taylor <-  auc.beta - auc.beta.hat
        ## observed <-  auc.hat.beta.hat - auc.beta
        ## taylor <-  auc.hat.beta - auc.beta.hat    
        ## observed <-  auc.hat.beta.hat - auc.beta.hat - auc.hat.beta #+ auc.beta.hat
        ## taylor <-  auc.hat.beta  - auc.beta - auc.hat.beta #+ auc.beta.hat
        observed <- auc.hat.beta.hat - auc.hat.gamma.hat
        taylor <- auc.beta.hat - auc.gamma.hat
        c(observed=observed,taylor=taylor)
        ## },error=function(e)c(NA,NA))
    })
    pairs['observed',]-pairs['taylor',]
})
matplot(ns*t(by.n),pch=1,col=1)



dd

## ## taylor [old]
## set.seed(2)
## require(mvtnorm)
## auc.continuous <- function(x,y)mean(outer(x,y,'<'))
## p <- 3
## n <-1e2
## beta <- runif(p)
## gamma <- runif(p)
## gamma <- beta;# gamma[p] <- 0
## beta.gamma <- c(beta,gamma)
## mu.x <- rep(1,p)
## mu.y <- rep(0,p)
## Sigma.x <- matrix(rnorm(p^2),nrow=p)
## Sigma.x <- Sigma.x%*%t(Sigma.x)
## Sigma.y <- matrix(rnorm(p^2),nrow=p)
## Sigma.y <- Sigma.y%*%t(Sigma.y)
## mu <- mu.x-mu.y
## Sigma <- Sigma.x+Sigma.y
## ## quad <- as.numeric(beta%*%Sigma%*%beta)
## ## D.prime.beta <- (Sigma%*%beta %*% (beta%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-beta%*%mu/sqrt(quad))
## ## D.prime.gamma <- (Sigma%*%beta %*% (beta%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-beta%*%mu/sqrt(quad))
## D.primes <- lapply(list(beta=beta,gamma=gamma), function(coefs) {
##     quad <- as.numeric(coefs%*%Sigma%*%coefs)
##     (Sigma%*%coefs %*% (coefs%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-coefs%*%mu/sqrt(quad))
## })
## D.prime <- c(D.primes$beta,D.primes$gamma)
## pairs <- replicate(1e3, {
##     beta.hat <- beta+rnorm(p)/sqrt(n)
##     gamma.hat <- gamma+rnorm(p)/sqrt(n)
##     beta.gamma.hat <- c(beta.hat,gamma.hat)
##     ## beta.gamma.hat <- c(beta+rnorm(p)/sqrt(n),gamma+rnorm(p)/sqrt(n))
##     ## beta.gamma.hat <- beta.gamma+rnorm(2*p)/sqrt(n)
##     ## delta <- (beta.gamma.hat-beta.gamma)
##     x <- rmvnorm(n,mu.x,Sigma.x)
##     y <- rmvnorm(n,mu.y,Sigma.y)
##     ## obs <- mean((x-y)%*%beta < 0)
##     ## fla <- pnorm(-(beta%*%mu)/sqrt(t(beta)%*%Sigma%*%beta))
##     ## auc.beta.hat <- auc.continuous(x%*%beta.hat,y%*%beta.hat)
##     aucs <- sapply(list(beta.hat=beta.hat,gamma.hat=gamma.hat), function(coefs)auc.continuous(x%*%coefs,y%*%coefs))
##     ## D <- aucs['beta']-aucs['gamma']
##     D.hat <- unname(aucs['beta.hat']-aucs['gamma.hat']    )
##     taylor <- D.hat + (beta.gamma.hat-beta.gamma)%*%D.prime
##     ## taylor <- D.hat + delta%*%D.prime
##     ## taylor <- D.hat + (rnorm(2*p)/sqrt(n))%*%D.prime
##     c(true=D.hat,taylor=taylor)
## })
## plot(pairs['true',],pairs['taylor',])
## abline(0,1)


## [old]
## ## variance of taylor apprx
## ## set.seed(2)
## require(mvtnorm)
## auc.continuous <- function(x,y)mean(outer(x,y,'<'))
## p <- 3
## n <-1e2
## ns <- round(seq(10,1e3,len=10))
## by.n <- sapply(ns, function(n) {
##     beta <- runif(p)
##     gamma <- runif(p)
##     gamma <- beta;# gamma[p] <- 0
##     beta.gamma <- c(beta,gamma)
##     mu.x <- rep(1,p)
##     mu.y <- rep(0,p)
##     Sigma.x <- matrix(rnorm(p^2),nrow=p)
##     Sigma.x <- Sigma.x%*%t(Sigma.x)
##     Sigma.y <- matrix(rnorm(p^2),nrow=p)
##     Sigma.y <- Sigma.y%*%t(Sigma.y)
##     mu <- mu.x-mu.y
##     Sigma <- Sigma.x+Sigma.y
##     ## quad <- as.numeric(beta%*%Sigma%*%beta)
##     ## D.prime.beta <- (Sigma%*%beta %*% (beta%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-beta%*%mu/sqrt(quad))
##     ## D.prime.gamma <- (Sigma%*%beta %*% (beta%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-beta%*%mu/sqrt(quad))
##     D.primes <- lapply(list(beta=beta,gamma=gamma), function(coefs) {
##         quad <- as.numeric(coefs%*%Sigma%*%coefs)
##         (Sigma%*%coefs %*% (coefs%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-coefs%*%mu/sqrt(quad))
##     })
##     D.prime <- c(D.primes$beta,D.primes$gamma)
##     pairs <- replicate(3e2, {
##         beta.hat <- beta+rnorm(p)/sqrt(n)
##         gamma.hat <- gamma+rnorm(p)/sqrt(n)
##         beta.gamma.hat <- c(beta.hat,gamma.hat)
##         ## beta.gamma.hat <- c(beta+rnorm(p)/sqrt(n),gamma+rnorm(p)/sqrt(n))
##         ## beta.gamma.hat <- beta.gamma+rnorm(2*p)/sqrt(n)
##         ## delta <- (beta.gamma.hat-beta.gamma)
##         x <- rmvnorm(n,mu.x,Sigma.x)
##         y <- rmvnorm(n,mu.y,Sigma.y)
##         ## obs <- mean((x-y)%*%beta < 0)
##         ## fla <- pnorm(-(beta%*%mu)/sqrt(t(beta)%*%Sigma%*%beta))
##         ## auc.beta.hat <- auc.continuous(x%*%beta.hat,y%*%beta.hat)
##         aucs <- sapply(list(beta.hat=beta.hat,gamma.hat=gamma.hat), function(coefs)auc.continuous(x%*%coefs,y%*%coefs))
##         ## D <- aucs['beta']-aucs['gamma']
##         D.hat <- unname(aucs['beta.hat']-aucs['gamma.hat']    )
##         taylor <- D.hat + (beta.gamma.hat-beta.gamma)%*%D.prime
##         ## taylor <- D.hat + delta%*%D.prime
##         ## taylor <- D.hat + (rnorm(2*p)/sqrt(n))%*%D.prime
##         c(true=D.hat,taylor=taylor)
##     })
##     c(observed=var(pairs['true',])*n,taylor=D.prime%*%D.prime)
## })
## plot(ns,by.n['observed',]-by.n['taylor',])
## plot(ns,by.n['observed',],ylim=range(by.n),type='l')
## lines(ns,by.n['taylor',],col=2)





## dd
## ## plot(pairs['true',],pairs['taylor',])
## ## abline(0,1)



## n*(X.bar+Y.bar)
n <- 1e2
obs <- replicate(5e3, {
    x <- rnorm(n)
    y <- rnorm(n)
    n^(1)*(mean(x)+mean(y))^2
})
plot(ecdf(obs))
curve(pchisq(x/2,1),add=TRUE,col=2)


## 2nd order taylor

## n*(X.bar+Y.bar) taylor
## obs <- replicate(1e3, {
##     ## x <- rnorm(n)
##     x.bar <- rnorm(1,sd=1/sqrt(n))
##     ## y <- rnorm(n)
##     ## n^(1)*(mean(x)+mean(y))^2
##     n*(x.bar^2 + rnorm(1,sd=1/sqrt(n))*2*x.bar)
## })
B <- 1e5
n <- 3e6
x.bar <- rnorm(B,sd=1/sqrt(n))
y.bar <- rnorm(B,sd=1/sqrt(n))
taylor <-    n*(x.bar^2 + y.bar*2*x.bar + x.bar^2)
plot(ecdf(taylor))
curve(pchisq(x/2,1),add=TRUE,col=2)



## degeneracy
require(mvtnorm)
auc.continuous <- function(x,y)mean(outer(x,y,'<'))
p <- 3
## n <- 1e2
beta <- runif(p)
gamma <- runif(p)
gamma <- beta#; gamma[p] <- 0
beta.gamma <- c(beta,gamma)
mu.x <- rep(1,p)
mu.y <- rep(1,p)
Sigma.x <- matrix(rnorm(p^2),nrow=p)
Sigma.x <- Sigma.x%*%t(Sigma.x)
Sigma.y <- matrix(rnorm(p^2),nrow=p)
Sigma.y <- Sigma.y%*%t(Sigma.y)
mu <- mu.x-mu.y
Sigma <- Sigma.x+Sigma.y
## D.prime <- do.call(rbind,lapply(list(beta=beta,gamma=gamma), function(coefs) {
##     quad <- as.numeric(coefs%*%Sigma%*%coefs)
##     (Sigma%*%coefs %*% (coefs%*%mu/quad) - mu)%*%(1/sqrt(quad)) %*% dnorm(-coefs%*%mu/sqrt(quad))
## }))

n <- 1e2
D.hats <- replicate(1e3, {        
    beta.hat <- beta+runif(p)/sqrt(n)
    gamma.hat <- gamma+runif(p)/sqrt(n)
    ## beta.gamma.hat <- c(beta.hat,gamma.hat)
    x <- rmvnorm(n,mu.x,Sigma.x)
    y <- rmvnorm(n,mu.y,Sigma.y)
    ## obs <- mean((x-y)%*%beta < 0)
    ## fla <- pnorm(-(beta%*%mu)/sqrt(t(beta)%*%Sigma%*%beta))
    ## auc.beta.hat <- auc.continuous(x%*%beta.hat,y%*%beta.hat)
    aucs <- sapply(list(beta.hat=beta.hat,gamma.hat=gamma.hat), function(coefs)auc.continuous(x%*%coefs,y%*%coefs))
    ## D <- aucs['beta']-aucs['gamma']
    D.hat <- unname(aucs['beta.hat']-aucs['gamma.hat']    )
    ## taylor <- D.hat + (beta.gamma.hat-beta.gamma)%*%D.prime
    ## c(true=D.hat,taylor=taylor)
})
op <- par(mfrow=c(1,2))
hist(n^(1/2)*D.hats)
hist(n^(1)*D.hats)
par(op)

plot(ecdf(n^(1/2)*D.hats))

ns <- round(seq(10,2e3,len=30))
by.n <- sapply(ns, function(n) {
    cat('.')
    D.hats <- replicate(1e2, {        
        beta.hat <- beta+runif(p)/sqrt(n)
        gamma.hat <- gamma+runif(p)/sqrt(n)
        ## beta.gamma.hat <- c(beta.hat,gamma.hat)
        x <- rmvnorm(n,mu.x,Sigma.x)
        y <- rmvnorm(n,mu.y,Sigma.y)
        ## obs <- mean((x-y)%*%beta < 0)
        ## fla <- pnorm(-(beta%*%mu)/sqrt(t(beta)%*%Sigma%*%beta))
        ## auc.beta.hat <- auc.continuous(x%*%beta.hat,y%*%beta.hat)
        aucs <- sapply(list(beta.hat=beta.hat,gamma.hat=gamma.hat), function(coefs)auc.continuous(x%*%coefs,y%*%coefs))
        ## D <- aucs['beta']-aucs['gamma']
        D.hat <- unname(aucs['beta.hat']-aucs['gamma.hat']    )
        ## taylor <- D.hat + (beta.gamma.hat-beta.gamma)%*%D.prime
        ## c(true=D.hat,taylor=taylor)
    })
    mean(D.hats)
})

op <- par(mfrow=c(1,2))
plot(ns,ns^(1/2)*by.n);abline(h=0)
plot(ns,ns^(1)*by.n);abline(h=0)
par(op)




## 5. "standard situation" asy variance formula from efron '75
## 5a. variance of beta_0. check using a random contrast.
source('misc.R')
## set.seed(2)
p <- 3
n <- 3e2
Delta <- runif(1) 
lambda <- rlogis(1) 
list2env(lda.sampler.init.std(Delta=Delta,lambda,p),globalenv())
asy.var <- diag(rep(1+Delta^2*pi.1*pi.0,p+1))
asy.var[1:2,1:2] <- matrix(c(1+Delta^2/4,-Delta/2*(pi.0-pi.1),-Delta/2*(pi.0-pi.1),1+2*Delta^2*pi.1*pi.0),byrow=TRUE,ncol=2)
asy.var <- asy.var/pi.0/pi.1
a <-runif(p+1)
vars <- replicate(1e2, {
    beta.hats <- replicate(1e1, {
        beta <- beta 
        ## gamma <- beta; gamma[p+1] <- 0
        ## pairs <- replicate(1e2, {
        ## start <- Sys.time()
        list2env(sample(n),globalenv())
        ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
        beta.hat <- lda.coefs(x.0,x.1)
    })
    ## n*var(hats)
    ## mean(hats)
    n*var(as.numeric(a%*%beta.hats))
})
hist(vars)
## with(beta.sampler, abline(v=(1+Delta^2/4)/(pi.0*pi.1),col=2))
abline(v=t(a)%*%asy.var%*%a,col=2)
abline(v=mean(vars),col=3)
(1+Delta^2/4)/(pi.0*pi.1)


dd

## ## variance of lambda_hat (sanity check)
## source('misc.R')
## ## set.seed(2)
## p <- 3
## n <- 1e2
## Delta <- runif(1) 
## lambda <- rlogis(1) 
## lda.sampler <- lda.sampler.init.std(Delta=Delta,lambda,p)
## beta <- lda.sampler$beta
## mean.vars <- replicate(1e2, {
##     hats <- replicate(1e2, {
##         ## gamma <- beta; gamma[p+1] <- 0
##         ## pairs <- replicate(1e2, {
##         ## start <- Sys.time()
##         list2env(lda.sampler$sample(n),globalenv())
##         ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
##         ## beta.hat <- lda.coefs(x.0,x.1)
##         x.0 <- x.0[,-1]; x.1 <- x.1[,-1]
##         n.0 <- nrow(x.0); n.1 <- nrow(x.1); n <- n.0+n.1
##         pi.0.hat <- n.0/n; pi.1.hat <- n.1/n
##         mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
##         Sigma.hat <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n
##         Sigma.hat.inv <- solve(Sigma.hat)
##         beta.hat <- c(log(pi.1.hat/pi.0.hat) - 1/2*(t(mu.1.hat)%*%Sigma.hat.inv%*%mu.1.hat - t(mu.0.hat)%*%Sigma.hat.inv%*%mu.0.hat), t(mu.1.hat-mu.0.hat)%*%Sigma.hat.inv)
##         log(pi.1.hat/pi.0.hat)
##     })
##     c(mean=mean(sqrt(n)*hats),var=var(sqrt(n)*hats))
## })
## ## mean(hats)
## ## })
## hist(mean.vars['mean',])
## abline(v=lambda,col=2)
## abline(v=mean(mean.vars['mean',]),col=3)
## hist(mean.vars['var',])
## with(lda.sampler, abline(v=1/pi.0/pi.1,col=2))
## abline(v=mean(mean.vars['var',]),col=3)
## with(lda.sampler, 1/pi.0/pi.1,col=2)



m <- matrix(c(1+Delta^2/4,-Delta/2*(pi.0-pi.1),-Delta/2*(pi.0-pi.1),1+2*Delta^2*pi.1*pi.0),byrow=TRUE,ncol=2)
t <- sqrt(sum(diag(m))+2*sqrt(det(m)))
m.root <- 1/t*(m+sqrt(det(m))*diag(2))


m.root <- matrix(c(2+Delta^2/2*(1/2+(Delta^2+6)*pi.0*pi.1),Delta/2*(pi.1-pi.0),Delta/2*(pi.1-pi.0),2+Delta^2*pi.0*pi.1*(2+1/2*(6+Delta^2))),byrow=TRUE,nrow=2) / sqrt(2+Delta^2*(1/4+2*pi.0*pi.1))



## 6. transformation

## for logistic regression as in lda, beta.hat -> A*beta.hat under x
## -> A*x when A is symmetric. Als0 when doing logistic regr with an
## intercept, beta.hat -> A*beta.hat for the non-intercept betas under
## x -> A*x+a when A is symmetric and a is constant.
p <- 3
n <- 20
A <- matrix(runif(p^2),nrow=p)
A <- A%*%t(A)
a <- runif(p)
x <- matrix(rnorm(n*p),ncol=p)
x.transformed <- t(A%*%t(x)+a)
y <- sample(0:1,n,replace=TRUE)
coef(glm(y ~ x, family=binomial))
A%*%coef(glm(y ~ x.transformed, family=binomial(logit)))[-1]
## A%*%coef(glm(y ~ I(x%*%A)-1, family=binomial(probit)))

## dbg. didnt notice A needs to be symmetric.
beta
## beta.bar
## A%*%beta.bar-beta
obj <- function(beta,x,y)-sum(log( exp(x%*%beta*y)/(1+exp(x%*%beta)) ))
optim(par=rep(0,p),obj,x=x,y=y)
## optim(par=rep(0,p),obj,x=x,y=y)$par-beta
## obj(beta,x,y)
x.bar <- t(A%*%t(x))
beta.bar <- coef(glm(y ~ x.bar-1, family=binomial))
beta.bar
optim(par=rep(0,p),obj,x=x.bar,y=y)
optim(par=rep(0,p),obj,x=x.bar,y=y)$par-beta.bar
## obj <- function(beta,x,y)-sum(log( exp(x%*%(solve(A)%*%beta)*y)/(1+exp(x%*%(solve(A)%*%beta))) ))
## optim(par=rep(0,p),obj,x=x.bar,y=y)$par
## solve(A)%*%optim(par=rep(0,p),obj,x=x.bar,y=y)$par
optim(par=rep(0,p),obj,x=x,y=y)$val
optim(par=rep(0,p),fn=function(beta,x,y)obj(A%*%beta,x,y),x=x.bar,y=y)$val

## x -> Ax leaves full and reduced model auc estimates the same when A is block diagonal
require(Matrix)
source('misc.R')
n <- 5e1
p.full <- 8
p.reduced <- 5
A <- bdiag(matrix(runif(p.reduced^2),nrow=p.reduced),matrix(runif((p.full-p.reduced)^2),nrow=p.full-p.reduced))
A <- as.matrix(A)
a <- runif(p.full)
## A <- matrix(runif(p.full^2),ncol=p.full)
## A <- A+t(A)
y <- rbinom(n,1,.7)
x <- matrix(runif(n*p.full),ncol=p.full)
sapply(list(x=x,transformed=t(A%*%t(x)+a)), function(x0) {
    sapply(list(full=x0,reduced=x0[,1:p.reduced]), function(x) {
        beta.hat <- coef(glm(y~x,family=binomial))
        scores.0 <- cbind(1,x[y==0,])%*%beta.hat
        scores.1 <- cbind(1,x[y==1,])%*%beta.hat
        auc.score.hat <- auc.hat(scores.0,scores.1)
        })
})

dd



## 6a
source('misc.R')
p <- 3
n <- 1e3
beta <- runif(p)
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
## auc.scores.deriv2.lda(beta,Sigma)
A <- matrix(runif(p^2),nrow=p)
beta.bar <- as.numeric(solve(t(A))%*%beta)
Sigma.bar <- A%*%Sigma%*%t(A)
A%*%auc.scores.deriv2.lda(beta,Sigma)%*%t(A)
auc.scores.deriv2.lda(beta.bar,Sigma.bar)


x.0 <- matrix(rnorm(n*p),ncol=p)
x.1 <- matrix(rnorm(n*p),ncol=p)
solve(t(A))%*%(lda.coefs(cbind(1,x.0),cbind(1,x.1))[-1])
## lda.coefs(cbind(1,x.0),cbind(1,x.1))
x.0.bar <- x.0%*%t(A)
x.1.bar <- x.1%*%t(A)
lda.coefs(cbind(1,x.0.bar),cbind(1,x.1.bar))


## 6b constant term invariance. P(x.0.bar%*%beta.bar <
## x.1.bar%*%beta.bar) = P(x.0%*%beta < x.1%*%beta) since the A
## factors in the .bar terms cancel out.
source('misc.R')
p <- 3
n <- 1e3
beta <- runif(p+1)
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
## auc.scores.deriv2.lda(beta,Sigma)
A <- matrix(runif(p^2),nrow=p)
list2env(lda.sampler.init(beta,Sigma),globalenv())
mu.diff <- c(0,mu.0-mu.1)
Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
auc.scores(beta,mu.diff,Sigma.diff)
mu.diff.bar <- c(0,A%*%(mu.0-mu.1))
Sigma.diff.bar <- Sigma.diff <- 2*cbind(0,rbind(0,A%*%Sigma%*%t(A)))
beta.bar <- c(beta[1],t(solve(A))%*%beta[-1])
auc.scores(beta.bar,mu.diff.bar,Sigma.diff.bar)

source('misc.R')
require(mvtnorm)
p <- 3
n <- 1e3
beta <- runif(p+1)
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
mu.0 <- runif(p); mu.1 <- runif(p)
## mean.mu <- (mu.0+mu.1)/2
## mu.0 <- mu.0-mean.mu; mu.1 <- mu.1-mean.mu
## lda.sampler <- lda.sampler.init(beta,Sigma)
aucs <- replicate(1e2, {
    x.0 <- x.0.bar <- cbind(1,rmvnorm(n,mean=mu.0,sigma=Sigma))
    x.1 <- x.1.bar <- cbind(1,rmvnorm(n,mean=mu.1,sigma=Sigma))
    ## x.0.bar[,2:(p+1)] <- x.0.bar[,2:(p+1)]%*%t(A); x.1.bar[,2:(p+1)] <- x.1.bar[,2:(p+1)]%*%t(A)    
    auc.hat.continuous(x.0%*%beta,x.1%*%beta)
    ## auc.hat.continuous(x.0.bar%*%beta.bar,x.1.bar%*%beta.bar)
})
hist(aucs)
A <- matrix(runif(p^2),nrow=p)
beta.bar <- c(beta[1],t(solve(A))%*%beta[-1])
mu.diff.bar <- c(0,A%*%(mu.1-mu.0))
Sigma.diff.bar <- 2*cbind(0,rbind(0,A%*%Sigma%*%t(A)))
abline(v=auc.scores(beta.bar,mu.diff.bar,Sigma.diff.bar),col=2)
abline(v=mean(aucs),col=3)
auc.scores(beta.bar,mu.diff=mu.diff.bar,Sigma.diff.bar)


## Efron transformation
source('misc.R')
require(mvtnorm)
p <- 3
n <- 1e3
beta <- runif(p+1)
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
mu.0 <- runif(p); mu.1 <- runif(p)
lda.sampler <- list2env(lda.sampler.init(beta,Sigma),globalenv())
aucs <- replicate(1e2, {
    ## x.0 <- x.0.bar <- cbind(1,rmvnorm(n,mean=mu.0,sigma=Sigma))
    ## x.1 <- x.1.bar <- cbind(1,rmvnorm(n,mean=mu.1,sigma=Sigma))
    ## x.0.bar[,2:(p+1)] <- x.0.bar[,2:(p+1)]%*%t(A); x.1.bar[,2:(p+1)] <- x.1.bar[,2:(p+1)]%*%t(A)
    list2env(lda.sampler$sample(n),globalenv())
    auc.hat.continuous(x.0%*%beta,x.1%*%beta)
    ## auc.hat.continuous(x.0.bar%*%beta.bar,x.1.bar%*%beta.bar)
})
hist(aucs)
## A <- matrix(runif(p^2),nrow=p)
## beta.bar <- c(beta[1],t(solve(A))%*%beta[-1])
## mu.diff.bar <- c(0,A%*%(mu.1-mu.0))
## Sigma.diff.bar <- 2*cbind(0,rbind(0,A%*%Sigma%*%t(A)))
Delta <- sqrt(t(mu.1-mu.0)%*%solve(Sigma)%*%(mu.1-mu.0))
beta.bar <- mu.diff.bar <- c(beta[1],Delta,rep(0,p-1))
mu.diff.bar[1] <- 0
Sigma.diff.bar <- 2*cbind(0,rbind(0,diag(p)))
abline(v=auc.scores(beta.bar,mu.diff.bar,Sigma.diff.bar),col=2)
abline(v=mean(aucs),col=3)
auc.scores(beta.bar,mu.diff=mu.diff.bar,Sigma.diff.bar)







## 7 hajek projection

## 7a checking formula
B <- 1e4
mu.0 <- runif(1); mu.1 <- runif(1)
sigma <- runif(1)
mean(pnorm(rnorm(B)*sigma+mu.0,mean=mu.1,sd=sigma))
pnorm((mu.0-mu.1)/sigma/sqrt(2))

source('misc.R')
require(mvtnorm)
p <- 3
n <- 1e1
beta <- runif(p)
mu.0 <- runif(p); mu.1 <- runif(p)
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
auc.hats <- replicate(1e3,{
    x.0 <- t(rmvnorm(n,mu.0,Sigma))
    x.1 <- t(rmvnorm(n,mu.1,Sigma))
    auc.hat.continuous(t(beta)%*%x.0,t(beta)%*%x.1)
})
hist(auc.hats)
quad <- as.numeric(sqrt(beta%*%Sigma%*%beta))
theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
abline(v=theta,col=2)
abline(v=mean(auc.hats),col=3)


## 7b projection versus observed
source('misc.R')
require(mvtnorm)
p <- 3
n <- 1e1
pairs <- replicate(1e2, {
    beta <- runif(p)
    mu.0 <- runif(p); mu.1 <- runif(p)
    Sigma <- matrix(runif(p^2),nrow=p)
    Sigma <- Sigma%*%t(Sigma)
    theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
    x.0 <- t(rmvnorm(n,mu.0,Sigma))
    x.1 <- t(rmvnorm(n,mu.1,Sigma))
    theta.hat <- auc.hat.continuous(t(beta)%*%x.0,t(beta)%*%x.1)
    observed <- theta.hat - theta
    quad <- as.numeric(beta%*%Sigma%*%beta)
    theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
    ## hajek <- -mean(pnorm(t(beta)%*%(x.0-mu.1) / sqrt(quad)) - as.numeric(pnorm(t(beta)%*%(mu.0-mu.1)/sqrt(2*quad)))) + mean(pnorm(t(beta)%*%(x.1-mu.0) / sqrt(quad)) - as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad))))
    hajek <- -mean(pnorm(t(beta)%*%(x.0-mu.1) / sqrt(quad)) - (1-theta)) + mean(pnorm(t(beta)%*%(x.1-mu.0) / sqrt(quad)) - theta)
    c(observed=observed,hajek=hajek)
})
plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)

## 7c check (U-proj.U)^2 should be O(1/n^2) 
source('misc.R')
require(mvtnorm)
p <- 2
## n <- 1e1
beta <- runif(p)
mu.0 <- runif(p); mu.1 <- runif(p)
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
quad <- as.numeric(beta%*%Sigma%*%beta)
auc.beta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
ns <- round(seq(1e2,5e2,len=20))
by.n <- sapply(ns, function(n) {
    cat('.')
    pairs <- replicate(1e2, {
        ## quad <- as.numeric(sqrt(beta%*%Sigma%*%beta))
        ## auc.beta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
        x.0 <- t(rmvnorm(n,mu.0,Sigma))
        x.1 <- t(rmvnorm(n,mu.1,Sigma))
        auc.hat.beta <- auc.hat.continuous(t(beta)%*%x.0,t(beta)%*%x.1)
        observed <- auc.hat.beta - auc.beta
        ## theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
        ## hajek <- -mean(pnorm(t(beta)%*%(x.0-mu.1) / sqrt(quad)) - as.numeric(pnorm(t(beta)%*%(mu.0-mu.1)/sqrt(2*quad)))) + mean(pnorm(t(beta)%*%(x.1-mu.0) / sqrt(quad)) - as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad))))
        hajek <- -mean(pnorm(t(beta)%*%(x.0-mu.1) / sqrt(quad)) - (1-auc.beta)) + mean(pnorm(t(beta)%*%(x.1-mu.0) / sqrt(quad)) - auc.beta)
        c(observed=observed,hajek=hajek)
    })
    ## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
    pairs['observed',]-pairs['hajek',]
})
## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
matplot(ns,ns^(2)*t(by.n)^2,pch='.',col=1,cex=3)
abline(h=0)

## ---------moved to misc.R
## ## encapsulate hajek projection calculation
## auc.scores.hajek <- function(x.0,x.1,mu.0,mu.1,beta,Sigma) {
##     x.0 <- t(x.0); x.1 <- t(x.1) 
##     quad <- as.numeric(beta%*%Sigma%*%beta)
##     theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
##     hajek <- -mean(pnorm(t(beta)%*%(x.0-mu.1) / sqrt(quad)) - (1-theta)) + mean(pnorm(t(beta)%*%(x.1-mu.0) / sqrt(quad)) - theta)
## }
## source('misc.R')
## require(mvtnorm)
## set.seed(1)
## p <- 2
## ## n <- 1e1
## beta <- runif(p)
## mu.0 <- runif(p); mu.1 <- runif(p)
## Sigma <- matrix(runif(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
## quad <- as.numeric(beta%*%Sigma%*%beta)
## auc.beta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
## ns <- round(seq(1e2,5e2,len=20))
## by.n <- sapply(ns, function(n) {
##     cat('.')
##     pairs <- replicate(1e2, {
##         ## quad <- as.numeric(sqrt(beta%*%Sigma%*%beta))
##         ## auc.beta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
##         x.0 <- rmvnorm(n,mu.0,Sigma)
##         x.1 <- rmvnorm(n,mu.1,Sigma)
##         auc.hat.beta <- auc.hat.continuous(x.0%*%beta,x.1%*%beta)
##         observed <- auc.hat.beta - auc.beta
##         ## theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
##         ## hajek <- -mean(pnorm(t(beta)%*%(x.0-mu.1) / sqrt(quad)) - as.numeric(pnorm(t(beta)%*%(mu.0-mu.1)/sqrt(2*quad)))) + mean(pnorm(t(beta)%*%(x.1-mu.0) / sqrt(quad)) - as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad))))
##         ## hajek <- -mean(pnorm(t(beta)%*%(x.0-mu.1) / sqrt(quad)) - (1-auc.beta)) + mean(pnorm(t(beta)%*%(x.1-mu.0) / sqrt(quad)) - auc.beta)
##         hajek <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta,Sigma)
##         c(observed=observed,hajek=hajek)
##     })
##     ## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
##     pairs['observed',]-pairs['hajek',]
## })
## ## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
## matplot(ns,ns^(2)*t(by.n)^2,pch='.',col=1,cex=3)
## abline(h=0)

## 7d check if (U-proj.U)^2 is O(1/n^2), now with perturbations to
## beta.
source('misc.R')
require(mvtnorm)
## set.seed(1)
p <- 2
## n <- 1e1
beta <- runif(p)
mu.0 <- runif(p); mu.1 <- runif(p)
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
## quad <- as.numeric(beta%*%Sigma%*%beta)
auc.beta <- auc.scores(beta,mu.1-mu.0,2*Sigma)
ns <- round(seq(1e2,5e2,len=20))
by.n <- sapply(ns, function(n) {
    cat('.')
    pairs <- replicate(1e2, {
        beta.hat <- beta+rnorm(length(beta))/sqrt(n)
        ## auc.beta <- auc.scores(beta,mu.1-mu.0,2*Sigma)
        auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
        x.0 <- rmvnorm(n,mu.0,Sigma)
        x.1 <- rmvnorm(n,mu.1,Sigma)
        auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
        ## observed <- auc.hat.beta - auc.beta
        observed <- auc.hat.beta.hat - auc.beta.hat
        hajek <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta.hat,Sigma)
        c(observed=observed,hajek=hajek)
    })
    ## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
    pairs['observed',]-pairs['hajek',]
})
## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
matplot(ns,ns^(2)*t(by.n)^2,pch='.',col=1,cex=3)
matplot(ns,ns*t(by.n),pch='.',col=1,cex=3)
abline(h=0)

dd
## OLD
## source('misc.R')
## require(mvtnorm)
## p <- 2
## ## n <- 1e1
## beta <- runif(p)
## mu.0 <- runif(p); mu.1 <- runif(p)
## Sigma <- matrix(runif(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
## quad <- as.numeric(beta%*%Sigma%*%beta)
## auc.beta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
## ns <- round(seq(1e2,5e2,len=20))
## by.n <- sapply(ns, function(n) {
##     cat('.')
##     pairs <- replicate(1e2, {
##         ## quad <- as.numeric(sqrt(beta%*%Sigma%*%beta))
##         ## auc.beta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
##         beta.hat <- beta+rnorm(length(beta))/sqrt(n)
##         x.0 <- t(rmvnorm(n,mu.0,Sigma))
##         x.1 <- t(rmvnorm(n,mu.1,Sigma))
##         auc.hat.beta <- auc.hat.continuous(t(beta.hat)%*%x.0,t(beta.hat)%*%x.1)
##         observed <- auc.hat.beta - auc.beta
##         ## theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
##         ## hajek <- -mean(pnorm(t(beta)%*%(x.0-mu.1) / sqrt(quad)) - as.numeric(pnorm(t(beta)%*%(mu.0-mu.1)/sqrt(2*quad)))) + mean(pnorm(t(beta)%*%(x.1-mu.0) / sqrt(quad)) - as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad))))
##         hajek <- -mean(pnorm(t(beta.hat)%*%(x.0-mu.1) / sqrt(quad)) - (1-auc.beta)) + mean(pnorm(t(beta.hat)%*%(x.1-mu.0) / sqrt(quad)) - auc.beta)
##         c(observed=observed,hajek=hajek)
##     })
##     ## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
##     pairs['observed',]-pairs['hajek',]
## })
## ## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
## matplot(ns,ns^(2)*t(by.n)^2,pch='.',col=1,cex=3)
## abline(h=0)



## 7e check if (U-proj.U)^2 is O(1/n^2), with lda data
source('misc.R')
require(mvtnorm)
## set.seed(1)
p <- 3
## n <- 1e1
## beta <- runif(p)
## mu.0 <- runif(p); mu.1 <- runif(p)
## Sigma <- matrix(runif(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
beta <- runif(p+1)
mu.0 <- c(1,runif(p)); mu.1 <- c(1,runif(p))
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
Sigma <- cbind(0,rbind(0,Sigma))
mu <- mu.1-mu.0
## quad <- as.numeric(beta%*%Sigma%*%beta)
auc.beta <- auc.scores(beta,mu.1-mu.0,2*Sigma)
ns <- round(seq(1e2,5e2,len=20))
by.n <- sapply(ns, function(n) {
    cat('.')
    pairs <- replicate(1e2, {
        x.0 <- rmvnorm(n,mu.0,Sigma)
        x.1 <- rmvnorm(n,mu.1,Sigma)
        beta.hat <- lda.coefs(x.0,x.1)#beta+rnorm(length(beta))/sqrt(n)
        auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
        auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
        observed <- auc.hat.beta.hat - auc.beta.hat
        hajek <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta.hat,Sigma)
        c(observed=observed,hajek=hajek)
    })
    ## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
    pairs['observed',]-pairs['hajek',]
})
## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
matplot(ns,ns^(2)*t(by.n)^2,pch='.',col=1,cex=3)
abline(h=0)

## #1 imbalance in n.0 and n.1 has a noticeable effect on the size of the
## error. but still same order.
source('misc.R')
require(mvtnorm)
set.seed(1)
p <- 3
## n <- 1e1
## beta <- runif(p)
## mu.0 <- runif(p); mu.1 <- runif(p)
## Sigma <- matrix(runif(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
beta <- runif(p+1)
## mu.0 <- c(1,runif(p)); mu.1 <- c(1,runif(p))
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
Sigma <- cbind(0,rbind(0,Sigma))
lda.sampler <- lda.sampler.init(beta,Sigma[-1,-1])
mu.1 <- c(1,lda.sampler$mu.1)
mu.0 <- c(1,lda.sampler$mu.0)
mu <- mu.1-mu.0
## quad <- as.numeric(beta%*%Sigma%*%beta)
auc.beta <- auc.scores(beta,mu.1-mu.0,2*Sigma)
ns <- round(seq(1e2,5e2,len=20))
by.n <- sapply(ns, function(n) {
    cat('.')
    pairs <- replicate(1e2, {
        n.0 <- rbinom(1,n,.1)#lda.sampler$pi.0) #1
        n.1 <- n-n.0
        x.0 <- rmvnorm(n.0,mu.0,Sigma)
        x.1 <- rmvnorm(n.1,mu.1,Sigma)
        ## x.0 <- rmvnorm(n,mu.0,Sigma)
        ## x.1 <- rmvnorm(n,mu.1,Sigma)
        ## rm(x.0);rm(x.1)
        ## samp <- lda.sampler$sample(n)
        ## x.0 <- samp$x.0; x.1 <- samp$x.1
        beta.hat <- lda.coefs(x.0,x.1)#beta+rnorm(length(beta))/sqrt(n)
        auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
        auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
        observed <- auc.hat.beta.hat - auc.beta.hat
        hajek <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta.hat,Sigma)
        c(observed=observed,hajek=hajek)
    })
    ## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
    pairs['observed',]-pairs['hajek',]
})
## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
matplot(ns,ns^(2)*t(by.n)^2,pch='.',col=1,cex=3)
abline(h=0)


dd


## sampling at beta=beta_LDA
source('misc.R')
require(mvtnorm)
## set.seed(1)
p <- 3
## n <- 1e1
## beta <- runif(p)
## mu.0 <- runif(p); mu.1 <- runif(p)
## Sigma <- matrix(runif(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
beta <- runif(p+1)
## mu.0 <- c(1,runif(p)); mu.1 <- c(1,runif(p))
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
Sigma <- cbind(0,rbind(0,Sigma))
lda.sampler <- lda.sampler.init(beta,Sigma[-1,-1])
mu.1 <- c(1,lda.sampler$mu.1)
mu.0 <- c(1,lda.sampler$mu.0)
mu <- mu.1-mu.0
## quad <- as.numeric(beta%*%Sigma%*%beta)
auc.beta <- auc.scores(beta,mu.1-mu.0,2*Sigma)
ns <- round(seq(1e2,5e2,len=20))
by.n <- sapply(ns, function(n) {
    cat('.')
    pairs <- replicate(1e2, {
        x.0 <- rmvnorm(n,mu.0,Sigma)
        x.1 <- rmvnorm(n,mu.1,Sigma)
        beta.hat <- lda.coefs(x.0,x.1)#beta+rnorm(length(beta))/sqrt(n)
        auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
        auc.hat.beta.hat <- auc.hat.continuous(x.0%*%beta.hat,x.1%*%beta.hat)
        observed <- auc.hat.beta.hat - auc.beta.hat
        hajek <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta.hat,Sigma)
        c(observed=observed,hajek=hajek)
    })
    ## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
    pairs['observed',]-pairs['hajek',]
})
## plot(pairs['observed',],pairs['hajek',]); abline(0,1,col=2)
matplot(ns,ns^(2)*t(by.n)^2,pch='.',col=1,cex=3)
abline(h=0)





n <- 1e2
y <- rbinom(n,1,1/2)
n.0 <- round(n*.9)
z <- c(y[1:n.0],rbinom(n-n.0,1,1/2))
x <- rnorm(n)/20
glm(y ~ z+x,family=binomial)


## trying a mean estimator otherwise similar to the difference of
## u-statistics. seems to be O(n^{1/2}). so 1/n rate comes from the
## u-statistic.  n <- 10
p <- 2
beta <- runif(p)
ns <- round(seq(10,2e4,len=40))
by.n <- sapply(ns, function(n) {
    diffs <- replicate(1e2, {
        x <- matrix(rnorm(n*p),ncol=p)
        epsilon <- runif(p)/sqrt(n)
        eta <- runif(p)/sqrt(n)
        beta.hat <- beta+epsilon
        gamma.hat <- beta+eta
        mean(x%*%beta.hat < 0) - mean(x%*%gamma.hat < 0)
    })
})
matplot(ns,ns*t(by.n),pch=1,col=1)
plot(ns,ns^.5*apply(by.n,2,sd),ylim=c(0,.5));abline(h=0)


## u-statistic version above is O(1/n)
p <- 2
beta <- runif(p)
ns <- round(seq(10,1e3,len=40))
by.n <- sapply(ns, function(n) {
    diffs <- replicate(1e2, {
        x.0 <- matrix(rnorm(n*p),ncol=p)
        x.1 <- matrix(rnorm(n*p),ncol=p)
        epsilon <- runif(p)/sqrt(n)
        eta <- runif(p)/sqrt(n)
        beta.hat <- beta+epsilon
        gamma.hat <- beta+eta
        ## mean(x%*%beta.hat < 0) - mean(x%*%gamma.hat < 0)
        mean(outer(x.0%*%beta.hat,x.1%*%beta.hat,'<')) - mean(outer(x.0%*%gamma.hat,x.1%*%gamma.hat,'<'))
    })
})
matplot(ns,ns*t(by.n),pch=1,col=1)
plot(ns,ns^1*apply(by.n,2,sd));abline(h=0)


## 8 von mises expansion

## 8a taylor expansion in distributions of x0,x1
require(mvtnorm)
source('misc.R')
taylor.term.0 <- function(mu.0,mu.1,Sigma,beta,gamma) pnorm(beta%*%(mu.1-mu.0)/sqrt(2*t(beta)%*%Sigma%*%beta)) - pnorm(gamma%*%(mu.1-mu.0)/sqrt(2*t(gamma)%*%Sigma%*%gamma))
p <- 3
n <- 1e1
pairs <- replicate(1e2, {
    mu.0 <- rep(0,p)
    mu.1 <- runif(p)
    Sigma <- matrix(runif(p^2),nrow=p)
    Sigma <- Sigma%*%t(Sigma)
    beta <- runif(p)
    gamma <- runif(p)
    approx <- taylor.term.0(mu.0,mu.1,Sigma,beta,gamma)# - taylor.term.0(mu.0,mu.1,Sigma,gamma)
    observed <- replicate(1e2, {
        x.0 <- t(rmvnorm(n,mu.0,Sigma))
        x.1 <- t(rmvnorm(n,mu.1,Sigma))
        ## observed <- mean(x.0%*%beta<x.1%*%beta) - mean(x.0%*%gamma<x.1%*%gamma)
        observed <- auc.hat.continuous(t(beta)%*%x.0,t(beta)%*%x.1) - auc.hat.continuous(gamma%*%x.0,gamma%*%x.1)
    })
    c(observed=mean(observed),approx=approx)
})
plot(pairs['observed',],pairs['approx',]);abline(0,1)

require(mvtnorm)
taylor.term.0 <- function(mu.0,mu.1,Sigma,beta,gamma) pnorm(beta%*%(mu.1-mu.0)/sqrt(2*t(beta)%*%Sigma%*%beta)) - pnorm(gamma%*%(mu.1-mu.0)/sqrt(2*t(gamma)%*%Sigma%*%gamma))
p <- 3
n <- 1e1
## pairs <- replicate(1e2, {
mu.0 <- rep(0,p)
mu.1 <- runif(p)
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
gamma <- runif(p)
approx <- taylor.term.0(mu.0,mu.1,Sigma,beta,gamma)
ns <- round(seq(10,1e3,len=20))
by.n <- sapply(ns, function(n) {
    observed <- replicate(1e2, {
        x.0 <- t(rmvnorm(n,mu.0,Sigma))
        x.1 <- t(rmvnorm(n,mu.1,Sigma))
        auc.hat.continuous(t(beta)%*%x.0,t(beta)%*%x.1) - auc.hat.continuous(gamma%*%x.0,gamma%*%x.1)
        ## mean(t(beta)%*%x.0<t(beta)%*%x.1) - mean(gamma%*%x.0<gamma%*%x.1)
    })
    ## c(observed=mean(observed),approx=approx)
    observed - c(approx)
})
plot(ns,ns^(1/2)*apply(by.n,2,sd))



taylor.term.0 <- function(mu.0,mu.1,Sigma,beta,gamma) pnorm(beta%*%(mu.1-mu.0)/sqrt(2*t(beta)%*%Sigma%*%beta)) - pnorm(gamma%*%(mu.1-mu.0)/sqrt(2*t(gamma)%*%Sigma%*%gamma))
taylor.term.1 <- function(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma) {
    partials.y <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.0-mu.1)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    partials.x <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.1-mu.0)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    unname(partials.y['gamma'] - partials.y['beta'] + partials.x['beta'] - partials.x['gamma']) - 2*taylor.term.0(mu.0,mu.1,Sigma,beta,gamma)
}
require(mvtnorm)
p <- 3
n <- 1e1
## pairs <- replicate(1e2, {
mu.0 <- rep(0,p)
mu.1 <- runif(p)
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
gamma <- runif(p)
ns <- round(seq(10,2e2,len=20))
const.term <- taylor.term.0(mu.0,mu.1,Sigma,beta,gamma)
by.n <- sapply(ns, function(n) {
    parts <- replicate(1e2, {
        x.0 <- t(rmvnorm(n,mu.0,Sigma))
        x.1 <- t(rmvnorm(n,mu.1,Sigma))
        auc.hat.beta <-  auc.hat.continuous(t(beta)%*%x.0,t(beta)%*%x.1) - auc.hat.continuous(gamma%*%x.0,gamma%*%x.1)
        ## truncated <- observed - const.term
        approx <- taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma)
        c(observed=auc.hat.beta - const.term,approx=approx)
        ## truncated-remainder
        ## partials.y <- sapply(list(beta=beta,gamma=gamma), function(beta)
        ##     mean(pnorm(t(beta)%*%(x.0-mu.1)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
        ## partials.x <- sapply(list(beta=beta,gamma=gamma), function(beta)
        ##     mean(pnorm(t(beta)%*%(x.1-mu.0)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
        ## observed + taylor.term.0(mu.0,mu.1,Sigma,beta,gamma) - (partials.y['gamma'] - partials.y['beta'] + partials.x['beta'] - partials.x['gamma'] )        
    })
    ## c(observed=mean(observed),approx=approx)
    parts['observed',] - parts['approx',]
})
plot(ns,ns^(1)*apply(by.n,2,sd))


taylor.term.0 <- function(mu.0,mu.1,Sigma,beta,gamma) pnorm(beta%*%(mu.1-mu.0)/sqrt(2*t(beta)%*%Sigma%*%beta)) - pnorm(gamma%*%(mu.1-mu.0)/sqrt(2*t(gamma)%*%Sigma%*%gamma))
taylor.term.1 <- function(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma) {
    partials.y <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.0-mu.1)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    partials.x <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.1-mu.0)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    unname(partials.y['gamma'] - partials.y['beta'] + partials.x['beta'] - partials.x['gamma']) - 2*taylor.term.0(mu.0,mu.1,Sigma,beta,gamma)
}
require(mvtnorm)
p <- 3
n <- 1e1
## pairs <- replicate(1e2, {
mu.0 <- rep(0,p)
mu.1 <- runif(p)
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
gamma <- runif(p)
ns <- round(seq(10,2e2,len=20))
const.term <- taylor.term.0(mu.0,mu.1,Sigma,beta,gamma)
by.n <- sapply(ns, function(n) {
    parts <- replicate(1e2, {
        x.0 <- t(rmvnorm(n,mu.0,Sigma))
        x.1 <- t(rmvnorm(n,mu.1,Sigma))
        ## auc.hat.beta <-  auc.hat.continuous(t(beta)%*%x.0,t(beta)%*%x.1) - auc.hat.continuous(gamma%*%x.0,gamma%*%x.1)
        ## truncated <- observed - const.term
        hajek <- taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma)
        ## c(observed=auc.hat.beta - const.term,approx=approx)
        ## truncated-remainder
        ## partials.y <- sapply(list(beta=beta,gamma=gamma), function(beta)
        ##     mean(pnorm(t(beta)%*%(x.0-mu.1)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
        ## partials.x <- sapply(list(beta=beta,gamma=gamma), function(beta)
        ##     mean(pnorm(t(beta)%*%(x.1-mu.0)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
        ## observed + taylor.term.0(mu.0,mu.1,Sigma,beta,gamma) - (partials.y['gamma'] - partials.y['beta'] + partials.x['beta'] - partials.x['gamma'] )        
    })
    ## c(observed=mean(observed),approx=approx)
    ## parts['observed',] - parts['approx',]
})
sds <- apply(by.n,2,sd)
lm0 <- lm(log(sds) ~ log(ns))
lm0
## plot(ns,log(sds))
## lines(ns,fitted(lm0),col=2)
## plot(ns,ns^(1)*apply(by.n,2,sd))


## 8b check difference of quad term evaluated at beta.hat/gamma.hat versus beta/gamma. seems o(1/n), maybe 1/n^1.25.
taylor.term.0 <- function(mu.0,mu.1,Sigma,beta,gamma) pnorm(beta%*%(mu.1-mu.0)/sqrt(2*t(beta)%*%Sigma%*%beta)) - pnorm(gamma%*%(mu.1-mu.0)/sqrt(2*t(gamma)%*%Sigma%*%gamma))
taylor.term.1 <- function(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma) {
    partials.y <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.0-mu.1)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    partials.x <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.1-mu.0)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    unname(partials.y['gamma'] - partials.y['beta'] + partials.x['beta'] - partials.x['gamma']) - 2*taylor.term.0(mu.0,mu.1,Sigma,beta,gamma)
}
start <- Sys.time()
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
p <- 3
n <- 1e1
## pairs <- replicate(1e2, {
mu.0 <- c(1,runif(p))
mu.1 <- c(1,runif(p))
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
Sigma <- rbind(0,cbind(0,Sigma))
beta <- runif(p+1)
gamma <- beta#runif(p+1)
## beta <- c(-1/2*(mu.1[-1]%*%solve(Sigma[-1,-1])%*%mu.1[-1]-mu.2[-1]%*%solve(Sigma[-1,-1])%*%mu.2[-1]),solve(Sigma[-1,-1])%*%(mu.1-mu.0))
ns <- round(seq(10,2e4,len=40))
## const.term <
by.n <- sapply(ns, function(n) {
## by.n <- mclapply(ns, mc.cores=detectCores()-3,FUN=function(n) {
    ## print(max(1,detectCores()-3))
    print(n)
    ## diffs <- mclapply(1:1e2, mc.cores=4,FUN=function(jj) {
    diffs <- replicate(1e2, {
        x.0 <- t(rmvnorm(n,mu.0,Sigma))
        x.1 <- t(rmvnorm(n,mu.1,Sigma))
        beta.hat <- beta+rnorm(p+1)/sqrt(n)
        gamma.hat <- gamma+rnorm(p+1)/sqrt(n)
        auc.hat.coefs <-  auc.hat.continuous(t(beta)%*%x.0,t(beta)%*%x.1) - auc.hat.continuous(gamma%*%x.0,gamma%*%x.1)
        ## truncated <- observed - const.term
        quad.beta <- auc.hat.coefs - taylor.term.0(mu.0,mu.1,Sigma,beta,gamma) - taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma) # quad.beta should be 0 so could speed this up by removing
        auc.hat.coefs.hat <-  auc.hat.continuous(t(beta.hat)%*%x.0,t(beta.hat)%*%x.1) - auc.hat.continuous(gamma.hat%*%x.0,gamma.hat%*%x.1)
        ## truncated <- observed - const.term
        quad.beta.hat <- auc.hat.coefs.hat - taylor.term.0(mu.0,mu.1,Sigma,beta.hat,gamma.hat) - taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta.hat,gamma.hat)
        quad.beta - quad.beta.hat
    })
    simplify2array(diffs)
})
## by.n <- simplify2array(by.n)
plot(ns,ns^(1)*apply(by.n,2,sd))
## save.image('221005.RData')
## save.image('221006.RData')
## save.image('221007.RData')
Sys.time() - start

## load('221006.RData')
## good.idx <- (1:40)%%5==1
## by.n <- by.n[good.idx]
## ns <- ns[good.idx]
## plot(ns,ns^(1)*sapply(by.n,sd))


## rate seems to be 1/n^(1.25)
load('221007.RData')
sds <- apply(by.n,2,sd)
lm0 <- lm(log(sds) ~ log(ns))
plot(ns,log(sds))
lines(ns,fitted(lm0),col=2)



## 8ba with lda beta.hat/gamma.hat. still seems o(1/n).
taylor.term.0 <- function(mu.0,mu.1,Sigma,beta,gamma) pnorm(beta%*%(mu.1-mu.0)/sqrt(2*t(beta)%*%Sigma%*%beta)) - pnorm(gamma%*%(mu.1-mu.0)/sqrt(2*t(gamma)%*%Sigma%*%gamma))
taylor.term.1 <- function(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma) {
    partials.y <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.0-mu.1)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    partials.x <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.1-mu.0)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    unname(partials.y['gamma'] - partials.y['beta'] + partials.x['beta'] - partials.x['gamma']) - 2*taylor.term.0(mu.0,mu.1,Sigma,beta,gamma)
}
start <- Sys.time()
require(mvtnorm)
require(parallel)
source('misc.R')
## set.seed(1)
p <- 3
n <- 1e1
## pairs <- replicate(1e2, {
## mu.0 <- c(1,runif(p))
## mu.1 <- c(1,runif(p))
## Sigma <- matrix(runif(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
## Sigma <- rbind(0,cbind(0,Sigma))
beta <- c(runif(p),0)
gamma <- beta
Sigma.lda <- matrix(rnorm(p^2),nrow=p)
Sigma.lda <- Sigma.lda%*%t(Sigma.lda)
beta.sampler <- lda.sampler.init(beta,Sigma.lda)        
## beta <- c(runif(p),0)
## gamma <- beta#runif(p+1)
## ## beta <- c(-1/2*(mu.1[-1]%*%solve(Sigma[-1,-1])%*%mu.1[-1]-mu.2[-1]%*%solve(Sigma[-1,-1])%*%mu.2[-1]),solve(Sigma[-1,-1])%*%(mu.1-mu.0))
ns <- round(seq(1e2,3e3,len=20))
## const.term <
by.n <- sapply(ns, function(n) {
## by.n <- mclapply(ns, mc.cores=detectCores()-3,FUN=function(n) {
    ## print(max(1,detectCores()-3))
    print(n)
    ## diffs <- mclapply(1:1e2, mc.cores=4,FUN=function(jj) {
    diffs <- replicate(1e2, {
        ## x.0 <- t(rmvnorm(n,mu.0,Sigma))
        ## x.1 <- t(rmvnorm(n,mu.1,Sigma))
        list2env(beta.sampler$sample(n),globalenv())
        ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
        ## gamma.hat <- gamma+rnorm(p+1)/sqrt(n)
        beta.hat <- lda.coefs(x.0,x.1)
        gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
        gamma.hat <- c(gamma.hat,0)
        x.0 <- t(x.0); x.1 <- t(x.1) ## decide row vs col oriented
        auc.hat.coefs <-  auc.hat.continuous(t(beta)%*%x.0,t(beta)%*%x.1) - auc.hat.continuous(gamma%*%x.0,gamma%*%x.1)
        mu.0 <- c(1,beta.sampler$mu.0); mu.1 <- c(1,beta.sampler$mu.1)
        Sigma <- cbind(0,rbind(0,Sigma.lda))
        quad.beta <- auc.hat.coefs - taylor.term.0(mu.0,mu.1,Sigma,beta,gamma) - taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma)
        auc.hat.coefs.hat <-  auc.hat.continuous(t(beta.hat)%*%x.0,t(beta.hat)%*%x.1) - auc.hat.continuous(gamma.hat%*%x.0,gamma.hat%*%x.1)
        quad.beta.hat <- auc.hat.coefs.hat - taylor.term.0(mu.0,mu.1,Sigma,beta.hat,gamma.hat) - taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta.hat,gamma.hat)
        quad.beta - quad.beta.hat
    })
    simplify2array(diffs)
})
plot(ns,ns^(1)*apply(by.n,2,sd))
Sys.time() - start

sds <- apply(by.n,2,sd)
lm0 <- lm(log(sds) ~ log(ns))
coef(lm0)
plot(ns,log(sds))
lines(ns,fitted(lm0),col=2)


## 8bb now looking also at rate of const and linear terms. as
## expected, const and linear terms are O(1/n) and quad term is o(1/n)
taylor.term.0 <- function(mu.0,mu.1,Sigma,beta,gamma) pnorm(beta%*%(mu.1-mu.0)/sqrt(2*t(beta)%*%Sigma%*%beta)) - pnorm(gamma%*%(mu.1-mu.0)/sqrt(2*t(gamma)%*%Sigma%*%gamma))
taylor.term.1 <- function(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma) {
    partials.y <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.0-mu.1)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    partials.x <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.1-mu.0)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    unname(partials.y['gamma'] - partials.y['beta'] + partials.x['beta'] - partials.x['gamma']) - 2*taylor.term.0(mu.0,mu.1,Sigma,beta,gamma)
}
start <- Sys.time()
require(mvtnorm)
require(parallel)
source('misc.R')
## set.seed(1)
p <- 3
n <- 1e1
## pairs <- replicate(1e2, {
## mu.0 <- c(1,runif(p))
## mu.1 <- c(1,runif(p))
## Sigma <- matrix(runif(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
## Sigma <- rbind(0,cbind(0,Sigma))
beta <- c(runif(p),0)
gamma <- beta
Sigma.lda <- matrix(rnorm(p^2),nrow=p)
Sigma.lda <- Sigma.lda%*%t(Sigma.lda)
beta.sampler <- lda.sampler.init(beta,Sigma.lda)        
## beta <- c(runif(p),0)
## gamma <- beta#runif(p+1)
## ## beta <- c(-1/2*(mu.1[-1]%*%solve(Sigma[-1,-1])%*%mu.1[-1]-mu.2[-1]%*%solve(Sigma[-1,-1])%*%mu.2[-1]),solve(Sigma[-1,-1])%*%(mu.1-mu.0))
ns <- round(seq(1e2,5e3,len=20))
## const.term <
by.n <- lapply(ns, function(n) {
## by.n <- mclapply(ns, mc.cores=detectCores()-3,FUN=function(n) {
    ## print(max(1,detectCores()-3))
    print(n)
    ## diffs <- mclapply(1:1e2, mc.cores=4,FUN=function(jj) {
    taylor.terms <- replicate(1e2, {
        ## x.0 <- t(rmvnorm(n,mu.0,Sigma))
        ## x.1 <- t(rmvnorm(n,mu.1,Sigma))
        list2env(beta.sampler$sample(n),globalenv())
        ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
        ## gamma.hat <- gamma+rnorm(p+1)/sqrt(n)
        beta.hat <- lda.coefs(x.0,x.1)
        gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
        gamma.hat <- c(gamma.hat,0)
        x.0 <- t(x.0); x.1 <- t(x.1) ## decide row vs col oriented
        ## auc.hat.coefs <-  auc.hat.continuous(t(beta)%*%x.0,t(beta)%*%x.1) - auc.hat.continuous(gamma%*%x.0,gamma%*%x.1)
        mu.0 <- c(1,beta.sampler$mu.0); mu.1 <- c(1,beta.sampler$mu.1)
        Sigma <- cbind(0,rbind(0,Sigma.lda))
        ## quad.beta <- auc.hat.coefs - taylor.term.0(mu.0,mu.1,Sigma,beta,gamma) - taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma)
        auc.hat.coefs.hat <-  auc.hat.continuous(t(beta.hat)%*%x.0,t(beta.hat)%*%x.1) - auc.hat.continuous(gamma.hat%*%x.0,gamma.hat%*%x.1)
        const.beta.hat <- taylor.term.0(mu.0,mu.1,Sigma,beta.hat,gamma.hat)
        linear.beta.hat <- taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta.hat,gamma.hat)
        quad.beta.hat <- auc.hat.coefs.hat - const.beta.hat - linear.beta.hat
        c(const=const.beta.hat,linear=linear.beta.hat,quad=quad.beta.hat)
    })
    ## simplify2array(taylor.terms)
})
by.n <- simplify2array(by.n)
Sys.time() - start
apply(by.n,1,function(m){
    sds <- apply(m,2,sd)
    lm0 <- lm(log(sds) ~ log(ns))
    coef(lm0)[2]
})

dd

## ## 8ba with lda beta.hat/gamma.hat. still seems o(1/n).
## taylor.term.0 <- function(mu.0,mu.1,Sigma,beta,gamma) pnorm(beta%*%(mu.1-mu.0)/sqrt(2*t(beta)%*%Sigma%*%beta)) - pnorm(gamma%*%(mu.1-mu.0)/sqrt(2*t(gamma)%*%Sigma%*%gamma))
## taylor.term.1 <- function(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma) {
##     partials.y <- sapply(list(beta=beta,gamma=gamma), function(beta)
##         mean(pnorm(t(beta)%*%(x.0-mu.1)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
##     partials.x <- sapply(list(beta=beta,gamma=gamma), function(beta)
##         mean(pnorm(t(beta)%*%(x.1-mu.0)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
##     unname(partials.y['gamma'] - partials.y['beta'] + partials.x['beta'] - partials.x['gamma']) - 2*taylor.term.0(mu.0,mu.1,Sigma,beta,gamma)
## }
## start <- Sys.time()
## require(mvtnorm)
## require(parallel)
## source('misc.R')
## ## set.seed(1)
## p <- 3
## n <- 1e1
## ## pairs <- replicate(1e2, {
## ## mu.0 <- c(1,runif(p))
## ## mu.1 <- c(1,runif(p))
## ## Sigma <- matrix(runif(p^2),nrow=p)
## ## Sigma <- Sigma%*%t(Sigma)
## ## Sigma <- rbind(0,cbind(0,Sigma))
## beta <- c(runif(p),0)
## gamma <- beta
## Sigma.lda <- matrix(rnorm(p^2),nrow=p)
## Sigma.lda <- Sigma.lda%*%t(Sigma.lda)
## beta.sampler <- lda.sampler.init(beta,Sigma.lda)        
## ## beta <- c(runif(p),0)
## ## gamma <- beta#runif(p+1)
## ## ## beta <- c(-1/2*(mu.1[-1]%*%solve(Sigma[-1,-1])%*%mu.1[-1]-mu.2[-1]%*%solve(Sigma[-1,-1])%*%mu.2[-1]),solve(Sigma[-1,-1])%*%(mu.1-mu.0))
## ns <- round(seq(1e2,3e3,len=20))
## ## const.term <
## ## by.n <- sapply(ns, function(n) {
## ## ## by.n <- mclapply(ns, mc.cores=detectCores()-3,FUN=function(n) {
## ##     ## print(max(1,detectCores()-3))
## ##     print(n)
## ##     ## diffs <- mclapply(1:1e2, mc.cores=4,FUN=function(jj) {
## ##     diffs <- replicate(1e2, {
## ##         ## x.0 <- t(rmvnorm(n,mu.0,Sigma))
## ##         ## x.1 <- t(rmvnorm(n,mu.1,Sigma))
## list2env(beta.sampler$sample(n),globalenv())
## ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
## ## gamma.hat <- gamma+rnorm(p+1)/sqrt(n)
## beta.hat <- lda.coefs(x.0,x.1)
## gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
## gamma.hat <- c(gamma.hat,0)
## x.0 <- t(x.0); x.1 <- t(x.1) ## decide row vs col oriented
## diff.hat.coefs <-  auc.hat.continuous(t(beta)%*%x.0,t(beta)%*%x.1) - auc.hat.continuous(gamma%*%x.0,gamma%*%x.1)
## mu.0 <- c(1,beta.sampler$mu.0); mu.1 <- c(1,beta.sampler$mu.1)
## Sigma <- cbind(0,rbind(0,Sigma.lda))
## quad.beta <- diff.hat.coefs - taylor.term.0(mu.0,mu.1,Sigma,beta,gamma) - taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma)
## diff.hat.coefs.hat <-  auc.hat.continuous(t(beta.hat)%*%x.0,t(beta.hat)%*%x.1) - auc.hat.continuous(gamma.hat%*%x.0,gamma.hat%*%x.1)
## quad.beta.hat <- diff.hat.coefs.hat - taylor.term.0(mu.0,mu.1,Sigma,beta.hat,gamma.hat) - taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta.hat,gamma.hat)
## quad.beta - quad.beta.hat


## })
## simplify2array(diffs)
## })
## plot(ns,ns^(1)*apply(by.n,2,sd))
## Sys.time() - start





## 8c growth rate of const term with hatted coefs, linear term with hatted coefs. both appear O(1/n).
taylor.term.0 <- function(mu.0,mu.1,Sigma,beta,gamma) pnorm(beta%*%(mu.1-mu.0)/sqrt(2*t(beta)%*%Sigma%*%beta)) - pnorm(gamma%*%(mu.1-mu.0)/sqrt(2*t(gamma)%*%Sigma%*%gamma))
taylor.term.1 <- function(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma) {
    partials.y <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.0-mu.1)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    partials.x <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.1-mu.0)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    unname(partials.y['gamma'] - partials.y['beta'] + partials.x['beta'] - partials.x['gamma']) - 2*taylor.term.0(mu.0,mu.1,Sigma,beta,gamma)
}
start <- Sys.time()
require(mvtnorm)
require(parallel)
source('misc.R')
## set.seed(1)
p <- 3
n <- 1e1
## pairs <- replicate(1e2, {
## mu.0 <- c(1,runif(p))
## mu.1 <- c(1,runif(p))
## Sigma <- matrix(runif(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
## Sigma <- rbind(0,cbind(0,Sigma))
beta <- c(runif(p),0)
gamma <- beta
Sigma.lda <- matrix(rnorm(p^2),nrow=p)
Sigma.lda <- Sigma.lda%*%t(Sigma.lda)
beta.sampler <- lda.sampler.init(beta,Sigma.lda)        
## beta <- c(runif(p),0)
## gamma <- beta#runif(p+1)
## ## beta <- c(-1/2*(mu.1[-1]%*%solve(Sigma[-1,-1])%*%mu.1[-1]-mu.2[-1]%*%solve(Sigma[-1,-1])%*%mu.2[-1]),solve(Sigma[-1,-1])%*%(mu.1-mu.0))
ns <- round(seq(1e2,1e3,len=20))
## const.term <
by.n <- sapply(ns, function(n) {
## by.n <- mclapply(ns, mc.cores=detectCores()-3,FUN=function(n) {
    ## print(max(1,detectCores()-3))
    print(n)
    ## diffs <- mclapply(1:1e2, mc.cores=4,FUN=function(jj) {
    terms <- replicate(1e2, {
        ## x.0 <- t(rmvnorm(n,mu.0,Sigma))
        ## x.1 <- t(rmvnorm(n,mu.1,Sigma))
        list2env(beta.sampler$sample(n),globalenv())
        ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
        ## gamma.hat <- gamma+rnorm(p+1)/sqrt(n)
        beta.hat <- lda.coefs(x.0,x.1)
        gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
        gamma.hat <- c(gamma.hat,0)
        x.0 <- t(x.0); x.1 <- t(x.1) ## decide row vs col oriented
        mu.0 <- c(1,beta.sampler$mu.0); mu.1 <- c(1,beta.sampler$mu.1)
        Sigma <- cbind(0,rbind(0,Sigma.lda))
        ## auc.hat.coefs <-  auc.hat.continuous(t(beta)%*%x.0,t(beta)%*%x.1) - auc.hat.continuous(gamma%*%x.0,gamma%*%x.1)
        ## quad.beta <- auc.hat.coefs - taylor.term.0(mu.0,mu.1,Sigma,beta,gamma) - taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma)
        ## auc.hat.coefs.hat <-  auc.hat.continuous(t(beta.hat)%*%x.0,t(beta.hat)%*%x.1) - auc.hat.continuous(gamma.hat%*%x.0,gamma.hat%*%x.1)
        ## quad.beta.hat <- auc.hat.coefs.hat - taylor.term.0(mu.0,mu.1,Sigma,beta.hat,gamma.hat) - taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta.hat,gamma.hat)
        ## quad.beta - quad.beta.hat
        c(const.hat=taylor.term.0(mu.0,mu.1,Sigma,beta.hat,gamma.hat),
        linear.hat=taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta.hat,gamma.hat))
    })    
},simplify=FALSE)
by.n <-  simplify2array(by.n)
## plot(ns,ns^(1)*apply(by.n,2,sd))
Sys.time() - start
op <- par(mfrow=c(1,2))
sds <- apply(by.n['const.hat',,],2,sd)
lm0 <- lm(log(sds) ~ log(ns))
coef(lm0)
plot(ns,log(sds))
lines(ns,fitted(lm0),col=2)
sds <- apply(by.n['linear.hat',,],2,sd)
lm0 <- lm(log(sds) ~ log(ns))
coef(lm0)
plot(ns,log(sds))
lines(ns,fitted(lm0),col=2)
par(op)



## 8d diff between diff.hat.coefs.hat and 1st order taylor
## approx--appears o(1/n) as expected from
## quad.beta.hat-quad.beta=quad.beta.hat being o(1/n)
taylor.term.0 <- function(mu.0,mu.1,Sigma,beta,gamma) pnorm(beta%*%(mu.1-mu.0)/sqrt(2*t(beta)%*%Sigma%*%beta)) - pnorm(gamma%*%(mu.1-mu.0)/sqrt(2*t(gamma)%*%Sigma%*%gamma))
taylor.term.1 <- function(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma) {
    partials.y <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.0-mu.1)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    partials.x <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.1-mu.0)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    unname(partials.y['gamma'] - partials.y['beta'] + partials.x['beta'] - partials.x['gamma']) - 2*taylor.term.0(mu.0,mu.1,Sigma,beta,gamma)
}
start <- Sys.time()
require(mvtnorm)
require(parallel)
source('misc.R')
## set.seed(1)
p <- 3
n <- 1e1
## pairs <- replicate(1e2, {
## mu.0 <- c(1,runif(p))
## mu.1 <- c(1,runif(p))
## Sigma <- matrix(runif(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
## Sigma <- rbind(0,cbind(0,Sigma))
beta <- c(runif(p),0)
gamma <- beta
Sigma.lda <- matrix(rnorm(p^2),nrow=p)
Sigma.lda <- Sigma.lda%*%t(Sigma.lda)
beta.sampler <- lda.sampler.init(beta,Sigma.lda)        
## beta <- c(runif(p),0)
## gamma <- beta#runif(p+1)
## ## beta <- c(-1/2*(mu.1[-1]%*%solve(Sigma[-1,-1])%*%mu.1[-1]-mu.2[-1]%*%solve(Sigma[-1,-1])%*%mu.2[-1]),solve(Sigma[-1,-1])%*%(mu.1-mu.0))
ns <- round(seq(1e2,3e3,len=20))
## const.term <
by.n <- sapply(ns, function(n) {
## by.n <- mclapply(ns, mc.cores=detectCores()-3,FUN=function(n) {
    ## print(max(1,detectCores()-3))
    print(n)
    ## diffs <- mclapply(1:1e2, mc.cores=4,FUN=function(jj) {
    diffs <- replicate(1e2, {
        ## x.0 <- t(rmvnorm(n,mu.0,Sigma))
        ## x.1 <- t(rmvnorm(n,mu.1,Sigma))
        list2env(beta.sampler$sample(n),globalenv())
        ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
        ## gamma.hat <- gamma+rnorm(p+1)/sqrt(n)
        beta.hat <- lda.coefs(x.0,x.1)
        gamma.hat <- lda.coefs(x.0[,-(p+1)],x.1[,-(p+1)])
        gamma.hat <- c(gamma.hat,0)
        x.0 <- t(x.0); x.1 <- t(x.1) ## decide row vs col oriented
        ## auc.hat.coefs <-  auc.hat.continuous(t(beta)%*%x.0,t(beta)%*%x.1) - auc.hat.continuous(gamma%*%x.0,gamma%*%x.1)
        mu.0 <- c(1,beta.sampler$mu.0); mu.1 <- c(1,beta.sampler$mu.1)
        Sigma <- cbind(0,rbind(0,Sigma.lda))
        ## quad.beta <- auc.hat.coefs - taylor.term.0(mu.0,mu.1,Sigma,beta,gamma) - taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma)
        auc.hat.coefs.hat <-  auc.hat.continuous(t(beta.hat)%*%x.0,t(beta.hat)%*%x.1) - auc.hat.continuous(gamma.hat%*%x.0,gamma.hat%*%x.1)
        ## quad.beta.hat <- auc.hat.coefs.hat - taylor.term.0(mu.0,mu.1,Sigma,beta.hat,gamma.hat) - taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta.hat,gamma.hat)
        ## quad.beta - quad.beta.hat
        const.coefs.hat <- taylor.term.0(mu.0,mu.1,Sigma,beta.hat,gamma.hat)
        linear.coefs.hat <- taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta.hat,gamma.hat)
        auc.hat.coefs.hat - const.coefs.hat - linear.coefs.hat
    })
    simplify2array(diffs)
})
plot(ns,ns^(1)*apply(by.n,2,sd))
Sys.time() - start
sds <- apply(by.n,2,sd)
lm0 <- lm(log(sds) ~ log(ns))
coef(lm0)
plot(ns,log(sds))
lines(ns,fitted(lm0),col=2)


## 8e expanding F-linear term in beta ("D.12"). is O(1/n).

## dbg--checking formulas
u <- function(beta,w)pnorm(t(beta)%*%w/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))
u.1 <- function(beta,w) {
    qform <- as.numeric(t(beta)%*%Sigma%*%beta)
    as.numeric(-dnorm(t(beta)%*%w/sqrt(qform))/sqrt(qform))*(as.numeric(t(beta)%*%w/qform)*Sigma%*%beta-w)
}
D.1 <- function(x.0,x.1,mu.0,mu.1,beta,gamma) {
    mean(apply(x.0-mu.1,2,function(w)u(gamma,w) - u(beta,w))) +  mean(apply(x.1-mu.0,2,function(w)u(beta,w) - u(gamma,w))) -2*(u(beta,(mu.1-mu.0)/sqrt(2)) - u(gamma,(mu.1-mu.0)/sqrt(2))) 
}
## 1. check u.1
## delta <- runif(length(beta)) 
## w <- runif(length(beta))
## ts <- seq(0,1,len=20)
## us <- sapply(ts, function(t)u(beta+t*delta,w))
## plot(ts,us,type='l')
## abline(a=us[1],b=t(u.1(beta,w))%*%delta,col=2) 
## 2. check u
taylor.term.1(x.0,x.1,mu.0,mu.1,Sigma,beta.hat,gamma.hat)
mean(apply(x.0-mu.1,2,function(w)u(gamma.hat,w) - u(beta.hat,w))) +  mean(apply(x.1-mu.0,2,function(w)u(beta.hat,w) - u(gamma.hat,w))) -2*(u(beta.hat,(mu.1-mu.0)/sqrt(2)) - u(gamma.hat,(mu.1-mu.0)/sqrt(2))) 
D.1(x.0,x.1,mu.0,mu.1,beta.hat,gamma.hat)
deriv <- rowMeans(rbind(-apply(x.0-mu.1,2,function(w)u.1(beta,w)),apply(x.0-mu.1,2,function(w)u.1(gamma,w)))) + rowMeans(rbind(apply(x.1-mu.0,2,function(w)u.1(beta,w)), -apply(x.1-mu.0,2,function(w)u.1(gamma,w)))) - 2*c(u.1(beta,(mu.1-mu.0)/sqrt(2)),-u.1(gamma,(mu.1-mu.0)/sqrt(2)))
delta.beta <- runif(length(beta)); delta.gamma <- runif(length(gamma))
ts <- seq(0,1,len=20)
D.1s <- sapply(ts, function(t)D.1(x.0,x.1,mu.0,mu.1,beta+t*delta.beta,gamma+t*delta.gamma))
plot(ts,D.1s,type='l')
abline(a=D.1s[1], b=t(deriv)%*%c(delta.beta,delta.gamma))

taylor.term.0 <- function(mu.0,mu.1,Sigma,beta,gamma) pnorm(beta%*%(mu.1-mu.0)/sqrt(2*t(beta)%*%Sigma%*%beta)) - pnorm(gamma%*%(mu.1-mu.0)/sqrt(2*t(gamma)%*%Sigma%*%gamma))
taylor.term.1 <- function(x.0,x.1,mu.0,mu.1,Sigma,beta,gamma) {
    partials.y <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.0-mu.1)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    partials.x <- sapply(list(beta=beta,gamma=gamma), function(beta)
        mean(pnorm(t(beta)%*%(x.1-mu.0)/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))) )
    unname(partials.y['gamma'] - partials.y['beta'] + partials.x['beta'] - partials.x['gamma']) - 2*taylor.term.0(mu.0,mu.1,Sigma,beta,gamma)
}
start <- Sys.time()
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
p <- 3
## n <- 1e1
## pairs <- replicate(1e2, {
mu.0 <- c(1,runif(p))
mu.1 <- c(1,runif(p))
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
Sigma <- rbind(0,cbind(0,Sigma))
beta <- runif(p+1)
gamma <- beta#runif(p+1)
## beta <- c(-1/2*(mu.1[-1]%*%solve(Sigma[-1,-1])%*%mu.1[-1]-mu.2[-1]%*%solve(Sigma[-1,-1])%*%mu.2[-1]),solve(Sigma[-1,-1])%*%(mu.1-mu.0))
ns <- round(seq(10,5e2,len=20))
by.n <- sapply(ns, function(n) {
## by.n <- mclapply(ns, mc.cores=detectCores()-3,FUN=function(n) {
    ## print(max(1,detectCores()-3))
    print(n)
    ## diffs <- mclapply(1:1e2, mc.cores=4,FUN=function(jj) {
    diffs <- replicate(1e2, {
        x.0 <- t(rmvnorm(n,mu.0,Sigma))
        x.1 <- t(rmvnorm(n,mu.1,Sigma))
        beta.hat <- beta+rnorm(p+1)/sqrt(n)
        gamma.hat <- gamma+rnorm(p+1)/sqrt(n)
        ## t(rowMeans(apply(x.1-mu.0,2,function(w)u.1(beta,w))) - rowMeans(apply(x.0-mu.1,2,function(w)u.1(beta,w))) - 2*u.1(beta,(mu.1-mu.0)/sqrt(2))) %*% (beta.hat-gamma.hat)
        qform <- as.numeric(t(beta)%*%Sigma%*%beta)
        deriv.1 <- colMeans(as.numeric(dnorm(t(beta)%*%(x.0-mu.1)/sqrt(qform)))*t((Sigma%*%beta)%*%t(beta)%*%(x.0-mu.1)/qform - x.0 + mu.1))    -            colMeans(as.numeric(dnorm(t(beta)%*%(x.1-mu.0)/sqrt(qform)))*t((Sigma%*%beta)%*%t(beta)%*%(x.1-mu.0)/qform - x.1 + mu.0))  + sqrt(2)*as.numeric(dnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*qform)))* t((Sigma%*%beta)%*%t(beta)%*%(mu.1-mu.0)/qform - (mu.1-mu.0))
        ## rowMeans(sapply(1:n, function(i)as.numeric(dnorm(t(beta)%*%(x.1[,i]-mu.0)/sqrt(qform)))*(as.numeric(t(beta)%*%(x.1[,i]-mu.0)/sqrt(qform))*Sigma%*%beta-x.1[,i]+mu.0)))
        deriv.1%*%(beta.hat-gamma.hat)/sqrt(qform)
    })
    simplify2array(diffs)
})
plot(ns,ns^(1)*apply(by.n,2,sd))
Sys.time() - start
sds <- apply(by.n,2,sd)
lm0 <- lm(log(sds) ~ log(ns))
coef(lm0)
plot(ns,log(sds))
lines(ns,fitted(lm0),col=2)




start <- Sys.time()
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
p <- 3
## n <- 1e1
## pairs <- replicate(1e2, {
mu.0 <- c(1,runif(p))
mu.1 <- c(1,runif(p))
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
Sigma <- rbind(0,cbind(0,Sigma))
beta <- runif(p+1)
gamma <- beta#runif(p+1)
## beta <- c(-1/2*(mu.1[-1]%*%solve(Sigma[-1,-1])%*%mu.1[-1]-mu.2[-1]%*%solve(Sigma[-1,-1])%*%mu.2[-1]),solve(Sigma[-1,-1])%*%(mu.1-mu.0))
ns <- round(seq(10,5e2,len=20))
by.n <- sapply(ns, function(n) {
## by.n <- mclapply(ns, mc.cores=detectCores()-3,FUN=function(n) {
    ## print(max(1,detectCores()-3))
    print(n)
    ## diffs <- mclapply(1:1e2, mc.cores=4,FUN=function(jj) {
    diffs <- replicate(1e2, {
        x.0 <- t(rmvnorm(n,mu.0,Sigma))
        x.1 <- t(rmvnorm(n,mu.1,Sigma))
        beta.hat <- beta+rnorm(p+1)/sqrt(n)
        gamma.hat <- gamma+rnorm(p+1)/sqrt(n)
        ## t(rowMeans(apply(x.1-mu.0,2,function(w)u.1(beta,w))) - rowMeans(apply(x.0-mu.1,2,function(w)u.1(beta,w))) - 2*u.1(beta,(mu.1-mu.0)/sqrt(2))) %*% (beta.hat-gamma.hat)
        qform <- as.numeric(t(beta)%*%Sigma%*%beta)
        deriv.1 <- colMeans(as.numeric(dnorm(t(beta)%*%(x.0-mu.1)/sqrt(qform)))*t((Sigma%*%beta)%*%t(beta)%*%(x.0-mu.1)/qform - x.0 + mu.1))    -            colMeans(as.numeric(dnorm(t(beta)%*%(x.1-mu.0)/sqrt(qform)))*t((Sigma%*%beta)%*%t(beta)%*%(x.1-mu.0)/qform - x.1 + mu.0))  + sqrt(2)*as.numeric(dnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*qform)))* t((Sigma%*%beta)%*%t(beta)%*%(mu.1-mu.0)/qform - (mu.1-mu.0))
        ## rowMeans(sapply(1:n, function(i)as.numeric(dnorm(t(beta)%*%(x.1[,i]-mu.0)/sqrt(qform)))*(as.numeric(t(beta)%*%(x.1[,i]-mu.0)/sqrt(qform))*Sigma%*%beta-x.1[,i]+mu.0)))
        deriv.1%*%(beta.hat-gamma.hat)/sqrt(qform)
    })
    simplify2array(diffs)
})
plot(ns,ns^(1)*apply(by.n,2,sd))
Sys.time() - start
sds <- apply(by.n,2,sd)
lm0 <- lm(log(sds) ~ log(ns))
coef(lm0)
plot(ns,log(sds))
lines(ns,fitted(lm0),col=2)


## 8ea formula for E(Xphi(X))
B <- 1e4
mu <- runif(1)
sigma <- runif(1)
x <- rnorm(B,mu,sigma)
mean(x*dnorm(x))
integrate(function(x)x*dnorm(x)*dnorm((x-mu)/sigma)/sigma,-10,10)$val
integrate(function(x)x/2/pi/sigma*exp(-1/2*((1+1/sigma^2)*x^2-2*x*mu/sigma^2+mu^2/sigma^2)),-10,10)$val
integrate(function(x)mu/2/pi/sigma^3/(1+1/sigma^2)*exp(-1/2*((1+1/sigma^2)*x^2-2*x*mu/sigma^2+mu^2/sigma^2)),-10,10)$val
## mu/4/sqrt(pi)/sigma^2*exp(-mu^2/4/sigma^2)
## mu/sqrt(2*pi)*sigma/(1+sigma^2)^(3/2)*exp(-mu^2/2/sigma^2/(1+sigma^2))
mu/sigma/2/pi/(1+sigma^2)*exp(-mu^2/2/(1+sigma^2))*sqrt(2*pi*sigma^2/(1+sigma^2))
mu/sqrt(2*pi)/(1+sigma^2)^(3/2)*exp(-mu^2/2/(1+sigma^2))



## 8eb
start <- Sys.time()
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
p <- 3
## n <- 1e1
## pairs <- replicate(1e2, {
mu.0 <- c(1,runif(p))
mu.1 <- c(1,runif(p))
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
Sigma <- rbind(0,cbind(0,Sigma))
beta <- runif(p+1)
gamma <- beta#runif(p+1)
## beta <- c(-1/2*(mu.1[-1]%*%solve(Sigma[-1,-1])%*%mu.1[-1]-mu.2[-1]%*%solve(Sigma[-1,-1])%*%mu.2[-1]),solve(Sigma[-1,-1])%*%(mu.1-mu.0))
a <- runif(p+1)
qform <- as.numeric(t(beta)%*%Sigma%*%beta)
plim.vec <- ((Sigma%*%beta)%*% (t(beta)%*%(mu.0-mu.1)/qform) + mu.1 - mu.0) %*% exp(-1/4*(t(beta)%*%(mu.0-mu.1))^2/qform)/4/sqrt(pi)
plim <- a%*%plim.vec
plim <- exp(-1/4*(t(beta)%*%(mu.0-mu.1))^2/qform)/4/sqrt(pi) * t(beta)%*%(mu.0-mu.1) #1
## plim.vec <- (Sigma%*%beta)%*% (t(beta)%*%(mu.0-mu.1)/qform/4/sqrt(pi)) %*% exp(-1/4*(t(beta)%*%(mu.0-mu.1))^2/qform) #2
## plim.vec <- c(exp(-1/4*(t(beta)%*%(mu.0-mu.1))^2/qform)/4/sqrt(pi)) * (mu.0-mu.1) #3
## plim <- a%*%plim.vec
ns <- round(seq(10,5e2,len=20))
by.n <- sapply(ns, function(n) {
## by.n <- mclapply(ns, mc.cores=detectCores()-3,FUN=function(n) {
    ## print(max(1,detectCores()-3))
    print(n)
    ## diffs <- mclapply(1:1e2, mc.cores=4,FUN=function(jj) {
    diffs <- replicate(1e2, {
        x.0 <- t(rmvnorm(n,mu.0,Sigma))
        x.1 <- t(rmvnorm(n,mu.1,Sigma))
        ## obs.vec <- colMeans(as.numeric(dnorm(t(beta)%*%(x.0-mu.1)/sqrt(qform)))*t((Sigma%*%beta)%*%t(beta)%*%(x.0-mu.1)/qform - x.0 + mu.1))
        ## obs <- a%*%obs.vec
        obs <- mean(dnorm(t(beta)%*%(x.0-mu.1)/sqrt(qform)) * t(beta)%*%(x.0-mu.1)) #1
        ## obs.vec <- mean(dnorm(t(beta)%*%(x.0-mu.1)/sqrt(qform))*t(beta)%*%(x.0-mu.1)/qform)*Sigma%*%beta #2
        ## obs.vec <- colMeans(as.numeric(dnorm(t(beta)%*%(x.0-mu.1)/sqrt(qform))) * t(x.0-mu.1)) #3
        ## obs <- a%*%obs.vec
        obs-plim
    })
    simplify2array(diffs)
})
matplot(t(by.n),pch=1,col=1)

plot(ns,colMeans(abs(by.n)))




require(mvtnorm)
p <- 3
mu <- runif(p)
## mu.1 <- runif(p)
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
u.1 <- function(beta,w) {
    qform <- as.numeric(t(beta)%*%Sigma%*%beta)
    as.numeric(-dnorm(t(beta)%*%w/sqrt(qform))/sqrt(qform))*(as.numeric(t(beta)%*%w/qform)*Sigma%*%beta-w)
}
B <- 1e5
a <- runif(p)
w <- t(rmvnorm(B,mu,Sigma))
obs <- apply(w,2,function(w.i)u.1(beta,w.i))
try <- u.1(beta,mu/sqrt(2))
hist(a%*%obs)
abline(v=a%*%try,col=2)



a <- runif(1)
c <- runif(1)
d <- runif(1)
integrate(function(u)dnorm(u+a)*dnorm(c*u+d),-7,7)$val
1/2/pi*exp(1/2*((a^2+2*a*c*d+c^2*d^2)/(c^2+1)-a^2-d^2))*sqrt(2*pi/(c^2+1))
1/sqrt(2*pi*(c^2+1))*exp(1/2*(-c^2/(c^2+1)*a^2+2*c*d*a/(c^2+1)-d^2/(c^2+1)))
1/sqrt(2*pi*(c^2+1))*integrate(function(a)exp(1/2*(-c^2/(c^2+1)*a^2+2*c*d*a/(c^2+1)-d^2/(c^2+1))),-7,7)$val
integrate(function(a)exp(-1/2*c^2/(c^2+1)*(a^2-2*d*a/c)),-7,7)$val*exp(-1/2*d^2/(c^2+1))
sqrt(2*pi*(1+1/c^2))


integrate(function(u)pnorm(u+a)*dnorm(c*u+d),-30,30)$val
## 1/sqrt(2*pi*(c^2+1))*integrate(function(a)exp(1/2*(-c^2/(c^2+1)*a^2+2*c*d*a/(c^2+1)-d^2/(c^2+1))),-7,7)$val
1/c*pnorm(c/sqrt(1+c^2)*(a-d/c))

## did not know this? integral of Phi(x) on [-a,a] is a. since Phi(-x)=1-Phi(x) so just integrating 1.
a <- runif(1,0,10); integrate(pnorm,-a,a)$val - a

## 8ec
require(mvtnorm)
B <- 1e6
mu <- runif(2)
Sigma <- matrix(runif(4),2,2)
Sigma <- Sigma%*%t(Sigma)
Bs <- round(seq(1e2,1e5,len=50))
try <- dnorm(mu[1]/sqrt(2*Sigma[1,1]))*(mu[1]/Sigma[1,1]*Sigma[1,2] - mu[2])/sqrt(2)
by.B <- sapply(Bs, function(B) {
    x <- t(rmvnorm(B,mu,Sigma))
    obs <- dnorm(x[1,]/sqrt(Sigma[1,1]))*(x[1,]/Sigma[1,1]*Sigma[1,2] - x[2,])
    obs - try
})
## hist(obs)
## abline(v=try,col=2)
## abline(v=mean(obs),col=3)
means <- sapply(by.B,mean)
plot(Bs,means)
abline(h=0)
means


require(mvtnorm)
B <- 1e6
mu <- runif(2)
Sigma <- matrix(runif(4),2,2)
Sigma <- Sigma%*%t(Sigma)
Bs <- round(seq(1e2,1e5,len=50))
try <- dnorm(mu[1]/sqrt(2))*(mu[1]/sqrt(Sigma[1,1])*Sigma[1,2] - mu[2])/sqrt(2)
by.B <- sapply(Bs, function(B) {
    x <- t(rmvnorm(B,mu,Sigma))
    obs <- dnorm(x[1,])*(x[1,]/sqrt(Sigma[1,1])*Sigma[1,2] - x[2,])
    obs - try
})
## hist(obs)
## abline(v=try,col=2)
## abline(v=mean(obs),col=3)
means <- sapply(by.B,mean)
plot(Bs,means)
abline(h=0)
means

## 8ed
a <- runif(1)
b <- runif(1)
c <- runif(1)
d <- runif(1)
integrate(function(u)dnorm(a*u+b)*dnorm(c*u+d),-10,10)$val
1/sqrt(2*pi)/sqrt(a^2+c^2)*exp(1/2*((a*b+c*d)^2/(a^2+c^2)-(b^2+d^2)))
1/2/pi*exp(-(b^2+d^2)/2)*integrate(function(u)exp(-1/2*(a^2+c^2)*(u^2+2*u*(a*b+c*d)/(a^2+c^2))),-10,10)$val
dnorm((a*d-b*c)/sqrt(a^2+c^2))/sqrt(a^2+c^2)

## 9 variance of linear term

## 9a
B <- 1e4
a <- runif(1)
b <- runif(1)
z <- rnorm(B)
mean(pnorm(a*z+b)*z)
a*dnorm(b/sqrt(a^2+1))/sqrt(a^2+1)
c1 <- runif(1)
c2 <- runif(1)
z1 <- rnorm(B)
z2 <- runif(1)
mean(pnorm(c1*z1+c2*z2)*c1*z1)
c1^2*dnorm(c2*z2/sqrt(c1^2+1))/sqrt(c1^2+1)

c1 <- runif(1)
c2 <- runif(1)
z1 <- rnorm(B)
z2 <- rnorm(B)
mean(pnorm(c1*z1+c2*z2)*c1*z1)
c1^2/sqrt(1+c1^2+c2^2)*dnorm(0)

B <- 1e5
c1 <- runif(1)
c2 <- runif(1)
z1 <- rnorm(B)
z2 <- rnorm(B)
mu1 <- runif(1)
mean(pnorm(c1*z1+c2*z2+mu1)*z1)
c1/sqrt(1+c1^2+c2^2)*dnorm(mu1/sqrt(1+c1^2+c2^2))
mean(pnorm(c1*z1+c2*z2+mu1))
pnorm(mu1/sqrt(1+c1^2+c2^2))

d <- runif(1);c <- runif(1)
u <- rnorm(B,-d/c,1/c)
mean(pnorm(a*u+b))
pnorm((b*c-a*d)/sqrt(a^2+c^2))
mean(dnorm(a*u+b))
c*dnorm((b*c-a*d)/sqrt(a^2+c^2))/sqrt(a^2+c^2)

mu <- runif(1)
sigma <- runif(1)
a <- runif(1)
b <- runif(1)
u <- rnorm(B,mu,sigma)
mean(pnorm(a*u+b))
pnorm((b+a*mu)/sqrt(1+a^2*sigma^2))
mean(dnorm(a*u+b))
dnorm((b+a*mu)/sqrt(1+a^2*sigma^2))/sqrt(1+a^2*sigma^2)

B <- 1e4
c1 <- runif(1)
d1 <- runif(1)
c2 <- runif(1)
d2 <- runif(1)
mu1 <- runif(1)
mu2 <- runif(1)
mean(pnorm(c1*z1+c2*z2+mu1)*(d1*z1+d2*z2+mu2))
(d1*c1+d2*c2)/sqrt(1+c1^2+c2^2)*dnorm(mu1/sqrt(1+c1^2+c2^2)) + mu2*pnorm(mu1/sqrt(1+c1^2+c2^2))

s <- sqrt(det(Sigma))
t <- sqrt(Sigma[1,1]+Sigma[2,2]+2*s)
p <- root <- 1/t*(Sigma+s*diag(2))
p[1,2]*sum(diag(p)) - Sigma[1,2]
1+p[1,1]^2+p[1,2]^2
1+Sigma[1,1]

require(mvtnorm)
B <- 1e4
mu <- runif(2)
Sigma <- matrix(runif(4),2,2)
Sigma <- Sigma%*%t(Sigma)
xy <- t(rmvnorm(B,mu,Sigma))
x <- xy[1,]; y <- xy[2,]
mean(pnorm(x)*y)
Sigma[1,2]/sqrt(1+Sigma[1,1])*dnorm(mu[1]/sqrt(1+Sigma[1,1])) + mu[2]*pnorm(mu1/sqrt(1+Sigma[1,1]))

require(mvtnorm)
B <- 1e4
p <- 4
a <- runif(p)
beta <- runif(p)
mu.w <- runif(p)
Sigma.w <- matrix(runif(p^2),nrow=p)
Sigma.w <- Sigma.w%*%t(Sigma.w)
mu <- rbind(beta,a)%*%mu.w
Sigma <- matrix(c(t(beta)%*%Sigma.w%*%beta,t(beta)%*%Sigma.w%*%a,t(beta)%*%Sigma.w%*%a,(a)%*%Sigma.w%*%a),nrow=2)
xy <- t(rmvnorm(B,mu,Sigma))
x <- xy[1,]; y <- xy[2,]
mean(pnorm(x)*y)
Sigma[1,2]/sqrt(1+Sigma[1,1])*dnorm(mu[1]/sqrt(1+Sigma[1,1])) + mu[2]*pnorm(mu[1]/sqrt(1+Sigma[1,1]))

require(mvtnorm)
B <- 1e3
p <- 4
mu <- runif(p)
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
a <- runif(p)
beta <- runif(p)
w <- t(rmvnorm(B,mu,Sigma))
mean(pnorm(t(beta)%*%w) * (t(a)%*%w))
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
t(a)%*%((Sigma%*%beta/sqrt(1+quad))%*%dnorm(t(beta)%*%mu/sqrt(1+quad)) +  mu*c(pnorm(t(beta)%*%mu/sqrt(1+quad))))


plot.matrix.deriv <- function(f,deriv,dim.x,dim.f) {
    x0 <- runif(dim.x)
    delta <- matrix(runif(dim.x),ncol=1)
    a <- matrix(runif(dim.f),ncol=1)
    ts <- seq(0,1,len=1e2)
    fs <- sapply(ts,function(t)t(a)%*%f(x0+t*delta))
    plot(ts,fs,type='l')
    deriv.scalar <- t(a)%*%deriv(x0)%*%delta
    abline(a=fs[1],b=deriv.scalar,col=2,lty=2)
}
dim.x <- sample(1:10,1)
dim.f <- sample(1:10,1)
A <- matrix(runif(dim.x*dim.f),ncol=dim.x)
f <- function(x)A%*%x
plot.matrix.deriv(f,function(x)A,dim.x,dim.f)

p <- 4
mu <- matrix(runif(p),ncol=1)
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
## a <- runif(p)
beta <- matrix(runif(p),ncol=1)
## w <- t(rmvnorm(B,mu,Sigma))
## mean(pnorm(t(beta)%*%w) * (t(a)%*%w))
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
f <- function(beta) {
        quad <- as.numeric(t(beta)%*%Sigma%*%beta)
        (Sigma%*%beta/sqrt(1+quad))%*%dnorm(t(beta)%*%mu/sqrt(1+quad)) +  mu*c(pnorm(t(beta)%*%mu/sqrt(1+quad)))
}
## f(beta)
deriv <- function(beta) {
    quad <- as.numeric(t(beta)%*%Sigma%*%beta)
    A <- mu-Sigma%*%beta%*%t(beta)%*%mu/(1+quad)
    B <- t(mu)/sqrt(1+quad)-(t(beta)%*%mu*(1+quad)^(-3/2))%*%(t(beta)%*%Sigma)
    C <- Sigma/sqrt(1+quad) - Sigma%*%beta%*%t(beta)%*%Sigma*(1+quad)^(-3/2)
    as.numeric(dnorm(t(beta)%*%mu/sqrt(1+quad))) * (A%*%B+C)
    ##     1/sqrt(1+quad)*mu%*%t(mu)   -
    ##     as.numeric(t(beta)%*%mu/(1+quad)^(3/2)) * (mu%*%t(beta)%*%Sigma + Sigma%*%beta%*%t(mu))   +
    ##     as.numeric((t(beta)%*%mu)^2/(1+quad)^(5/2) - 1/(1+quad)^(3/2))*(Sigma%*%beta)%*%t(Sigma%*%beta)    -
    ##     Sigma/sqrt(1+quad)
    ## )
}
plot.matrix.deriv(f,deriv,p,p)
ddd




## f <- function(beta) {
##         quad <- as.numeric(t(beta)%*%Sigma%*%beta)
##         mu%*%pnorm(t(beta)%*%mu/sqrt(1+quad))
##         ## (Sigma%*%beta/sqrt(1+quad))%*%dnorm(t(beta)%*%mu/sqrt(1+quad)) +  mu*c(pnorm(t(beta)%*%mu/sqrt(1+quad)))
## }
## deriv <- function(beta) {
##         quad <- as.numeric(t(beta)%*%Sigma%*%beta)
##         A <- mu%*%dnorm(t(beta)%*%mu/sqrt(1+quad))
##         B <- (1/sqrt(1+quad))%*%t(mu) - (1+quad)^(-3/2)*(t(beta)%*%mu)%*%(t(beta)%*%Sigma)
##         A%*%B    
## }
## plot.matrix.deriv(f,deriv,p,p)
## ## f <- function(beta) {
## ##         quad <- as.numeric(t(beta)%*%Sigma%*%beta)
## ##         mu%*%pnorm(t(beta)%*%mu * 1/sqrt(1+quad))
## ##         ## (Sigma%*%beta/sqrt(1+quad))%*%dnorm(t(beta)%*%mu/sqrt(1+quad)) +  mu*c(pnorm(t(beta)%*%mu/sqrt(1+quad)))
## ## }
## ## deriv <- function(beta) {
## ##         quad <- as.numeric(t(beta)%*%Sigma%*%beta)
## ##         A <- mu%*%dnorm(t(beta)%*%mu/sqrt(1+quad))
## ##         B <-  (1/sqrt(1+quad))%*%t(mu) -  (1+quad)^(-3/2)*(t(beta)%*%mu)%*%(t(beta)%*%Sigma)
## ##         A%*%B    
## ## }
## ## plot.matrix.deriv(f,deriv,p,p)
## dd

## ## ## f <- function(beta)t(beta)%*%Sigma%*%beta
## ## ## deriv <- function(beta)2*t(beta)%*%Sigma
## ## ## plot.matrix.deriv(f,deriv,p,1)
## ## f <- function(beta) {
## ##         quad <- as.numeric(t(beta)%*%Sigma%*%beta)
## ##         mu%*%pnorm(t(beta)%*%mu)
## ##         ## (Sigma%*%beta/sqrt(1+quad))%*%dnorm(t(beta)%*%mu/sqrt(1+quad)) +  mu*c(pnorm(t(beta)%*%mu/sqrt(1+quad)))
## ## }
## ## deriv <- function(beta) {
## ##         quad <- as.numeric(t(beta)%*%Sigma%*%beta)
## ##         A <- mu%*%dnorm(t(beta)%*%mu)
## ##         B <- t(mu)# - (1+quad)^(-3/2)*(t(beta)%*%mu)%*%(t(beta)%*%Sigma)
## ##         A%*%B    
## ## }
## ## plot.matrix.deriv(f,deriv,p,p)
## f <- function(beta) {
##         quad <- as.numeric(t(beta)%*%Sigma%*%beta)
##         ## mu%*%pnorm(t(beta)%*%mu + 1/sqrt(1+quad))
##         ## (Sigma%*%beta/sqrt(1+quad))%*%dnorm(t(beta)%*%mu/sqrt(1+quad)) +  mu*c(pnorm(t(beta)%*%mu/sqrt(1+quad)))
##         Sigma%*%beta/sqrt(1+quad)
## }
## deriv <- function(beta) {
##         quad <- as.numeric(t(beta)%*%Sigma%*%beta)
##         A <- Sigma/sqrt(1+quad)
##         B <- -Sigma%*%beta%*%t(beta)%*%Sigma * as.numeric((1+quad)^(-3/2))
##         A+B    
## }
## plot.matrix.deriv(f,deriv,p,p)


## f <- function(beta) {
##         quad <- as.numeric(t(beta)%*%Sigma%*%beta)
##         dnorm(t(beta)%*%mu/sqrt(1+quad))
##         ## (Sigma%*%beta/sqrt(1+quad))%*%dnorm(t(beta)%*%mu/sqrt(1+quad)) +  mu*c(pnorm(t(beta)%*%mu/sqrt(1+quad)))
## }
## deriv <- function(beta) {
##     quad <- as.numeric(t(beta)%*%Sigma%*%beta)
##     A <-   -dnorm(t(beta)%*%mu/sqrt(1+quad))*t(beta)%*%mu/sqrt(1+quad)
##     B <- (1/sqrt(1+quad))%*%t(mu) - (1+quad)^(-3/2)*(t(beta)%*%mu)%*%(t(beta)%*%Sigma)
##     A%*%B    
## }
## plot.matrix.deriv(f,deriv,p,1)
## f <- function(beta) {
##         quad <- as.numeric(t(beta)%*%Sigma%*%beta)
##                 (Sigma%*%beta/sqrt(1+quad)) %*% dnorm(t(beta)%*%mu/sqrt(1+quad))
##         ## (Sigma%*%beta/sqrt(1+quad))%*%dnorm(t(beta)%*%mu/sqrt(1+quad)) +  mu*c(pnorm(t(beta)%*%mu/sqrt(1+quad)))
## }
## deriv <- function(beta) {
##     quad <- as.numeric(t(beta)%*%Sigma%*%beta)
##     A <-   -dnorm(t(beta)%*%mu/sqrt(1+quad))*t(beta)%*%mu/sqrt(1+quad)
##     B <- (1/sqrt(1+quad))%*%t(mu) - (1+quad)^(-3/2)*(t(beta)%*%mu)%*%(t(beta)%*%Sigma)
##     C <- Sigma/sqrt(1+quad)
##     D <- -Sigma%*%beta%*%t(beta)%*%Sigma * as.numeric((1+quad)^(-3/2))
##     ## A+B    
##     (Sigma%*%beta/sqrt(1+quad)) %*% (A%*%B) +
##         as.numeric(dnorm(t(beta)%*%mu/sqrt(1+quad)))*(C+D)
## }
## plot.matrix.deriv(f,deriv,p,p)

## 26 Ciy Hall Mall 02155
## for Dr. Yang
## 781 292 7700
## 1177 Providence Hwy Norwood MA 02062

require(mvtnorm)
p <- 4
B <- 1e3
mu <- matrix(runif(p),ncol=1)
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
a <- runif(p)
beta <- matrix(runif(p),ncol=1)
w <- t(rmvnorm(B,mu,Sigma))
## mean(pnorm(t(beta)%*%w) * (t(a)%*%w))
observed <- dnorm(t(beta)%*%w)*(t(a)%*%w)^2
## quad <- as.numeric(t(beta)%*%Sigma%*%beta)
## f <- function(beta) {
##         quad <- as.numeric(t(beta)%*%Sigma%*%beta)
##         (Sigma%*%beta/sqrt(1+quad))%*%dnorm(t(beta)%*%mu/sqrt(1+quad)) +  mu*c(pnorm(t(beta)%*%mu/sqrt(1+quad)))
## }
## f(beta)
deriv <- function(beta) {
    quad <- as.numeric(t(beta)%*%Sigma%*%beta)
    A <- mu-Sigma%*%beta%*%t(beta)%*%mu/(1+quad)
    B <- t(mu)/sqrt(1+quad)-(t(beta)%*%mu*(1+quad)^(-3/2))%*%(t(beta)%*%Sigma)
    C <- Sigma/sqrt(1+quad) - Sigma%*%beta%*%t(beta)%*%Sigma*(1+quad)^(-3/2)
    as.numeric(dnorm(t(beta)%*%mu/sqrt(1+quad))) * (A%*%B+C)
    ##     1/sqrt(1+quad)*mu%*%t(mu)   -
    ##     as.numeric(t(beta)%*%mu/(1+quad)^(3/2)) * (mu%*%t(beta)%*%Sigma + Sigma%*%beta%*%t(mu))   +
    ##     as.numeric((t(beta)%*%mu)^2/(1+quad)^(5/2) - 1/(1+quad)^(3/2))*(Sigma%*%beta)%*%t(Sigma%*%beta)    -
    ##     Sigma/sqrt(1+quad)
    ## )
}
hist(observed)
abline(v=t(a)%*%deriv(beta)%*%a,col=2)
abline(v=mean(observed),col=3)

## 9aa
require(mvtnorm)
p <- 4
B <- 1e3
mu.bar <- matrix(runif(p),ncol=1)
Sigma.bar <- matrix(runif(p^2),nrow=p)
Sigma.bar <- Sigma.bar%*%t(Sigma.bar)
a <- runif(p)
beta <- matrix(runif(p),ncol=1)
quad.bar <- as.numeric(t(beta)%*%Sigma.bar%*%beta)
mu <- sqrt(2)*mu.bar/sqrt(quad.bar)
Sigma <- 2*Sigma.bar/quad.bar
v <- t(rmvnorm(B,mu,Sigma))
## mean(pnorm(t(beta)%*%w) * (t(a)%*%w))
observed <- dnorm(t(beta)%*%v)*(t(a)%*%v)^2
## quad <- as.numeric(t(beta)%*%Sigma%*%beta)
## f <- function(beta) {
##         quad <- as.numeric(t(beta)%*%Sigma%*%beta)
##         (Sigma%*%beta/sqrt(1+quad))%*%dnorm(t(beta)%*%mu/sqrt(1+quad)) +  mu*c(pnorm(t(beta)%*%mu/sqrt(1+quad)))
## }
## f(beta)
deriv <- function(beta) {
    quad <- as.numeric(t(beta)%*%Sigma%*%beta)
    A <- mu-Sigma%*%beta%*%t(beta)%*%mu/(1+quad)
    B <- t(mu)/sqrt(1+quad)-(t(beta)%*%mu*(1+quad)^(-3/2))%*%(t(beta)%*%Sigma)
    C <- Sigma/sqrt(1+quad) - Sigma%*%beta%*%t(beta)%*%Sigma*(1+quad)^(-3/2)
    as.numeric(dnorm(t(beta)%*%mu/sqrt(1+quad))) * (A%*%B+C)
    ##     1/sqrt(1+quad)*mu%*%t(mu)   -
    ##     as.numeric(t(beta)%*%mu/(1+quad)^(3/2)) * (mu%*%t(beta)%*%Sigma + Sigma%*%beta%*%t(mu))   +
    ##     as.numeric((t(beta)%*%mu)^2/(1+quad)^(5/2) - 1/(1+quad)^(3/2))*(Sigma%*%beta)%*%t(Sigma%*%beta)    -
    ##     Sigma/sqrt(1+quad)
    ## )
}
hist(observed)
abline(v=t(a)%*%deriv(beta)%*%a,col=2)
abline(v=mean(observed),col=3)

deriv(beta)
A <- sqrt(2/quad.bar)*mu.bar-Sigma.bar%*%beta%*%(sqrt(2)/3*t(beta)%*%mu.bar/sqrt(quad.bar)*2/quad.bar)
B <- sqrt(2/3/quad.bar)*t(mu.bar)-(sqrt(2/quad.bar)*t(beta)%*%mu.bar/sqrt(27)*2/quad.bar)%*%(t(beta)%*%Sigma.bar)
C <- as.numeric(2/sqrt(3)/quad.bar)*Sigma.bar - (2/quad.bar)^2/sqrt(27)*Sigma.bar%*%beta%*%t(beta)%*%Sigma.bar
as.numeric(dnorm(sqrt(2/3/quad.bar)*t(beta)%*%mu.bar))*(A%*%B+C)
A <- mu.bar  - (Sigma.bar%*%beta)%*%(2/3*t(beta)%*%mu.bar/quad.bar)
C <- Sigma.bar - (2/3/quad.bar)*Sigma.bar%*%beta%*%t(beta)%*%Sigma.bar
2*sqrt(2*pi)*as.numeric(1/quad.bar/sqrt(6*pi)*dnorm(sqrt(2/3/quad.bar)*t(beta)%*%mu.bar))*(A%*%t(A) + C)



## 9ab
require(mvtnorm)
## derivative of E(pnorm(t(beta)%*%w)*w) when w is N(mu,Sigma)
deriv <- function(beta,mu,Sigma) {
    quad <- as.numeric(t(beta)%*%Sigma%*%beta)
    A <- mu-Sigma%*%beta%*%t(beta)%*%mu/(1+quad)
    B <- t(mu)/sqrt(1+quad)-(t(beta)%*%mu*(1+quad)^(-3/2))%*%(t(beta)%*%Sigma)
    C <- Sigma/sqrt(1+quad) - Sigma%*%beta%*%t(beta)%*%Sigma*(1+quad)^(-3/2)
    as.numeric(dnorm(t(beta)%*%mu/sqrt(1+quad))) * (A%*%B+C)
}
p <- 4
B <- 1e4
mu <- mu.bar <- matrix(runif(p),ncol=1)
Sigma.bar <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma.bar <- Sigma.bar%*%t(Sigma.bar)
a <- runif(p)
beta <- matrix(runif(p),ncol=1)
quad.bar <- as.numeric(t(beta)%*%Sigma.bar%*%beta)
## mu <- sqrt(2)*mu.bar/sqrt(quad.bar)
## Sigma <- 2*Sigma.bar/quad.bar
w <- t(rmvnorm(B,mu.bar,Sigma.bar))
## mean(pnorm(t(beta)%*%w) * (t(a)%*%w))
observed <- dnorm(t(beta)%*%w/sqrt(quad.bar))^2*(t(a)%*%w/sqrt(quad.bar))^2
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
hist(observed)
abline(v=1/(2*sqrt(2*pi))*t(a)%*%deriv(beta,mu=sqrt(2)*mu/sqrt(quad),Sigma=2*Sigma/quad)%*%a,col=2)
abline(v=mean(observed),col=3,lty=2)

1/(2*sqrt(2*pi))*deriv(beta,mu=sqrt(2)*mu/sqrt(quad),Sigma=2*Sigma/quad)
A <- mu.bar  - (Sigma.bar%*%beta)%*%(2/3*t(beta)%*%mu.bar/quad.bar)
C <- Sigma.bar - (2/3/quad.bar)*Sigma.bar%*%beta%*%t(beta)%*%Sigma.bar
as.numeric(1/quad.bar/sqrt(6*pi)*dnorm(sqrt(2/3/quad.bar)*t(beta)%*%mu.bar))*(A%*%t(A) + C)



## 9b

## 9ba
## derivative of E(dnorm(a*x)) when x is N(mu,sigma^2)
deriv <- function(a,mu,sigma)-a*f(a)/(1+a^2*sigma^2)*(mu^2/(1+a^2*sigma^2)+sigma^2)
mu <- runif(1)
sigma <- runif(1)
f <- function(a)dnorm(a*mu/sqrt(1+a^2*sigma^2))/sqrt(1+a^2*sigma^2)
as <- seq(0,4,len=30)
fs <- f(as)
plot(as,fs,type='l')
## curve(f)
## deriv <- function(a)-a*f(a)/(1+a^2*sigma^2)*(mu^2/(1+a^2*sigma^2)+sigma^2)
i <- sample(length(as),1)
m <- deriv(as[i],mu,sigma)
abline(a=fs[i]-m*as[i],b=m,lty=2,col=2)
abline(v=as[i])


## 9bb
## derivative of E(dnorm(a*x)) when x is N(mu,sigma^2)
deriv <- function(a,mu,sigma)-a*f(a)/(1+a^2*sigma^2)*(mu^2/(1+a^2*sigma^2)+sigma^2)
B <- 1e3
mu <- runif(1)
sigma <- runif(1)
a <- runif(1)
x <- rnorm(B,mu,sigma)
observed <- dnorm(a*x)*x^2
hist(observed)
try <- -1/a*deriv(a,mu,sigma)
abline(v=try,col=2)
abline(v=mean(observed),col=3,lty=2)

## 9bc
## deriv <- function(a,mu,sigma)-a*f(a)/(1+a^2*sigma^2)*(mu^2/(1+a^2*sigma^2)+sigma^2)
require(mvtnorm)
B <- 1e3
p <- 4
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
a <- runif(p)
w <- t(rmvnorm(B,mu,Sigma))
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
observed <- dnorm(sqrt(2)*t(beta)%*%w/sqrt(quad))/quad * (t(beta)%*%w/quad)^2 * as.numeric(t(a)%*%(Sigma%*%beta%*%t(beta)%*%Sigma)%*%a)
try <- dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(quad))/3/sqrt(3)*(1/3*(t(beta)%*%mu)^2/quad+1)/quad^2* as.numeric(t(a)%*%(Sigma%*%beta%*%t(beta)%*%Sigma)%*%a)
hist(observed)
abline(v=try,col=2)
abline(v=mean(observed),col=3,lty=2)


## 9c

## 9ca
require(mvtnorm)
B <- 1e3
p <- 4
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
a <- runif(p)
x <- t(rmvnorm(B,mu,Sigma))
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
observed <- dnorm(t(beta)%*%x)*(t(beta)%*%x) * (t(a)%*%x)
try <-dnorm(t(beta)%*%mu/sqrt(1+quad))*(1+quad)^(-3/2)*( (t(beta)%*%mu)*(t(mu)%*%a) + (1-(t(beta)%*%mu)^2/(1+quad))*(t(beta)%*%Sigma%*%a))
hist(observed)
abline(v=try,col=2)
abline(v=mean(observed),col=3,lty=2)




## 9cb
require(mvtnorm)
B <- 1e3
p <- 4
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
a <- runif(p)
w <- t(rmvnorm(B,mu,Sigma))
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
observed <- dnorm(sqrt(2)*t(beta)%*%w/sqrt(quad))/quad * (t(beta)%*%w/quad) * as.numeric(t(a)%*%(Sigma%*%beta)) * t(a)%*%w
try <- t(a)%*%Sigma%*%beta * dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(quad))*1/3/sqrt(3)/quad^2*(   t(beta)%*%mu * t(mu)%*%a + (1-2*(t(beta)%*%mu)^2/3/quad)*(t(beta)%*%Sigma%*%a))
hist(observed)
abline(v=try,col=2)
abline(v=mean(observed),col=3,lty=2)



## 9d combining
p <- 4
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
A.1 <- mu  - (Sigma%*%beta)%*%(2/3*t(beta)%*%mu/quad)
A.2 <- Sigma - (2/3/quad)*Sigma%*%beta%*%t(beta)%*%Sigma
A <- sqrt(2*pi)*as.numeric(1/quad/sqrt(6*pi)*dnorm(sqrt(2/3/quad)*t(beta)%*%mu))*(A.1%*%t(A.1) + A.2)
B <- as.numeric(dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(quad))/3/sqrt(3)*(1/3*(t(beta)%*%mu)^2/quad+1)/quad^2)*Sigma%*%beta%*%t(beta)%*%Sigma
C <-     as.numeric(dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(quad))*1/3/sqrt(3)/quad^2)*Sigma%*%beta %*%(   (t(beta)%*%mu) %*% t(mu)+ (1-2*(t(beta)%*%mu)^2/3/quad)%*%(t(beta)%*%Sigma))
A+B-C-t(C)
as.numeric(dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(quad))/sqrt(3)/quad) * (
    mu%*%t(mu) + Sigma - as.numeric(t(beta)%*%mu/quad) * (mu%*%t(beta)%*%Sigma + Sigma%*%beta%*%t(mu)) + as.numeric((t(beta)%*%mu/quad)^2-1/quad)*Sigma%*%beta%*%t(beta)%*%Sigma
)


## dbg
as.numeric(1/dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(quad))*3*sqrt(3)*quad) * A
A.1 <- sqrt(3)*mu  - (Sigma%*%beta)%*%(2/sqrt(3)*t(beta)%*%mu/quad)
A.2 <- 3*Sigma - (2/quad)*Sigma%*%beta%*%t(beta)%*%Sigma
A.1%*%t(A.1)+A.2

as.numeric(1/dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(quad))*3*sqrt(3)*quad) * B
Sigma%*%beta%*%t(beta)%*%Sigma/quad*as.numeric(1/3*(t(beta)%*%mu)^2/quad+1)

as.numeric(1/dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(quad))*3*sqrt(3)*quad) * C
(Sigma%*%beta/quad)%*% (t(beta)%*%mu%*%t(mu)  + (1-2*(t(beta)%*%mu)^2/3/quad)%*%t(beta)%*%Sigma  )



## 9e
require(mvtnorm)
p <- 3
mu <- runif(p)
## mu.1 <- runif(p)
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
u.1 <- function(beta,w) {
    qform <- as.numeric(t(beta)%*%Sigma%*%beta)
    as.numeric(-dnorm(t(beta)%*%w/sqrt(qform))/sqrt(qform))*(as.numeric(t(beta)%*%w/qform)*Sigma%*%beta-w)
}
B <- 1e3
a <- runif(p)
w <- t(rmvnorm(B,mu,Sigma))
u.1s <- apply(w,2,function(w.i)u.1(beta,w.i))
observed <- (t(a)%*%u.1s)^2
## try <- u.1(beta,mu/sqrt(2))
## hist(a%*%obs)
## abline(v=a%*%try,col=2)
moment.2 <- 1/sqrt(2*pi)*as.numeric(dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(quad))/sqrt(3)/quad) * (
    mu%*%t(mu) + Sigma - as.numeric(t(beta)%*%mu/quad) * (mu%*%t(beta)%*%Sigma + Sigma%*%beta%*%t(mu)) + as.numeric((t(beta)%*%mu/quad)^2-1/quad)*Sigma%*%beta%*%t(beta)%*%Sigma
)
try <- t(a)%*%moment.2%*%a
hist(observed)
abline(v=try,col=2)
abline(v=mean(observed),col=3,lty=2)

moment.1 <- u.1(beta,mu/sqrt(2))
try <- t(a)%*%moment.1
observed <- t(a)%*%u.1s
hist(observed)
abline(v=try,col=2)
abline(v=mean(observed),col=3,lty=2)


## 9f variance of u.1
require(mvtnorm)
p <- 3
mu <- runif(p)
## mu.1 <- runif(p)
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
u.1 <- function(beta,w) {
    qform <- as.numeric(t(beta)%*%Sigma%*%beta)
    as.numeric(-dnorm(t(beta)%*%w/sqrt(qform))/sqrt(qform))*(as.numeric(t(beta)%*%w/qform)*Sigma%*%beta-w)
}
moment.1 <- u.1(beta,mu/sqrt(2))
moment.2 <- 1/sqrt(2*pi)*as.numeric(dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(quad))/sqrt(3)/quad) * (
    mu%*%t(mu) + Sigma - as.numeric(t(beta)%*%mu/quad) * (mu%*%t(beta)%*%Sigma + Sigma%*%beta%*%t(mu)) + as.numeric((t(beta)%*%mu/quad)^2-1/quad)*Sigma%*%beta%*%t(beta)%*%Sigma
)
var.u.1 <- moment.2 - moment.1%*%t(moment.1)
aa <- (
    mu%*%t(mu)  - as.numeric(t(beta)%*%mu/quad) * (mu%*%t(beta)%*%Sigma + Sigma%*%beta%*%t(mu)) + as.numeric((t(beta)%*%mu/quad)^2)*Sigma%*%beta%*%t(beta)%*%Sigma
)
## as.numeric(1/2/sqrt(2*pi)/quad*dnorm(t(beta)%*%mu/sqrt(quad)))*aa
## moment.1%*%t(moment.1)
## as.numeric(-dnorm(t(beta)%*%mu/sqrt(2)/sqrt(quad))/sqrt(2*quad))*(as.numeric(t(beta)%*%mu/quad)*Sigma%*%beta-mu)
## 1/sqrt(2*pi)/2/quad*dnorm(t(beta)%*%mu/sqrt(quad))
as.numeric(1/sqrt(2*pi)/quad*(dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(quad))/sqrt(3) - dnorm(t(beta)%*%mu/sqrt(quad))/2)) * aa + as.numeric(1/sqrt(6*pi)/quad*dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(quad)))*(Sigma - 1/quad*Sigma%*%beta%*%t(beta)%*%Sigma)
var.u.1

a <- runif(p)
B <- 1e3
Bs <- round(seq(10,5e3,len=30))
by.B <- sapply(Bs, function(B) {
    w <- t(rmvnorm(B,mu,Sigma))
    u.1s <- apply(w,2,function(w.i)u.1(beta,w.i))
    observed <- as.numeric(t(a)%*%u.1s)
    ## try <- u.1(beta,mu/sqrt(2))
    ## hist(a%*%obs)
    ## abline(v=a%*%try,col=2)
    var(observed) - t(a)%*%var.u.1%*%a
})
plot(Bs,by.B)
abline(h=0)






## 10 beta\neq gamma

## 10a check o(1/sqrt(n)) convergence of hajek projection #1 using noise for beta.hat and gamma.hat
require(mvtnorm)
start <- Sys.time()
source('misc.R')
set.seed(1)
p <- 3
## n <- 1e3
ns <- round(seq(1e2,1e3,len=4e1))
## by.n <- sapply(ns, function(n) {
## terms <- replicate(2e1, {
## tryCatch({
beta <- runif(p)/4
gamma <- runif(p)/4
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
mu.1 <- runif(p)
mu.0 <- runif(p)
## beta.sampler <- lda.sampler.init(beta,Sigma)        
## mu.diff <- c(0,mu.0-mu.1)
Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
by.n <- sapply(ns, function(n) {
    cat('.')
    pairs <- replicate(1e2, {
        ## list2env(beta.sampler$sample(n),globalenv())
        ## data <- beta.sampler$sample(n)
        ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
        x.0 <- rmvnorm(n,mu.0,Sigma)
        x.1 <- rmvnorm(n,mu.1,Sigma)
        beta.hat <- beta+rnorm(p)/sqrt(n)#lda.coefs(x.0,x.1)
        gamma.hat <- gamma+rnorm(p)/sqrt(n)
        auc.hat.beta.hat <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat(x.0%*%gamma.hat,x.1%*%gamma.hat)
        diff.hat.coef.hat <- auc.hat.beta.hat-auc.hat.gamma.hat
        auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
        auc.gamma.hat <- auc.scores(gamma.hat,mu.1-mu.0,2*Sigma)
        diff.coef.hat <- auc.beta.hat-auc.gamma.hat
        auc.hajek.beta.hat <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta.hat,Sigma)
        auc.hajek.gamma.hat <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,gamma.hat,Sigma)
        diff.hajek.coef.hat <- auc.hajek.beta.hat - auc.hajek.gamma.hat
        diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef.hat
    })
})
Sys.time() - start
sds <- apply(by.n,2,sd)
lm0 <- lm(log(sds) ~ log(ns))
lm0


## and #2 using lda coefs. for reduced.len>0, get faster than O(1/n)
## convergence, as noted previously (allows hatted quadratic von mises
## term to be ignored in the beta==gamma case)
require(mvtnorm)
start <- Sys.time()
source('misc.R')
set.seed(1)
p <- 10
## n <- 1e3
ns <- round(seq(1e2,1e3,len=4e1))
## by.n <- sapply(ns, function(n) {
## terms <- replicate(2e1, {
## tryCatch({
reduced.len <- 1
beta <- runif(p)/4
beta <- c(beta[1:(p-reduced.len)],rep(0,reduced.len))
## gamma <- runif(p)/4
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
## mu.1 <- runif(p)
## mu.0 <- runif(p)
## beta.sampler <- lda.sampler.init(beta,Sigma)
mu.0 <- rep(0,p)
mu.1 <- mu.diff <- as.numeric(Sigma%*%beta)
## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
by.n <- sapply(ns, function(n) {
    cat('.')
    pairs <- replicate(1e2, {
        ## list2env(beta.sampler$sample(n),globalenv())
        ## data <- beta.sampler$sample(n)
        ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
        x.0 <- rmvnorm(n,mu.0,Sigma)
        x.1 <- rmvnorm(n,mu.1,Sigma)
        beta.hat <- lda.coefs(cbind(1,x.0),cbind(1,x.1))[-1]
        gamma.hat <- lda.coefs(cbind(1,x.0[,1:(p-reduced.len)]),cbind(1,x.1[,1:(p-reduced.len)]))[-1]
        gamma.hat <- c(gamma.hat,rep(0,length(beta.hat)-length(gamma.hat)))
        auc.hat.beta.hat <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat(x.0%*%gamma.hat,x.1%*%gamma.hat)
        diff.hat.coef.hat <- auc.hat.beta.hat-auc.hat.gamma.hat
        auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
        auc.gamma.hat <- auc.scores(gamma.hat,mu.1-mu.0,2*Sigma)
        diff.coef.hat <- auc.beta.hat-auc.gamma.hat
        auc.hajek.beta.hat <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta.hat,Sigma)
        auc.hajek.gamma.hat <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,gamma.hat,Sigma)
        diff.hajek.coef.hat <- auc.hajek.beta.hat - auc.hajek.gamma.hat
        diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef.hat
    })
})
Sys.time() - start
sds <- apply(by.n,2,sd)
lm0 <- lm(log(sds) ~ log(ns))
lm0



## 10b logistic model. auc(beta^t x) seems maximized where beta is the
## true logistic coefficients
source('misc.R')
p <- 4
n <- 1e2
beta <- runif(p+1)
delta <- runif(length(beta))
ts <- seq(-.1,1,len=40)*2
by.t <- sapply(ts, function(t) {
    aucs <- replicate(1e3, {
        x <- cbind(1,matrix(runif(n*p),ncol=p))
        risk <- plogis(x%*%beta)
        d <- rbinom(n,1,prob=risk)
        x0 <- x[d==0,]; x1 <- x[d==1,]
        auc.hat(x0%*%(beta+t*delta), x1%*%(beta+t*delta))
    })
    mean(aucs)
})       
plot(ts,by.t)    
abline(v=0)

## 10c logistic model.
source('misc.R')
p <- 4
n <- 1e3
beta <- runif(p+1)
diffs <- replicate(2e2, {
    x <- cbind(1,matrix(rnorm(n*p,mean=0),ncol=p))
    risk <- plogis(x%*%beta)
    d <- rbinom(n,1,prob=risk)
    ## x0 <- x[d==0,]; x1 <- x[d==1,]
    glm.full <- glm(d ~ x,family=binomial)
    glm.reduced <- glm(d ~ x[,-ncol(x)],family=binomial)
    x.0.full <- fitted(glm.full)[d==0]
    x.1.full <- fitted(glm.full)[d==1]
    x.0.reduced <- fitted(glm.reduced)[d==0]
    x.1.reduced <- fitted(glm.reduced)[d==1]
    auc.full <- auc.hat(x.0.full,x.1.full)
    auc.reduced <- auc.hat(x.0.reduced,x.1.reduced)
    auc.full - auc.reduced
})
hist(diffs)
qqnorm((diffs-mean(diffs))/sd(diffs));abline(0,1)


## 10c checking formulas
require(mvtnorm)
start <- Sys.time()
source('misc.R')
set.seed(1)
p <- 3
## n <- 1e3
ns <- round(seq(1e2,5e2,len=2e1))
## by.n <- sapply(ns, function(n) {
## terms <- replicate(2e1, {
## tryCatch({
beta <- runif(p)/4
gamma <- runif(p)/4
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
mu.1 <- runif(p)
mu.0 <- runif(p)
## beta.sampler <- lda.sampler.init(beta,Sigma)        
## mu.diff <- c(0,mu.0-mu.1)
## Sigma.diff <- 2*Sigma#cbind(0,rbind(0,Sigma))
u <- function(beta,w)pnorm(t(beta)%*%w/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))
u.1 <- function(beta,w) {
    qform <- as.numeric(t(beta)%*%Sigma%*%beta)
    as.numeric(-dnorm(t(beta)%*%w/sqrt(qform))/sqrt(qform))*(as.numeric(t(beta)%*%w/qform)*Sigma%*%beta-w)
}
## version of u.1 to take matrix of w's, each w in a column
u.1 <- function(beta,w) {
    qform <- as.numeric(t(beta)%*%Sigma%*%beta)
    t(as.numeric(-dnorm(t(beta)%*%w/sqrt(qform))/sqrt(qform)) *  t(kronecker(t(beta)%*%w/qform,Sigma%*%beta) - w))
}
## w <- matrix(runif(n*p),nrow=p)
## apply(w,2,function(w.i)u.1(beta,w.i)) - u.1a(beta,w)
auc.scores.hajek <- function(x.0,x.1,mu.0,mu.1,beta,Sigma) {
    x.0 <- t(x.0); x.1 <- t(x.1) 
    quad <- as.numeric(beta%*%Sigma%*%beta)
    theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
    hajek <- -mean(pnorm(t(beta)%*%(x.0-mu.1) / sqrt(quad)) - (1-theta)) + mean(pnorm(t(beta)%*%(x.1-mu.0) / sqrt(quad)) - theta)
}
## theta.bar <- function(x.0,x.1,beta) auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta,Sigma)
theta.bar <- function(beta,x.0,x.1) {
    x.0 <- t(x.0); x.1 <- t(x.1) 
    ## mean(apply(x.1,1,function(x)u(beta,x-mu.0))) +    mean(apply(x.0,1,function(x)u(beta,mu.1-x))) - 2*u(beta,(mu.1-mu.0)/sqrt(2))
    mean(u(beta,x.1-mu.0)) + mean(u(beta,-(x.0-mu.1))) - 2*u(beta,(mu.1-mu.0)/sqrt(2))
}
theta.bar.1 <- function(beta,x.0,x.1) rowMeans(u.1(beta,t(x.1)-mu.0)) + rowMeans(u.1(beta,-(t(x.0)-mu.1))) - 2*u.1(beta,(mu.1-mu.0)/sqrt(2))
by.n <- sapply(ns, function(n) {
    cat('.')
    pairs <- replicate(1e2, {
        ## list2env(beta.sampler$sample(n),globalenv())
        ## data <- beta.sampler$sample(n)
        ## beta.hat <- beta+rnorm(p+1)/sqrt(n)
        x.0 <- rmvnorm(n,mu.0,Sigma)
        x.1 <- rmvnorm(n,mu.1,Sigma)
        beta.hat <- beta+rnorm(p)/sqrt(n)#lda.coefs(x.0,x.1)
        gamma.hat <- gamma+rnorm(p)/sqrt(n)
        auc.hat.beta.hat <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat(x.0%*%gamma.hat,x.1%*%gamma.hat)
        diff.hat.coef.hat <- auc.hat.beta.hat-auc.hat.gamma.hat
        auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
        auc.gamma.hat <- auc.scores(gamma.hat,mu.1-mu.0,2*Sigma)
        diff.coef.hat <- auc.beta.hat-auc.gamma.hat
        ## auc.hajek.beta.hat <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta.hat,Sigma)
        ## auc.hajek.gamma.hat <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,gamma.hat,Sigma)
        auc.hajek.beta.hat <- theta.bar(beta.hat,x.0,x.1)
        auc.hajek.gamma.hat <- theta.bar(gamma.hat,x.0,x.1)
        diff.hajek.coef.hat <- auc.hajek.beta.hat - auc.hajek.gamma.hat
        ## diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef.hat
        diff <- auc.scores(beta,mu.1-mu.0,2*Sigma) - auc.scores(gamma,mu.1-mu.0,2*Sigma)
        ## ## part.1 <- diff.hat.coef.hat-diff.coef.hat
        ## ## part.1 <- diff.hajek.coef.hat
        ## part.1 <- theta.bar(beta,x.0,x.1) + t(beta.hat-beta)%*%theta.bar.1(beta,x.0,x.1) - theta.bar(gamma,x.0,x.1) - t(gamma.hat-gamma)%*%theta.bar.1(gamma,x.0,x.1) #1
        ## ## part.2 <- diff.coef.hat - diff
        ## part.2 <- (beta.hat-beta)%*%u.1(beta,(mu.1-mu.0)/sqrt(2)) -  (gamma.hat-gamma)%*%u.1(gamma,(mu.1-mu.0)/sqrt(2)) #1
        ## diff.hat.coef.hat - diff - (part.1 + part.2)
        auc.hajek.beta <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta,Sigma)
        auc.hajek.gamma <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,gamma,Sigma)
        diff.hajek.coef <- auc.hajek.beta - auc.hajek.gamma
        diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef
        ## ## A <- (beta.hat-beta)%*%(rowMeans(apply(x.1,1, function(x.1i)u.1(beta,x.1i-mu.0))) +  rowMeans(apply(x.0,1, function(x.0i)u.1(beta,mu.1-x.0i))) - u.1(beta,(mu.1-mu.0)/sqrt(2))) #1
        ## ## B <- (gamma.hat-gamma)%*%(rowMeans(apply(x.1,1, function(x.1i)u.1(gamma,x.1i-mu.0))) +  rowMeans(apply(x.0,1, function(x.0i)u.1(gamma,mu.1-x.0i))) - u.1(gamma,(mu.1-mu.0)/sqrt(2))) #1
        ## ## A <- (beta.hat-beta)%*%(rowMeans(u.1(beta,t(x.1)-mu.0)) + rowMeans(u.1(beta,-(t(x.0)-mu.1))) - u.1(beta,(mu.1-mu.0)/sqrt(2))) #1
        ## ## B <- (gamma.hat-gamma)%*%(rowMeans(u.1(gamma,t(x.1)-mu.0)) + rowMeans(u.1(gamma,-(t(x.0)-mu.1))) - u.1(gamma,(mu.1-mu.0)/sqrt(2))) #1
        ## ## approx <- diff.hajek.coef + A - B #1
        ## ## browser()
        ## C <- c(beta.hat-beta,gamma.hat-gamma)%*%rbind(rowMeans(u.1(beta,t(x.1)-mu.0)) - u.1(beta,(mu.1-mu.0)/sqrt(2)), -rowMeans(u.1(gamma,t(x.1)-mu.0))+ u.1(gamma,(mu.1-mu.0)/sqrt(2)))
        ## D <- c(beta.hat-beta,gamma.hat-gamma)%*%c(rowMeans(u.1(beta,-(t(x.0)-mu.1))), -rowMeans(u.1(gamma,-(t(x.0)-mu.1))))
        ## approx <- diff.hajek.coef + C + D
        ## obs <- diff.hat.coef.hat - diff
        ## obs-approx
    })
})
Sys.time() - start
sds <- apply(by.n,2,sd)
lm0 <- lm(log(sds) ~ log(ns))
lm0

## 10d  diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef is o(n^(-1/2))
require(mvtnorm)
start <- Sys.time()
source('misc.R')
set.seed(3)
p <- 6
## n <- 1e3
ns <- round(seq(1e2,5e2,len=2e1))
beta <- runif(p)/4
gamma <- runif(p)/4
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
mu.1 <- runif(p)
mu.0 <- runif(p)
by.n <- sapply(ns, function(n) {
    cat('.')
    pairs <- replicate(1e2, {
        x.0 <- rmvnorm(n,mu.0,Sigma)
        x.1 <- rmvnorm(n,mu.1,Sigma)
        beta.hat <- beta+rnorm(p)/sqrt(n)#lda.coefs(x.0,x.1)
        gamma.hat <- gamma+rnorm(p)/sqrt(n)
        auc.hat.beta.hat <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat(x.0%*%gamma.hat,x.1%*%gamma.hat)
        diff.hat.coef.hat <- auc.hat.beta.hat-auc.hat.gamma.hat
        auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
        auc.gamma.hat <- auc.scores(gamma.hat,mu.1-mu.0,2*Sigma)
        diff.coef.hat <- auc.beta.hat-auc.gamma.hat
        diff <- auc.scores(beta,mu.1-mu.0,2*Sigma) - auc.scores(gamma,mu.1-mu.0,2*Sigma)
        auc.hajek.beta <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta,Sigma)
        auc.hajek.gamma <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,gamma,Sigma)
        diff.hajek.coef <- auc.hajek.beta - auc.hajek.gamma
        diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef
    })
})
Sys.time() - start
sds <- apply(by.n,2,sd)
lm0 <- lm(log(sds) ~ log(ns))
lm0

## using estimated coefs in diff.hajek.coef, still O(1/n)
require(mvtnorm)
start <- Sys.time()
source('misc.R')
set.seed(3)
p <- 6
## n <- 1e3
ns <- round(seq(1e2,5e2,len=2e1))
beta <- runif(p)/4
gamma <- runif(p)/4
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
mu.1 <- runif(p)
mu.0 <- runif(p)
by.n <- sapply(ns, function(n) {
    cat('.')
    pairs <- replicate(1e2, {
        x.0 <- rmvnorm(n,mu.0,Sigma)
        x.1 <- rmvnorm(n,mu.1,Sigma)
        beta.hat <- beta+rnorm(p)/sqrt(n)#lda.coefs(x.0,x.1)
        gamma.hat <- gamma+rnorm(p)/sqrt(n)
        auc.hat.beta.hat <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat(x.0%*%gamma.hat,x.1%*%gamma.hat)
        diff.hat.coef.hat <- auc.hat.beta.hat-auc.hat.gamma.hat
        auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
        auc.gamma.hat <- auc.scores(gamma.hat,mu.1-mu.0,2*Sigma)
        diff.coef.hat <- auc.beta.hat-auc.gamma.hat
        diff <- auc.scores(beta,mu.1-mu.0,2*Sigma) - auc.scores(gamma,mu.1-mu.0,2*Sigma)
        auc.hajek.beta <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta.hat,Sigma)
        auc.hajek.gamma <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,gamma.hat,Sigma)
        diff.hajek.coef <- auc.hajek.beta - auc.hajek.gamma
        diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef
    })
})
Sys.time() - start
sds <- apply(by.n,2,sd)
lm0 <- lm(log(sds) ~ log(ns))
lm0

## with lda coef estimates, still O(1/n)
require(mvtnorm)
start <- Sys.time()
source('misc.R')
set.seed(3)
p <- 6
reduced.len <- 1
## n <- 1e3
ns <- round(seq(1e2,5e2,len=2e1))
beta <- runif(p)/4
gamma <- runif(p)/4
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
mu.1 <- runif(p)
mu.0 <- runif(p)
by.n <- sapply(ns, function(n) {
    cat('.')
    pairs <- replicate(3e2, {
        x.0 <- rmvnorm(n,mu.0,Sigma)
        x.1 <- rmvnorm(n,mu.1,Sigma)
        beta.hat <- lda.coefs(x.0,x.1)
        gamma.hat <- lda.coefs(cbind(1,x.0[,1:(p-reduced.len)]),cbind(1,x.1[,1:(p-reduced.len)]))[-1]
        gamma.hat <- c(gamma.hat,rep(0,length(beta.hat)-length(gamma.hat)))
        auc.hat.beta.hat <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat(x.0%*%gamma.hat,x.1%*%gamma.hat)
        diff.hat.coef.hat <- auc.hat.beta.hat-auc.hat.gamma.hat
        auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
        auc.gamma.hat <- auc.scores(gamma.hat,mu.1-mu.0,2*Sigma)
        diff.coef.hat <- auc.beta.hat-auc.gamma.hat
        diff <- auc.scores(beta,mu.1-mu.0,2*Sigma) - auc.scores(gamma,mu.1-mu.0,2*Sigma)
        auc.hajek.beta <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta.hat,Sigma)
        auc.hajek.gamma <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,gamma.hat,Sigma)
        diff.hajek.coef <- auc.hajek.beta - auc.hajek.gamma
        diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef
    })
})
Sys.time() - start
sds <- apply(by.n,2,sd)
lm0 <- lm(log(sds) ~ log(ns))
lm0

dd

## 10e variance of delong part 

delong.var <- function(x=NULL,y=NULL,xy=NULL,g=NULL) {
    if(is.null(x)||is.null(y)) {
        if(is.null(xy)||is.null(g)) stop()
        g <- factor(g,labels=0:1)
        x <- xy[g==0,]; y <- xy[g==1,]
    }
    m <- nrow(x); n <- nrow(y)
    theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
    V.10.1 <- sapply(x[,1], function(x.i)mean(x.i<y[,1]))
    V.10.2 <- sapply(x[,2], function(x.i)mean(x.i<y[,2]))
    V.01.1 <- sapply(y[,1], function(y.i)mean(x[,1]<y.i))
    V.01.2 <- sapply(y[,2], function(y.i)mean(x[,2]<y.i))
    S.10 <- cov(cbind(V.10.1,V.10.2))
    S.01 <- cov(cbind(V.01.1,V.01.2))
    S <- S.10/m + S.01/n
    contrast <- matrix(c(1,-1),ncol=1)
    ## browser()
    t(contrast)%*%S%*%contrast
    mean((V.01.1-V.01.2)^2)/n + mean((V.10.1-V.10.2)^2)/m
}
require(mvtnorm)
start <- Sys.time()
source('misc.R')
set.seed(3)
p <- 6
## n <- 1e3
ns <- round(seq(1e2,4e2,len=2e1))
beta <- runif(p)/4
gamma <- runif(p)/4
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
mu.1 <- runif(p)
mu.0 <- runif(p)
by.n <- sapply(ns, function(n) {
    cat('.')
    stats <- replicate(1e2, {
        x.0 <- rmvnorm(n,mu.0,Sigma)
        x.1 <- rmvnorm(n,mu.1,Sigma)
        beta.hat <- beta+rnorm(p)/sqrt(n)#lda.coefs(x.0,x.1)
        gamma.hat <- gamma+rnorm(p)/sqrt(n)
        auc.hat.beta.hat <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat(x.0%*%gamma.hat,x.1%*%gamma.hat)
        diff.hat.coef.hat <- auc.hat.beta.hat-auc.hat.gamma.hat
        auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
        auc.gamma.hat <- auc.scores(gamma.hat,mu.1-mu.0,2*Sigma)
        diff.coef.hat <- auc.beta.hat-auc.gamma.hat
        ## diff <- auc.scores(beta,mu.1-mu.0,2*Sigma) - auc.scores(gamma,mu.1-mu.0,2*Sigma)
        ## auc.hajek.beta <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta,Sigma)
        ## auc.hajek.gamma <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,gamma,Sigma)
        ## diff.hajek.coef <- auc.hajek.beta - auc.hajek.gamma
        stat <- diff.hat.coef.hat-diff.coef.hat
        ## x <- 
        ## y <- cbind(x.1%*%beta,x.1%*%gamma)
        ## var.hat <- delong.var(x=cbind(x.0%*%beta,x.0%*%gamma),y=cbind(x.1%*%beta,x.1%*%gamma))
        var.hat <- delong.var(x=cbind(x.0%*%beta.hat,x.0%*%gamma.hat),y=cbind(x.1%*%beta.hat,x.1%*%gamma.hat))
        c(stat=stat,var.hat=var.hat) 
    })
    c(var.obs=var(stats['stat',]),var.fla=mean(stats['var.hat',]))
})
plot(ns,by.n['var.obs',])
points(ns,by.n['var.fla',],col=2)
plot(ns,by.n['var.obs',]-by.n['var.fla',])
lm(log(abs(by.n['var.obs',]-by.n['var.fla',])) ~ log(ns))




## 10e-1 variance of delong part, using coefs estimated by
delong.var <- function(x=NULL,y=NULL,xy=NULL,g=NULL) {
    if(is.null(x)||is.null(y)) {
        if(is.null(xy)||is.null(g)) stop()
        g <- factor(g,labels=0:1)
        x <- xy[g==0,]; y <- xy[g==1,]
    }
    m <- nrow(x); n <- nrow(y)
    theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
    V.10.1 <- sapply(x[,1], function(x.i)mean(x.i<y[,1]))
    V.10.2 <- sapply(x[,2], function(x.i)mean(x.i<y[,2]))
    V.01.1 <- sapply(y[,1], function(y.i)mean(x[,1]<y.i))
    V.01.2 <- sapply(y[,2], function(y.i)mean(x[,2]<y.i))
    S.10 <- cov(cbind(V.10.1,V.10.2))
    S.01 <- cov(cbind(V.01.1,V.01.2))
    S <- S.10/m + S.01/n
    contrast <- matrix(c(1,-1),ncol=1)
    ## browser()
    t(contrast)%*%S%*%contrast
    mean((V.01.1-V.01.2)^2)/n + mean((V.10.1-V.10.2)^2)/m
}
require(mvtnorm)
start <- Sys.time()
source('misc.R')
set.seed(1)
p <- 4
## n <- 1e3
ns <- round(seq(1e2,7e2,len=2e1))
beta <- runif(p)/4
## gamma <- runif(p)/4
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
mu.1 <- runif(p)
mu.0 <- runif(p)
by.n <- sapply(ns, function(n) {
    cat('.')
    stats <- replicate(1e2, {
        x.0 <- rmvnorm(n,mu.0,Sigma)
        x.1 <- rmvnorm(n,mu.1,Sigma)
        ## beta.hat <- beta+rnorm(p)/sqrt(n)#lda.coefs(x.0,x.1)
        ## gamma.hat <- gamma+rnorm(p)/sqrt(n)
        beta.hat <- lda.coefs(x.0,x.1)
        gamma.hat <- lda.coefs(x.0[,-p],x.1[,-p])
        gamma.hat <- c(gamma.hat,0)
        auc.hat.beta.hat <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat)
        auc.hat.gamma.hat <- auc.hat(x.0%*%gamma.hat,x.1%*%gamma.hat)
        diff.hat.coef.hat <- auc.hat.beta.hat-auc.hat.gamma.hat
        auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
        auc.gamma.hat <- auc.scores(gamma.hat,mu.1-mu.0,2*Sigma)
        diff.coef.hat <- auc.beta.hat-auc.gamma.hat
        ## diff <- auc.scores(beta,mu.1-mu.0,2*Sigma) - auc.scores(gamma,mu.1-mu.0,2*Sigma)
        ## auc.hajek.beta <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta,Sigma)
        ## auc.hajek.gamma <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,gamma,Sigma)
        ## diff.hajek.coef <- auc.hajek.beta - auc.hajek.gamma
        stat <- diff.hat.coef.hat-diff.coef.hat
        ## x <- 
        ## y <- cbind(x.1%*%beta,x.1%*%gamma)
        ## var.hat <- delong.var(x=cbind(x.0%*%beta,x.0%*%gamma),y=cbind(x.1%*%beta,x.1%*%gamma))
        var.hat <- delong.var(x=cbind(x.0%*%beta.hat,x.0%*%gamma.hat),y=cbind(x.1%*%beta.hat,x.1%*%gamma.hat))
        c(stat=stat,var.hat=var.hat) 
    })
    c(var.obs=var(stats['stat',]),var.fla=mean(stats['var.hat',]))
})
plot(ns,by.n['var.obs',])
points(ns,by.n['var.fla',],col=2)
plot(ns,by.n['var.obs',]-by.n['var.fla',])
lm(log(abs(by.n['var.obs',]-by.n['var.fla',])) ~ log(ns))



## 10f simplify delong var formula instead of using contrasts and full
## covariance matrix. looks close but a little off for smaller
## values. probably an m vs m-1 or n vs n-1 somewhere. TODO theta.hat^2 terms
delong.var.old <- function(x=NULL,y=NULL,xy=NULL,g=NULL) {
    if(is.null(x)||is.null(y)) {
        if(is.null(xy)||is.null(g)) stop()
        g <- factor(g,labels=0:1)
        x <- xy[g==0,]; y <- xy[g==1,]
    }
    m <- nrow(x); n <- nrow(y)
    theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
    V.10.1 <- sapply(x[,1], function(x.i)mean(x.i<y[,1]))
    V.10.2 <- sapply(x[,2], function(x.i)mean(x.i<y[,2]))
    V.01.1 <- sapply(y[,1], function(y.i)mean(x[,1]<y.i))
    V.01.2 <- sapply(y[,2], function(y.i)mean(x[,2]<y.i))
    S.10 <- cov(cbind(V.10.1,V.10.2))
    S.01 <- cov(cbind(V.01.1,V.01.2))
    S <- S.10/m + S.01/n
    contrast <- matrix(c(1,-1),ncol=1)
    t(contrast)%*%S%*%contrast
}
delong.var.new <- function(x=NULL,y=NULL,xy=NULL,g=NULL) {
    browser()
    if(is.null(x)||is.null(y)) {
        if(is.null(xy)||is.null(g)) stop()
        g <- factor(g,labels=0:1)
        x <- xy[g==0,]; y <- xy[g==1,]
    }
    m <- nrow(x); n <- nrow(y)
    theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
    V.10.1 <- sapply(x[,1], function(x.i)mean(x.i<y[,1]))
    V.10.2 <- sapply(x[,2], function(x.i)mean(x.i<y[,2]))
    V.01.1 <- sapply(y[,1], function(y.i)mean(x[,1]<y.i))
    V.01.2 <- sapply(y[,2], function(y.i)mean(x[,2]<y.i))
    var(V.10.1-V.10.2)/m + var(V.01.1-V.01.2)/n
    ## S.10 <- cov(cbind(V.10.1,V.10.2))
    ## S.01 <- cov(cbind(V.01.1,V.01.2))
    ## 1/m*(S.10[1,1]+S.10[2,2]-2*S.10[1,2]) + 1/n*(S.01[1,1]+S.01[2,2]-2*S.01[1,2])
}
delong.var.new <- function(x=NULL,y=NULL,xy=NULL,g=NULL) {
    ## browser()
    if(is.null(x)||is.null(y)) {
        if(is.null(xy)||is.null(g)) stop()
        g <- factor(g,labels=0:1)
        x <- xy[g==0,]; y <- xy[g==1,]
    }
    ## m <- nrow(x); n <- nrow(y)
    ## theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
    ## V.10.1 <- sapply(x[,1], function(x.i)mean(x.i<y[,1]))
    ## V.10.2 <- sapply(x[,2], function(x.i)mean(x.i<y[,2]))
    ## V.01.1 <- sapply(y[,1], function(y.i)mean(x[,1]<y.i))
    ## V.01.2 <- sapply(y[,2], function(y.i)mean(x[,2]<y.i))
    F.hat.1 <- ecdf(x[,1])
    F.hat.2 <- ecdf(x[,2])
    G.hat.1 <- ecdf(y[,1])
    G.hat.2 <- ecdf(y[,2])
    ## V.01.1 - F.hat.1(y[,1])
    ## V.10.1 - (1-G.hat.1(x[,1]))
    ## V.01.2 - F.hat.2(y[,2])
    ## V.10.2 - (1-G.hat.2(x[,2]))
    ## diff.hajek.hat <- mean(G.hat.2(x[,2]) - G.hat.1(x[,1])) + mean(F.hat.1(y[,1]) - F.hat.2(y[,2])) + 2*diff(theta.hats)
    var(G.hat.2(x[,2]) - G.hat.1(x[,1]))/m + var(F.hat.1(y[,1]) - F.hat.2(y[,2]))/n
    ## var(V.10.1-V.10.2)/m + var(V.01.1-V.01.2)/n
}
## set.seed(1)
m <- 10
n <- 20
x.0 <- matrix(runif(2*m),ncol=2)
x.1 <- matrix(runif(2*n),ncol=2)
delong.var.old(x=x.0,y=x.1)
delong.var.new(x=x.0,y=x.1)


# 10g refactoring hajek calculation. auc.hajek is now a function
# of F and G (van der vaart notation). #1 Using the true F and G in
# the normal case, replace auc.scores.hajek. #2 Using the estimated F
# and G, replace delong.var().

##1 
auc.hajek <- function(F,G,x,y,auc,terms.only=TRUE,IID=FALSE) {
    ## -G(x) + F(y) + 1-2*auc
    terms <- list(control = -(G(x)-(1-auc)), case = F(y)-auc)
    if(terms.only) return(terms)
    if(IID) {
        ## stopifnot(terms.only)
        n.0 <- length(x); n.1 <- length(y)
        n <- n.0+n.1
        ## g <- c(rep(0,n.0),rep(1,n.1))
        return(c(terms$control*n/n.0,terms$case*n/n.1))
        }
    return(sum(sapply(terms,mean)))
}
## m <- 10
## n <- 20
## x <- matrix(runif(2*m),ncol=2)
## y <- matrix(runif(2*n),ncol=2)
## F.hat.1 <- ecdf(x[,1])
## F.hat.2 <- ecdf(x[,2])
## G.hat.1 <- ecdf(y[,1])
## G.hat.2 <- ecdf(y[,2])
## theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
## hajek.diff.terms <- mapply('-',auc.hajek(F.hat.1,G.hat.1,x[,1],y[,1],theta.hats[1]),auc.hajek(F.hat.2,G.hat.2,x[,2],y[,2],theta.hats[2]))
## sum(sapply(hajek.diff.terms,mean))
source('misc.R')
p <- 3
n <- 10
mu.0 <- runif(p)
mu.1 <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%Sigma
x.0 <- matrix(runif(n*p),ncol=p)
x.1 <- matrix(runif(n*p),ncol=p)
beta <- runif(p)
quad <- as.numeric(beta%*%Sigma%*%beta)
F <- function(u)pnorm((u-as.numeric(beta%*%mu.0))/sqrt(quad)) #ecdf of beta^Tx.0
G <- function(u)pnorm((u-as.numeric(beta%*%mu.1))/sqrt(quad)) #ecdf of beta^Tx.1
theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
auc.hajek(F,G,x.0%*%beta,x.1%*%beta,auc=theta,FALSE)
auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta,Sigma)
auc.hajek.lda <- function(x.0,x.1,mu.0,mu.1,beta,Sigma,terms.only=FALSE) {
    quad <- as.numeric(beta%*%Sigma%*%beta)
    F <- function(u)pnorm((u-as.numeric(beta%*%mu.0))/sqrt(quad)) #ecdf of beta^Tx.0
    G <- function(u)pnorm((u-as.numeric(beta%*%mu.1))/sqrt(quad)) #ecdf of beta^Tx.1
    theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
    auc.hajek(F,G,x.0%*%beta,x.1%*%beta,auc=theta,terms.only)
} 


#2
delong.var.old <- function(x=NULL,y=NULL,xy=NULL,g=NULL) {
    if(is.null(x)||is.null(y)) {
        if(is.null(xy)||is.null(g)) stop()
        g <- factor(g,labels=0:1)
        x <- xy[g==0,]; y <- xy[g==1,]
    }
    m <- nrow(x); n <- nrow(y)
    theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
    V.10.1 <- sapply(x[,1], function(x.i)mean(x.i<y[,1]))
    V.10.2 <- sapply(x[,2], function(x.i)mean(x.i<y[,2]))
    V.01.1 <- sapply(y[,1], function(y.i)mean(x[,1]<y.i))
    V.01.2 <- sapply(y[,2], function(y.i)mean(x[,2]<y.i))
    S.10 <- cov(cbind(V.10.1,V.10.2))
    S.01 <- cov(cbind(V.01.1,V.01.2))
    S <- S.10/m + S.01/n
    contrast <- matrix(c(1,-1),ncol=1)
    t(contrast)%*%S%*%contrast
}
delong.var.new <- function(x=NULL,y=NULL,xy=NULL,g=NULL) {
    if(is.null(x)||is.null(y)) {
        if(is.null(xy)||is.null(g)) stop()
        g <- factor(g,labels=0:1)
        x <- xy[g==0,]; y <- xy[g==1,]
    }
    m <- nrow(x); n <- nrow(y)
    theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
    F.hat.1 <- ecdf(x[,1])
    F.hat.2 <- ecdf(x[,2])
    G.hat.1 <- ecdf(y[,1])
    G.hat.2 <- ecdf(y[,2])
    hajek.diff.terms <- mapply('-',auc.hajek(F.hat.1,G.hat.1,x[,1],y[,1],0),auc.hajek(F.hat.2,G.hat.2,x[,2],y[,2],0))
    sum(sapply(hajek.diff.terms,var) / c(m,n))
}
m <- 10
n <- 20
x <- matrix(runif(2*m),ncol=2)
y <- matrix(runif(2*n),ncol=2)
## F.hat.1 <- ecdf(x[,1])
## F.hat.2 <- ecdf(x[,2])
## G.hat.1 <- ecdf(y[,1])
## G.hat.2 <- ecdf(y[,2])
## theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
## hajek.diff.terms <- mapply('-',auc.hajek(F.hat.1,G.hat.1,x[,1],y[,1],0),auc.hajek(F.hat.2,G.hat.2,x[,2],y[,2],0))
## sum(sapply(hajek.diff.terms,var) / c(m,n))
delong.var.old(x,y)
delong.var.new(x,y)



m <- nrow(x); n <- nrow(y)
theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
F.hat.1 <- ecdf(x[,1])
F.hat.2 <- ecdf(x[,2])
G.hat.1 <- ecdf(y[,1])
G.hat.2 <- ecdf(y[,2])
## Below are the IID terms for the delong part, control in the first
## list item and case in the second. The control and case terms for
## the derivative part will be added to these, before taking the
## sample variance as in the last line.
hajek.diff.terms <- mapply('-',auc.hajek(F.hat.1,G.hat.1,x[,1],y[,1],0),auc.hajek(F.hat.2,G.hat.2,x[,2],y[,2],0))
sum(sapply(hajek.diff.terms,var) / c(m,n))

source('misc.r')
n <- 1e2
n.0 <- sample(n,1)
n.1 <- sample(n,1)
x.a0 <- rnorm(n.0)
x.a1 <- rnorm(n.1)
x.b0 <- rnorm(n.0)
x.b1 <- rnorm(n.1)
hajek.a <- auc.hajek(x=x.a0,y=x.a1)
hajek.b <- auc.hajek(x=x.b0,y=x.b1)
var(hajek.a$control-hajek.b$control)/n.0 + var(hajek.a$case-hajek.b$case)/n.1
delong.var.new(cbind(x.a0,x.b0),cbind(x.a1,x.b1))

dd

## checking formulas
## x <- x.0; y <- x.1
## if(is.null(x)||is.null(y)) {
##     if(is.null(xy)||is.null(g)) stop()
##     g <- factor(g,labels=0:1)
##     x <- xy[g==0,]; y <- xy[g==1,]
## }
## m <- nrow(x); n <- nrow(y)
## theta.hat <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
## V.10.1 <- sapply(x[,1], function(x.i)mean(x.i<y[,1]))
## V.10.2 <- sapply(x[,2], function(x.i)mean(x.i<y[,2]))
## V.01.1 <- sapply(y[,1], function(y.i)mean(x[,1]<y.i))
## V.01.2 <- sapply(y[,2], function(y.i)mean(x[,2]<y.i))
## S.10 <- cov(cbind(V.10.1,V.10.2))
## S.01 <- cov(cbind(V.01.1,V.01.2))
## 1/m*(S.10[1,1]+S.10[2,2]-2*S.10[1,2]) + 1/n*(S.01[1,1]+S.01[2,2]-2*S.01[1,2])
## 1/m/(m-1)*( sum((V.10.1-theta.hat[1])^2) + sum((V.10.2-theta.hat[2])^2) - 2*sum((V.10.1-theta.hat[1])*(V.10.2-theta.hat[2])) ) +
##     1/n/(n-1)*( sum((V.01.1-theta.hat[1])^2) + sum((V.01.2-theta.hat[2])^2) - 2*sum((V.01.1-theta.hat[1])*(V.01.2-theta.hat[2])) )
## var(V.10.1-V.10.2)/m + var(V.01.1-V.01.2)/n


## 11 logistic part

## Flaw in many of these is that the replicate reps should be increasing as well as n.
## Also: I've been using fitted() when I should have been using predict(). The
## difference is that fitted applies the link function.



## 11a misspecified logistic model. sandwich se doesnt seem too
## different from regular se. both seem very noisy. #1 non random x
## and #2 random x.
set.seed(2)
require(sandwich)
var.logistic <- function(glm0) {
    X <- cbind(1,glm0$model[[2]])
    Y <- glm0$model[[1]]
    fits <- fitted(glm0)
    W.hat <- diag(fits*(1-fits))
    W.tilde <- diag((Y-fits)^2)
    ## cov.beta <-  solve(t(X)%*%W%*%X)
    ## diag(cov.beta)
    ## browser()
    bread  <-  solve(t(X)%*%W.hat%*%X)
    meat <- t(X)%*%W.tilde%*%X
    cov.beta <- bread%*%meat%*%bread
    diag(cov.beta)
}
n <- 5e2
## set.seed(1)
p <- 4
## beta0 <- runif(1)
ns <- round(seq(1e2,1e4,len=70))
by.n <- sapply(ns, function(n) {
beta <- runif(p+1)*2
    ## x <- matrix(rnorm(n*p),ncol=p)  #1
    hats <- replicate(1e2, {
        x <- matrix(rnorm(n*p),ncol=p)
        ## risk <- plogis(x%*%beta) #2
        risk <- plogis(cbind(1,x)%*%beta)
        g <- rbinom(n,1,risk)
        ## coef(glm(g~x-1))-beta
        ## sqrt(diag(vcov(glm(g~x-1))))
        glm0 <- glm(g~x[,1],family=binomial)
        p0 <- length(coef(glm0))
        coef  <- unname(coef(glm0)[p0])
        ## browser()
        se <- sqrt(vcov(glm0)[p0,p0])
        ## coef.se <- sqrt(var.logistic(glm0))[p]
        se.robust <- unname(sqrt(diag(sandwich(glm0)))[p0])
        c(coef=coef,se=se,se.robust=se.robust) 
    })
    ## hist(hats['se',])
    ## abline(v=sd(hats['coef',]),col=2)
    se <- mean(hats['se',])
    se.robust <- mean(hats['se.robust',])
    obs <- sd(hats['coef',])
    c(obs=obs,se=se,se.robust=se.robust)
})
## plot(ns,by.n[1,])
## points(ns,by.n[2,],col=2)
## diffs <- abs(by.n[1,]-by.n[2,])
## ## hist(hats['var',])
## ##  abline(v=var(hats['coef',]),col=2)
## ##  abline(v=var.true,col=3)
## plot(ns,sqrt(ns)*diffs,ylim=range(sqrt(ns)*diffs))
diffs <- by.n['obs',]-by.n['se',]
diffs.robust <- by.n['obs',]-by.n['se.robust',]
plot(ns,sqrt(ns)*diffs,ylim=range(c(sqrt(ns)*diffs,sqrt(ns)*diffs.robust)))
points(ns,sqrt(ns)*diffs.robust,col=2)
abline(h=0)
lm(log(abs(diffs))~log(ns))
lm(log(abs(diffs.robust))~log(ns))

n <- 5e2
## set.seed(1)
p <- 4
## beta0 <- runif(1)
ns <- round(seq(1e2,1e3,len=20))
by.n <- sapply(ns, function(n) {
    beta <- runif(p+1)*2
    x <- matrix(rnorm(n*p),ncol=p)
    hats <- replicate(1e2, {
        ## risk <- plogis(x%*%beta)
        risk <- plogis(cbind(1,x)%*%beta)
        g <- rbinom(n,1,risk)
        ## coef(glm(g~x-1))-beta
        ## sqrt(diag(vcov(glm(g~x-1))))
        glm0 <- glm(g~x[,1],family=binomial)
        p0 <- length(coef(glm0))
        coef  <- coef(glm0)[p0]
        ## browser()
        var.beta <- vcov(glm0)[p0,p0]
        ## coef.se <- sqrt(var.logistic(glm0))[p]
        var.beta.robust <- diag(sandwich(glm0))[p0]
        c(coef=unname(coef),var.beta=var.beta,var.beta.robust=unname(var.beta.robust))
    })
    ## hist(hats['se',])
    ## abline(v=sd(hats['coef',]),col=2)
    var.beta <- mean(hats['var.beta',])
    var.beta.robust <- mean(hats['var.beta.robust',])
    obs <- var(hats['coef',])
    c(obs=obs,var.beta=var.beta,var.beta.robust=var.beta.robust)
})
plot(ns,ns*by.n['obs',] - ns*by.n['var.beta.robust',])

dd

## dbg. 1. sanity check with lm insead of glm. based on histogram,
## seems reported vars's are centered quite far from observed var/true
## var, as with glm. difference goes to 0 at a O(1/n) rate, but so
## does the variance. 2. Is var.hat(beta) as reported by lm converge
## to the true var(beta) faster than 1/n, so that the variance
## estimate of sqrt(n)*beta converges to the true variance of
## sqrt(n)beta? Or is it 1/n, and n * the variance estimate is an
## approximately unbiased estimate of n * the true variance?
## 3. shouldnt be looking at convergence rates, should be looking at
## bias. normalized by sqrt(n) the true and estimated sd do not seem
## to come together, but the estimate is still unbiased around the
## true.
set.seed(4)
n <- 1e2
## p <- 1
beta <- runif(2)
ns <- round(seq(1e2,5e3,len=40))
by.n <- sapply(ns, function(n) {
    ## hats <- replicate(1e2, {
    x <- rnorm(n)
    sigma <- .1
    var.true <-  sigma^2/sum((x-mean(x))^2)
    hats <- replicate(1e2, {
        y <- beta[1] + x*beta[2] + rnorm(n,sd=sigma)
        lm0 <- lm(y~x)
        coef  <- unname(coef(lm0)[2])
        coef.var <- unname(diag(vcov(lm0))[2])
        ## coef.se <- unname(summary(lm0)$coef[p,2])
        c(coef=coef,var=coef.var)
    })
    c(obs=var(hats['coef',]),fla=var.true)
})
op <- par(mfrow=c(1,2))
plot(ns,ns*by.n[1,])
points(ns,ns*by.n[2,],col=2)
diffs <- abs(by.n[1,]-by.n[2,])
## hist(hats['var',])
##  abline(v=var(hats['coef',]),col=2)
##  abline(v=var.true,col=3)
plot(ns,ns*diffs,ylim=range(ns*diffs))
par(op)
lm(log(diffs)~log(ns))

diffs <- by.n[1,]-by.n[2,]
plot(ns,ns*diffs)
abline(h=0)




## as above but with one fixed n
library(sandwich)
n <- 1e4
## set.seed(1)
p <- 4
beta <- runif(p+1)*2
## x <- matrix(rnorm(n*p),ncol=p) #1
sim <- replicate(1e3, {
    x <- matrix(rnorm(n*p),ncol=p) #2
    risk <- plogis(cbind(1,x)%*%beta)
    y <- rbinom(n,1,risk)
    glm0 <- glm(y~x[,1],family=binomial)
    p0 <- 2
    coef  <- unname(coef(glm0)[p0])
    se <- sqrt(vcov(glm0)[p0,p0])
    se.robust <- sqrt(diag(sandwich(glm0)))[p0]
    c(coef=coef,se=se,se.robust=unname(se.robust) )
})
op <- par(mfrow=c(1,2))
observed.se <- sd(sim['coef',])
hist(sim['se',]-observed.se)
hist(sim['se.robust',]-observed.se)
par(op)
hist(sqrt(n)*(sim['se',]-sim['se.robust',])) # se and robust se are not quite

dd

## ## 11b IID representation of  diff.hat.coef.hat-diff.coef
## require(mvtnorm)
## ## start <- Sys.time()
## source('misc.R')
## set.seed(1)
## p <- 3
## ## n <- 1e3
## ns <- round(seq(1e2,5e2,len=2e1))
## beta <- runif(p)/4
## gamma <- runif(p)/4
## Sigma <- matrix(rnorm(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
## mu.1 <- runif(p)
## mu.0 <- runif(p)
## n.0 <- 1e2; n.1 <- 1e2
## n <- n.0+n.1
## u <- function(beta,w)pnorm(t(beta)%*%w/as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))
## u.1 <- function(beta,w) {
##     qform <- as.numeric(t(beta)%*%Sigma%*%beta)
##     as.numeric(-dnorm(t(beta)%*%w/sqrt(qform))/sqrt(qform))*(as.numeric(t(beta)%*%w/qform)*Sigma%*%beta-w)
## }
## ## version of u.1 to take matrix of w's, each w in a column
## u.1 <- function(beta,w) {
##     qform <- as.numeric(t(beta)%*%Sigma%*%beta)
##     t(as.numeric(-dnorm(t(beta)%*%w/sqrt(qform))/sqrt(qform)) *  t(kronecker(t(beta)%*%w/qform,Sigma%*%beta) - w))
## }
## by.n <- sapply(ns, function(n) {
##     cat('.')
##     pairs <- replicate(1e2, {
##         x.0 <- rmvnorm(n,mu.0,Sigma)
##         x.1 <- rmvnorm(n,mu.1,Sigma)
##         beta.hat <- beta+rnorm(p)/sqrt(n)#lda.coefs(x.0,x.1)
##         gamma.hat <- gamma+rnorm(p)/sqrt(n)
##         auc.hat.beta.hat <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat)
##         auc.hat.gamma.hat <- auc.hat(x.0%*%gamma.hat,x.1%*%gamma.hat)
##         diff.hat.coef.hat <- auc.hat.beta.hat-auc.hat.gamma.hat
##         diff.coef <- auc.scores(beta,mu.1-mu.0,2*Sigma) - auc.scores(gamma,mu.1-mu.0,2*Sigma)
##         ## quad.beta <- as.numeric(t(beta)%*%Sigma%*%beta)
##         ## quad.gamma <- as.numeric(t(gamma)%*%Sigma%*%gamma)
##         ## infl.beta <- t(X)%*%(
##         approx <- mean(u(beta,t(x.1)-mu.0) - u(gamma,t(x.1)-mu.0)) + 
##             mean(u(gamma,t(x.0)-mu.1) - u(beta,t(x.0)-mu.1)) +
##             t(auc.scores.deriv(beta,mu.1-mu.0,2*Sigma))%*%(beta.hat-beta)/n -
##             t(auc.scores.deriv(gamma,mu.1-mu.0,2*Sigma))%*%(gamma.hat-gamma)/n        
##         auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
##         auc.gamma.hat <- auc.scores(gamma.hat,mu.1-mu.0,2*Sigma)
##         diff.coef.hat <- auc.beta.hat-auc.gamma.hat
##         ## diff <- auc.scores(beta,mu.1-mu.0,2*Sigma) - auc.scores(gamma,mu.1-mu.0,2*Sigma)
##         auc.hajek.beta <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta,Sigma)
##         auc.hajek.gamma <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,gamma,Sigma)
##         diff.hajek.coef <- auc.hajek.beta - auc.hajek.gamma
##         diff.hat.coef.hat-diff.coef - approx
##         diff.hat.coef.hat-diff.coef.hat 
##         mean(u(beta,t(x.1)-mu.0) - u(gamma,t(x.1)-mu.0)) + 
##             mean(u(gamma,t(x.0)-mu.1) - u(beta,t(x.0)-mu.1))
##         diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef
##         ## x.0 <- rmvnorm(n,mu.0,Sigma)
##         ## x.1 <- rmvnorm(n,mu.1,Sigma)
##         ## beta.hat <- beta+rnorm(p)/sqrt(n)#lda.coefs(x.0,x.1)
##         ## gamma.hat <- gamma+rnorm(p)/sqrt(n)
##         ## auc.hat.beta.hat <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat)
##         ## auc.hat.gamma.hat <- auc.hat(x.0%*%gamma.hat,x.1%*%gamma.hat)
##         ## diff.hat.coef.hat <- auc.hat.beta.hat-auc.hat.gamma.hat
##         ## auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
##         ## auc.gamma.hat <- auc.scores(gamma.hat,mu.1-mu.0,2*Sigma)
##         ## diff.coef.hat <- auc.beta.hat-auc.gamma.hat
##         ## diff <- auc.scores(beta,mu.1-mu.0,2*Sigma) - auc.scores(gamma,mu.1-mu.0,2*Sigma)
##         ## auc.hajek.beta <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta,Sigma)
##         ## auc.hajek.gamma <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,gamma,Sigma)
##         ## diff.hajek.coef <- auc.hajek.beta - auc.hajek.gamma
##         ## diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef
##     })
## })
## sds <- apply(by.n,2,sd)
## lm0 <- lm(log(sds) ~ log(ns))
## lm0


## 11b logistic regr influence function. #1 at true beta, approx error
## is O(1/n). Approx error is O(1/sqrt(n)) using hatted coefs, though
## (using random or non-random x).
p <- 4
n <- 3e2
beta <- runif(p+1)
ns <- round(seq(1e2,1e3,len=20))
by.n <- sapply(ns, function(n) {
    diffs <- replicate(1e2, {
        x <- cbind(1,matrix(runif(n*p),ncol=p))
        risk <- as.numeric(plogis(x%*%beta))
        y <- rbinom(n,1,risk)
        glm0 <- glm(y~x-1,family=binomial)
        fits <- fitted(glm0)
        beta.hat <- coef(glm0)
        W <- diag(risk*(1-risk))  #1
        W <- diag(fits*(1-fits)) #2
        l.prime.prime <- -t(x)%*%W%*%x 
        l.prime <- t(x)%*%(y-risk) #1
        l.prime <- t(x)%*%(y-fits) #2
        ## W.hat <- diag(fits*(1-fits))
        ## l.prime.prime <- -t(x)%*%W.hat%*%x
        ## l.prime <- t(x)%*%(y-fits)
        infl <- -solve(l.prime.prime)%*%l.prime
        beta.hat-beta - infl
    })
    apply(diffs,1,sd)
})
coef(lm(log(colSums(by.n))~log(ns)))

## rates for sds. sqrt(n)*infl.hat has different limit that
## sqrt(n)(beta.hat-beta)
p <- 4
n <- 3e2
beta <- runif(p+1)
ns <- round(seq(1e2,3e3,len=20))
by.n <- sapply(ns, function(n) {
    diffs <- replicate(1e2, {
        x <- cbind(1,matrix(runif(n*p),ncol=p))
        risk <- as.numeric(plogis(x%*%beta))
        y <- rbinom(n,1,risk)
        glm0 <- glm(y~x-1,family=binomial)
        fits <- fitted(glm0)
        beta.hat <- coef(glm0)
        W <- diag(risk*(1-risk))  #1
        l.prime.prime <- -t(x)%*%W%*%x 
        l.prime <- t(x)%*%(y-risk) #1
        W.hat <- diag(fits*(1-fits)) #2
        l.prime.prime.hat <- -t(x)%*%W.hat%*%x 
        l.prime.hat <- t(x)%*%(y-fits) #2
        infl <- -solve(l.prime.prime)%*%l.prime
        infl.hat <- -solve(l.prime.prime.hat)%*%l.prime.hat
        p0 <- 3
        c(obs=unname(beta.hat-beta)[p0],infl=infl[p0],infl.hat=infl.hat[p0])
    })
    apply(diffs,1,sd)
})
coef(lm(log(abs(by.n['obs',]-by.n['infl',]))~log(ns)))
coef(lm(log(abs(by.n['obs',]-by.n['infl.hat',]))~log(ns)))
plot(sqrt(ns)*(by.n['obs',]-by.n['infl.hat',]))

## V.sandwich^{-1/2}*beta.hat is nearly the identity. clearest from
## the heatmap. looking at the matrix itself it doesn't look very much
## like hte identity but then in #2 even in the correclty psecified
## model with the usual se's it doesn't look. asymptotics must be
## slow.
require(sandwich)
p <- 4
n <- 5e3
beta <- runif(p+1)
## ns <- round(seq(1e2,3e3,len=20))
## by.n <- sapply(ns, function(n) {
z.stats <- replicate(1e2, {
    x <- cbind(1,matrix(runif(n*p),ncol=p))
    risk <- as.numeric(plogis(x%*%beta))
    y <- rbinom(n,1,risk)
    glm0 <- glm(y~x[,2:3]-1,family=binomial)
    beta.hat <- coef(glm0)
    cov.beta.robust <- sandwich(glm0)
    solve(expm::sqrtm(cov.beta.robust))%*%beta.hat
    ## cov.beta <- vcov(glm0) #2
    ## solve(expm::sqrtm(cov.beta))%*%beta.hat #2
})
z.stats <- z.stats[,,,drop=TRUE]
abs(cov(t(z.stats)))
heatmap(abs(cov(t(z.stats))))







## 12 probit model

## 12a checking formulas
require(mvtnorm)
n <- 1e4
Sigma <- matrix(runif(4),nrow=2)
Sigma <- Sigma%*%t(Sigma)
obs <- replicate(1e2, {
    mu <- runif(2)*4
    xy <- rmvnorm(n,mu,Sigma)
    x <- xy[,1]; y <- xy[,2]
    proj <- mu[2] + (x-mu[1])*Sigma[1,2]/Sigma[1,1]
    ## var(proj)
    ## var(y-proj)
    var(x+proj)
})
## try <- Sigma[1,2]^2/Sigma[1,1]
try <- Sigma[2,2]-Sigma[1,2]^2/Sigma[1,1]
try <- Sigma[1,1]+2*Sigma[1,2]+Sigma[1,2]^2/Sigma[1,1]
## try <- Sigma[1,1]*(1+Sigma[1,2]/Sigma[2,2])^2
hist(obs)
abline(v=try,col=2)
abline(v=mean(obs),col=3)

a <- runif(1)
b <- runif(1)
c <- runif(1)
x <- runif(1)
integrate(function(u)pnorm(a*u+b*x+c)*dnorm(u),-5,5)$val
pnorm((b*x+c)/sqrt(1+a^2))

## routine to visually check conditional expectation on a continuous
## variable
## cond.exp.viz <- function(x,y,f,n.quantiles=1e2) {
##     bin <- cut(x,breaks=c(-Inf,quantile(x,probs=seq(0,1,len=n.quantiles)),Inf))
##     y.means <- aggregate(y,list(bin=bin),mean)$x
##     bin.centers <- aggregate(x,list(bin=bin),median)$x
##     plot(bin.centers,y.means,xlab='',asp=1)
##     curve(f,add=TRUE,col=2)
## }
cond.exp.viz <- function(x,y,f,n.quantiles=1e2) {
    ## browser()
    bin <- cut(x,breaks=c(-Inf,quantile(x,probs=seq(0,1,len=n.quantiles)),Inf))
    y.means <- aggregate(y,list(bin=bin),mean)[,2]
    bin.centers <- aggregate(x,list(bin=bin),median)[,2]
    plot(bin.centers,y.means,xlab='',asp=1)
    curve(f,add=TRUE,col=2)
}
require(mvtnorm)
n <- 1e4
Sigma <- matrix(runif(4),nrow=2)
Sigma <- Sigma%*%t(Sigma)
## Sigma <- diag(c(1,1))
mu <- runif(2)*4
xy <- rmvnorm(n,mu,Sigma)
x <- xy[,1]; y <- xy[,2]
## bin <- cut(x,breaks=c(-Inf,quantile(x,probs=seq(0,1,len=1000)),Inf))
## y.means <- aggregate(y,list(bin=bin),mean)$x
## bin.centers <- aggregate(x,list(bin=bin),median)$x
## plot(bin.centers,y.means,xlab='',asp=1)
## abline(a=mu[2]-Sigma[1,2]/Sigma[1,1]*mu[1],b=Sigma[1,2]/Sigma[1,1])
cond.exp <- function(x)mu[2]+Sigma[1,2]/Sigma[1,1]*(x-mu[1])
cond.exp.viz(x,y,f=cond.exp)

hist.paired <-  function(x,y,n.quantiles=1e2,new.plot=TRUE,...) {
    bin <- cut(x,breaks=c(-Inf,quantile(x,probs=seq(0,1,len=n.quantiles)),Inf))
    y.means <- aggregate(y,list(bin=bin),mean)$x
    bin.centers <- aggregate(x,list(bin=bin),median)$x
    if(new.plot) plot(bin.centers,y.means,xlab='',asp=1,...) else points(bin.centers,y.means) 
    ## curve(f,add=TRUE,col=2)
}

require(mvtnorm)
n <- 1e4
Sigma <- matrix(runif(4),nrow=2)
Sigma <- Sigma%*%t(Sigma)
## Sigma <- diag(c(1,1))
mu <- runif(2)
xy <- rmvnorm(n,mu,Sigma)
x <- xy[,1]; y <- xy[,2]
beta0 <- runif(1)
cond.exp <- function(x,mu,Sigma)pnorm((x+beta0+mu[2]+Sigma[1,2]/Sigma[1,1]*(x-mu[1]))/sqrt(1+Sigma[2,2]-Sigma[1,2]^2/Sigma[1,1]))
cond.exp.viz(x=x,y=pnorm(x+y+beta0),f=function(x)cond.exp(x,mu,Sigma))

cond.exp.viz <- function(x,y,f,n.quantiles=1e2,...) {
    ## browser()
    bin <- cut(x,breaks=c(-Inf,quantile(x,probs=seq(0,1,len=n.quantiles)),Inf))
    y.means <- aggregate(y,list(bin=bin),mean)[,2]
    bin.centers <- aggregate(x,list(bin=bin),median)[,2]
    plot(bin.centers,y.means,xlab='',asp=1,...)
    curve(f,add=TRUE,col=2)
}
require(mvtnorm)
n <- 1e4
p.x <- 3
p.y <- 3
Sigma <- matrix(runif((p.x+p.y)^2),nrow=p.x+p.y)
Sigma <- Sigma%*%t(Sigma)
Sigma.xx <- Sigma[1:p.x,1:p.x]
Sigma.yy <- Sigma[(p.x+1):(p.x+p.y),(p.x+1):(p.x+p.y)]
Sigma.xy <- Sigma[1:p.x,(p.x+1):(p.x+p.y)]
mu.x <- runif(p.x); mu.y <- runif(p.y)
mu <- c(mu.x,mu.y)
xy <- rmvnorm(n,mu,Sigma)
x <- xy[,1:p.x]; y <- xy[,(p.x+1):(p.x+p.y)]
beta.0 <- runif(1)
beta.x <- runif(p.x); beta.y <- runif(p.y)
beta <- c(beta.x,beta.y)
cond.exp <- function(u) pnorm ( ( u+beta.0+mu.y%*%beta.y + (t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)*(u-mu.x%*%beta.x)  ) / sqrt(  1+t(beta.y)%*%Sigma.yy%*%beta.y-(t(beta.x)%*%Sigma.xy%*%beta.y)^2/(t(beta.x)%*%Sigma.xx%*%beta.x) ) )
denom <- sqrt(1+t(beta.y)%*%(Sigma.yy-Sigma.xy%*%solve(Sigma.xx)%*%Sigma.xy)%*%beta.y)
cond.exp <- function(u)beta.0+t(beta.x)%*%x+t(beta.y)%*%mu.y+t(beta.y)%*%Sigma.xy%*%solve(Sigma.xx)%*%(x-mu.x)
cond.exp.viz(x=x%*%beta.x, y=pnorm(beta.0+x%*%beta.x+y%*%beta.y),f=cond.exp,1e3,cex=1/2)


dd

## ## cov(x%*%beta.x,y%*%beta.y)
## ## t(beta.x)%*%Sigma.xy%*%beta.y
## var(x%*%beta.x)
## t(beta.x)%*%Sigma.xx%*%beta.x
## var(y%*%beta.y)
## t(beta.y)%*%Sigma.yy%*%beta.y
## mean(y%*%beta.y)
## beta.y%*%mu.y
## mean(x%*%beta.x)
## beta.x%*%mu.x
## cond.exp <- function(x,mu,Sigma)pnorm((x+beta0+mu[2]+Sigma[1,2]/Sigma[1,1]*(x-mu[1]))/sqrt(1+Sigma[2,2]-Sigma[1,2]^2/Sigma[1,1]))
## mu.bar <- c(beta.x%*%mu.x,beta.y%*%mu.y)
## Sigma.bar <- matrix(c(t(beta.x)%*%Sigma.xx%*%beta.x,t(beta.x)%*%Sigma.xy%*%beta.y,t(beta.x)%*%Sigma.xy%*%beta.y,t(beta.y)%*%Sigma.yy%*%beta.y),nrow=2)
## cond.exp.viz(x=x%*%beta.x, y=pnorm(beta.0+x%*%beta.x+y%*%beta.y),f=function(x)cond.exp(x,mu.bar,Sigma.bar),n=1e3,cex=1/2)



## 12b formula for reduced model coefficients
require(mvtnorm)
set.seed(1)
n <- 1e2
p.x <- 3
p.y <- 3
Sigma <- matrix(runif((p.x+p.y)^2),nrow=p.x+p.y)
Sigma <- Sigma%*%t(Sigma)
Sigma.xx <- Sigma[1:p.x,1:p.x]
Sigma.yy <- Sigma[(p.x+1):(p.x+p.y),(p.x+1):(p.x+p.y)]
Sigma.xy <- Sigma[1:p.x,(p.x+1):(p.x+p.y)]
Sigma.yx <- t(Sigma.xy)
mu.x <- runif(p.x); mu.y <- runif(p.y)
mu <- c(mu.x,mu.y)
beta.0 <- runif(1)/10
beta.x <- runif(p.x)/10; beta.y <- runif(p.y)/10
beta <- c(beta.x,beta.y)
## denom <- sqrt(  1+t(beta.y)%*%Sigma.yy%*%beta.y-(t(beta.x)%*%Sigma.xy%*%beta.y)^2/(t(beta.x)%*%Sigma.xx%*%beta.x) ) 
## beta.0.reduced <- ( beta.0 + t(beta.y)%*%mu.y - (t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)*t(beta.x)%*%mu.x ) / denom
## beta.x.reduced <- (1+(t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)) / denom * beta.x
## beta.reduced <- c(beta.0.reduced,beta.x.reduced)
denom <- sqrt(1+t(beta.y)%*%(Sigma.yy-Sigma.yx%*%solve(Sigma.xx)%*%Sigma.xy)%*%beta.y)
beta.0.reduced <- beta.0 + t(beta.y)%*%mu.y - t(beta.y)%*%Sigma.yx%*%solve(Sigma.xx)%*%mu.x
beta.x.reduced <- t(beta.x)+t(beta.y)%*%Sigma.yx%*%solve(Sigma.xx)
beta.reduced <- c(beta.0.reduced,beta.x.reduced) / c(denom)
ns <- round(seq(1e2,1e3,len=10))
by.n <- sapply(ns, function(n) {
    coef.hats <- replicate(1e2, {
        xy <- rmvnorm(n,mu,Sigma)
        x <- xy[,1:p.x]; y <- xy[,(p.x+1):(p.x+p.y)]
        risk <- pnorm(beta.0 + xy%*%beta)
        g <- rbinom(n,1,risk)
        glm0 <- glm(g~x, family=binomial(link='probit'))
        coef(glm0)
    })
    ## hist((coef.hats-beta.reduced)[2,])
    ## ## hist(coef.hats)
    ## denom <- sqrt(  1+t(beta.y)%*%Sigma.yy%*%beta.y-(t(beta.x)%*%Sigma.xy%*%beta.y)^2/(t(beta.x)%*%Sigma.xx%*%beta.x) ) 
    ## beta.0.reduced <- ( beta.0 + t(beta.y)%*%mu.y - (t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)*t(beta.x)%*%mu.x ) / denom
    ## ## abline(v=beta.0.reduced,col=2)
    ## ## abline(v=mean(coef.hats),col=3)
    ## ## ( beta.0+mu.y%*%beta.y + (t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)*(-mu.x%*%beta.x)  ) / denom
    ## beta.x.reduced <- ( beta.x+(t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)*beta.x  ) / denom
    ## beta.reduced <- c(beta.0.reduced,beta.x.reduced)
    ## hist(coef.hats-beta.reduced)
    sum(apply(coef.hats-beta.reduced,1,sd))
})
plot(ns,by.n)



## dd



## encapsulating. thought about not distinguishing the intercept term
## \beta_0, but the marginalizing creates an intercept
## term. [update:shouldnt matter, bc auc will ignore it anyway]
set.seed(1)
require(mvtnorm)
n <- 1e2
p.x <- 3
p.y <- 3
Sigma <- matrix(runif((p.x+p.y)^2),nrow=p.x+p.y)
Sigma <- Sigma%*%t(Sigma)
Sigma.xx <- Sigma[1:p.x,1:p.x]
Sigma.yy <- Sigma[(p.x+1):(p.x+p.y),(p.x+1):(p.x+p.y)]
Sigma.xy <- Sigma[1:p.x,(p.x+1):(p.x+p.y)]
Sigma.yx <- t(Sigma.xy)
mu.x <- runif(p.x); mu.y <- runif(p.y)
mu <- c(mu.x,mu.y)
beta.0 <- runif(1)/10
beta.x <- runif(p.x)/10; beta.y <- runif(p.y)/10
beta <- c(beta.x,beta.y)
## returns the coefficients gamma in Phi(gamma.0+gamma^Tx) = E(
## Phi(beta.0+beta.x^Tx+beta.y^Ty) | x ) where (x,y) are multivariate
## normal (mu,Sigma), beta=c(beta.0,beta.x,beta.y), and p.reduced is
## the length of beta.x
probit.coef.reduced <- function(mu,Sigma,beta,p.reduced) {
    ## browser()
    p.full <- length(mu)
    stopifnot(dim(Sigma)==p.full)
    p.x <- p.reduced;  p.y <- p.full-p.x
    Sigma.xx <- Sigma[1:p.x,1:p.x]
    Sigma.yy <- Sigma[(p.x+1):(p.x+p.y),(p.x+1):(p.x+p.y)]
    Sigma.xy <- Sigma[1:p.x,(p.x+1):(p.x+p.y)]
    Sigma.yx <- t(Sigma.xy)
    mu.x <- mu[1:p.x]; mu.y <- mu[(p.x+1):(p.x+p.y)]
    beta.0 <- beta[1]
    beta.x <- beta[2:(p.x+1)]; beta.y <- beta[(p.x+2):(p.x+p.y+1)]
    ## denom <- sqrt(  1+t(beta.y)%*%Sigma.yy%*%beta.y-(t(beta.x)%*%Sigma.xy%*%beta.y)^2/(t(beta.x)%*%Sigma.xx%*%beta.x) )
    ## beta.0.reduced <- ( beta.0 + t(beta.y)%*%mu.y - (t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)*t(beta.x)%*%mu.x ) / denom
    ## beta.x.reduced <- (1+(t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)) / denom * beta.x
    ## beta.reduced <- c(beta.0.reduced,beta.x.reduced)
    denom <- sqrt(1+t(beta.y)%*%(Sigma.yy-Sigma.yx%*%solve(Sigma.xx)%*%Sigma.xy)%*%beta.y)
    beta.0.reduced <- beta.0 + t(beta.y)%*%mu.y - t(beta.y)%*%Sigma.yx%*%solve(Sigma.xx)%*%mu.x
    beta.x.reduced <- t(beta.x)+t(beta.y)%*%Sigma.yx%*%solve(Sigma.xx)
    beta.reduced <- c(beta.0.reduced,beta.x.reduced) / c(denom)
}
beta.reduced <- probit.coef.reduced(mu,Sigma,c(beta.0,beta),p.reduced=p.x)
ns <- round(seq(1e2,6e3,len=10))
by.n <- sapply(ns, function(n) {
    coef.hats <- replicate(1e2, {
        xy <- rmvnorm(n,mu,Sigma)
        x <- xy[,1:p.x]; y <- xy[,(p.x+1):(p.x+p.y)]
        risk <- pnorm(beta.0 + xy%*%beta)
        g <- rbinom(n,1,risk)
        glm0 <- glm(g~x, family=binomial(link='probit'))
        coef(glm0)
    })
    ## sum(apply(coef.hats-beta.reduced,1,sd))
    sum(rowMeans(abs(coef.hats-beta.reduced)))
})
plot(ns,by.n)
lm0 <- lm(log(by.n)~log(ns))
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)



## check that z-statistics look OK
require(mvtnorm)
source('misc.R')
set.seed(2)
n <-2e3
p.x <- 6
p.y <- 2
Sigma <- matrix(runif((p.x+p.y)^2),nrow=p.x+p.y)
Sigma <- Sigma%*%t(Sigma)
## Sigma <- diag(p.x+p.y)
Sigma.xx <- Sigma[1:p.x,1:p.x]
Sigma.yy <- Sigma[(p.x+1):(p.x+p.y),(p.x+1):(p.x+p.y)]
Sigma.xy <- Sigma[1:p.x,(p.x+1):(p.x+p.y)]
mu.x <- runif(p.x); mu.y <- runif(p.y)
mu <- c(mu.x,mu.y)
beta.0 <- runif(1)/5
beta.x <- runif(p.x)/5; beta.y <- runif(p.y)/5
beta.xy <- c(beta.x,beta.y)
beta.reduced <- probit.coef.reduced(mu,Sigma,c(beta.0,beta.xy),p.x)
## ns <- round(seq(1e2,1e4,len=30))
## by.n <- sapply(ns, function(n) {
z.stats <- replicate(1e3, {
    xy <- rmvnorm(n,mu,Sigma)
    x <- xy[,1:p.x]#; y <- xy[,(p.x+1):(p.x+p.y)]
    risk <- pnorm(beta.0 + xy%*%beta.xy)
    g <- rbinom(n,1,risk)
    ## glm0 <- glm(g~xy, family=binomial(link='probit'))
    ## solve(cov.rt)%*%(coef(glm0) - c(beta.0,beta.xy))
    glm0 <- glm(g~x, family=binomial(link='probit'))
    cov.rt <- expm::sqrtm(vcov(glm0))
    solve(cov.rt)%*%(coef(glm0) - beta.reduced)
})
## rowMeans(abs(coefs-c(beta.0,beta.xy)))
z.stats <- z.stats[,,,drop=TRUE]
qqnorm(z.stats)
abline(0,1,col=2)
apply(z.stats,1,var)
apply(z.stats,1,mean)
## hist(z.stats)

## ## debugging. see notes.tex.

## ## checking collapsibility for linear regression. requires
## ## marginalized covariates to be orthogonal to remaining covariates
## ## for right se.
## require(mvtnorm)
## ## set.seed(1)
## n <- 1e2
## beta <- runif(3)
## Sigma <- matrix(runif(4),2)
## Sigma <- Sigma%*%t(Sigma)
## ## Sigma <- diag(2)
## mu <- runif(2)
## z.stats <- replicate(1e3, {
##     x <- rmvnorm(n,mu,Sigma)
##     y <- cbind(1,x)%*%beta + rnorm(n)/5
##     beta.reduced <- beta[1:2] + beta[3]*c(mu[2]-Sigma[1,2]/Sigma[2,2]*mu[1], Sigma[1,2]/Sigma[2,2])
##     lm0 <- lm(y ~ x[,1])
##     vcov.root <- expm::sqrtm(vcov(lm0))
##     solve(vcov.root)%*%(coef(lm0)-beta.reduced)
## })
## qqnorm(z.stats); abline(0,1)



## ## collapsing probit when p.x=p.y=1. looks fine.
## ## set.seed(1)
## mu <- runif(2)
## Sigma <- matrix(runif(4),2)
## Sigma <- Sigma%*%t(Sigma)
## beta <- runif(3)
## x <- runif(1)
## var.cond <- (1-Sigma[1,2]^2/Sigma[1,1]/Sigma[2,2])*Sigma[2,2]
## ## mu.cond <- 
## integrate(function(y)pnorm(beta[1]+beta[2]*x+beta[3]*y)*dnorm(y,mu[2]+Sigma[1,2]/Sigma[1,1]*(x-mu[1]),sd=sqrt(var.cond)),-10,10)$val
## denom <- sqrt(1+beta[3]^2*var.cond )
## beta.reduced <- c(beta[1]+beta[3]*mu[2]-beta[3]*Sigma[1,2]/Sigma[1,1]*mu[1], beta[2]+beta[3]*Sigma[1,2]/Sigma[1,1]) / denom
## pnorm(beta.reduced%*%c(1,x))
## fn <- Vectorize(function(x)integrate(function(y)pnorm(beta[1]+beta[2]*x+beta[3]*y)*dnorm(y,mu[2]+Sigma[1,2]/Sigma[1,1]*(x-mu[1]),sd=sqrt(var.cond)),-10,10)$val)
## fn.try <- Vectorize(function(x)pnorm(beta.reduced%*%c(1,x)))
## curve(fn)
## curve(fn.try,add=TRUE,col=2)


## ## ## SO. 2 issues here--sigma[2,2] should be sigma[1,1] in cond mean formula. also issue discussed in notes about wrong conidtioning variable.
## ## require(mvtnorm)
## ## set.seed(1)
## ## n <- 2e3
## ## mu <- c(0,0)
## ## Sigma <- matrix(runif(4),2)
## ## Sigma <- Sigma%*%t(Sigma)
## ## beta <- c(-.3,.6,.7)
## ## var.cond <- 1-Sigma[1,2]^2/Sigma[1,1]/Sigma[2,2]
## ## denom <- sqrt(1+beta[3]^2*var.cond)
## ## beta.prime <- c(beta[1]+beta[3]*mu[2]-beta[3]*Sigma[1,2]/Sigma[2,2]*mu[1], beta[2]+beta[3]*Sigma[1,2]/Sigma[2,2]) / denom 
## ## z.stats <- replicate(100, {
## ##     xy <- rmvnorm(n,mu,Sigma)
## ##     x <- xy[,1]
## ##     risk <- pnorm(cbind(1,xy)%*%beta)
## ##     g <- rbinom(n,1,risk)
## ##     glm0 <- glm(g~x,family=binomial(link='probit'))
## ##     var.root <- expm::sqrtm(vcov(glm0)) ## matrix square root
## ##     solve(var.root)%*%(coef(glm0)-beta.prime)
## ## })
## ## qqnorm(z.stats)
## ## abline(0,1)


## require(mvtnorm)
## ## set.seed(1)
## n <- 2e3
## mu <- runif(2)
## Sigma <- matrix(runif(4),2)
## Sigma <- Sigma%*%t(Sigma)
## ## Sigma <- diag(2)
## beta <- runif(3)
## var.cond <- 1-Sigma[1,2]^2/Sigma[1,1]/Sigma[2,2]
## denom <- sqrt(1+beta[3]^2*var.cond )
## beta.reduced <- c(beta[1]+beta[3]*mu[2]-beta[3]*Sigma[1,2]/Sigma[1,1]*mu[1], beta[2]+beta[3]*Sigma[1,2]/Sigma[1,1]) / denom
## z.stats <- replicate(5e2, {
##     xy <- rmvnorm(n,mu,Sigma)
##     x <- xy[,1]
##     risk <- pnorm(cbind(1,xy)%*%beta)
##     g <- rbinom(n,1,risk)
##     glm0 <- glm(g~x,family=binomial(link='probit'))
##     var.root <- expm::sqrtm(vcov(glm0))
##     solve(var.root)%*%(coef(glm0)-beta.reduced)
## })
## z.stats <- z.stats[,,,drop=TRUE]
## qqnorm(z.stats)
## abline(0,1)
## ## apply(z.stats,1,var)
## ## apply(z.stats,1,mean)
## ## hist(z.stats)

## source('misc.R')
## require(mvtnorm)
## ## set.seed(1)
## n <- 2e3
## mu <- runif(2)
## Sigma <- matrix(runif(4),2)
## Sigma <- Sigma%*%t(Sigma)
## ## Sigma <- diag(2)
## beta <- runif(3)
## var.cond <- 1-Sigma[1,2]^2/Sigma[1,1]/Sigma[2,2]
## denom <- sqrt(1+beta[3]^2*var.cond )
## beta.reduced <- c(beta[1]+beta[3]*mu[2]-beta[3]*Sigma[1,2]/Sigma[1,1]*mu[1], beta[2]+beta[3]*Sigma[1,2]/Sigma[1,1]) / denom
## xy <- rmvnorm(n,mu,Sigma)
## x <- xy[,1]
## risk <- pnorm(cbind(1,xy)%*%beta)
## fn.try <- Vectorize(function(x)pnorm(beta.reduced%*%c(1,x)))
## cond.exp.viz(x=x,y=risk,f=fn.try)

## require(mvtnorm)
## ## set.seed(1)
## n <- 2e3
## mu <- runif(2)
## Sigma <- matrix(runif(4),2)
## Sigma <- Sigma%*%t(Sigma)
## ## Sigma <- diag(2)
## beta <- runif(3)
## var.cond <- 1-Sigma[1,2]^2/Sigma[1,1]/Sigma[2,2]
## denom <- sqrt(1+beta[3]^2*var.cond )
## beta.reduced <- c(beta[1]+beta[3]*mu[2]-beta[3]*Sigma[1,2]/Sigma[1,1]*mu[1], beta[2]+beta[3]*Sigma[1,2]/Sigma[1,1]) / denom
## ns <- round(seq(1e2,5e3,len=20))
## by.n <- sapply(ns, function(n) {
##     diffs <- replicate(1e2, {
##         xy <- rmvnorm(n,mu,Sigma)
##         x <- xy[,1]
##         risk <- pnorm(cbind(1,xy)%*%beta)
##         g <- rbinom(n,1,risk)
##         glm0 <- glm(g~x,family=binomial(link='probit'))
##         ## var.root <- expm::sqrtm(vcov(glm0))
##         coef(glm0)-beta.reduced
##     })
##     diffs[1,]
## })
## mad <- apply(abs(by.n),2,mean)
## plot(ns,mad)


## 12c influence function and sandwich var
set.seed(1)
require(mvtnorm)
require(sandwich)
n <- 5e3
p <- 3
beta <- runif(p+1)/5
ns <- round(seq(1e2,1e3,len=30))
by.n <- sapply(ns, function(n) {
    hats <- replicate(1e2, {
        x <- cbind(1,matrix(runif(n*p),ncol=p))
        risk <- pnorm(x%*%beta)
        y <- rbinom(n,1,risk)
        glm0 <- glm(y~x-1,family=binomial(link='probit'))
        risk.hat <- predict(glm0)
        scores <- dnorm(risk.hat)*(y/pnorm(risk.hat)/(1-pnorm(risk.hat)) - 1/(1-pnorm(risk.hat))) * x
        W <- diag(dnorm(risk.hat)^2/pnorm(risk.hat)/(1-pnorm(risk.hat)))
        hessian <- -t(x)%*%W%*%x / n
        infl <- solve(hessian)%*%t(scores) / sqrt(n)
        ## infl <- solve(bread(glm0))%*%t(scores) / sqrt(n)
        var.sandwich <- diag(cov(t(infl)))
        ## var.sandwich <- solve(hessian)%*%meat(glm0)%*%solve(hessian)
        ## var.sandwich <- diag(sandwich(glm0))
        ## var.sandwich <- diag(vcov(glm0)) # !!
        coefs <- coef(glm0)
        list(var.sandwich=var.sandwich,coefs=coefs)
    },simplify=FALSE)
    coefs <- sapply(hats, function(lst)lst$coefs)
    var.observed <- apply(coefs,1,var)
    var.sandwich <- rowMeans(sapply(hats,function(lst)lst$var.sandwich))
    var.observed - var.sandwich
})
diffs <- colSums(by.n)
coef(lm(log(abs(diffs))~log(ns)))
plot(ns,ns*diffs)

## bread/hessian is almost the same as sandwich library's value,
## but meat/outer product is a little off (1e-4 discrepancy with
## these settings)
set.seed(1)
require(mvtnorm)
require(sandwich)
n <- 5e1
p <- 3
beta <- runif(p+1)/5
## ns <- round(seq(1e2,1e3,len=30))
## by.n <- sapply(ns, function(n) {
##     hats <- replicate(1e2, {
x <- cbind(1,matrix(runif(n*p),ncol=p))
risk <- pnorm(x%*%beta)
y <- rbinom(n,1,risk)
glm0 <- glm(y~x-1,family=binomial(link='probit'))
risk.hat <- predict(glm0)
scores <- dnorm(risk.hat)*(y/pnorm(risk.hat)/(1-pnorm(risk.hat)) - 1/(1-pnorm(risk.hat))) * x
W <- diag(dnorm(risk.hat)^2/pnorm(risk.hat)/(1-pnorm(risk.hat)))
hessian <- -t(x)%*%W%*%x
bread.try <- -solve(hessian)
## infl <- solve(hessian)%*%t(scores) / sqrt(n)
bread(glm0) - n*vcov(glm0)
vcov(glm0)-bread.try
cov(scores) - meat(glm0)



set.seed(1)
require(mvtnorm)
require(sandwich)
n <- 5e2
p <- 3
beta <- runif(p+1)/5
ns <- round(seq(1e2,1e3,len=30))
by.n <- sapply(ns, function(n) {
    cat('.')
    coef.hats <- replicate(n, {
        x <- cbind(1,matrix(runif(n*p),ncol=p))
        risk <- pnorm(x%*%beta)
        y <- rbinom(n,1,risk)
        glm0 <- glm(y~x-1,family=binomial(link='probit'))
        coef(glm0)
    })
    var.observed <- apply(sqrt(n)*coef.hats,1,var)
    var.sandwiches <- replicate(1e2, {
        x <- cbind(1,matrix(runif(n*p),ncol=p))
        risk <- pnorm(x%*%beta)
        y <- rbinom(n,1,risk)
        glm0 <- glm(y~x-1,family=binomial(link='probit'))
        ## risk.hat <- predict(glm0)
        ## scores <- dnorm(risk.hat)*(y/pnorm(risk.hat)/(1-pnorm(risk.hat)) - 1/(1-pnorm(risk.hat))) * x
        ## W <- diag(dnorm(risk.hat)^2/pnorm(risk.hat)/(1-pnorm(risk.hat)))
        ## hessian <- -t(x)%*%W%*%x / n
        ## infl <- solve(hessian)%*%t(scores) / sqrt(n)
        ## infl <- solve(bread(glm0))%*%t(scores) / sqrt(n)
        ## diag(cov(t(infl)))
        ## var.sandwich <- solve(hessian)%*%meat(glm0)%*%solve(hessian)
        ## var.sandwich <- diag(sandwich(glm0))
        diag(vcov(glm0)) # !!
        ## coefs <- coef(glm0)
        ## list(var.sandwich=var.sandwich,coefs=coefs)
    })
    ## coefs <- sapply(hats, function(lst)lst$coefs)
    ## var.observed <- apply(coefs,1,var)
    var.sandwich <- n*rowMeans(var.sandwiches)
    ## var.observed - var.sandwich
    cbind(var.observed=var.observed,var.sandwich=var.sandwich)
},simplify=FALSE)
by.n <- simplify2array(by.n)
plot(ns, (by.n[,'var.observed',] - by.n[,'var.sandwich',])[1,])
plot(ns,by.n[1,'var.observed',])
points(ns,by.n[1,'var.sandwich',],col=2)


## checking coverage
options(warn=-1)
set.seed(1)
n <- 1e2
p <- 3
alpha <- .1
beta <- runif(p+1)/5
contrast <- runif(length(beta))
## ns <- round(seq(1e2,1e3,len=30))
## by.n <- sapply(ns, function(n) {
##     hats <- replicate(1e2, {
cover <- replicate(1e3, {
    x <- cbind(1,matrix(runif(n*p),ncol=p))
    risk <- pnorm(x%*%beta)
    y <- rbinom(n,1,risk)
    glm0 <- glm(y~x-1,family=binomial(link='probit'))
    ## risk.hat <- predict(glm0)
    ## scores <- dnorm(risk.hat)*(y/pnorm(risk.hat)/(1-pnorm(risk.hat)) - 1/(1-pnorm(risk.hat))) * x
    ## W <- diag(dnorm(risk.hat)^2/pnorm(risk.hat)/(1-pnorm(risk.hat)))
    ## hessian <- -t(x)%*%W%*%x
    ## bread.try <- -solve(hessian)
    ## ## infl <- solve(hessian)%*%t(scores) / sqrt(n)
    ## bread(glm0) - n*vcov(glm0)
    ## vcov(glm0)-bread.try
    ## cov(scores) - meat(glm0)
    vcov.beta <- sandwich(glm0)
    var.contrast <- t(contrast)%*%vcov.beta%*%contrast
    beta.hat <- coef(glm0)
    CI <- beta.hat%*%contrast + c(-1,1)*sqrt(var.contrast)*qnorm(alpha/2)
    prod(CI - beta%*%contrast) < 0
})
mean(cover)


## checking z-stats. also encapsulating routine to get influence
## function from probit glm object.
options(warn=-1)
set.seed(1)
infl.hat.probit <- function(glm0) {
    x <- model.matrix(glm0)
    y <- glm0$y
    risk.hat <- predict(glm0)
    scores <- dnorm(risk.hat)*(y/pnorm(risk.hat)/(1-pnorm(risk.hat)) - 1/(1-pnorm(risk.hat))) * x
    W <- diag(dnorm(risk.hat)^2/pnorm(risk.hat)/(1-pnorm(risk.hat)))
    hessian <- -t(x)%*%W%*%x / n
    bread <- -solve(hessian)
    infl <- bread%*%t(scores)
}    
n <- 1e2
p <- 6
## alpha <- .1
beta <- runif(p+1)/5
contrast <- runif(length(beta))
## ns <- round(seq(1e2,1e3,len=30))
## by.n <- sapply(ns, function(n) {
##     hats <- replicate(1e2, {
z.stats <- replicate(1e3, {
    x <- cbind(1,matrix(runif(n*p),ncol=p))
    risk <- pnorm(x%*%beta)
    y <- rbinom(n,1,risk)
    glm0 <- glm(y~x-1,family=binomial(link='probit'))
    ## risk.hat <- predict(glm0)
    ## scores <- dnorm(risk.hat)*(y/pnorm(risk.hat)/(1-pnorm(risk.hat)) - 1/(1-pnorm(risk.hat))) * x
    ## W <- diag(dnorm(risk.hat)^2/pnorm(risk.hat)/(1-pnorm(risk.hat)))
    ## hessian <- -t(x)%*%W%*%x / n
    ## bread.try <- -solve(hessian)
    ## ## ## infl <- solve(hessian)%*%t(scores) / sqrt(n)
    ## ## bread(glm0) - n*vcov(glm0)
    ## ## vcov(glm0)-bread.try
    ## ## cov(scores) - meat(glm0)
    ## ## infl <- bread(glm0)%*%t(scores)
    ## infl <- bread.try%*%t(scores)
    ## ## vcov.beta <- sandwich(glm0)
    infl <- infl.hat.probit(glm0)
    vcov.beta <- cov(t(infl))/n
    var.contrast <- t(contrast)%*%vcov.beta%*%contrast
    beta.hat <- coef(glm0)
    (beta.hat-beta)%*%contrast / sqrt(var.contrast)
})
qqnorm(z.stats)
abline(0,1)


## var matrix of full and reduced coef estimates jointly. looks
## fine. also, doesnt seem to make a difference whether or not
## conditioning on x.

require(mvtnorm)
source('misc.R')
set.seed(1)
n <- 3e3
## p.x <- 3
## p.y <- 3
p.full <- 6
p.reduced <- 3
Sigma <- matrix(runif((p.full)^2),nrow=p.full)
Sigma <- Sigma%*%t(Sigma)
mu <- runif(p.full)#c(mu.x,mu.y)
beta.full <- runif(p.full+1)/5
beta.reduced <- probit.coef.reduced(mu,Sigma,beta.full,p.reduced)
    ## x <- rmvnorm(n,mu,Sigma)
    ## x.reduced <- x[,1:p.reduced]
    ## risk <- pnorm(cbind(1,x)%*%beta.full)
    ## g <- rbinom(n,1,risk)
z.stats <- replicate(1e2, {
    ## browser()
    x <- rmvnorm(n,mu,Sigma)
    x.reduced <- x[,1:p.reduced]
    risk <- pnorm(cbind(1,x)%*%beta.full)
    g <- rbinom(n,1,risk)
    glm.full <- glm(g~x, family=binomial(link='probit'))
    ## infl.full <- infl.hat.probit(glm.full)
    infl.full <- infl.probit(model.matrix(glm.full),g,coef(glm.full))
    ## cov.root.full <- expm::sqrtm(cov(t(infl.full))/n)
    ## cov.root.full <- expm::sqrtm(vcov(glm.full))
    ## solve(cov.root.full)%*%(coef(glm.full) - beta.full)
    glm.reduced <- glm(g~x.reduced, family=binomial(link='probit'))
    ## infl.reduced <- infl.hat.probit(glm.reduced)
    infl.reduced <- infl.probit(model.matrix(glm.reduced),g,coef(glm.reduced))
    ## cov.root.reduced <- expm::sqrtm(vcov(glm.reduced))
    ## cov.root.reduced <- expm::sqrtm(cov(t(infl.reduced))/n)
    ## solve(cov.root.reduced)%*%(coef(glm.reduced) - beta.reduced)
    infl.combined <- rbind(infl.full,infl.reduced)
    cov.root.combined <- expm::sqrtm(cov(t(infl.combined))/n)
    solve(cov.root.combined)%*%(c(coef(glm.full),coef(glm.reduced)) - c(beta.full,beta.reduced))
})
## rowMeans(abs(coefs-c(beta.0,beta.xy)))
qqnorm(z.stats)
abline(0,1,col=2)

## checking influence function for probit beta.hat is actually
## o(1/sqrt(n))
start <- Sys.time()
require(mvtnorm)
source('misc.R')
set.seed(1)
n <- 1.5e2
p.full <- 6
p.reduced <- 3
Sigma <- matrix(runif((p.full-1)^2),nrow=p.full-1)
Sigma <- Sigma%*%t(Sigma)
Sigma <- rbind(0,cbind(0,Sigma))
Sigma.reduced <- Sigma[1:p.reduced,1:p.reduced]
mu <- c(1,runif(p.full-1))
mu.reduced <- mu[1:p.reduced]
beta.full <- runif(p.full)/5
beta.reduced <- probit.coef.reduced(mu[-1],Sigma[-1,-1],beta.full,p.reduced-1)
beta.combined <- c(beta.full,beta.reduced)
auc.full <- auc.probit(mu,Sigma,beta.full)
auc.reduced <- auc.probit(mu[1:(p.reduced)],Sigma[1:(p.reduced),1:(p.reduced)],beta.reduced)
diff.coef <- auc.full-auc.reduced
ns <- round(seq(1.5e2,1e3,len=20))
by.n <- sapply(ns, function(n) {
    cat('.')
    diffs <- replicate(1e2, {
        x <- rmvnorm(n,mu,Sigma)
        x <- matrix(runif(n*p.full),ncol=p.full)
        ## x.reduced <- x[,1:p.reduced]
        risk <- pnorm(x%*%beta.full)
        g <- rbinom(n,1,risk)
        glm.full <- glm(g~x-1, family=binomial(link='probit'))
        ## glm.reduced <- glm(g~x.reduced-1, family=binomial(link='probit'))
        beta.hat.full <- coef(glm.full)
        obs <- beta.hat.full - beta.full
        try <- rowMeans(infl.probit(x,g,beta.full))
        try - obs
    })
    diffs[2,]
    ## diffs
})
Sys.time()-start
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)


dd


## debugging. probably OK to delete.

## require(mvtnorm)
## ## options(warn=-1)
## probit.coef.reduced <- function(mu,Sigma,beta,p.reduced) {
##     ## browser()
##     p.full <- length(mu)
##     stopifnot(dim(Sigma)==p.full)
##     p.x <- p.reduced;  p.y <- p.full-p.x
##     Sigma.xx <- Sigma[1:p.x,1:p.x]
##     Sigma.yy <- Sigma[(p.x+1):(p.x+p.y),(p.x+1):(p.x+p.y)]
##     Sigma.xy <- Sigma[1:p.x,(p.x+1):(p.x+p.y)]
##     mu.x <- mu[1:p.x]; mu.y <- mu[(p.x+1):(p.x+p.y)]
##     beta.0 <- beta[1]
##     beta.x <- beta[2:(p.x+1)]; beta.y <- beta[(p.x+2):(p.x+p.y+1)]
##     denom <- sqrt(  1+t(beta.y)%*%Sigma.yy%*%beta.y-(t(beta.x)%*%Sigma.xy%*%beta.y)^2/(t(beta.x)%*%Sigma.xx%*%beta.x) )
##     beta.0.reduced <- ( beta.0 + t(beta.y)%*%mu.y - (t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)*t(beta.x)%*%mu.x ) / denom
##     beta.x.reduced <- (1+(t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)) / denom * beta.x
##     beta.reduced <- c(beta.0.reduced,beta.x.reduced)
## }
## ## source('misc.R')
## n <- 4e3
## p.full <- 6 # no. of covariates other than intercept in full model
## p.reduced <- 3
## beta.full <- runif(p.full+1)/5
## mu <- runif(p.full)
## Sigma <- matrix(runif(p.full^2),nrow=p.full)
## Sigma <- Sigma%*%t(Sigma)
## beta.reduced <- probit.coef.reduced(mu,Sigma,beta.full,p.reduced)
## ns <- round(seq(1e2,3e3,len=20))
## ## by.n <- sapply(ns, function(n) {
## ## cat('.')
## diffs <- replicate(1e2, {
##     x <- cbind(1,rmvnorm(n,mu,Sigma))
##     risk <- pnorm(x%*%beta.full)
##     y <- rbinom(n,1,risk)
##     x.reduced <- x[,2:(1+p.reduced)]
##     glm0 <- glm(y~x.reduced,family=binomial(link='probit'))
##     ## x.reduced <- x[,1:(1+p.reduced)]
##     ## glm0 <- glm(y~x.reduced-1,family=binomial(link='probit'))
##     coef(glm0)-beta.reduced
## })
## hist(diffs[1,])
## abline(v=mean(diffs[1,]),col=2)

## dd
## ## sum(rowMeans(abs(diffs)))
## ## })
## ## plot(ns,by.n)
## ## coef(lm(log(by.n)~log(ns)))

## require(mvtnorm)
## n <- 3e3
## p.x <- 3
## p.y <- 3
## p.full <- p.x+p.y
## Sigma <- matrix(runif((p.x+p.y)^2),nrow=p.x+p.y)
## Sigma <- Sigma%*%t(Sigma)
## ## mu.x <- runif(p.x); mu.y <- runif(p.y)
## ## mu <- c(mu.x,mu.y)
## mu <- runif(p.x+p.y)
## ## beta.0 <- runif(1)/10
## ## beta.x <- runif(p.x)/10; beta.y <- runif(p.y)/10
## ## beta <- c(beta.x,beta.y)
## beta <- beta.full <- runif(p.full+1)/10
## beta.0 <- beta.full[1]
## beta.x <- beta.full[2:(1+p.x)]
## probit.coef.reduced <- function(mu,Sigma,beta,p.reduced) {
##     ## browser()
##     p.full <- length(mu)
##     stopifnot(dim(Sigma)==p.full)
##     p.x <- p.reduced;  p.y <- p.full-p.x
##     Sigma.xx <- Sigma[1:p.x,1:p.x]
##     Sigma.yy <- Sigma[(p.x+1):(p.x+p.y),(p.x+1):(p.x+p.y)]
##     Sigma.xy <- Sigma[1:p.x,(p.x+1):(p.x+p.y)]
##     mu.x <- mu[1:p.x]; mu.y <- mu[(p.x+1):(p.x+p.y)]
##     beta.0 <- beta[1]
##     beta.x <- beta[2:(p.x+1)]; beta.y <- beta[(p.x+2):(p.x+p.y+1)]
##     denom <- sqrt(  1+t(beta.y)%*%Sigma.yy%*%beta.y-(t(beta.x)%*%Sigma.xy%*%beta.y)^2/(t(beta.x)%*%Sigma.xx%*%beta.x) )
##     beta.0.reduced <- ( beta.0 + t(beta.y)%*%mu.y - (t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)*t(beta.x)%*%mu.x ) / denom
##     beta.x.reduced <- (1+(t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)) / denom * beta.x
##     beta.reduced <- c(beta.0.reduced,beta.x.reduced)
## }
## beta.reduced <- probit.coef.reduced(mu,Sigma,c(beta.0,beta),p.reduced=p.x)
## ns <- round(seq(1e2,3e3,len=10))
## diffs <- replicate(1e2, {
##     ## xy <- rmvnorm(n,mu,Sigma)
##     ## x <- xy[,1:p.x]; y <- xy[,(p.x+1):(p.x+p.y)]
##     ## risk <- pnorm(beta.0 + xy%*%beta)
##     ## g <- rbinom(n,1,risk)    
##     ## glm0 <- glm(g~x, family=binomial(link='probit'))
##     ## x <- rmvnorm(n,mu,Sigma)
##     ## risk <- pnorm(beta.0+x%*%beta)
##     ## y <- rbinom(n,1,risk)
##     ## x.reduced <- x[,1:(p.reduced)]
##     ## glm0 <- glm(y~x.reduced,family=binomial(link='probit'))
##     x <- cbind(1,rmvnorm(n,mu,Sigma))
##     risk <- pnorm(x%*%c(beta))
##     y <- rbinom(n,1,risk)
##     x.reduced <- x[,2:(p.reduced+1)]
##     glm0 <- glm(y~x.reduced,family=binomial(link='probit'))
##     coef(glm0)-beta.reduced
## })
## hist(diffs[1,])
## abline(v=mean(diffs[1,]),col=2)

## require(mvtnorm)
## source('misc.R')
## n <- 1e2
## p.x <- 3
## p.y <- 3
## Sigma <- matrix(runif((p.x+p.y)^2),nrow=p.x+p.y)
## Sigma <- Sigma%*%t(Sigma)
## Sigma.xx <- Sigma[1:p.x,1:p.x]
## Sigma.yy <- Sigma[(p.x+1):(p.x+p.y),(p.x+1):(p.x+p.y)]
## Sigma.xy <- Sigma[1:p.x,(p.x+1):(p.x+p.y)]
## mu.x <- runif(p.x); mu.y <- runif(p.y)
## mu <- c(mu.x,mu.y)
## beta.0 <- runif(1)/10
## beta.x <- runif(p.x)/10; beta.y <- runif(p.y)/10
## beta <- c(beta.x,beta.y)
## ## beta.reduced <- probit.coef.reduced(mu,Sigma,c(beta.0,beta),p.reduced=p.x)
## denom <- sqrt(  1+t(beta.y)%*%Sigma.yy%*%beta.y-(t(beta.x)%*%Sigma.xy%*%beta.y)^2/(t(beta.x)%*%Sigma.xx%*%beta.x) )
## beta.0.reduced <- ( beta.0 + t(beta.y)%*%mu.y - (t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)*t(beta.x)%*%mu.x ) / denom 
## beta.x.reduced.factor <- (1+(t(beta.x)%*%Sigma.xy%*%beta.y)/(t(beta.x)%*%Sigma.xx%*%beta.x)) / denom 
## xy <- rmvnorm(n,mu,Sigma)
## x <- xy[,1:p.x]; y <- xy[,(p.x+1):(p.x+p.y)]
## risk <- pnorm(beta.0 + xy%*%beta)
## g <- rbinom(n,1,risk)
## risk.reduced <- x%*%beta.x
## cond.exp <- function(u)pnorm(beta.0.reduced+beta.x.reduced.factor*u)
## cond.exp.viz(x=risk.reduced,y=risk,f=cond.exp,1e3,cex=1/2)
## coef(glm(g~xy,family=binomial(link='probit'))) - c(beta.0,beta)

## source('misc.R')
## set.seed(1)
## p.x <- 3
## p.y <- 3
## Sigma <- matrix(runif((p.x+p.y)^2),nrow=p.x+p.y)
## Sigma <- Sigma%*%t(Sigma)
## Sigma.xx <- Sigma[1:p.x,1:p.x]
## Sigma.yy <- Sigma[(p.x+1):(p.x+p.y),(p.x+1):(p.x+p.y)]
## Sigma.xy <- Sigma[1:p.x,(p.x+1):(p.x+p.y)]
## mu.x <- runif(p.x); mu.y <- runif(p.y)
## mu <- c(mu.x,mu.y)
## beta.0 <- runif(1)/5
## beta.x <- runif(p.x)/5; beta.y <- runif(p.y)/5
## beta.xy <- c(beta.x,beta.y)
## beta.reduced <- probit.coef.reduced(mu,Sigma,c(beta.0,beta.xy),p.x)
## ns <- round(seq(1e2,1e4,len=30))
## by.n <- sapply(ns, function(n) {
##     coefs <- replicate(1e2, {
##         xy <- rmvnorm(n,mu,Sigma)
##         x <- xy[,1:p.x]; y <- xy[,(p.x+1):(p.x+p.y)]
##         risk <- pnorm(beta.0 + xy%*%beta.xy)
##         g <- rbinom(n,1,risk)
##         ## glm0 <- glm(g~xy, family=binomial(link='probit'))
##         glm0 <- glm(g~x, family=binomial(link='probit'))
##         coef(glm0)
##     })
##     ## rowMeans(abs(coefs-c(beta.0,beta.xy)))
##     rowMeans(abs(coefs-beta.reduced))
## })
## diffs <- colSums(by.n)
## diffs <- by.n[4,]
## plot(ns,diffs)
## coef(lm(log(diffs)~log(ns)))
## coefs <- coef(lm(log(diffs)~log(ns)))
## curve(exp(coefs[1]+coefs[2]*log(x)),add=TRUE)

## source('misc.R')
## set.seed(1)
## n <- 1e2
## p.x <- 3
## p.y <- 3
## Sigma <- matrix(runif((p.x+p.y)^2),nrow=p.x+p.y)
## Sigma <- Sigma%*%t(Sigma)
## Sigma.xx <- Sigma[1:p.x,1:p.x]
## Sigma.yy <- Sigma[(p.x+1):(p.x+p.y),(p.x+1):(p.x+p.y)]
## Sigma.xy <- Sigma[1:p.x,(p.x+1):(p.x+p.y)]
## mu.x <- runif(p.x); mu.y <- runif(p.y)
## mu <- c(mu.x,mu.y)
## beta.0 <- runif(1)/5
## beta.x <- runif(p.x)/5; beta.y <- runif(p.y)/5
## beta.xy <- c(beta.x,beta.y)
## beta.reduced <- probit.coef.reduced(mu,Sigma,c(beta.0,beta.xy),p.x)
## ## ns <- round(seq(1e2,1e4,len=30))
## ## by.n <- sapply(ns, function(n) {
## z.stats <- replicate(1e2, {
##     xy <- rmvnorm(n,mu,Sigma)
##     x <- xy[,1:p.x]; y <- xy[,(p.x+1):(p.x+p.y)]
##     risk <- pnorm(beta.0 + xy%*%beta.xy)
##     g <- rbinom(n,1,risk)
##     ## glm0 <- glm(g~xy, family=binomial(link='probit'))
##     glm0 <- glm(g~x, family=binomial(link='probit'))
##     cov.rt <- expm::sqrtm(vcov(glm0))
##     solve(cov.rt)%*%(coef(glm0) - beta.reduced)
## })
## ## rowMeans(abs(coefs-c(beta.0,beta.xy)))
## qqnorm(z.stats)
## abline(0,1,col=2)


## 12d formula for variance of #1 full and #2 reduced infl function,
## and #3 covariance of full and reduced influence functions, and #4
## combined full and reduced infl functions. conditioning on x, using
## true influence function.

require(mvtnorm)
source('misc.R')
set.seed(1)
n <- 1.5e3
## p.x <- 3
## p.y <- 3
p.full <- 6
p.reduced <- 3
Sigma <- matrix(runif((p.full)^2),nrow=p.full)
Sigma <- Sigma%*%t(Sigma)
mu <- runif(p.full)#c(mu.x,mu.y)
## mu.x <- mu[1:p.reduced]; mu.y <- 
## beta.0 <- runif(1)/5
## beta.x <- runif(p.x)/5; beta.y <- runif(p.y)/5
## beta.xy <- c(beta.x,beta.y)
## beta.reduced <- probit.coef.reduced(mu,Sigma,c(beta.0,beta.xy),p.x)
beta.full <- runif(p.full+1)/5
beta.reduced <- probit.coef.reduced(mu,Sigma,beta.full,p.reduced)
beta.combined <- c(beta.full,beta.reduced)
a <- runif(p.full+1)
b <- runif(p.reduced+1)
## ns <- round(seq(1e2,1e4,len=30))
## by.n <- sapply(ns, function(n) {
infl.probit <- function(x,g,beta) { 
    eta <- as.numeric(x%*%beta)
    scores <- dnorm(eta)*(g/pnorm(eta)/(1-pnorm(eta)) - 1/(1-pnorm(eta))) * x
    W <- diag(dnorm(eta)^2/pnorm(eta)/(1-pnorm(eta)))
    hessian <- -t(x)%*%W%*%x / n
    bread <- -solve(hessian)
    infl <- bread%*%t(scores)
}
x <- cbind(1,rmvnorm(n,mu,Sigma))
x.reduced <- x[,1:(p.reduced+1)]
risk <- pnorm(x%*%beta.full)
obs <- replicate(3e2, {
    ## browser()
    g <- rbinom(n,1,risk)
    infl.full <- infl.probit(x,g,beta.full)
    infl.reduced <- infl.probit(x.reduced,g,beta.reduced)
    t(a)%*%infl.full%*%t(infl.full)%*%a / n #1
    t(b)%*%infl.reduced%*%t(infl.reduced)%*%b / n #2
    t(a)%*%(infl.full%*%t(infl.reduced))%*%b / n #3
    infl.combined <- rbind(infl.full,infl.reduced) #4
    t(c(a,b))%*%(infl.combined%*%t(infl.combined))%*%c(a,b) / n #4
    glm.full <- glm(g~x-1, family=binomial(link='probit'))#5
    glm.reduced <- glm(g~x.reduced-1, family=binomial(link='probit')) #5
    beta.combined.hat <- c(coef(glm.full),coef(glm.reduced))  #5
    (t(c(a,b))%*%(sqrt(n)*(beta.combined.hat - beta.combined) ))^2  #5
    })
eta.full <- as.numeric(x%*%beta.full)
eta.reduced <- as.numeric(x.reduced%*%beta.reduced)
W.full <- diag(dnorm(eta.full)^2/pnorm(eta.full)/(1-pnorm(eta.full))) #1
vcov.full <- solve(t(x)%*%W.full%*%x / n) #1
cov.try <- t(a)%*%vcov.full%*%a #1
W.reduced <- diag(dnorm(eta.reduced)^2/pnorm(eta.reduced)/(1-pnorm(eta.reduced))) #2
vcov.reduced <- solve(t(x.reduced)%*%W.reduced%*%x.reduced / n) #2
cov.try <- t(b)%*%vcov.reduced%*%b #2
W.mixed <- diag(dnorm(eta.full)*dnorm(eta.reduced)/pnorm(eta.reduced)/(1-pnorm(eta.reduced))) #3
vcov.mixed <- vcov.full%*%t(x)%*%W.mixed%*%x.reduced %*% vcov.reduced / n #3
cov.try <- t(a)%*%vcov.mixed%*%b #3 
vcov.combined <- rbind(cbind(vcov.full,vcov.mixed),cbind(t(vcov.mixed),vcov.reduced)) #4
cov.try <- t(c(a,b))%*%vcov.combined%*%c(a,b) #4 and #5
hist(obs)
abline(v=cov.try,col=2)
abline(v=mean(obs),col=3)
c(cov.try,mean(obs))




##12e formula for AUC under probit model and its derivative wrt beta

## 12e-1 formula for P(w<w'|g=0,g'=1), P(g=1|w)=pnorm(w), w ~ N(mu,sigma^2)
source('misc.R')
set.seed(1)
n <- 3e2
bw <- 8
mu <- runif(1)
sigma <- runif(1)
auc.hats <- replicate(3e3, {
    w <- rnorm(n,mu,sigma)
    risk <- pnorm(w)
    g <- rbinom(n,1,risk)
    w.0 <- w[g==0]; w.1 <- w[g==1]
    auc.hat(w.0,w.1)
})
auc.probit <- function(mu,sigma,bw=8) {
    E.g <- pnorm(mu/sqrt(1+sigma^2))
    f.case <- function(w)pnorm(w)*dnorm((w-c(mu))/sigma)/c(sigma*E.g)
    f.control <- function(w)(1-pnorm(w))*dnorm((w-c(mu))/sigma)/c(sigma*(1-E.g))
    int.inner <- Vectorize(function(v)integrate(function(u)f.control(u),mu-bw*sigma,v)$val)
    integrate(function(v)f.case(v)*int.inner(v),mu-bw*sigma,mu+bw*sigma)$val
}
hist(auc.hats)
abline(v=auc.probit(mu,sigma),col=2)

## 12e-2 with w an index
require(mvtnorm)
source('misc.R')
set.seed(1)
p <- 4
n <- 1e2
beta <- runif(p+1)/2
mu <- c(1,runif(p))
Sigma <- matrix(runif((p+1)^2),p+1)
Sigma <- Sigma%*%t(Sigma)
Sigma[,1] <- Sigma[1,] <- 0
auc.hats <- replicate(1e3, {
    x <- rmvnorm(n,mu,Sigma)
    risk <- pnorm(x%*%beta)
    g <- rbinom(n,1,risk)
    x.0 <- x[g==0,]; x.1 <- x[g==1,]
    auc.hat(x.0%*%beta,x.1%*%beta)
})
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
## bw <- 8
## int.inner <- Vectorize(function(v)integrate(function(u)(1-pnorm(u))*dnorm((u-c(t(beta)%*%mu))/quad),t(beta)%*%mu-bw*sqrt(quad),v)$val)
## int.outer <- integrate(function(v) pnorm(v)*dnorm((v-c(t(beta)%*%mu))/sqrt(quad))*int.inner(v),t(beta)%*%mu-bw*sqrt(quad),t(beta)%*%mu+bw*sqrt(quad))$val
## E.g <- pnorm((t(beta)%*%mu)/sqrt(1+quad))
try <- auc.probit(beta%*%mu,quad)
hist(auc.hats)
abline(v=try,col=2)
abline(v=mean(auc.hats),col=3)



v <- 2
a <- runif(1)
t <- runif(1)
a <- t <- 1
integrate(function(u)dnorm(a*u)*dnorm(t*u),-6,v)$val
1/2/pi*integrate(function(u)exp(-1/2*(t^2+a^2)*u^2),-6,v)$val
1/sqrt(2*pi*(t^2+a^2))*pnorm(v)
1/sqrt(2*pi)* integrate(function(u)exp(-1/2*(sqrt(t^2+a^2)*u)^2),-10,v)$val
1/sqrt(t^2+a^2)*pnorm(v*(sqrt(t^2+a^2)))

## 12e-3 derivative of auc.probit(0,sigma) wrt sigma
source('misc.R')
bw <- 8
f <- Vectorize(function(sigma)auc.probit(0,sigma^2))
curve(f,.1,5)
sigma <- runif(1,.1,5)
E.g <- 1/2
f.case <- function(w)pnorm(w)*dnorm(w/sigma)*2/sigma
f.control <- function(w)(1-pnorm(w))*dnorm(w/sigma)*2/sigma
A <- -2/sigma*f(sigma)
int.inner <- Vectorize(function(v)integrate(function(u)f.control(u),-bw*sigma,v)$val)
B <- int.outer <- integrate(function(v)f.case(v)*v^2/sigma^3*int.inner(v),-bw*sigma,+bw*sigma)$val
int.inner <- Vectorize(function(v)integrate(function(u)f.control(u)*u^2/sigma^3,-bw*sigma,v)$val)
C <- int.outer <- integrate(function(v)f.case(v)*int.inner(v),-bw*sigma,+bw*sigma)$val
slope <- A+B+C
abline(b=slope,a=f(sigma)-slope*sigma,col=2)


## 12e-4 derivative of auc.probit(0,sqrt(t(beta)%*%Sigma%*%beta)) wrt beta
source('misc.R')
p <- 3
bw <- 8
## beta <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
f <- function(beta)auc.probit(0,t(beta)%*%Sigma%*%beta)
## f(runif(3))
auc.probit.prime <- function(beta,mu,Sigma) {
    f <- function(sigma)auc.probit(0,sigma)
    sigma <- as.numeric(sqrt(t(beta)%*%Sigma%*%beta))
    E.g <- 1/2
    f.case <- function(w)pnorm(w)*dnorm(w/sigma)*2/sigma
    f.control <- function(w)(1-pnorm(w))*dnorm(w/sigma)*2/sigma
    A <- -2/sigma*f(sigma)
    int.inner <- Vectorize(function(v)integrate(function(u)f.control(u),-bw*sigma,v)$val)
    B <- int.outer <- integrate(function(v)f.case(v)*v^2/sigma^3*int.inner(v),-bw*sigma,+bw*sigma)$val
    int.inner <- Vectorize(function(v)integrate(function(u)f.control(u)*u^2/sigma^3,-bw*sigma,v)$val)
    C <- int.outer <- integrate(function(v)f.case(v)*int.inner(v),-bw*sigma,+bw*sigma)$val
    scalar.part <- A+B+C
    scalar.part*t(beta)%*%Sigma/sigma
}
f.prime <- function(beta)auc.probit.prime(beta,mu,Sigma)
plot.matrix.deriv(f,f.prime,p,1)

## 12e-5 derivative of auc.probit(mu,sigma) wrt sigma
source('misc.R')
mu <- 0
sigma <- 3
bw <- 7
f.sigma <- Vectorize(function(sigma)auc.probit(mu,sigma^2))
E.g <- function(sigma)pnorm(mu/sqrt(1+sigma^2))
A <- function(sigma)sigma^2*E.g(sigma)*(1-E.g(sigma))
B <- function(v,sigma)pnorm(v)*dnorm((v-mu)/sigma)
C <- function(u,sigma)(1-pnorm(u))*dnorm((u-mu)/sigma)
A.prime <- function(sigma)2*sigma*E.g(sigma)*(1-E.g(sigma)) - sigma^2*dnorm(mu/sqrt(1+sigma^2))*mu*sigma/(1+sigma^2)^(3/2)*(1-E.g(sigma)) + sigma^2*E.g(sigma)*dnorm(mu/sqrt(1+sigma^2))*mu*sigma/(1+sigma^2)^(3/2)
## -A(sigma)^(-2)*(  2*sigma*E.g(sigma)*(1-E.g(sigma)) - sigma^2*dnorm(mu/sqrt(1+sigma^2))*mu*sigma/(1+sigma^2)^(3/2)*(1-E.g(sigma)) + sigma^2*E.g(sigma)*dnorm(mu/sqrt(1+sigma^2))*mu*sigma/(1+sigma^2)^(3/2) )
viz.deriv(f=function(sigma)1/A(sigma),deriv=function(sigma)-1/A(sigma)^2*A.prime(sigma),1,1)


## f.sigma.prime <- function(mu,sigma) {
##     E.g <- function(sigma)pnorm(mu/sqrt(1+sigma^2))
##     A <- function(sigma)sigma^2*E.g(sigma)*(1-E.g(sigma))
##     B <- function(v,sigma)pnorm(v)*dnorm((v-mu)/sigma)
##     C <- function(u,sigma)(1-pnorm(u))*dnorm((u-mu)/sigma)
##     A.prime <- function(sigma)2*sigma*E.g(sigma)*(1-E.g(sigma)) - sigma^2*dnorm(mu/sqrt(1+sigma^2))*mu*sigma/(1+sigma^2)^(3/2)*(1-E.g(sigma)) + sigma^2*E.g(sigma)*dnorm(mu/sqrt(1+sigma^2))*mu*sigma/(1+sigma^2)^(3/2)
##     int.inner <- Vectorize(function(v)integrate(function(u)C(u,sigma),mu-bw*sigma,v)$val)
##     int.outer <- integrate(function(v)B(v,sigma)*int.inner(v),mu-bw*sigma,mu+bw*sigma)$val
##     line.1 <- -1/A(sigma)^2*A.prime(sigma)*int.outer
##     int.inner <- Vectorize(function(v)integrate(function(u)pnorm(v)*dnorm((v-mu)/sigma)*(v-mu)^2/sigma^3*C(u,sigma) + B(v,sigma)*(1-pnorm(u))*dnorm((u-mu)/sigma)*(u-mu)^2/sigma^3,mu-bw*sigma,v)$val)
##     int.outer <- integrate(int.inner,mu-bw*sigma,mu+bw*sigma)$val
##     line.2 <- int.outer/A(sigma)
##     line.1+line.2
## }
f.sigma.prime <- function(mu,sigma,bw=7) {
    E.g <- function(sigma)pnorm(mu/sqrt(1+sigma^2))
    A <- function(mu,sigma)sigma^2*E.g(sigma)*(1-E.g(sigma))
    ## B <- function(v,sigma)pnorm(v)*dnorm((v-mu)/sigma)
    ## C <- function(u,sigma)(1-pnorm(u))*dnorm((u-mu)/sigma)
    A.sigma.prime <- function(sigma)2*sigma*E.g(sigma)*(1-E.g(sigma)) - sigma^2*dnorm(mu/sqrt(1+sigma^2))*mu*sigma/(1+sigma^2)^(3/2)*(1-E.g(sigma)) + sigma^2*E.g(sigma)*dnorm(mu/sqrt(1+sigma^2))*mu*sigma/(1+sigma^2)^(3/2)
    ## int.inner <- Vectorize(function(v)integrate(function(u)C(u,sigma),mu-bw*sigma,v)$val)
    ## int.outer <- integrate(function(v)B(v,sigma)*int.inner(v),mu-bw*sigma,mu+bw*sigma)$val
    ## line.1 <- -1/A(sigma)^2*A.prime(sigma)*int.outer
    ## int.inner <- Vectorize(function(v)integrate(function(u)pnorm(v)*dnorm((v-mu)/sigma)*(v-mu)^2/sigma^3*C(u,sigma) + B(v,sigma)*(1-pnorm(u))*dnorm((u-mu)/sigma)*(u-mu)^2/sigma^3,mu-bw*sigma,v)$val)
    ## int.outer <- integrate(int.inner,mu-bw*sigma,mu+bw*sigma)$val
    ## line.2 <- int.outer/A(sigma)
    ## line.1+line.2
    int.inner <- Vectorize(function(v)integrate(function(u)pnorm(v)*(1-pnorm(u))*dnorm((v-mu)/sigma)*dnorm((u-mu)/sigma) * (-A.sigma.prime(sigma)/A(mu,sigma)^2 + ((v-mu)^2+(u-mu)^2)/A(mu,sigma)/sigma^3   ),mu-bw*sigma,v)$val)
    int.outer <- integrate(int.inner,mu-bw*sigma,mu+bw*sigma)$val
    int.outer
}
## mu <- runif(1); sigma <- runif(1)
## f.sigma.prime(mu,sigma)
## f.sigma.prime.old(mu,sigma)
viz.deriv(f.sigma,function(sigma)f.sigma.prime(mu,sigma),1,1)

dd

## 12e-6 derivative wrt mu

f.mu <- Vectorize(function(mu)auc.probit(mu,sigma^2))
E.g <- function(mu)pnorm(mu/sqrt(1+sigma^2))
A <- function(mu)sigma^2*E.g(mu)*(1-E.g(mu))
B <- function(v,mu)pnorm(v)*dnorm((v-mu)/sigma)
C <- function(u,mu)(1-pnorm(u))*dnorm((u-mu)/sigma)
A.prime <- function(mu)sigma^2/sqrt(1+sigma^2)*dnorm(mu/sqrt(1+sigma^2))*(1-2*E.g(mu))
viz.deriv(f=function(mu)1/A(mu),deriv=function(mu)-1/A(mu)^2*A.prime(mu),1,1)

dd

## v <- 1
source('misc.R')
## bw <- 7
## f.mu.prime.old <- function(mu,sigma) {
##     E.g <- function(mu)pnorm(mu/sqrt(1+sigma^2))
##     A <- function(mu)sigma^2*E.g(mu)*(1-E.g(mu))
##     B <- function(v,mu)pnorm(v)*dnorm((v-mu)/sigma)
##     C <- function(u,mu)(1-pnorm(u))*dnorm((u-mu)/sigma)
##     A.prime <- function(mu)sigma^2/sqrt(1+sigma^2)*dnorm(mu/sqrt(1+sigma^2))*(1-2*E.g(mu))
##     int.inner <- Vectorize(function(v)integrate(function(u)C(u,mu),mu-bw*sigma,v)$val)
##     int.outer <- integrate(function(v)B(v,mu)*int.inner(v),mu-bw*sigma,mu+bw*sigma)$val
##     line.1 <- -1/A(mu)^2*A.prime(mu)*int.outer
##     int.inner <- Vectorize(function(v)integrate(function(u)pnorm(v)*dnorm((v-mu)/sigma)*(v-mu)/sigma^2*C(u,mu) + B(v,mu)*(1-pnorm(u))*dnorm((u-mu)/sigma)*(u-mu)/sigma^2,mu-bw*sigma,v)$val)
##     int.outer <- integrate(int.inner,mu-bw*sigma,mu+bw*sigma)$val
##     line.2 <- int.outer/A(mu)
##     line.1+line.2
## }
f.mu.prime <- function(mu,sigma,bw=7) {
    E.g <- function(mu)pnorm(mu/sqrt(1+sigma^2))
    A <- function(mu,sigma)sigma^2*E.g(mu)*(1-E.g(mu))
    A.mu.prime <- function(mu)sigma^2/sqrt(1+sigma^2)*dnorm(mu/sqrt(1+sigma^2))*(1-2*E.g(mu))
    int.inner <- Vectorize(function(v)integrate(function(u)pnorm(v)*(1-pnorm(u))*dnorm((v-mu)/sigma)*dnorm((u-mu)/sigma)*(-A.mu.prime(mu)/A(mu,sigma)^2+(v+u-2*mu)/A(mu,sigma)/sigma^2),mu-bw*sigma,v)$val)
    int.outer <- integrate(int.inner,mu-bw*sigma,mu+bw*sigma)$val
    int.outer
}
## mu <- runif(1); sigma <- runif(1)
## f.mu.prime(mu,sigma)
## f.mu.prime.old(mu,sigma)
viz.deriv(f.mu,function(mu)f.mu.prime(mu,sigma),1,1)

dd

## 12e-4 derivative of
## auc.probit(t(beta)%*%mu,sqrt(t(beta)%*%Sigma%*%beta)) wrt beta

## source('misc.R')
## p <- 3
## bw <- 8
## beta <- runif(p)
## mu <- runif(p)
## Sigma <- matrix(runif(p^2),p)
## Sigma <- Sigma%*%t(Sigma)
## f <- function(beta)auc.probit(t(beta)%*%mu,t(beta)%*%Sigma%*%beta)
## auc.probit.prime <- function(beta,mu,Sigma) {
##     ## f <- function(mu,sigma)auc.probit(mu,sigma)
##     ## sigma <- as.numeric(sqrt(t(beta)%*%Sigma%*%beta))
##     ## mu <- as.numeric(t(beta)%*%mu)
##     f.sigma.prime(as.numeric(t(beta)%*%mu),as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))*t(beta)%*%Sigma/sigma +
##     f.mu.prime(as.numeric(t(beta)%*%mu),as.numeric(sqrt(t(beta)%*%Sigma%*%beta)))*t(mu)
## }
## f.prime <- function(beta)auc.probit.prime(beta,mu,Sigma)
## viz.deriv(f,f.prime,p,1)


source('misc.R')
p <- 3
bw <- 8
beta <- runif(p)
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
f <- function(beta)auc.probit(t(beta)%*%mu,t(beta)%*%Sigma%*%beta)
auc.probit.prime <- function(beta,mu,Sigma) {
    f.mu.prime <- function(mu,sigma,bw=7) {
        E.g <- function(mu)pnorm(mu/sqrt(1+sigma^2))
        A <- function(mu,sigma)sigma^2*E.g(mu)*(1-E.g(mu))
        A.mu.prime <- function(mu)sigma^2/sqrt(1+sigma^2)*dnorm(mu/sqrt(1+sigma^2))*(1-2*E.g(mu))
        int.inner <- Vectorize(function(v)integrate(function(u)pnorm(v)*(1-pnorm(u))*dnorm((v-mu)/sigma)*dnorm((u-mu)/sigma)*(-A.mu.prime(mu)/A(mu,sigma)^2+(v+u-2*mu)/A(mu,sigma)/sigma^2),mu-bw*sigma,v)$val)
        int.outer <- integrate(int.inner,mu-bw*sigma,mu+bw*sigma)$val
        int.outer
    }
    f.sigma.prime <- function(mu,sigma,bw=7) {
        E.g <- function(sigma)pnorm(mu/sqrt(1+sigma^2))
        A <- function(mu,sigma)sigma^2*E.g(sigma)*(1-E.g(sigma))
        A.sigma.prime <- function(sigma)2*sigma*E.g(sigma)*(1-E.g(sigma)) - sigma^2*dnorm(mu/sqrt(1+sigma^2))*mu*sigma/(1+sigma^2)^(3/2)*(1-E.g(sigma)) + sigma^2*E.g(sigma)*dnorm(mu/sqrt(1+sigma^2))*mu*sigma/(1+sigma^2)^(3/2)
        int.inner <- Vectorize(function(v)integrate(function(u)pnorm(v)*(1-pnorm(u))*dnorm((v-mu)/sigma)*dnorm((u-mu)/sigma) * (-A.sigma.prime(sigma)/A(mu,sigma)^2 + ((v-mu)^2+(u-mu)^2)/A(mu,sigma)/sigma^3   ),mu-bw*sigma,v)$val)
        int.outer <- integrate(int.inner,mu-bw*sigma,mu+bw*sigma)$val
        int.outer
    }
    index.mu <- as.numeric(t(beta)%*%mu)
    index.sigma <- as.numeric(sqrt(t(beta)%*%Sigma%*%beta))
    f.sigma.prime(index.mu,index.sigma)*t(beta)%*%Sigma/index.sigma +
    f.mu.prime(index.mu,index.sigma)*t(mu)
}
f.prime <- function(beta)auc.probit.prime(beta,mu,Sigma)
viz.deriv(f,f.prime,p,1)






## 12f diff.coef.hat - diff.coef

## probably OK to delete below commented out code
## ## require(mvtnorm)
## ## start <- Sys.time()
## ## source('misc.R')
## ## set.seed(3)
## ## p <- 6
## ## ## n <- 1e3
## ## ns <- round(seq(1e2,5e2,len=2e1))
## ## beta <- runif(p)/4
## ## gamma <- runif(p)/4
## ## Sigma <- matrix(rnorm(p^2),nrow=p)
## ## Sigma <- Sigma%*%t(Sigma)
## ## mu.1 <- runif(p)
## ## mu.0 <- runif(p)
## ## by.n <- sapply(ns, function(n) {
## ##     cat('.')
## ##     pairs <- replicate(1e2, {
## ##         x.0 <- rmvnorm(n,mu.0,Sigma)
## ##         x.1 <- rmvnorm(n,mu.1,Sigma)
## ##         beta.hat <- beta+rnorm(p)/sqrt(n)#lda.coefs(x.0,x.1)
## ##         gamma.hat <- gamma+rnorm(p)/sqrt(n)
## ##         auc.hat.beta.hat <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat)
## ##         auc.hat.gamma.hat <- auc.hat(x.0%*%gamma.hat,x.1%*%gamma.hat)
## ##         diff.hat.coef.hat <- auc.hat.beta.hat-auc.hat.gamma.hat
## ##         auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
## ##         auc.gamma.hat <- auc.scores(gamma.hat,mu.1-mu.0,2*Sigma)
## ##         diff.coef.hat <- auc.beta.hat-auc.gamma.hat
## ##         diff <- auc.scores(beta,mu.1-mu.0,2*Sigma) - auc.scores(gamma,mu.1-mu.0,2*Sigma)
## ##         auc.hajek.beta <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta.hat,Sigma)
## ##         auc.hajek.gamma <- auc.scores.hajek(x.0,x.1,mu.0,mu.1,gamma.hat,Sigma)
## ##         diff.hajek.coef <- auc.hajek.beta - auc.hajek.gamma
## ##         diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef
## ##     })
## ## })
## ## Sys.time() - start
## ## sds <- apply(by.n,2,sd)
## ## lm0 <- lm(log(sds) ~ log(ns))
## ## lm0


## require(mvtnorm)
## source('misc.R')
## set.seed(1)
## n <- 1.5e3
## ## p.x <- 3
## ## p.y <- 3
## p.full <- 6
## p.reduced <- 3
## Sigma <- matrix(runif((p.full)^2),nrow=p.full)
## Sigma <- Sigma%*%t(Sigma)
## Sigma <- rbind(0,cbind(0,Sigma))
## mu <- c(1,runif(p.full))
## beta.full <- runif(p.full+1)/5
## beta.reduced <- probit.coef.reduced(mu[-1],Sigma[-1,-1],beta.full,p.reduced)
## beta.combined <- c(beta.full,beta.reduced)
## a <- runif(p.full+1)
## b <- runif(p.reduced+1)
## ## ns <- round(seq(1e2,1e4,len=30))
## ## by.n <- sapply(ns, function(n) {
## ## infl.probit <- function(x,g,beta) { 
## ##     eta <- as.numeric(x%*%beta)
## ##     scores <- dnorm(eta)*(g/pnorm(eta)/(1-pnorm(eta)) - 1/(1-pnorm(eta))) * x
## ##     W <- diag(dnorm(eta)^2/pnorm(eta)/(1-pnorm(eta)))
## ##     hessian <- -t(x)%*%W%*%x / n
## ##     bread <- -solve(hessian)
## ##     infl <- bread%*%t(scores)
## ## }
## auc.full <- auc.probit(mu,Sigma,beta.full)
## auc.reduced <- auc.probit(mu[1:(p.reduced+1)],Sigma[1:(p.reduced+1),1:(p.reduced+1)],beta.reduced)
## obs <- replicate(3e2, {
##     x <- rmvnorm(n,mu,Sigma)
##     x.reduced <- x[,1:(p.reduced+1)]
##     risk <- pnorm(x%*%beta.full)
##     ## browser()
##     g <- rbinom(n,1,risk)
##     glm.full <- glm(g~x-1, family=binomial(link='probit'))
##     glm.reduced <- glm(g~x.reduced-1, family=binomial(link='probit'))
##     beta.hat.full <- coef(glm.full)
##     beta.hat.reduced <- coef(glm.reduced)
##     ## auc.beta.hat <- auc.scores(beta.hat,mu.1-mu.0,2*Sigma)
##     ## auc.gamma.hat <- auc.scores(gamma.hat,mu.1-mu.0,2*Sigma)
##     ## diff.coef.hat <- auc.beta.hat-auc.gamma.hat
##     ## diff.coef <- auc.scores(beta,mu.1-mu.0,2*Sigma) - auc.scores(gamma,mu.1-mu.0,2*Sigma)
##     diff.coef <- auc.full-auc.reduced
##     auc.full.hat <- auc.probit(mu,Sigma,beta.hat.full)
##     auc.reduced.hat <- auc.probit(mu[1:(p.reduced+1)],Sigma[1:(p.reduced+1),1:(p.reduced+1)],beta.hat.reduced)
##     diff.coef.hat <- auc.full.hat-auc.reduced.hat
##     diff.coef - diff.coef.hat
## })


## x <- rmvnorm(n,mu,Sigma)
## x.reduced <- x[,1:p.reduced]
## risk <- pnorm(cbind(1,x)%*%beta.full)
## g <- rbinom(n,1,risk)
## glm.full <- glm(g~x, family=binomial(link='probit'))
## infl.full <- infl.hat.probit(glm.full)
## ## cov.root.full <- expm::sqrtm(cov(t(infl.full))/n)
## ## cov.root.full <- expm::sqrtm(vcov(glm.full))
## ## solve(cov.root.full)%*%(coef(glm.full) - beta.full)
## glm.reduced <- glm(g~x.reduced, family=binomial(link='probit'))
## infl.reduced <- infl.hat.probit(glm.reduced)
## ## cov.root.reduced <- expm::sqrtm(vcov(glm.reduced))
## ## cov.root.reduced <- expm::sqrtm(cov(t(infl.reduced))/n)
## ## solve(cov.root.reduced)%*%(coef(glm.reduced) - beta.reduced)
## infl.combined <- rbind(infl.full,infl.reduced)



## lda case--taylor error is O(1/n). just using beta.hats, no influence function computed yet.
require(mvtnorm)
source('misc.R')
set.seed(1)
n <- 1.5e2
## p.x <- 3
## p.y <- 3
p.full <- 6
p.reduced <- 3
Sigma <- matrix(runif((p.full-1)^2),nrow=p.full-1)
Sigma <- Sigma%*%t(Sigma)
Sigma <- rbind(0,cbind(0,Sigma))
Sigma.reduced <- Sigma[1:p.reduced,1:p.reduced]
mu.0 <- c(1,rep(0,p.full-1))
mu.1 <- c(1,runif(p.full-1))
mu.diff <- mu.1-mu.0
mu.0.reduced <- mu.0[1:p.reduced]; mu.1.reduced <- mu.1[1:p.reduced]
mu.diff.reduced <- mu.1.reduced - mu.0.reduced
pi.0 <- 1/3; pi.1 <- 1-pi.0
beta.full <- c(log(pi.1/pi.0) - 1/2*(t(mu.1[-1])%*%solve(Sigma[-1,-1])%*%mu.1[-1]-t(mu.0[-1])%*%solve(Sigma[-1,-1])%*%mu.0[-1] ) , mu.diff[-1]%*%solve(Sigma[-1,-1])) # could parametrize lda model by beta rather than mu
beta.reduced <- c(log(pi.1/pi.0) - 1/2*(t(mu.1.reduced[-1])%*%solve(Sigma.reduced[-1,-1])%*%mu.1.reduced[-1]-t(mu.0.reduced[-1])%*%solve(Sigma.reduced[-1,-1])%*%mu.0.reduced[-1] ) , mu.diff.reduced[-1]%*%solve(Sigma.reduced[-1,-1]))
beta.combined <- c(beta.full,beta.reduced)
auc.full <- auc.lda(mu.diff,Sigma,beta.full)
auc.reduced <- auc.lda(mu.diff.reduced,Sigma.reduced,beta.reduced)
diff.coef <- auc.full-auc.reduced
ns <- round(seq(1.5e2,1e3,len=20))
by.n <- sapply(ns, function(n) {
    diffs <- replicate(3e2, {
        n.0 <- rbinom(1,n,pi.0)
        n.1 <- n-n.0    
        x.0 <- rmvnorm(n.0,mu.0,Sigma)
        x.1 <- rmvnorm(n.1,mu.1,Sigma)
        x <- rbind(x.0,x.1)
        g <- c(rep(0,n.0),rep(1,n.1))
        x.0.reduced <- x.0[,1:p.reduced]
        x.1.reduced <- x.1[,1:p.reduced]
        beta.hat.full <- lda.coefs(x.0,x.1)
        beta.hat.reduced <- lda.coefs(x.0.reduced,x.1.reduced)
        auc.full.hat <- auc.lda(mu.diff,Sigma,beta.hat.full)
        auc.reduced.hat <- auc.lda(mu.diff.reduced,Sigma.reduced,beta.hat.reduced)
        diff.coef.hat <- auc.full.hat-auc.reduced.hat
        ## obs <- auc.full.hat - auc.full
        ## try <- (beta.hat.full - beta.full)%*%auc.lda.deriv(beta.full,mu.diff,Sigma)
        ## obs <- auc.reduced.hat - auc.reduced
        ## try <- (beta.hat.reduced - beta.reduced)%*%auc.lda.deriv(beta.reduced,mu.diff.reduced,Sigma.reduced)
        ## obs <- diff.coef.hat - diff.coef
        ## try <- (beta.hat.full - beta.full)%*%auc.lda.deriv(beta.full,mu.diff,Sigma) - (beta.hat.reduced - beta.reduced)%*%auc.lda.deriv(beta.hat.reduced,mu.diff.reduced,Sigma.reduced)
        obs <- diff.coef.hat - diff.coef
        try <- (beta.hat.full - beta.full)%*%auc.lda.deriv(beta.full,mu.diff,Sigma) - (beta.hat.reduced - beta.reduced)%*%auc.lda.deriv(beta.hat.reduced,mu.diff.reduced,Sigma.reduced)
        try - obs
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)



## probit case--taylor error is O(1/n). #1 just full auc error, #2
## just reduced auc error, #3 diff.coef.hat - diff.coef, #4 also using
## hatted coefs in the derivative, #5 also using influence function
## for coef.hat-coef. #6 when further using estimated coefs in the
## infl functions, taylor error is O(1/sqrt(n)).
start <- Sys.time()
require(mvtnorm)
source('misc.R')
set.seed(1)
n <- 1.5e2
## p.x <- 3
## p.y <- 3
p.full <- 6
p.reduced <- 3
Sigma <- matrix(runif((p.full-1)^2),nrow=p.full-1)
Sigma <- Sigma%*%t(Sigma)
Sigma <- rbind(0,cbind(0,Sigma))
Sigma.reduced <- Sigma[1:p.reduced,1:p.reduced]
## mu.0 <- c(1,rep(0,p.full-1))
## mu.1 <- c(1,runif(p.full-1))
## mu.diff <- mu.1-mu.0
## mu.0.reduced <- mu.0[1:p.reduced]; mu.1.reduced <- mu.1[1:p.reduced]
## mu.diff.reduced <- mu.1.reduced - mu.0.reduced
## pi.0 <- 1/3; pi.1 <- 1-pi.0
mu <- c(1,runif(p.full-1))
mu.reduced <- mu[1:p.reduced]
beta.full <- runif(p.full)/5
beta.reduced <- probit.coef.reduced(mu[-1],Sigma[-1,-1],beta.full,p.reduced-1)
beta.combined <- c(beta.full,beta.reduced)
auc.full <- auc.probit(mu,Sigma,beta.full)
## auc.reduced <- auc.probit(mu[1:(p.reduced+1)],Sigma[1:(p.reduced+1),1:(p.reduced+1)],beta.reduced)
auc.reduced <- auc.probit(mu[1:(p.reduced)],Sigma[1:(p.reduced),1:(p.reduced)],beta.reduced)
diff.coef <- auc.full-auc.reduced
ns <- round(seq(1.5e2,1e3,len=20))
## deriv.0 <- auc.probit.deriv(beta.full,mu,Sigma)
by.n <- sapply(ns, function(n) {
    cat('.')
    diffs <- replicate(1e2, {
        x <- rmvnorm(n,mu,Sigma)
        x.reduced <- x[,1:p.reduced]
        risk <- pnorm(x%*%beta.full)
        g <- rbinom(n,1,risk)
        glm.full <- glm(g~x-1, family=binomial(link='probit'))
        glm.reduced <- glm(g~x.reduced-1, family=binomial(link='probit'))
        beta.hat.full <- coef(glm.full)
        beta.hat.reduced <- coef(glm.reduced)
        auc.full.hat <- auc.probit(mu,Sigma,beta.hat.full)
        auc.reduced.hat <- auc.probit(mu.reduced,Sigma.reduced,beta.hat.reduced)
        diff.coef.hat <- auc.full.hat-auc.reduced.hat
        ## obs <- auc.full.hat - auc.full #1
        ## try <- t(beta.hat.full - beta.full)%*%auc.probit.deriv(beta.full,mu,Sigma) #1
        ## obs <- auc.reduced.hat - auc.reduced #2
        ## try <- t(beta.hat.reduced - beta.reduced)%*%auc.probit.deriv(beta.reduced,mu.reduced,Sigma.reduced) #2
        ## obs <- diff.coef.hat - diff.coef #3
        ## try <- t(beta.hat.full - beta.full)%*%auc.probit.deriv(beta.full,mu,Sigma) - (beta.hat.reduced - beta.reduced)%*%auc.probit.deriv(beta.reduced,mu.reduced,Sigma.reduced) #3
        ## obs <- diff.coef.hat - diff.coef #4
        ## try <- t(beta.hat.full - beta.full)%*%auc.probit.deriv(beta.hat.full,mu,Sigma) - (beta.hat.reduced - beta.reduced)%*%auc.probit.deriv(beta.hat.reduced,mu.reduced,Sigma.reduced) #4
        ## browser()
        obs <- diff.coef.hat - diff.coef #5
        try <- t(rowMeans(infl.probit(x,g,beta.full)))%*%auc.probit.deriv(beta.hat.full,mu,Sigma) - t(rowMeans(infl.probit(x.reduced,g,beta.reduced)))%*%auc.probit.deriv(beta.hat.reduced,mu.reduced,Sigma.reduced) #5
        ## obs <- diff.coef.hat - diff.coef #6
        ## try <- t(rowMeans(infl.hat.probit(glm.full)))%*%auc.probit.deriv(beta.hat.full,mu,Sigma) - t(rowMeans(infl.hat.probit(glm.reduced)))%*%auc.probit.deriv(beta.hat.reduced,mu.reduced,Sigma.reduced) #6
        try - obs
    })
})
Sys.time()-start
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)

## debugging
## p <- 4
## mu <- runif(p)
## Sigma <- matrix(runif(p^2),p)
## Sigma <- Sigma%*%t(Sigma)
## f <- function(beta)auc.probit(t(beta)%*%mu,t(beta)%*%Sigma%*%beta)
## f.prime <- function(beta)auc.probit.prime(beta,mu,Sigma)
## ## viz.deriv(f,f.prime,p,1)
## beta0 <- runif(p)
## ts <- seq(0.001,1,len=20)
## delta <- runif(p)
## deriv0 <- as.numeric(f.prime(beta0))
## error <- sapply(ts, function(t) f(beta0+t*delta) - f(beta0) - t*(delta%*%deriv0))
## lm(log(abs(error))~log(ts))







## 12g   check hajek projection in probit case is o(n^(-1/2))
require(mvtnorm)
start <- Sys.time()
source('misc.R')
set.seed(3)
## op <- options(warn=-1)
## auc.hajek.probit <- function(x.0,x.1,mu,Sigma,beta,bw=7) {    
##     mu <- t(beta)%*%mu
##     sigma <- sqrt(t(beta)%*%Sigma%*%beta)
##     f.0 <- function(w)(1-pnorm(w))*dnorm((w-mu)/sigma)/sigma/(1-pnorm(mu/sqrt(1+sigma^2)))
##     f.1 <- function(w)pnorm(w)*dnorm((w-mu)/sigma)/sigma/pnorm(mu/sqrt(1+sigma^2))
##     F.0 <- Vectorize(function(w)integrate(f.0,min(w,mu-bw*sigma),w)$val)
##     F.1 <- Vectorize(function(w)integrate(f.1,min(w,mu-bw*sigma),w)$val)
##     ## auc <- integrate(function(x)F.1(x)*f.0(x),mu-bw*sigma,mu+bw*sigma)$val
##     auc <- auc.probit(mu,sigma^2)
##     auc.hajek(F.0,F.1,x.0%*%beta,x.1%*%beta,auc,terms.only=FALSE)
## }
p <- 4
## n <- 1e3
ns <- round(seq(1e2,3e2,len=2e1))
beta <- runif(p)/4
Sigma <- matrix(runif(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
mu <- runif(p)
auc <- auc.probit(mu,Sigma,beta)
by.n <- sapply(ns, function(n) {
    cat('.')
    stats <- replicate(1e2, {
        x <- rmvnorm(n,mu,Sigma)
        ## x.reduced <- x[,1:p.reduced]
        risk <- pnorm(x%*%beta)
        g <- rbinom(n,1,risk)
        x.0 <- x[g==0,]; x.1 <- x[g==1,]
        obs <- auc.hat(x.0%*%beta,x.1%*%beta) - auc
        try <- auc.hajek.probit(x.0,x.1,mu,Sigma,beta)
        obs-try
    })
    ## hist(stats)
    ## abline(v=auc,col=2)
## abline(v=mean(stats),col=3) 
})
Sys.time() - start
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)




## ## 12h true variance of hajek term. abandoned when i realized i need
## ## the joint distribution of beta.full%*%x.0 |D=i and
## ## beta.reduced%*%x.0.reduced |D=i.
## mu.a <- runif(1)
## mu.b <- runif(1)
## sigma.a <- runif(1)
## sigma.b <- runif(1)
## bw <- 7
## ## mu <- mu.a; sigma <- sigma.a
## ab <- lapply(list(a=list(mu=mu.a,sigma=sigma.a),b=list(mu=mu.b,sigma=sigma.b)), function(params) {
##     with(params, {
##         f.0 <- function(w)(1-pnorm(w))*dnorm((w-mu)/sigma)/sigma/(1-pnorm(mu/sqrt(1+sigma^2)))
##         f.1 <- function(w)pnorm(w)*dnorm((w-mu)/sigma)/sigma/pnorm(mu/sqrt(1+sigma^2))
##         F.0 <- Vectorize(function(w)integrate(f.0,min(w,mu-bw*sigma),w)$val)
##         F.1 <- Vectorize(function(w)integrate(f.1,min(w,mu-bw*sigma),w)$val)
##         ## auc <- integrate(function(x)F.1(x)*f.0(x),mu-bw*sigma,mu+bw*sigma)$val
##         auc <- auc.probit(mu,sigma^2)
##         list(mu=mu, sigma=sigma, auc=auc,f.0=f.0,f.1=f.1,F.0=F.0,F.1=F.1)
##     })
## })
## a <- ab$a; b <- ab$b
## ## x.b0 <- 1
## inner.int <- Vectorize(function(x.b0)integrate(function(x.a0)(a$F.1(x.a0)-b$F.1(x.b0)-(a$auc-b$auc))^2*a$f.0(x.a0)*b$f.0(x.b0), a$mu-bw*a$sigma,a$mu+bw*a$sigma)$val)
## outer.int <- integrate(inner.int, b$mu-bw*a$sigma,b$mu+bw*b$sigma)$val
## control.var <- outer.int
## inner.int <- Vectorize(function(x.b1)integrate(function(x.a1)(a$F.0(x.a1)-b$F.0(x.b1)-(a$auc-b$auc))^2*a$f.1(x.a1)*b$f.1(x.b1), a$mu-bw*a$sigma,a$mu+bw*a$sigma)$val)
## outer.int <- integrate(inner.int, b$mu-bw*a$sigma,b$mu+bw*b$sigma)$val
## case.var <- outer.int

## n <- 1e2
## stats <- replicate(1e3, {
##     x.a <- rnorm(n,mu.a,sigma.a)
##     risk.a <- pnorm(x.a)
##     g.a <- rbinom(n,1,risk.a)
##     x.b <- rnorm(n,mu.b,sigma.b)
##     risk.b <- pnorm(x.b)
##     g.b <- rbinom(n,1,risk.b)
##     x.a0 <- x.a[g.a==0]; x.a1 <- x.a[g.a==1]
##     x.b0 <- x.b[g.b==0]; x.b1 <- x.b[g.b==1]
##     mean((a$F.1(x.a0)-b$F.1(x.b0)-(a$auc-b$auc))^2)
##     ## auc.hat(x.a0,x.a1)
## })
## hist(stats)
## abline(v=control.var,col=2)
## abline(v=mean(stats),col=3)



## inner.int <- Vectorize(function(x.b0)integrate(function(x.a0)(a$F.1(x.a0)-b$F.1(x.b0)-(a$auc-b$auc))^2*a$f.0(x.a0)*b$f.0(x.b0), a$mu-bw*a$sigma,a$mu+bw*a$sigma)$val)
## outer.int <- integrate(inner.int, b$mu-bw*a$sigma,b$mu+bw*b$sigma)$val
## control.var <- outer.int
## inner.int <- Vectorize(function(x.b1)integrate(function(x.a1)(a$F.0(x.a1)-b$F.0(x.b1)-(a$auc-b$auc))^2*a$f.1(x.a1)*b$f.1(x.b1), a$mu-bw*a$sigma,a$mu+bw*a$sigma)$val)
## outer.int <- integrate(inner.int, b$mu-bw*a$sigma,b$mu+bw*b$sigma)$val
## case.var <- outer.int


## require(mvtnorm)
## source('misc.R')
## set.seed(1)
## n <- 1.5e2
## ## p.x <- 3
## ## p.y <- 3
## p.full <- 6
## p.reduced <- 3
## Sigma.full <- matrix(runif((p.full-1)^2),nrow=p.full-1)
## Sigma.full <- Sigma.full%*%t(Sigma.full)
## Sigma.full <- rbind(0,cbind(0,Sigma.full))
## Sigma.reduced <- Sigma.full[1:p.reduced,1:p.reduced]
## mu.full <- c(1,runif(p.full-1))
## mu.reduced <- mu.full[1:p.reduced]
## beta.full <- runif(p.full)/5
## beta.reduced <- probit.coef.reduced(mu.full[-1],Sigma.full[-1,-1],beta.full,p.reduced-1)
## ## beta.combined <- c(beta.full,beta.reduced)
## auc.full <- auc.probit(mu.full,Sigma.full,beta.full)
## ## auc.reduced <- auc.probit(mu[1:(p.reduced+1)],Sigma[1:(p.reduced+1),1:(p.reduced+1)],beta.reduced)
## auc.reduced <- auc.probit(mu.reduced,Sigma.reduced,beta.reduced)
## ## diff.coef <- auc.full-auc.reduced
## ns <- round(seq(1.5e2,1e3,len=20))
## ## deriv.0 <- auc.probit.deriv(beta.full,mu,Sigma)
## ## by.n <- sapply(ns, function(n) {
## ##     cat('.')
## ##     diffs <- replicate(1e2, {
## hajeks <- replicate(1e2, {
##     x <- rmvnorm(n,mu.full,Sigma.full)
##     x.reduced <- x[,1:p.reduced]
##     risk <- pnorm(x%*%beta.full)
##     g <- rbinom(n,1,risk)
##     x.0 <- x[g==0,]; x.1 <- x[g==1,]
##     x.0.reduced <- x.reduced[g==0,]; x.1.reduced <- x.reduced[g==1,]
##     ## glm.full <- glm(g~x-1, family=binomial(link='probit'))
##     hajek.full <- auc.hajek.probit(x.0,x.1,mu.full,Sigma.full,beta.full,terms.only=TRUE)
##     hajek.reduced <- auc.hajek.probit(x.0.reduced,x.1.reduced,mu.reduced,Sigma.reduced,beta.reduced,terms.only=TRUE)
##     n.0 <- nrow(x.0); n.1 <- nrow(x.1)
##     sqrt(n.0)*mean(hajek.full$control - hajek.reduced$control)
## })
## var(hajeks)

## bw <- 7
## mu.a <- mu.full%*%beta.full; sigma.a <- sqrt(t(beta.full)%*%Sigma.full%*%beta.full)
## mu.b <- mu.reduced%*%beta.reduced; sigma.b <- sqrt(t(beta.reduced)%*%Sigma.reduced%*%beta.reduced)
## ab <- lapply(list(a=list(mu=mu.a,sigma=sigma.a),b=list(mu=mu.b,sigma=sigma.b)), function(params) {
##     with(params, {
##         f.0 <- function(w)(1-pnorm(w))*dnorm((w-mu)/sigma)/sigma/(1-pnorm(mu/sqrt(1+sigma^2)))
##         f.1 <- function(w)pnorm(w)*dnorm((w-mu)/sigma)/sigma/pnorm(mu/sqrt(1+sigma^2))
##         F.0 <- Vectorize(function(w)integrate(f.0,min(w,mu-bw*sigma),w)$val)
##         F.1 <- Vectorize(function(w)integrate(f.1,min(w,mu-bw*sigma),w)$val)
##         ## auc <- integrate(function(x)F.1(x)*f.0(x),mu-bw*sigma,mu+bw*sigma)$val
##         auc <- auc.probit(mu,sigma^2)
##         list(mu=mu, sigma=sigma, auc=auc,f.0=f.0,f.1=f.1,F.0=F.0,F.1=F.1)
##     })
## })
## a <- ab$a; b <- ab$b
## inner.int <- Vectorize(function(x.b0)integrate(function(x.a0)(a$F.1(x.a0)-b$F.1(x.b0)-(a$auc-b$auc))^2*a$f.0(x.a0)*b$f.0(x.b0), a$mu-bw*a$sigma,a$mu+bw*a$sigma)$val)
## outer.int <- integrate(inner.int, b$mu-bw*a$sigma,b$mu+bw*b$sigma)$val
## control.var <- outer.int
## inner.int <- Vectorize(function(x.b1)integrate(function(x.a1)(a$F.0(x.a1)-b$F.0(x.b1)-(a$auc-b$auc))^2*a$f.1(x.a1)*b$f.1(x.b1), a$mu-bw*a$sigma,a$mu+bw*a$sigma)$val)
## outer.int <- integrate(inner.int, b$mu-bw*a$sigma,b$mu+bw*b$sigma)$val
## case.var <- outer.int
## control.var




## ## debugging

## require(mvtnorm)
## source('misc.R')
## set.seed(1)
## n <- 1e2
## ## p.x <- 3
## ## p.y <- 3
## p.full <- 6
## p.reduced <- 3
## Sigma.full <- matrix(runif((p.full-1)^2),nrow=p.full-1)
## Sigma.full <- Sigma.full%*%t(Sigma.full)
## Sigma.full <- rbind(0,cbind(0,Sigma.full))
## Sigma.reduced <- Sigma.full[1:p.reduced,1:p.reduced]
## mu.full <- c(1,runif(p.full-1))
## mu.reduced <- mu.full[1:p.reduced]
## beta.full <- runif(p.full)/5
## beta.reduced <- probit.coef.reduced(mu.full[-1],Sigma.full[-1,-1],beta.full,p.reduced-1)
## ## beta.combined <- c(beta.full,beta.reduced)
## auc.full <- auc.probit(mu.full,Sigma.full,beta.full)
## ## auc.reduced <- auc.probit(mu[1:(p.reduced+1)],Sigma[1:(p.reduced+1),1:(p.reduced+1)],beta.reduced)
## auc.reduced <- auc.probit(mu.reduced,Sigma.reduced,beta.reduced)
## ## diff.coef <- auc.full-auc.reduced
## ns <- round(seq(1.5e2,1e3,len=20))
## ## deriv.0 <- auc.probit.deriv(beta.full,mu,Sigma)
## ## by.n <- sapply(ns, function(n) {
## ##     cat('.')
## ##     diffs <- replicate(1e2, {
## ## hajeks <- replicate(1e2, {
## stats <- replicate(1e2, {
## x <- rmvnorm(n,mu.full,Sigma.full)
## x.reduced <- x[,1:p.reduced]
## risk <- pnorm(x%*%beta.full)
## g <- rbinom(n,1,risk)
## x.0 <- x[g==0,]; x.1 <- x[g==1,]
## x.0.reduced <- x.reduced[g==0,]; x.1.reduced <- x.reduced[g==1,]
## ## glm.full <- glm(g~x-1, family=binomial(link='probit'))
## ## stat <- auc.hajek.probit(x.0,x.1,mu.full,Sigma.full,beta.full,terms.only=TRUE)$control # F.1(x.0)    
## ## stat <- auc.hajek.probit(x.0,x.1,mu.full,Sigma.full,beta.full,terms.only=TRUE)$case # F.0(x.1)    
## ## stat <- auc.hajek.probit(x.0.reduced,x.1.reduced,mu.reduced,Sigma.reduced,beta.reduced,terms.only=TRUE)$case # F.0(x.1)    
## ## stat <- auc.hajek.probit(x.0,x.1,mu.full,Sigma.full,beta.full,terms.only=TRUE)$control # F.1(x.0)
## hajek.full <- auc.hajek.probit(x.0,x.1,mu.full,Sigma.full,beta.full,terms.only=TRUE)
## hajek.reduced <- auc.hajek.probit(x.0.reduced,x.1.reduced,mu.reduced,Sigma.reduced,beta.reduced,terms.only=TRUE)
## stat <- hajek.full$control 
## stat <- hajek.reduced$control
## stat <- hajek.full$control - hajek.reduced$control
## ## hajek.reduced <- auc.hajek.probit(x.0.reduced,x.1.reduced,mu.reduced,Sigma.reduced,beta.reduced,terms.only=TRUE)
## ## n.0 <- nrow(x.0); n.1 <- nrow(x.1)
## ## sqrt(n.0)*mean(hajek.full$control - hajek.reduced$control)
## ## })
## ## var(hajeks)
## mean(stat^2)
## })
## ## with(list(mu=mu.full,Sigma=Sigma.full,beta=beta.full,auc=auc.full), {
## ##     mu <- t(beta)%*%mu
## ##     sigma <- sqrt(t(beta)%*%Sigma%*%beta)
## ##     f.0 <- function(w)(1-pnorm(w))*dnorm((w-mu)/sigma)/sigma/(1-pnorm(mu/sqrt(1+sigma^2)))
## ##     f.1 <- function(w)pnorm(w)*dnorm((w-mu)/sigma)/sigma/pnorm(mu/sqrt(1+sigma^2))
## ##     F.0 <- Vectorize(function(w)integrate(f.0,min(w,mu-bw*sigma),w)$val)
## ##     F.1 <- Vectorize(function(w)integrate(f.1,min(w,mu-bw*sigma),w)$val)
## ##     ## auc <- integrate(function(x)F.1(x)*f.0(x),mu-bw*sigma,mu+bw*sigma)$val
## ##     integrate(function(x.0)(F.1(x.0)-(1-auc))^2*f.0(x.0),mu-bw*sigma,mu+bw*sigma)$val
## ##     ## integrate(function(x.1)auc*f.1(x.1),mu-bw*sigma,mu+bw*sigma)$val
## ## })
## bw <- 7
## ab <- lapply(list(a=list(mu=mu.full,Sigma=Sigma.full,beta=beta.full),b=list(mu=mu.reduced,Sigma=Sigma.reduced,beta=beta.reduced)), function(params) {
##     with(params, {
##         mu <- t(beta)%*%mu
##         sigma <- sqrt(t(beta)%*%Sigma%*%beta)
##         f.0 <- function(w)(1-pnorm(w))*dnorm((w-mu)/sigma)/sigma/(1-pnorm(mu/sqrt(1+sigma^2)))
##         f.1 <- function(w)pnorm(w)*dnorm((w-mu)/sigma)/sigma/pnorm(mu/sqrt(1+sigma^2))
##         F.0 <- Vectorize(function(w)integrate(f.0,min(w,mu-bw*sigma),w)$val)
##         F.1 <- Vectorize(function(w)integrate(f.1,min(w,mu-bw*sigma),w)$val)
##         ## auc <- integrate(function(x)F.1(x)*f.0(x),mu-bw*sigma,mu+bw*sigma)$val
##         auc <- auc.probit(mu,sigma^2)
##         list(mu=mu, sigma=sigma, auc=auc,f.0=f.0,f.1=f.1,F.0=F.0,F.1=F.1)
##     })
## })
## a <- ab$a; b <- ab$b
## hist(stats)
## ## try <- with(a,    integrate(function(x.0)(F.1(x.0)-(1-auc))^2*f.0(x.0),mu-bw*sigma,mu+bw*sigma)$val)
## ## try <-  integrate(function(x.a0)(a$F.1(x.a0)-(1-a$auc))^2*a$f.0(x.a0),a$mu-bw*a$sigma,a$mu+bw*a$sigma)$val
## ## try <-  integrate(function(x.b0)(b$F.1(x.b0)-(1-b$auc))^2*b$f.0(x.b0),b$mu-bw*b$sigma,b$mu+bw*b$sigma)$val

## inner.int <- Vectorize(function(x.b0)integrate(function(x.a0)(a$F.1(x.a0)-(1-a$auc)- b$F.1(x.b0)+(1-b$auc))^2*a$f.0(x.a0),-5,5)$val)
## try <- integrate(function(x.b0)inner.int(x.b0)*b$f.0(x.b0),-5,5)$val
## try
## abline(v=try,col=2)
## abline(v=mean(stats),col=3)


## 12h   check in probit case that diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef is o(n^(-1/2))

require(mvtnorm)
source('misc.R')
set.seed(1)
start <- Sys.time()
n <- 1.5e2
## p.x <- 3
## p.y <- 3
p.full <- 6
p.reduced <- 3
Sigma.full <- matrix(runif((p.full-1)^2),nrow=p.full-1)
Sigma.full <- Sigma.full%*%t(Sigma.full)
Sigma.full <- rbind(0,cbind(0,Sigma.full))
Sigma.reduced <- Sigma.full[1:p.reduced,1:p.reduced]
mu.full <- c(1,runif(p.full-1))
mu.reduced <- mu.full[1:p.reduced]
beta.full <- runif(p.full)/5
beta.reduced <- probit.coef.reduced(mu.full[-1],Sigma.full[-1,-1],beta.full,p.reduced-1)
## beta.combined <- c(beta.full,beta.reduced)
## auc.full <- auc.probit(mu.full,Sigma.full,beta.full)
## ## auc.reduced <- auc.probit(mu[1:(p.reduced+1)],Sigma[1:(p.reduced+1),1:(p.reduced+1)],beta.reduced)
## auc.reduced <- auc.probit(mu.reduced,Sigma.reduced,beta.reduced)
ns <- round(seq(1.5e2,3e3,len=20))
by.n <- sapply(ns, function(n) {
    diffs <-  replicate(10, {
        x <- rmvnorm(n,mu.full,Sigma.full)
        x.reduced <- x[,1:p.reduced]
        risk <- pnorm(x%*%beta.full)
        g <- rbinom(n,1,risk)
        x.0 <- x[g==0,]; x.1 <- x[g==1,]
        x.0.reduced <- x.reduced[g==0,]; x.1.reduced <- x.reduced[g==1,]
        glm.full <- glm(g~x-1, family=binomial(link='probit'))
        glm.reduced <- glm(g~x.reduced-1, family=binomial(link='probit'))
        beta.hat.full <- coef(glm.full)
        beta.hat.reduced <- coef(glm.reduced)
        ## n.0 <- nrow(x.0); n.1 <- nrow(x.1)
        ## sqrt(n.0)*mean(hajek.full$control - hajek.reduced$control)
        auc.hat.full.hat <- auc.hat(x.0%*%beta.full,x.1%*%beta.full)
        auc.hat.reduced.hat <- auc.hat(x.0.reduced%*%beta.reduced,x.1.reduced%*%beta.reduced)
        diff.hat.coef.hat <- auc.hat.full.hat - auc.hat.reduced.hat
        auc.full.hat <- auc.probit(mu.full,Sigma.full,beta.hat.full)
        ## auc.reduced <- auc.probit(mu[1:(p.reduced+1)],Sigma[1:(p.reduced+1),1:(p.reduced+1)],beta.reduced)
        auc.reduced.hat <- auc.probit(mu.reduced,Sigma.reduced,beta.hat.reduced)
        diff.coef.hat <- auc.full.hat - auc.reduced.hat
        auc.hajek.full <- auc.hajek.probit(x.0,x.1,mu.full,Sigma.full,beta.hat.full,terms.only=FALSE)
        auc.hajek.reduced <- auc.hajek.probit(x.0.reduced,x.1.reduced,mu.reduced,Sigma.reduced,beta.hat.reduced,terms.only=FALSE)
        diff.hajek.coef <- auc.hajek.full - auc.hajek.reduced
        diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef
        ## auc.hat.full.hat - auc.full.hat - auc.hajek.full
    })
})        
Sys.time()-start
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)



## 12i hajek as iid sum

## 12i-1 estimating variance after splitting sample. doesn't seem to affect
## the variance of the variance estimate.
ns <- round(seq(10,1e2,len=20))
by.n <- sapply(ns, function(n) {
    pairs <- replicate(1e2, {
        x <- rnorm(n)
        c( var(x)/n,var(x[1:(n/2)])/2/n+var(x[(n/2+1):n])/2/n )
    })
    apply(pairs,1,var)
})
matplot(t(by.n))


## 12i-2 direct computation of delong variance
source('misc.r')
n <- 1e2
n.0 <- sample(n,1)
n.1 <- sample(n,1)
x.a0 <- rnorm(n.0)
x.a1 <- rnorm(n.1)
x.b0 <- rnorm(n.0)
x.b1 <- rnorm(n.1)
hajek.a <- auc.hajek(x=x.a0,y=x.a1)
hajek.b <- auc.hajek(x=x.b0,y=x.b1)
var(hajek.a$control-hajek.b$control)/n.0 + var(hajek.a$case-hajek.b$case)/n.1
delong.var.new(cbind(x.a0,x.b0),cbind(x.a1,x.b1))

## 12i-3 variance of iid sum , for a single auc estimate
source('misc.r')
n <- 1e2
mu <- runif(1)
sigma <- runif(1)
x <- rnorm(n,mu,sigma)
risk <- pnorm(x)
g <- rbinom(n,1,risk)
x.0 <- x[g==0]; x.1 <- x[g==1]
n.0 <- length(x.0); n.1 <- length(x.1)
## auc.hat(x.0,x.1)
hajek <- auc.hajek(x=x.0,y=x.1)
var.hat <- var(hajek$control)/n.0 + var(hajek$case)/n.1
iid <- auc.hajek(x=x.0,y=x.1,IID=TRUE)
var.hat.try <- var(iid)/length(iid)
c(var.hat,var.hat.try)







## 12j #1 combining delong and deriv part approximations. #2 writing the approximation as an IID sum.

require(mvtnorm)
source('misc.R')
set.seed(1)
start <- Sys.time()
## n <- 1.5e2
p.full <- 6
p.reduced <- 3
Sigma.full <- matrix(runif((p.full-1)^2),nrow=p.full-1)
Sigma.full <- Sigma.full%*%t(Sigma.full)
Sigma.full <- rbind(0,cbind(0,Sigma.full))
Sigma.reduced <- Sigma.full[1:p.reduced,1:p.reduced]
mu.full <- c(1,runif(p.full-1))
mu.reduced <- mu.full[1:p.reduced]
beta.full <- runif(p.full)/5
beta.reduced <- probit.coef.reduced(mu.full[-1],Sigma.full[-1,-1],beta.full,p.reduced-1)
auc.full <- auc.probit(mu.full,Sigma.full,beta.full)
auc.reduced <- auc.probit(mu.reduced,Sigma.reduced,beta.reduced)
diff.coef <- auc.full-auc.reduced
ns <- round(seq(1.5e2,3e3,len=20))
by.n <- sapply(ns, function(n) {
    diffs <-  replicate(10, {
        x <- rmvnorm(n,mu.full,Sigma.full)
        x.reduced <- x[,1:p.reduced]
        risk <- pnorm(x%*%beta.full)
        g <- rbinom(n,1,risk)
        x.0 <- x[g==0,]; x.1 <- x[g==1,]
        x.0.reduced <- x.reduced[g==0,]; x.1.reduced <- x.reduced[g==1,]
        glm.full <- glm(g~x-1, family=binomial(link='probit'))
        glm.reduced <- glm(g~x.reduced-1, family=binomial(link='probit'))
        beta.hat.full <- coef(glm.full)
        beta.hat.reduced <- coef(glm.reduced)
        auc.hat.full.hat <- auc.hat(x.0%*%beta.full,x.1%*%beta.full)
        auc.hat.reduced.hat <- auc.hat(x.0.reduced%*%beta.reduced,x.1.reduced%*%beta.reduced)
        diff.hat.coef.hat <- auc.hat.full.hat - auc.hat.reduced.hat
        auc.full.hat <- auc.probit(mu.full,Sigma.full,beta.hat.full)
        auc.reduced.hat <- auc.probit(mu.reduced,Sigma.reduced,beta.hat.reduced)
        diff.coef.hat <- auc.full.hat - auc.reduced.hat
        ## auc.hajek.full <- auc.hajek.probit(x.0,x.1,mu.full,Sigma.full,beta.hat.full,terms.only=FALSE)
        ## auc.hajek.reduced <- auc.hajek.probit(x.0.reduced,x.1.reduced,mu.reduced,Sigma.reduced,beta.hat.reduced,terms.only=FALSE) #1
        ## ## delong.part <- diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef #1
        ## taylor <- t(rowMeans(infl.probit(x,g,beta.full)))%*%auc.probit.deriv(beta.hat.full,mu.full,Sigma.full) - t(rowMeans(infl.probit(x.reduced,g,beta.reduced)))%*%auc.probit.deriv(beta.hat.reduced,mu.reduced,Sigma.reduced) #1
        ## ## deriv.part <- diff.coef.hat - diff.coef - taylor #1
        ## ## delong.part + deriv.part #1
        ## diff.hajek.coef <- auc.hajek.full - auc.hajek.reduced #1
        ## target <- diff.hat.coef.hat - diff.coef #1
        ## try <- diff.hajek.coef + taylor #1
        auc.hajek.full <- auc.hajek.probit(x.0,x.1,mu.full,Sigma.full,beta.hat.full,terms.only=TRUE,IID=TRUE) #2
        auc.hajek.reduced <- auc.hajek.probit(x.0.reduced,x.1.reduced,mu.reduced,Sigma.reduced,beta.hat.reduced,terms.only=TRUE,IID=TRUE) #2
        taylor <-  t(auc.probit.deriv(beta.hat.full,mu.full,Sigma.full))%*%infl.probit(x,g,beta.full) - t(auc.probit.deriv(beta.hat.reduced,mu.reduced,Sigma.reduced))%*%infl.probit(x.reduced,g,beta.reduced) #2
        diff.hajek.coef <- auc.hajek.full - auc.hajek.reduced #2
        target <- diff.hat.coef.hat - diff.coef #2
        try <- mean(diff.hajek.coef + taylor) #2
        target - try
    })
})        
Sys.time()-start
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)




## 12k #1 relative error of variance estimate. vanishes nicely. #2
## using estimated influence function. #3 using estimated hajek
## projection.

require(mvtnorm)
source('misc.R')
set.seed(1)
start <- Sys.time()
## n <- 1.5e2
p.full <- 6
p.reduced <- 3
Sigma.full <- matrix(runif((p.full-1)^2),nrow=p.full-1)
Sigma.full <- Sigma.full%*%t(Sigma.full)
Sigma.full <- rbind(0,cbind(0,Sigma.full))
Sigma.reduced <- Sigma.full[1:p.reduced,1:p.reduced]
mu.full <- c(1,runif(p.full-1))
mu.reduced <- mu.full[1:p.reduced]
beta.full <- runif(p.full)/5
beta.reduced <- probit.coef.reduced(mu.full[-1],Sigma.full[-1,-1],beta.full,p.reduced-1)
auc.full <- auc.probit(mu.full,Sigma.full,beta.full)
auc.reduced <- auc.probit(mu.reduced,Sigma.reduced,beta.reduced)
diff.coef <- auc.full-auc.reduced
ns <- round(seq(1.5e2,2e3,len=20))
by.n <- sapply(ns, function(n) {
    stats <-  replicate(10, {
        x <- rmvnorm(n,mu.full,Sigma.full)
        x.reduced <- x[,1:p.reduced]
        risk <- pnorm(x%*%beta.full)
        g <- rbinom(n,1,risk)
        x.0 <- x[g==0,]; x.1 <- x[g==1,]
        x.0.reduced <- x.reduced[g==0,]; x.1.reduced <- x.reduced[g==1,]
        glm.full <- glm(g~x-1, family=binomial(link='probit'))
        glm.reduced <- glm(g~x.reduced-1, family=binomial(link='probit'))
        beta.hat.full <- coef(glm.full)
        beta.hat.reduced <- coef(glm.reduced)
        auc.hat.full.hat <- auc.hat(x.0%*%beta.full,x.1%*%beta.full)
        auc.hat.reduced.hat <- auc.hat(x.0.reduced%*%beta.reduced,x.1.reduced%*%beta.reduced)
        diff.hat.coef.hat <- auc.hat.full.hat - auc.hat.reduced.hat
        ## auc.full.hat <- auc.probit(mu.full,Sigma.full,beta.hat.full)
        ## auc.reduced.hat <- auc.probit(mu.reduced,Sigma.reduced,beta.hat.reduced)
        ## diff.coef.hat <- auc.full.hat - auc.reduced.hat
        ## auc.hajek.full <- auc.hajek.probit(x.0,x.1,mu.full,Sigma.full,beta.hat.full,terms.only=TRUE,IID=TRUE) #1
        ## auc.hajek.reduced <- auc.hajek.probit(x.0.reduced,x.1.reduced,mu.reduced,Sigma.reduced,beta.hat.reduced,terms.only=TRUE,IID=TRUE) #1
        auc.hajek.full <- auc.hajek(x=x.0%*%beta.hat.full,y=x.1%*%beta.hat.full,terms.only=TRUE,IID=TRUE) #3
        auc.hajek.reduced <- auc.hajek(x=x.0.reduced%*%beta.hat.reduced,y=x.1.reduced%*%beta.hat.reduced,terms.only=TRUE,IID=TRUE) #3
        ## taylor <-  t(auc.probit.deriv(beta.hat.full,mu.full,Sigma.full))%*%infl.probit(x,g,beta.full) - t(auc.probit.deriv(beta.hat.reduced,mu.reduced,Sigma.reduced))%*%infl.probit(x.reduced,g,beta.reduced) #1
        taylor <-  t(auc.probit.deriv(beta.hat.full,mu.full,Sigma.full))%*%infl.probit(x,g,beta.hat.full) - t(auc.probit.deriv(beta.hat.reduced,mu.reduced,Sigma.reduced))%*%infl.probit(x.reduced,g,beta.hat.reduced) #2
        diff.hajek.coef <- auc.hajek.full - auc.hajek.reduced 
        target <- diff.hat.coef.hat - diff.coef
        approx <- as.numeric(diff.hajek.coef + taylor) 
        c(estimate=target, var.hat=var(approx))
    })    
    B <- length(stats['estimate',])
    ## hist(stats['var.hat',])
    ## abline(v=B*var(stats['estimate',]),col=2)
    ## abline(v=mean(stats['var.hat',]),col=3)
    (stats['var.hat',] - var(sqrt(B)*stats['estimate',])) / var(sqrt(B)*stats['estimate',])
    apply(stats,1,var)
})        
Sys.time()-start
matplot(ns,t(by.n),pch=1,col=1)
abline(h=0)


## chekcing z-stats. they're off. issue with auc.probit?
require(mvtnorm)
source('misc.R')
set.seed(1)
start <- Sys.time()
## n <- 1.5e2
p.full <- 6
p.reduced <- 3
Sigma.full <- matrix(runif((p.full-1)^2),nrow=p.full-1)
Sigma.full <- Sigma.full%*%t(Sigma.full)
Sigma.full <- rbind(0,cbind(0,Sigma.full))
Sigma.reduced <- Sigma.full[1:p.reduced,1:p.reduced]
mu.full <- c(1,runif(p.full-1))
mu.reduced <- mu.full[1:p.reduced]
beta.full <- runif(p.full)/5
beta.reduced <- probit.coef.reduced(mu.full[-1],Sigma.full[-1,-1],beta.full,p.reduced-1)
auc.full <- auc.probit(mu.full,Sigma.full,beta.full)
auc.reduced <- auc.probit(mu.reduced,Sigma.reduced,beta.reduced)
diff.coef <- auc.full-auc.reduced
ns <- round(seq(1.5e2,2e3,len=20))
by.n <- sapply(ns, function(n) {

    stats <-  replicate(1e2, {
        x <- rmvnorm(n,mu.full,Sigma.full)
        x.reduced <- x[,1:p.reduced]
        risk <- pnorm(x%*%beta.full)
        g <- rbinom(n,1,risk)
        x.0 <- x[g==0,]; x.1 <- x[g==1,]
        x.0.reduced <- x.reduced[g==0,]; x.1.reduced <- x.reduced[g==1,]
        glm.full <- glm(g~x-1, family=binomial(link='probit'))
        glm.reduced <- glm(g~x.reduced-1, family=binomial(link='probit'))
        beta.hat.full <- coef(glm.full)
        beta.hat.reduced <- coef(glm.reduced)
        auc.hat.full.hat <- auc.hat(x.0%*%beta.full,x.1%*%beta.full)
        auc.hat.reduced.hat <- auc.hat(x.0.reduced%*%beta.reduced,x.1.reduced%*%beta.reduced)
        diff.hat.coef.hat <- auc.hat.full.hat - auc.hat.reduced.hat
        ## auc.full.hat <- auc.probit(mu.full,Sigma.full,beta.hat.full)
        ## auc.reduced.hat <- auc.probit(mu.reduced,Sigma.reduced,beta.hat.reduced)
        ## diff.coef.hat <- auc.full.hat - auc.reduced.hat
        ## auc.hajek.full <- auc.hajek.probit(x.0,x.1,mu.full,Sigma.full,beta.hat.full,terms.only=TRUE,IID=TRUE) #1
        ## auc.hajek.reduced <- auc.hajek.probit(x.0.reduced,x.1.reduced,mu.reduced,Sigma.reduced,beta.hat.reduced,terms.only=TRUE,IID=TRUE) #1
        auc.hajek.full <- auc.hajek(x=x.0%*%beta.hat.full,y=x.1%*%beta.hat.full,terms.only=TRUE,IID=TRUE) #3
        auc.hajek.reduced <- auc.hajek(x=x.0.reduced%*%beta.hat.reduced,y=x.1.reduced%*%beta.hat.reduced,terms.only=TRUE,IID=TRUE) #3
        ## taylor <-  t(auc.probit.deriv(beta.hat.full,mu.full,Sigma.full))%*%infl.probit(x,g,beta.full) - t(auc.probit.deriv(beta.hat.reduced,mu.reduced,Sigma.reduced))%*%infl.probit(x.reduced,g,beta.reduced) #1
        taylor <-  t(auc.probit.deriv(beta.hat.full,mu.full,Sigma.full))%*%infl.probit(x,g,beta.hat.full) - t(auc.probit.deriv(beta.hat.reduced,mu.reduced,Sigma.reduced))%*%infl.probit(x.reduced,g,beta.hat.reduced) #2
        diff.hajek.coef <- auc.hajek.full - auc.hajek.reduced 
        target <- diff.hat.coef.hat - diff.coef
        approx <- as.numeric(diff.hajek.coef + taylor) 
        target / sd(approx)
    })
    qqnorm(stats)
    abline(0,1)

})        
Sys.time()-start
matplot(ns,t(by.n),pch=1,col=1)
abline(h=0)



## 12l fixing auc.probit

## 12l-1 densities of gamma^t x | g=i
require(mvtnorm)
set.seed(1)
p <- 3
n <- 5e3
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
gamma <- runif(p)
x <- rmvnorm(n,mu,Sigma)
risk <- pnorm(x%*%beta)
g <- rbinom(n,1,risk)
hist((x%*%gamma)[g==1],prob=TRUE)
beta.quad <- as.numeric(t(beta)%*%Sigma%*%beta)
gamma.quad <- as.numeric(t(gamma)%*%Sigma%*%gamma)
E.g <- pnorm((beta%*%mu)/sqrt(1+beta.quad))
rho <- (t(gamma)%*%Sigma%*%beta) / sqrt(beta.quad) / sqrt(gamma.quad)
sigma.cond <- (1-rho^2)*gamma.quad
## v <- runif(mu)
f <- Vectorize(function(v)integrate(function(w)dnorm((v-gamma%*%mu-(t(gamma)%*%Sigma%*%beta)/beta.quad*(w-t(beta)%*%mu))/sqrt(sigma.cond))/sqrt(sigma.cond)*pnorm(w)*dnorm((w-beta%*%mu)/sqrt(beta.quad))/sqrt(beta.quad)/E.g,-10,10)$val)
a <- -t(gamma)%*%Sigma%*%beta/beta.quad/sqrt(sigma.cond)
b <- function(v)(v-gamma%*%mu)/sqrt(sigma.cond)+t(gamma)%*%Sigma%*%beta/beta.quad/sqrt(sigma.cond)*(beta%*%mu)
c <- 1/sqrt(beta.quad)
d <- -t(beta)%*%mu/sqrt(beta.quad)
f.try <- Vectorize(function(v)integrate(function(w)dnorm(a*w + b(v))      /sqrt(sigma.cond)*pnorm(w)*dnorm((w-beta%*%mu)/sqrt(beta.quad))/sqrt(beta.quad)/E.g,-10,10)$val)
## f.try <- Vectorize(function(v)1/sqrt(beta.quad)/E.g/sqrt(sigma.cond)*integrate(function(w)pnorm(w)*dnorm(a*w+b(v))*dnorm(c*w+d),-10,10)$val)
## integrate(f,-8,8)
## a <- runif(1);b <- runif(1);c <- runif(1);d <- runif(1);w <- runif(1)
## dnorm(a*w+b)*dnorm(c*w+d)
## 1/sqrt(2*pi)*dnorm(w*sqrt(a^2+c^2) + (a*b+c*d)/sqrt(a^2+c^2))* exp(1/2*(a*b+c*d)^2/(a^2+c^2)-(b^2+d^2)/2)
## f.try <- function(v)1/sqrt(a^2+c^2)*pnorm(-(a*b(v)+c*d)/sqrt(a^2+c^2)/sqrt(1+a^2+c^2))*dnorm(sqrt(b(v)^2+d^2-(a*b(v)+c*d)^2/(a^2+c^2)))/sqrt(beta.quad)/E.g/sqrt(sigma.cond)
## curve(f)
curve(f.try,add=TRUE,col=2)


require(mvtnorm)
set.seed(1)
p <- 3
n <- 5e3
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
gamma <- runif(p)
x <- rmvnorm(n,mu,Sigma)
risk <- pnorm(x%*%beta)
g <- rbinom(n,1,risk)
hist((x%*%gamma)[g==0],prob=TRUE)
beta.quad <- as.numeric(t(beta)%*%Sigma%*%beta)
gamma.quad <- as.numeric(t(gamma)%*%Sigma%*%gamma)
E.g <- pnorm((beta%*%mu)/sqrt(1+beta.quad))
rho <- (t(gamma)%*%Sigma%*%beta) / sqrt(beta.quad) / sqrt(gamma.quad)
sigma.cond <- (1-rho^2)*gamma.quad
a <- -t(gamma)%*%Sigma%*%beta/beta.quad/sqrt(sigma.cond)
b <- function(v)(v-gamma%*%mu)/sqrt(sigma.cond)+t(gamma)%*%Sigma%*%beta/beta.quad/sqrt(sigma.cond)*(beta%*%mu)
c <- 1/sqrt(beta.quad)
d <- -t(beta)%*%mu/sqrt(beta.quad)
## f <- function(x)1/sqrt(sigma.cond)/abs(a) - E.g/(1-E.g)*f.gamma.1(x,mu,Sigma,gamma,beta)
## f <- Vectorize(function(x)integrate(function(w)dnorm(a*w+b(x))/sqrt(sigma.cond)*(1-pnorm(w))*dnorm(c*w+d)/sqrt(beta.quad)/(1-E.g),-10,10))
f <- Vectorize(function(x)1/sqrt(beta.quad)/(1-E.g)/sqrt(sigma.cond)*integrate(function(w)dnorm(a*w+b(x))*dnorm(c*w+d),-10,10)$val - 1/sqrt(beta.quad)/(1-E.g)/sqrt(sigma.cond)*integrate(function(w)pnorm(w)*dnorm(a*w+b(x))*dnorm(c*w+d),-10,10)$val)
## x <- 1
## integrate(function(w)dnorm(a*w+b(x))*dnorm(c*w+d),-10,10)$val
## 1/sqrt(a^2+c^2)/sqrt(2*pi)*exp(1/2*(a*b(x)+c*d)^2/(a^2+c^2)-(b(x)^2+d^2)/2)
f <- Vectorize(function(x)1/sqrt(beta.quad)/(1-E.g)/sqrt(sigma.cond)*1/sqrt(a^2+c^2)/sqrt(2*pi)*exp(1/2*(a*b(x)+c*d)^2/(a^2+c^2)-(b(x)^2+d^2)/2) - E.g/(1-E.g)*f.gamma.1(x,mu,Sigma,gamma,beta))
curve(f(x),col=2,add=TRUE)
curve(f.gamma.1(x,mu,Sigma,gamma,beta),add=TRUE,col=2)



require(mvtnorm)
set.seed(1)
p <- 3
n <- 5e3
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
gamma <- runif(p)
x <- rmvnorm(n,mu,Sigma)
risk <- pnorm(x%*%beta)
g <- rbinom(n,1,risk)
f.gamma.1 <- function(x,mu,Sigma,gamma,beta) {
    beta.quad <- as.numeric(t(beta)%*%Sigma%*%beta)
    gamma.quad <- as.numeric(t(gamma)%*%Sigma%*%gamma)
    E.g <- pnorm((beta%*%mu)/sqrt(1+beta.quad))
    rho <- (t(gamma)%*%Sigma%*%beta) / sqrt(beta.quad) / sqrt(gamma.quad)
    sigma.cond <- (1-rho^2)*gamma.quad
    a <- -t(gamma)%*%Sigma%*%beta/beta.quad/sqrt(sigma.cond)
    b <- function(v)(v-gamma%*%mu)/sqrt(sigma.cond)+t(gamma)%*%Sigma%*%beta/beta.quad/sqrt(sigma.cond)*(beta%*%mu)
    c <- 1/sqrt(beta.quad)
    d <- -t(beta)%*%mu/sqrt(beta.quad)
    1/sqrt(a^2+c^2)*pnorm(-(a*b(x)+c*d)/sqrt(a^2+c^2)/sqrt(1+a^2+c^2))*dnorm(sqrt(b(x)^2+d^2-(a*b(x)+c*d)^2/(a^2+c^2)))/sqrt(beta.quad)/E.g/sqrt(sigma.cond)
}
f.gamma.0 <- function(x,mu,Sigma,gamma,beta) {
    beta.quad <- as.numeric(t(beta)%*%Sigma%*%beta)
    gamma.quad <- as.numeric(t(gamma)%*%Sigma%*%gamma)
    E.g <- pnorm((beta%*%mu)/sqrt(1+beta.quad))
    rho <- (t(gamma)%*%Sigma%*%beta) / sqrt(beta.quad) / sqrt(gamma.quad)
    sigma.cond <- (1-rho^2)*gamma.quad
    a <- -t(gamma)%*%Sigma%*%beta/beta.quad/sqrt(sigma.cond)
    b <- function(v)(v-gamma%*%mu)/sqrt(sigma.cond)+t(gamma)%*%Sigma%*%beta/beta.quad/sqrt(sigma.cond)*(beta%*%mu)
    c <- 1/sqrt(beta.quad)
    d <- -t(beta)%*%mu/sqrt(beta.quad)
    1/sqrt(beta.quad)/(1-E.g)/sqrt(sigma.cond)*1/sqrt(a^2+c^2)/sqrt(2*pi)*exp(1/2*(a*b(x)+c*d)^2/(a^2+c^2)-(b(x)^2+d^2)/2) - E.g/(1-E.g)*f.gamma.1(x,mu,Sigma,gamma,beta)
}
pdf.index.probit <- function(x,g,mu,Sigma,gamma,beta) {
    beta.quad <- as.numeric(t(beta)%*%Sigma%*%beta)
    gamma.quad <- as.numeric(t(gamma)%*%Sigma%*%gamma)
    E.g <- pnorm((beta%*%mu)/sqrt(1+beta.quad))
    rho <- (t(gamma)%*%Sigma%*%beta) / sqrt(beta.quad) / sqrt(gamma.quad)
    sigma.cond <- (1-rho^2)*gamma.quad
    a <- -t(gamma)%*%Sigma%*%beta/beta.quad/sqrt(sigma.cond)
    b <- function(v)(v-gamma%*%mu)/sqrt(sigma.cond)+t(gamma)%*%Sigma%*%beta/beta.quad/sqrt(sigma.cond)*(beta%*%mu)
    c <- 1/sqrt(beta.quad)
    d <- -t(beta)%*%mu/sqrt(beta.quad)
    g.1 <-     1/sqrt(a^2+c^2)*pnorm(-(a*b(x)+c*d)/sqrt(a^2+c^2)/sqrt(1+a^2+c^2))*dnorm(sqrt(b(x)^2+d^2-(a*b(x)+c*d)^2/(a^2+c^2)))/sqrt(beta.quad)/E.g/sqrt(sigma.cond)    
    if(g==0)return ( 1/sqrt(beta.quad)/(1-E.g)/sqrt(sigma.cond)*1/sqrt(a^2+c^2)/sqrt(2*pi)*exp(1/2*(a*b(x)+c*d)^2/(a^2+c^2)-(b(x)^2+d^2)/2) - E.g/(1-E.g)*g.1 )
    if(g==1)return(g.1)
    stop()
}
op <- par(mfrow=c(1,2))
hist((x%*%gamma)[g==1],prob=TRUE)
## curve(f.gamma.1(x,mu,Sigma,gamma,beta),add=TRUE,col=2)
curve(pdf.index.probit(x,1,mu,Sigma,gamma,beta),add=TRUE,col=2)
hist((x%*%gamma)[g==0],prob=TRUE)
## curve(f.gamma.0(x,mu,Sigma,gamma,beta),add=TRUE,col=2)
curve(pdf.index.probit(x,0,mu,Sigma,gamma,beta),add=TRUE,col=2)
par(op)

## 12l-2 auc
require(mvtnorm)
source('misc.R')
set.seed(1)
p <- 3
n <- 5e2
auc.probit <- function(mu,Sigma,beta,gamma,bw=5){
    lb <- gamma%*%mu-bw*t(gamma)%*%Sigma%*%gamma
    ub <- gamma%*%mu+bw*t(gamma)%*%Sigma%*%gamma
    ## inner.int <- Vectorize(function(v)integrate(function(w)f.gamma.0(w,mu,Sigma,gamma,beta),lb,v)$val)
    ## try <- integrate(function(v)inner.int(v)*f.gamma.1(v,mu,Sigma,gamma,beta),round(lb),round(ub))$val
    inner.int <- Vectorize(function(v)integrate(function(w)pdf.index.probit(w,0,mu,Sigma,gamma,beta),lb,v)$val)
    integrate(function(v)inner.int(v)*pdf.index.probit(v,1,mu,Sigma,gamma,beta),round(lb),round(ub))$val
}
pairs <- replicate(3e1,{
    cat('.')
    mu <- runif(p)
    Sigma <- matrix(runif(p^2),p)
    Sigma <- Sigma%*%t(Sigma)
    beta <- runif(p)
    gamma <- runif(p)
    x <- rmvnorm(n,mu,Sigma)
    risk <- pnorm(x%*%beta)
    g <- rbinom(n,1,risk)
    obs <- auc.hat((x%*%gamma)[g==0],(x%*%gamma)[g==1])
    c(obs=obs,try=auc.probit(mu,Sigma,beta,gamma))
})
plot(pairs['obs',],pairs['try',]);abline(0,1)



## 12l-3 auc maximized when gamma==beta
source('misc.R')
set.seed(3)
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
beta <- runif(p)
delta <- runif(p)
ts <- seq(-1/2,1/2,len=20)
aucs <- sapply(ts, function(t)auc.probit(mu,Sigma,beta,beta+t*delta))
plot(ts,aucs)
abline(v=0)


## 12l-4 therefore, with probit data, no need to account for estimated beta.hats

require(mvtnorm)
source('misc.R')
set.seed(1)
## using the old auc.probit that assumes beta==gamma, because the
## "repaired" one is currently crashing when beta==gamma
auc.probit <- function(mean.x,var.x,beta=NULL,bw=8) {
    if(!is.null(beta)) {
        stopifnot(length(beta)==length(mean.x))
        mean.x <- beta%*%mean.x
        var.x <- t(beta)%*%var.x%*%beta
    } else     stopifnot(length(mean.x)==1 && length(var.x)==1)
    E.g <- pnorm(mean.x/sqrt(1+var.x))
    mu <- mean.x; sigma <- sqrt(var.x) # clean up
    f.case <- function(w)pnorm(w)*dnorm((w-c(mu))/c(sigma))/c(sigma*E.g)
    f.control <- function(w)(1-pnorm(w))*dnorm((w-c(mu))/c(sigma))/c(sigma*(1-E.g))
    int.inner <- Vectorize(function(v)integrate(function(u)f.control(u),mu-bw*sigma,v)$val)
    integrate(function(v)f.case(v)*int.inner(v),mu-bw*sigma,mu+bw*sigma)$val
}
delong.var <- function(x=NULL,y=NULL,xy=NULL,g=NULL) {
    if(is.null(x)||is.null(y)) {
        if(is.null(xy)||is.null(g)) stop()
        g <- factor(g,labels=0:1)
        x <- xy[g==0,]; y <- xy[g==1,]
    }
    m <- nrow(x); n <- nrow(y)
    theta.hats <- sapply(1:2, function(i) mean(outer(x[,i],y[,i],'<')))
    F.hat.1 <- ecdf(x[,1])
    F.hat.2 <- ecdf(x[,2])
    G.hat.1 <- ecdf(y[,1])
    G.hat.2 <- ecdf(y[,2])
    hajek.diff.terms <- mapply('-',auc.hajek(F.hat.1,G.hat.1,x[,1],y[,1],0),auc.hajek(F.hat.2,G.hat.2,x[,2],y[,2],0))
    sum(sapply(hajek.diff.terms,var) / c(m,n))
}
## start <- Sys.time()
n <- 2e3
p.full <- 6
p.reduced <- 3
Sigma.full <- matrix(runif((p.full-1)^2),nrow=p.full-1)
Sigma.full <- Sigma.full%*%t(Sigma.full)
Sigma.full <- rbind(0,cbind(0,Sigma.full))
Sigma.reduced <- Sigma.full[1:p.reduced,1:p.reduced]
mu.full <- c(1,runif(p.full-1))
mu.reduced <- mu.full[1:p.reduced]
beta.full <- runif(p.full)/5
beta.reduced <- probit.coef.reduced(mu.full[-1],Sigma.full[-1,-1],beta.full,p.reduced-1)
auc.full <- auc.probit(mu.full,Sigma.full,beta.full)
auc.reduced <- auc.probit(mu.reduced,Sigma.reduced,beta.reduced)
diff.coef <- auc.full-auc.reduced
## ns <- round(seq(1.5e2,2e3,len=20))
## by.n <- sapply(ns, function(n) {
stats <-  replicate(3e2, {
    x <- rmvnorm(n,mu.full,Sigma.full)
    x.reduced <- x[,1:p.reduced]
    risk <- pnorm(x%*%beta.full)
    g <- rbinom(n,1,risk)
    x.0 <- x[g==0,]; x.1 <- x[g==1,]
    x.0.reduced <- x.reduced[g==0,]; x.1.reduced <- x.reduced[g==1,]
    glm.full <- glm(g~x-1, family=binomial(link='probit'))
    glm.reduced <- glm(g~x.reduced-1, family=binomial(link='probit'))
    beta.hat.full <- coef(glm.full)
    beta.hat.reduced <- coef(glm.reduced)
    auc.hat.full.hat <- auc.hat(x.0%*%beta.hat.full,x.1%*%beta.hat.full)
    auc.hat.reduced.hat <- auc.hat(x.0.reduced%*%beta.hat.reduced,x.1.reduced%*%beta.hat.reduced)
    diff.hat.coef.hat <- auc.hat.full.hat - auc.hat.reduced.hat
    ## auc.full.hat <- auc.probit(mu.full,Sigma.full,beta.hat.full)
    ## auc.reduced.hat <- auc.probit(mu.reduced,Sigma.reduced,beta.hat.reduced)
    ## diff.coef.hat <- auc.full.hat - auc.reduced.hat
    ## auc.hajek.full <- auc.hajek.probit(x.0,x.1,mu.full,Sigma.full,beta.hat.full,terms.only=TRUE,IID=TRUE) #1
    ## auc.hajek.reduced <- auc.hajek.probit(x.0.reduced,x.1.reduced,mu.reduced,Sigma.reduced,beta.hat.reduced,terms.only=TRUE,IID=TRUE) #1
    ## auc.hajek.full <- auc.hajek(x=x.0%*%beta.hat.full,y=x.1%*%beta.hat.full,terms.only=TRUE,IID=TRUE) #3
    ## auc.hajek.reduced <- auc.hajek(x=x.0.reduced%*%beta.hat.reduced,y=x.1.reduced%*%beta.hat.reduced,terms.only=TRUE,IID=TRUE) #3
    ## ## taylor <-  t(auc.probit.deriv(beta.hat.full,mu.full,Sigma.full))%*%infl.probit(x,g,beta.full) - t(auc.probit.deriv(beta.hat.reduced,mu.reduced,Sigma.reduced))%*%infl.probit(x.reduced,g,beta.reduced) #1
    ## taylor <-  t(auc.probit.deriv(beta.hat.full,mu.full,Sigma.full))%*%infl.probit(x,g,beta.hat.full) - t(auc.probit.deriv(beta.hat.reduced,mu.reduced,Sigma.reduced))%*%infl.probit(x.reduced,g,beta.hat.reduced) #2
    ## diff.hajek.coef <- auc.hajek.full - auc.hajek.reduced 
    ## target <- diff.hat.coef.hat - diff.coef
    ## approx <- as.numeric(diff.hajek.coef + taylor) 
    ## target / sd(approx)
    var.hat <- delong.var(x=cbind(x.0%*%beta.hat.full,x.0.reduced%*%beta.hat.reduced),y=cbind(x.1%*%beta.hat.full,x.1.reduced%*%beta.hat.reduced))
    (diff.hat.coef.hat - diff.coef) / sqrt(var.hat)
})
qqnorm(stats)
abline(0,1)



## re-did to get rid of intercept-related code, since the intercept
## doesn't affect the auc diff. can probably delete below.

## 11 LDA
## data

## ## 11a [from 10a] check o(1/sqrt(n)) convergence of hajek projection 
## require(mvtnorm)
## start <- Sys.time()
## source('misc.R')
## set.seed(1)
## p.full <- 6
## p.red <- 3
## ## n <- 1e3
## ns <- round(seq(1e2,1e3,len=2e1))
## ## by.n <- sapply(ns, function(n) {
## ## terms <- replicate(2e1, {
## ## tryCatch({
## ## reduced.len <- 1
## ## beta.full <- runif(p)/4
## ## gamma <- runif(p)/4
## Sigma.full <- matrix(rnorm(p.full^2),nrow=p.full)
## Sigma.full <- Sigma.full%*%t(Sigma.full)
## Sigma.full[1,] <- Sigma.full[,1] <- 0
## Sigma.red <- Sigma.full[1:p.red,1:p.red]
## beta.full <- runif(p.full)
## mu.1.full <- mu.diff.full <- Sigma.full%*%beta.full
## mu.1.red <- mu.diff.red <- mu.diff.full[1:p.red]
## mu.0.full <- 0*mu.diff.full; mu.0.red <- 0*mu.diff.red
## beta.red <- c(0,solve(Sigma.red[-1,-1])%*%mu.diff.red[-1])
## ## params.full <- list(mu.0.full,mu.1.full,Sigma.full,beta.full)
## ## params.red <- list(mu.0.red,mu.1.red,Sigma.red,beta.red)
## by.n <- sapply(ns, function(n) {
##     cat('.')
##     pairs <- replicate(1e2, {
##         x.0 <- rmvnorm(n,mu.0.full,Sigma.full)
##         x.1 <- rmvnorm(n,mu.1.full,Sigma.full)
##         x.0.red <- x.0[,1:p.red]; x.1.red <- x.1[,1:p.red]
##         beta.hat.full <- lda.coefs(x.0,x.1)
##         beta.hat.red <- lda.coefs(x.0.red,x.1.red)
##         auc.hat.full.hat <- auc.hat(x.0%*%beta.hat.full,x.1%*%beta.hat.full)
##         auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
##         diff.hat.coef.hat <- auc.hat.full.hat-auc.hat.red.hat
##         auc.full.hat <- auc.lda(mu.diff.full,Sigma.full,beta.hat.full)
##         auc.red.hat <- auc.lda(mu.diff.red,Sigma.red,beta.hat.red)
##         diff.coef.hat <- auc.full.hat-auc.red.hat
##         auc.hajek.full.hat <- auc.hajek.lda(x.0,x.1,mu.0.full,mu.1.full,beta.hat.full,Sigma.full)
##         auc.hajek.red.hat <- auc.hajek.lda(x.0.red,x.1.red,mu.0.red,mu.1.red,beta.hat.red,Sigma.red)
##         diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
##         diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef.hat
##     })
## })
## Sys.time() - start
## sds <- apply(by.n,2,sd)
## lm0 <- lm(log(sds) ~ log(ns))
## lm0


## ## refactoring auc.lda, auc.hajek
## require(mvtnorm)
## start <- Sys.time()
## source('misc.R')
## set.seed(1)
## p.full <- 6
## p.red <- 3
## ## n <- 1e3
## ns <- round(seq(1e2,1e3,len=2e1))
## Sigma.full <- matrix(rnorm(p.full^2),nrow=p.full)
## Sigma.full <- Sigma.full%*%t(Sigma.full)
## Sigma.full[1,] <- Sigma.full[,1] <- 0
## Sigma.red <- Sigma.full[1:p.red,1:p.red]
## beta.full <- runif(p.full)
## mu.1.full <- mu.diff.full <- Sigma.full%*%beta.full
## mu.1.red <- mu.diff.red <- mu.diff.full[1:p.red]
## mu.0.full <- 0*mu.diff.full; mu.0.red <- 0*mu.diff.red
## beta.red <- c(0,solve(Sigma.red[-1,-1])%*%mu.diff.red[-1])
## params.full <- list(mu.0=mu.0.full,mu.1=mu.1.full,Sigma=Sigma.full,beta=beta.full)
## params.red <- list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red,beta=beta.red)
## by.n <- sapply(ns, function(n) {
##     cat('.')
##     pairs <- replicate(1e2, {
##         x.0 <- rmvnorm(n,mu.0.full,Sigma.full)
##         x.1 <- rmvnorm(n,mu.1.full,Sigma.full)
##         x.0.red <- x.0[,1:p.red]; x.1.red <- x.1[,1:p.red]
##         beta.hat.full <- lda.coefs(x.0,x.1)
##         beta.hat.red <- lda.coefs(x.0.red,x.1.red)
##         auc.hat.full.hat <- auc.hat(x.0%*%beta.hat.full,x.1%*%beta.hat.full)
##         auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
##         diff.hat.coef.hat <- auc.hat.full.hat-auc.hat.red.hat
##         ## auc.full.hat <- auc.lda(mu.diff.full,Sigma.full,beta.hat.full)
##         auc.full.hat <- auc.lda(beta.hat.full,params.full)
##         ## auc.red.hat <- auc.lda(mu.diff.red,Sigma.red,beta.hat.red)
##         auc.red.hat <- auc.lda(beta.hat.red,params.red)
##         diff.coef.hat <- auc.full.hat-auc.red.hat
##         ## auc.hajek.full.hat <- auc.hajek.lda(x.0,x.1,mu.0.full,mu.1.full,beta.hat.full,Sigma.full)
##         auc.hajek.full.hat <- auc.hajek.lda(x.0,x.1,beta.hat.full,params.full,terms.only=TRUE)
##         auc.hajek.red.hat <- auc.hajek.lda(x.0.red,x.1.red,beta.hat.red,params.red,terms.only=TRUE)
##         diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
##         diff.hat.coef.hat-diff.coef.hat - mean(diff.hajek.coef.hat)
##     })
## })
## Sys.time() - start
## sds <- apply(by.n,2,sd)
## lm0 <- lm(log(sds) ~ log(ns))
## lm0


## ## checking var(hajek) is the delong var
## require(mvtnorm)
## source('misc.R')
## start <- Sys.time()
## set.seed(1)
## p.full <- 6
## p.red <- 3
## n <- 5e2
## ns <- round(seq(1e2,1e3,len=2e1))
## Sigma.full <- matrix(rnorm(p.full^2),nrow=p.full)
## Sigma.full <- Sigma.full%*%t(Sigma.full)
## Sigma.full[1,] <- Sigma.full[,1] <- 0
## Sigma.red <- Sigma.full[1:p.red,1:p.red]
## beta.full <- runif(p.full)
## mu.1.full <- mu.diff.full <- Sigma.full%*%beta.full
## mu.1.red <- mu.diff.red <- mu.diff.full[1:p.red]
## mu.0.full <- 0*mu.diff.full; mu.0.red <- 0*mu.diff.red
## beta.red <- c(0,solve(Sigma.red[-1,-1])%*%mu.diff.red[-1])
## params.full <- list(mu.0=mu.0.full,mu.1=mu.1.full,Sigma=Sigma.full,beta=beta.full)
## params.red <- list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red,beta=beta.red)
## ## by.n <- sapply(ns, function(n) {
##     ## cat('.')
##     pairs <- replicate(1e2, {
##         x.0 <- rmvnorm(n,mu.0.full,Sigma.full)
##         x.1 <- rmvnorm(n,mu.1.full,Sigma.full)
##         x.0.red <- x.0[,1:p.red]; x.1.red <- x.1[,1:p.red]
##         beta.hat.full <- lda.coefs(x.0,x.1)
##         beta.hat.red <- lda.coefs(x.0.red,x.1.red)
##         ## auc.hat.full.hat <- auc.hat(x.0%*%beta.hat.full,x.1%*%beta.hat.full)
##         ## auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
##         ## diff.hat.coef.hat <- auc.hat.full.hat-auc.hat.red.hat
##         ## auc.full.hat <- auc.lda(mu.diff.full,Sigma.full,beta.hat.full)
##         ## auc.full.hat <- auc.lda(beta.hat.full,params.full)
##         ## auc.red.hat <- auc.lda(mu.diff.red,Sigma.red,beta.hat.red)
##         ## auc.red.hat <- auc.lda(beta.hat.red,params.red)
##         ## diff.coef.hat <- auc.full.hat-auc.red.hat
##         ## auc.hajek.full.hat <- auc.hajek.lda(x.0,x.1,mu.0.full,mu.1.full,beta.hat.full,Sigma.full)
##         auc.hajek.full.hat <- auc.hajek.lda(x.0,x.1,beta.hat.full,params.full,terms.only=TRUE)
##         auc.hajek.red.hat <- auc.hajek.lda(x.0.red,x.1.red,beta.hat.red,params.red,terms.only=TRUE)
##         diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
##         ## diff.hat.coef.hat-diff.coef.hat - mean(diff.hajek.coef.hat)
##         c(old=delong.var(x=cbind(x.0%*%beta.hat.full,x.0.red%*%beta.hat.red),y=cbind(x.1%*%beta.hat.full,x.1.red%*%beta.hat.red)), try=var(diff.hajek.coef.hat)/length(diff.hajek.coef.hat))
##     })
## plot(pairs['old',],pairs['try',]); abline(0,1)



## ## 11b influence function


## pi.1 <- runif(1,1/4,3/4)
## ns <- round(seq(1e2,1e3,len=30))
## by.n <- sapply(ns, function(n) {
##     replicate(1e2, {
##         n.1 <- rbinom(1,n,pi.1)
##         pi.1.hat <- n.1/n
##         obs <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1))
##         infl <- pi.1.hat/pi.1/(1-pi.1)-1/(1-pi.1)
##         obs-infl
##     })
## })
## mad <- colMeans(abs(by.n))
## plot(ns,mad)
## lm0 <- lm(log(mad)~log(ns))
## coef(lm0)
## curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)


## p <- 5
## mu.0 <- runif(p)
## mu.1 <- runif(p)
## Sigma <- matrix(runif(p^2),p)
## Sigma <- Sigma%*%t(Sigma)
## ns <- round(seq(1e2,1e3,len=30))
## by.n <- sapply(ns, function(n) {
##     replicate(1e2, {
##         x.0 <- rmvnorm(n,mu.0,Sigma)
##         x.1 <- rmvnorm(n,mu.1,Sigma)
##         mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
##         Sigma.hat <- ((nrow(x.0)-1)*cov(x.0) + (nrow(x.1)-1)*cov(x.1))/(nrow(x.0)+nrow(x.1))
##         obs <- t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat - t(mu.0)%*%solve(Sigma)%*%mu.0
##         infl <- 2*t(mu.0)%*%solve(Sigma)%*%(mu.0.hat-mu.0)
##         obs-infl
##     })
## })
## mad <- colMeans(abs(by.n))
## plot(ns,mad)
## lm0 <- lm(log(mad)~log(ns))
## coef(lm0)
## curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)



## p <- 1
## mu <- runif(p)
## Sigma <- matrix(runif(p^2),p)
## Sigma <- Sigma%*%t(Sigma)
## ns <- round(seq(1e2,1e3,len=30))
## a <- matrix(runif(p),ncol=1)
## by.n <- sapply(ns, function(n) {
##     replicate(1e2, {
##         x <- rmvnorm(n,mu,Sigma)
##         ## x.1 <- rmvnorm(n,mu.1,Sigma)
##         ## mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
##         Sigma.hat <- cov(x)
##         ## obs <- t(mu.0.hat)%*%solve(Sigma.hat)%*%mu.0.hat - t(mu.0)%*%solve(Sigma)%*%mu.0
##         ## infl <- 2*t(mu.0)%*%solve(Sigma)%*%(mu.0.hat-mu.0)
##         ## obs-infl
##         obs <- t(a)%*%(solve(Sigma.hat)-solve(Sigma))%*%a
##         try <- -t(a)%*%solve(Sigma)%*%(Sigma.hat - Sigma)%*%solve(Sigma)%*%a
##         obs-try
##     })
## })
## mad <- colMeans(abs(by.n))
## plot(ns,mad)
## lm0 <- lm(log(mad)~log(ns))
## coef(lm0)
## curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)


## require(mvtnorm)
## set.seed(1)
## p <- 5
## mu.0 <- runif(p)
## mu.1 <- runif(p)
## Sigma <- matrix(runif(p^2),p)
## Sigma <- Sigma%*%t(Sigma)
## ns <- round(seq(1e2,1e3,len=30))
## pi.1 <- 2/3
## by.n <- sapply(ns, function(n) {
##     diffs <- replicate(1e2, {
##         n.1 <- rbinom(1,n,pi.1)
##         n.0 <- n-n.1
##         x.0 <- t(rmvnorm(n.0,mu.0,Sigma))
##         x.1 <- t(rmvnorm(n.1,mu.1,Sigma))
##         mu.0.hat <- rowMeans(x.0);  mu.1.hat <- rowMeans(x.1)
##         ## n.0 <- ncol(x.0); n.1 <- ncol(x.1)
##         g <- c(rep(0,n.0),rep(1,n.1))
##         pi.1.hat <- n.1/(n.0+n.1)
##         ## obs <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1))
##         ## infl <- pi.1.hat/pi.1/(1-pi.1)-1/(1-pi.1)
##         ## Sigma.hat <- ((nrow(x.0)-1)*cov(x.0) + (nrow(x.1)-1)*cov(x.1))/(nrow(x.0)+nrow(x.1))
##         ## obs.0 <- t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat - t(mu.0)%*%solve(Sigma)%*%mu.0
##         ## obs.1 <- t(mu.1.hat)%*%solve(Sigma)%*%mu.1.hat - t(mu.1)%*%solve(Sigma)%*%mu.1
##         ## infl.0 <- 2*t(mu.0)%*%solve(Sigma)%*%(mu.0.hat-mu.0)
##         ## infl.1 <- 2*t(mu.1)%*%solve(Sigma)%*%(mu.1.hat-mu.1)
##         ## obs.mu <- t(mu.0)%*%solve(Sigma)%*%mu.0 - t(mu.1)%*%solve(Sigma)%*%mu.1 - (t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat - t(mu.1.hat)%*%solve(Sigma)%*%mu.1.hat)
##         ## ## obs <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1)) + t(mu.0)%*%solve(Sigma)%*%mu.0 - t(mu.1)%*%solve(Sigma)%*%mu.1 - (t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat - t(mu.1.hat)%*%solve(Sigma)%*%mu.1.hat)
##         ## infl.mu <- c(2*t(mu.0)%*%solve(Sigma)%*%(x.0-mu.0) * -(n.0+n.1)/n.0, 2*t(mu.1)%*%solve(Sigma)%*%(x.1-mu.1) * (n.0+n.1)/n.1)
##         ## obs.pi <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1))
##         ## ## infl.pi <- pi.1.hat/pi.1/(1-pi.1)-1/(1-pi.1)
##         ## infl.pi <- g/pi.1/(1-pi.1)-1/(1-pi.1)
##         obs.0 <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1)) + 1/2*(t(mu.1.hat)%*%solve(Sigma)%*%mu.1.hat - t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat) - 1/2*(t(mu.1)%*%solve(Sigma)%*%mu.1 - t(mu.0)%*%solve(Sigma)%*%mu.0)
##         ## infl <- mean(infl.pi+infl.mu)
##         infl.0 <- g/pi.1/(1-pi.1)-1/(1-pi.1) + 1/2*c(2*t(mu.0)%*%solve(Sigma)%*%(x.0-mu.0) * -n/n.0, 2*t(mu.1)%*%solve(Sigma)%*%(x.1-mu.1) * n/n.1)
##         obs.0 - mean(infl.0)
##         ## ## obs.pi - mean(infl.pi)
##         ## ## obs-infl
##         ## obs <- solve(Sigma)%*%(mu.1.hat - mu.0.hat) - solve(Sigma)%*%(mu.1 - mu.0)
##         ## ## infl <- cbind(solve(Sigma)%*%(x.1-mu.1)*n/n.1,-solve(Sigma)%*%(x.0-mu.0)*n/n.0)
##         ## infl <- cbind(-solve(Sigma)%*%(x.0-mu.0)*n/n.0,solve(Sigma)%*%(x.1-mu.1)*n/n.1)
##         ## obs-rowMeans(infl)
##         ## obs <- rbind(obs.0,obs)
##         ## infl <- rbind(infl.0,infl)
##         ## obs - rowMeans(infl)
##     })
##     ## diffs <- diffs[,,]
##     ## colSums(abs(diffs))
## })
## mad <- colMeans(abs(by.n))
## plot(ns,mad)
## lm0 <- lm(log(mad)~log(ns))
## coef(lm0)
## curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)


## ## encapsulate
## ## expects x in model.matrix format
## source('misc.R')
## infl.lda <- function(x,g,params,terms.only=TRUE) {
##     ## browser()
##     x <- t(x)
##     mu.0 <- params$mu.0
##     mu.1 <- params$mu.1
##     Sigma <- params$Sigma
##     n <- length(g)
##     n.1 <- sum(g); n.0 <- n-n.1
##     x.0 <- t(x[,g==0]); x.1 <- t(x[,g==1])
##     g <- sort(g)
##     mu.0.hat <- rowMeans(x.0);  mu.1.hat <- rowMeans(x.1)
##     pi.1.hat <- n.1/(n.0+n.1)
##     infl.0 <-   g/pi.1/(1-pi.1)-1/(1-pi.1) - t(mu.1)%*%solve(Sigma)%*%(x-mu.1) * n/n.1 * g  + t(mu.0)%*%solve(Sigma)%*%(x-mu.0) * n/n.0 * (1-g)
##     infl <- t(g*t(solve(Sigma)%*%(x-mu.1)*n/n.1)  -  (1-g)*t(solve(Sigma)%*%(x-mu.0)*n/n.0))
##     infl <- rbind(infl.0,infl)
##     if(terms.only) return(infl) else return(rowMeans(infl))
## }
## require(mvtnorm)
## set.seed(1)
## p <- 5
## mu.0 <- runif(p)
## mu.1 <- runif(p)
## Sigma <- matrix(runif(p^2),p)
## Sigma <- Sigma%*%t(Sigma)
## beta <- solve(Sigma)%*%(mu.1-mu.0)
## beta.0 <- log(pi.1/(1-pi.1)) - 1/2*(t(mu.1)%*%solve(Sigma)%*%mu.1-t(mu.0)%*%solve(Sigma)%*%mu.0)
## ns <- round(seq(1e2,1e3,len=30))
## pi.1 <- 2/3
## by.n <- sapply(ns, function(n) {
##     diffs <- replicate(1e2, {
##         n.1 <- rbinom(1,n,pi.1)
##         n.0 <- n-n.1
##         x.0 <- cbind(1,rmvnorm(n.0,mu.0,Sigma))
##         x.1 <- cbind(1,rmvnorm(n.1,mu.1,Sigma))
##         x <- t(rbind(x.0,x.1)[,-1])
##         mu.0.hat <- colMeans(x.0)[-1];  mu.1.hat <- colMeans(x.1)[-1]
##         ## n.0 <- ncol(x.0); n.1 <- ncol(x.1)
##         g <- c(rep(0,n.0),rep(1,n.1))
##         ## pi.1.hat <- n.1/(n.0+n.1)
##         ## obs <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1))
##         ## infl <- pi.1.hat/pi.1/(1-pi.1)-1/(1-pi.1)
##         ## Sigma.hat <- ((nrow(x.0)-1)*cov(x.0) + (nrow(x.1)-1)*cov(x.1))/(nrow(x.0)+nrow(x.1))
##         ## obs.0 <- t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat - t(mu.0)%*%solve(Sigma)%*%mu.0
##         ## obs.1 <- t(mu.1.hat)%*%solve(Sigma)%*%mu.1.hat - t(mu.1)%*%solve(Sigma)%*%mu.1
##         ## infl.0 <- 2*t(mu.0)%*%solve(Sigma)%*%(mu.0.hat-mu.0)
##         ## infl.1 <- 2*t(mu.1)%*%solve(Sigma)%*%(mu.1.hat-mu.1)
##         ## obs.mu <- t(mu.0)%*%solve(Sigma)%*%mu.0 - t(mu.1)%*%solve(Sigma)%*%mu.1 - (t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat - t(mu.1.hat)%*%solve(Sigma)%*%mu.1.hat)
##         ## ## obs <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1)) + t(mu.0)%*%solve(Sigma)%*%mu.0 - t(mu.1)%*%solve(Sigma)%*%mu.1 - (t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat - t(mu.1.hat)%*%solve(Sigma)%*%mu.1.hat)
##         ## infl.mu <- c(2*t(mu.0)%*%solve(Sigma)%*%(x.0-mu.0) * -(n.0+n.1)/n.0, 2*t(mu.1)%*%solve(Sigma)%*%(x.1-mu.1) * (n.0+n.1)/n.1)
##         ## obs.pi <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1))
##         ## ## infl.pi <- pi.1.hat/pi.1/(1-pi.1)-1/(1-pi.1)
##         ## infl.pi <- g/pi.1/(1-pi.1)-1/(1-pi.1)
##         beta.hat <- coefs.lda(x.0,x.1,Sigma)
##         ## beta.hat.plim <- c(log(pi.1/(1-pi.1)) + t(mu.0)%*%solve(Sigma)%*%mu.0/2 - t(mu.1)%*%solve(Sigma)%*%mu.1/2, solve(Sigma)%*%(mu.1 - mu.0))
##         ## obs.0 <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1)) -t(mu.1.hat)%*%solve(Sigma)%*%mu.1.hat/2  + t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat/2  + t(mu.1)%*%solve(Sigma)%*%mu.1  / 2 -  t(mu.0)%*%solve(Sigma)%*%mu.0/2 
##         ## obs.0 <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1))
##         ## infl <- mean(infl.pi+infl.mu)
##         ## infl.0 <- g/pi.1/(1-pi.1)-1/(1-pi.1) + c(2*t(mu.0)%*%solve(Sigma)%*%(x.0-mu.0) * -n/n.0, 2*t(mu.1)%*%solve(Sigma)%*%(x.1-mu.1) * n/n.1)
##         ## infl.0 <-  g/pi.1/(1-pi.1)-1/(1-pi.1) - t(mu.1)%*%solve(Sigma)%*%(x-mu.1) * n/n.1 * g  + t(mu.0)%*%solve(Sigma)%*%(x-mu.0) * n/n.0 * (1-g)
##         ## obs.0 - mean(infl.0)
##         ## obs.pi - mean(infl.pi)
##         ## obs-infl
##         ## obs <- solve(Sigma)%*%(mu.1.hat - mu.0.hat) - solve(Sigma)%*%(mu.1 - mu.0)
##         ## infl <- cbind(solve(Sigma)%*%(x.1-mu.1)*n/n.1,-solve(Sigma)%*%(x.0-mu.0)*n/n.0)
##         ## infl <- cbind(-solve(Sigma)%*%(x.0-mu.0)*n/n.0,solve(Sigma)%*%(x.1-mu.1)*n/n.1)
##         ## infl <- t(g*t(solve(Sigma)%*%(x-mu.1)*n/n.1)  -  (1-g)*t(solve(Sigma)%*%(x-mu.0)*n/n.0))
##         ## obs-rowMeans(infl)
##         ## obs <- rbind(obs.0,obs)
##         ## infl <- rbind(infl.0,infl)
##         ## obs - rowMeans(infl)
##         beta.hat - c(beta.0,beta) - rowMeans(infl.lda(t(x),g,list(mu.0=mu.0,mu.1=mu.1,Sigma=Sigma)))
##         ## obs.0-mean(infl.0)
##     })
##     ## diffs <- diffs[,,]
##     ## colSums(abs(diffs))
## })
## mad <- colMeans(abs(by.n))
## plot(ns,mad)
## lm0 <- lm(log(mad)~log(ns))
## coef(lm0)
## curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)



## ## 11c #1 taylor expansion. #2 using influence function. #3 using IID
## ## form of infl. #4 using combined influence function and gradient.

## require(mvtnorm)
## source('misc.R')
## set.seed(1)
## n <- 1.5e2
## p.full <- 6
## p.red <- 3
## Sigma.full <- matrix(runif(p.full^2),nrow=p.full)
## Sigma.full <- Sigma.full%*%t(Sigma.full)
## ## Sigma.full <- rbind(0,cbind(0,Sigma.full))
## Sigma.full[1,] <- Sigma.full[,1] <- 0
## Sigma.red <- Sigma.full[1:p.red,1:p.red]
## mu.0.full <- rep(0,p.full)
## mu.1.full <- runif(p.full)
## mu.0.full[1] <- mu.1.full[1] <- 1
## mu.diff.full <- mu.1.full-mu.0.full
## mu.0.red <- mu.0.full[1:p.red]; mu.1.red <- mu.1.full[1:p.red]
## mu.diff.red <- mu.1.red - mu.0.red
## pi.0 <- 1/3; pi.1 <- 1-pi.0
## beta.full <- c(log(pi.1/pi.0) - 1/2*(t(mu.1.full[-1])%*%solve(Sigma.full[-1,-1])%*%mu.1.full[-1]-t(mu.0.full[-1])%*%solve(Sigma.full[-1,-1])%*%mu.0.full[-1] ) , mu.diff.full[-1]%*%solve(Sigma.full[-1,-1])) # could parametrize lda model by beta rather than mu
## beta.red <- c(log(pi.1/pi.0) - 1/2*(t(mu.1.red[-1])%*%solve(Sigma.red[-1,-1])%*%mu.1.red[-1]-t(mu.0.red[-1])%*%solve(Sigma.red[-1,-1])%*%mu.0.red[-1] ) , mu.diff.red[-1]%*%solve(Sigma.red[-1,-1]))
## beta.combined <- c(beta.full,beta.red)
## ## params.full <- list(mu.0=mu.0.full[,mu.1=mu.1.full,Sigma=Sigma.full,beta=beta.full)
## ## params.red <- list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red,beta=beta.red)
## params.full <- list(mu.0=mu.0.full[-1],mu.1=mu.1.full[-1],Sigma=Sigma.full[-1,-1],beta=beta.full)
## params.red <- list(mu.0=mu.0.red[-1],mu.1=mu.1.red[-1],Sigma=Sigma.red[-1,-1],beta=beta.red)
## auc.full <- auc.lda(beta.full,params.full)
## auc.red <- auc.lda(beta.red,params.red)
## diff.coef <- auc.full-auc.red
## ns <- round(seq(1.5e2,1e3,len=20))
## by.n <- sapply(ns, function(n) {
##     diffs <- replicate(3e2, {
##         n.0 <- rbinom(1,n,pi.0)
##         n.1 <- n-n.0    
##         x.0 <- rmvnorm(n.0,mu.0.full,Sigma.full)
##         x.1 <- rmvnorm(n.1,mu.1.full,Sigma.full)
##         x.full <- rbind(x.0,x.1)
##         g <- c(rep(0,n.0),rep(1,n.1))
##         x.0.red <- x.0[,1:p.red]
##         x.1.red <- x.1[,1:p.red]
##         x.red <- rbind(x.0.red,x.1.red)
##         beta.hat.full <- coefs.lda(x.0,x.1,Sigma.full[-1,-1])
##         beta.hat.red <- coefs.lda(x.0.red,x.1.red,Sigma.red[-1,-1])
##         auc.full.hat <- auc.lda(beta.hat.full,params.full)
##         auc.red.hat <- auc.lda(beta.hat.red,params.red)
##         diff.coef.hat <- auc.full.hat-auc.red.hat
##         obs <- diff.coef.hat - diff.coef
##         ## try <- (beta.hat.full - beta.full)%*%auc.lda.deriv(beta.full,params.full) - (beta.hat.red - beta.red)%*%auc.lda.deriv(beta.hat.red,params.red) #1
##         ## infl.full <- infl.lda(x.full,g,params.full,terms.only=FALSE) #2
##         ## infl.red <- infl.lda(x.red,g,params.red,terms.only=FALSE) #2
##         ## try <- infl.full%*%auc.lda.deriv(beta.full,params.full) - infl.red%*%auc.lda.deriv(beta.hat.red,params.red)
##         ## infl.full <- infl.lda(x.full,g,params.full,terms.only=TRUE) #3
##         ## infl.red <- infl.lda(x.red,g,params.red,terms.only=TRUE) #3
##         ## try <- mean(t(auc.lda.deriv(beta.full,params.full))%*%infl.full - t(auc.lda.deriv(beta.hat.red,params.red))%*%infl.red)
##         infl.full <- infl.lda(x.full,g,params.full,terms.only=TRUE) #4
##         infl.red <- infl.lda(x.red,g,params.red,terms.only=TRUE) #4
##         infl <- rbind(infl.full,infl.red) #4
##         deriv <- t(c(auc.lda.deriv(beta.full,params.full),-auc.lda.deriv(beta.red,params.red))) #4
##         try <- mean(deriv%*%infl) #4
##         try - obs
##     })
## })
## mad <- colMeans(abs(by.n))
## plot(ns,mad)
## lm0 <- lm(log(mad)~log(ns))
## coef(lm0)
## curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)





## 11 LDA data. redo-ing ignoring intercept term, which doesnt affect the index auc.

## ## 11a [from 10a] check o(1/sqrt(n)) convergence of hajek projection 
## require(mvtnorm)
## start <- Sys.time()
## source('misc.R')
## set.seed(1)
## p.full <- 6
## p.red <- 3
## ## n <- 1e3
## ns <- round(seq(1e2,1e3,len=2e1))
## ## by.n <- sapply(ns, function(n) {
## ## terms <- replicate(2e1, {
## ## tryCatch({
## ## reduced.len <- 1
## ## beta.full <- runif(p)/4
## ## gamma <- runif(p)/4
## Sigma.full <- matrix(rnorm(p.full^2),nrow=p.full)
## Sigma.full <- Sigma.full%*%t(Sigma.full)
## Sigma.full[1,] <- Sigma.full[,1] <- 0
## Sigma.red <- Sigma.full[1:p.red,1:p.red]
## beta.full <- runif(p.full)
## mu.1.full <- mu.diff.full <- Sigma.full%*%beta.full
## mu.1.red <- mu.diff.red <- mu.diff.full[1:p.red]
## mu.0.full <- 0*mu.diff.full; mu.0.red <- 0*mu.diff.red
## beta.red <- c(0,solve(Sigma.red[-1,-1])%*%mu.diff.red[-1])
## ## params.full <- list(mu.0.full,mu.1.full,Sigma.full,beta.full)
## ## params.red <- list(mu.0.red,mu.1.red,Sigma.red,beta.red)
## by.n <- sapply(ns, function(n) {
##     cat('.')
##     pairs <- replicate(1e2, {
##         x.0 <- rmvnorm(n,mu.0.full,Sigma.full)
##         x.1 <- rmvnorm(n,mu.1.full,Sigma.full)
##         x.0.red <- x.0[,1:p.red]; x.1.red <- x.1[,1:p.red]
##         beta.hat.full <- lda.coefs(x.0,x.1)
##         beta.hat.red <- lda.coefs(x.0.red,x.1.red)
##         auc.hat.full.hat <- auc.hat(x.0%*%beta.hat.full,x.1%*%beta.hat.full)
##         auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
##         diff.hat.coef.hat <- auc.hat.full.hat-auc.hat.red.hat
##         auc.full.hat <- auc.lda(mu.diff.full,Sigma.full,beta.hat.full)
##         auc.red.hat <- auc.lda(mu.diff.red,Sigma.red,beta.hat.red)
##         diff.coef.hat <- auc.full.hat-auc.red.hat
##         auc.hajek.full.hat <- auc.hajek.lda(x.0,x.1,mu.0.full,mu.1.full,beta.hat.full,Sigma.full)
##         auc.hajek.red.hat <- auc.hajek.lda(x.0.red,x.1.red,mu.0.red,mu.1.red,beta.hat.red,Sigma.red)
##         diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
##         diff.hat.coef.hat-diff.coef.hat - diff.hajek.coef.hat
##     })
## })
## Sys.time() - start
## sds <- apply(by.n,2,sd)
## lm0 <- lm(log(sds) ~ log(ns))
## lm0


## refactoring auc.lda, auc.hajek
require(mvtnorm)
## start <- Sys.time()
source('misc.R')
set.seed(1)
p.full <- 6
p.red <- 3
## n <- 1e3
ns <- round(seq(1e2,1e3,len=2e1))
Sigma.full <- matrix(rnorm(p.full^2),nrow=p.full)
Sigma.full <- Sigma.full%*%t(Sigma.full)
## Sigma.full[1,] <- Sigma.full[,1] <- 0
Sigma.red <- Sigma.full[1:p.red,1:p.red]
beta.full <- runif(p.full)
mu.1.full <- mu.diff.full <- Sigma.full%*%beta.full
mu.1.red <- mu.diff.red <- mu.diff.full[1:p.red]
mu.0.full <- 0*mu.diff.full; mu.0.red <- 0*mu.diff.red
beta.red <- solve(Sigma.red)%*%mu.diff.red
params.full <- list(mu.0=mu.0.full,mu.1=mu.1.full,Sigma=Sigma.full)#,beta=beta.full)
params.red <- list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red)#,beta=beta.red)
by.n <- sapply(ns, function(n) {
    cat('.')
    pairs <- replicate(1e2, {
        x.0 <- rmvnorm(n,mu.0.full,Sigma.full)
        x.1 <- rmvnorm(n,mu.1.full,Sigma.full)
        x.0.red <- x.0[,1:p.red]; x.1.red <- x.1[,1:p.red]
        beta.hat.full <- coefs.lda(x.0,x.1)
        beta.hat.red <- coefs.lda(x.0.red,x.1.red)
        ## beta.hat.full <- beta.full; beta.hat.red <- beta.red ##!!!!
        auc.hat.full.hat <- auc.hat(x.0%*%beta.hat.full,x.1%*%beta.hat.full)
        auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
        diff.hat.coef.hat <- auc.hat.full.hat-auc.hat.red.hat
        ## auc.full.hat <- auc.lda(mu.diff.full,Sigma.full,beta.hat.full)
        auc.full.hat <- auc.lda(beta.hat.full,params.full)
        ## auc.red.hat <- auc.lda(mu.diff.red,Sigma.red,beta.hat.red)
        auc.red.hat <- auc.lda(beta.hat.red,params.red)
        diff.coef.hat <- auc.full.hat-auc.red.hat
        ## auc.hajek.full.hat <- auc.hajek.lda(x.0,x.1,mu.0.full,mu.1.full,beta.hat.full,Sigma.full)
        auc.hajek.full.hat <- auc.hajek.lda(x.0,x.1,beta.hat.full,params.full,terms.only=TRUE)
        auc.hajek.red.hat <- auc.hajek.lda(x.0.red,x.1.red,beta.hat.red,params.red,terms.only=TRUE)
        diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
        diff.hat.coef.hat-diff.coef.hat - mean(diff.hajek.coef.hat)
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)

## checking auc.hajek

## refactoring auc.lda, auc.hajek
require(mvtnorm)
## start <- Sys.time()
source('misc.R')
set.seed(1)
p <- 6
## n <- 1e3
beta <- runif(p)
Sigma <- matrix(rnorm(p^2),nrow=p)
Sigma <- Sigma%*%t(Sigma)
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
## Sigma[1,] <- Sigma[,1] <- 0
mu.1 <- mu.diff <- Sigma%*%beta
mu.0 <- 0*mu.diff#; mu.0.red <- 0*mu.diff.red
ns <- round(seq(1e2,1e3,len=2e1))
by.n <- sapply(ns, function(n) {
    cat('.')
    pairs <- replicate(1e2, {
        x.0 <- rmvnorm(n,mu.0,Sigma)
        x.1 <- rmvnorm(n,mu.1,Sigma)
        x <- as.numeric(x.0%*%beta)
        y <- as.numeric(x.1%*%beta)
        ## F <- function(u)pnorm((u-t(beta)%*%mu.0)/sqrt(quad))
        ## G <- function(u)pnorm((u-t(beta)%*%mu.1)/sqrt(quad))
        ## theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
        ## hajek <- auc.hajek(F,G,x,y,theta,terms.only=FALSE)
        hajek <- auc.hajek.lda(x.0,x.1,beta,list(mu.0=mu.0,mu.1=mu.1,Sigma=Sigma,beta=beta),terms.only=FALSE)
        obs <- auc.hat(x,y) - theta
        obs - hajek
    })
})
## Sys.time() - start
## sds <- apply(by.n,2,sd)
## lm0 <- lm(log(sds) ~ log(ns))
## lm0
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)

dd

## checking var(hajek) is the delong var
require(mvtnorm)
source('misc.R')
n <- 5e2
set.seed(1)
p.full <- 6
p.red <- 3
## n <- 1e3
ns <- round(seq(1e2,1e3,len=2e1))
Sigma.full <- matrix(rnorm(p.full^2),nrow=p.full)
Sigma.full <- Sigma.full%*%t(Sigma.full)
## Sigma.full[1,] <- Sigma.full[,1] <- 0
Sigma.red <- Sigma.full[1:p.red,1:p.red]
beta.full <- runif(p.full)
mu.1.full <- mu.diff.full <- Sigma.full%*%beta.full
mu.1.red <- mu.diff.red <- mu.diff.full[1:p.red]
mu.0.full <- 0*mu.diff.full; mu.0.red <- 0*mu.diff.red
beta.red <- solve(Sigma.red)%*%mu.diff.red
params.full <- list(mu.0=mu.0.full,mu.1=mu.1.full,Sigma=Sigma.full,beta=beta.full)
params.red <- list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red,beta=beta.red)
## by.n <- sapply(ns, function(n) {
    ## cat('.')
    pairs <- replicate(1e2, {
        x.0 <- rmvnorm(n,mu.0.full,Sigma.full)
        x.1 <- rmvnorm(n,mu.1.full,Sigma.full)
        x.0.red <- x.0[,1:p.red]; x.1.red <- x.1[,1:p.red]
        beta.hat.full <- coefs.lda(x.0,x.1)
        beta.hat.red <- coefs.lda(x.0.red,x.1.red)
        auc.hajek.full.hat <- auc.hajek.lda(x.0,x.1,beta.hat.full,params.full,terms.only=TRUE)
        auc.hajek.red.hat <- auc.hajek.lda(x.0.red,x.1.red,beta.hat.red,params.red,terms.only=TRUE)
        diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
        ## diff.hat.coef.hat-diff.coef.hat - mean(diff.hajek.coef.hat)
        c(old=delong.var(x=cbind(x.0%*%beta.hat.full,x.0.red%*%beta.hat.red),y=cbind(x.1%*%beta.hat.full,x.1.red%*%beta.hat.red)), try=var(diff.hajek.coef.hat)/length(diff.hajek.coef.hat))
    })
plot(pairs['old',],pairs['try',]); abline(0,1)






## 11b influence function


pi.1 <- runif(1,1/4,3/4)
ns <- round(seq(1e2,1e3,len=30))
by.n <- sapply(ns, function(n) {
    replicate(1e2, {
        n.1 <- rbinom(1,n,pi.1)
        pi.1.hat <- n.1/n
        obs <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1))
        infl <- pi.1.hat/pi.1/(1-pi.1)-1/(1-pi.1)
        obs-infl
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)


p <- 5
mu.0 <- runif(p)
mu.1 <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
ns <- round(seq(1e2,1e3,len=30))
by.n <- sapply(ns, function(n) {
    replicate(1e2, {
        x.0 <- rmvnorm(n,mu.0,Sigma)
        x.1 <- rmvnorm(n,mu.1,Sigma)
        mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
        Sigma.hat <- ((nrow(x.0)-1)*cov(x.0) + (nrow(x.1)-1)*cov(x.1))/(nrow(x.0)+nrow(x.1))
        obs <- t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat - t(mu.0)%*%solve(Sigma)%*%mu.0
        infl <- 2*t(mu.0)%*%solve(Sigma)%*%(mu.0.hat-mu.0)
        obs-infl
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)



p <- 1
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
ns <- round(seq(1e2,1e3,len=30))
a <- matrix(runif(p),ncol=1)
by.n <- sapply(ns, function(n) {
    replicate(1e2, {
        x <- rmvnorm(n,mu,Sigma)
        ## x.1 <- rmvnorm(n,mu.1,Sigma)
        ## mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
        Sigma.hat <- cov(x)
        ## obs <- t(mu.0.hat)%*%solve(Sigma.hat)%*%mu.0.hat - t(mu.0)%*%solve(Sigma)%*%mu.0
        ## infl <- 2*t(mu.0)%*%solve(Sigma)%*%(mu.0.hat-mu.0)
        ## obs-infl
        obs <- t(a)%*%(solve(Sigma.hat)-solve(Sigma))%*%a
        try <- -t(a)%*%solve(Sigma)%*%(Sigma.hat - Sigma)%*%solve(Sigma)%*%a
        obs-try
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)


require(mvtnorm)
set.seed(1)
p <- 5
mu.0 <- runif(p)
mu.1 <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
ns <- round(seq(1e2,1e3,len=30))
pi.1 <- 2/3
by.n <- sapply(ns, function(n) {
    diffs <- replicate(1e2, {
        n.1 <- rbinom(1,n,pi.1)
        n.0 <- n-n.1
        x.0 <- t(rmvnorm(n.0,mu.0,Sigma))
        x.1 <- t(rmvnorm(n.1,mu.1,Sigma))
        mu.0.hat <- rowMeans(x.0);  mu.1.hat <- rowMeans(x.1)
        ## n.0 <- ncol(x.0); n.1 <- ncol(x.1)
        g <- c(rep(0,n.0),rep(1,n.1))
        pi.1.hat <- n.1/(n.0+n.1)
        ## obs <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1))
        ## infl <- pi.1.hat/pi.1/(1-pi.1)-1/(1-pi.1)
        ## Sigma.hat <- ((nrow(x.0)-1)*cov(x.0) + (nrow(x.1)-1)*cov(x.1))/(nrow(x.0)+nrow(x.1))
        ## obs.0 <- t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat - t(mu.0)%*%solve(Sigma)%*%mu.0
        ## obs.1 <- t(mu.1.hat)%*%solve(Sigma)%*%mu.1.hat - t(mu.1)%*%solve(Sigma)%*%mu.1
        ## infl.0 <- 2*t(mu.0)%*%solve(Sigma)%*%(mu.0.hat-mu.0)
        ## infl.1 <- 2*t(mu.1)%*%solve(Sigma)%*%(mu.1.hat-mu.1)
        ## obs.mu <- t(mu.0)%*%solve(Sigma)%*%mu.0 - t(mu.1)%*%solve(Sigma)%*%mu.1 - (t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat - t(mu.1.hat)%*%solve(Sigma)%*%mu.1.hat)
        ## ## obs <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1)) + t(mu.0)%*%solve(Sigma)%*%mu.0 - t(mu.1)%*%solve(Sigma)%*%mu.1 - (t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat - t(mu.1.hat)%*%solve(Sigma)%*%mu.1.hat)
        ## infl.mu <- c(2*t(mu.0)%*%solve(Sigma)%*%(x.0-mu.0) * -(n.0+n.1)/n.0, 2*t(mu.1)%*%solve(Sigma)%*%(x.1-mu.1) * (n.0+n.1)/n.1)
        ## obs.pi <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1))
        ## ## infl.pi <- pi.1.hat/pi.1/(1-pi.1)-1/(1-pi.1)
        ## infl.pi <- g/pi.1/(1-pi.1)-1/(1-pi.1)
        obs.0 <- log(pi.1.hat/(1-pi.1.hat)) - log(pi.1/(1-pi.1)) + 1/2*(t(mu.1.hat)%*%solve(Sigma)%*%mu.1.hat - t(mu.0.hat)%*%solve(Sigma)%*%mu.0.hat) - 1/2*(t(mu.1)%*%solve(Sigma)%*%mu.1 - t(mu.0)%*%solve(Sigma)%*%mu.0)
        ## infl <- mean(infl.pi+infl.mu)
        infl.0 <- g/pi.1/(1-pi.1)-1/(1-pi.1) + 1/2*c(2*t(mu.0)%*%solve(Sigma)%*%(x.0-mu.0) * -n/n.0, 2*t(mu.1)%*%solve(Sigma)%*%(x.1-mu.1) * n/n.1)
        obs.0 - mean(infl.0)
        ## ## obs.pi - mean(infl.pi)
        ## ## obs-infl
        ## obs <- solve(Sigma)%*%(mu.1.hat - mu.0.hat) - solve(Sigma)%*%(mu.1 - mu.0)
        ## ## infl <- cbind(solve(Sigma)%*%(x.1-mu.1)*n/n.1,-solve(Sigma)%*%(x.0-mu.0)*n/n.0)
        ## infl <- cbind(-solve(Sigma)%*%(x.0-mu.0)*n/n.0,solve(Sigma)%*%(x.1-mu.1)*n/n.1)
        ## obs-rowMeans(infl)
        ## obs <- rbind(obs.0,obs)
        ## infl <- rbind(infl.0,infl)
        ## obs - rowMeans(infl)
    })
    ## diffs <- diffs[,,]
    ## colSums(abs(diffs))
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)


## encapsulate
## expects x in model.matrix format
source('misc.R')
infl.lda <- function(x,g,params,terms.only=TRUE) {
    ## browser()
    x <- t(x)
    mu.0 <- params$mu.0
    mu.1 <- params$mu.1
    Sigma <- params$Sigma
    n <- length(g)
    n.1 <- sum(g); n.0 <- n-n.1
    x.0 <- t(x[,g==0]); x.1 <- t(x[,g==1])
    g <- sort(g)
    mu.0.hat <- rowMeans(x.0);  mu.1.hat <- rowMeans(x.1)
    infl <- t(g*t(solve(Sigma)%*%(x-mu.1)*n/n.1)  -  (1-g)*t(solve(Sigma)%*%(x-mu.0)*n/n.0))
    ## infl <- rbind(infl.0,infl)
    if(terms.only) return(infl) else return(rowMeans(infl))
}
require(mvtnorm)
set.seed(1)
p <- 5
mu.0 <- runif(p)
mu.1 <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
beta <- solve(Sigma)%*%(mu.1-mu.0)
ns <- round(seq(1e2,1e3,len=30))
pi.1 <- 2/3
by.n <- sapply(ns, function(n) {
    diffs <- replicate(1e2, {
        n.1 <- rbinom(1,n,pi.1)
        n.0 <- n-n.1
        x.0 <- rmvnorm(n.0,mu.0,Sigma)
        x.1 <- rmvnorm(n.1,mu.1,Sigma)
        x <- rbind(x.0,x.1)
        mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
        ## n.0 <- ncol(x.0); n.1 <- ncol(x.1)
        g <- c(rep(0,n.0),rep(1,n.1))
        beta.hat <- coefs.lda(x.0,x.1,Sigma)
        beta.hat - beta - rowMeans(infl.lda(x,g,list(mu.0=mu.0,mu.1=mu.1,Sigma=Sigma)))
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)



## 11c #1 taylor expansion. #2 using influence function. #3 using IID
## form of infl. #4 using combined influence function and gradient.

require(mvtnorm)
source('misc.R')
set.seed(1)
## n <- 1.5e2
p.full <- 6
p.red <- 3
Sigma.full <- matrix(runif(p.full^2),nrow=p.full)
Sigma.full <- Sigma.full%*%t(Sigma.full)
## Sigma.full <- rbind(0,cbind(0,Sigma.full))
## Sigma.full[1,] <- Sigma.full[,1] <- 0
Sigma.red <- Sigma.full[1:p.red,1:p.red]
mu.0.full <- rep(0,p.full)
mu.1.full <- runif(p.full)
## mu.0.full[1] <- mu.1.full[1] <- 1
mu.diff.full <- mu.1.full-mu.0.full
mu.0.red <- mu.0.full[1:p.red]; mu.1.red <- mu.1.full[1:p.red]
mu.diff.red <- mu.1.red - mu.0.red
pi.0 <- 1/3; pi.1 <- 1-pi.0
## beta.full <- c(log(pi.1/pi.0) - 1/2*(t(mu.1.full[-1])%*%solve(Sigma.full[-1,-1])%*%mu.1.full[-1]-t(mu.0.full[-1])%*%solve(Sigma.full[-1,-1])%*%mu.0.full[-1] ) , mu.diff.full[-1]%*%solve(Sigma.full[-1,-1])) # could parametrize lda model by beta rather than mu
## beta.red <- c(log(pi.1/pi.0) - 1/2*(t(mu.1.red[-1])%*%solve(Sigma.red[-1,-1])%*%mu.1.red[-1]-t(mu.0.red[-1])%*%solve(Sigma.red[-1,-1])%*%mu.0.red[-1] ) , mu.diff.red[-1]%*%solve(Sigma.red[-1,-1]))
beta.full <-  solve(Sigma.full)%*%mu.diff.full
beta.red <- solve(Sigma.red)%*%mu.diff.red
beta.combined <- c(beta.full,beta.red)
## params.full <- list(mu.0=mu.0.full[,mu.1=mu.1.full,Sigma=Sigma.full,beta=beta.full)
## params.red <- list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red,beta=beta.red)
params.full <- list(mu.0=mu.0.full,mu.1=mu.1.full,Sigma=Sigma.full,beta=beta.full)
params.red <- list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red,beta=beta.red)
auc.full <- auc.lda(beta.full,params.full)
auc.red <- auc.lda(beta.red,params.red)
diff.coef <- auc.full-auc.red
ns <- round(seq(1.5e2,1e3,len=20))
by.n <- sapply(ns, function(n) {
    diffs <- replicate(3e2, {
        n.0 <- rbinom(1,n,pi.0)
        n.1 <- n-n.0    
        x.0 <- rmvnorm(n.0,mu.0.full,Sigma.full)
        x.1 <- rmvnorm(n.1,mu.1.full,Sigma.full)
        x.full <- rbind(x.0,x.1)
        g <- c(rep(0,n.0),rep(1,n.1))
        x.0.red <- x.0[,1:p.red]
        x.1.red <- x.1[,1:p.red]
        x.red <- rbind(x.0.red,x.1.red)
        beta.hat.full <- coefs.lda(x.0,x.1,list(Sigma=Sigma.full))
        beta.hat.red <- coefs.lda(x.0.red,x.1.red,list(Sigma=Sigma.red))
        auc.full.hat <- auc.lda(beta.hat.full,params.full)
        auc.red.hat <- auc.lda(beta.hat.red,params.red)
        diff.coef.hat <- auc.full.hat-auc.red.hat
        obs <- diff.coef.hat - diff.coef
        ## try <- t(beta.hat.full - beta.full)%*%auc.lda.deriv(beta.full,params.full) - t(beta.hat.red - beta.red)%*%auc.lda.deriv(beta.hat.red,params.red) #1
        ## infl.full <- infl.lda(x.full,g,params.full,terms.only=FALSE) #2
        ## infl.red <- infl.lda(x.red,g,params.red,terms.only=FALSE) #2
        ## try <- infl.full%*%auc.lda.deriv(beta.full,params.full) - infl.red%*%auc.lda.deriv(beta.hat.red,params.red) #2
        ## infl.full <- infl.lda(x.full,g,params.full,terms.only=TRUE) #3
        ## infl.red <- infl.lda(x.red,g,params.red,terms.only=TRUE) #3
        ## try <- mean(t(auc.lda.deriv(beta.full,params.full))%*%infl.full - t(auc.lda.deriv(beta.hat.red,params.red))%*%infl.red) #3
        infl.full <- infl.lda(x.full,g,params.full,terms.only=TRUE) #4
        infl.red <- infl.lda(x.red,g,params.red,terms.only=TRUE) #4
        infl <- rbind(infl.full,infl.red) #4
        deriv <- t(c(auc.lda.deriv(beta.full,params.full),-auc.lda.deriv(beta.red,params.red))) #4
        try <- mean(deriv%*%infl) #4
        try - obs
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)




 ## 11d combine the two parts
require(mvtnorm)
source('misc.R')
set.seed(1)
## n <- 1.5e2
p.full <- 6
p.red <- 3
Sigma.full <- matrix(runif(p.full^2),nrow=p.full)
Sigma.full <- Sigma.full%*%t(Sigma.full)
## Sigma.full <- rbind(0,cbind(0,Sigma.full))
## Sigma.full[1,] <- Sigma.full[,1] <- 0
Sigma.red <- Sigma.full[1:p.red,1:p.red]
mu.0.full <- rep(0,p.full)
mu.1.full <- runif(p.full)
## mu.0.full[1] <- mu.1.full[1] <- 1
mu.diff.full <- mu.1.full-mu.0.full
mu.0.red <- mu.0.full[1:p.red]; mu.1.red <- mu.1.full[1:p.red]
mu.diff.red <- mu.1.red - mu.0.red
pi.0 <- 1/3; pi.1 <- 1-pi.0
## beta.full <- c(log(pi.1/pi.0) - 1/2*(t(mu.1.full[-1])%*%solve(Sigma.full[-1,-1])%*%mu.1.full[-1]-t(mu.0.full[-1])%*%solve(Sigma.full[-1,-1])%*%mu.0.full[-1] ) , mu.diff.full[-1]%*%solve(Sigma.full[-1,-1])) # could parametrize lda model by beta rather than mu
## beta.red <- c(log(pi.1/pi.0) - 1/2*(t(mu.1.red[-1])%*%solve(Sigma.red[-1,-1])%*%mu.1.red[-1]-t(mu.0.red[-1])%*%solve(Sigma.red[-1,-1])%*%mu.0.red[-1] ) , mu.diff.red[-1]%*%solve(Sigma.red[-1,-1]))
beta.full <-  solve(Sigma.full)%*%mu.diff.full
beta.red <- solve(Sigma.red)%*%mu.diff.red
beta.combined <- c(beta.full,beta.red)
## params.full <- list(mu.0=mu.0.full[,mu.1=mu.1.full,Sigma=Sigma.full,beta=beta.full)
## params.red <- list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red,beta=beta.red)
params.full <- list(mu.0=mu.0.full,mu.1=mu.1.full,Sigma=Sigma.full,beta=beta.full)
params.red <- list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red,beta=beta.red)
auc.full <- auc.lda(beta.full,params.full)
auc.red <- auc.lda(beta.red,params.red)
diff.coef <- auc.full-auc.red
ns <- round(seq(1.5e2,5e2,len=20))
by.n <- sapply(ns, function(n) {
    cat('.')
    diffs <- replicate(3e2, {
        n.0 <- rbinom(1,n,pi.0)
        n.1 <- n-n.0    
        x.0 <- rmvnorm(n.0,mu.0.full,Sigma.full)
        x.1 <- rmvnorm(n.1,mu.1.full,Sigma.full)
        x.full <- rbind(x.0,x.1)
        g <- c(rep(0,n.0),rep(1,n.1))
        x.0.red <- x.0[,1:p.red]
        x.1.red <- x.1[,1:p.red]
        x.red <- rbind(x.0.red,x.1.red)
        beta.hat.full <- coefs.lda(x.0,x.1,params.full)#Sigma.full)
        beta.hat.red <- coefs.lda(x.0.red,x.1.red,params.red)#Sigma.red)
        auc.full.hat <- auc.lda(beta.hat.full,params.full)
        auc.red.hat <- auc.lda(beta.hat.red,params.red)
        diff.coef.hat <- auc.full.hat-auc.red.hat
        obs.deriv <- diff.coef.hat - diff.coef
        infl.full <- infl.lda(x.full,g,params.full,terms.only=TRUE) 
        infl.red <- infl.lda(x.red,g,params.red,terms.only=TRUE) 
        infl <- rbind(infl.full,infl.red) 
        deriv <- t(c(auc.lda.deriv(beta.full,params.full),-auc.lda.deriv(beta.red,params.red))) 
        try.deriv <- mean(deriv%*%infl) 
        try.deriv - obs.deriv
        auc.hat.full.hat <- auc.hat(x.0%*%beta.hat.full,x.1%*%beta.hat.full)
        auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
        diff.hat.coef.hat <- auc.hat.full.hat-auc.hat.red.hat
        auc.hajek.full.hat <- auc.hajek.lda(x.0,x.1,beta.hat.full,params.full,terms.only=TRUE)
        auc.hajek.red.hat <- auc.hajek.lda(x.0.red,x.1.red,beta.hat.red,params.red,terms.only=TRUE)
        diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
        obs.hajek <- diff.hat.coef.hat-diff.coef.hat
        try.hajek <- mean(diff.hajek.coef.hat)
        obs.hajek - try.hajek
        obs <- diff.hat.coef.hat - diff.coef
        try <- diff.hajek.coef.hat + deriv%*%infl
        obs - mean(try)
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)







## 11e #1 check z-stats. #2 with estimating F, G in hajek part. #3
## with estimated sigma in coefficient estimtes. #4 all params
## estimated other than derivative.
require(mvtnorm)
source('misc.R')
set.seed(1)
## n <- 1.5e2
p.full <- 6
p.red <- 3
Sigma.full <- matrix(runif(p.full^2),nrow=p.full)
Sigma.full <- Sigma.full%*%t(Sigma.full)
## Sigma.full <- rbind(0,cbind(0,Sigma.full))
## Sigma.full[1,] <- Sigma.full[,1] <- 0
Sigma.red <- Sigma.full[1:p.red,1:p.red]
mu.0.full <- rep(0,p.full)
mu.1.full <- runif(p.full)
## mu.0.full[1] <- mu.1.full[1] <- 1
mu.diff.full <- mu.1.full-mu.0.full
mu.0.red <- mu.0.full[1:p.red]; mu.1.red <- mu.1.full[1:p.red]
mu.diff.red <- mu.1.red - mu.0.red
pi.0 <- 1/3; pi.1 <- 1-pi.0
## beta.full <- c(log(pi.1/pi.0) - 1/2*(t(mu.1.full[-1])%*%solve(Sigma.full[-1,-1])%*%mu.1.full[-1]-t(mu.0.full[-1])%*%solve(Sigma.full[-1,-1])%*%mu.0.full[-1] ) , mu.diff.full[-1]%*%solve(Sigma.full[-1,-1])) # could parametrize lda model by beta rather than mu
## beta.red <- c(log(pi.1/pi.0) - 1/2*(t(mu.1.red[-1])%*%solve(Sigma.red[-1,-1])%*%mu.1.red[-1]-t(mu.0.red[-1])%*%solve(Sigma.red[-1,-1])%*%mu.0.red[-1] ) , mu.diff.red[-1]%*%solve(Sigma.red[-1,-1]))
beta.full <-  solve(Sigma.full)%*%mu.diff.full
beta.red <- solve(Sigma.red)%*%mu.diff.red
beta.combined <- c(beta.full,beta.red)
## params.full <- list(mu.0=mu.0.full[,mu.1=mu.1.full,Sigma=Sigma.full,beta=beta.full)
## params.red <- list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red,beta=beta.red)
## params.full <- list(mu.0=mu.0.full,mu.1=mu.1.full,Sigma=Sigma.full)
## params.red <- list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red)
auc.full <- auc.lda(beta.full,list(mu.0=mu.0.full,mu.1=mu.1.full,Sigma=Sigma.full))
auc.red <- auc.lda(beta.red,list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red))
diff.coef <- auc.full-auc.red
ns <- round(seq(1.5e2,5e2,len=20))
## by.n <- sapply(ns, function(n) {
## cat('.')
n <- 3e2
z.stats <- replicate(3e2, {
    n.0 <- rbinom(1,n,pi.0)
    n.1 <- n-n.0    
    x.0 <- rmvnorm(n.0,mu.0.full,Sigma.full)
    x.1 <- rmvnorm(n.1,mu.1.full,Sigma.full)
    x.full <- rbind(x.0,x.1)
    g <- c(rep(0,n.0),rep(1,n.1))
    x.0.red <- x.0[,1:p.red]
    x.1.red <- x.1[,1:p.red]
    x.red <- rbind(x.0.red,x.1.red)
    Sigma.hat.full <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n #4 
    Sigma.hat.red <- (  with(list(x.scaled=scale(x.0.red,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.red,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n #4
    params.full <- list(mu.0=colMeans(x.0),mu.1=colMeans(x.1),Sigma=Sigma.hat.full) #4
    params.red <- list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.hat.red) #4
    ## beta.hat.full <- coefs.lda(x.0,x.1,Sigma.full) #1,#2
    ## beta.hat.red <- coefs.lda(x.0.red,x.1.red,Sigma.red) #1,#2
    beta.hat.full <- coefs.lda(x.0,x.1,params.full)
    beta.hat.red <- coefs.lda(x.0.red,x.1.red,params.red)
    ## derivative part
    ## auc.full.hat <- auc.lda(beta.hat.full,params.full)
    ## auc.red.hat <- auc.lda(beta.hat.red,params.red)
    ## diff.coef.hat <- auc.full.hat-auc.red.hat
    ## obs.deriv <- diff.coef.hat - diff.coef
    infl.full <- infl.lda(x.full,g,params.full,terms.only=TRUE) 
    infl.red <- infl.lda(x.red,g,params.red,terms.only=TRUE) 
    infl <- rbind(infl.full,infl.red) 
    deriv <- t(c(auc.lda.deriv(beta.full,params.full),-auc.lda.deriv(beta.red,params.red))) 
    auc.hat.full.hat <- auc.hat(x.0%*%beta.hat.full,x.1%*%beta.hat.full)
    auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
    diff.hat.coef.hat <- auc.hat.full.hat-auc.hat.red.hat
    auc.hajek.full.hat <- auc.hajek(x=x.0%*%beta.hat.full,y=x.1%*%beta.hat.full,terms.only=TRUE,IID=TRUE) #2
    auc.hajek.red.hat <- auc.hajek(x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red,terms.only=TRUE,IID=TRUE) #2
    diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
    obs <- diff.hat.coef.hat - diff.coef
    try <- as.numeric(diff.hajek.coef.hat + deriv%*%infl)
    obs - mean(try)
    obs / sqrt(var(try) / length(try))
})
qqnorm(z.stats)
abline(0,1)



## 13 trying lda with non-normal data. first make derivative
## calculation nonprametric.

## viz--derivative estimation
set.seed(1)
source('misc.R')
p <- 4
beta <- runif(p)
n <- 1e3
x <- matrix(runif(n*p),ncol=p)
y <- matrix(runif(n*p),ncol=p)
delta <- rep(0,length(beta)); delta[2] <- 1
ts <- seq(-1,1,len=100)
aucs <- sapply(ts, function(t)auc.hat(x%*%(beta+t*delta),y%*%(beta+t*delta)))
plot(ts,aucs)



## viz--lda data--gradient should be 0
require(mvtnorm)
source('misc.R')
set.seed(1)
n <- 1e3
p <- 4
mu.0 <- rep(0,p)
mu.1 <- runif(p)/10
mu.diff <- mu.1-mu.0
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
beta <- solve(Sigma)%*%mu.diff
## beta.0 <- -1/2*(t(mu.1)%*%solve(Sigma)%*%mu.1 - t(mu.0)%*%solve(Sigma)%*%mu.0)

x.0 <- rmvnorm(n,mu.0,Sigma)
x.1 <- rmvnorm(n,mu.1,Sigma)
delta <- rep(0,length(beta)); delta[2] <- 1
ts <- seq(-1,1,len=50)/1e+1
## aucs <- sapply(ts, function(t)auc.hat(beta.0+as.numeric(x.0%*%(beta+t*delta)),beta.0+as.numeric(x.1%*%(beta+t*delta))))
## aucs <- sapply(ts, function(t)auc.hat(x.0%*%(beta+t*delta),x.1%*%(beta+t*delta)))
beta.hat <- lda.coefs(cbind(1,x.0),cbind(1,x.1))[-1]
aucs <- sapply(ts, function(t)auc.hat(x.0%*%(beta.hat+t*delta),x.1%*%(beta.hat+t*delta)))
plot(ts,aucs)
abline(v=0)
abline(v=c(-1,1)*1/sqrt(n),lty=2)


## estimating derivative--lda data
require(mvtnorm)
source('misc.R')
set.seed(1)
n <- 1e3
p <- 4
mu.0 <- rep(0,p)
mu.1 <- runif(p)/10
mu.diff <- mu.1-mu.0
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
x.0 <- rmvnorm(n,mu.0,Sigma)
x.1 <- rmvnorm(n,mu.1,Sigma)
beta.hat <- lda.coefs(cbind(1,x.0),cbind(1,x.1))[-1]
bw <- 1/n^(1/4)
delta <- rep(0,length(beta)); delta[2] <- 1
(  auc.hat(x.0%*%(beta.hat+bw*delta),x.1%*%(beta.hat+bw*delta)) - auc.hat(x.0%*%(beta.hat-bw*delta),x.1%*%(beta.hat-bw*delta))  ) / (2*bw)

## encapsulate
est.gradient <- function(f,u0,bw) {
    p <- length(u0)
    sapply(1:p, function(i) {
        delta <- rep(0,p); delta[i] <- 1
        (f(u0+bw*delta) - f(u0-bw*delta)) / (2*bw)
    })
}
bw <- 1/n^(1/4)
u0 <- beta.hat
f <- function(u)auc.hat(x.0%*%u,x.1%*%u)
est.gradient(f,u0,bw)

## auc.hat looks smooth as a function of beta for n~100, trying numDeriv
require(mvtnorm)
source('misc.R')
set.seed(1)
n <- 5e2
p <- 4
mu.0 <- rep(0,p)
mu.1 <- runif(p)/10
mu.diff <- mu.1-mu.0
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
beta <- solve(Sigma)%*%mu.diff
## beta.0 <- -1/2*(t(mu.1)%*%solve(Sigma)%*%mu.1 - t(mu.0)%*%solve(Sigma)%*%mu.0)
x.0 <- rmvnorm(n,mu.0,Sigma)
x.1 <- rmvnorm(n,mu.1,Sigma)
delta <- rep(0,length(beta)); delta[2] <- 1
f <- function(u)auc.hat(x.0%*%(beta.hat+u*delta),x.1%*%(beta.hat+u*delta))
f <- Vectorize(f)
curve(f)

require(numDeriv)
f <- function(u)auc.hat(x.0%*%u,x.1%*%u)
## f <- Vectorize(f)
grad(f,beta.hat)


## 13b redo 11e with numeric derivative. looked a little skew at first
## but better and better with larger n.

## 11e #1 check z-stats. #2 with estimating F, G in hajek part. #3
## with estimated sigma in coefficient estimtes. #4 all params
## estimated other than derivative.
require(mvtnorm)
require(numDeriv)
source('misc.R')
set.seed(1)
n <- 3e3
p.full <- 6
p.red <- 3
Sigma.full <- matrix(runif(p.full^2),nrow=p.full)
Sigma.full <- Sigma.full%*%t(Sigma.full)
Sigma.red <- Sigma.full[1:p.red,1:p.red]
mu.0.full <- rep(0,p.full)
mu.1.full <- runif(p.full)
mu.diff.full <- mu.1.full-mu.0.full
mu.0.red <- mu.0.full[1:p.red]; mu.1.red <- mu.1.full[1:p.red]
mu.diff.red <- mu.1.red - mu.0.red
pi.0 <- 1/3; pi.1 <- 1-pi.0
beta.full <-  solve(Sigma.full)%*%mu.diff.full
beta.red <- solve(Sigma.red)%*%mu.diff.red
beta.combined <- c(beta.full,beta.red)
auc.full <- auc.lda(beta.full,list(mu.0=mu.0.full,mu.1=mu.1.full,Sigma=Sigma.full))
auc.red <- auc.lda(beta.red,list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red))
diff.coef <- auc.full-auc.red
## ns <- round(seq(1.5e2,5e2,len=20))
## by.n <- sapply(ns, function(n) {
## cat('.')
z.stats <- replicate(3e2, {
    n.0 <- rbinom(1,n,pi.0)
    n.1 <- n-n.0    
    x.0 <- rmvnorm(n.0,mu.0.full,Sigma.full)
    x.1 <- rmvnorm(n.1,mu.1.full,Sigma.full)
    x.full <- rbind(x.0,x.1)
    g <- c(rep(0,n.0),rep(1,n.1))
    x.0.red <- x.0[,1:p.red]
    x.1.red <- x.1[,1:p.red]
    x.red <- rbind(x.0.red,x.1.red)
    Sigma.hat.full <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n #4 
    Sigma.hat.red <- (  with(list(x.scaled=scale(x.0.red,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.red,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n #4
    params.full <- list(mu.0=colMeans(x.0),mu.1=colMeans(x.1),Sigma=Sigma.hat.full) #4
    params.red <- list(mu.0=colMeans(x.0.red),mu.1=colMeans(x.1.red),Sigma=Sigma.hat.red) #4
    ## params.red <- list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.hat.red) #4
    beta.hat.full <- coefs.lda(x.0,x.1,params.full)
    beta.hat.red <- coefs.lda(x.0.red,x.1.red,params.red)
    infl.full <- infl.lda(x.full,g,params.full,terms.only=TRUE) 
    infl.red <- infl.lda(x.red,g,params.red,terms.only=TRUE) 
    infl <- rbind(infl.full,infl.red)
    ## deriv <- t(c(auc.lda.deriv(beta.full,params.full),-auc.lda.deriv(beta.red,params.red))) 
    ## deriv.full <- auc.lda.deriv(beta.full,params.full)
    ## deriv.red <- auc.lda.deriv(beta.red,params.red)
    deriv.full <- grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat.full)
    deriv.red <- grad(function(u)auc.hat(x.0.red%*%u,x.1.red%*%u),beta.hat.red)    
    deriv <- t(c(deriv.full,-deriv.red))
    auc.hat.full.hat <- auc.hat(x.0%*%beta.hat.full,x.1%*%beta.hat.full)
    auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
    diff.hat.coef.hat <- auc.hat.full.hat-auc.hat.red.hat
    auc.hajek.full.hat <- auc.hajek(x=x.0%*%beta.hat.full,y=x.1%*%beta.hat.full,terms.only=TRUE,IID=TRUE) #2
    auc.hajek.red.hat <- auc.hajek(x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red,terms.only=TRUE,IID=TRUE) #2
    diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
    obs <- diff.hat.coef.hat - diff.coef
    try <- as.numeric(diff.hajek.coef.hat + deriv%*%infl)
    obs - mean(try)
    obs / sqrt(var(try) / length(try))
})
## save.image('sessions/13b.RData')
qqnorm(z.stats)
abline(0,1)

hist(z.stats,prob=TRUE)
curve(dnorm,add=TRUE)



## 13c trying sum of exponentials

## 13c-1 formula for auc of gammas with common rate
plot(ecdf(replicate(1e3,sum(rexp(5,2)))))
curve(pgamma(x,5,2),add=TRUE,col=2)

source('misc.R')
p <- 4
auc.gamma <- function(rate.0,rate.1,shape)integrate(function(x)pgamma(x,shape=p,rate=rate.0)*dgamma(x,shape=p,rate=rate.1),0,p/rate.1+4*sqrt(p)/rate.1)$val
pairs <- replicate(1e2, {
    rate.0 <- runif(1,1,3)
    rate.1 <- runif(1,1,3)
    try <- auc.gamma(rate.0,rate.1,shape)
    x <- replicate(1e2,sum(rexp(p,rate.0)))
    y <- replicate(1e2,sum(rexp(p,rate.1)))
    obs <- auc.hat(x,y)
    c(obs,try)
})
plot(pairs[1,],pairs[2,])
abline(0,1)




df <- read.csv('/mnt/c/Users/haben/Desktop/logs_STATISTC516_145047_FA22_20221209-2313.csv')
df <- subset(df,select=c(Time,Event.name))
df <- subset(df,Event.name=='A submission has been submitted.')
df$Time <- as.POSIXct(df$Time, format='%m/%d/%y, %H:%M')
df <- df[order(df$Time),]
df$cumul <- 1:nrow(df)
plot(cumul~Time,data=df,type='l')

df <- read.csv('/mnt/c/Users/haben/OneDrive - University of Massachusetts/umass/teaching/516/STATISTC516_145047_FA22 Grades-20221219_2210-comma_separated.csv')


online <- as.numeric(with(df,pmin(quiz..1..Real.,quiz..5..Real.)))
inclass <- as.numeric(with(df,pmax(quiz..2..Real.,quiz..3..Real.,quiz..4..Real.)))
df$Last.name[order(online - inclass)]


ratios <- with(df,as.numeric(quiz.subtotal..Real.) / as.numeric(homework.subtotal..Real.))
ratios[43] <- NA # stopped attending
plot(ratios,rep(0,length(ratios)),xlab='(quiz subtotal) / (hw subtotal)',ylab='')
points(ratios[which(ratios>3)],0,col=2)


x <- c(1.8,1.5,2,2.5,1.8,2.5,1.6,1.5)
y <- c(51,51,115,150,126,150,118,106)

## 13c-2. z-stats dont look good. debugging...

require(numDeriv)
source('misc.R')
auc.gamma <- function(rate.0,rate.1,shape)integrate(function(x)pgamma(x,shape,rate=rate.0)*dgamma(x,shape,rate=rate.1),0,shape/rate.1+4*sqrt(shape)/rate.1)$val
set.seed(1)
n <- 1e3
p.full <- 3
p.red <- 2
## Sigma.full <- matrix(runif(p.full^2),nrow=p.full)
## Sigma.full <- Sigma.full%*%t(Sigma.full)
## Sigma.red <- Sigma.full[1:p.red,1:p.red]
## mu.0.full <- rep(0,p.full)
## mu.1.full <- runif(p.full)
## mu.diff.full <- mu.1.full-mu.0.full
## mu.0.red <- mu.0.full[1:p.red]; mu.1.red <- mu.1.full[1:p.red]
## mu.diff.red <- mu.1.red - mu.0.red
## pi.0 <- 1/3; pi.1 <- 1-pi.0
beta.full <-  rep(1,p.full)
beta.red <- rep(1,p.red)
rate.0 <- runif(1,1,3)
rate.1 <- runif(1,1,3)
## beta.combined <- c(beta.full,beta.red)
auc.full <- auc.gamma(rate.0,rate.1,p.full)
auc.red <- auc.gamma(rate.0,rate.1,p.red)
## auc.full <- auc.lda(beta.full,list(mu.0=mu.0.full,mu.1=mu.1.full,Sigma=Sigma.full))
## auc.red <- auc.lda(beta.red,list(mu.0=mu.0.red,mu.1=mu.1.red,Sigma=Sigma.red))
diff.coef <- auc.full-auc.red
## ns <- round(seq(1.5e2,5e2,len=20))
## by.n <- sapply(ns, function(n) {
## cat('.')
z.stats <- replicate(3e2, {
    n.0 <- n.1 <- round(n/2)
    x.0 <- matrix(rexp(n.0*p.full,rate.0),ncol=p.full)
    x.1 <- matrix(rexp(n.1*p.full,rate.1),ncol=p.full)
    ## x.1 <- replicate(1e2,sum(rexp(p.full,rate.1)))
    ## x.0 <- rmvnorm(n.0,mu.0.full,Sigma.full)
    ## x.1 <- rmvnorm(n.1,mu.1.full,Sigma.full)
    x.full <- rbind(x.0,x.1)
    g <- c(rep(0,n.0),rep(1,n.1))
    x.0.red <- x.0[,1:p.red]
    x.1.red <- x.1[,1:p.red]
    x.red <- rbind(x.0.red,x.1.red)
    Sigma.hat.full <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n #4 
    Sigma.hat.red <- (  with(list(x.scaled=scale(x.0.red,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.red,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n #4
    params.full <- list(mu.0=colMeans(x.0),mu.1=colMeans(x.1),Sigma=Sigma.hat.full) #4
    params.red <- list(mu.0=colMeans(x.0.red),mu.1=colMeans(x.1.red),Sigma=Sigma.hat.red) #4
    beta.hat.full <- coefs.lda(x.0,x.1,params.full)
    beta.hat.red <- coefs.lda(x.0.red,x.1.red,params.red)
    infl.full <- infl.lda(x.full,g,params.full,terms.only=TRUE) 
    infl.red <- infl.lda(x.red,g,params.red,terms.only=TRUE) 
    infl <- rbind(infl.full,infl.red)
    ## deriv <- t(c(auc.lda.deriv(beta.full,params.full),-auc.lda.deriv(beta.red,params.red))) 
    ## deriv.full <- auc.lda.deriv(beta.full,params.full)
    ## deriv.red <- auc.lda.deriv(beta.red,params.red)
    deriv.full <- grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat.full)
    deriv.red <- grad(function(u)auc.hat(x.0.red%*%u,x.1.red%*%u),beta.hat.red)    
    deriv <- t(c(deriv.full,-deriv.red))
    auc.hat.full.hat <- auc.hat(x.0%*%beta.hat.full,x.1%*%beta.hat.full)
    auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
    diff.hat.coef.hat <- auc.hat.full.hat-auc.hat.red.hat
    auc.hajek.full.hat <- auc.hajek(x=x.0%*%beta.hat.full,y=x.1%*%beta.hat.full,terms.only=TRUE,IID=TRUE) #2
    auc.hajek.red.hat <- auc.hajek(x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red,terms.only=TRUE,IID=TRUE) #2
    diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
    obs <- diff.hat.coef.hat - diff.coef
    try <- as.numeric(diff.hajek.coef.hat + deriv%*%infl)
    obs - mean(try)
    obs / sqrt(var(try) / length(try))
})
## save.image('sessions/13b.RData')
qqnorm(z.stats)
abline(0,1)

hist(z.stats,prob=TRUE)
curve(dnorm,add=TRUE)







require(numDeriv)
source('misc.R')
auc.gamma <- function(rate.0,rate.1,shape)integrate(function(x)pgamma(x,shape,rate=rate.0)*dgamma(x,shape,rate=rate.1),0,shape/rate.1+4*sqrt(shape)/rate.1)$val
set.seed(1)
n <- 1e3
p.full <- 3
p.red <- 2
beta.full <-  rep(1,p.full)
beta.red <- rep(1,p.red)
rate.0 <- runif(1,1,3)
rate.1 <- runif(1,1,3)
auc.full <- auc.gamma(rate.0,rate.1,p.full)
auc.red <- auc.gamma(rate.0,rate.1,p.red)
diff.coef <- auc.full-auc.red
z.stats <- replicate(3e2, {
    n.0 <- n.1 <- round(n/2)
    x.0 <- matrix(rexp(n.0*p.full,rate.0),ncol=p.full)
    x.1 <- matrix(rexp(n.1*p.full,rate.1),ncol=p.full)
    x.full <- rbind(x.0,x.1)
    g <- c(rep(0,n.0),rep(1,n.1))
    x.0.red <- x.0[,1:p.red]
    x.1.red <- x.1[,1:p.red]
    x.red <- rbind(x.0.red,x.1.red)
    Sigma.hat.full <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n #4 
    Sigma.hat.red <- (  with(list(x.scaled=scale(x.0.red,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.red,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n #4
    params.full <- list(mu.0=colMeans(x.0),mu.1=colMeans(x.1),Sigma=Sigma.hat.full) #4
    params.red <- list(mu.0=colMeans(x.0.red),mu.1=colMeans(x.1.red),Sigma=Sigma.hat.red) #4
    beta.hat.full <- coefs.lda(x.0,x.1,params.full)
    ## beta.hat.red <- coefs.lda(x.0.red,x.1.red,params.red)
    infl.full <- infl.lda(x.full,g,params.full,terms.only=TRUE) 
    ## infl.red <- infl.lda(x.red,g,params.red,terms.only=TRUE) 
    ## infl <- rbind(infl.full,infl.red)
    deriv.full <- grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat.full)
    ## deriv.red <- grad(function(u)auc.hat(x.0.red%*%u,x.1.red%*%u),beta.hat.red)    
    ## deriv <- t(c(deriv.full,-deriv.red))
    auc.hat.full.hat <- auc.hat(x.0%*%beta.hat.full,x.1%*%beta.hat.full)
    ## auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
    ## diff.hat.coef.hat <- auc.hat.full.hat-auc.hat.red.hat
    auc.hajek.full.hat <- auc.hajek(x=x.0%*%beta.hat.full,y=x.1%*%beta.hat.full,terms.only=TRUE,IID=TRUE) #2
    ## auc.hajek.red.hat <- auc.hajek(x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red,terms.only=TRUE,IID=TRUE) #2
    ## diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
    ## obs <- diff.hat.coef.hat - diff.coef
    ## try <- as.numeric(diff.hajek.coef.hat + deriv%*%infl)
    obs <- auc.hat.full.hat - auc.full
    try <- as.numeric(auc.hajek.full.hat + deriv.full%*%infl.full)
    obs - mean(try)
    obs / sqrt(var(try) / length(try))
})
## save.image('sessions/13b.RData')
qqnorm(z.stats-mean(z.stats))
abline(0,1)


## 13c-3 checking mu, Sigma, beta.star formulas
require(numDeriv)
source('misc.R')
auc.gamma <- function(rate.0,rate.1,shape)integrate(function(x)pgamma(x,shape,rate=rate.0)*dgamma(x,shape,rate=rate.1),0,shape/rate.1+4*sqrt(shape)/rate.1)$val
set.seed(1)
n <- 1e3
n.0 <- n.1 <- round(n/2)
p.full <- 3
p.red <- 2
beta.full <-  rep(1,p.full)
beta.red <- rep(1,p.red)
rate.0 <- runif(1,1,3)
rate.1 <- runif(1,1,3)
## beta.combined <- c(beta.full,beta.red)
auc.full <- auc.gamma(rate.0,rate.1,p.full)
auc.red <- auc.gamma(rate.0,rate.1,p.red)
diff.coef <- auc.full-auc.red
Sigma.full <- diag(p.full)*(n.0/n*1/rate.0^2 + n.1/n*1/rate.1^2)
Sigma.red <- Sigma.full[1:p.red,1:p.red]
mu.diff.full <- rep(1/rate.1-1/rate.0,p.full)
## mu.diff.red <- rep(1/rate.1-1/rate.0,p.red)
beta.star.full <- solve(Sigma.full)%*%mu.diff.full
beta.star.red <- beta.star.full[1:p.red]
## ns <- round(seq(1.5e2,5e2,len=20))
## by.n <- sapply(ns, function(n) {
## cat('.')
a.full <- runif(p.full); a.red <- runif(p.red)
stats <- replicate(3e2, {
    x.0.full <- matrix(rexp(n.0*p.full,rate.0),ncol=p.full)
    x.1.full <- matrix(rexp(n.1*p.full,rate.1),ncol=p.full)
    x.full <- rbind(x.0.full,x.1.full)
    g <- c(rep(0,n.0),rep(1,n.1))
    x.0.red <- x.0.full[,1:p.red]
    x.1.red <- x.1.full[,1:p.red]
    x.red <- rbind(x.0.red,x.1.red)
    Sigma.hat.full <- (  with(list(x.scaled=scale(x.0.full,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.full,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n #4 
    Sigma.hat.red <- (  with(list(x.scaled=scale(x.0.red,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.red,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n #4
    params.full <- list(mu.0=colMeans(x.0.full),mu.1=colMeans(x.1.full),Sigma=Sigma.hat.full) #4
    params.red <- list(mu.0=colMeans(x.0.red),mu.1=colMeans(x.1.red),Sigma=Sigma.hat.red) #4
    beta.hat.full <- coefs.lda(x.0.full,x.1.full,params.full)
    beta.hat.red <- coefs.lda(x.0.red,x.1.red,params.red)
    c(full=a.full%*%beta.hat.full, red=a.red%*%beta.hat.red)
})
op <- par(mfrow=c(1,2))
hist(stats['full',])
abline(v=a.full%*%beta.star.full,col=2)
hist(stats['red',])
abline(v=a.red%*%beta.star.red,col=2)
par(op)









## 13c-4 check auc.lda.gamma
require(sdprisk)
require(numDeriv)
source('misc.R')
## auc.gamma <- function(rate.0,rate.1,shape)integrate(function(x)pgamma(x,shape,rate=rate.0)*dgamma(x,shape,rate=rate.1),0,shape/rate.1+4*sqrt(shape)/rate.1)$val
set.seed(1)
n <- 1e2
## n.0 <- n.1 <- round(n/2)
p.full <- 3
p.red <- 2
pi.0 <- 1/3
beta.full <-  rep(1,p.full)
beta.red <- rep(1,p.red)
rate.0 <- runif(1,1,3)
rate.1 <- runif(1,1,rate.0)
## beta.combined <- c(beta.full,beta.red)
## auc.full <- auc.gamma(rate.0,rate.1,p.full)
## auc.red <- auc.gamma(rate.0,rate.1,p.red)
diff.coef <- auc.full-auc.red
Sigma.full <- diag(p.full)*(pi.0*1/rate.0^2 + (1-pi.0)*1/rate.1^2)
Sigma.red <- Sigma.full[1:p.red,1:p.red]
mu.diff.full <- rep(1/rate.1-1/rate.0,p.full)
## mu.diff.red <- rep(1/rate.1-1/rate.0,p.red)
beta.star.full <- solve(Sigma.full)%*%mu.diff.full
beta.star.red <- beta.star.full[1:p.red]
F.0.full <- function(q)pgamma(q,shape=p.full,rate=rate.0/unique(beta.star.full))
F.1.full <- function(q)pgamma(q,shape=p.full,rate=rate.1/unique(beta.star.full))
F.0.red <- function(q)pgamma(q,shape=p.red,rate=rate.0/unique(beta.star.red))
F.1.red <- function(q)pgamma(q,shape=p.red,rate=rate.1/unique(beta.star.red))
ns <- round(seq(1.5e2,5e2,len=20))
## P(x<y) where x ~ sum_i beta_i x_i ~ hypoexp(rate.0/beta), y ~ sum_i
## beta_i y_i, x_i ~ Exp(rate.0) ~ hypoexp(rate.1/beta), x_i ~
## Exp(rate.0), y_i ~ Exp(rate.1).
## auc.gamma <- function(beta)integrate(function(x)phypoexp(x,rate=rate.0/beta)*dhypoexp(x,rate=rate.1/beta),0,sum(beta/rate.1) + sqrt(sum((beta / rate.1)^2))*5)$val
## auc.lda.gamma <- function(beta) {
##     if(length(unique(beta))!=1) return( integrate(function(x)phypoexp(x,rate=rate.0/beta)*dhypoexp(x,rate=rate.1/beta),0,sum(beta/rate.1) + sqrt(sum((beta / rate.1)^2))*5)$val )
##     if(length(unique(beta))==1) return( with(list(shape=length(beta)),
##                                              integrate(function(x)pgamma(x,shape,rate=rate.0)*dgamma(x,shape,rate=rate.1),0,shape/rate.1+4*sqrt(shape)/rate.1)$val ) )
##     stop()
## }    
## by.n <- sapply(ns, function(n) {
    ## cat('.')    
beta <- runif(p.full)
n.0 <- round(pi.0*n); n.1 <- n-n.0
stats <- replicate(3e2, {
    x.0.full <- matrix(rexp(n.0*p.full,rate.0),ncol=p.full)
    x.1.full <- matrix(rexp(n.1*p.full,rate.1),ncol=p.full)
    auc.hat(x.0.full%*%beta,x.1.full%*%beta)
    ## plot(ecdf(x.0.full%*%beta))
    ## curve(phypoexp(x,rate=rate.0/beta),add=TRUE,col=2)
    ## hist(x.1.full%*%beta,prob=TRUE)
    ## curve(dhypoexp(x,rate=rate.1/beta),add=TRUE,col=2)
})
hist(stats)
abline(v=auc.lda.gamma(beta,rate.0,rate.1),col=2)
abline(v=mean(stats),col=3)



## need functionality for a sum of pair of gammas, needed to get
## partials. since hypoexponential distr requires all rate params to
## be distinct.
n <- 1e4
shape1 <- runif(1); rate1 <- runif(1)
shape2 <- runif(1); rate2 <- runif(1)
hypergeom1F1 <- function(a,b,z)
    gamma(b)/gamma(a)/gamma(b-a)*integrate(function(u)exp(z*u)*u^(a-1)*(1-u)^(b-a-1),0,1)$val
dgammasum <- Vectorize(function(x,shape1,rate1,shape2,rate2) {
    if(x<=0) {
        return(0)
    } else {
        return( rate1^shape1*rate2^shape2/gamma(shape1+shape2)*exp(-rate1*x)*x^(shape1+shape2-1)*hypergeom1F1(shape2,shape1+shape2,(rate1-rate2)*x) )
    }
}, vectorize.args='x')  
hist(rgamma(n,shape1,rate1)+rgamma(n,shape2,rate2),prob=TRUE)
curve(dgammasum(x,shape1,rate1,shape2,rate2),add=TRUE,col=2)
pgammasum <- Vectorize(function(q,shape1,rate1,shape2,rate2)integrate(function(x)dgammasum(x,shape1,rate1,shape2,rate2),0,q)$val,, vectorize.args='q')
plot(ecdf(rgamma(n,shape1,rate1)+rgamma(n,shape2,rate2)))
## curve(pgammasum(x,shape1,rate1,shape2,rate2),add=TRUE,col=2)
## for some reason curve() not working with pgammasum
x <- seq(0,10,len=1e2)
lines(x,pgammasum(x,shape1,rate1,shape2,rate2),col=2)

## check continuity of auc.lda.gamma
source('misc.R')
require(sdprisk)
set.seed(1)
rate.0 <- runif(1); rate.1 <- runif(1)
ts <- seq(-1,1,len=101)/1
beta <- c(1,2,3)
aucs <- sapply(ts, function(t)auc.lda.gamma(beta+c(0,0,t),rate.0,rate.1))
plot(ts,aucs)
points(0,auc.lda.gamma(beta,rate.0,rate.1),col=2)

source('misc.R')
require(sdprisk)
set.seed(1)
rate.0 <- runif(1); rate.1 <- runif(1)
ts <- seq(-1,1,len=11)/10
beta <- c(1,1,1)
aucs <- sapply(ts, function(t)auc.lda.gamma(beta+c(0,0,t),rate.0,rate.1))
plot(ts,aucs)
## points(0,auc.lda.gamma(beta,rate.0,rate.1),col=2)

n <- 1e3
beta <- c(1,2,2 + ts[1])
auc.hats <- replicate(3e2, {
    x.0 <- matrix(rexp(3*n,rate.0),ncol=3)
    x.1 <- matrix(rexp(3*n,rate.1),ncol=3)
    auc.hat(x.0%*%beta,x.1%*%beta)
})
hist(auc.hats)
abline(v=auc.lda.gamma(beta,rate.0,rate.1),col=2)
abline(v=mean(auc.hats),col=3)

x.0 <- matrix(rexp(3*n,rate.0),ncol=3)
x.1 <- matrix(rexp(3*n,rate.1),ncol=3)
rle0 <- rle(sort(beta))
count <- rle0$len; coef <- rle0$val
xs <- seq(0,80,len=30)
plot(ecdf(x.0%*%beta))
lines(xs, sapply(xs,function(x)pgammasum(x,count[1],rate.0/coef[1],count[2],rate.0/coef[2])),col=2)
hist(x.1%*%beta,prob=TRUE)
curve(dgammasum(x,count[1],rate.1/coef[1],count[2],rate.1/coef[2]),add=TRUE,col=2)
plot(ecdf(x.1%*%beta))
curve(pgammasum(x,count[1],rate.1/coef[1],count[2],rate.1/coef[2]),add=TRUE,col=2)

## check pdf/cdf are continuous
y <- sapply(ts, function(t)dhypoexp(2,c(1,2,2+t)))
plot(ts,y)
points(0,dgammasum(2,1,1,2,2),col=2)
y <- sapply(ts, function(t)phypoexp(2,c(1,2,2+t)))
plot(ts,y)
points(0,pgammasum(2,1,1,2,2),col=2)

dd








## 13c-5 check hajek. fine but seems a little slower than 1/n rate. needed to clip
## beta.hats at 0 otherwise cant compute auc.full.hat, auc.red.hat
## since hypoexp rates>0. verified instead when adding noise to beta.star to
## get beta.hat, convergence is 1/n. also distance from 1/n rate seems
## to go down with sample size as there should be less clipping. n^-.9
## with n>2000.
require(sdprisk)
require(numDeriv)
source('misc.R')
## auc.gamma <- function(rate.0,rate.1,shape)integrate(function(x)pgamma(x,shape,rate=rate.0)*dgamma(x,shape,rate=rate.1),0,shape/rate.1+4*sqrt(shape)/rate.1)$val
set.seed(1)
## n <- 1e2
## n.0 <- n.1 <- round(n/2)
p.full <- 3
p.red <- 2
pi.0 <- 1/3
beta.full <-  rep(1,p.full)
beta.red <- rep(1,p.red)
rate.0 <- runif(1,1,3)
rate.1 <- runif(1,1,rate.0)
## beta.combined <- c(beta.full,beta.red)
diff.coef <- auc.full-auc.red
Sigma.full <- diag(p.full)*(pi.0*1/rate.0^2 + (1-pi.0)*1/rate.1^2)
Sigma.red <- Sigma.full[1:p.red,1:p.red]
mu.diff.full <- rep(1/rate.1-1/rate.0,p.full)
## mu.diff.red <- rep(1/rate.1-1/rate.0,p.red)
beta.star.full <- as.numeric(solve(Sigma.full)%*%mu.diff.full)
beta.star.red <- beta.star.full[1:p.red]
F.0.full <- function(q)pgamma(q,shape=p.full,rate=rate.0/unique(beta.star.full))
F.1.full <- function(q)pgamma(q,shape=p.full,rate=rate.1/unique(beta.star.full))
F.0.red <- function(q)pgamma(q,shape=p.red,rate=rate.0/unique(beta.star.red))
F.1.red <- function(q)pgamma(q,shape=p.red,rate=rate.1/unique(beta.star.red))
ns <- round(seq(1.5e2,1e3,len=20))
## P(x<y) where x ~ sum_i beta_i x_i ~ hypoexp(rate.0/beta), y ~ sum_i
## beta_i y_i, x_i ~ Exp(rate.0) ~ hypoexp(rate.1/beta), x_i ~
## Exp(rate.0), y_i ~ Exp(rate.1).
## auc.lda.gamma <- function(beta) {
##     beta <- as.numeric(beta)
##     if(length(unique(beta))!=1) return( integrate(function(x)phypoexp(x,rate=rate.0/beta)*dhypoexp(x,rate=rate.1/beta),0,sum(beta/rate.1) + sqrt(sum((beta / rate.1)^2))*5)$val )
##     if(length(unique(beta))==1) return( with(list(shape=length(beta)),
##                                              integrate(function(x)pgamma(x,shape,rate=rate.0)*dgamma(x,shape,rate=rate.1),0,shape/rate.1+4*sqrt(shape)/rate.1)$val ) )
##     stop()
## }    
auc.full <- auc.lda.gamma(beta.star.full,rate.0,rate.1)
auc.red <- auc.lda.gamma(beta.star.red,rate.0,rate.1)
diff.coef <- auc.full - auc.red
by.n <- sapply(ns, function(n) {
    cat('.')    
    n.0 <- round(pi.0*n); n.1 <- n-n.0
    ## err <- 0
    stats <- replicate(3e2, {
        x.0.full <- matrix(rexp(n.0*p.full,rate.0),ncol=p.full)
        x.1.full <- matrix(rexp(n.1*p.full,rate.1),ncol=p.full)
        x.full <- rbind(x.0.full,x.1.full)
        g <- c(rep(0,n.0),rep(1,n.1))
        x.0.red <- x.0.full[,1:p.red]
        x.1.red <- x.1.full[,1:p.red]
        x.red <- rbind(x.0.red,x.1.red)
        Sigma.hat.full <- (  with(list(x.scaled=scale(x.0.full,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.full,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n  
        Sigma.hat.red <- (  with(list(x.scaled=scale(x.0.red,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.red,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n 
        params.full <- list(mu.0=colMeans(x.0.full),mu.1=colMeans(x.1.full),Sigma=Sigma.hat.full) 
        params.red <- list(mu.0=colMeans(x.0.red),mu.1=colMeans(x.1.red),Sigma=Sigma.hat.red) 
        beta.hat.full <- coefs.lda(x.0.full,x.1.full,params.full)
        beta.hat.red <- coefs.lda(x.0.red,x.1.red,params.red)
        ## if(sum(c(beta.hat.full,beta.hat.red)<0)>1) err <<- err+1
        ## beta.hat.full <- pmax(0,beta.hat.full) + rnorm(length(beta.hat.full))/n
        ## beta.hat.red <- pmax(0,beta.hat.red) + rnorm(length(beta.hat.red))/m
        beta.hat.full <- abs(beta.hat.full)
        beta.hat.red <- abs(beta.hat.red) 
        ## beta.hat.full <- abs(beta.star.full + rnorm(p.full)/sqrt(n))
        ## beta.hat.red <- abs(beta.star.red + rnorm(p.red)/sqrt(n))
        ## plot(ecdf(x.0.full%*%beta.star.full)) #1
        ## curve(F.0.full,add=TRUE,col=2)
        ## plot(ecdf(x.1.full%*%beta.star.full))
        ## curve(F.1.full,add=TRUE,col=2)
        ## plot(ecdf(x.0.full%*%beta.hat.full))
        ## curve(F.0.full,add=TRUE,col=2)
        ## plot(ecdf(x.0.red%*%beta.hat.red))
        ## curve(F.0.red,add=TRUE,col=2)
        ## plot(ecdf(x.1.red%*%beta.hat.red))
        ## curve(F.1.red,add=TRUE,col=2)
        auc.hat.full.hat <- auc.hat(x.0.full%*%beta.hat.full,x.1.full%*%beta.hat.full)
        auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
        diff.hat.coef.hat <- auc.hat.full.hat-auc.hat.red.hat
        auc.full.hat <- auc.lda.gamma(beta.hat.full,rate.0,rate.1)
        auc.red.hat <- auc.lda.gamma(beta.hat.red,rate.0,rate.1)
        diff.coef.hat <- auc.full.hat-auc.red.hat
        auc.hajek.full.hat <- auc.hajek(F=F.0.full,G=F.1.full,x=x.0.full%*%beta.hat.full,y=x.1.full%*%beta.hat.full,auc=auc.full.hat,terms.only=TRUE,IID=TRUE)
        auc.hajek.red.hat <- auc.hajek(F=F.0.red,G=F.1.red,x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red,auc=auc.red.hat,terms.only=TRUE,IID=TRUE)
        diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
        obs <- diff.hat.coef.hat-diff.coef.hat
        approx <- mean(diff.hajek.coef.hat)
        obs - approx
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)
## save.image('sessions/13c-5.RData')




## 13c-6 differentiating auc.lda.gamma

source('misc.R')
rate.0 <- runif(1)
rate.1 <- runif(1)
beta <- runif(3)
f <- function(beta)auc.lda.gamma(beta,rate.0,rate.1)
f.prime <- function(beta) numDeriv::grad(f,beta)
viz.deriv(f,f.prime,dim.x=length(beta),dim.f=1)


## 13c-7 #1 taylor expansion

require(sdprisk)
require(numDeriv)
source('misc.R')
## auc.gamma <- function(rate.0,rate.1,shape)integrate(function(x)pgamma(x,shape,rate=rate.0)*dgamma(x,shape,rate=rate.1),0,shape/rate.1+4*sqrt(shape)/rate.1)$val
set.seed(1)
## n <- 1e2
## n.0 <- n.1 <- round(n/2)
p.full <- 3
p.red <- 2
pi.0 <- 1/3
beta.full <-  rep(1,p.full)
beta.red <- rep(1,p.red)
rate.0 <- runif(1,1,3)
rate.1 <- runif(1,1,rate.0)
## beta.combined <- c(beta.full,beta.red)
## diff.coef <- auc.full-auc.red
Sigma.full <- diag(p.full)*(pi.0*1/rate.0^2 + (1-pi.0)*1/rate.1^2)
Sigma.red <- Sigma.full[1:p.red,1:p.red]
mu.diff.full <- rep(1/rate.1-1/rate.0,p.full)
## mu.diff.red <- rep(1/rate.1-1/rate.0,p.red)
beta.star.full <- as.numeric(solve(Sigma.full)%*%mu.diff.full)
beta.star.red <- beta.star.full[1:p.red]
F.0.full <- function(q)pgamma(q,shape=p.full,rate=rate.0/unique(beta.star.full))
F.1.full <- function(q)pgamma(q,shape=p.full,rate=rate.1/unique(beta.star.full))
F.0.red <- function(q)pgamma(q,shape=p.red,rate=rate.0/unique(beta.star.red))
F.1.red <- function(q)pgamma(q,shape=p.red,rate=rate.1/unique(beta.star.red))
ns <- round(seq(2e2,1e3,len=20))
## P(x<y) where x ~ sum_i beta_i x_i ~ hypoexp(rate.0/beta), y ~ sum_i
## beta_i y_i, x_i ~ Exp(rate.0) ~ hypoexp(rate.1/beta), x_i ~
## Exp(rate.0), y_i ~ Exp(rate.1).
auc.full <- auc.lda.gamma(beta.star.full,rate.0,rate.1)
auc.red <- auc.lda.gamma(beta.star.red,rate.0,rate.1)
diff.coef <- auc.full - auc.red
deriv.full <- grad(function(beta)auc.lda.gamma(beta,rate.0,rate.1),beta.star.full)
deriv.red <- grad(function(beta)auc.lda.gamma(beta,rate.0,rate.1),beta.star.red)
deriv <- c(deriv.full, -deriv.red)
## source('misc.R')
## f <- function(beta)auc.lda.gamma(beta,rate.0,rate.1)
## ## f.prime <- function(beta) numDeriv::grad(f,beta)
## ## viz.deriv(f,f.prime,dim.x=length(beta),dim.f=1,x0=beta.star.full)
## xs <- seq(-1,1,len=20)/10
## xs <- sort(c(1/seq(1e3,1e6,len=20)^1,0))
## a <- runif(2)
## aucs <- sapply(xs,function(x)f(beta.star.full+c(0,1,2)*x))
## plot(xs,aucs)
by.n <- sapply(ns, function(n) {
    cat('.')    
    n.0 <- round(pi.0*n); n.1 <- n-n.0
    ## err <- 0
    stats <- replicate(3e2, {
        x.0.full <- matrix(rexp(n.0*p.full,rate.0),ncol=p.full)
        x.1.full <- matrix(rexp(n.1*p.full,rate.1),ncol=p.full)
        x.full <- rbind(x.0.full,x.1.full)
        g <- c(rep(0,n.0),rep(1,n.1))
        x.0.red <- x.0.full[,1:p.red]
        x.1.red <- x.1.full[,1:p.red]
        x.red <- rbind(x.0.red,x.1.red)
        Sigma.hat.full <- (  with(list(x.scaled=scale(x.0.full,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.full,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n  
        Sigma.hat.red <- (  with(list(x.scaled=scale(x.0.red,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.red,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n 
        params.hat.full <- list(mu.0=colMeans(x.0.full),mu.1=colMeans(x.1.full),Sigma=Sigma.full) 
        params.hat.red <- list(mu.0=colMeans(x.0.red),mu.1=colMeans(x.1.red),Sigma=Sigma.red) 
        ## beta.hat.full <- beta.star.full+rnorm(p.full)/sqrt(n)#coefs.lda(x.0.full,x.1.full,params.full)
        ## beta.hat.red <- beta.star.red+rnorm(p.red)/sqrt(n)#coefs.lda(x.0.red,x.1.red,params.red)
        beta.hat.full <- coefs.lda(x.0.full,x.1.full,params.hat.full)
        beta.hat.red <- coefs.lda(x.0.red,x.1.red,params.hat.red)
        ## if(sum(c(beta.hat.full,beta.hat.red)<0)>1) err <<- err+1
        beta.hat.full <- abs(beta.hat.full)
        beta.hat.red <- abs(beta.hat.red)
        auc.full.hat <- auc.lda.gamma(beta.hat.full,rate.0,rate.1)
        auc.red.hat <- auc.lda.gamma(beta.hat.red,rate.0,rate.1)
        diff.coef.hat <- auc.full.hat-auc.red.hat
        obs <- diff.coef.hat - diff.coef
        ## approx <- deriv.full%*%(beta.hat.full - beta.star.full) - deriv.red%*%(beta.hat.red - beta.star.red)
        ## approx <- deriv %*% rbind(beta.hat.full - beta.star.full, beta.hat.red - beta.star.red)
        infl.full <- infl.lda(x.full,g,params.hat.full,terms.only=TRUE) 
        infl.red <- infl.lda(x.red,g,params.hat.red,terms.only=TRUE) 
        infl <- rbind(infl.full,infl.red)
        approx <- mean(deriv %*% infl)
        obs - approx
    })
    ## print(err)
    stats
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)
## save.image('sessions/13c-5.RData')



## 13c-8 combining hajek and taylor parts

require(sdprisk)
require(numDeriv)
source('misc.R')
set.seed(1)
## n <- 1e2
## n.0 <- n.1 <- round(n/2)
p.full <- 3
p.red <- 2
pi.0 <- 1/3
beta.full <-  rep(1,p.full)
beta.red <- rep(1,p.red)
rate.0 <- runif(1,1,3)
rate.1 <- runif(1,1,rate.0)
Sigma.full <- diag(p.full)*(pi.0*1/rate.0^2 + (1-pi.0)*1/rate.1^2)
Sigma.red <- Sigma.full[1:p.red,1:p.red]
mu.diff.full <- rep(1/rate.1-1/rate.0,p.full)
## mu.diff.red <- rep(1/rate.1-1/rate.0,p.red)
beta.star.full <- as.numeric(solve(Sigma.full)%*%mu.diff.full)
beta.star.red <- beta.star.full[1:p.red]
F.0.full <- function(q)pgamma(q,shape=p.full,rate=rate.0/unique(beta.star.full))
F.1.full <- function(q)pgamma(q,shape=p.full,rate=rate.1/unique(beta.star.full))
F.0.red <- function(q)pgamma(q,shape=p.red,rate=rate.0/unique(beta.star.red))
F.1.red <- function(q)pgamma(q,shape=p.red,rate=rate.1/unique(beta.star.red))
ns <- round(seq(2e2,1e3,len=20))
auc.full <- auc.lda.gamma(beta.star.full,rate.0,rate.1)
auc.red <- auc.lda.gamma(beta.star.red,rate.0,rate.1)
diff.coef <- auc.full - auc.red
deriv.full <- grad(function(beta)auc.lda.gamma(beta,rate.0,rate.1),beta.star.full)
deriv.red <- grad(function(beta)auc.lda.gamma(beta,rate.0,rate.1),beta.star.red)
deriv <- c(deriv.full, -deriv.red)
by.n <- sapply(ns, function(n) {
    cat('.')    
    n.0 <- round(pi.0*n); n.1 <- n-n.0
    ## err <- 0
    stats <- replicate(3e2, {
        x.0.full <- matrix(rexp(n.0*p.full,rate.0),ncol=p.full)
        x.1.full <- matrix(rexp(n.1*p.full,rate.1),ncol=p.full)
        x.full <- rbind(x.0.full,x.1.full)
        g <- c(rep(0,n.0),rep(1,n.1))
        x.0.red <- x.0.full[,1:p.red]
        x.1.red <- x.1.full[,1:p.red]
        x.red <- rbind(x.0.red,x.1.red)
        Sigma.hat.full <- (  with(list(x.scaled=scale(x.0.full,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.full,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n  
        Sigma.hat.red <- (  with(list(x.scaled=scale(x.0.red,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.red,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n 
        params.hat.full <- list(mu.0=colMeans(x.0.full),mu.1=colMeans(x.1.full),Sigma=Sigma.full) 
        params.hat.red <- list(mu.0=colMeans(x.0.red),mu.1=colMeans(x.1.red),Sigma=Sigma.red) 
        beta.hat.full <- coefs.lda(x.0.full,x.1.full,params.hat.full)
        beta.hat.red <- coefs.lda(x.0.red,x.1.red,params.hat.red)
        beta.hat.full <- abs(beta.hat.full)
        beta.hat.red <- abs(beta.hat.red)
        auc.full.hat <- auc.lda.gamma(beta.hat.full,rate.0,rate.1)
        auc.red.hat <- auc.lda.gamma(beta.hat.red,rate.0,rate.1)
        diff.coef.hat <- auc.full.hat-auc.red.hat
        infl.full <- infl.lda(x.full,g,params.hat.full,terms.only=TRUE) 
        infl.red <- infl.lda(x.red,g,params.hat.red,terms.only=TRUE) 
        infl <- rbind(infl.full,infl.red)
        ## obs <- diff.coef.hat - diff.coef
        ## approx <- mean(deriv %*% infl)
        ## obs - approx
        auc.hat.full.hat <- auc.hat(x.0.full%*%beta.hat.full,x.1.full%*%beta.hat.full)
        auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
        diff.hat.coef.hat <- auc.hat.full.hat-auc.hat.red.hat
        ## auc.full.hat <- auc.lda.gamma(beta.hat.full,rate.0,rate.1)
        ## auc.red.hat <- auc.lda.gamma(beta.hat.red,rate.0,rate.1)
        ## diff.coef.hat <- auc.full.hat-auc.red.hat
        auc.hajek.full.hat <- auc.hajek(F=F.0.full,G=F.1.full,x=x.0.full%*%beta.hat.full,y=x.1.full%*%beta.hat.full,auc=auc.full.hat,terms.only=TRUE,IID=TRUE)
        auc.hajek.red.hat <- auc.hajek(F=F.0.red,G=F.1.red,x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red,auc=auc.red.hat,terms.only=TRUE,IID=TRUE)
        diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
        ## obs <- diff.hat.coef.hat-diff.coef.hat
        ## approx <- mean(diff.hajek.coef.hat)
        ## obs - approx
        obs <- diff.hat.coef.hat - diff.coef
        approx <- mean(diff.hajek.coef.hat + as.numeric(deriv %*% infl))
        obs - approx
    })
    ## print(err)
    stats
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)
## save.image('sessions/13c-5.RData')




## 13c-9 #1 checking z-stats. #2 with estimated hajek. #3 w estimated
## influence function #4 estimated deriv at non estimated beta.star. a
## little to light tailed/concentrated at 0. #4 all estimated
## params. light tails again. but improving with n.

## update: looks more normal when omitting derivative adjustment?
## actually no. serious outliers (~40 sds) not going away with larger
## n even getting worse. looked normal because of the scale on the
## vertical that results. since derivative is 0 though (see notes) the
## outliers are unexpected when the deriv approx is set to 0: indeed
## set the derivatives to 0 and comment out the derivative
## calculations in the loop and set n > 1e4 to see quite
## normal-looking z-stats.

require(parallel)
require(sdprisk)
require(numDeriv)
source('misc.R')
start <- Sys.time()
set.seed(2)
n <- 1.3e4
## ns <- round(seq(2e2,1e3,len=20))
## n.0 <- n.1 <- round(n/2)
p.full <- 3
p.red <- 2
pi.0 <- 1/3
beta.full <-  rep(1,p.full)
beta.red <- rep(1,p.red)
rate.0 <- runif(1,1,3)
rate.1 <- runif(1,1,rate.0)
Sigma.full <- diag(p.full)*(pi.0*1/rate.0^2 + (1-pi.0)*1/rate.1^2)
Sigma.red <- Sigma.full[1:p.red,1:p.red]
mu.diff.full <- rep(1/rate.1-1/rate.0,p.full)
## mu.diff.red <- rep(1/rate.1-1/rate.0,p.red)
beta.star.full <- as.numeric(solve(Sigma.full)%*%mu.diff.full)
beta.star.red <- beta.star.full[1:p.red]
## F.0.full <- function(q)pgamma(q,shape=p.full,rate=rate.0/unique(beta.star.full)) #1
## F.1.full <- function(q)pgamma(q,shape=p.full,rate=rate.1/unique(beta.star.full)) #1
## F.0.red <- function(q)pgamma(q,shape=p.red,rate=rate.0/unique(beta.star.red)) #1
## F.1.red <- function(q)pgamma(q,shape=p.red,rate=rate.1/unique(beta.star.red)) #1
auc.full <- auc.lda.gamma(beta.star.full,rate.0,rate.1)
auc.red <- auc.lda.gamma(beta.star.red,rate.0,rate.1)
diff.coef <- auc.full - auc.red
## deriv.full <- grad(function(beta)auc.lda.gamma(beta,rate.0,rate.1),beta.star.full) #1,2,3
## deriv.red <- grad(function(beta)auc.lda.gamma(beta,rate.0,rate.1),beta.star.red) #1,2,3
## deriv <- c(deriv.full, -deriv.red) #1,2,3
## by.n <- sapply(ns, function(n) {
## cat('.')    
n.0 <- round(pi.0*n); n.1 <- n-n.0
## err <- 0
z.stats <- mclapply(1:3e2, mc.cores=detectCores()-3, FUN=function(dd) {
## z.stats <- replicate(3e2, {
    x.0.full <- matrix(rexp(n.0*p.full,rate.0),ncol=p.full)
    x.1.full <- matrix(rexp(n.1*p.full,rate.1),ncol=p.full)
    x.full <- rbind(x.0.full,x.1.full)
    g <- c(rep(0,n.0),rep(1,n.1))
    x.0.red <- x.0.full[,1:p.red]
    x.1.red <- x.1.full[,1:p.red]
    x.red <- rbind(x.0.red,x.1.red)
    Sigma.hat.full <- (  with(list(x.scaled=scale(x.0.full,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.full,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n  
    Sigma.hat.red <- (  with(list(x.scaled=scale(x.0.red,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1.red,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n 
    ## params.hat.full <- list(mu.0=colMeans(x.0.full),mu.1=colMeans(x.1.full),Sigma=Sigma.full) #1,#2 
    ## params.hat.red <- list(mu.0=colMeans(x.0.red),mu.1=colMeans(x.1.red),Sigma=Sigma.red) #1,#2
    params.hat.full <- list(mu.0=colMeans(x.0.full),mu.1=colMeans(x.1.full),Sigma=Sigma.hat.full) #3
    params.hat.red <- list(mu.0=colMeans(x.0.red),mu.1=colMeans(x.1.red),Sigma=Sigma.hat.red) #3
    beta.hat.full <- coefs.lda(x.0.full,x.1.full,params.hat.full)
    beta.hat.red <- coefs.lda(x.0.red,x.1.red,params.hat.red)
    ## beta.hat.full <- abs(beta.hat.full)
    ## beta.hat.red <- abs(beta.hat.red)
    ## auc.full.hat <- auc.lda.gamma(beta.hat.full,rate.0,rate.1)
    ## auc.red.hat <- auc.lda.gamma(beta.hat.red,rate.0,rate.1)
    ## diff.coef.hat <- auc.full.hat-auc.red.hat
    infl.full <- infl.lda(x.full,g,params.hat.full,terms.only=TRUE) 
    infl.red <- infl.lda(x.red,g,params.hat.red,terms.only=TRUE) 
    infl <- rbind(infl.full,infl.red)
    deriv.full <- grad(function(beta)auc.hat(x.0.full%*%beta,x.1.full%*%beta),beta.hat.full)#!!!!!!!!!!!
    deriv.red <- grad(function(beta)auc.hat(x.0.red%*%beta,x.1.red%*%beta),beta.hat.red)
    deriv <- c(deriv.full, -deriv.red)
    ## obs <- diff.coef.hat - diff.coef
    ## approx <- mean(deriv %*% infl)
    ## obs - approx
    auc.hat.full.hat <- auc.hat(x.0.full%*%beta.hat.full,x.1.full%*%beta.hat.full)
    auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
    diff.hat.coef.hat <- auc.hat.full.hat-auc.hat.red.hat
    ## auc.hajek.full.hat <- auc.hajek(F=F.0.full,G=F.1.full,x=x.0.full%*%beta.hat.full,y=x.1.full%*%beta.hat.full,auc=auc.full.hat,terms.only=TRUE,IID=TRUE) #1
    ## auc.hajek.red.hat <- auc.hajek(F=F.0.red,G=F.1.red,x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red,auc=auc.red.hat,terms.only=TRUE,IID=TRUE) #1
    auc.hajek.full.hat <- auc.hajek(F=NULL,G=NULL,x=x.0.full%*%beta.hat.full,y=x.1.full%*%beta.hat.full,auc=auc.hat.full.hat,terms.only=TRUE,IID=TRUE) #2
    auc.hajek.red.hat <- auc.hajek(F=NULL,G=NULL,x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red,auc=auc.hat.red.hat,terms.only=TRUE,IID=TRUE) #2
    diff.hajek.coef.hat <- auc.hajek.full.hat - auc.hajek.red.hat
    ## obs <- diff.hat.coef.hat-diff.coef.hat
    ## approx <- mean(diff.hajek.coef.hat)
    ## obs - approx
    obs <- diff.hat.coef.hat - diff.coef
    approx <- diff.hajek.coef.hat + as.numeric(deriv %*% infl)
    obs / sqrt(var(approx) / length(approx))
})
Sys.time() - start
z.stats <- unlist(z.stats)
qqnorm(z.stats)
abline(0,1)

hist(z.stats,prob=TRUE)
curve(dnorm,add=TRUE)
## save.image('sessions/13c-9.RData')


r0 <- rate.0
r1 <- rate.1
pi.1 <- 1-pi.0
beta.star.full
(r1-r0)/(pi.0*r1/r0+pi.1*r0/r1)
1/(r1/r0/(r1-r0) - pi.1*(1/r1+1/r0))
r0*r1/(r1+r0) / (r1^2/(r1^2-r0^2) -  pi.1)

pairs <- replicate(1e3, {
r0 <- runif(1)
r1 <- runif(1)
c((r1-r0)/(pi.0*r1/r0+pi.1*r0/r1),  1/r0-1/r1)
})
plot(pairs[1,],pairs[2,],xlim=c(-.01,.01))






## 14 general glm + gaussian

## 14a plim of beta.hat under misspecified binary response model

require(mvtnorm)
p <- 4
n <- 1e3
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta.0 <- runif(p)
beta <- runif(p)
x <- rmvnorm(n,mu,Sigma)
obs <- h(x%*%beta.0)*log(h(x%*%beta))
w <- runif(1)
cond.mean <- function(w) beta.0%*%mu + t(beta.0)%*%Sigma%*%beta / (t(beta)%*%Sigma%*%beta) * (w- t(beta)%*%mu)
cond.var <- (1- (t(beta.0)%*%Sigma%*%beta)^2 / (t(beta)%*%Sigma%*%beta) / (t(beta.0)%*%Sigma%*%beta.0)) * t(beta.0)%*%Sigma%*%beta.0
inner <- Vectorize(function(w)integrate(function(w.0)h(w.0)*dnorm(w.0,cond.mean(w),sqrt(cond.var)),cond.mean(w) - 4*sqrt(cond.var), cond.mean(w)+4*sqrt(cond.var))$val)
try <- integrate(function(w)log(h(w))*dnorm(w,t(beta)%*%mu,sqrt(t(beta)%*%Sigma%*%beta))*inner(w),t(beta)%*%mu - 5*sqrt(t(beta)%*%Sigma%*%beta), t(beta)%*%mu + 5*sqrt(t(beta)%*%Sigma%*%beta))$val
hist(obs)
abline(v=try,col=2)
abline(v=mean(obs),col=3)



require(mvtnorm)
set.seed(1)
p <- 4
n <- 1e3
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta.0 <- runif(p)/5
lim <- 5
cond.mean <- function(w,beta) beta.0%*%mu + t(beta.0)%*%Sigma%*%beta / (t(beta)%*%Sigma%*%beta) * (w- t(beta)%*%mu)
cond.var <- function(beta)(1- (t(beta.0)%*%Sigma%*%beta)^2 / (t(beta)%*%Sigma%*%beta) / (t(beta.0)%*%Sigma%*%beta.0)) * t(beta.0)%*%Sigma%*%beta.0
inner <- Vectorize(function(w,beta)integrate(function(w.0)  (h(w.0)*log(h(w)) + (1-h(w.0))*log(1-h(w)))    *dnorm(w,t(beta)%*%mu,sqrt(t(beta)%*%Sigma%*%beta))*dnorm(w.0,cond.mean(w,beta),sqrt(cond.var(beta))),-lim,lim)$val,vectorize.args='w')
criterion <- function(beta)integrate(function(w)inner(w,beta),-lim,lim,subdivisions=1e3)$val
beta <- runif(p)
x <- rmvnorm(n,mu,Sigma)
obs <- h(x%*%beta.0)*log(h(x%*%beta)) + (1-h(x%*%beta.0))*log(1-h(x%*%beta))
hist(obs)
abline(v=criterion(beta),col=2)
abline(v=mean(obs),col=3)

## alternate form of criterion
require(mvtnorm)
## set.seed(1)
p <- 4
n <- 1e3
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta.0 <- runif(p)/5
beta <- runif(p)/2
lim <- 5
criterion.old <- function(beta){
    cond.mean <- function(w,beta) beta.0%*%mu + t(beta.0)%*%Sigma%*%beta / (t(beta)%*%Sigma%*%beta) * (w- t(beta)%*%mu)
    cond.var <- function(beta)(1- (t(beta.0)%*%Sigma%*%beta)^2 / (t(beta)%*%Sigma%*%beta) / (t(beta.0)%*%Sigma%*%beta.0)) * t(beta.0)%*%Sigma%*%beta.0
    inner <- Vectorize(function(w,beta)integrate(function(w.0)  (h(w.0)*log(h(w)) + (1-h(w.0))*log(1-h(w)))    *dnorm(w,t(beta)%*%mu,sqrt(t(beta)%*%Sigma%*%beta))*dnorm(w.0,cond.mean(w,beta),sqrt(cond.var(beta))),-lim,lim)$val,vectorize.args='w')
    integrate(function(w)inner(w,beta),-lim,lim,subdivisions=1e3)$val
}
criterion.old(beta)
criterion <- function(beta){
    cond.mean <- function(w,beta) beta.0%*%mu + t(beta.0)%*%Sigma%*%beta / (t(beta)%*%Sigma%*%beta) * (w- t(beta)%*%mu)
    cond.var <- function(beta)(1- (t(beta.0)%*%Sigma%*%beta)^2 / (t(beta)%*%Sigma%*%beta) / (t(beta.0)%*%Sigma%*%beta.0)) * t(beta.0)%*%Sigma%*%beta.0
    f <- function(w,beta)dnorm(w,t(beta)%*%mu,sqrt(t(beta)%*%Sigma%*%beta))
    f.cond <- function(w,w.0,beta)dnorm(w.0,cond.mean(w,beta),sqrt(cond.var(beta)))
    inner <- Vectorize(function(w,beta)integrate(function(w.0)w*f(w,beta)*h(w.0)*f.cond(w,w.0,beta),-lim,lim)$val,vectorize.args='w')
    integrate(function(w)inner(w,beta),-lim,lim,subdivisions=1e3)$val + integrate(function(w)log(1-h(w))*f(w,beta),-lim,lim)$val
    }
criterion(beta)
pairs <- replicate(1e1, {
    beta <- runif(p)/5
    c(criterion.old(beta),criterion(beta))
})
plot(pairs[1,],pairs[2,]); abline(0,1)



require(mvtnorm)
## set.seed(1)
p <- 4
n <- 1e3
mu <- runif(p)*0
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta.0 <- runif(p)/5
beta <- runif(p)/2
lim <- 10
cond.mean <- function(w,beta) beta.0%*%mu + t(beta.0)%*%Sigma%*%beta / (t(beta)%*%Sigma%*%beta) * (w- t(beta)%*%mu)
cond.var <- function(beta)(1- (t(beta.0)%*%Sigma%*%beta)^2 / (t(beta)%*%Sigma%*%beta) / (t(beta.0)%*%Sigma%*%beta.0)) * t(beta.0)%*%Sigma%*%beta.0
f <- function(w,beta)dnorm(w,t(beta)%*%mu,sqrt(t(beta)%*%Sigma%*%beta))
f.cond <- function(w,w.0,beta)dnorm(w.0,cond.mean(w,beta),sqrt(cond.var(beta)))
inner <- Vectorize(function(w,beta)integrate(function(w.0)w*f(w,beta)*h(w.0)*f.cond(w,w.0,beta),-lim,lim)$val,vectorize.args='w')
c(integrate(function(w)inner(w,beta),-lim,lim,subdivisions=1e3)$val, + integrate(function(w)log(1-h(w))*f(w,beta),-lim,lim)$val)


p.red <- 2
beta.star.red <- optim(runif(p),function(beta)-criterion(c(beta[1:p.red],rep(0,p-p.red))))$par[1:p.red]

n <- 1e3
a <- runif(p.red)
obs <- replicate(5e2, {
    x <- rmvnorm(n,mu,Sigma)
    x.red <- x[,1:p.red]
    risk <- h(x%*%beta.0)
    g <- rbinom(n,1,risk)
    ## beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link='probit')))
    beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link='logit')))
    beta.hat.red%*%a
})
hist(obs)
abline(v=beta.star.red%*%a,col=2)
abline(v=mean(obs),col=3)


## encapsulate
require(mvtnorm)
coef.reduced.glm.gaussian <- function(p.red,params,lim=Inf) {
    h <- params$link
    mu <- params$mu
    Sigma <- params$Sigma
    beta.0 <- params$beta
    p.full <- length(beta.0)
    if(p.red==p.full)return(beta.0)   
    cond.mean <- function(w,beta) beta.0%*%mu + t(beta.0)%*%Sigma%*%beta / (t(beta)%*%Sigma%*%beta) * (w- t(beta)%*%mu)
    cond.var <- function(beta)(1- (t(beta.0)%*%Sigma%*%beta)^2 / (t(beta)%*%Sigma%*%beta) / (t(beta.0)%*%Sigma%*%beta.0)) * t(beta.0)%*%Sigma%*%beta.0
    inner <- Vectorize(function(w,beta)integrate(function(w.0)  (h(w.0)*log(h(w)) + (1-h(w.0))*log(1-h(w)))    *dnorm(w,t(beta)%*%mu,sqrt(t(beta)%*%Sigma%*%beta))*dnorm(w.0,cond.mean(w,beta),sqrt(cond.var(beta))),-lim,lim)$val,vectorize.args='w')
    criterion <- function(beta)integrate(function(w)inner(w,beta),-lim,lim,subdivisions=1e3)$val
    optim(runif(p.full),function(beta)-criterion(c(beta[1:p.red],rep(0,p.full-p.red))))$par[1:p.red]
}
## source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 1e3
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
beta.0 <- runif(p.full)/2
h <- plogis
beta.star.red <- coef.reduced.glm.gaussian(p.red=2,params=list(mu=mu,Sigma=Sigma,beta=beta.0,link=h),lim=10)

n <- 1e3
a <- runif(p.red)
obs <- replicate(1e3, {
    x <- rmvnorm(n,mu,Sigma)
    x.red <- x[,1:p.red]
    risk <- h(x%*%beta.0)
    g <- rbinom(n,1,risk)
    ## beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link='probit')))
    beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link='logit')))
    beta.hat.red%*%a
})
hist(obs)
abline(v=beta.star.red%*%a,col=2)
abline(v=mean(obs),col=3)

## check on adversarial model
set.seed(1)
p.full <- 4
p.red <- 2
n <- 1e3
a <- 5
b <- 5
mu <- rep(0,p.full)
Sigma <- diag(p.full)*sqrt(a)
beta.0 <- c(rep(0,p.red),rep(1,p.full-p.red))/a 
h <- plogis
beta.star.red <- coef.reduced.glm.gaussian(p.red=2,params=list(mu=mu,Sigma=Sigma,beta=beta.0,link=h),lim=10) ## should be 0
n <- 1e3
pairs <- replicate(1e2, {
    a <- runif(p.red)
    obs <- replicate(1e2, {
        x <- rmvnorm(n,mu,Sigma)
        x.red <- x[,1:p.red]
        risk <- h(x%*%beta.0)
        g <- rbinom(n,1,risk)
        ## beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link='probit')))
        beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link='logit')))
        beta.hat.red%*%a
    })
    c(beta.star.red%*%a,mean(obs))
})
plot(pairs[1,],pairs[2,]); abline(0,1)
plot(pairs[1,] - pairs[2,]); abline(h=0)

## 14b pdf of beta^x | g=i assuming g satisfies P(g=1)=h(beta.0%*%x), and auc calculation
require(mvtnorm)
p.full <- 4
p.red <- 2
n <- 1e3
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta.0 <- runif(p.full)/2
beta <- runif(p.full)
f.w <- function(w)dnorm(w,t(beta)%*%mu, sd=sqrt(t(beta)%*%Sigma%*%beta))
f.w0 <- function(w0)dnorm(w0,t(beta.0)%*%mu, sd=sqrt(t(beta.0)%*%Sigma%*%beta.0))
E.g <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
x <- rmvnorm(n,mu,Sigma)
risk <- h(x%*%beta.0)
g <- rbinom(n,1,risk)
## mean(g)
## E.g
lim <- Inf
quad.ww <- as.numeric(t(beta)%*%Sigma%*%beta)
quad.ww0 <- as.numeric(t(beta.0)%*%Sigma%*%beta)
quad.w0w0 <- as.numeric(t(beta.0)%*%Sigma%*%beta.0)
f.w0.cond <- function(w0,w)dnorm(w0, mean=t(beta.0)%*%mu + quad.ww0/quad.ww*(w-t(beta)%*%mu), sd=sqrt((1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0))
f.w.w0.cond <- function(w,w0,g) h(w0)^g*(1-h(w0))^(1-g)*f.w0.cond(w0,w)*f.w(w)/(E.g^g*(1-E.g)^(1-g))
f.w.cond <- Vectorize( function(w,g) integrate(function(w0)f.w.w0.cond(w,w0,g),-lim,lim)$val, vectorize.args='w')
g0 <- sample(0:1,1)
hist((x%*%beta)[g==g0],prob=TRUE)
curve(f.w.cond(x,g0),add=TRUE)
integrate(function(x)f.w.cond(x,g0),-10,10)


## using logitnorm package
require(logitnorm)
w <- runif(1)
g <- rbinom(1,1,.5)
integrate(function(w0)f.w.w0.cond(w,w0,g),-lim,lim)$val
integrate(function(w0)h(w0)^g*(1-h(w0))^(1-g)*f.w0.cond(w0,w),-lim,lim)$val*f.w(w)/(E.g^g*(1-E.g)^(1-g))
integrate(function(w0)h(w0)*f.w0.cond(w0,w),-lim,lim)$val
ln.mean <- unname(momentsLogitnorm(mu=t(beta.0)%*%mu + quad.ww0/quad.ww*(w-t(beta)%*%mu), sigma=sqrt((1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0))['mean'])
ln.mean^g*(1-ln.mean)^(1-g)*f.w(w)/(E.g^g*(1-E.g)^(1-g))



## auc computation using logitnorm mean
ln.mean <- function(x,beta,params,lim=Inf) {
    mu <- params$mu; Sigma <- params$Sigma; beta.0 <- params$beta
    quad.ww <- as.numeric(t(beta)%*%Sigma%*%beta)
    quad.ww0 <- as.numeric(t(beta.0)%*%Sigma%*%beta)
    quad.w0w0 <- as.numeric(t(beta.0)%*%Sigma%*%beta.0)
    Sigma.cond <- sqrt((1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0) 
    ln.mean <-  sapply(x, function(x)integrate(function(w0)plogis(w0)*dnorm(w0,mean=sum(beta.0*mu) + quad.ww0/quad.ww*(x-sum(beta*mu)), sd=Sigma.cond),-lim,lim)$val)
}
## ln.mean^g*(1-ln.mean)^(1-g)*f.w(w)/(E.g^g*(1-E.g)^(1-g))

pdf.index.glm.gaussian.new <- function(w,g,beta,params,lim=Inf) {
    ## browser()
    h <- params$link
    f.w <- function(w)dnorm(w,t(beta)%*%mu, sd=sqrt(t(beta)%*%Sigma%*%beta))
    f.w0 <- function(w0)dnorm(w0,t(beta.0)%*%mu, sd=sqrt(t(beta.0)%*%Sigma%*%beta.0))
    E.g <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    ln.means <- ln.mean(w,beta,params,lim)
    ln.means^g*(1-ln.means)^(1-g)*f.w(w)/(E.g^g*(1-E.g)^(1-g))
}


## check new pdf routine using logitnorm is the same
require(mvtnorm)
source('misc.R')
set.seed(1)
p <- 4
## p.red <- 2
n <- 1e3
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta.0 <- runif(p)/2
beta <- runif(p)
## beta <- beta.0
params <- list(mu=mu,Sigma=Sigma,beta=beta.0)
## x <- rmvnorm(n,mu,Sigma)
## risk <- h(x%*%beta.0)
## g <- rbinom(n,1,risk)
g0 <- sample(0:1,1)
## hist((x%*%beta)[g==g0],prob=TRUE)
curve(pdf.index.glm.gaussian(x,g0,beta,params),-2,2)
curve(pdf.index.glm.gaussian.new(x,g0,beta,params),add=TRUE,col=2,lty=2)


## ## save(w,g,beta,mu,Sigma,beta.0,file='sessions/14b.RData')
## load('sessions/14b.RData')
## require(logitnorm)
## h <- plogis
## lim <- Inf
## f.w <- function(w)dnorm(w,t(beta)%*%mu, sd=sqrt(t(beta)%*%Sigma%*%beta))
## f.w0 <- function(w0)dnorm(w0,t(beta.0)%*%mu, sd=sqrt(t(beta.0)%*%Sigma%*%beta.0))
## quad.ww <- as.numeric(t(beta)%*%Sigma%*%beta)
## quad.ww0 <- as.numeric(t(beta.0)%*%Sigma%*%beta)
## quad.w0w0 <- as.numeric(t(beta.0)%*%Sigma%*%beta.0)
## f.w0.cond <- function(w0,w)dnorm(w0, mean=t(beta.0)%*%mu + quad.ww0/quad.ww*(w-t(beta)%*%mu), sd=sqrt((1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0))
## integrate(function(w0)h(w0)*f.w0.cond(w0,w),-lim,lim)$val
## momentsLogitnorm(mu=t(beta.0)%*%mu + quad.ww0/quad.ww*(w-t(beta)%*%mu), sigma=sqrt((1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0))['mean']


## encapsulate
require(mvtnorm)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 1e3
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta.0 <- runif(p.full)/2
beta <- runif(p.full)
## beta <- beta.0
params <- list(mu=mu,Sigma=Sigma,beta=beta.0)
## f.w <- function(w)dnorm(w,t(beta)%*%mu, sd=sqrt(t(beta)%*%Sigma%*%beta))
## f.w0 <- function(w0)dnorm(w0,t(beta.0)%*%mu, sd=sqrt(t(beta.0)%*%Sigma%*%beta.0))
## E.g <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
x <- rmvnorm(n,mu,Sigma)
risk <- h(x%*%beta.0)
g <- rbinom(n,1,risk)
g0 <- sample(0:1,1)
hist((x%*%beta)[g==g0],prob=TRUE)
curve(pdf.index.glm.gaussian(x,g0,beta,params),add=TRUE)
integrate(function(x)pdf.index.glm.gaussian(x,g0,beta,params),-10,10)

## make sure logitnorm version gives the same results
summary(as.numeric(pdf.index.glm.gaussian.new(x%*%beta,g0,beta,params)-pdf.index.glm.gaussian(x%*%beta,g0,beta,params)))

## newer version takes about a 1/3 the time
set.seed(1)
n <- 1e4
x <- rnorm(n)
start <- Sys.time()
invisible(pdf.index.glm.gaussian.new(x,g0,beta,params))
print(Sys.time() - start)
start <- Sys.time()
invisible(pdf.index.glm.gaussian(x,g0,beta,params))
print(Sys.time() - start)

## 14c auc = P(beta^tx0 < beta^tx1 | g0=0,g1=1)

require(mvtnorm)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 1e2
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- pnorm
beta.0 <- runif(p.full)/2
beta <- runif(p.full)
## beta <- beta.0
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,h=h)
f.0 <- function(x)pdf.index.glm.gaussian(x,g=0,beta=beta.0,params=params)
f.1 <- function(x)pdf.index.glm.gaussian(x,g=1,beta=beta.0,params=params)
F.0 <- Vectorize(function(x)integrate(function(w)f.0(w),-Inf,x)$val)
## x <- rmvnorm(n,mu,Sigma)
## risk <- h(x%*%beta.0)
## g <- rbinom(n,1,risk)
## plot(ecdf((x%*%beta.0)[g==0]))
## curve(F.0(x),add=TRUE,col=2)
auc <- integrate(function(x)F.0(x)*f.1(x),-Inf,Inf)$val
## auc <- auc.glm.gaussian(beta.0,params)
obs <- replicate(1e3, {
    x <- rmvnorm(n,mu,Sigma)
    risk <- h(x%*%beta.0)
    g <- rbinom(n,1,risk)
    auc.hat(x[g==0,]%*%beta.0,x[g==1,]%*%beta.0)
})
hist(obs)
abline(v=auc,col=2)
abline(v=mean(obs),col=3)


## auc on reduced coefs

require(mvtnorm)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 3e2
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- pnorm
beta.0 <- runif(p.full)/2
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,h=h)
beta.star.red <- coef.reduced.glm.gaussian(p.red=2,params,lim=5)
## beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
## beta <- beta.0
f.0 <- function(x)pdf.index.glm.gaussian(x,g=0,beta=beta.star.red,params=params)
f.1 <- function(x)pdf.index.glm.gaussian(x,g=1,beta=beta.star.red,params=params)
F.0 <- Vectorize(function(x)integrate(function(w)f.0(w),-Inf,x)$val)
## x <- rmvnorm(n,mu,Sigma)
## risk <- h(x%*%beta.0)
## g <- rbinom(n,1,risk)
## plot(ecdf((x%*%beta.star.red)[g==0]))
## curve(F.0(x),add=TRUE,col=2)
auc <- integrate(function(x)F.0(x)*f.1(x),-Inf,Inf)$val
obs <- replicate(1e3, {
    x <- rmvnorm(n,mu,Sigma)
    risk <- h(x%*%beta.0)
    g <- rbinom(n,1,risk)
    auc.hat(x[g==0,]%*%beta.star.red,x[g==1,]%*%beta.star.red)
})
hist(obs)
abline(v=auc,col=2)
abline(v=mean(obs),col=3)


## encapsulate

require(mvtnorm)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 3e2
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta.0 <- runif(p.full)/2
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h)
## beta <- beta.0
beta <- runif(p.full)
beta <- coef.reduced.glm.gaussian(p.red=2,params,lim=6)
## beta <- c(beta,rep(0,p.full-p.red))
obs <- replicate(3e2, {
    x <- rmvnorm(n,mu,Sigma)
    risk <- h(x%*%beta.0)
    g <- rbinom(n,1,risk)
    ## auc.hat(x[g==0,]%*%beta,x[g==1,]%*%beta)
    auc.hat(x[g==0,1:p.red]%*%beta,x[g==1,1:p.red]%*%beta)
})
auc.try <- auc.glm.gaussian(beta,params)
## auc.try <- auc.new(beta,params,5)
hist(obs)
abline(v=auc.try,col=2)
abline(v=mean(obs),col=3)


## verify grad is 0 at true beta.0 / true risk
numDeriv::grad(function(beta)auc.glm.gaussian(beta,params),x=beta.0)


start <- Sys.time()
auc.try <- auc.glm.gaussian(beta,params)
print(Sys.time() - start)


## check  auc routine using logitnorm is the same (assumes h=logit)
p <- 4
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta.0 <- runif(p)/2
beta <- runif(p)
## beta <- beta.0
params <- list(mu=mu,Sigma=Sigma,beta=beta.0)
auc.glm.gaussian.new <- function(beta,params,lim=Inf) {
    mu <- params$mu; Sigma <- params$Sigma; beta.0 <- params$beta
    f.w <- function(w)dnorm(w,t(beta)%*%mu, sd=sqrt(t(beta)%*%Sigma%*%beta))
    f.w0 <- function(w0)dnorm(w0,t(beta.0)%*%mu, sd=sqrt(t(beta.0)%*%Sigma%*%beta.0))
    E.g <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    inner <- Vectorize(function(v,beta)integrate(function(w)(1-ln.mean(w,beta,params))*ln.mean(v,beta,params)*f.w(w)*f.w(v),-lim,v)$val,vectorize.args='v')
    integrate(function(v)inner(v,beta),-lim,lim)$val / (E.g*(1-E.g))
}
c(auc.glm.gaussian.new(beta,params),     auc.glm.gaussian(beta,params))



## 14d hajek. seems O(1/n) like expected but routine is very slow so
## not much data. and many NAs. [update: redo below in 16f]

start <- Sys.time()
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 3e2
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta.0 <- runif(p.full)/2
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link.name='logit')
beta.star.full <- coef.reduced.glm.gaussian(p.full,h,params)
## beta.star.red <- coef.reduced.glm.gaussian(p.red=2,link=h,params)
## beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
## beta <- beta.0
lim <- 5
ns <- round(seq(10,2e2,len=40))
auc.full <- auc.glm.gaussian(beta.star.full,h,params,lim=5)
## by.n <- sapply(ns, function(n) {
by.n <- mclapply(ns, mc.cores=detectCores()-4, FUN=function(n) {
    cat('.')
    replicate(10, {
        tryCatch({
            f.0 <- function(x)pdf.index.glm.gaussian(x,g=0,beta=beta.hat,link=h,params=params,lim=lim)
            f.1 <- function(x)pdf.index.glm.gaussian(x,g=1,beta=beta.hat,link=h,params=params,lim=lim)
            F.0 <- Vectorize(function(x)integrate(function(w)f.0(w),-lim,x)$val)
            F.1 <- Vectorize(function(x)integrate(function(w)f.1(w),-lim,x)$val)
            ## F.0(x.1%*%beta.hat)
            x <- rmvnorm(n,mu,Sigma)
            risk <- h(x%*%beta.0)
            g <- rbinom(n,1,risk)
            x.0 <- x[g==0,]
            x.1 <- x[g==1,]
            ## beta.hat <- beta.star.full + rnorm(p.full)/sqrt(n)
            beta.hat <- coef(glm(g~x-1,family=binomial(link=params$link.name)))
            auc.hajek.hat <- auc.hajek(x=x.0%*%beta.hat,y=x.1%*%beta.hat,F=F.0,G=F.1,auc=auc.full,terms.only=FALSE)
            auc.hat.full.hat <- auc.hat(x=x.0%*%beta.hat,y=x.1%*%beta.hat)
            auc.full.hat <- auc.glm.gaussian(beta.hat,h,params,lim=5)
            obs <- auc.hat.full.hat - auc.full.hat
            try <- auc.hajek.hat
            obs - try},
            error=function(e){print(e); NA})
    })
})
by.n <- simplify2array(by.n)
Sys.time() - start
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)
## save.image('sessions/14d.RData')


## 14e influence function


## checking score formula against known probit case
score <- function(x,g,params) { # x in model matrix format
    h <- params$link; h.1 <- params$link.deriv
    eta <- as.numeric(x%*%params$beta)
    t( (g-h(eta))/h(eta)/(1-h(eta))*h.1(eta) * x )
}
fi <- function(x,g,params) { # x in model matrix format
    h <- params$link; h.1 <- params$lnk.deriv; h.2 <- params$link.deriv2
    eta <- as.numeric(x%*%params$beta)
    denom <- h(eta)*(1-h(eta))
    t((h.2(eta)*(g-h(eta))/denom - h.1(eta)/denom^2 * (g*(h.1(eta)-2*h(eta)*h.1(eta))+h.1(eta)*h(eta)^2))*x) %*% x
}
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 1e2
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- pnorm
h.1 <- dnorm
h.2 <- function(x)-x*dnorm(x)
beta.0 <- runif(p.full)/2
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2)
x <- rmvnorm(n,mu,Sigma)
eta <- as.numeric(x%*%beta.0)
risk <- h(eta)
g <- rbinom(n,1,risk)
sum(abs(score(x,g,params) - t(dnorm(eta)*(g/pnorm(eta) - (1-g)/(1-pnorm(eta)))*x)))

## checking fisher inf formula against known logit case
score <- function(x,g,params) { # x in model matrix format
    h <- params$link; h.1 <- params$lnk.deriv
    eta <- as.numeric(x%*%params$beta)
    t( (g-h(eta))/h(eta)/(1-h(eta))*h.1(eta) * x )
}
fi <- function(x,g,params) { # x in model matrix format
    h <- params$link; h.1 <- params$lnk.deriv; h.2 <- params$link.deriv2
    eta <- as.numeric(x%*%params$beta)
    denom <- h(eta)*(1-h(eta))
    -t((h.2(eta)*(g-h(eta))/denom - h.1(eta)/denom^2 * (g*(h.1(eta)-2*h(eta)*h.1(eta))+h.1(eta)*h(eta)^2))*x) %*% x
}
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 1e2
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
beta.0 <- runif(p.full)/2
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2)
x <- rmvnorm(n,mu,Sigma)
eta <- as.numeric(x%*%beta.0)
risk <- h(eta)
g <- rbinom(n,1,risk)
## sum(abs(score(x,g,params) - t(dnorm(eta)*(g/pnorm(eta) - (1-g)/(1-pnorm(eta)))*x)))
fi(x,g,params)
t(x)%*%diag(exp(eta)/(1+exp(eta))^2)%*%x


## checking infl fun has right rate of convergence
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 1e2
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
beta.0 <- runif(p.full)/2
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.full)
ns <- round(seq(1e2,1e3,len=20))
by.n <- sapply(ns, function(n) {
    replicate(1e2, {
        x <- rmvnorm(n,mu,Sigma)
        ## x <- matrix(rgamma(n*p.full,1,5),ncol=p.full)
        eta <- as.numeric(x%*%beta.0)
        risk <- h(eta)
        g <- rbinom(n,1,risk)
        beta.hat <- coef(glm(g~x-1,family=binomial(link=link.name)))
        obs <- beta.hat-beta.0
        ## infl <- solve(fi(x,g,params))%*%score(x,g,params)
        ## approx <- rowSums(infl)
        approx <- infl.glm(x,g,params,terms.only=FALSE)
        obs-approx
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)


## check influence function standardizes beta.hat
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 5e3
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
beta.0 <- runif(p.full)/2
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.full)
## ns <- round(seq(1e2,1e3,len=20))
z.stats <-     replicate(4e2, {
    x <- rmvnorm(n,mu,Sigma)
    ## x <- matrix(runif(n*p.full),nrow=n)
    eta <- as.numeric(x%*%beta.0)
    risk <- h(eta)
    g <- rbinom(n,1,risk)
    beta.hat <- coef(glm(g~x-1,family=binomial(link=link.name)))
    obs <- beta.hat-beta.0
    infl <- infl.glm.gaussian(x,g,params,terms.only=TRUE)
    var.hat <- var(t(infl))/n
    expm::sqrtm(solve(var.hat))%*%obs 
})
z.stats <- z.stats[,,]
qqnorm(z.stats)
abline(0,1)

hist(z.stats,prob=TRUE)
curve(dnorm,add=TRUE)


dd

## ## checking score formula against known probit case--reduced coefs
## score <- function(x,g,params) { # x in model matrix format
##     h <- params$link; h.1 <- params$link.deriv
##     eta <- as.numeric(x%*%params$beta)
##     t( (g-h(eta))/h(eta)/(1-h(eta))*h.1(eta) * x )
## }
## fi <- function(x,g,params) { # x in model matrix format
##     h <- params$link; h.1 <- params$lnk.deriv; h.2 <- params$link.deriv2
##     eta <- as.numeric(x%*%params$beta)
##     denom <- h(eta)*(1-h(eta))
##     t((h.2(eta)*(g-h(eta))/denom - h.1(eta)/denom^2 * (g*(h.1(eta)-2*h(eta)*h.1(eta))+h.1(eta)*h(eta)^2))*x) %*% x
## }
## require(mvtnorm)
## require(parallel)
## source('misc.R')
## set.seed(1)
## p.full <- 4
## p.red <- 2
## n <- 1e2
## mu <- runif(p.full)
## Sigma <- matrix(runif(p.full^2),p.full)
## Sigma <- Sigma%*%t(Sigma)
## h <- pnorm
## h.1 <- dnorm
## h.2 <- function(x)-x*dnorm(x)
## beta.0 <- runif(p.full)/2
## params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2)
## beta.star.red <- coef.reduced.glm.gadssian(p.red,params,5)

## x <- rmvnorm(n,mu,Sigma)
## eta <- as.numeric(x%*%beta.0)
## risk <- h(eta)
## g <- rbinom(n,1,risk)
## sum(abs(score(x,g,params) - t(dnorm(eta)*(g/pnorm(eta) - (1-g)/(1-pnorm(eta)))*x)))


## infl function conv rate with reduced coefs
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 1e2
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
beta.0 <- runif(p.full)/2
params.full <- list(mu=mu,Sigma=Sigma,beta=beta.0,p=p.full,link=h,link.deriv=h.1,link.deriv2=h.2)
beta.star.red <- coef.reduced.glm.gaussian(p.red,params.full,lim=5)
params.red <- list(mu=mu[1:p.red],Sigma=Sigma[1:p.red,1:p.red],beta=beta.star.red,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
## beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
ns <- round(seq(1e2,1e3,len=20))
by.n <- sapply(ns, function(n) {
    replicate(1e2, {
        x <- rmvnorm(n,mu,Sigma)
        ## x <- matrix(rexp(n*p.full)
        eta <- as.numeric(x%*%beta.0)
        risk <- h(eta)
        g <- rbinom(n,1,risk)
        ## beta.hat <- coef(glm(g~x-1,family=binomial(link=link.name)))
        ## obs <- beta.hat-beta.0
        ## approx <- infl.glm.gaussian(x,g,params.full,terms.only=FALSE)
        ## obs-approx
        x.red <- x[,1:p.red]
        beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link=link.name)))
        obs <- beta.hat.red - beta.star.red
        approx <- infl.glm.gaussian(x.red,g,params.red,terms.only=FALSE)
        obs - approx
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)

## check variance from infl function on redcued coef
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 5e3
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
beta.0 <- runif(p.full)/2
params.full <- list(mu=mu,Sigma=Sigma,beta=beta.0,p=p.full,link=h,link.deriv=h.1,link.deriv2=h.2)
beta.star.red <- coef.reduced.glm.gaussian(p.red,params.full,15)
params.red <- list(mu=mu[1:p.red],Sigma=Sigma[1:p.red,1:p.red],beta=beta.star.red,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
## beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
## ns <- round(seq(1e2,1e3,len=20))
## by.n <- sapply(ns, function(n) {
z.stats <-     replicate(4e2, {
    x <- rmvnorm(n,mu,Sigma)
    eta <- as.numeric(x%*%beta.0)
    risk <- h(eta)
    g <- rbinom(n,1,risk)
    ## beta.hat <- coef(glm(g~x-1,family=binomial(link=link.name)))
    ## obs <- beta.hat-beta.0
    ## approx <- infl.glm.gaussian(x,g,params.full,terms.only=FALSE)
    ## obs-approx
    x.red <- x[,1:p.red]
    beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link=link.name)))
    obs <- beta.hat.red - beta.star.red
    infl <- infl.glm.gaussian(x.red,g,params.red,terms.only=TRUE)
    var.hat <- var(t(infl))/n
    expm::sqrtm(solve(var.hat))%*%obs 
})
z.stats <- z.stats[,,]
op <- par(mfrow=c(1,2))
qqnorm(z.stats)
abline(0,1)
hist(z.stats,prob=TRUE)
curve(dnorm,add=TRUE)
par(op)

dd


require(mvtnorm)
require(numDeriv)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 5e3
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
beta.0 <- runif(p.full)/2
## bb <- runif(p.full)
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2)
beta.star.full <- coef.reduced.glm.gaussian(p.full,params)
auc.full <- auc.glm.gaussian(beta.star.full,params)
## ns <- round(seq(1e2,1e3,len=20))
## by.n <- sapply(ns, function(n) {
z.stats <- replicate(3e2, {
    x <- rmvnorm(n,mu,Sigma)
    eta <- as.numeric(x%*%beta.0)
    risk <- h(eta)
    g <- rbinom(n,1,risk)
    x.0.full <- x[g==0,]; x.1.full <- x[g==1,]
    beta.hat.full <- coef(glm(g~x-1,family=binomial(link=link.name)))
    auc.hat.full.hat <- auc.hat(x.0.full%*%beta.hat.full,x.1.full%*%beta.hat.full)
    obs <- auc.hat.full.hat - auc.full
    auc.hajek.hat <- auc.hajek(x=x.0.full%*%beta.hat.full,y=x.1.full%*%beta.hat.full,F=NULL,G=NULL,auc=auc.hat.full.hat,terms.only=TRUE,IID=TRUE)
    params.hat <- list(mu=colMeans(x),Sigma=cov(x),beta=beta.hat.full,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.full)
    infl.hat <- infl.glm.gaussian(x,g,params.hat)
    ## deriv.full <- rep(0,p.full)
    deriv.full <- grad(function(beta)auc.hat(x.0.full%*%beta,x.1.full%*%beta),beta.hat.full,method.args=list(r=6))
    approx <- auc.hajek.hat + as.numeric(deriv.full%*%infl.hat)
    obs / sqrt(var(approx) / length(approx))
    deriv.full
})
qqnorm(z.stats)
abline(0,1)



## checking z-stats on reduced data
require(mvtnorm)
require(numDeriv)
source('misc.R')
set.seed(2)
p.full <- 8
p.red <- 3
n <- 3e3
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
beta.0 <- runif(p.full)/2
params.full <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2)
beta.star.red <- as.numeric(coef.reduced.glm.gaussian(p.red,params.full,lim=5))
## beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
params.red <- list(mu=mu[1:p.red],Sigma=Sigma[1:p.red,1:p.red],beta=beta.star.red,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
auc.red <- auc.glm.gaussian(beta.star.red,params.full)
## ns <- round(seq(1e2,1e3,len=20))
## by.n <- sapply(ns, function(n) {
## set.seed(2)
z.stats <- replicate(1e2, {
    cat('.')
    x <- rmvnorm(n,mu,Sigma)
    ## x <- matrix(runif(n*p.full),ncol=p.full)#!!!
    risk <- h(x%*%beta.0)
    ## risk <- h((x%*%beta.0)^2)##!!!
    ## risk <- h((x^2)%*%beta.0)##!!!
    g <- rbinom(n,1,risk)
    x.red <- x[,1:p.red]
    x.0 <- x[g==0,]; x.1<- x[g==1,]    
    x.0.red <- x.red[g==0,]; x.1.red <- x.red[g==1,]    
    beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link=link.name)))
    ## beta.hat.red <- c(beta.hat.red,rep(0,p.full-p.red))
    auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
    obs <- auc.hat.red.hat - auc.red
    auc.hajek.red <- auc.hajek(x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red,F=NULL,G=NULL,auc=auc.hat.red.hat,terms.only=TRUE,IID=TRUE)
    params.hat.red <- list(mu=colMeans(x.red),Sigma=cov(x.red),beta=beta.hat.red,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
    infl.hat <- infl.glm.gaussian(x.red,g,params.hat.red)
    deriv.red <- grad(function(beta)auc.hat(x.0.red%*%beta,x.1.red%*%beta),beta.hat.red[1:p.red])
    approx <- auc.hajek.red + as.numeric(deriv.red%*%infl.hat) 
    obs / sqrt(var(approx) / length(approx))
    ## deriv.red[3]
})
## mean(z.stats)
qqnorm(z.stats)
abline(0,1)


n<- 1e6
x <- rmvnorm(n,mu,Sigma)
risk <- h((x%*%beta.0)^1)##!!!
g <- rbinom(n,1,risk)
x.red <- x[,1:p.red]
beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link=link.name)))
beta.hat.red - beta.star.red
rbind(beta.hat.red, beta.star.red)



seed <- seed+1
require(mvtnorm)
require(numDeriv)
source('misc.R')
set.seed(seed)
start <- Sys.time()
p.full <- 4
p.red <- 2
n <- 3e3
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
beta.0 <- runif(p.full)
params.full <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2)
beta.star.red <- as.numeric(coef.reduced.glm.gaussian(p.red,params.full,lim=10))
beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
## params.red <- list(mu=mu[1:p.red],Sigma=Sigma[1:p.red,1:p.red],beta=beta.star.red,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
auc.red <- auc.glm.gaussian(beta.star.red,params.full)
grad(function(beta)auc.glm.gaussian(c(beta[1:2],0,0),params.full,5),beta.star.red)
Sys.time() - start

a <- runif(2)
ts <- seq(-1,1,len=10)/2
vals <- sapply(ts, function(t)auc.glm.gaussian(c(beta.star.red[1:2]+t*a,0,0),params.full,5))
plot(ts,vals)
abline(v=0)


seed <- 1

seed <- seed+1
set.seed(seed)
p <- 4
n <- 1e4
gamma <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
## x <- matrix(rexp(n*p),ncol=p)%*%expm::sqrtm(Sigma)
x <- matrix(rbeta(n*p,runif(1),runif(1)),ncol=p)%*%expm::sqrtm(Sigma)
x <- matrix(rbeta(n*p,2,2),ncol=p)%*%expm::sqrtm(Sigma)
risk <- plogis(x%*%gamma)
y <- rbinom(n,1,risk)
beta.lda <- solve(Sigma)%*%(colMeans(x[y==1,])-colMeans(x[y==0,]))
plot(x%*%beta.lda,risk)


p <- 4
n <- 1e4
gamma <- matrix(runif(p),ncol=1)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
x <- matrix(runif(p),ncol=1)
sum(diag((t(gamma)%*%x) * (t(gamma)%*%Sigma%*%x)))
sum(diag((t(gamma)%*%x) * (t(x)%*%Sigma%*%gamma)))
sum(diag(gamma%*%t(gamma)%*%x%*%t(x)%*%Sigma))


p <- 4
n <- 1e4
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
pairs <- replicate(1e2,{
x <- matrix(rnorm(p),ncol=1)
gamma <- matrix(rnorm(p),ncol=1)
c(t(gamma)%*%x,t(gamma)%*%Sigma%*%x)
if(prod(c(t(gamma)%*%x,t(gamma)%*%Sigma%*%x))<0)browser()
cat('.')
})


sum(diag(t(gamma)%*%x * t(x)%*%Sigma%*%gamma))
sum(diag((gamma%*%t(gamma))%*%(x %*% t(x))%*%Sigma))

A <- (gamma%*%t(gamma))%*%(x %*% t(x))
eigen(A)$val



p <- 4
n <- 1e2
params <- replicate(2,{
    Sigma <- matrix(runif(p^2),p)
    Sigma <- Sigma%*%t(Sigma)
    mu <- runif(p)
    list(mu=mu,Sigma=Sigma)
},simplify=FALSE)
pi.0 <- runif(1); pi.1 <- 1-pi.0
f0 <- params[[1]]; f1 <- params[[2]]
x <- t(matrix(rnorm(n*p),n))
risk <- lapply(list(f0,f1), function(f) with(f, diag(t(x-mu)%*%solve(Sigma)%*%(x-mu))  + 2*t(mu)%*%solve(Sigma)%*%x))
risk <- -(risk[[1]]+risk[[2]])
lda <- t(f1$mu-f0$mu)%*%solve(pi.0*f0$Sigma+pi.1*f1$Sigma)%*%x


## adversarial parameters
require(mvtnorm)
require(numDeriv)
source('misc.R')
seed <- 1

## component at .2 with seed==21. update: no, set lim=10 and it becomes much smaller
seed <- seed+1
set.seed(seed)
start <- Sys.time()
p.full <- 4
p.red <- 2
n <- 3e3
## mu <- runif(p.full)
mu <- rep(0,p.full)
## Sigma <- matrix(runif(p.full^2),p.full)
## Sigma <- Sigma%*%t(Sigma)
Sigma <- diag(p.full)*1.5+2
## Sigma <- diag(c(rep(0.1,p.full/2),rep(1,p.full/2)))*4
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
beta.0 <- runif(p.full)
params.full <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2)
beta.star.red <- as.numeric(coef.reduced.glm.gaussian(p.red,params.full,lim=10))
beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
## params.red <- list(mu=mu[1:p.red],Sigma=Sigma[1:p.red,1:p.red],beta=beta.star.red,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
auc.red <- auc.glm.gaussian(beta.star.red,params.full)
grad0 <- grad(function(beta)auc.glm.gaussian(c(beta[1:2],0,0),params.full,5),beta.star.red)
print(grad0)
print(Sys.time() - start)
## a <- runif(2)
## a <- c(1,0)
## ts <- seq(-1,1,len=10)/10
## vals <- sapply(ts, function(t)auc.glm.gaussian(c(beta.star.red[1:2]+t*a,0,0),params.full,5))
plot(ts,vals,asp=1)
abline(v=0)

dd




## 15 lda + gaussian covariates but different Sigma's (misspecified lda)

## 15a checking formulas

## check auc formula
require(mvtnorm)
source('misc.R')
p <- 4
n <- 1e2
mu.1 <- runif(p)
mu.0 <- runif(p)
Sigma.0 <- matrix(runif(p^2),p)
Sigma.0 <- Sigma.0%*%t(Sigma.0)
Sigma.1 <- matrix(runif(p^2),p)
Sigma.1 <- Sigma.1%*%t(Sigma.1)
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1)
beta <- runif(p)
pi.0 <- 1/3; pi.1 <- 1-pi.0
n.0 <- round(n*pi.0); n.1 <- round(n*pi.1)
aucs <- replicate(1e3, {
    x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
    x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
    auc.hat(x.0%*%beta,x.1%*%beta)
})
hist(aucs)
try <- auc.lda.gaussian(beta,params)
abline(v=try,col=2)
abline(v=mean(aucs),col=3)

## check beta_{LDA} formula
require(mvtnorm)
source('misc.R')
set.seed(1)
p <- 4
n <- 1e3
mu.1 <- runif(p)
mu.0 <- runif(p)
Sigma.0 <- matrix(runif(p^2),p)
Sigma.0 <- Sigma.0%*%t(Sigma.0)
Sigma.1 <- matrix(runif(p^2),p)
Sigma.1 <- Sigma.1%*%t(Sigma.1)
## Sigma.1 <- Sigma.0
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1)
beta <- runif(p)
pi.0 <- 1/3; pi.1 <- 1-pi.0
n.0 <- round(n*pi.0); n.1 <- round(n*pi.1)
a <- runif(p)
obs <- replicate(1e3, {
    x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
    x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
    coefs.lda(x.0,x.1,params=list())
})
obs <- obs[,,]
hist(a%*%obs)
try <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
abline(v=a%*%try,col=2)
abline(v=mean(a%*%obs),col=3)


## check formula for derivative of auc
source('misc.R')
p <- 4
mu.1 <- runif(p)
mu.0 <- runif(p)
Sigma.0 <- matrix(runif(p^2),p)
Sigma.0 <- Sigma.0%*%t(Sigma.0)
Sigma.1 <- matrix(runif(p^2),p)
Sigma.1 <- Sigma.1%*%t(Sigma.1)
## Sigma.1 <- Sigma.0
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1)
viz.deriv(f=function(beta)auc.lda.gaussian(beta,params), deriv=function(beta)auc.deriv.lda.gaussian(beta,params),p,1)

## checking derivative of auc at beta.star is nonzero
source('misc.R')
p <- 4
mu.0 <- rep(0,p)
mu.1 <- runif(p)/10
Sigma.0 <- matrix(runif(p^2),p)
Sigma.0 <- Sigma.0%*%t(Sigma.0)
Sigma.1 <- matrix(runif(p^2),p)
Sigma.1 <- Sigma.1%*%t(Sigma.1)
pi.0 <- 1/3#runif(1)
pi.1 <- 1-pi.0
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
beta <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
auc.deriv.lda.gaussian(beta,params) 


## how big does the deriv get?
derivs <- replicate(1e5, {
    mu.0 <- rep(0,p)
    mu.1 <- runif(p) * 10
    Sigma.0 <- matrix(runif(p^2),p)
    Sigma.0 <- Sigma.0%*%t(Sigma.0) * 10
    Sigma.1 <- matrix(runif(p^2),p) 
    Sigma.1 <- Sigma.1%*%t(Sigma.1) * 10
    pi.0 <- runif(1)
    pi.1 <- 1-pi.0
    params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
    beta <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
    deriv <- auc.deriv.lda.gaussian(beta,params)
    max(abs(deriv))
})
max(derivs) ## ~ .8


## 15b taylor expansion is O(1/n). 

require(mvtnorm)
source('misc.R')
p <- 4
n <- 1e3
mu.1 <- runif(p)
mu.0 <- runif(p)
Sigma.0 <- matrix(runif(p^2),p)
Sigma.0 <- Sigma.0%*%t(Sigma.0)
Sigma.1 <- matrix(runif(p^2),p)
Sigma.1 <- Sigma.1%*%t(Sigma.1)
## Sigma.1 <- Sigma.0
pi.0 <- 1/3; pi.1 <- 1-pi.0
## load('sessions/15d-2.RData')
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
## beta <- runif(p)
## a <- runif(p)
beta.star <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
ns <- round(seq(1e3,7e3,len=20))
by.n <- sapply(ns, function(n) {
    n.0 <- round(n*pi.0); n.1 <- round(n*pi.1)
    obs <- replicate(1e2, {
        x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
        x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
        beta.hat <- coefs.lda(x.0,x.1,params=list(Sigma=NULL))
        obs <- auc.lda.gaussian(beta.hat,params) - auc.lda.gaussian(beta.star,params)
        approx <- deriv.star%*%(beta.hat - beta.star)
        obs - approx
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)

## taylor remainder seems O(1/sqrt(n)) with adversarial parameters
require(mvtnorm)
source('misc.R')
p <- 4
n <- 1e3
## load('sessions/15d-2.RData')
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
beta.star <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
ns <- round(seq(1e3,7e3,len=20))
by.n <- sapply(ns, function(n) {
    n.0 <- round(n*pi.0); n.1 <- round(n*pi.1)
    obs <- replicate(1e2, {
        x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
        x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
        beta.hat <- coefs.lda(x.0,x.1,params=list(Sigma=NULL))
        obs <- auc.lda.gaussian(beta.hat,params) - auc.lda.gaussian(beta.star,params)
        approx <- deriv.star%*%(beta.hat - beta.star)
        obs - approx
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)

load('sessions/15d-2.RData')
source('misc.R')
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
beta.star <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
a <- runif(length(beta.star))
ts <- seq(-1,1,len=1e4)/1e3
aucs <- sapply(ts, function(t)auc.lda.gaussian(beta.star+t*a,params))
plot(ts,aucs,type='l')

dd


## 15c influence function. 0 as before when Sigma is provided, 1/sqrt(n) rate otherwise.
require(mvtnorm)
source('misc.R')
p <- 4
n <- 1e3
mu.1 <- runif(p)
mu.0 <- runif(p)
Sigma.0 <- matrix(runif(p^2),p)
Sigma.0 <- Sigma.0%*%t(Sigma.0)
Sigma.1 <- matrix(runif(p^2),p)
Sigma.1 <- Sigma.1%*%t(Sigma.1)
## Sigma.1 <- Sigma.0
pi.0 <- 1/3; pi.1 <- 1-pi.0
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
## beta <- runif(p)
## a <- runif(p)
beta.star <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
ns <- round(seq(1e2,1e3,len=20))
by.n <- sapply(ns, function(n) {
    n.0 <- round(n*pi.0); n.1 <- round(n*pi.1)
    obs <- replicate(1e2, {
        x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
        x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
        beta.hat <- coefs.lda(x.0,x.1,params=list())#Sigma=pi.0*Sigma.0+pi.1*Sigma.1))
        x <- rbind(x.0,x.1)
        g <- c(rep(0,n.0),rep(1,n.1))
        approx <- infl.lda(x,g,params,terms.only=FALSE)
        obs <- beta.hat - beta.star
        sum(abs(obs - approx))
        ## sum(abs(beta.hat - beta.star))
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)



## 15c hajek
require(mvtnorm)
source('misc.R')
p <- 4
n <- 1e3
mu.1 <- runif(p)
mu.0 <- runif(p)
Sigma.0 <- matrix(runif(p^2),p)
Sigma.0 <- Sigma.0%*%t(Sigma.0)
Sigma.1 <- matrix(runif(p^2),p)
Sigma.1 <- Sigma.1%*%t(Sigma.1)
## Sigma.1 <- Sigma.0
pi.0 <- 1/3; pi.1 <- 1-pi.0
## load('sessions/15d-2.RData')
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
## beta <- runif(p)
## a <- runif(p)
beta.star <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
ns <- round(seq(1e2,1e3,len=20))
by.n <- sapply(ns, function(n) {
    n.0 <- round(n*pi.0); n.1 <- round(n*pi.1)
    obs <- replicate(1e2, {
        x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
        x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
        beta.hat <- coefs.lda(x.0,x.1,params=list())
        x <- rbind(x.0,x.1)
        g <- c(rep(0,n.0),rep(1,n.1))
        obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc.lda.gaussian(beta.hat,params)
        ## F <- function(x)pnorm((x-t(beta.hat)%*%mu.0)/sqrt(t(beta.hat)%*%Sigma.0%*%beta.hat))
        ## G <- function(y)pnorm((y-t(beta.hat)%*%mu.1)/sqrt(t(beta.hat)%*%Sigma.1%*%beta.hat))
        ## approx <- auc.hajek((x.0%*%beta.hat)[,],(x.1%*%beta.hat)[,],F=F,G=G,auc=auc.lda.gaussian(beta.hat,params),terms.only=FALSE)
        approx <- auc.hajek.lda.gaussian(x,g,beta.hat,params,terms.only=FALSE)
        ## approx <- auc.hajek.lda.gaussian(x,g,beta.star,params,terms.only=FALSE)
        obs - approx
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)

sds <- apply(by.n,2,sd)
plot(ns,sds)
lm0 <- lm(log(sds)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)


## 15d combining
require(mvtnorm)
source('misc.R')
set.seed(1)
p <- 4
n <- 1e3
## load('sessions/15d-2.RData')
mu.1 <- runif(p)
mu.0 <- runif(p)
Sigma.0 <- matrix(runif(p^2),p)
Sigma.0 <- Sigma.0%*%t(Sigma.0)
Sigma.1 <- matrix(runif(p^2),p)
Sigma.1 <- Sigma.1%*%t(Sigma.1)
## Sigma.1 <- Sigma.0
pi.0 <- 1/3; pi.1 <- 1-pi.0
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
## beta <- runif(p)
## a <- runif(p)
beta.star <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
ns <- round(seq(1e2,1e3,len=20))
by.n <- sapply(ns, function(n) {
    n.0 <- round(n*pi.0); n.1 <- round(n*pi.1)
    obs <- replicate(1e2, {
        x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
        x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
        beta.hat <- coefs.lda(x.0,x.1,params=params)
        x <- rbind(x.0,x.1)
        g <- c(rep(0,n.0),rep(1,n.1))
        ## obs.1 <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc.lda.gaussian(beta.hat,params)
        ## approx.1 <- auc.hajek.lda.gaussian(x,g,beta.hat,params,terms.only=FALSE)
        approx.1 <- auc.hajek.lda.gaussian(x,g,beta.hat,params,terms.only=TRUE,IID=TRUE)
        ## obs - approx
        ## infl <- infl.lda(x,g,params,terms.only=FALSE)
        infl <- infl.lda(x,g,params,terms.only=TRUE)
        ## obs.2 <- auc.lda.gaussian(beta.hat,params) - auc.lda.gaussian(beta.star,params)
        ## approx.2 <- deriv.star%*%(beta.hat - beta.star)
        approx.2 <- deriv.star%*%infl
        ## obs.1+obs.2 - (approx.1+approx.2)
        obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc.lda.gaussian(beta.star,params)
        approx <- approx.1+approx.2
        obs - mean(approx)
        ## obs.2 - approx.2
        ## approx <- infl.lda(x,g,params,terms.only=FALSE)
        ## sum(abs(obs - approx))
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)


## adversarial data parameters
source('misc.R')
p <- 4
while(1) {
    mu.0 <- rep(0,p)
    mu.1 <- runif(p) * 10
    Sigma.0 <- matrix(runif(p^2),p)
    Sigma.0 <- Sigma.0%*%t(Sigma.0) * 10
    Sigma.1 <- matrix(runif(p^2),p) 
    Sigma.1 <- Sigma.1%*%t(Sigma.1) * 10
    pi.0 <- runif(1)
    pi.1 <- 1-pi.0
    params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
    beta <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
    deriv <- auc.deriv.lda.gaussian(beta,params)
    if(max(abs(deriv))>1)break
}


save(mu.0,mu.1,Sigma.0,Sigma.1,pi.0,pi.1,file='sessions/15d-1.RData')


## adversarial parameters. derivative can be > 20, scale of Sigma controls.
source('misc.R')
p <- 2
m <- 1e3
max.deriv <- -Inf
while(1) {
    ## Sigma.0 <- matrix(runif(p^2),p)
    ## Sigma.0 <- Sigma.0%*%t(Sigma.0)  *m
    ## Sigma.1 <- matrix(runif(p^2),p) 
    ## Sigma.1 <- Sigma.1%*%t(Sigma.1)  *m
Sigma.0 <- diag(runif(p)) * m
Sigma.0 <- diag(p) * m
Sigma.1 <- diag(runif(p)) * m
    pi.0 <- runif(1)
    pi.1 <- 1-pi.0
    beta <- runif(p) /m
    beta <- rep(1,p )/m
    ## mu.0 <- rep(0,p)
    ## mu.1 <- (pi.0*Sigma.0+pi.1*Sigma.1)%*%beta
    ## params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
    ## auc.deriv.lda.gaussian(beta,params)
    Sigma.pi <- pi.0*Sigma.0+pi.1*Sigma.1
    Sigma <- Sigma.0+Sigma.1
    quad <- as.numeric(t(beta)%*%Sigma%*%beta)
    quad.pi <- as.numeric(t(beta)%*%Sigma.pi%*%beta)
    ## deriv <- max(abs(dnorm(quad.pi/sqrt(quad))/sqrt(quad)*t(beta)%*%(Sigma.pi-quad.pi/quad*Sigma)))
    deriv <- max(abs(1/sqrt(quad)*t(beta)%*%(Sigma.pi-quad.pi/quad*Sigma)))
    if(deriv>max.deriv){
        max.deriv <<- deriv
        print(max.deriv)
    }
    if(deriv>4)break
}
mu.0 <- rep(0,p)
mu.1 <- (pi.0*Sigma.0+pi.1*Sigma.1)%*%beta
## params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
save(mu.0,mu.1,Sigma.0,Sigma.1,pi.0,pi.1,file='sessions/15d-2.RData')


p <- 4
pi.0 <- 3/4; pi.1 <- 1-pi.0
Sigma.0 <- matrix(runif(p^2),p)
Sigma.0 <- Sigma.0%*%t(Sigma.0)  /10
Sigma.1 <- matrix(runif(p^2),p) 
Sigma.1 <- Sigma.1%*%t(Sigma.1)  /10
Sigma.pi <- pi.0*Sigma.0+pi.1*Sigma.1
Sigma <- Sigma.0+Sigma.1
mu <- eigen(Sigma%*%solve(Sigma.pi))$vec[,1]


A <- matrix(runif(p^2),p)
A <- A%*%t(A)  
B <- matrix(runif(p^2),p) 
B <- B%*%t(B)  
mu <- eigen(A%*%solve(B))$vec[,1]
(A%*%solve(B)%*%mu) / mu
(t(mu)%*%solve(B)%*%A%*%solve(B)%*%mu) / (t(mu)%*%solve(B)%*%mu)

require(expm)
Sigma%*%solve(Sigma.pi)
(Sigma.0+Sigma.1)%*%solve(diag(p)+pi.1/pi.0*solve(Sigma.0)%*%Sigma.1)%*%solve(Sigma.0)/pi.0

solve(diag(p)+pi.1/pi.0*solve(Sigma.0)%*%Sigma.1)


dd

## require(expm)
## p <- 4
## m <- 1e4
## A <- matrix(runif(p^2),p)
## A <- A%*%t(A) 
## B <- matrix(runif(p^2),p) 
## B <- B%*%t(B) 
## ev <- eigen(solve(B)%*%A)
## colSums((sqrtm(A)%*%ev$vectors)^2)/colSums((sqrtm(B)%*%ev$vectors)^2) - ev$values


## 15d checking z-stats, using true derivative at beta.star
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
load('sessions/15d-2.RData') ## very adversarial (auc almost discontinuous at beta.star)
load('sessions/15d-1.RData')
pi.0 <- .95; pi.1 <- 1-pi.0
p <- 4
n <- 2e4
## mu.1 <- runif(p)
## mu.0 <- runif(p)
## Sigma.0 <- matrix(runif(p^2),p)
## Sigma.0 <- Sigma.0%*%t(Sigma.0)
## Sigma.1 <- matrix(runif(p^2),p)
## Sigma.1 <- Sigma.1%*%t(Sigma.1)
## Sigma.1 <- Sigma.0
## pi.0 <- 1/3; pi.1 <- 1-pi.0
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
beta.star <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
## ns <- round(seq(1e2,1e3,len=20))
## by.n <- sapply(ns, function(n) {
n.0 <- round(n*pi.0); n.1 <- n-n.0
## z.stats <- replicate(5e2, {
z.stats <- mclapply(1:5e2, mc.cores=detectCores()-3, FUN= function(dd){
    x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
    x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
    beta.hat <- coefs.lda(x.0,x.1,params=list())
    params.hat <- list(mu.0=colMeans(x.0),mu.1=colMeans(x.1),pi.0=nrow(x.0)/(nrow(x.0)+nrow(x.1)),Sigma.0=var(x.0),Sigma.1=var(x.1))
    x <- rbind(x.0,x.1)
    g <- c(rep(0,n.0),rep(1,n.1))
    approx.1 <- auc.hajek.lda.gaussian(x,g,beta.hat,params.hat,terms.only=TRUE,IID=TRUE)
    infl <- infl.lda(x,g,params.hat,terms.only=TRUE)
    approx.2 <- deriv.star%*%infl
    obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc.lda.gaussian(beta.star,params)
    approx <- as.numeric(approx.1+approx.2)
    ## obs - mean(approx)
    ## obs.2 - approx.2
    ## approx <- infl.lda(x,g,params,terms.only=FALSE)
    ## obs <- diff.hat.coef.hat - diff.coef
    ## approx <- diff.hajek.coef.hat + as.numeric(deriv %*% infl)
    obs / sqrt(var(approx) / length(approx))
})
z.stats <- unlist(z.stats)
## })
## mad <- colMeans(abs(by.n))
## plot(ns,mad)
## lm0 <- lm(log(mad)~log(ns))
## coef(lm0)
## curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)
op <- par(mfrow=c(1,2))
qqnorm(z.stats)
abline(0,1)
hist(z.stats,prob=TRUE)
curve(dnorm,add=TRUE)
par(op)



## checking z-stats, estimating derivative at beta.star
require(mvtnorm)
require(parallel)
require(numDeriv)
source('misc.R')
start <- Sys.time()
set.seed(1)
## load('sessions/15d-2.RData') ## very adversarial (auc almost discontinuous at beta.star)
load('sessions/15d-1.RData')
pi.0 <- .95; pi.1 <- 1-pi.0 ## to avoid rounding issues with the adversarial pi.0
p <- 4
n <- 5e4
## mu.1 <- runif(p)
## mu.0 <- runif(p)
## Sigma.0 <- matrix(runif(p^2),p)
## Sigma.0 <- Sigma.0%*%t(Sigma.0)
## Sigma.1 <- matrix(runif(p^2),p)
## Sigma.1 <- Sigma.1%*%t(Sigma.1)
## Sigma.1 <- Sigma.0
## pi.0 <- 1/3; pi.1 <- 1-pi.0
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
beta.star <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
## deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
## ns <- round(seq(1e2,1e3,len=20))
## by.n <- sapply(ns, function(n) {
n.0 <- round(n*pi.0); n.1 <- n-n.0
## z.stats <- replicate(5e2, {
z.stats <- mclapply(1:3e2, mc.cores=detectCores()-3, FUN= function(dd){
    x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
    x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
    beta.hat <- coefs.lda(x.0,x.1,params=list())
    params.hat <- list(mu.0=colMeans(x.0),mu.1=colMeans(x.1),pi.0=nrow(x.0)/(nrow(x.0)+nrow(x.1)),Sigma.0=var(x.0),Sigma.1=var(x.1))
    x <- rbind(x.0,x.1)
    g <- c(rep(0,n.0),rep(1,n.1))
    approx.1 <- auc.hajek.lda.gaussian(x,g,beta.hat,params.hat,terms.only=TRUE,IID=TRUE)
    infl <- infl.lda(x,g,params.hat,terms.only=TRUE)
    ## deriv.star.hat <- rep(0,p)
    deriv.star.hat <- grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat)
    approx.2 <- deriv.star.hat%*%infl
    obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc.lda.gaussian(beta.star,params)
    approx <- as.numeric(approx.1+approx.2)
    ## obs - mean(approx)
    ## obs.2 - approx.2
    ## approx <- infl.lda(x,g,params,terms.only=FALSE)
    ## obs <- diff.hat.coef.hat - diff.coef
    ## approx <- diff.hajek.coef.hat + as.numeric(deriv %*% infl)
    obs / sqrt(var(approx) / length(approx))
})
z.stats <- unlist(z.stats)
Sys.time() - start
## })
## mad <- colMeans(abs(by.n))
## plot(ns,mad)
## lm0 <- lm(log(mad)~log(ns))
## coef(lm0)
## curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)
op <- par(mfrow=c(1,2))
qqnorm(z.stats)
abline(0,1)
hist(z.stats,prob=TRUE)
curve(dnorm,add=TRUE)
par(op)
alpha <- .05
mean(abs(z.stats)<qnorm(1-alpha/2))
## save.image('sessions/15d.RData')
## n=30k, using 15-1.RData. fpr 83% when nominal rate is 95%. without derivative term though it is just 54%.
## save.image('sessions/15d-3.RData') # n=50k, 3hrs to run, fpr 88%. without deriv term 52%.

dd










## checking derivative estimate

require(mvtnorm)
require(parallel)
require(numDeriv)
source('misc.R')
set.seed(1)
## load('sessions/15d-2.RData') ## very adversarial (auc almost discontinuous at beta.star)
load('sessions/15d-1.RData')
pi.0 <- .95; pi.1 <- 1-pi.0 ## to avoid rounding issues with the adversarial pi.0
p <- 4
n <- 2e3
## mu.1 <- runif(p)
## mu.0 <- runif(p)
## Sigma.0 <- matrix(runif(p^2),p)
## Sigma.0 <- Sigma.0%*%t(Sigma.0)
## Sigma.1 <- matrix(runif(p^2),p)
## Sigma.1 <- Sigma.1%*%t(Sigma.1)
## Sigma.1 <- Sigma.0
## pi.0 <- 1/3; pi.1 <- 1-pi.0
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
beta.star <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
## ns <- round(seq(1e2,1e3,len=20))
## by.n <- sapply(ns, function(n) {
n.0 <- round(n*pi.0); n.1 <- n-n.0
## z.stats <- replicate(5e2, {
derivs <- mclapply(1:1e3, mc.cores=detectCores()-3, FUN= function(dd){
    x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
    x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
    beta.hat <- coefs.lda(x.0,x.1,params=list())
    deriv.star.hat <- grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat)
})
derivs <- simplify2array(derivs)
a <- runif(p)
hist(t(a)%*%derivs)
abline(v=deriv.star%*%a,col=2)
abline(v=mean(a%*%derivs),col=3)


## load('sessions/15d-1.RData')
## source('misc.R')
## params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
## beta.star <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
a <- runif(length(beta.star))
ts <- seq(-1,1,len=1e4+1)
stopifnot(ts[which.min(abs(ts))]==0)
aucs <- sapply(ts, function(t)auc.lda.gaussian(beta.star+t*a,params))
plot(ts,aucs,type='l')
abline(v=0)
abline(a=aucs[which.min(abs(ts))],b=deriv.star%*%a)


p <- 4
Sigma.0 <- matrix(runif(p^2),p)
Sigma.0 <- Sigma.0%*%t(Sigma.0)  *m
Sigma.1 <- matrix(runif(p^2),p) 
Sigma.1 <- Sigma.1%*%t(Sigma.1)  /m
Sigma.0 <- diag(runif(p)) * m
Sigma.1 <- diag(runif(p))
pi.0 <- runif(1)
pi.1 <- 1-pi.0
beta <- runif(p) /m
## mu.0 <- rep(0,p)
## mu.1 <- (pi.0*Sigma.0+pi.1*Sigma.1)%*%beta
## params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
## auc.deriv.lda.gaussian(beta,params)
Sigma.pi <- pi.0*Sigma.0+pi.1*Sigma.1
Sigma <- Sigma.0+Sigma.1
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
quad.pi <- as.numeric(t(beta)%*%Sigma.pi%*%beta)
deriv <- max(abs(dnorm(quad.pi/sqrt(quad))/sqrt(quad)*t(beta)%*%(Sigma.pi-quad.pi/quad*Sigma)))
deriv


## 15e simulation using parametrization giving adversarial derivative


## assumes p even, pi.0>pi.1. update: this parameterization doesn't work well; it is O(1/beta) so after multiplying by the influence function is bounded
p <- 8
target.deriv <- 1e1
pi.0 <- runif(1); pi.0 <- max(pi.0,1-pi.0)
pi.1 <- 1-pi.0
a <- p*(4*sqrt(pi*exp(1))/(pi.0-pi.1)*target.deriv)^2
b <- 0#((pi.0-pi.1)/3/sqrt(p))^(2/3)-2
c <- b+2
d <- sqrt(1+(b+c)/2)/sqrt(a*p)/(pi.0+pi.1*(b+c)/2)
Sigma.0 <- diag(p) * a
Sigma.1 <- diag(c(rep(b,p/2),rep(c,p/2))) * a
beta <- rep(1,p) * d
Sigma.pi <- pi.0*Sigma.0+pi.1*Sigma.1
Sigma <- Sigma.0+Sigma.1
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
quad.pi <- as.numeric(t(beta)%*%Sigma.pi%*%beta)
## deriv <- max(abs(dnorm(quad.pi/sqrt(quad))/sqrt(quad)*t(beta)%*%(Sigma.pi-quad.pi/quad*Sigma)))
## 1/sqrt(quad)*t(beta)%*%(Sigma.pi-quad.pi/quad*Sigma)
## t(beta)%*%(Sigma.pi-quad.pi/quad*Sigma)
## sqrt(a)*3
dnorm(quad.pi/sqrt(quad))/sqrt(quad)*t(beta)%*%(Sigma.pi-quad.pi/quad*Sigma)
1/4/sqrt(pi*exp(1))*sqrt(a/p)*(pi.0-pi.1)
## dnorm(quad.pi/sqrt(quad))
## dnorm(d*sqrt(a*p)*(pi.0+pi.1*(b+c)/2)/sqrt(1+(b+c)/2))
## 1/sqrt(2*pi*exp(1))

dd


## another adversarial parameterization, also considering infl function = O(|Sigma.pi|^.5)
require(expm)
source('misc.R')
set.seed(1)
p <- 6
m <- 1
max.deriv <- -Inf
b <- runif(1)
c <- runif(1)
while(1) {
    ## Sigma.0 <- matrix(runif(p^2),p)
    ## Sigma.0 <- Sigma.0%*%t(Sigma.0)  *m
    ## Sigma.1 <- matrix(runif(p^2),p) 
    ## Sigma.1 <- Sigma.1%*%t(Sigma.1)  *m
    ## Sigma.0 <- diag(runif(p)) * m
    Sigma.1 <- diag(p) /m
    ## Sigma.1 <- diag(runif(p)) * m
    ## Sigma.0 <- diag(c(rep(0,p/2),rep(1,p/2)) )*m
    Sigma.0 <- diag(c(rep(runif(1),p/2),rep(runif(1),p/2)) )*3
    pi.0 <- runif(1); pi.0 <- max(pi.0,1-pi.0)
    pi.1 <- 1-pi.0
    ## beta <- runif(p) 
    beta <- rep(1,p )
    ## mu.0 <- rep(0,p)
    ## mu.1 <- (pi.0*Sigma.0+pi.1*Sigma.1)%*%beta
    ## params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
    ## auc.deriv.lda.gaussian(beta,params)
    Sigma.pi <- pi.0*Sigma.0+pi.1*Sigma.1
    Sigma <- Sigma.0+Sigma.1
    quad <- as.numeric(t(beta)%*%Sigma%*%beta)
    quad.pi <- as.numeric(t(beta)%*%Sigma.pi%*%beta)
    deriv <- 1/sqrt(quad)*t(beta)%*%(Sigma.pi-quad.pi/quad*Sigma)%*%solve(sqrtm(Sigma.pi))
    ## deriv <- max(abs(1/sqrt(quad)*t(beta)%*%(Sigma.pi-quad.pi/quad*Sigma)))
    if(max(abs(deriv))>max.deriv){
        max.deriv <<- max(abs(deriv))
        print(max.deriv)
    }
    if(max.deriv>2)break
}


require(expm)
p <- 6
a <- runif(1)
b <- runif(1)
c <- runif(1)
d <- runif(1)
eps <- runif(1)
## b <- eps
## c <- 2-eps
Sigma.1 <- diag(p) 
Sigma.0 <- diag(c(rep(b,p/2),rep(c,p/2)) )
Sigma.0 <- diag(c(rep(eps,p/2),rep(2-eps,p/2)) )
pi.0 <- runif(1); pi.0 <- max(pi.0,1-pi.0)
pi.1 <- 1-pi.0
## beta <- runif(p) 
beta <- rep(1,p )*d
## mu.0 <- rep(0,p)
## mu.1 <- (pi.0*Sigma.0+pi.1*Sigma.1)%*%beta
## params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
## auc.deriv.lda.gaussian(beta,params)
Sigma.pi <- pi.0*Sigma.0+pi.1*Sigma.1
Sigma <- Sigma.0+Sigma.1
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
quad.pi <- as.numeric(t(beta)%*%Sigma.pi%*%beta)
## 1/sqrt(quad)*t(beta)%*%(Sigma.pi-quad.pi/quad*Sigma)%*%solve(sqrtm(Sigma.pi))
1/sqrt(quad)*t(beta)%*%(Sigma.pi-quad.pi/quad*Sigma)%*%solve(sqrtm(Sigma.pi))
## 1/sqrt(p)*(1+(b+c)/2)^(-3/2)*(pi.1-pi.0)*(c-b)/2*1/sqrt(pi.1+pi.0*b)
1/sqrt(p)*2^(-3/2)*(pi.1-pi.0)*(1-eps)*(1+(eps-1)*pi.0)^(-.5)

p <- 6
pi.0 <- runif(1); pi.0 <- max(pi.0,1-pi.0)
pi.1 <- 1-pi.0
hist(4*pi.0^3/27/(pi.0-pi.1)^2 - replicate(1e6,{
a <- runif(1)
b <- runif(1)
c <- runif(1)
(c-b)^2/(1+b+c)^3/(pi.0+pi.1*b)
}))

curve(4*x^3/27/(2*x-1),.75,1)

curve(-1+3*x-2*x^2,.5,.75)

require(expm)
p <- 4
g <- .1
uu <- 4*p*pi*exp(1)*g^2
pi.0.min <- 1/2*(1-uu+sqrt(uu*(uu+2)))
pi.0 <- runif(1,pi.0.min)
d <- sqrt(2/p)
a <- (1-2*pi.0)^2/16/p/pi/exp(1)
b <- -2*a-pi.0*g^2
c <- a+(pi.0-1)*g^2
f <- function(x)a*x^2+b*x+c
## eps <- uniroot(f,c(0,.1))$root
eps <- min(1/2/a*(-b+c(-1,1)*sqrt(b^2-4*a*c)))
stopifnot(eps>0)
Sigma.1 <- diag(p) 
Sigma.0 <- diag(c(rep(eps,p/2),rep(2-eps,p/2)) )
pi.1 <- 1-pi.0
## beta <- runif(p) 
beta <- rep(1,p )*d
## mu.0 <- rep(0,p)
## mu.1 <- (pi.0*Sigma.0+pi.1*Sigma.1)%*%beta
## params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
## auc.deriv.lda.gaussian(beta,params)
Sigma.pi <- pi.0*Sigma.0+pi.1*Sigma.1
Sigma <- Sigma.0+Sigma.1
quad <- as.numeric(t(beta)%*%Sigma%*%beta)
quad.pi <- as.numeric(t(beta)%*%Sigma.pi%*%beta)
## 1/sqrt(quad)*t(beta)%*%(Sigma.pi-quad.pi/quad*Sigma)%*%solve(sqrtm(Sigma.pi))
dnorm(quad.pi/sqrt(quad))*1/sqrt(quad)*t(beta)%*%(Sigma.pi-quad.pi/quad*Sigma)%*%solve(sqrtm(Sigma.pi))

dd


## checking z-stats with adversarial parameters
require(mvtnorm)
require(parallel)
require(numDeriv)
source('misc.R')
start <- Sys.time()
set.seed(1)
p <- 8
target <- 1/10 ## derivative part
## pi.0 <- runif(1); pi.0 <- max(pi.0,1-pi.0)
## pi.1 <- 1-pi.0
## a <- p*(4*sqrt(pi*exp(1))/(pi.0-pi.1)*target)^2
## b <- 0#((pi.0-pi.1)/3/sqrt(p))^(2/3)-2
## c <- b+2
## d <- sqrt(1+(b+c)/2)/sqrt(a*p)/(pi.0+pi.1*(b+c)/2)
## Sigma.0 <- diag(p) * a
## Sigma.1 <- diag(c(rep(b,p/2),rep(c,p/2))) * a
## beta.star <- rep(1,p) * d
uu <- 4*p*pi*exp(1)*target^2
pi.0.min <- 1/2*(1-uu+sqrt(uu*(uu+2)))
stopifnot(pi.0.min>1/2)
pi.0 <- runif(1,pi.0.min)
pi.1 <- 1-pi.0
d <- sqrt(2/p)
a <- (1-2*pi.0)^2/16/p/pi/exp(1)
b <- -2*a-pi.0*target^2
c <- a+(pi.0-1)*target^2
## f <- function(x)a*x^2+b*x+c
eps <- min(1/2/a*(-b+c(-1,1)*sqrt(b^2-4*a*c)))
stopifnot(eps>0)
beta.star <- rep(1,p )*d
## mu.0 <- rep(0,p)
## mu.1 <- (pi.0*Sigma.0+pi.1*Sigma.1)%*%beta
## params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
## auc.deriv.lda.gaussian(beta,params)
Sigma.1 <- diag(p) 
Sigma.0 <- diag(c(rep(eps,p/2),rep(2-eps,p/2)) )
Sigma.pi <- pi.0*Sigma.0+pi.1*Sigma.1
Sigma <- Sigma.0+Sigma.1
mu.0 <- rep(0,p)
mu.1 <- (pi.0*Sigma.0+pi.1*Sigma.1)%*%beta.star
params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
## beta.star <- solve(pi.0*Sigma.0+pi.1*Sigma.1)%*%(mu.1-mu.0)
deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
## deriv.star%*%solve(expm::sqrtm(Sigma.pi)) ## approx 2
## ns <- round(seq(1e2,1e3,len=20))
## by.n <- sapply(ns, function(n) {
auc <- auc.lda.gaussian(beta.star,params)
n <- 5e3
n.0 <- round(n*pi.0); n.1 <- n-n.0
## z.stats <- replicate(5e2, {
z.stats <- mclapply(1:1e2, mc.cores=detectCores()-3, FUN= function(dd){
    x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
    x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
    beta.hat <- coefs.lda(x.0,x.1,params=list())
    params.hat <- list(mu.0=colMeans(x.0),mu.1=colMeans(x.1),pi.0=nrow(x.0)/(nrow(x.0)+nrow(x.1)),Sigma.0=var(x.0),Sigma.1=var(x.1))
    x <- rbind(x.0,x.1)
    g <- c(rep(0,n.0),rep(1,n.1))
    ## approx.1 <- auc.hajek.lda.gaussian(x,g,beta.hat,params.hat,terms.only=TRUE,IID=TRUE)
    approx.1 <- auc.hajek(x.0%*%beta.hat,x.1%*%beta.hat,terms.only=TRUE,IID=TRUE)
    infl <- infl.lda(x,g,params.hat,terms.only=TRUE)
    ## infl <- infl.lda(x,g,params,terms.only=TRUE)##!!
    ## deriv.star.hat <- rep(0,p)
    deriv.star.hat <- deriv.star ##!!!
    ## deriv.star.hat <- grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat)
    approx.2 <- deriv.star.hat%*%infl
    obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc
    approx <- as.numeric(approx.1+approx.2)
    obs / sqrt(var(approx) / length(approx))
})
z.stats <- unlist(z.stats)
Sys.time() - start
op <- par(mfrow=c(1,2))
qqnorm(z.stats)
abline(0,1)
hist(z.stats,prob=TRUE)
curve(dnorm,add=TRUE)
par(op)
## fpr
alpha <- .05
mean(abs(z.stats)<qnorm(1-alpha/2))
## L_inf distance from standard normal cdf
z.stats <- sort(z.stats)
distances <- pnorm(z.stats) - 1:length(z.stats)/length(z.stats)
max.dist <- max(abs(distances))
max.dist


dd



## as above, but comparing derivatives: 0, oracle, estimated. TODO:
## make this a comparison with a true oracle estimator, right now only
## the derivative is oracle.

require(mvtnorm)
require(parallel)
require(numDeriv)
source('misc.R')
sim <- function(n,p,taylor.part,alpha=.05,reps=3e2) {
    ## browser()
    uu <- 4*p*pi*exp(1)*taylor.part^2
    pi.0.min <- 1/2*(1-uu+sqrt(uu*(uu+2)))
    stopifnot(pi.0.min>1/2)
    pi.0 <- runif(1,pi.0.min)
    pi.1 <- 1-pi.0
    d <- sqrt(2/p)
    a <- (1-2*pi.0)^2/16/p/pi/exp(1)
    b <- -2*a-pi.0*taylor.part^2
    c <- a+(pi.0-1)*taylor.part^2
    eps <- min(1/2/a*(-b+c(-1,1)*sqrt(b^2-4*a*c)))
    stopifnot(eps>0)
    beta.star <- rep(1,p )*d
    Sigma.1 <- diag(p) 
    Sigma.0 <- diag(c(rep(eps,p/2),rep(2-eps,p/2)) )
    Sigma.pi <- pi.0*Sigma.0+pi.1*Sigma.1
    Sigma <- Sigma.0+Sigma.1
    mu.0 <- rep(0,p)
    mu.1 <- (pi.0*Sigma.0+pi.1*Sigma.1)%*%beta.star
    params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0)
    deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
    auc <- auc.lda.gaussian(beta.star,params)
    ## n <- 5e3
    n.0 <- round(n*pi.0); n.1 <- n-n.0
    ## z.stats <- replicate(5e2, {
    z.stats <- mclapply(1:reps, mc.cores=detectCores()-3, FUN= function(dd){
        x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
        x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
        beta.hat <- coefs.lda(x.0,x.1,params=list())
        params.hat <- list(mu.0=colMeans(x.0),mu.1=colMeans(x.1),pi.0=nrow(x.0)/(nrow(x.0)+nrow(x.1)),Sigma.0=var(x.0),Sigma.1=var(x.1))
        x <- rbind(x.0,x.1)
        g <- c(rep(0,n.0),rep(1,n.1))
        ## approx.1 <- auc.hajek.lda.gaussian(x,g,beta.hat,params.hat,terms.only=TRUE,IID=TRUE)
        approx.1 <- auc.hajek(x.0%*%beta.hat,x.1%*%beta.hat,terms.only=TRUE,IID=TRUE)
        infl <- infl.lda(x,g,params.hat,terms.only=TRUE)
        deriv.star.hat.1 <- rep(0,p)
        deriv.star.hat.2 <- deriv.star
        deriv.star.hat.3 <- grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat)
        sapply(list(deriv.star.hat.1=deriv.star.hat.1,deriv.star.hat.2=deriv.star.hat.2,deriv.star.hat.3=deriv.star.hat.3), function(deriv.star.hat) {
            approx.2 <- deriv.star.hat%*%infl
            obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc
            approx <- as.numeric(approx.1+approx.2)
            obs / sqrt(var(approx) / length(approx))        
        })    
    })
    z.stats <- simplify2array(z.stats)
    fpr <- rowMeans(abs(z.stats)<qnorm(1-alpha/2))
    ## L_inf distance from standard normal cdf
    max.dist <- apply(z.stats,1,function(knots) {
        knots <- sort(knots)
        distances <- pnorm(knots) - 1:length(knots)/length(knots)
        max(abs(distances))
    })
    cbind(fpr=fpr,max.dist=max.dist)
}
start <- Sys.time()
set.seed(1)
p <- 8
## n <- 5e3
taylor.part <- 1/10
ns <- round(seq(3e2,5e3,len=30))
ns <- round(seq(3e2,3e3,len=10))
reps <- 1e2#1e3
alpha <- .05
sim0 <- lapply(ns, function(n) sim(n,p,taylor.part,alpha=alpha,reps=reps))
fpr <- do.call(rbind,lapply(sim0,function(m)m[,'fpr']))
max.dist <- do.call(rbind,lapply(sim0,function(m)m[,'max.dist']))
Sys.time() - start
matplot(x=ns,fpr,pch=1,type='l',col=1:3,lty=1)
legend('bottomleft',legend=c('ignored','oracle','estimated'),col=1:3,lty=1)
abline(h=1-alpha,lty=2)
## matplot(max.dist,pch=1,type='l',col=1:3,lty=1)
## legend('bottomleft',legend=c('ignored','oracle','estimated'),col=1:3,lty=1)
## abline(h=1-alpha,lty=2)
## save.image('sessions/15e.RData')




## 16 redoing 14e, logistic + guassian

## 16a checking z-stats on reduced data
require(mvtnorm)
require(numDeriv)
source('misc.R')
## set.seed(21)
p.full <- 5
p.red <- 3
n <- 1e3
mu <- runif(p.full)
Sigma <- matrix(runif(p.full^2),p.full)*4
Sigma <- Sigma%*%t(Sigma)
## p.full <- 4
## p.red <- 2
## n <- 3e3
## mu <- rep(0,p.full)
Sigma <- diag(p.full)*1.5+2
beta.0 <- runif(p.full)
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
params.full <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2)
beta.star <- as.numeric(coef.reduced.glm.gaussian(p.red,params.full,lim=5))
## beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
params.red <- list(mu=mu[1:p.red],Sigma=Sigma[1:p.red,1:p.red],beta=beta.star,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
auc.red <- auc.glm.gaussian(beta.star,params.full)
deriv.star <- grad(function(beta)auc.glm.gaussian(c(beta[1:2],0,0),params.full,10),beta.star)
## ns <- round(seq(1e2,1e3,len=20))
## by.n <- sapply(ns, function(n) {
## set.seed(2)
z.stats <- replicate(1e2, {
    cat('.')
    x <- rmvnorm(n,mu,Sigma)
    risk <- h(x%*%beta.0)
    g <- rbinom(n,1,risk)
    x.red <- x[,1:p.red]
    x.0 <- x[g==0,]; x.1<- x[g==1,]    
    x.0.red <- x.red[g==0,]; x.1.red <- x.red[g==1,]    
    beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link=link.name)))
    ## beta.hat.red <- c(beta.hat.red,rep(0,p.full-p.red))
    auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
    obs <- auc.hat.red.hat - auc.red
    auc.hajek.red <- auc.hajek(x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red,F=NULL,G=NULL,auc=auc.hat.red.hat,terms.only=TRUE,IID=TRUE)
    params.hat.red <- list(mu=colMeans(x.red),Sigma=cov(x.red),beta=beta.hat.red,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
    infl.hat <- infl.glm.gaussian(x.red,g,params.hat.red)
    deriv.red <- grad(function(beta)auc.hat(x.0.red%*%beta,x.1.red%*%beta),beta.hat.red[1:p.red])
    deriv.red <- deriv.star ##!!!!
    deriv.red <- rep(0,p.red)
    approx <- auc.hajek.red + as.numeric(deriv.red%*%infl.hat) 
    obs / sqrt(var(approx) / length(approx))
    ## deriv.red[3]
})
op <- par(mfrow=c(1,2))
qqnorm(z.stats)
abline(0,1)
hist(z.stats,prob=TRUE)
curve(dnorm,add=TRUE)
par(op)
## fpr
alpha <- .05
mean(abs(z.stats)<qnorm(1-alpha/2))
## L_inf distance from standard normal cdf
z.stats <- sort(z.stats)
distances <- pnorm(z.stats) - 1:length(z.stats)/length(z.stats)
max.dist <- max(abs(distances))
max.dist


## ts <- seq(-1,1,len=200)/10
## a <- runif(2)
## vals <- sapply(ts,function(t)auc.hat(x.0.red%*%(beta.star.red[1:2]+t*a),x.1.red%*%(beta.star.red[1:2]+t*a)))
## plot(ts,vals)
## abline(v=0)





## 16b adversarial param search
source('misc.R')
require(numDeriv)
seed <- seed+1
set.seed(seed)
start <- Sys.time()
p.full <- 4
p.red <- 3
## n <- 3e3
## mu <- runif(p.full)
mu <- rep(0,p.full)+1
## Sigma <- matrix(runif(p.full^2),p.full)
## Sigma <- Sigma%*%t(Sigma)*7
## Sigma <- diag(p.full)*3+3
## Sigma <- diag(runif(p.full))+4
## Sigma <- diag(c(rep(0,p.full/2),rep(3,p.full/2)))+3
Sigma <- diag(c(3,rep(1,p.full-1)))
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
## beta.0 <- runif(p.full)
beta.0 <- rep(1,p.full)
params.full <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2)
beta.star.red <- as.numeric(coef.reduced.glm.gaussian(p.red,params.full,lim=10))
beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
## params.red <- list(mu=mu[1:p.red],Sigma=Sigma[1:p.red,1:p.red],beta=beta.star.red,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
auc.red <- auc.glm.gaussian(beta.star.red,params.full)
deriv.star <- grad(function(beta)auc.glm.gaussian(c(beta[1:p.red],rep(0,p.full-p.red)),params.full,10),beta.star.red)
print(deriv.star[1:p.red])
print(Sys.time() - start)
## a <- runif(2)
## a <- c(1,0)
## ts <- seq(-1,1,len=10)/10
## vals <- sapply(ts, function(t)auc.glm.gaussian(c(beta.star.red[1:2]+t*a,0,0),params.full,5))
## plot(ts,vals,asp=1)
## abline(v=0)


x <- rmvnorm(1e4,mu,Sigma)
mean((as.numeric(h(x%*%beta.0) - h(x%*%beta.star.red))*x)[,4])

Sigma.11 <- runif(1)
mu <- c(1,2)
## Sigma <- matrix(c(Sigma.11,mu[2]/mu[1]*Sigma.11,mu[2]/mu[1]*Sigma.11,(mu[2]/mu[1])^2*Sigma.11),2)
gamma <- runif(2)
Sigma <- matrix(c(Sigma.11,mu[2]/mu[1]*Sigma.11,mu[2]/mu[1]*Sigma.11,(mu[2]/mu[1])^2*Sigma.11),2)
Sigma.11^2*(2*gamma[1]/gamma[2]*mu[2]/mu[1]*(1-gamma[2])+(mu[2]/mu[1])^2) - ((gamma[2]*mu[2]/mu[1])*Sigma.11)^2

beta <- gamma[1]+gamma[2]/mu[1]*mu[2]
beta[1]*mu[1] - sum(gamma*mu)
beta*Sigma.11 - (gamma[1]*Sigma.11 + gamma[2]*Sigma[1,2])
beta^2*Sigma.11 - t(gamma)%*%Sigma%*%gamma

x <- rmvnorm(1e4,mu,Sigma)
h(beta*x[,1])-h(x%*%gamma)

risk <- h(x%*%gamma)



## require(mvtnorm)
## require(numDeriv)
## source('misc.R')
## set.seed(1)
## p <- 4
## mu <- rep(0,p)
## Sigma <- matrix(runif(p^2),p)
## Sigma <- Sigma%*%t(Sigma)
## h <- plogis
## beta <- runif(p)
## quad <- as.numeric(t(beta)%*%Sigma%*%beta)
## ## inner <- Vectorize(function(u)integrate(function(v)h(-v)*dnorm(v,mean=0,sd=sqrt(quad))*h(u)*dnorm(u,0,sqrt(quad)),-Inf,u)$val)
## ## auc <- 4 * integrate(function(u)inner(u),-Inf,Inf)$val
## E.g <- function(beta)integrate(function(u)h(u)*dnorm(u,0,sqrt(quad(beta))),-Inf,Inf)$val
## quad <- function(beta)as.numeric(t(beta)%*%Sigma%*%beta)
## inner <- Vectorize(function(u,beta)integrate(function(v)h(-v)*dnorm(v,mean=0,sd=sqrt(quad(beta)))*h(u)*dnorm(u,0,sqrt(quad(beta))),-Inf,u)$val,vectorize.args='u')
## auc <- function(beta) 4 * integrate(function(u)inner(u,beta),-Inf,Inf)$val
## ## check auc formula
## n <- 1e2
## auc.hats <- replicate(1e3, {
##     x <- rmvnorm(n,mu,Sigma)
##     risk <- h(x%*%beta)
##     y <- rbinom(n,1,risk)
##     auc.hat(x[y==0,]%*%beta, x[y==1,]%*%beta)
## })
## hist(auc.hats)
## abline(v=auc(beta),col=2)
## abline(v=mean(auc.hats),col=3)

## inner <- Vectorize(function(u)integrate(function(v)h(-v)*h(u)*dnorm(v,mean=0,sd=sqrt(quad))*dnorm(u,0,sqrt(quad))*((u^2+v^2)/2/quad^2 - 1/quad),-Inf,u)$val)
## 4 * integrate(function(u)inner(u),-Inf,Inf)$val

## grad(auc,x=beta)

## ts <- seq(-1,1,len=30)
## a <- runif(p)
## vals <- sapply(ts, function(t)auc(beta+t*a))
## plot(ts,vals); abline(v=0)


## hessian

require(mvtnorm)
require(numDeriv)
source('misc.R')
set.seed(1)
p <- 4
## p.red <- 2
## n <- 3e2
mu <- rep(0,p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta <- runif(p)/2
params <- list(mu=mu,Sigma=Sigma,beta=beta,link=h)
## hessian(function(beta)auc.glm.gaussian(beta,params,lim=10),x=beta)
ts <- seq(-1,1,len=30)
a <- runif(p)
vals <- sapply(ts, function(t)auc.glm.gaussian(beta+t*a,params))
plot(ts,vals); abline(v=0)


require(mvtnorm)
p <- 4
## p.red <- 2
n <- 1e3
## beta <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
mu.0 <- rep(0,p)
mu.1 <- runif(p)/3#Sigma%*%beta
pi.0 <- .3; pi.1 <- 1-pi.0
## mu <- pi.0*mu.0+pi.1*mu.1
beta <- as.numeric(solve(Sigma)%*%(mu.1-mu.0))
beta.hats <- replicate(1e3, {
    n.0 <- round(rbinom(1,n,pi.0))
    n.1 <- n-n.0
    x.0 <- rmvnorm(n.0,mu.0,Sigma)
    x.1 <- rmvnorm(n.1,mu.1,Sigma)
    x <- rbind(x.0,x.1)
    y <- rep(c(0,1),c(n.0,n.1))
    beta.hat <- coef(glm(y~x,family=binomial('logit')))[-1]
})
a <- runif(p)
hist(a%*%beta.hats)
abline(v=beta%*%a,col=2)
abline(v=mean(a%*%beta.hats),col=3)



require(mvtnorm)
p <- 4
p.red <- 2
n <- 1e3
## beta <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
mu.0 <- rep(0,p)
mu.1 <- runif(p)/3#Sigma%*%beta
pi.0 <- .1; pi.1 <- 1-pi.0
## mu <- pi.0*mu.0+pi.1*mu.1
beta <- as.numeric(solve(Sigma)%*%(mu.1-mu.0))
beta.red <- as.numeric(solve(Sigma[1:p.red,1:p.red])%*%((mu.1-mu.0)[1:p.red]))
beta.0 <- log(pi.1/pi.0)-1/2*(t(mu.1)%*%solve(Sigma)%*%mu.1 - t(mu.0)%*%solve(Sigma)%*%mu.0)
beta.hats <- replicate(1e3, {
    n.0 <- round(rbinom(1,n,pi.0))
    n.1 <- n-n.0
    x.0 <- rmvnorm(n.0,mu.0,Sigma)
    x.1 <- rmvnorm(n.1,mu.1,Sigma)
    x <- rbind(x.0,x.1)
    x.red <- x[,1:p.red]
    y <- rep(c(0,1),c(n.0,n.1))
    beta.hat <- coef(glm(y~x.red,family=binomial('logit')))[-1]
})
a <- runif(p.red)
hist(a%*%beta.hats)
abline(v=beta.red%*%a,col=2)
abline(v=mean(a%*%beta.hats),col=3)


## formulas--likelihood of lda has the same lik as logistic regression
## with a gaussian mixture for x
y <- rbinom(1,1,.5)
x <- rmvnorm(1,(mu.0+mu.1)/2,Sigma)
dmvnorm(x,mu.0,Sigma)*pi.0*exp(y*(beta.0+sum(beta*x)))
(dmvnorm(x,mu.1,Sigma)*pi.1)^y*(dmvnorm(x,mu.0,Sigma)*pi.0)^(1-y) 
y
dmvnorm(x,mu.1,Sigma) / dmvnorm(x,mu.0,Sigma) * (pi.1/pi.0)
exp(beta.0+sum(beta*x))




p <- 4
p.red <- 2
n <- 1e3
## beta <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
mu.0 <- rep(0,p)
mu.1 <- runif(p)/3#Sigma%*%beta
pi.0 <- .1; pi.1 <- 1-pi.0
## mu <- pi.0*mu.0+pi.1*mu.1
beta <- as.numeric(solve(Sigma)%*%(mu.1-mu.0))
beta.red <- as.numeric(solve(Sigma[1:p.red,1:p.red])%*%((mu.1-mu.0)[1:p.red]))
beta.0 <- log(pi.1/pi.0)-1/2*(t(mu.1)%*%solve(Sigma)%*%mu.1 - t(mu.0)%*%solve(Sigma)%*%mu.0)
x <- runif(p)
f <- function(x)pi.0*dmvnorm(x,mu.0,Sigma)+pi.1*dmvnorm(x,mu.1,Sigma)
risk <- plogis(beta.0+sum(beta*x))
pi.0/pi.1*risk*dmvnorm(x,mu.0,Sigma)+risk*dmvnorm(x,mu.1,Sigma)
dmvnorm(x,mu.1,Sigma)
## dmvnorm(x,mu.1,Sigma)^2 / (1+exp(beta.0+sum(beta*x)))
1 / (1+exp(beta.0+sum(beta*x)))*(2*pi)^(-p/2)/sqrt(det(Sigma))*exp(-1/2*( t(x-mu.0)%*%solve(Sigma)%*%(x-mu.0)+t(x-mu.1)%*%solve(Sigma)%*%(x-mu.1)) + log(pi.0/pi.1)+beta.0+sum(beta*x))#t(mu.1-mu.0)%*%solve(Sigma)%*%x)





## logitnorm mean when mu=0; then ln.mean(-x)=1-ln.mean(x)
ln.mean <- function(x,beta,params,lim=Inf) {
    mu <- params$mu; Sigma <- params$Sigma; beta.0 <- params$beta
    h <- params$h
    quad.ww <- as.numeric(t(beta)%*%Sigma%*%beta)
    quad.ww0 <- as.numeric(t(beta.0)%*%Sigma%*%beta)
    quad.w0w0 <- as.numeric(t(beta.0)%*%Sigma%*%beta.0)
    Sigma.cond <- sqrt((1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0) 
    sapply(x, function(x)integrate(function(w0)h(w0)*dnorm(w0,mean=sum(beta.0*mu) + quad.ww0/quad.ww*(x-sum(beta*mu)), sd=Sigma.cond),-lim,lim)$val)
}
require(mvtnorm)
source('misc.R')
set.seed(1)
p <- 4
## p.red <- 2
n <- 1e3
mu <- rep(0,p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta.0 <- runif(p)/2
beta <- runif(p)
## beta <- beta.0
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,h=h)
curve(ln.mean(-x,beta,params),-10,10)
curve(1-ln.mean(x,beta,params),add=TRUE,lty=2)



## 16c logitnorm mean--changing convolution in integral, more stable numerically
ln.mean <- function(w,beta,params,lim=Inf) {
    ## browser()
    mu <- params$mu; Sigma <- params$Sigma; beta.0 <- params$beta
    h <- params$h
    if(all.equal(beta,beta.0)==TRUE)return(h(w))
    quad.ww <- as.numeric(t(beta)%*%Sigma%*%beta)
    quad.ww0 <- as.numeric(t(beta.0)%*%Sigma%*%beta)
    quad.w0w0 <- as.numeric(t(beta.0)%*%Sigma%*%beta.0)
    mu.cond <- function(w)beta.0%*%mu + quad.ww0/quad.ww*(w-beta%*%mu)
    sigma.cond <- sqrt((1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0) 
    ## sapply(w, function(w)integrate(function(w0)plogis(w0)*dnorm(w0,mean=mu.cond(w), sd=sigma.cond),-lim,lim)$val)
    sapply(w, function(w)integrate(function(w0)h(mu.cond(w)-w0)*dnorm(w0/sigma.cond)/sigma.cond,-lim,lim)$val)
}
require(mvtnorm)
source('misc.R')
set.seed(1)
p <- 4
## p.red <- 2
n <- 1e3
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta.0 <- runif(p)/2
beta <- runif(p)
## beta <- beta.0
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h)
curve(ln.mean(-x,beta,params),-10,10)

## ln(beta) -> plogis as beta -> beta.0
plot(0,type='n',ylim=c(0,1),xlim=c(-1,1))
ts <- seq(0,3,len=10)/5
a <- runif(p)
for(i in 1:length(ts))curve(ln.mean(x,beta.0+ts[i]*a,params),add=TRUE,col=terrain.colors(length(ts))[i])
## curve(h,add=TRUE,col=terrain.colors(length(ts))[1],lty=2)
legend('bottomright',lty=1,col=terrain.colors(length(ts)),legend=ts)

## 16c-1 sensitivity of ln.mean to mu. seems possible to use mu to separate
## out h terms, and sigma then to concentrate the densities at
## different h terms
require(mvtnorm)
source('misc.R')
set.seed(1)
p <- 4
## p.red <- 2
n <- 1e3
mu <- runif(p)*0
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)/10
h <- plogis
beta.0 <- runif(p)
## beta <- runif(p)
a <- rep(1,p)
t <- .2
beta <- beta.0+a*t
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h)
## curve(ln.mean(x,beta.0,params)*(1-ln.mean(x,beta.0,params)),-5,5) # wrong fla
## f <- function(x)ln.mean(x,beta.0+a*t,params)
## curve(f(x)*(1-f(x)),add=TRUE,col=2)
## abline(v=0,lty=2)
## abline(h=1/2,lty=2)
## just one term of product
curve(ln.mean(x,beta.0,params),-5,5)
curve(dnorm(x,mean=beta.0%*%mu,sd=sqrt(t(beta.0)%*%Sigma%*%beta.0)),add=TRUE,lty=2)
f <- function(x)ln.mean(x,beta.0+a*t,params)
curve(f(x),add=TRUE,col=2)
curve(dnorm(x,mean=beta%*%mu,sd=sqrt(t(beta)%*%Sigma%*%beta)),add=TRUE,lty=2,col=2)
abline(v=0,lty=2)
abline(h=1/2,lty=2)
## abline(v=beta.0%*%mu)
## abline(v=beta%*%mu,col=2)

## 16c-2 see if translates to auc. [doesnt seem to]
ln.mean <- function(w,beta,params,lim=Inf) {
    ## browser()
    mu <- params$mu; Sigma <- params$Sigma; beta.0 <- params$beta
    h <- params$h
    if(all.equal(beta,beta.0)==TRUE)return(h(w))
    quad.ww <- as.numeric(t(beta)%*%Sigma%*%beta)
    quad.ww0 <- as.numeric(t(beta.0)%*%Sigma%*%beta)
    quad.w0w0 <- as.numeric(t(beta.0)%*%Sigma%*%beta.0)
    mu.cond <- function(w)beta.0%*%mu + quad.ww0/quad.ww*(w-beta%*%mu)
    sigma.cond <- sqrt((1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0) 
    ## sapply(w, function(w)integrate(function(w0)plogis(w0)*dnorm(w0,mean=mu.cond(w), sd=sigma.cond),-lim,lim)$val)
    sapply(w, function(w)integrate(function(w0)h(mu.cond(w)-w0)*dnorm(w0/sigma.cond)/sigma.cond,-lim,lim)$val)
}
auc <- function(beta,params,lim=Inf) {
    ## pdf.index.glm.gaussian <- function(x,g,beta,params,lim=Inf) {
    ## browser()
    h <- params$link
    mu <- params$mu
    Sigma <- params$Sigma
    beta.0 <- params$beta
    f.w <- function(w)dnorm(w,t(beta)%*%mu, sd=sqrt(t(beta)%*%Sigma%*%beta))
    f.w0 <- function(w0)dnorm(w0,t(beta.0)%*%mu, sd=sqrt(t(beta.0)%*%Sigma%*%beta.0))
    pi.1 <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    inner <- Vectorize(function(v) integrate(function(w)(1-ln.mean(w,beta,params))*f.w(w),-lim,v)$val)
    integrate(function(v)inner(v)*ln.mean(v,beta,params)*f.w(v),-lim,lim)$val / (pi.1*(1-pi.1))
}
p <- 4
## p.red <- 2
n <- 1e3
mu <- runif(p)*3
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
h <- plogis
beta.0 <- runif(p)
beta <- runif(p)/30
## beta <- beta.0
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h)
## ## make sure above version of auc looks right on data
## try <- auc(beta,params)
## auc.hats <- replicate(1e2, {
##     x <- rmvnorm(n,mu,Sigma)
##     risk <- h(x%*%beta.0)
##     d <- rbinom(n,1,risk)
##     x.0 <- x[d==0,]; x.1 <- x[d==1,]
##     auc.hat(x.0%*%beta.0,x.1%*%beta.0)
## })
## hist(auc.hats)
## abline(v=try,col=2)
## abline(v=mean(auc.hats),col=3)
## check auc curve as beta varies around beta.0
## a <- rep(1,p)
ts <- seq(0,1,len=5)
aucs <- sapply(ts, function(t)auc(beta.0*t+beta*(1-t),params))
plot(ts,aucs,type='l',asp=1)
## save.image('sessions/16c-2.RData') ##randomly found with large slope


dd


## 16c-3 see why auc doesnt change much even when ln.means do
require(mvtnorm)
source('misc.R')
set.seed(1)
p <- 4
## p.red <- 2
n <- 1e3
mu <- runif(p)*1
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)/20
h <- plogis
beta.0 <- runif(p)
## beta <- beta.0
a <- runif(p)-.4
t <- .2
beta <- runif(p)
## beta <- beta.0+a*t
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h)
f.w <- function(w,beta)dnorm(w,t(beta)%*%mu, sd=sqrt(t(beta)%*%Sigma%*%beta))
f <- function(w,v,beta)(1-ln.mean(w,beta,params))*ln.mean(v,beta,params)*f.w(w,beta)*f.w(v,beta)
## f <- function(w,v,beta)f.w(w,beta)*f.w(v,beta)
## f <- function(w,v,beta)(1-ln.mean(w,beta,params))*ln.mean(v,beta,params)
## f <- function(w,v,beta)(1-h(w))*h(v)
grid <- seq(-1,1,len=20)*3
z <- outer(grid,grid,f,beta)
contour(grid,grid,z,asp=1)
z.0 <- outer(grid,grid,f,beta.0)
contour(grid,grid,z.0,add=TRUE,col=2)
abline(0,1)

f <- function(w,v,beta)f.w(w,beta)*f.w(v,beta)
z <- outer(grid,grid,f,beta.0+a*t)
contour(grid,grid,z,add=TRUE)
z.0 <- outer(grid,grid,f,beta.0)
contour(z,add=TRUE,col=2)
abline(0,1)

c(auc(beta.0,params),auc(beta,params))


## side issue noticed: seems auc \neq 1/2 when beta.0=0? works with
## noise. so must be related to estimation of the parameter? but
## derivative term is 0 as the model is well specified. Seems possible
## it is just a finite-sample issue.
require(mvtnorm)
n <- 1e4
p <- 4
mu <- rep(0,p)
Sigma <- diag(p)
h <- plogis
ns <- round(seq(1e3,2e4,len=10))
by.n <- sapply(ns, function(n) {
    auc.hats <- replicate(1e1, {
        x <- matrix(rnorm(n*p),ncol=p)#rmvnorm(n,mu,Sigma)
        d <- rbinom(n,1,1/2)
        beta.hat <- coef(glm(d~x-1,family=binomial('logit')))
        ## beta.hat <- rnorm(p)
        auc.hat(x[d==0,]%*%beta.hat,x[d==1,]%*%beta.hat)
        ## beta.hat[1]
    })
    ## hist(auc.hats)
    ## abline(v=mean(auc.hats),col=2)
    mean(auc.hats)
})
plot(ns,by.n,type='o')
lm0 <- lm(log(by.n)~log(ns))
coef(lm0)
lines(ns,exp(fitted(lm0)),col=2)






## 16d taking beta out of the picture, just mu and sigma
## logit normal mean of a normal rv with mean and variance given by y | x
h.cond <- function(x,mu,Sigma,h,lim=Inf) {
    mu.cond <- function(x)mu[2]+Sigma[1,2]/Sigma[1,1]*(x-mu[1])
    sigma.cond <- sqrt((1-Sigma[1,2]^2/prod(diag(Sigma)))*Sigma[2,2])
    sapply(x, function(x)integrate(function(y)h(mu.cond(x)-y)*dnorm(y/sigma.cond)/sigma.cond,-lim,lim)$val)
}
auc <- function(mu,Sigma,lim=Inf) {
    ## browser()
    f.x <- function(x,mu,Sigma)dnorm(x,mu[1],sd=sqrt(Sigma[1,1]))
    f <- function(w,v,mu,Sigma)(1-h.cond(w,mu,Sigma,h,lim))*h.cond(v,mu,Sigma,h,lim)*f.x(w,mu,Sigma)*f.x(v,mu,Sigma)
    ## pi.1 <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    f.y <- function(y,mu,Sigma)dnorm(y,mu[2],sd=sqrt(Sigma[2,2]))
    pi.1 <- integrate(function(x)h(x)*f.y(x,mu,Sigma),-Inf,Inf)$val
    inner <- Vectorize(function(v)integrate(function(w)f(w,v,mu,Sigma),-lim,v)$val)
    integrate(function(v)inner(v),-lim,lim)$val / (pi.1*(1-pi.1))
}
require(mvtnorm)
source('misc.R')
set.seed(1)
## n <- 1e3
mu <- runif(2)*0
Sigma <- matrix(runif(4),2)
Sigma <- Sigma%*%t(Sigma)
sigma.12 <- .5
sigma.22 <- 1
Sigma <- matrix(c(1,sigma.12,sigma.12,sigma.22^2),2)
#### match with beta-generated coefficients
## mu <- runif(p)*1
## Sigma <- matrix(runif(p^2),p)
## Sigma <- Sigma%*%t(Sigma)/20
## h <- plogis
## beta.0 <- runif(p)
a <- c(0,1)#runif(2)
t <- -3
## beta <- runif(p)
## mu <- c(beta%*%mu,beta.0%*%mu)
## Sigma <- matrix(c(t(beta)%*%Sigma%*%beta,t(beta)%*%Sigma%*%beta.0,t(beta)%*%Sigma%*%beta.0,t(beta.0)%*%Sigma%*%beta.0),2)
###
h <- plogis
f.x <- function(x,mu,Sigma)dnorm(x,mu[1],sd=sqrt(Sigma[1,1]))
f <- function(w,v,mu,Sigma)(1-h.cond(w,mu,Sigma,h))*h.cond(v,mu,Sigma,h)*f.x(w,mu,Sigma)*f.x(v,mu,Sigma)
f <- function(w,v,mu,Sigma)f.x(w,mu,Sigma)*f.x(v,mu,Sigma)
f <- function(w,v,mu,Sigma)(1-h.cond(w,mu,Sigma,h))*h.cond(v,mu,Sigma,h)
grid <- seq(-1,1,len=50)*10
z <- outer(grid,grid,f,mu,Sigma)
contour(grid,grid,z,asp=1)#,xlim=range(grid),ylim=range(grid))
## contour(grid,grid,z,xlim=range(grid),ylim=range(grid))
## print(auc(mu,Sigma))
## Sigma[2,2] <- Sigma[2,2]+3
Sigma[1,2] <- Sigma[2,1] <- sigma.12+.3
## Sigma[2,2] <- Sigma[1,2]^2+.01
## mu[2] <- mu[2]+2
z <- outer(grid,grid,f,mu,Sigma/1)
contour(grid,grid,z,add=TRUE,col=2)
abline(0,1)
## print(auc(mu,Sigma))

deltas <- seq(.001,.1,len=4)
aucs <- sapply(deltas, function(delta) {
    Sigma[1,2]  <- Sigma[2,1] <- sigma.22-delta
    auc(mu,Sigma)
})
plot(deltas,aucs,asp=1,type='l')

auc <- function(beta,params,lim=Inf) {
    ## pdf.index.glm.gaussian <- function(x,g,beta,params,lim=Inf) {
    ## browser()
    h <- params$link
    mu <- params$mu
    Sigma <- params$Sigma
    beta.0 <- params$beta
    f.w <- function(w)dnorm(w,t(beta)%*%mu, sd=sqrt(t(beta)%*%Sigma%*%beta))
    f.w0 <- function(w0)dnorm(w0,t(beta.0)%*%mu, sd=sqrt(t(beta.0)%*%Sigma%*%beta.0))
    pi.1 <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    inner <- Vectorize(function(v) integrate(function(w)(1-ln.mean(w,beta,params))*f.w(w),-lim,v)$val)
    integrate(function(v)inner(v)*ln.mean(v,beta,params)*f.w(v),-lim,lim)$val / (pi.1*(1-pi.1))
}
p <- 4
mu.xy <- mu ## ran after previous section
mu <- rep(0,p)
beta <- c(sqrt(sigma.12),sqrt(1-sigma.12),rep(0,p-2))
beta.0 <- c(sqrt(sigma.12),0,sqrt(sigma.22-sigma.12),rep(0,p-3))
Sigma <- diag(p)
## mu <- rep(0,p)
## beta <- runif(p/2)
## beta <- beta / sqrt(sum(beta^2))
## beta <- c(beta,rep(0,p/2))
## ## unit <- runif(p/2); unit <- unit/sqrt(sum(unit^2))
## beta.0 <- runif(p)
## beta.0 <- beta.0/sqrt(sum(beta.0^2))*sqrt(2)
## Sigma <- diag(p)
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=plogis)
auc(beta,params)
## mu.2 <- 2
## mu[(p/2+1):p] <- beta.0[(p/2+1):p]/sum(beta.0[(p/2+1):p]^2)*mu.2
## params$mu <- mu
sigma.12 <- sigma.12 + .1
beta <- c(sqrt(sigma.12),sqrt(1-sigma.12),rep(0,p-2))
beta.0 <- c(sqrt(sigma.12),0,sqrt(sigma.22-sigma.12),rep(0,p-3))
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=plogis)
auc(beta,params)

## numDeriv::grad(function(beta)auc(beta,params),beta)
## a <- rep(1,p)
ts <- seq(0,1,len=10)
aucs <- sapply(ts, function(t)auc(t*beta.0+(1-t)*beta,params))




## 16e forcing beta==0
require(mvtnorm)
p <- 4
n <- 3e2
mu <- rep(0,p)
Sigma <- diag(p)
h <- plogis
## beta.0 <- c(rep(0,p/2),rep(1,p/2))
ratios <- seq(1,10,len=10)
by.ratio <- sapply(ratios, function(ratio) {
    beta.0 <- c(rep(1,p/2),rep(1,p/2)*ratio)/max(ratios)
    beta.0 <- c(rep(1,p/2)/ratio,rep(1,p/2))
    coefs <- replicate(1e3, {
        x <- rmvnorm(n,mu,Sigma)
        ## beta.0 <- rep(1,p)
        risk <- h(x%*%beta.0)
        d <- rbinom(n,1,risk)
        x.red <- x[,1:(p/2)]
        coef(glm(d~x.red-1,family=binomial('logit')))
    })
    mean(coefs[1,])
    ## hist(coefs[1,])
    ## abline(v=mean(coefs[1,]),col=2)
    ## abline(v=0,col=3)
})
plot(ratios,by.ratio)
lm(log(by.ratio) ~ log(ratios))

ln.mean <- function(w,beta,params,lim=Inf) {
    ## browser()
    mu <- params$mu; Sigma <- params$Sigma; beta.0 <- params$beta
    h <- params$h
    if(all.equal(beta,beta.0)==TRUE)return(h(w))
    quad.ww <- as.numeric(t(beta)%*%Sigma%*%beta)
    quad.ww0 <- as.numeric(t(beta.0)%*%Sigma%*%beta)
    quad.w0w0 <- as.numeric(t(beta.0)%*%Sigma%*%beta.0)
    mu.cond <- function(w)beta.0%*%mu + quad.ww0/quad.ww*(w-beta%*%mu)
    sigma.cond <- sqrt((1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0) 
    ## sapply(w, function(w)integrate(function(w0)plogis(w0)*dnorm(w0,mean=mu.cond(w), sd=sigma.cond),-lim,lim)$val)
    sapply(w, function(w)integrate(function(w0)h(mu.cond(w)-w0)*dnorm(w0/sigma.cond)/sigma.cond,-lim,lim)$val)
}
auc <- function(beta,params,lim=Inf) {
    ## pdf.index.glm.gaussian <- function(x,g,beta,params,lim=Inf) {
    ## browser()
    h <- params$link
    mu <- params$mu
    Sigma <- params$Sigma
    beta.0 <- params$beta
    f.w <- function(w)dnorm(w,t(beta)%*%mu, sd=sqrt(t(beta)%*%Sigma%*%beta))
    f.w0 <- function(w0)dnorm(w0,t(beta.0)%*%mu, sd=sqrt(t(beta.0)%*%Sigma%*%beta.0))
    pi.1 <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    inner <- Vectorize(function(v) integrate(function(w)(1-ln.mean(w,beta,params))*f.w(w),-lim,v)$val)
    integrate(function(v)inner(v)*ln.mean(v,beta,params)*f.w(v),-lim,lim)$val / (pi.1*(1-pi.1))
}
require(mvtnorm)
source('misc.R')
set.seed(1)
p <- 4
n <- 1e3
mu <- rep(0,p)
a <- 5
Sigma <- diag(p)*sqrt(a)
h <- plogis
## beta.0 <- c(rep(1,p/2),rep(1,p/2)*ratio)/max(ratios)
bs <- seq(1,30,len=20)
by.b <- sapply(bs, function(b) {
    ## beta.0 <- c(rep(1,p/2),rep(1,p/2)*b)#/max(bs)
    ## beta.0 <- c(rep(1,p/2)/b,rep(1,p/2))
    beta.0 <- c(rep(1,p/2)/b,rep(1,p/2))/a
    ## params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h)
    aucs <- replicate(3e1, {
        ## diag(Sigma)[(p/2+1):p] <- 10
        x <- rmvnorm(n,mu,Sigma)
        risk <- h(x%*%beta.0)
        d <- rbinom(n,1,risk)
        x.red <- x[,1:(p/2)]
        beta.hat <- coef(glm(d~x.red-1,family=binomial('logit')))
        beta.hat <- c(beta.hat,rep(0,p/2))
        ## d <- rbinom(n,1,1/2)
        ## beta.hat <- coef(glm(d~x-1,family=binomial('logit')))
        auc.hat(x[d==0,]%*%beta.hat,x[d==1,]%*%beta.hat)
        ## numDeriv::grad(function(beta)auc.hat(x[d==0,]%*%beta,x[d==1,]%*%beta),x=beta.hat)       
        ## beta.hat[2]
    })
    mean(aucs)
})
plot(bs,by.b)
lm(log(abs(by.b)) ~ log(bs))



## 16e-1 auc surface in p=2 case
require(mvtnorm)
source('misc.R')
set.seed(1)
p <- 2
n <- 1e3
mu <- rep(0,p)
a <- 5
Sigma <- diag(p)*sqrt(a)
h <- plogis
## beta.0 <- c(rep(1,p/2),rep(1,p/2)*ratio)/max(ratios)
b <- 5#seq(1,30,len=20)
## by.b <- sapply(bs, function(b) {
## beta.0 <- c(rep(1,p/2),rep(1,p/2)*b)#/max(bs)
## beta.0 <- c(rep(1,p/2)/b,rep(1,p/2))
beta.0 <- c(rep(1,p/2)/b,rep(1,p/2))/a
grid<- seq(-1,1,len=30)
## grid <- outer(grid.pts,grid.pts)
## beta <- runif(2)
zs <- replicate(1e1, {
    x <- rmvnorm(n,mu,Sigma)
    risk <- h(x%*%beta.0)
    d <- rbinom(n,1,risk)
    f <- Vectorize(function(beta1,beta2) {
        beta <- c(beta1,beta2)
        auc.hat(x[d==0,]%*%beta,x[d==1,]%*%beta)
    })
    outer(grid,grid,f)
},simplify=FALSE)
z <- Reduce('+',zs) / length(zs)
contour(grid,grid,z)
image(grid,grid,z)
contour(grid,grid,z,add=TRUE)
points(beta.0[1],beta.0[2],col=1,cex=3,pch=20)

## check against numerically calculated auc. roughly similar.
ln.mean <- function(x,beta,params,lim=Inf) {
    mu <- params$mu; Sigma <- params$Sigma; beta.0 <- params$beta
    quad.ww <- as.numeric(t(beta)%*%Sigma%*%beta)
    quad.ww0 <- as.numeric(t(beta.0)%*%Sigma%*%beta)
    quad.w0w0 <- as.numeric(t(beta.0)%*%Sigma%*%beta.0)
    Sigma.cond <- sqrt((1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0) 
    ln.mean <-  sapply(x, function(x)integrate(function(w0)plogis(w0)*dnorm(w0,mean=sum(beta.0*mu) + quad.ww0/quad.ww*(x-sum(beta*mu)), sd=Sigma.cond),-lim,lim)$val)
}
## ln.mean^g*(1-ln.mean)^(1-g)*f.w(w)/(E.g^g*(1-E.g)^(1-g))
auc.logit.gaussian <- function(beta,params,lim=Inf) { # ==auc.glm.gaussian.new from 14c
    mu <- params$mu; Sigma <- params$Sigma; beta.0 <- params$beta
    f.w <- function(w)dnorm(w,t(beta)%*%mu, sd=sqrt(t(beta)%*%Sigma%*%beta))
    f.w0 <- function(w0)dnorm(w0,t(beta.0)%*%mu, sd=sqrt(t(beta.0)%*%Sigma%*%beta.0))
    E.g <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    inner <- Vectorize(function(v,beta)integrate(function(w)(1-ln.mean(w,beta,params))*ln.mean(v,beta,params)*f.w(w)*f.w(v),-lim,v)$val,vectorize.args='v')
    integrate(function(v)inner(v,beta),-lim,lim)$val / (E.g*(1-E.g))
}
require(mvtnorm)
source('misc.R')
set.seed(1)
p <- 2
## n <- 1e3
mu <- rep(0,p)
a <- 5
Sigma <- diag(p)*sqrt(a)
h <- plogis
## beta.0 <- c(rep(1,p/2),rep(1,p/2)*ratio)/max(ratios)
b <- 5#seq(1,30,len=20)
beta.0 <- c(rep(1,p/2)/b,rep(1,p/2))/a
params <- list(mu=mu,Sigma=Sigma,beta=beta.0)
grid<- seq(-1,1,len=15)
f <- Vectorize(function(beta1,beta2) {
    cat('.')
    beta <- c(beta1,beta2)
    ## auc.hat(x[d==0,]%*%beta,x[d==1,]%*%beta)
    tryCatch(auc.logit.gaussian(beta,params,lim=10),error=function(e)NA)
})
z <- outer(grid,grid,f)
contour(grid,grid,z)
image(grid,grid,z)
points(beta.0[1],beta.0[2],col=1,cex=3,pch=20)
## save.image('sessions/16e-1.RData')

dd


## 16f checking z-stats on reduced data with the adversarial
## model. unadjusted delong model performs much better... update: main
## culprit was numDeriv::grad setting. now looks much closer to the
## unadjusted delong (though not better than it)
require(mvtnorm)
require(numDeriv)
require(parallel)
source('misc.R')
start <- Sys.time()
lim <- 15
set.seed(1)
p.full <- 4
p.red <- 2
mu <- rep(0,p.full)
a <- 10
b <- 10
Sigma <- diag(p.full)*sqrt(a)
beta.0 <- c(rep(1,p.red)/b,rep(1,p.full-p.red))/a 
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
params.full <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2)
beta.star <- as.numeric(coef.reduced.glm.gaussian(p.red,params.full,lim=lim))
## beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
## beta.star.hat <- rowMeans(replicate(1e3, {
##     x <- rmvnorm(n,mu,Sigma)
##     risk <- h(x%*%beta.0)
##     g <- rbinom(n,1,risk)
##     x.red <- x[,1:p.red]
##     x.0 <- x[g==0,]; x.1<- x[g==1,]    
##     x.0.red <- x.red[g==0,]; x.1.red <- x.red[g==1,]    
##     beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link=link.name)))
## }))
params.red <- list(mu=mu[1:p.red],Sigma=Sigma[1:p.red,1:p.red],beta=beta.star,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
## auc.red <- auc.glm.gaussian(beta.star,params.full,lim=lim)
## auc.red.hats <- replicate(1e4, {
##     x <- rmvnorm(n,mu,Sigma)
##     risk <- h(x%*%beta.0)
##     g <- rbinom(n,1,risk)
##     x.red <- x[,1:p.red]
##     x.0 <- x[g==0,]; x.1<- x[g==1,]    
##     x.0.red <- x.red[g==0,]; x.1.red <- x.red[g==1,]    
##     beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link=link.name)))
##     auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
## })
n <- 5e5
x <- rmvnorm(n,mu,Sigma)
risk <- h(x%*%beta.0)
g <- rbinom(n,1,risk)
x.red <- x[,1:p.red]
x.0 <- x[g==0,]; x.1<- x[g==1,]    
x.0.red <- x.red[g==0,]; x.1.red <- x.red[g==1,]    
beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link=link.name)))
auc.red <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
## auc.red <- mean(auc.red.hats) # maybe numerical issues. beta.star.red OK but auc.red seems off. estimating.
## save.image('sessions/16f-4.RData')
load('sessions/16f-4.RData')
set.seed(1)
start <- Sys.time()
## deriv.star <- grad(function(beta)auc.glm.gaussian(c(beta[1:p.red],rep(0,p.full-p.red)),params.full,10),beta.star)
## ns <- round(seq(1e2,1e3,len=20))
## by.n <- sapply(ns, function(n) {
## set.seed(2)
## z.stats <- replicate(1e2, {
n <- 1.5e3
z.stats <- mclapply(1:5e2, mc.cores=detectCores()-3, FUN=function(dd){
    cat('.')
    x <- rmvnorm(n,mu,Sigma)
    risk <- h(x%*%beta.0)
    g <- rbinom(n,1,risk)
    x.red <- x[,1:p.red]
    x.0 <- x[g==0,]; x.1<- x[g==1,]    
    x.0.red <- x.red[g==0,]; x.1.red <- x.red[g==1,]    
    beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link=link.name)))
    ## beta.hat.red <- c(beta.hat.red,rep(0,p.full-p.red))
    auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
    obs <- auc.hat.red.hat - auc.red
    auc.hajek.red <- auc.hajek(x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red,F=NULL,G=NULL,auc=auc.hat.red.hat,terms.only=TRUE,IID=TRUE)
    params.hat.red <- list(mu=colMeans(x.red),Sigma=cov(x.red),beta=beta.hat.red,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
    infl.hat <- infl.glm.gaussian(x.red,g,params.hat.red)
    ## deriv.red <- grad(function(beta)auc.hat(x.0.red%*%beta,x.1.red%*%beta),beta.hat.red[1:p.red])
    ## deriv.red <- grad(function(beta)auc.hat(x.0.red%*%beta,x.1.red%*%beta),beta.star)
    ## deriv.red <- deriv.star 
    deriv.red <- grad(function(beta)auc.hat(x.0.red%*%beta,x.1.red%*%beta),beta.hat.red[1:p.red],method='simple',method.args=list(eps=1/n^(.5)))
    ## deriv.red <- rep(0,p.red)
    approx <- auc.hajek.red + as.numeric(deriv.red%*%infl.hat) 
    obs / sqrt(var(approx) / length(approx))
})
z.stats <- simplify2array(z.stats)
Sys.time() - start
z.stats <- z.stats-mean(z.stats)#!!! centering seems off
op <- par(mfrow=c(1,2))
qqnorm(z.stats)
abline(0,1)
hist(z.stats,prob=TRUE)
curve(dnorm,add=TRUE)
par(op)
## fpr
alpha <- .05
print(mean(abs(z.stats)<qnorm(1-alpha/2)))
## print(mean(abs(z.stats-mean(z.stats))<qnorm(1-alpha/2)))#!!!
## L_inf distance from standard normal cdf
## z.stats <- sort(z.stats)
## distances <- pnorm(z.stats) - 1:length(z.stats)/length(z.stats)
## max.dist <- max(abs(distances))
## max.dist
var(z.stats)

dd



## below trying to track down problem with convergence. turned out to
## be mainly numDeriv::grad

## 16f-1 check influence function on adversarial model
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
n <- 1.5e4
mu <- rep(0,p.full)
a <- 1
b <- 1
Sigma <- diag(p.full)*sqrt(a)
beta.0 <- c(rep(1,p.red)/b,rep(1,p.full-p.red))/a
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
## beta.0 <- runif(p.full)/2
params.full <- list(mu=mu,Sigma=Sigma,beta=beta.0,p=p.full,link=h,link.deriv=h.1,link.deriv2=h.2)
beta.star.red <- coef.reduced.glm.gaussian(p.red,params.full,15)
params.red <- list(mu=mu[1:p.red],Sigma=Sigma[1:p.red,1:p.red],beta=beta.star.red,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
## beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
## ns <- round(seq(1e2,1e3,len=20))
## by.n <- sapply(ns, function(n) {
z.stats <-     replicate(4e2, {
    x <- rmvnorm(n,mu,Sigma)
    eta <- as.numeric(x%*%beta.0)
    risk <- h(eta)
    g <- rbinom(n,1,risk)
    ## beta.hat <- coef(glm(g~x-1,family=binomial(link=link.name)))
    ## obs <- beta.hat-beta.0
    ## approx <- infl.glm.gaussian(x,g,params.full,terms.only=FALSE)
    ## obs-approx
    x.red <- x[,1:p.red]
    beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link=link.name)))
    obs <- beta.hat.red - beta.star.red
    infl <- infl.glm.gaussian(x.red,g,params.red,terms.only=TRUE)
    var.hat <- var(t(infl))/n
    expm::sqrtm(solve(var.hat))%*%obs 
})
z.stats <- z.stats[,,]
op <- par(mfrow=c(1,2))
qqnorm(z.stats)
abline(0,1)
hist(z.stats,prob=TRUE)
curve(dnorm,add=TRUE)
par(op)




## 16f-2 check sd of unadjusted estimator with adversrial data. is too
## large, resulting in too small z-stats, doesn't seem to improve with
## n possibly getting worse.
auc.var <- function(x,y) {
    terms <- auc.hajek(x,y,terms.only=TRUE,IID=TRUE)
    var(terms) / length(terms)
}
## n <- 1e2
## standardized <- replicate(3e3, {
##     x <- runif(n); y <- runif(n)
##     theta.hat <- auc.hat(x,y)
##     var.hat <- auc.var(x,y)
##     theta.hat / sqrt(var.hat)
## })
## var(standardized)
require(mvtnorm)
require(numDeriv)
require(parallel)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
## n <- 5e3
mu <- rep(0,p.full)
a <- 1
b <- 1
h <- plogis
Sigma <- diag(p.full)*sqrt(a)
beta.0 <- c(rep(1,p.red)/b,rep(1,p.full-p.red))/a 
ns <- round(seq(1e2,5e3,len=30))
by.n <- sapply(ns, function(n) {
    cat('.')
    standardized <- mclapply(1:3e2, mc.cores=detectCores()-3, FUN=function(dd){
        ## cat('.')
        x <- rmvnorm(n,mu,Sigma)
        risk <- h(x%*%beta.0)
        g <- rbinom(n,1,risk)
        x.red <- x[,1:p.red]
        x.0 <- x[g==0,]; x.1<- x[g==1,]    
        x.0.red <- x.red[g==0,]; x.1.red <- x.red[g==1,]    
        beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link='logit')))   
        var.hat <- auc.var(x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red)
        auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red) / sqrt(var.hat)
    })
    standardized <- simplify2array(standardized)
    sd(standardized)
})
plot(ns,by.n)
## save.image('sessions/16f-1.RData')

dd



## 16f-3 check hajek part convergence with adversarial data. also check
## convergence to auc hajek at beta star, besides betahat.
start <- Sys.time()
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
p.full <- 4
p.red <- 2
## n <- 3e2
mu <- rep(0,p.full)
a <- 5
h <- plogis
Sigma <- diag(p.full)*sqrt(a)
beta.0 <- c(rep(0,p.red),rep(1,p.full-p.red))/a 
h <- plogis
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link.name='logit')
beta.star.full <- coef.reduced.glm.gaussian(p.full,params)
## beta.star.red <- coef.reduced.glm.gaussian(p.red=2,link=h,params)
## beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
## beta <- beta.0
lim <- 5
ns <- round(seq(10,3e2,len=20))
auc.full <- auc.glm.gaussian(beta.star.full,params,lim=5)
## by.n <- sapply(ns, function(n) {
by.n <- mclapply(ns, mc.cores=detectCores()-4, FUN=function(n) {
    cat('.')
    replicate(10, {
        tryCatch({
            x <- rmvnorm(n,mu,Sigma)
            risk <- h(x%*%beta.0)
            g <- rbinom(n,1,risk)
            x.0 <- x[g==0,]
            x.1 <- x[g==1,]
            beta.hat <- coef(glm(g~x-1,family=binomial(link=params$link.name)))
            ## f.0 <- function(x)pdf.index.glm.gaussian(x,g=0,beta=beta.hat,params=params,lim=lim) #1 start
            ## f.1 <- function(x)pdf.index.glm.gaussian(x,g=1,beta=beta.hat,params=params,lim=lim) 
            ## F.0 <- Vectorize(function(x)integrate(function(w)f.0(w),-lim,x)$val) 
            ## F.1 <- Vectorize(function(x)integrate(function(w)f.1(w),-lim,x)$val)
            ## auc.hajek.hat <- auc.hajek(x=x.0%*%beta.hat,y=x.1%*%beta.hat,F=F.0,G=F.1,auc=NULL,terms.only=FALSE)
            ## try <- auc.hajek.hat #1 end
            f.0 <- function(x)pdf.index.glm.gaussian(x,g=0,beta=beta.star.full,params=params,lim=lim) #2 start
            f.1 <- function(x)pdf.index.glm.gaussian(x,g=1,beta=beta.star.full,params=params,lim=lim)
            F.0 <- Vectorize(function(x)integrate(function(w)f.0(w),-lim,x)$val)
            F.1 <- Vectorize(function(x)integrate(function(w)f.1(w),-lim,x)$val)
            auc.hajek.star <- auc.hajek(x=x.0%*%beta.star.full,y=x.1%*%beta.star.full,F=F.0,G=F.1,auc=auc.full,terms.only=FALSE)
            try <- auc.hajek.star #2 end
            auc.hat.full.hat <- auc.hat(x=x.0%*%beta.hat,y=x.1%*%beta.hat)
            auc.full.hat <- auc.glm.gaussian(beta.hat,params,lim=5)
            obs <- auc.hat.full.hat - auc.full.hat
            obs - try
        }, error=function(e){print(e); NA})
    })
})
by.n <- simplify2array(by.n)
Sys.time() - start
## mad <- colMeans(abs(by.n))
## plot(ns,mad)
## lm0 <- lm(log(mad)~log(ns))
## coef(lm0)
## curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)
vars <- apply(by.n,2,var,na.rm=TRUE)
plot(ns,vars)
lm0 <- lm(log(vars)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)

sds <- apply(by.n,2,sd,na.rm=TRUE)
plot(ns,sds)
lm0 <- lm(log(sds)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)


## speed improvement for F.0, F.1
source('misc.R')
p <- 4
lim <- 5
beta <- runif(p)
params <- list(mu=runif(p),Sigma=diag(runif(p)),beta=runif(p),link.name='logit',link=plogis)
f.0 <- function(x)pdf.index.logit.gaussian(x,g=0,beta=beta,params=params,lim=lim) #1 start
F.0 <- Vectorize(function(x)integrate(function(w)f.0(w),-lim,x)$val) 
x <- runif(100)*3
start <- Sys.time()
F.0(x)
Sys.time() - start
F.0 <- function(x) {
    old <- rank(x)
    y <- c(-lim,sort(x))
    segments <- sapply(2:length(y), function(i) integrate(function(w)f.0(w),y[i-1],y[i])$val)
    return(cumsum(segments)[old])
}
start <- Sys.time()
F.0(x)
Sys.time() - start

## 16f-4 check hajek part convergence with adversarial data at reduced
## parameter. 1. fine 1/n rate as expected--looking at the
## nonestimated parameters though: in the hajak part, F.0,F.1, and auc
## are computed numerically. beta.hat is estimated of course. also
## needed that lim>=10, with lim=5 doesnt even seem to converge to
## 0. Also when a,b are larger, the rate drops. 1/n rate is at
## a=b=1. at a=b=5, more like n^-.78. at a=b=10 breaks down to
## n^-.15. tried increasing lim to 20, little improvement.

## 2. with estimated parameters (either the F.i or the auc) converges at
## sqrt(n) rate, as expected.

start <- Sys.time()
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(2)
p.full <- 4
p.red <- 2
## n <- 3e2
mu <- rep(0,p.full)
a <- 10
b <- 10
Sigma <- diag(p.full)*sqrt(a)
beta.0 <- c(rep(1,p.red)/b,rep(1,p.full-p.red))/a
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link.name='logit')
lim <- 20
beta.star.red <- coef.reduced.glm.gaussian(p.red,params,lim=lim)
beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
## beta.star.red <- rep(0,p.full)
## params.red <- list(mu=mu[1:p.red],Sigma=Sigma[1:p.red,1:p.red],beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
auc.red <- auc.glm.gaussian(beta.star.red,params)
ns <- round(seq(5e2,3e3,len=20))
## params$auc <- auc.glm.gaussian(beta.0,params,lim=5)
## by.n <- sapply(ns, function(n) {
by.n <- mclapply(ns, mc.cores=detectCores()-4, FUN=function(n) {
    cat('.')
    replicate(30, {
        tryCatch({
            x <- rmvnorm(n,mu,Sigma)
            risk <- h(x%*%beta.0)
            g <- rbinom(n,1,risk)
            ## x.red <- x[,1:p.red]
            x.0 <- x[g==0,]; x.1<- x[g==1,]    
            ## x.0.red <- x.red[g==0,]; x.1.red <- x.red[g==1,]    
            beta.hat.red <- coef(glm(g~x[,1:p.red]-1,family=binomial(link='logit')))
            beta.hat.red <- c(beta.hat.red,rep(0,p.full-p.red))
            ## ## beta.hat.red <- beta.star.red + c(rnorm(p.red)/n^.5,rep(0,p.full-p.red))
            ## ## beta.hat.red <- beta.star.red
            f.0 <- function(x)pdf.index.logit.gaussian(x,g=0,beta=beta.hat.red,params=params,lim=lim) #1 start
            f.1 <- function(x)pdf.index.logit.gaussian(x,g=1,beta=beta.hat.red,params=params,lim=lim) 
            ## F.0 <- Vectorize(function(x)integrate(function(w)f.0(w),-lim,x)$val) 
            ## F.1 <- Vectorize(function(x)integrate(function(w)f.1(w),-lim,x)$val)
            F.0 <- function(x) {
                old <- rank(x)
                y <- c(-lim,sort(x))
                segments <- sapply(2:length(y), function(i) integrate(function(w)f.0(w),y[i-1],y[i])$val)
                return(cumsum(segments)[old])
            }
            F.1 <- function(x) {
                old <- rank(x)
                y <- c(-lim,sort(x))
                 segments <- sapply(2:length(y), function(i) integrate(function(w)f.1(w),y[i-1],y[i])$val)
                return(cumsum(segments)[old])
            }
            auc.hajek.hat <- auc.hajek(x=x.0%*%beta.hat.red,y=x.1%*%beta.hat.red,F=F.0,G=F.1,auc=auc.red,terms.only=FALSE) #1 end
            ## auc.hajek.hat <- 0#auc.hajek(x=x.0%*%beta.hat.red,y=x.1%*%beta.hat.red,terms.only=FALSE)
            try <- auc.hajek.hat 
            auc.hat.red.hat <- auc.hat(x=x.0%*%beta.hat.red,y=x.1%*%beta.hat.red)
            auc.red.hat <- auc.glm.gaussian(beta.hat.red,params,lim=lim)
            obs <- auc.hat.red.hat - auc.red.hat
            obs - try
        }, error=function(e){print(e); NA})
    })
})
by.n <- simplify2array(by.n)
Sys.time() - start
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)
## save.image('sessions/16f-2.RData')


## 16f-5 check infl function


## infl function conv rate with reduced coefs. checks out at O(1/n)
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(1)
mu <- rep(0,p.full)
a <- 10
b <- 10
Sigma <- diag(p.full)*sqrt(a)
beta.0 <- c(rep(1,p.red)/b,rep(1,p.full-p.red))/a
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
## link.name='logit'
## params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link.name='logit')
lim <- 10
params.full <- list(mu=mu,Sigma=Sigma,beta=beta.0,p=p.full,link=h,link.deriv=h.1,link.deriv2=h.2)
beta.star.red <- coef.reduced.glm.gaussian(p.red,params.full,lim=5)
params.red <- list(mu=mu[1:p.red],Sigma=Sigma[1:p.red,1:p.red],beta=beta.star.red,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
## beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
ns <- round(seq(1e2,1e3,len=20))
by.n <- sapply(ns, function(n) {
    replicate(1e2, {
        x <- rmvnorm(n,mu,Sigma)
        eta <- as.numeric(x%*%beta.0)
        risk <- h(eta)
        g <- rbinom(n,1,risk)
        ## beta.hat <- coef(glm(g~x-1,family=binomial(link=link.name)))
        ## obs <- beta.hat-beta.0
        ## approx <- infl.glm.gaussian(x,g,params.full,terms.only=FALSE)
        ## obs-approx
        x.red <- x[,1:p.red]
        beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link=link.name)))
        obs <- beta.hat.red - beta.star.red
        approx <- infl.glm.gaussian(x.red,g,params.red,terms.only=FALSE)
        obs - approx
    })
})
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)


## 16f-6 check taylor expansion

start <- Sys.time()
require(mvtnorm)
require(parallel)
source('misc.R')
set.seed(2)
p.full <- 4
p.red <- 2
## n <- 3e2
mu <- rep(0,p.full)
a <- 5
b <- 5
Sigma <- diag(p.full)*sqrt(a)
beta.0 <- c(rep(1,p.red)/b,rep(1,p.full-p.red))/a
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link.name='logit')
lim <- 20
beta.star.red <- coef.reduced.glm.gaussian(p.red,params,lim=lim)
beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
## beta.star.red <- rep(0,p.full)
## params.red <- list(mu=mu[1:p.red],Sigma=Sigma[1:p.red,1:p.red],beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
auc.red <- auc.logit.gaussian(beta.star.red,params)
deriv <- numDeriv::grad(function(beta)auc.logit.gaussian(beta,params),x=beta.star.red)
ns <- round(seq(2e2,3e3,len=50))
## params$auc <- auc.glm.gaussian(beta.0,params,lim=5)
## by.n <- sapply(ns, function(n) {
by.n <- mclapply(ns, mc.cores=detectCores()-4, FUN=function(n) {
    cat('.')
    replicate(1e2, {
        tryCatch({
            x <- rmvnorm(n,mu,Sigma)
            risk <- h(x%*%beta.0)
            g <- rbinom(n,1,risk)
            ## x.red <- x[,1:p.red]
            x.0 <- x[g==0,]; x.1<- x[g==1,]    
            ## x.0.red <- x.red[g==0,]; x.1.red <- x.red[g==1,]    
            beta.hat.red <- coef(glm(g~x[,1:p.red]-1,family=binomial(link='logit')))
            beta.hat.red <- c(beta.hat.red,rep(0,p.full-p.red))
            ## ## ## beta.hat.red <- beta.star.red + c(rnorm(p.red)/n^.5,rep(0,p.full-p.red))
            ## ## ## beta.hat.red <- beta.star.red
            ## deriv <- 
            try <- deriv%*%(beta.hat.red-beta.star.red) 
            ## auc.hat.red.hat <- auc.hat(x=x.0%*%beta.hat.red,y=x.1%*%beta.hat.red)
            auc.red.hat <- auc.logit.gaussian(beta.hat.red,params,lim=lim)
            obs <-  auc.red.hat - auc.red
            obs - try
        }, error=function(e){print(e); NA})
    })
})
by.n <- simplify2array(by.n)
Sys.time() - start
mad <- colMeans(abs(by.n))
plot(ns,mad)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)
curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)
## save.image('sessions/16f-3.RData') overnight sim, bad rate n^-.4




## 16g non-gaussian covariates

## redo 16f checking z-stats using non-gaussian data
require(mvtnorm)
require(numDeriv)
require(parallel)
source('misc.R')
## set.seed(21)
start <- Sys.time()
p.full <- 4
p.red <- 2
n <- 5e2
mu <- rep(0,p.full)
a <- 1
b <- 1
Sigma <- diag(p.full)*sqrt(a)
beta.0 <- c(0*rep(1,p.red)/b,rep(1,p.full-p.red))/a ##!!!!
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
params.full <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2)
auc.red <- 0
z.stats <- mclapply(1:1e3, mc.cores=detectCores()-3, FUN=function(dd){
    ## cat('.')
    ## x <- rmvnorm(n,mu,Sigma)
    x <- matrix(rexp(n*p.full,rate=1:p.full),ncol=p.full)%*%diag(1:p.full)
    ## x <- matrix(rbeta(n*p.full,4,.4),ncol=p.full)
    ## x <- matrix(rcauchy(n*p.full),ncol=p.full)
    ## x <- matrix(rgamma(n*p.full,1,5),ncol=p.full)
    risk <- h(x%*%beta.0)
    g <- rbinom(n,1,risk)
    x.red <- x[,1:p.red]
    x.0 <- x[g==0,]; x.1<- x[g==1,]    
    x.0.red <- x.red[g==0,]; x.1.red <- x.red[g==1,]    
    beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link=link.name)))
    ## beta.hat.red <- c(beta.hat.red,rep(0,p.full-p.red))
    auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
    obs <- auc.hat.red.hat - auc.red
    auc.hajek.red <- auc.hajek(x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red,F=NULL,G=NULL,auc=auc.hat.red.hat,terms.only=TRUE,IID=TRUE)
    params.hat.red <- list(mu=colMeans(x.red),Sigma=cov(x.red),beta=beta.hat.red,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
    infl.hat <- infl.glm(x.red,g,params.hat.red)
    deriv.red <- grad(function(beta)auc.hat(x.0.red%*%beta,x.1.red%*%beta),beta.hat.red[1:p.red],method='simple',method.args=list(eps=1/n^(.5)))
    ## deriv.red <- rep(0,p.red)
    approx <- auc.hajek.red + as.numeric(deriv.red%*%infl.hat) 
    obs / sqrt(var(approx) / length(approx))
})
Sys.time()-start
z.stats <- simplify2array(z.stats)
var(z.stats,na.rm=TRUE)
hist((z.stats-mean(z.stats,na.rm=TRUE))/sd(z.stats,na.rm=TRUE),prob=TRUE)
curve(dnorm,add=TRUE)


curve(3*x + log(1-plogis(x)),-50,50)

dd
## op <- par(mfrow=c(1,2))
## qqnorm(z.stats)
## abline(0,1)
## hist(z.stats,prob=TRUE)
## curve(dnorm,add=TRUE)
## par(op)
## ## fpr
## alpha <- .05
## print(mean(abs(z.stats)<qnorm(1-alpha/2)))
## ## L_inf distance from standard normal cdf
## z.stats <- sort(z.stats)
## distances <- pnorm(z.stats) - 1:length(z.stats)/length(z.stats)
## max.dist <- max(abs(distances))
## max.dist
## var(z.stats)


## redo 16f checking variance of z-stats using non-gaussian data

## even the hajek/delong term with nonestimated beta seems biased
## above--investigating. [update: it was because i was
## centering with auc=0 rather than the real auc (obtaine
## numerically). i guess even though youre subtracting a constant
## because the denominator of the z-stat is random, it affects the
## variance of the z-stats.

## strange critical fpr behavior at b=9.9973. must be related to numeric routine?
require(mvtnorm)
require(numDeriv)
require(parallel)
source('misc.R')
set.seed(2)
start <- Sys.time()
p.full <- 4
p.red <- 2
## n <- 5e2
mu <- rep(0,p.full)
a <- 5
b <- 9.9973
## a <- b <- 10
Sigma <- diag(p.full)*sqrt(a)
beta.0 <- c(rep(1,p.red)/b,rep(1,p.full-p.red))/a 
## Sigma <- matrix(runif(p.full^2),p.full)
## Sigma <- Sigma%*%t(Sigma)
## beta.0 <- runif(p.full)
## mu <- runif(p.full)
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
link.name='logit'
params.full <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=h,link.deriv=h.1,link.deriv2=h.2)
## auc.red <- 0
beta.star.red <- coef.reduced.glm.gaussian(p.red,params.full,lim=15)
beta.star.red <- c(beta.star.red,rep(0,p.full-p.red))
auc.red <- auc.glm.gaussian(beta.star.red,params.full)
alpha <- .05
ns <- round(seq(1e3,1.5e3,len=2))
## ns <- c(1e3,5e3)
by.n <- sapply(ns, function(n) {
    ## z.stats <- mclapply(1:1e1, mc.cores=detectCores()-4, mc.preschedule=FALSE, FUN=function(dd){
    z.stats <- lapply(1:3e2,  FUN=function(dd){
        ## cat('.')
        x <- rmvnorm(n,mu,Sigma)
        ## x <- matrix(rnorm(n*p.full),ncol=p.full)
        ## x <- pnorm(x)
        ## x <- qexp(x)
        ## x <- matrix(rexp(n*p.full,rate=1:p.full),ncol=p.full)%*%diag(1:p.full)
        ## x <- matrix(rbeta(n*p.full,4,.4),ncol=p.full)
        ## x <- matrix(rcauchy(n*p.full),ncol=p.full)
        ## x <- matrix(rgamma(n*p.full,1,5),ncol=p.full)
        risk <- h(x%*%beta.0)
        g <- rbinom(n,1,risk)
        x.red <- x[,1:p.red]
        x.0 <- x[g==0,]; x.1<- x[g==1,]    
        x.0.red <- x.red[g==0,]; x.1.red <- x.red[g==1,]    
        beta.hat.red <- coef(glm(g~x.red-1,family=binomial(link=link.name)))
        ## beta.hat.red <- c(beta.hat.red,rep(0,p.full-p.red))
        auc.hat.red.hat <- auc.hat(x.0.red%*%beta.hat.red,x.1.red%*%beta.hat.red)
        obs <- auc.hat.red.hat - auc.red
        ## auc.hajek.red <- auc.hajek(x=x.0.red%*%beta.star.red,y=x.1.red%*%beta.star.red,F=NULL,G=NULL,auc=auc.hat.red.hat,terms.only=TRUE,IID=TRUE)
        auc.hajek.red <- auc.hajek(x=x.0.red%*%beta.hat.red,y=x.1.red%*%beta.hat.red,F=NULL,G=NULL,auc=auc.hat.red.hat,terms.only=TRUE,IID=TRUE)
        params.hat.red <- list(mu=colMeans(x.red),Sigma=cov(x.red),beta=beta.hat.red,link=h,link.deriv=h.1,link.deriv2=h.2,p=p.red)
        infl.hat <- infl.glm(x.red,g,params.hat.red)
        pair <- sapply(list(estimated=grad(function(beta)auc.hat(x.0.red%*%beta,x.1.red%*%beta),beta.hat.red[1:p.red],method='simple',method.args=list(eps=1/n^(.5))),zero=rep(0,p.red)), function(deriv.red) {
            approx <- auc.hajek.red + as.numeric(deriv.red%*%infl.hat) 
            obs / sqrt(var(approx) / length(approx))
        })
        ## approx <- auc.hajek.red #+ as.numeric(deriv.red%*%infl.hat) 
        ## (obs / sqrt(var(approx) / length(approx))   )*c(1,1) ##!!!!!!!
    })
    ## Sys.time()-start
    z.stats <- simplify2array(z.stats)
    ## apply(z.stats,1,var,na.rm=TRUE)
    apply(z.stats,1,function(x)mean(abs(x)<qnorm(1-alpha/2)))
})
print(Sys.time()-start)
matplot(ns,t(by.n),pch=1,col=1,type='o',lty=1:2)
legend('bottomleft',lty=1:2,legend=rownames(by.n))
by.n[1,]-by.n[2,]
## save.image('sessions/16f-5.RData') ## adj and unadjusted seem to converge



dd



dd

## z.stats <- simplify2array(z.stats)
## apply(z.stats,1,var,na.rm=TRUE)
## })
## print(Sys.time()-start)
## matplot(ns,t(by.n),pch=1,col=1,type='o',lty=1:2)
## legend('bottomleft',lty=1:2,legend=rownames(by.n))
## ## save.image('sessions/16f-5.RData')


## 16g-1 check derivative performance

## check how gradient estimate converges in cdf case

## n=1. rate seems n^-.25
require(numDeriv)
n <- 1e4
rF <- rnorm
dF <- dnorm
## p <- 2
## mu <- rep(0,p)
## Sigma <- diag(p)
x0 <- runif(1)
ns <- round(seq(1e2,1e4,len=30))
by.n <- sapply(ns, function(n) {
    deriv.hat <- replicate(5e2, {
        x <- rF(n)
        f <- ecdf(x)
        grad(f,x=x0,method='simple',method.args=list(eps=1/n^(.5)))
    })
})
true <- dF(x0)
## hist(deriv.hat)
## abline(v=dF(x0),col=2)
## abline(v=mean(deriv.hat),col=3)
mad <- apply(abs(by.n-true),2,mean)
plot(mad)
## sds <- apply(by.n-true,2,sd)
## plot(sds)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)




## n=2. same settings to grad, rate stays the same at about n^-.25 
require(numDeriv)
require(mvtnorm)
p <- 2
n <- 1e6
mu <- runif(p)
Sigma <- matrix(runif(p^2),p)
Sigma <- Sigma%*%t(Sigma)
pF <- function(q)pmvnorm(lower=-Inf,upper=q,mean=mu,sigma=Sigma)
rF <- function(n)rmvnorm(n,mu,Sigma)
r=c(1,2)
pF <- function(q)pexp(q[1],r[1])*pexp(q[2],r[2])
rF <- function(n)cbind(rexp(n,rate[1]),rexp(n,rate[2]))
pF.hat <- function(q,x)mean(x[,1]<q[1] & x[,2]<q[2])
q0 <- runif(2)
a <- runif(2) # random direction for gradient
true <- as.numeric(a%*%grad(function(q)pF(q),q0))
ns <- round(seq(1e3,2e4,len=30))
by.n <- sapply(ns, function(n) {
    deriv.hat <- replicate(1e2, {
        x <- rF(n)
        ## f(q0,x)
        ## as.numeric(pmvnorm(lower=-Inf,upper=q0,mean=mu,sigma=Sigma))
        grad(function(q)pF.hat(q,x),q0,method='simple',method.args=list(eps=1/n^(.5)))
        ## grad(function(q)pF.hat(q,x),q0)
    })
    a%*%deriv.hat
})
## hist(deriv.hat)
## abline(v=true,col=2)
## abline(v=mean(deriv.hat),col=3)
mad <- apply(abs((by.n-true)),2,mean)
plot(mad)
## sds <- apply(by.n-true,2,sd)
## plot(sds)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)






## curve(h(1,mu.2,sigma.12=x,sigma.2=2),.1,1.8,ylim=c(-1,1))
## x0 <- runif(1,.1,2)
## m <- h.prime(1,mu.2,sigma.12=x0,sigma.2=2)
## abline(b=m,a=h(1,mu.2,sigma.12=x0,sigma.2=2)-m*x0,col=2)


h.cond <- function(x,mu,Sigma,h,lim=Inf) {
    mu.cond <- function(x)mu[2]+Sigma[1,2]/Sigma[1,1]*(x-mu[1])
    sigma.cond <- sqrt((1-Sigma[1,2]^2/prod(diag(Sigma)))*Sigma[2,2])
    sapply(x, function(x)integrate(function(y)h(mu.cond(x)-y)*dnorm(y/sigma.cond)/sigma.cond,-lim,lim)$val)
}
auc <- function(mu,Sigma,h=plogis,lim=Inf) {
    ## browser()
    f.x <- function(x,mu,Sigma)dnorm(x,mu[1],sd=sqrt(Sigma[1,1]))
    f <- function(w,v,mu,Sigma)(1-h.cond(w,mu,Sigma,h,lim))*h.cond(v,mu,Sigma,h,lim)*f.x(w,mu,Sigma)*f.x(v,mu,Sigma)
    ## pi.1 <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    f.y <- function(y,mu,Sigma)dnorm(y,mu[2],sd=sqrt(Sigma[2,2]))
    pi.1 <- integrate(function(x)h(x)*f.y(x,mu,Sigma),-Inf,Inf)$val
    inner <- Vectorize(function(v)integrate(function(w)f(w,v,mu,Sigma),-lim,v)$val)
    integrate(function(v)inner(v),-lim,lim)$val / (pi.1*(1-pi.1))
}
auc.try <- function(mu,Sigma,h=plogis,lim=Inf) {
    ## browser()
    Sigma[1,2] <- Sigma[1,2]/sqrt(Sigma[1,1])
    f.y <- function(y,mu,Sigma)dnorm(y,mu[2],sd=sqrt(Sigma[2,2]))
    pi.1 <- integrate(function(x)h(x)*f.y(x,mu,Sigma),-Inf,Inf)$val
    h.cond <- Vectorize(function(u)integrate(function(xi)h(xi*sqrt(Sigma[2,2]-Sigma[1,2]^2)+mu[2]+Sigma[1,2]*u)*dnorm(xi),-lim,lim)$val)
    ## f.x <- function(x,mu,Sigma)dnorm(x,mu[1],sd=sqrt(Sigma[1,1]))
    f <- function(u,v,mu,Sigma)(1-h.cond(u))*h.cond(v)*dnorm(u)*dnorm(v)
    ## pi.1 <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    inner <- Vectorize(function(v)integrate(function(u)f(u,v,mu,Sigma),-lim,v)$val)
    integrate(function(v)inner(v),-lim,lim)$val / (pi.1*(1-pi.1))
    ## browser()
    ## f.y <- function(y,mu,Sigma)dnorm(y,mu[2],sd=sqrt(Sigma[2,2]))
    ## pi.1 <- integrate(function(x)h(x)*f.y(x,mu,Sigma),-Inf,Inf)$val
    ## ## sigma.cond <- sqrt((1-Sigma[1,2]^2/prod(diag(Sigma)))*Sigma[2,2])
    ## sigma.cond <- sqrt(Sigma[2,2]-Sigma[1,2]^2)
    ## h.cond <- Vectorize(function(u)integrate(function(xi)h(xi)*dnorm((xi-(mu[2]+Sigma[1,2]*u))/sigma.cond)/sigma.cond,-lim,lim)$val,vectorize.args='u')
    ## h.cond <- Vectorize(function(u)integrate(function(xi)h(xi*sigma.cond+mu[2]+Sigma[1,2]*u)*dnorm(xi),-lim,lim)$val)
    ## h.cond <- Vectorize(function(u)integrate(function(xi)h(xi*sqrt(Sigma[2,2]-Sigma[1,2]^2)+mu[2]+Sigma[1,2]*u)*dnorm(xi),-lim,lim)$val)
    ## ## f.x <- function(x,mu,Sigma)dnorm(x,0,sd=sqrt(Sigma[1,1]))
    ## f <- function(w,v,mu,Sigma)(1-h.cond(w))*h.cond(v)*dnorm(w)*dnorm(v)
    ## ## pi.1 <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    ## inner <- Vectorize(function(v)integrate(function(w)f(w,v,mu,Sigma),-lim,v)$val)
    ## integrate(function(v)inner(v),-lim,lim)$val / (pi.1*(1-pi.1))
}
require(mvtnorm)
## source('misc.R')
set.seed(1)
## n <- 1e3
mu <- runif(2)
Sigma <- matrix(runif(4),2)
Sigma <- Sigma%*%t(Sigma)
sigma.12 <- 1
sigma.22 <- 1.3
Sigma <- matrix(c(1,sigma.12,sigma.12,sigma.22^2),2)
## a <- c(0,1)#runif(2)
## t <- -3
h <- pnorm
auc(mu,Sigma,h=h)
auc.try(mu,Sigma,h=h)



## verify auc has a maximum at reduced coefs under probit model
source('misc.R')
auc <- function(mu,Sigma,h=plogis,lim=Inf) {
    ## browser()
    Sigma[1,2] <- Sigma[1,2]/sqrt(Sigma[1,1])
    f.y <- function(y,mu,Sigma)dnorm(y,mu[2],sd=sqrt(Sigma[2,2]))
    ## pi.1 <- integrate(function(x)h(x)*f.y(x,mu,Sigma),-Inf,Inf)$val
    h.cond <- Vectorize(function(u)integrate(function(xi)h(xi*sqrt(Sigma[2,2]-Sigma[1,2]^2)+mu[2]+Sigma[1,2]*u)*dnorm(xi),-lim,lim)$val)
    ## f.x <- function(x,mu,Sigma)dnorm(x,mu[1],sd=sqrt(Sigma[1,1]))
    f <- function(u,v,mu,Sigma)(1-h.cond(u))*h.cond(v)*dnorm(u)*dnorm(v)
    ## pi.1 <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    inner <- Vectorize(function(v)integrate(function(u)f(u,v,mu,Sigma),-lim,v)$val)
    integrate(function(v)inner(v),-lim,lim)$val# / (pi.1*(1-pi.1))#!!!!!!!
}
## set.seed(1)
## pairs <- replicate(5,{
n <- 1e2
h <- pnorm
p.full <- 6
p.red <- 3
Sigma <- matrix(runif(p.full^2),nrow=p.full)
Sigma <- Sigma%*%t(Sigma)
mu <- runif(p.full)
beta.0 <- runif(p.full)
beta.x.reduced <- probit.coef.reduced(mu,Sigma,c(0,beta.0),p.reduced=p.red)[-1]
ts <- seq(-1,1,len=10)/10
a <- c(runif(length(beta.x.reduced)),rep(0,p.full-p.red))
aucs <- sapply(ts, function(t) {
    beta <- c(beta.x.reduced,rep(0,p.full-p.red))
    beta <- beta+a*t
    var.1 <- t(beta)%*%Sigma%*%beta
    var.12 <- t(beta)%*%Sigma%*%beta.0
    var.2 <- t(beta.0)%*%Sigma%*%beta.0
    Sigma.simple <-  matrix(c(var.1,var.12,var.12,var.2),2)
    ## params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=pnorm)
    auc(mu,Sigma.simple,h=h)
})
## auc.glm.gaussian(beta,params)
## )
## })
## plot(pairs[1,],pairs[2,]); abline(0,1)
## (t(beta)%*%Sigma%*%beta)*beta.0 - (t(beta)%*%Sigma%*%beta.0)*beta
## ts <- seq(-.001,.001,len=4)
## aucs <- sapply(ts,function(t){
##     print(t)
##     var.12 <- t(beta)%*%Sigma%*%beta.0 + t
##     Sigma.simple <-  matrix(c(var.1,var.12,var.12,var.2),2)
##     auc(mu,Sigma.simple,h=h)
## })
## plot(ts,aucs)
plot(ts,aucs,type='l')
abline(v=0)


## parameterizing the bivariate normal directly, no betas/quad forms
## involved
## 1. realized not the d/d(sigma_{12}(AUC) that is 0 but d/d(beta)(sigma_{12})
h <- pnorm
var.1 <- 1
var.2 <- 2
mu <- runif(2)
ts <- seq(0,var.2-var.1-.1,len=5)/1
aucs <- sapply(ts, function(t) {
    var.1.t <- var.1+t
    var.12.t <- var.1.t*sqrt((var.2+1)/(var.1.t+1))
    Sigma <-  matrix(c(var.1.t,var.12.t,var.12.t,var.2),2)
    ## print(det(Sigma.simple))
    ## params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=pnorm)
    auc(mu,Sigma,h=h) 
})
plot(ts,aucs,type='l')

auc <- function(mu,Sigma,h=plogis,lim=Inf) {
    ## browser()
    Sigma[1,2] <- Sigma[1,2]/sqrt(Sigma[1,1])
    f.y <- function(y,mu,Sigma)dnorm(y,mu[2],sd=sqrt(Sigma[2,2]))
    ## pi.1 <- integrate(function(x)h(x)*f.y(x,mu,Sigma),-Inf,Inf)$val
    h.cond <- Vectorize(function(u)integrate(function(xi)h(xi*sqrt(Sigma[2,2]-Sigma[1,2]^2)+mu[2]+Sigma[1,2]*u)*dnorm(xi),-lim,lim)$val)
    ## f.x <- function(x,mu,Sigma)dnorm(x,mu[1],sd=sqrt(Sigma[1,1]))
    f <- function(u,v)(1-h.cond(u))*h.cond(v)*dnorm(u)*dnorm(v)
    ## pi.1 <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    inner <- Vectorize(function(v)integrate(function(u)f(u,v),-lim,v)$val)
    integrate(function(v)inner(v),-lim,lim)$val# / (pi.1*(1-pi.1))#!!!!!!!
}
## set.seed(1)
## pairs <- replicate(5,{
n <- 1e2
h <- pnorm
p.full <- 6
p.red <- 3
Sigma <- matrix(runif(p.full^2),nrow=p.full)
Sigma <- Sigma%*%t(Sigma)
mu <- runif(p.full)
beta.0 <- runif(p.full)
beta.x.reduced <- probit.coef.reduced(mu,Sigma,c(0,beta.0),p.reduced=p.red)[-1]
ts <- seq(-1,1,len=10)/10
a <- c(runif(length(beta.x.reduced)),rep(0,p.full-p.red))
aucs <- sapply(ts, function(t) {
    beta <- c(beta.x.reduced,rep(0,p.full-p.red))
    beta <- beta+a*t
    var.1 <- t(beta)%*%Sigma%*%beta
    var.12 <- t(beta)%*%Sigma%*%beta.0
    var.2 <- t(beta.0)%*%Sigma%*%beta.0
    Sigma.simple <-  matrix(c(1,var.12/sqrt(var.1),var.12/sqrt(var.1),var.2),2)
    params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=pnorm)
    ## auc(mu,Sigma.simple,h=h)
    var.12/sqrt(var.1) ## 1
})
plot(ts,aucs,type='l')
abline(v=0)

numDeriv::grad(function(beta) {
    beta <- c(beta,rep(0,p.full-p.red))
    var.1 <- t(beta)%*%Sigma%*%beta
    var.12 <- t(beta)%*%Sigma%*%beta.0
    var.12/sqrt(var.1)
}, x=beta.x.reduced)

dd


## glm derivative with h=probit. should vanish.
## set.seed(1)
p.x <- 3
p.y <- 3
Sigma <- matrix(runif((p.x+p.y)^2),nrow=p.x+p.y)
Sigma <- Sigma%*%t(Sigma)
## Sigma.xx <- Sigma[1:p.x,1:p.x]
## Sigma.yy <- Sigma[(p.x+1):(p.x+p.y),(p.x+1):(p.x+p.y)]
## Sigma.xy <- Sigma[1:p.x,(p.x+1):(p.x+p.y)]
## Sigma.yx <- t(Sigma.xy)
mu.x <- runif(p.x); mu.y <- runif(p.y)
mu <- c(mu.x,mu.y)
beta.x <- runif(p.x)/10; beta.y <- runif(p.y)/10
beta <- c(beta.x,beta.y)

source('misc.R')
set.seed(1)
p.full <- 6
p.red <- 3
Sigma <- matrix(runif(p.full^2),nrow=p.full)
Sigma <- Sigma%*%t(Sigma)
mu <- runif(p.full)
beta.0 <- runif(p.full)
beta.x.reduced <- probit.coef.reduced(mu,Sigma,c(0,beta.0),p.reduced=p.red)[-1]
## ts <- seq(-1,1,len=10)/10
## a <- c(runif(length(beta.x.reduced)),rep(0,p.full-p.red))
beta <- c(beta.x.reduced,rep(0,p.full-p.red))
mu.2 <- mu%*%beta.0
sigma.12 <- (t(beta)%*%Sigma%*%beta.0) / sqrt(t(beta)%*%Sigma%*%beta)
sigma.2 <- sqrt( (t(beta.0)%*%Sigma%*%beta.0))
h0 <- pnorm
h0.prime <- dnorm
## mu.2 <- runif(1)
## sigma.12 <- runif(1)
## sigma.2 <- runif(1,sigma.12,2)
h <- Vectorize(function(u,mu.2,sigma.12,sigma.2) integrate(function(xi)h0(xi*sqrt(sigma.2^2-sigma.12^2)+mu.2+sigma.12*u)*dnorm(xi),-Inf,Inf)$val,vectorize.args=c('u','sigma.12'))
h.prime <- Vectorize(function(u,mu.2,sigma.12,sigma.2) integrate(function(xi)h0.prime(xi*sqrt(sigma.2^2-sigma.12^2)+mu.2+sigma.12*u)*(u-sigma.12/sqrt(sigma.2^2-sigma.12^2)*xi)*dnorm(xi),-Inf,Inf)$val,vectorize.args=c('u','sigma.12'))
inner <- Vectorize(function(v)integrate(function(u)dnorm(u)*dnorm(v)*(-h.prime(u,mu.2,sigma.12,sigma.2)*h(v,mu.2,sigma.12,sigma.2)+(1-h(u,mu.2,sigma.12,sigma.2))*h.prime(v,mu.2,sigma.12,sigma.2)),-Inf,v)$val)
integrate(function(v)inner(v),-Inf,Inf)$val

## omits pi.1(1-pi.1) denominator
auc <- function(mu,Sigma,h0=pnorm,lim=Inf) {
    ## browser()
    Sigma[1,2] <- Sigma[1,2]/sqrt(Sigma[1,1])
    ## f.y <- function(y,mu,Sigma)dnorm(y,mu[2],sd=sqrt(Sigma[2,2]))
    ## pi.1 <- integrate(function(x)h0(x)*f.y(x,mu,Sigma),-Inf,Inf)$val
    h.cond <- Vectorize(function(u)integrate(function(xi)h0(xi*sqrt(Sigma[2,2]-Sigma[1,2]^2)+mu[2]+Sigma[1,2]*u)*dnorm(xi),-lim,lim)$val)
    ## f.x <- function(x,mu,Sigma)dnorm(x,mu[1],sd=sqrt(Sigma[1,1]))
    f <- function(u,v,mu,Sigma)(1-h.cond(u))*h.cond(v)*dnorm(u)*dnorm(v)
    ## pi.1 <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    inner <- Vectorize(function(v)integrate(function(u)f(u,v,mu,Sigma),-lim,v)$val)
    integrate(function(v)inner(v),-lim,lim)$val #/ (pi.1*(1-pi.1))
}
auc.prime <- function(mu,Sigma,h0=pnorm,h0.prime=dnorm) {
    mu.2 <- mu[2]
    sigma.12 <- Sigma[1,2]; sigma.2 <- sqrt(Sigma[2,2])
    ## h <- Vectorize(function(u) integrate(function(xi)h0(xi*sqrt(sigma.2^2-sigma.12^2)+mu.2+sigma.12*u)*dnorm(xi),-Inf,Inf)$val)
    h <- Vectorize(function(u) integrate(function(xi)h0(xi*sqrt(Sigma[2,2]-Sigma[1,2]^2)+mu[2]+Sigma[1,2]*u)*dnorm(xi),-Inf,Inf)$val)
    h.prime <- Vectorize(function(u,mu.2,sigma.12,sigma.2) integrate(function(xi)h0.prime(xi*sqrt(sigma.2^2-sigma.12^2)+mu.2+sigma.12*u)*(u-sigma.12/sqrt(sigma.2^2-sigma.12^2)*xi)*dnorm(xi),-Inf,Inf)$val,vectorize.args=c('u','sigma.12'))
    inner <- Vectorize(function(v)integrate(function(u)dnorm(u)*dnorm(v)*(-h.prime(u,mu.2,sigma.12,sigma.2)*h(v)+(1-h(u))*h.prime(v,mu.2,sigma.12,sigma.2)),-Inf,v)$val)
    integrate(function(v)inner(v),-Inf,Inf)$val
}
lim <- Inf
mu <- runif(2)
sigma.2 <- 2
sigma.12 <- 1
Sigma <- matrix(c(1,sigma.12,sigma.12,sigma.2^2),2)
ts <- seq(0,sqrt(Sigma[2,2])-Sigma[1,2]-.1,len=10)
aucs <- sapply(ts, function(t) {
    Sigma.t <- Sigma; Sigma.t[1,2] <- Sigma.t[2,1] <- Sigma[1,2]+t
    ## print(Sigma.t[2,2]-Sigma.t[1,2]^2)
    auc(mu,Sigma.t,h=pnorm)
    })
deriv <- auc.prime(mu,Sigma)
plot(ts,aucs)
abline(b=deriv,a=aucs[1],col=2)


n <- 1e2
h <- pnorm
p.full <- 6
p.red <- 3
Sigma <- matrix(runif(p.full^2),nrow=p.full)
Sigma <- Sigma%*%t(Sigma)
mu <- runif(p.full)
beta.0 <- runif(p.full)
beta.x.reduced <- probit.coef.reduced(mu,Sigma,c(0,beta.0),p.reduced=p.red)[-1]
ts <- seq(-1,1,len=5)/10
a <- c(runif(length(beta.x.reduced)),rep(0,p.full-p.red))
aucs <- sapply(ts, function(t) {
    beta <- c(beta.x.reduced,rep(0,p.full-p.red))
    beta <- beta+a*t
    var.1 <- t(beta)%*%Sigma%*%beta
    var.12 <- t(beta)%*%Sigma%*%beta.0
    var.2 <- t(beta.0)%*%Sigma%*%beta.0
    Sigma.simple <-  matrix(c(var.1,var.12,var.12,var.2),2)
    ## params <- list(mu=mu,Sigma=Sigma,beta=beta.0,link=pnorm)
    auc(c(.3,.4),Sigma.simple,h=h)
})
plot(ts,aucs,type='l')


Sigma.xx <- Sigma[1:p.red,1:p.red]
Sigma.yy <- Sigma[(p.red+1):(p.full),(p.red+1):(p.full)]
Sigma.xy <- Sigma[1:p.red,(p.red+1):(p.full)]
Sigma.yx <- t(Sigma.xy)
beta.x <- beta.0[1:p.red]; beta.y <- beta.0[(p.red+1):p.full]
denom <- sqrt(1+t(beta.y)%*%(Sigma.yy-Sigma.yx%*%solve(Sigma.xx)%*%Sigma.xy)%*%beta.y)
## beta.x.reduced / (beta.x+solve(Sigma.xx)%*%Sigma.xy%*%beta.y)
t(beta.x)%*%Sigma.xx%*%beta.x + t(beta.y)%*%Sigma.yx%*%solve(Sigma.xx)%*%Sigma.xy%*%beta.y+2*t(beta.x)%*%Sigma.xy%*%beta.y
 

mu <- c(mu%*%beta,mu%*%beta.0)
beta <- c(beta.x.reduced,rep(0,p.full-p.red))
var.1 <- t(beta)%*%Sigma%*%beta
var.12 <- t(beta)%*%Sigma%*%beta.0
var.2 <- t(beta.0)%*%Sigma%*%beta.0
Sigma.simple <-  matrix(c(1,var.12/sqrt(var.1),var.12/sqrt(var.1),var.2),2)
deriv <- auc.prime(mu,Sigma.simple)

abline(a=aucs[1],b=deriv,col=2)




## 17 data analysis

fhs <- read.csv('framingham.csv')
## fhs <- within(fhs, cigsPerDay[cigsPerDay==0] <- .1) # bump for log transform
fhs$smoker <- fhs$cigsPerDay>0
for(col in c('age','cigsPerDay','totChol','sysBP','diaBP','BMI','heartRate','glucose')) fhs[,col] <- log(fhs[,col])
x <- subset(fhs,select=c(TenYearCHD,age,totChol,sysBP,diaBP,diabetes,male,smoker))
x <- na.omit(x)
d <- x$TenYearCHD
x <- subset(x,select=-TenYearCHD)
x <- as.matrix(cbind(x))
glm0 <- glm(d~x-1,family=binomial('logit'))
summary(glm0)
source('misc.R')
auc.index.linearize <- function(x,d,beta,infl.fn,deriv.fn=NULL) {
    ## browser()
    n <- nrow(x)
    x.0 <- x[d==0,]; x.1 <- x[d==1,] # clean up
    if(is.null(deriv.fn)) deriv.fn <- function(x,d,beta)numDeriv::grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta,method='simple',method.args=list(eps=1/n^(.5)))
    approx.1 <- auc.hajek(x.0%*%beta,x.1%*%beta,terms.only=TRUE,IID=TRUE)
    infl.hat <- infl.fn(x,d)
    deriv.hat <- deriv.fn(x,d,beta)#grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat)
    ## print(deriv.hat)
    ## print(head(infl.hat))
    approx.2 <- deriv.hat%*%infl.hat
    ## obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc
    approx <- as.numeric(approx.1+approx.2)
    ## var(approx) / length(approx)
}
## params.full$x <- xd$x; params.full$d <- xd$d
## params.red$x <- xd$x[,1:p.red]; params.red$d <- xd$d
x.full <- x
x.red <- subset(x,select=-diabetes)
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
h.2 <- function(x)h.1(x)*(1-2*h(x))
out <- lapply(list(full=x.full,red=x.red), function(x) {
    ## with(params, {
    ## browser()
    x.0 <- x[d==0,]; x.1 <- x[d==1,]
    ## beta.hat <- coefs.lda(x.0,x.1)
    beta.hat <- coef(glm(d~x-1,family=binomial('logit')))
    obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat)
    params <- list(p=ncol(x),beta=beta.hat,link=h,link.deriv=h.1,link.deriv2=h.2)
    linearized <- lapply(list(delong= function(x,d,beta)rep(0,ncol(x)), proposed=NULL), function(deriv.fn)
        auc.index.linearize(x,d,beta.hat,infl.fn=function(x,d)infl.glm(x,d,params,terms.only=TRUE),deriv.fn=deriv.fn))
    linearized <- simplify2array(linearized)
    list(obs=obs,linearized=linearized)
})
diff.linearized <- out$full$linearized - out$red$linearized
diff.obs <- out$full$obs - out$red$obs
## diff.linearized <-  out$red$linearized
## diff.obs <-  out$red$obs
diff.var.hat <- apply(diff.linearized,2,function(iid)var(iid) / length(iid))
c(diff.obs) / sqrt(diff.var.hat)
