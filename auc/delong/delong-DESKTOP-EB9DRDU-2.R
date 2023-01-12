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
## save.image('220930.RData')
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
ns <- round(seq(1e2,3e4,len=2e1))
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
## save.image('221001.RData')
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
ns <- round(seq(1e2,1e3,len=20))
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
        Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
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
        observed <- auc.hat.beta.hat - auc.hat.beta - (auc.hat.gamma.hat - auc.hat.beta)
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
deriv.2 <- -t(delta)%*%auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)%*%delta
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
