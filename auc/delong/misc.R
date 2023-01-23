delong.test <- function(x=NULL,y=NULL,xy=NULL,g=NULL) {
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
    z.stat <- diff(rev(theta.hats)) / sqrt(t(contrast)%*%S%*%contrast)
    return(z.stat)
}

delong.var <- function(x=NULL,y=NULL,xy=NULL,g=NULL) {
    ## browser()
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
    hajek.diff.terms <- mapply('-',auc.hajek(F.hat.1,G.hat.1,x[,1],y[,1],0),auc.hajek(F.hat.2,G.hat.2,x[,2],y[,2],0),SIMPLIFY=FALSE)
    sum(sapply(hajek.diff.terms,var) / c(m,n))
}

## assumes x.0 and x.1 have intercept terms in first column
lda.coefs <- function(x.0,x.1) {
    x.0 <- x.0[,-1]; x.1 <- x.1[,-1]
    n.0 <- nrow(x.0); n.1 <- nrow(x.1); n <- n.0+n.1
    pi.0.hat <- n.0/n; pi.1.hat <- n.1/n
    mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
    Sigma.hat <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n
    Sigma.hat.inv <- solve(Sigma.hat)
    beta.hat <- c(log(pi.1.hat/pi.0.hat) - 1/2*(t(mu.1.hat)%*%Sigma.hat.inv%*%mu.1.hat - t(mu.0.hat)%*%Sigma.hat.inv%*%mu.0.hat), t(mu.1.hat-mu.0.hat)%*%Sigma.hat.inv)
    return(beta.hat)
}


## coefs.lda <- function(x.0,x.1,Sigma=NULL) {
##     ## browser()
##     x.0 <- x.0[,-1]; x.1 <- x.1[,-1]
##     n.0 <- nrow(x.0); n.1 <- nrow(x.1); n <- n.0+n.1
##     pi.0.hat <- n.0/n; pi.1.hat <- n.1/n
##     mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
##     if(is.null(Sigma)) Sigma <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n
##     Sigma.inv <- solve(Sigma)
##     beta.hat <- c(log(pi.1.hat/pi.0.hat) - 1/2*(t(mu.1.hat)%*%Sigma.inv%*%mu.1.hat - t(mu.0.hat)%*%Sigma.inv%*%mu.0.hat), t(mu.1.hat-mu.0.hat)%*%Sigma.inv)
##     return(beta.hat)
## }

## ## ignore intercept term
## coefs.lda <- function(x.0,x.1,params) {
##     ## browser()
##     ## x.0 <- x.0[,-1]; x.1 <- x.1[,-1]
##     n.0 <- nrow(x.0); n.1 <- nrow(x.1); n <- n.0+n.1
##     ## pi.0.hat <- n.0/n; pi.1.hat <- n.1/n
##     mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
##     Sigma <- params$Sigma
##     if(is.null(Sigma)) Sigma <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n
##     Sigma.inv <- solve(Sigma)
##     beta.hat <- Sigma.inv%*%(mu.1.hat-mu.0.hat)
##     return(beta.hat)
## }

## add support for heteroskedasticity
coefs.lda <- function(x.0,x.1,params=NULL) {
    ## x.0 <- x.0[,-1]; x.1 <- x.1[,-1]
    ## browser()
    n.0 <- nrow(x.0); n.1 <- nrow(x.1); n <- n.0+n.1
    mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
    Sigma <- params$Sigma
    if(is.null(Sigma)) Sigma <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n
    pi.0 <- params$pi.0 ## for heteroskedasticity
    ## browser()
    Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
    if(is.null(pi.0))pi.0 <- 1/2 
    if(is.null(Sigma.0))Sigma.0 <- Sigma
    if(is.null(Sigma.1))Sigma.1 <- Sigma
    Sigma <- pi.0*Sigma.0 + (1-pi.0)*Sigma.1
    Sigma.inv <- solve(Sigma)
    ## if(max(abs(Sigma-ss))>.003)browser()
    beta.hat <- Sigma.inv%*%(mu.1.hat-mu.0.hat)
    return(beta.hat)
}

## TODO case with ties
auc.hat <- function(x,y,ties=FALSE)if(!ties)mean(outer(x,y,'<')) else NA
## deprecate
auc.hat.continuous <- function(x,y)mean(outer(x,y,'<'))


lda.sampler.init <- function(beta,Sigma) {
    p <- length(beta)-1
    stopifnot(nrow(Sigma)==p)
    Sigma.inv <- solve(Sigma)
    mu.0 <- rep(0,p)
    mu.1 <- mu.0+Sigma%*%beta[-1]
    pi.1 <- plogis(beta[1]+1/2*(t(mu.1)%*%Sigma.inv%*%mu.1 - t(mu.0)%*%Sigma.inv%*%mu.0))
    pi.0 <- 1-pi.1
    return(function(n) {
        n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
        stopifnot(n.0 > 0 && n.1 > 0)
        x.0 <- mvtnorm::rmvnorm(n.0,mean=mu.0,sigma=Sigma)
        x.1 <- mvtnorm::rmvnorm(n-n.0,mean=mu.1,sigma=Sigma)
        return(list(n.0=n.0,n.1=n.1,x.0=x.0,x.1=x.1))
        })    
}


lda.sampler.init <- function(beta,Sigma) {
    p <- length(beta)-1
    stopifnot(nrow(Sigma)==p)
    Sigma.inv <- solve(Sigma)
    mu.0 <- rep(0,p)
    mu.1 <- mu.0+Sigma%*%beta[-1]
    pi.1 <- plogis(beta[1]+1/2*(t(mu.1)%*%Sigma.inv%*%mu.1 - t(mu.0)%*%Sigma.inv%*%mu.0))
    pi.0 <- 1-pi.1
    return(list(
        mu.0=mu.0,
        mu.1=mu.1,
        pi.0=pi.0,
        pi.1=pi.1,
        sample=function(n) {
            n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
            stopifnot(n.0 > 0 && n.1 > 0)
            x.0 <- mvtnorm::rmvnorm(n.0,mean=mu.0,sigma=Sigma)
            x.1 <- mvtnorm::rmvnorm(n.1,mean=mu.1,sigma=Sigma)
            return(list(n.0=n.0,n.1=n.1,x.0=cbind(1,x.0),x.1=cbind(1,x.1)))
        }
        ))
    }

lda.sampler.init.std <- function(Delta,lambda,p) {
    ## p <- length(beta)-1
    ## stopifnot(nrow(Sigma)==p)
    ## Sigma.inv <- solve(Sigma)
    ## mu.1 <- rep(0,p)
    ## mu.0 <- mu.0+Sigma%*%beta[-1]
    ## pi.1 <- plogis(beta[1]+1/2*(t(mu.1)%*%Sigma.inv%*%mu.1 - t(mu.0)%*%Sigma.inv%*%mu.0))
    ## pi.0 <- 1-pi.1
    pi.1 <- plogis(lambda)
    pi.0 <- 1-pi.1
    e.1 <- c(1,rep(0,p-1))
    mu.1 <- Delta/2 * e.1
    mu.0 <- -mu.1
    Sigma <- Sigma.inv <- diag(p)
    beta <- c(lambda,Delta*e.1)
    return(list(
        mu.0=mu.0,
        mu.1=mu.1,
        pi.0=pi.0,
        pi.1=pi.1,
        beta=beta,
        Sigma=Sigma,
        mu.diff=c(0,mu.1-mu.0),
        Sigma.diff= 2*cbind(0,rbind(0,Sigma)),
        sample=function(n) {
            n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
            stopifnot(n.0 > 0 && n.1 > 0)
            x.0 <- matrix(rnorm(n.0*p),ncol=p)
            x.1 <- matrix(rnorm(n.1*p),ncol=p)
            x.1[,1] <- x.1[,1] + Delta/2
            x.0[,1] <- x.0[,1] - Delta/2
            ## x.0 <- mvtnorm::rmvnorm(n.0,mean=mu.0,sigma=Sigma)
            ## x.1 <- mvtnorm::rmvnorm(n-n.0,mean=mu.1,sigma=Sigma)
            return(list(n.0=n.0,n.1=n.1,x.0=cbind(1,x.0),x.1=cbind(1,x.1)))
        }
        ))
    }
## auc.scores <- function(coefs,mu.diff,Sigma.diff) {
##     ## mu.diff <- c(0,mu.diff)
##     ## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
##     pnorm(-(coefs%*%mu.diff)/sqrt(t(coefs)%*%(Sigma.diff)%*%coefs))
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
## ## mu.diff <- c(0,mu.0-mu.1)
## ## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
## abline(v=auc.scores(beta,mu.diff=mu.diff,Sigma.diff),col=2)
## abline(v=mean(aucs),col=3)



## calling "auc.scores" by "auc.lda" now...deprecate

## P(coefs%*%x.0 < coefs%*%x.1) where x.0-x.1 is multivariate normal
## with mean mu.diff and variance Sigma.diff
auc.scores <- function(coefs,mu.diff,Sigma.diff) {
    ## mu.diff <- c(0,mu.diff)
    ## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
    pnorm((coefs%*%mu.diff)/sqrt(t(coefs)%*%(Sigma.diff)%*%coefs))
}
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
## mu.diff <- c(0,lda.sampler$mu.0-lda.sampler$mu.1)
## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
## abline(v=auc.scores(beta,mu.diff=mu.diff,Sigma.diff),col=2)
## abline(v=mean(aucs),col=3)

## derivative of auc.scores
auc.scores.deriv <-  function(coefs,mu.diff,Sigma.diff) {
    quad <- as.numeric(coefs%*%Sigma.diff%*%coefs)
    -(Sigma.diff%*%coefs %*% (coefs%*%mu.diff/quad) - mu.diff)%*%(1/sqrt(quad)) %*% dnorm(coefs%*%mu.diff/sqrt(quad))
}
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


## second derivative of auc(beta), but only valid when evaluated at beta==the lda coefs
auc.scores.deriv2.lda <- function(coefs,Sigma.diff) {
      quad <- as.numeric(coefs%*%Sigma.diff%*%coefs)
      -dnorm(sqrt(quad)/2)/sqrt(quad)*(Sigma.diff - (Sigma.diff%*%coefs)%*%t(Sigma.diff%*%coefs)/quad)/2
      }
## ## source('misc.R')
## p <- 3
## n <- 1e3
## beta <- runif(p+1)
## Sigma <- matrix(rnorm(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
## lda.sampler <- lda.sampler.init(beta,Sigma)
## delta <- runif(p+1)
## a <- runif(p+1)
## ts <- seq(0,1,len=2000)
## mu.diff <- c(0,lda.sampler$mu.1-lda.sampler$mu.0)
## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
## ## means <- sapply(ts, function(t)auc.scores(coefs=beta+t*delta,mu.diff=mu.diff,Sigma.diff=Sigma.diff))
## ## plot(ts,means,type='l')
## derivs <- sapply(ts, function(t)a%*%auc.scores.deriv(coefs=beta+t*delta,mu.diff=mu.diff,Sigma.diff=Sigma.diff))
## plot(ts,derivs,type='l')
## ## curve(means[1]+x*deriv.1,add=TRUE,col=2)      
## ## abline(a=means[1],b=slope,col=2)
## deriv.2 <- t(delta)%*%auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)%*%a
## curve(derivs[1]+x*deriv.2,add=TRUE,col=2)      

## p <- 3
## n <- 1e3
## beta <- runif(p+1)
## Sigma <- matrix(rnorm(p^2),nrow=p)
## Sigma <- Sigma%*%t(Sigma)
## lda.sampler <- lda.sampler.init(beta,Sigma)
## delta <- runif(p+1)
## ## a <- runif(p+1)
## ts <- seq(0,1,len=2000)
## mu.diff <- c(0,lda.sampler$mu.1-lda.sampler$mu.0)
## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
## means <- sapply(ts, function(t)auc.scores(coefs=beta+t*delta,mu.diff=mu.diff,Sigma.diff=Sigma.diff))
## plot(ts,means,type='l')
## ## derivs <- sapply(ts, function(t)a%*%auc.scores.deriv(coefs=beta+t*delta,mu.diff=mu.diff,Sigma.diff=Sigma.diff))
## ## plot(ts,derivs,type='l')
## ## curve(means[1]+x*deriv.1,add=TRUE,col=2)      
## ## abline(a=means[1],b=slope,col=2)
## deriv.1 <- t(delta)%*%auc.scores.deriv(coefs=beta,mu.diff=mu.diff,Sigma.diff=Sigma.diff)
## curve(means[1]+deriv.1*x,add=TRUE,col=2)
## deriv.2 <- t(delta)%*%auc.scores.deriv2.lda(coefs=beta,Sigma.diff=Sigma.diff)%*%delta
## ## curve(derivs[1]+x*deriv.2,add=TRUE,col=2)      
## curve(means[1]+deriv.1*x+deriv.2*x^2/2,add=TRUE,col=3)


## auc.hajek <- function(F,G,x,y,auc,terms.only=TRUE) {
##     ## -G(x) + F(y) + 1-2*auc
##     terms <- list(control = -(G(x)-(1-auc)), case = F(y)-auc)
##     if(terms.only) return(terms) else return(sum(sapply(terms,mean)))
## }

## auc.hajek <- function(F=NULL,G=NULL,x,y,auc,terms.only=TRUE,IID=FALSE) {
##     ## -G(x) + F(y) + 1-2*auc
##     terms <- list(control = -(G(x)-(1-auc)), case = F(y)-auc)
##     if(terms.only) return(terms)
##     if(IID) {
##         ## stopifnot(terms.only)
##         n.0 <- length(x); n.1 <- length(y)
##         n <- n.0+n.1
##         ## g <- c(rep(0,n.0),rep(1,n.1))
##         return(c(terms$control*n/n.0,terms$case*n/n.1))
##     }
##     return(sum(sapply(terms,mean)))
## }

## support for empirically estimated parameters
auc.hajek <- function(x,y,F=NULL,G=NULL,auc=NULL,terms.only=TRUE,IID=FALSE) {
    ## -G(x) + F(y) + 1-2*auc
    ## browser()
    if(is.null(F)) F <- ecdf(x)
    if(is.null(G)) G <- ecdf(y)
    if(is.null(auc)) auc <- auc.hat(x,y)
    terms <- list(control = -(G(x)-(1-auc)), case = F(y)-auc)
    if(terms.only) {
        if(IID) {
            ## stopifnot(terms.only)
            n.0 <- length(x); n.1 <- length(y)
            n <- n.0+n.1
            ## g <- c(rep(0,n.0),rep(1,n.1))
            return(c(terms$control*n/n.0,terms$case*n/n.1))
        } else  return(terms)
    }
    return(sum(sapply(terms,mean)))
}


## auc.hajek.lda <- function(x.0,x.1,mu.0,mu.1,beta,Sigma) {
##     quad <- as.numeric(beta%*%Sigma%*%beta)
##     F <- function(u)pnorm((u-as.numeric(beta%*%mu.0))/sqrt(quad)) #cdf of beta^Tx.0
##     G <- function(u)pnorm((u-as.numeric(beta%*%mu.1))/sqrt(quad)) #cdf of beta^Tx.1
##     theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
##     auc.hajek(F,G,x.0%*%beta,x.1%*%beta,auc=theta,FALSE)
## } 


auc.hajek.lda <- function(x.0,x.1,beta,params,terms.only=TRUE,IID=terms.only){
    beta.index <- beta
    with(params, {
        ## browser()
        quad <- as.numeric(t(beta)%*%Sigma%*%beta)
        F <- function(u)pnorm((u-as.numeric(t(beta)%*%mu.0))/sqrt(quad)) #cdf of beta^Tx.0
        G <- function(u)pnorm((u-as.numeric(t(beta)%*%mu.1))/sqrt(quad)) #cdf of beta^Tx.1
        theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
        auc.hajek(F,G,x.0%*%beta.index,x.1%*%beta.index,auc=theta,terms.only,IID)
    })
} 

## auc.scores.hajek <- function(x.0,x.1,mu.0,mu.1,beta,Sigma) {
##     x.0 <- t(x.0); x.1 <- t(x.1) 
##     quad <- as.numeric(beta%*%Sigma%*%beta)
##     theta <- as.numeric(pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*quad)))
##     hajek <- -mean(pnorm(t(beta)%*%(x.0-mu.1) / sqrt(quad)) - (1-theta)) + mean(pnorm(t(beta)%*%(x.1-mu.0) / sqrt(quad)) - theta)
##     return(hajek)
## }

auc.scores.hajek <- auc.hajek.lda # deprecate

## ## test for auc.hajek.gaussian, the rewritten auc.scores.hajek (#10g)
## cmp <- replicate(1e2, {
##     p <- 3
##     n <- 10
##     beta <- runif(p)
##     mu.0 <- runif(p)
##     mu.1 <- runif(p)
##     Sigma <- matrix(runif(p^2),p)
##     Sigma <- Sigma%*%Sigma
##     x.0 <- matrix(runif(n*p),ncol=p)
##     x.1 <- matrix(runif(n*p),ncol=p)
##     auc.scores.hajek(x.0,x.1,mu.0,mu.1,beta,Sigma) - auc.hajek.gaussian(x.0,x.1,mu.0,mu.1,beta,Sigma)
## })
## sum(abs(cmp))



## test for auc.scores.hajek
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

## unit test at 12g
auc.hajek.probit <- function(x.0,x.1,mu,Sigma,beta,bw=7,terms.only=FALSE,IID=FALSE) {    
    mu <- t(beta)%*%mu
    sigma <- sqrt(t(beta)%*%Sigma%*%beta)
    f.0 <- function(w)(1-pnorm(w))*dnorm((w-mu)/sigma)/sigma/(1-pnorm(mu/sqrt(1+sigma^2)))
    f.1 <- function(w)pnorm(w)*dnorm((w-mu)/sigma)/sigma/pnorm(mu/sqrt(1+sigma^2))
    F.0 <- Vectorize(function(w)integrate(f.0,min(w,mu-bw*sigma),w)$val)
    F.1 <- Vectorize(function(w)integrate(f.1,min(w,mu-bw*sigma),w)$val)
    ## auc <- integrate(function(x)F.0(x)*f.1(x),mu-bw*sigma,mu+bw*sigma)$val
    auc <- auc.probit(mu,sigma^2)
    auc.hajek(F.0,F.1,x.0%*%beta,x.1%*%beta,auc,terms.only,IID)
}
    

## returns the coefficients gamma in Phi(gamma.0+gamma^Tx) = E(
## Phi(beta.0+beta.x^Tx+beta.y^Ty) | x ) where (x,y) are multivariate
## normal (mu,Sigma), beta=c(beta.0,beta.x,beta.y), and p.reduced is
## the length of beta.x
probit.coef.reduced <- function(mu,Sigma,beta,p.reduced) {
    ## browser()
    ## op <- options(warn=-1)
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
    denom <- sqrt(1+t(beta.y)%*%(Sigma.yy-Sigma.yx%*%solve(Sigma.xx)%*%Sigma.xy)%*%beta.y)
    beta.0.reduced <- beta.0 + t(beta.y)%*%mu.y - t(beta.y)%*%Sigma.yx%*%solve(Sigma.xx)%*%mu.x
    beta.x.reduced <- t(beta.x)+t(beta.y)%*%Sigma.yx%*%solve(Sigma.xx)
    beta.reduced <- c(beta.0.reduced,beta.x.reduced) / c(denom)
    return(beta.reduced)
    ## on.exit(options(op))
}

## x assumed in model matrix format
infl.probit <- function(x,g,beta=NULL) {
    if(is.null(beta)) beta <- coef(glm(g~x-1,family=binomial('probit')))
    n <- nrow(x)
    eta <- as.numeric(x%*%beta)
    scores <- dnorm(eta)*(g/pnorm(eta)/(1-pnorm(eta)) - 1/(1-pnorm(eta))) * x
    W <- diag(dnorm(eta)^2/pnorm(eta)/(1-pnorm(eta)))
    hessian <- -t(x)%*%W%*%x / n
    bread <- -solve(hessian)
    infl <- bread%*%t(scores)
}

## ## deprecate
## infl.hat.probit <- function(glm0) {
##     x <- model.matrix(glm0)
##     y <- glm0$y
##     n <- nrow(x)
##     eta.hat <- predict(glm0)
##     scores <- dnorm(eta.hat)*(y/pnorm(eta.hat)/(1-pnorm(eta.hat)) - 1/(1-pnorm(eta.hat))) * x
##     W <- diag(dnorm(eta.hat)^2/pnorm(eta.hat)/(1-pnorm(eta.hat)))
##     hessian <- -t(x)%*%W%*%x / n
##     bread <- -solve(hessian)
##     infl <- bread%*%t(scores)
## }    
## ## ## infl.hat.probit is just infl.probit evaluated at beta=beta.hat
## ## set.seed(1)
## ## n <- 30
## ## p <- 6
## ## ## alpha <- .1
## ## beta <- runif(p+1)/5
## ## x <- cbind(1,matrix(runif(n*p),ncol=p))
## ## risk <- pnorm(x%*%beta)
## ## y <- rbinom(n,1,risk)
## ## glm0 <- glm(y~x-1,family=binomial(link='probit'))
## ## old <- infl.hat.probit(glm0)
## ## new <- infl.probit(x,y,coef(glm0))
## ## old-new


## ## expects x in model matrix format
## infl.lda <- function(x,g,params,terms.only=TRUE) {
##     ## browser()
##     x <- t(x)#[-1,]
##     mu.0 <- params$mu.0
##     mu.1 <- params$mu.1
##     Sigma <- params$Sigma
##     ## mu.0 <- c(1,params$mu.0)
##     ## mu.1 <- c(1,params$mu.1)
##     ## Sigma <- cbind(0,rbind(0,params$Sigma))
##     n <- length(g)
##     n.1 <- sum(g); n.0 <- n-n.1
##     x.0 <- t(x[,g==0]); x.1 <- t(x[,g==1])
##     g <- sort(g)
##     mu.0.hat <- rowMeans(x.0);  mu.1.hat <- rowMeans(x.1)
##     ## pi.1.hat <- n.1/(n.0+n.1)
##     ## infl.0 <-   g/pi.1/(1-pi.1)-1/(1-pi.1) - t(mu.1)%*%solve(Sigma)%*%(x-mu.1) * n/n.1 * g  + t(mu.0)%*%solve(Sigma)%*%(x-mu.0) * n/n.0 * (1-g)
##     infl <- t(g*t(solve(Sigma)%*%(x-mu.1)*n/n.1)  -  (1-g)*t(solve(Sigma)%*%(x-mu.0)*n/n.0))
##     ## infl <- rbind(infl.0,infl)
##     if(terms.only) return(infl) else return(rowMeans(infl))
## }


## expects x in model matrix format. adding support for heteroskedasticity.
infl.lda <- function(x,g,params,terms.only=TRUE) {
    ## browser()
    x <- t(x)#[-1,]
    mu.0 <- params$mu.0
    mu.1 <- params$mu.1
    pi.0 <- params$pi.0
    pi.1 <- 1-pi.0
    Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
    if(is.null(Sigma.0))Sigma.0 <- params$Sigma
    if(is.null(Sigma.1))Sigma.1 <- params$Sigma
    Sigma <- pi.0*Sigma.0 + pi.1*Sigma.1
    ## Sigma <- params$Sigma
    ## mu.0 <- c(1,params$mu.0)
    ## mu.1 <- c(1,params$mu.1)
    ## Sigma <- cbind(0,rbind(0,params$Sigma))
    n <- length(g)
    n.1 <- sum(g); n.0 <- n-n.1
    x.0 <- t(x[,g==0]); x.1 <- t(x[,g==1])
    g <- sort(g)
    mu.0.hat <- rowMeans(x.0);  mu.1.hat <- rowMeans(x.1) # clean up
    ## pi.1.hat <- n.1/(n.0+n.1)
    ## infl.0 <-   g/pi.1/(1-pi.1)-1/(1-pi.1) - t(mu.1)%*%solve(Sigma)%*%(x-mu.1) * n/n.1 * g  + t(mu.0)%*%solve(Sigma)%*%(x-mu.0) * n/n.0 * (1-g)
    infl <- t(g*t(solve(Sigma)%*%(x-mu.1)*n/n.1)  -  (1-g)*t(solve(Sigma)%*%(x-mu.0)*n/n.0))
    ## infl <- rbind(infl.0,infl)
    if(terms.only) return(infl) else return(rowMeans(infl))
}

## (23/1/13) expects x in model matrix format. making hetero/homoskedasticity
## explicit. adding empirical estimates to allow for params=NULL. 
infl.lda <- function(x,d,params,var.equal=TRUE,terms.only=TRUE) {
    ## browser()
    g <- d # clean up
    mu.0 <- params$mu.0
    mu.1 <- params$mu.1
    pi.0 <- params$pi.0
    x.0 <- x[g==0,]; x.1 <- x[g==1,]
    if(is.null(mu.0))mu.0 <- colMeans(x.0)
    if(is.null(mu.1))mu.1 <- colMeans(x.1)
    if(is.null(pi.0))pi.0 <- nrow(x.0)/(nrow(x.0)+nrow(x.1))
    pi.1 <- 1-pi.0
    if(var.equal) {
        Sigma <- params$Sigma
        if(is.null(Sigma))...
    } else {
        Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
        if(is.null(Sigma.0))Sigma.0 <- var(x.0)
        if(is.null(Sigma.1))Sigma.1 <- var(x.1)
        Sigma <- pi.0*Sigma.0 + pi.1*Sigma.1
    }
    x <- t(x)#[-1,]
    n <- length(g) # clean up redundancies here
    n.1 <- sum(g); n.0 <- n-n.1
    ## x.0 <- t(x[,g==0]); x.1 <- t(x[,g==1])
    g <- sort(g)
    ## mu.0.hat <- rowMeans(x.0);  mu.1.hat <- rowMeans(x.1)
    ## pi.1.hat <- n.1/(n.0+n.1)
    ## infl.0 <-   g/pi.1/(1-pi.1)-1/(1-pi.1) - t(mu.1)%*%solve(Sigma)%*%(x-mu.1) * n/n.1 * g  + t(mu.0)%*%solve(Sigma)%*%(x-mu.0) * n/n.0 * (1-g)
    infl <- t(g*t(solve(Sigma)%*%(x-mu.1)*n/n.1)  -  (1-g)*t(solve(Sigma)%*%(x-mu.0)*n/n.0))
    ## infl <- rbind(infl.0,infl)
    if(terms.only) return(infl) else return(rowMeans(infl))
}

## below auc.probit and derivative routines were not useful. they
## assumed the beta used to generate the risks is the same beta used
## for the indexes i the auc.

## ## deprecate
## auc.probit <- function(mu,sigma,bw=8) {
##     E.g <- pnorm(mu/sqrt(1+sigma^2))
##     f.case <- function(w)pnorm(w)*dnorm((w-c(mu))/c(sigma))/c(sigma*E.g)
##     f.control <- function(w)(1-pnorm(w))*dnorm((w-c(mu))/c(sigma))/c(sigma*(1-E.g))
##     int.inner <- Vectorize(function(v)integrate(function(u)f.control(u),mu-bw*sigma,v)$val)
##     integrate(function(v)f.case(v)*int.inner(v),mu-bw*sigma,mu+bw*sigma)$val
## }


## ## if beta==NULL, return P(x.0<x.1) where x ~ N(mean.x,var.x),
## ## g=bernoulli(pnorm(x)), x.i ~ x|g=i. If beta!=NULL, first transform
## ## mean.x by t(beta)%*%mean.x and var.x by t(beta)%*%var.x%*%beta.
## auc.probit <- function(mean.x,var.x,beta=NULL,bw=8) {
##     if(!is.null(beta)) {
##         stopifnot(length(beta)==length(mean.x))
##         mean.x <- beta%*%mean.x
##         var.x <- t(beta)%*%var.x%*%beta
##     } else     stopifnot(length(mean.x)==1 && length(var.x)==1)
##     E.g <- pnorm(mean.x/sqrt(1+var.x))
##     mu <- mean.x; sigma <- sqrt(var.x) # clean up
##     f.case <- function(w)pnorm(w)*dnorm((w-c(mu))/c(sigma))/c(sigma*E.g)
##     f.control <- function(w)(1-pnorm(w))*dnorm((w-c(mu))/c(sigma))/c(sigma*(1-E.g))
##     int.inner <- Vectorize(function(v)integrate(function(u)f.control(u),mu-bw*sigma,v)$val)
##     integrate(function(v)f.case(v)*int.inner(v),mu-bw*sigma,mu+bw*sigma)$val
## }


## ## ## check auc.probit.
## ## require(mvtnorm)
## ## start <- Sys.time()
## ## source('misc.R')
## ## set.seed(3)
## ## p <- 2
## ## ## n <- 1e3
## ## ns <- round(seq(1e2,7e2,len=2e1))
## ## beta <- runif(p)/4
## ## Sigma <- matrix(runif(p^2),nrow=p)
## ## Sigma <- Sigma%*%t(Sigma)
## ## mu <- runif(p)
## ## auc <- auc.probit(mu,Sigma,beta)
## ## by.n <- sapply(ns, function(n) {
## ##     cat('.')
## ##     diffs <- replicate(1e2, {
## ##         x <- rmvnorm(n,mu,Sigma)
## ##         risk <- pnorm(x%*%beta)
## ##         g <- rbinom(n,1,risk)
## ##         x.0 <- x[g==0,]; x.1 <- x[g==1,]
## ##         auc.hat(x.0,x.1) - auc
## ##     })
## ## })
## ## Sys.time() - start
## ## mad <- colMeans(abs(by.n))
## ## plot(ns,mad)
## ## lm0 <- lm(log(mad)~log(ns))
## ## coef(lm0)
## ## curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)



## ## derivative of auc.probit(0,sqrt(t(beta)%*%Sigma%*%beta)) wrt
## ## beta. not useful if probit model has an intercept term.
## auc.probit.deriv <- function(beta,mu,Sigma,bw=8) {
##     f <- function(sigma)auc.probit(0,sigma)
##     sigma <- as.numeric(sqrt(t(beta)%*%Sigma%*%beta))
##     E.g <- 1/2
##     f.case <- function(w)pnorm(w)*dnorm(w/sigma)*2/sigma
##     f.control <- function(w)(1-pnorm(w))*dnorm(w/sigma)*2/sigma
##     A <- -2/sigma*f(sigma)
##     int.inner <- Vectorize(function(v)integrate(function(u)f.control(u),-bw*sigma,v)$val)
##     B <- int.outer <- integrate(function(v)f.case(v)*v^2/sigma^3*int.inner(v),-bw*sigma,+bw*sigma)$val
##     int.inner <- Vectorize(function(v)integrate(function(u)f.control(u)*u^2/sigma^3,-bw*sigma,v)$val)
##     C <- int.outer <- integrate(function(v)f.case(v)*int.inner(v),-bw*sigma,+bw*sigma)$val
##     scalar.part <- A+B+C
##     scalar.part*t(beta)%*%Sigma/sigma
## }


## ## derivative of
## ## auc.probit(t(beta)%*%mu,sqrt(t(beta)%*%Sigma%*%beta)) wrt beta
## auc.probit.deriv <- function(beta,mu,Sigma) {
##     f.mu.prime <- function(mu,sigma,bw=7) {
##         E.g <- function(mu)pnorm(mu/sqrt(1+sigma^2))
##         A <- function(mu,sigma)sigma^2*E.g(mu)*(1-E.g(mu))
##         A.mu.prime <- function(mu)sigma^2/sqrt(1+sigma^2)*dnorm(mu/sqrt(1+sigma^2))*(1-2*E.g(mu))
##         int.inner <- Vectorize(function(v)integrate(function(u)pnorm(v)*(1-pnorm(u))*dnorm((v-mu)/sigma)*dnorm((u-mu)/sigma)*(-A.mu.prime(mu)/A(mu,sigma)^2+(v+u-2*mu)/A(mu,sigma)/sigma^2),mu-bw*sigma,v)$val)
##         int.outer <- integrate(int.inner,mu-bw*sigma,mu+bw*sigma)$val
##         int.outer
##     }
##     f.sigma.prime <- function(mu,sigma,bw=7) {
##         E.g <- function(sigma)pnorm(mu/sqrt(1+sigma^2))
##         A <- function(mu,sigma)sigma^2*E.g(sigma)*(1-E.g(sigma))
##         A.sigma.prime <- function(sigma)2*sigma*E.g(sigma)*(1-E.g(sigma)) - sigma^2*dnorm(mu/sqrt(1+sigma^2))*mu*sigma/(1+sigma^2)^(3/2)*(1-E.g(sigma)) + sigma^2*E.g(sigma)*dnorm(mu/sqrt(1+sigma^2))*mu*sigma/(1+sigma^2)^(3/2)
##         int.inner <- Vectorize(function(v)integrate(function(u)pnorm(v)*(1-pnorm(u))*dnorm((v-mu)/sigma)*dnorm((u-mu)/sigma) * (-A.sigma.prime(sigma)/A(mu,sigma)^2 + ((v-mu)^2+(u-mu)^2)/A(mu,sigma)/sigma^3   ),mu-bw*sigma,v)$val)
##         int.outer <- integrate(int.inner,mu-bw*sigma,mu+bw*sigma)$val
##         int.outer
##     }
##     index.mu <- as.numeric(t(beta)%*%mu)
##     index.sigma <- as.numeric(sqrt(t(beta)%*%Sigma%*%beta))
##     mu.deriv <- f.sigma.prime(index.mu,index.sigma)*t(beta)%*%Sigma/index.sigma
##     sigma.deriv <- f.mu.prime(index.mu,index.sigma)*t(mu)
##     matrix(mu.deriv+sigma.deriv,ncol=1)
## }
## ## mu <- runif(p)
## ## Sigma <- matrix(runif(p^2),p)
## ## Sigma <- Sigma%*%t(Sigma)
## ## f <- function(beta)auc.probit(t(beta)%*%mu,t(beta)%*%Sigma%*%beta)
## ## f.prime <- function(beta)auc.probit.prime(beta,mu,Sigma)
## ## viz.deriv(f,f.prime,p,1)



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
    ## stop()
}

auc.probit <- function(mu,Sigma,beta,gamma,bw=5){
    lb <- gamma%*%mu-bw*t(gamma)%*%Sigma%*%gamma
    ub <- gamma%*%mu+bw*t(gamma)%*%Sigma%*%gamma
    ## inner.int <- Vectorize(function(v)integrate(function(w)f.gamma.0(w,mu,Sigma,gamma,beta),lb,v)$val)
    ## try <- integrate(function(v)inner.int(v)*f.gamma.1(v,mu,Sigma,gamma,beta),round(lb),round(ub))$val
    inner.int <- Vectorize(function(v)integrate(function(w)pdf.index.probit(w,0,mu,Sigma,gamma,beta),lb,v)$val)
    try <- integrate(function(v)inner.int(v)*pdf.index.probit(v,1,mu,Sigma,gamma,beta),round(lb),round(ub))$val
}


## ## P(x.0<x.1) where x.0~N(mu.0,var.x), x.1~N(mu.1,var.x),
## ## mean.x.diff=mu.1-mu.0. if beta is not null, then first transform x
## ## -> t(beta)%*%x.
## auc.lda <- function(mean.x.diff,var.x,beta=NULL) {
##     ## mu.diff <- c(0,mu.diff)
##     ## Sigma.diff <- 2*cbind(0,rbind(0,Sigma))
##     ## pnorm((coefs%*%mu.diff)/sqrt(t(coefs)%*%(2*Sigma)%*%coefs))
##     if(!is.null(beta)) {
##         stopifnot(length(beta)==length(mean.x.diff))
##         mean.x.diff <- beta%*%mean.x.diff
##         var.x <- t(beta)%*%var.x%*%beta
##     } else     stopifnot(length(mean.x.diff)==1 && length(var.x)==1)
                  
##     pnorm(mean.x.diff/sqrt(var.x)/2)
## }

## P(x.0<x.1) where x.0~N(mu.0,var.x), x.1~N(mu.1,var.x),
## mean.x.diff=mu.1-mu.0. if beta is not null, then first transform x
## -> t(beta)%*%x. test: 11a
auc.lda <- function(beta=NULL,params){
    ## mu.diff <- with(params, c(0,mu.1 - mu.0))
    ## Sigma <- cbind(0,rbind(0,params$Sigma))
    mu.diff <- with(params, mu.1 - mu.0)
    Sigma <- params$Sigma
    if(!is.null(beta)) {
        stopifnot(length(beta)==length(mu.diff))
        mu.diff <- t(beta)%*%mu.diff
        Sigma <- t(beta)%*%Sigma%*%beta
    } else     stopifnot(length(mu.diff)==1 && length(Sigma)==1)                  
    pnorm(mu.diff/sqrt(Sigma*2))
}



## ## ignoring interceptt term
## auc.lda <- function(beta=NULL,params){
##     mu.diff <- with(params, mu.1 - mu.0)
##     Sigma <- params$Sigma
##     if(!is.null(beta)) {
##         stopifnot(length(beta)==length(mu.diff))
##         mu.diff <- beta%*%mu.diff
##         Sigma <- t(beta)%*%Sigma%*%beta
##     } else     stopifnot(length(mu.diff)==1 && length(Sigma)==1)                  
##     pnorm(mu.diff/sqrt(Sigma)/2)
## }


## auc.scores.deriv <-  function(coefs,mu.diff,Sigma.diff) {
##     quad <- as.numeric(coefs%*%Sigma.diff%*%coefs)
##     -(Sigma.diff%*%coefs %*% (coefs%*%mu.diff/quad) - mu.diff)%*%(1/sqrt(quad)) %*% dnorm(coefs%*%mu.diff/sqrt(quad))
## }

## auc.lda.deriv <-  function(coefs,mean.x.diff,var.x) {
##     var.x <- var.x*2
##     quad <- as.numeric(coefs%*%var.x%*%coefs)
##     -(var.x%*%coefs %*% (coefs%*%mean.x.diff/quad) - mean.x.diff)%*%(1/sqrt(quad)) %*% dnorm(coefs%*%mean.x.diff/sqrt(quad))
## }

## auc.deriv.lda <-  function(beta,params) {
##     ## browser()
##     mu.diff <- with(params, mu.1 - mu.0)
##     Sigma <- params$Sigma
##     Sigma <- Sigma*2
##     quad <- as.numeric(t(beta)%*%Sigma%*%beta)
##     -(Sigma%*%beta %*% (t(beta)%*%mu.diff/quad) - mu.diff)%*%(1/sqrt(quad)) %*% dnorm(t(beta)%*%mu.diff/sqrt(quad))
## }

## allow for heteroscedasticity (misspecified lda)
auc.deriv.lda.gaussian <-  function(beta,params) {
    ## browser()
    mu.diff <- with(params, mu.1 - mu.0)
    Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
    if(is.null(Sigma.0))Sigma.0 <- params$Sigma
    if(is.null(Sigma.1))Sigma.1 <- params$Sigma
    Sigma <- Sigma.0 + Sigma.1
    ## Sigma <- Sigma*2
    quad <- as.numeric(t(beta)%*%Sigma%*%beta)
    t(  -(Sigma%*%beta %*% (t(beta)%*%mu.diff/quad) - mu.diff)%*%(1/sqrt(quad)) %*% dnorm(t(beta)%*%mu.diff/sqrt(quad))  )
}



## deprecate. make naming consisntent.
auc.lda.deriv <- auc.deriv.lda.gaussian

################ dev ops

## routine to visually check conditional expectation on a continuous
## variable
viz.cond.exp <- function(x,y,f,n.quantiles=1e2,...) {
    ## browser()
    bin <- cut(x,breaks=c(-Inf,quantile(x,probs=seq(0,1,len=n.quantiles)),Inf))
    y.means <- aggregate(y,list(bin=bin),mean)[,2]
    bin.centers <- aggregate(x,list(bin=bin),median)[,2]
    plot(bin.centers,y.means,xlab='',asp=1,...)
    curve(f,add=TRUE,col=2)
}


## routine to check derivative of a vector-valued f wrt a vector
## x. using a random direction and random starting point if none
## passed.
viz.deriv <- function(f,deriv,dim.x,dim.f,x0=NULL,delta=NULL) {
    if(is.null(x0)) x0 <- runif(dim.x)
    if(is.null(delta)) delta <- matrix(runif(dim.x),ncol=1)
    a <- matrix(runif(dim.f),ncol=1)
    ts <- seq(0,1,len=1e2)
    fs <- sapply(ts,function(t)t(a)%*%f(x0+t*delta))
    plot(ts,fs,type='l')
    deriv.scalar <- t(a)%*%deriv(x0)%*%delta
    abline(a=fs[1],b=deriv.scalar,col=2,lty=2)
}
## dim.x <- sample(1:10,1)
## dim.f <- sample(1:10,1)
## A <- matrix(runif(dim.x*dim.f),ncol=dim.x)
## f <- function(x)A%*%x
## plot.matrix.deriv(f,function(x)A,dim.x,dim.f)

## deprecate
plot.matrix.deriv <- viz.deriv
cond.exp.viz <- viz.cond.exp

################




hypergeom1F1 <- function(a,b,z)
    gamma(b)/gamma(a)/gamma(b-a)*integrate(function(u)exp(z*u)*u^(a-1)*(1-u)^(b-a-1),0,1)$val

dgammasum <- Vectorize(function(x,shape1,rate1,shape2,rate2) {
    if(x<=0) {
        return(0)
    } else {
        return( rate1^shape1*rate2^shape2/gamma(shape1+shape2)*exp(-rate1*x)*x^(shape1+shape2-1)*hypergeom1F1(shape2,shape1+shape2,(rate1-rate2)*x) )
    }
}, vectorize.args='x')  

pgammasum <- Vectorize(function(q,shape1,rate1,shape2,rate2)integrate(function(x)dgammasum(x,shape1,rate1,shape2,rate2),0,q)$val,, vectorize.args='q')


## P(x<y) where x ~ sum_i beta_i x_i ~ hypoexp(rate.0/beta), y ~ sum_i
## beta_i y_i, x_i ~ Exp(rate.0) ~ hypoexp(rate.1/beta), x_i ~
## Exp(rate.0), y_i ~ Exp(rate.1).
## TODO should take a params for model params 
auc.lda.gamma <- function(beta,rate.0,rate.1) {
    beta <- as.numeric(beta)
    n.unique <- length(unique(beta))
    if(n.unique==length(beta)) {
        ## cat('3')
        ul <- Inf#sum(beta/rate.1) + sqrt(sum((beta / rate.1)^2))*5
        return( integrate(function(x)sdprisk::phypoexp(x,rate=rate.0/beta)*sdprisk::dhypoexp(x,rate=rate.1/beta),0,ul)$val )
        }
    if(n.unique==1) {
        ## cat('1')
        ul <- Inf
        return( with(list(shape=length(beta)),
                     integrate(function(x)pgamma(x,shape,rate=rate.0)*dgamma(x,shape,rate=rate.1),0,ul)$val ) )
        }
    if(n.unique==2) {
        ## cat('2')
        rle0 <- rle(sort(beta))
        count <- rle0$len; coef <- rle0$val
        ul <- Inf#sum(count*coef/rate.1) + 5*sum(sqrt(count*coef)/rate.1)
        ## print(ul)
        return( integrate(function(x)pgammasum(x,count[1],rate.0/coef[1],count[2],rate.0/coef[2])*dgammasum(x,count[1],rate.1/coef[1],count[2],rate.1/coef[2]),0,ul)$val )#sum(count*coef/rate.1) + 3*sum(sqrt(count*coef)/rate.1))$val )
        }
    stop()
}    





## returns the mle coefficients in the model P(g=1)=link(beta^t x.red)
## assuming g has been generated under the model
## P(g=1)=link(params$beta^t x). numerically maximizes the expected
## log likelihood given x, which is itself obtain by numerical
## integration. test: 14a
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


## inner <- Vectorize(function(w)integrate(function(w0)h(w0)*dnorm(w.0,cond.mean(w,beta),sqrt(cond.var(beta)))/sqrt(cond.var(beta)),-lim,lim)$val)

## conditional pdf of the index, beta^tx | g=i, assuming model P(g=i)
## = link(params$beta^t x). test: 14b.
pdf.index.glm.gaussian <- function(x,g,beta,params,lim=Inf) {
    ## browser()
    h <- params$link
    w <- x # clean up
    mu <- params$mu
    Sigma <- params$Sigma
    beta.0 <- params$beta
    ## p.full <- length(beta.0)
    f.w <- function(w)dnorm(w,t(beta)%*%mu, sd=sqrt(t(beta)%*%Sigma%*%beta))
    f.w0 <- function(w0)dnorm(w0,t(beta.0)%*%mu, sd=sqrt(t(beta.0)%*%Sigma%*%beta.0))
    E.g <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    if(isTRUE(all.equal(beta,beta.0))) return(h(x)^g*(1-h(x))^(1-g)*f.w0(x)/(E.g^g*(1-E.g)^(1-g)))
    quad.ww <- as.numeric(t(beta)%*%Sigma%*%beta)
    quad.ww0 <- as.numeric(t(beta.0)%*%Sigma%*%beta)
    quad.w0w0 <- as.numeric(t(beta.0)%*%Sigma%*%beta.0)
    f.w0.cond <- function(w0,w)dnorm(w0, mean=t(beta.0)%*%mu + quad.ww0/quad.ww*(w-t(beta)%*%mu), sd=sqrt((1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0))
    f.w.w0.cond <- function(w,w0,g) h(w0)^g*(1-h(w0))^(1-g)*f.w0.cond(w0,w)*f.w(w)/(E.g^g*(1-E.g)^(1-g))
    ## f.w.cond <- Vectorize( function(w,g) integrate(function(w0)f.w.w0.cond(w,w0,g),-lim,lim)$val, vectorize.args='w')
    ## sapply(w, function(w) integrate(function(w0)f.w.w0.cond(w,w0,g),-lim,lim)$val)
    sapply(w, function(w) integrate(function(w0)f.w.w0.cond(w,w0,g),-lim,lim)$val)
    ## Sigma.cond <- (1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0
    ## f.w0.cond <- function(w0,w)dnorm(w0, mean=t(beta.0)%*%mu + quad.ww0/quad.ww*(w-t(beta)%*%mu), sd=sqrt(Sigma.cond))
    ## f.w.w0.cond <- function(w,w0,g) h(w0)^g*(1-h(w0))^(1-g)*f.w0.cond(w0,w)*f.w(w)/(E.g^g*(1-E.g)^(1-g))
    ## sapply(w, function(w) integrate(function(w0)h(w0)^g*(1-h(w0))^(1-g)*f.w0.cond(w0,w)*f.w(w)/(E.g^g*(1-E.g)^(1-g)),-lim,lim)$val)
}
## pdf.index.glm.gaussian <- Vectorize(pdf.index.glm.gaussian,vectorize.args='x')


## conditional pdf of the index, beta^tx | g=i, assuming model P(g=i)
## = link(params$beta^t x). test: 14b.
## !!! this version assumes h is the loistic link, but is faster than the general version
pdf.index.logit.gaussian <- function(x,g,beta,params,lim=Inf) {
    ## browser()
    h <- params$link
    w <- x # clean up
    mu <- params$mu
    Sigma <- params$Sigma
    beta.0 <- params$beta
    f.w <- function(w)dnorm(w,t(beta)%*%mu, sd=sqrt(t(beta)%*%Sigma%*%beta))
    f.w0 <- function(w0)dnorm(w0,t(beta.0)%*%mu, sd=sqrt(t(beta.0)%*%Sigma%*%beta.0))
    E.g <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
    if(isTRUE(all.equal(beta,beta.0))) return(h(x)^g*(1-h(x))^(1-g)*f.w0(x)/(E.g^g*(1-E.g)^(1-g)))
    quad.ww <- as.numeric(t(beta)%*%Sigma%*%beta)
    quad.ww0 <- as.numeric(t(beta.0)%*%Sigma%*%beta)
    quad.w0w0 <- as.numeric(t(beta.0)%*%Sigma%*%beta.0)
    Sigma.cond <- sqrt((1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0) 
    ## fExp <- function(x) plogis(x) * dnorm(x, mean = mu, sd = sigma)
    ## .exp <- integrate(fExp, -Inf, Inf, abs.tol = abs.tol, ...)$value    
    ## ln.mean <-  sapply(w, function(w)logitnorm::momentsLogitnorm(mu=t(beta.0)%*%mu + quad.ww0/quad.ww*(w-t(beta)%*%mu), sigma=Sigma.cond))
    ln.mean <-  sapply(w, function(w)integrate(function(w0)plogis(w0)*dnorm(w0,mean=sum(beta.0*mu) + quad.ww0/quad.ww*(w-sum(beta*mu)), sd=Sigma.cond),-lim,lim)$val)
    ln.mean^g*(1-ln.mean)^(1-g)*f.w(w)/(E.g^g*(1-E.g)^(1-g))
}



## dd



## infl.glm.gaussian <- function(x,g,params,terms.only=TRUE) {
##     p <- params$p
##     score <- function(x,g,params) { # x in model matrix format
##         h <- params$link; h.1 <- params$lnk.deriv
##         eta <- as.numeric(x%*%params$beta)
##         t( (g-h(eta))/h(eta)/(1-h(eta))*h.1(eta) * x )
##     }
##     fi <- function(x,g,params) { # x in model matrix format
##         h <- params$link; h.1 <- params$lnk.deriv; h.2 <- params$link.deriv2
##         eta <- as.numeric(x%*%params$beta)
##         denom <- h(eta)*(1-h(eta))
##         -t((h.2(eta)*(g-h(eta))/denom - h.1(eta)/denom^2 * (g*(h.1(eta)-2*h(eta)*h.1(eta))+h.1(eta)*h(eta)^2))*x) %*% x
##     }
##     terms <- solve(fi(x,g,params)[1:p,1:p])%*%score(x,g,params)[1:p,]
##     if (terms.only) return(terms) else return(rowSums(terms))
## }


## properly normalizing score. in the old version, the terms are the
## terms of the mean (divided by 1/n). routines were expecting just
## the terms, so infl=mean(terms)
## TODO: p.red and beta currently stored inside params but
## not really part of the data generation.

infl.glm <- function(x,g,params,terms.only=TRUE) {
    ## browser()
    p <- params$p
    score <- function(x,g,params) { # x in model matrix format
        h <- params$link; h.1 <- params$lnk.deriv
        eta <- as.numeric(x%*%params$beta)
        t( (g-h(eta))/h(eta)/(1-h(eta))*h.1(eta) * x )
    }
    fi <- function(x,g,params) { # x in model matrix format
        h <- params$link; h.1 <- params$lnk.deriv; h.2 <- params$link.deriv2
        n <- nrow(x)
        eta <- as.numeric(x%*%params$beta)
        denom <- h(eta)*(1-h(eta))
        -t((h.2(eta)*(g-h(eta))/denom - h.1(eta)/denom^2 * (g*(h.1(eta)-2*h(eta)*h.1(eta))+h.1(eta)*h(eta)^2))*x) %*% x / n
    }
    ## terms <- solve(fi(x,g,params)[1:p,1:p])%*%score(x,g,params)[1:p,]
    terms <- solve(fi(x,g,params))%*%score(x,g,params)
    if (terms.only) return(terms) else return(rowMeans(terms))
}

## deprecate--no assumption made about gaussian covariates 
infl.glm.gaussian <- infl.glm


## P(x.0%*%beta < x.1%*%beta | g.0=0,g.1=1), where P(g=i|x) = h(params$beta %*% x). test: 14c
auc.glm.gaussian <- function(beta,params,lim=Inf) {
    ## h <- link
    ## mu <- params$mu
    ## Sigma <- params$Sigma
    ## beta.0 <- params$beta
    ## browser()
    if(length(beta)<length(params$beta)) beta <- c(beta,rep(0,length(params$beta)-length(beta)))
    f.0 <- function(x)pdf.index.glm.gaussian(x,g=0,beta=beta,params=params,lim=lim)
    f.1 <- function(x)pdf.index.glm.gaussian(x,g=1,beta=beta,params=params,lim=lim)
    F.0 <- Vectorize(function(x)integrate(function(w)f.0(w),-lim,x)$val)
    integrate(function(x)F.0(x)*f.1(x),-lim,lim)$val
}



## P(x.0%*%beta < x.1%*%beta | g.0=0,g.1=1), where P(g=i|x) = h(params$beta %*% x). test: 14c
auc.logit.gaussian <- function(beta,params,lim=Inf) {
    ## h <- link
    ## mu <- params$mu
    ## Sigma <- params$Sigma
    ## beta.0 <- params$beta
    ## browser()
    if(length(beta)<length(params$beta)) beta <- c(beta,rep(0,length(params$beta)-length(beta)))
    f.0 <- function(x)pdf.index.logit.gaussian(x,g=0,beta=beta,params=params,lim=lim)
    f.1 <- function(x)pdf.index.logit.gaussian(x,g=1,beta=beta,params=params,lim=lim)
    F.0 <- Vectorize(function(x)integrate(function(w)f.0(w),-lim,x)$val)
    integrate(function(x)F.0(x)*f.1(x),-lim,lim)$val
}



## h <- link
## mu <- params$mu
## Sigma <- params$Sigma
## beta.0 <- params$beta
## quad.ww <- as.numeric(t(beta)%*%Sigma%*%beta)
## quad.ww0 <- as.numeric(t(beta.0)%*%Sigma%*%beta)
## quad.w0w0 <- as.numeric(t(beta.0)%*%Sigma%*%beta.0)
## Sigma.cond <- (1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0
## mu.1 <- sum(beta*mu)
## sigma.1 <- sqrt(quad.ww)
## mu.2 <- sum(beta.0*mu)+quad.ww0/quad.ww*(w-sum(beta*mu))
## sigma.2 <- sqrt(Sigma.cond)
## integrate(function(u)h(u)*dnorm(u-


## auc.new <- function(beta,params,lim=Inf) {
##     h <- params$link
##     mu <- params$mu
##     Sigma <- params$Sigma
##     beta.0 <- params$beta
##     quad.ww <- as.numeric(t(beta)%*%Sigma%*%beta)
##     quad.ww0 <- as.numeric(t(beta.0)%*%Sigma%*%beta)
##     quad.w0w0 <- as.numeric(t(beta.0)%*%Sigma%*%beta.0)
##     f.w <- function(w)dnorm(w,t(beta)%*%mu, sd=sqrt(quad.ww))
##     f.w0 <- function(w0)dnorm(w0,t(beta.0)%*%mu, sd=sqrt(quad.w0w0))
##     E.g <- integrate(function(w0)h(w0)*f.w0(w0),-Inf,Inf)$val
##     if(all.equal(beta,beta.0)==TRUE) return(h(x)^g*(1-h(x))^(1-g)*f.w0(x)/(E.g^g*(1-E.g)^(1-g)))
##     f.w0.cond <- function(w0,w)dnorm(w0, mean=t(beta.0)%*%mu + quad.ww0/quad.ww*(w-t(beta)%*%mu), sd=sqrt((1-quad.ww0^2/quad.ww/quad.w0w0)*quad.w0w0))
##     f.w.w0.cond <- function(w,w0,g) h(w0)^g*(1-h(w0))^(1-g)*f.w0.cond(w0,w)*f.w(w)/(E.g^g*(1-E.g)^(1-g))
##     int.1 <- Vectorize(function(v)integrate(function(w0)f.w0.cond(w0,v),-lim,lim)$val)
##     int.2 <- Vectorize(function(w)integrate(function(v)(1-h(v))*f.w(v)*int.1(v),-lim,w)$val)
##     integrate(function(w)h(w)*f.w(w)*int.1(w)*int.2(w),-lim,lim)$val / E.g / (1-E.g)
## }




## P(x.0<x.1) where x.0~N(mu.0,Sigma.0), x.1~N(mu.1,Sigma.1). Support
## for different Sigmas for testing though heteroscedasticity would
## violate lda assumptions. test:15a
auc.lda.gaussian <- function(beta,params){
    mu.diff <- with(params, mu.1 - mu.0)
    Sigma <- with(params, Sigma.0+Sigma.1)
    pnorm(t(beta)%*%mu.diff / sqrt(t(beta)%*%Sigma%*%beta))
}


auc.hajek.lda.gaussian <- function(x,g,beta,params,terms.only=TRUE,IID=FALSE) {
    x.0 <- x[g==0,]; x.1 <- x[g==1,]
    with(params, {
        F <- function(x)pnorm((x-t(beta)%*%mu.0)/sqrt(t(beta)%*%Sigma.0%*%beta))
        G <- function(y)pnorm((y-t(beta)%*%mu.1)/sqrt(t(beta)%*%Sigma.1%*%beta))
        auc.hajek((x.0%*%beta)[,],(x.1%*%beta)[,],F=F,G=G,auc=auc.lda.gaussian(beta,params),terms.only=terms.only,IID=IID)
    })
    }
    
auc.hajek.logit.gaussian <- function(x,g,beta,params,lim=lim,terms.only=TRUE,IID=FALSE) {
    x.0 <- x[g==0,]; x.1 <- x[g==1,]
    f.0 <- function(x)pdf.index.logit.gaussian(x,g=0,beta=beta.,params=params,lim=lim) 
    f.1 <- function(x)pdf.index.logit.gaussian(x,g=1,beta=beta.hat.red,params=params,lim=lim) 
    F.0 <- Vectorize(function(x)integrate(function(w)f.0(w),-lim,x)$val) 
    F.1 <- Vectorize(function(x)integrate(function(w)f.1(w),-lim,x)$val)
    ## F.0 <- function(x) {
    ##     old <- rank(x)
    ##     y <- c(-lim,sort(x))
    ##     segments <- sapply(2:length(y), function(i) integrate(function(w)f.0(w),y[i-1],y[i])$val)
    ##     return(cumsum(segments)[old])
    ## }
    ## F.1 <- function(x) {
    ##     old <- rank(x)
    ##     y <- c(-lim,sort(x))
    ##      segments <- sapply(2:length(y), function(i) integrate(function(w)f.1(w),y[i-1],y[i])$val)
    ##     return(cumsum(segments)[old])
    ## }
    auc.hajek.hat <- auc.hajek(x=x.0%*%beta,y=x.1%*%beta,F=F.0,G=F.1,auc=NULL,terms.only=terms.only,IID=IID)
}
