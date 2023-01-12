## replication of tables in binormal vs lehmann manuscript 
## adapted from postdoc/musie/vivek/sim.R
require(survival)
require(parallel)
## require(dplyr)
## require(ggplot2)
require(xtable)

auc.cox <- function(y,d,alpha=.05) {
    fit <-         coxph(Surv(y) ~ d)
    log.theta.hat <- unname(coef(fit))
    se.hat <- sqrt(vcov(fit))
    q <- qnorm(1-alpha/2)
    CI.lower <- log.theta.hat-q*se.hat
    CI.upper <- log.theta.hat+q*se.hat
    auc.hat=1/(exp(log.theta.hat)+1)
    est.cox <- c(auc.hat=auc.hat,se.hat=se.hat,CI.lower=1/(exp(CI.upper)+1),CI.upper=1/(exp(CI.lower)+1))
}

auc.binormal <- function(y,d,alpha=.05) {
    n <- sum(d) # assumes same # of diseased as non-diseased
    y.0 <- y[d==0]; y.1 <- y[d==1]
    s0 <- sd(y.0); s1 <- sd(y.1)
    auc.hat <- pnorm((mean(y.1)-mean(y.0))/sqrt(s0^2+s1^2))        
    var.hat <- 1/n + (mean(y.1)-mean(y.0))^2/2/(s0^2+s1^2)^3*(s0^4+s1^4)/(n-1)
    q <- qnorm(1-alpha/2)
    CI.lower <- qnorm(auc.hat)-q*sqrt(var.hat)
    CI.upper <- qnorm(auc.hat)+q*sqrt(var.hat)
    est.binormal <- c(auc.hat=auc.hat,se.hat=sqrt(var.hat),CI.lower=pnorm(CI.lower),CI.upper=pnorm(CI.upper))
}


reps <- 1e2
ns <- c(5,10,15,30,60,100,200)#seq(10,40,by=10)
aucs <- AUCs <- c(.55,.7,.9)#seq(.5,.9,by=.1)


scale.0 <- 1
shape <- 1
by.auc <- lapply(AUCs, function(auc) {
    scale.1 <- scale.0*(auc/(1-auc))^(1/shape)
    by.n <- mclapply(ns, mc.cores=detectCores()-4,function(n) {
        ## print(n)
        ## print(scale.1)
        auc.hats <- replicate(reps, {
            ## print('.')
            y.0 <- rweibull(n,shape=shape,scale=scale.0)
            y.1 <- rweibull(n,shape=shape,scale=scale.1)
            y <- c(y.0,y.1)
            d <- rep(c(0,1),each=n)
            rbind(cox=auc.cox(y=y,d=d),binormal=auc.binormal(y=y,d=d))
        })
        auc.hats <- simplify2array(auc.hats)
        auc.mean <- rowMeans(auc.hats[,'auc.hat',])
        auc.var <-  rowMeans(auc.hats[,'se.hat',]^2)
        summary <- rbind(mean=auc.mean,bias=auc.mean-auc,mse=(auc.mean-auc)^2+auc.var)
    })
    by.n <- simplify2array(by.n)
})
by.auc <- simplify2array(by.auc)
dimnames(by.auc)[[3]] <- ns
dimnames(by.auc)[[4]] <- aucs


auc <- aucs[1]
sd.0 <- sd.1 <- 1
mu.0 <- 0
## mean(rnorm(1e4,mu.0,sd.0) < rnorm(1e4,mu.1,sd.1))
by.auc <- lapply(aucs, function(auc) {
    ## scale.1 <- scale.0*(auc/(1-auc))^(1/shape)
    mu.1 <- mu.0 + sd.1*qnorm(auc)*sqrt(1+(sd.0/sd.1)^2)
    by.n <- mclapply(ns, mc.cores=detectCores()-4,function(n) {
        ## print(n)
        ## print(scale.1)
        auc.hats <- replicate(reps, {
            ## print('.')
            y.0 <- rnorm(n,mu.0,sd.0)#rweibull(n,shape=shape,scale=scale.0)
            y.1 <- rnorm(n,mu.1,sd.1)#rweibull(n,shape=shape,scale=scale.1)
            y <- c(y.0,y.1)
            d <- rep(c(0,1),each=n)
            rbind(cox=auc.cox(y=y,d=d),binormal=auc.binormal(y=y,d=d))
        })
        auc.hats <- simplify2array(auc.hats)
        auc.mean <- rowMeans(auc.hats[,'auc.hat',])
        auc.var <-  rowMeans(auc.hats[,'se.hat',]^2)
        summary <- rbind(mean=auc.mean,bias=auc.mean-auc,mse=(auc.mean-auc)^2+auc.var)
    })
    by.n <- simplify2array(by.n)
})
by.auc <- simplify2array(by.auc)
dimnames(by.auc)[[3]] <- ns
dimnames(by.auc)[[4]] <- aucs


## for(j in 1:10) {
options(scipen=999)
ftbl <- ftable(by.auc,row.vars=c(3,1),col.vars=c(4,2))
xftbl <- xtableFtable(ftbl,method='compact',digits=3)
print.xtableFtable(xftbl, booktabs = TRUE)









auc.cox.mle <- function(y,d,alpha=.05) {

    est.cox <- c(auc.hat=auc.hat,se.hat=se.hat,CI.lower=1/(exp(CI.upper)+1),CI.upper=1/(exp(CI.lower)+1))
}

auc <- .6
scale.0 <- 1
scale.1 <- scale.0*(auc/(1-auc))^(1/shape)


## max of log likelihood--looks OK
m <- n <- 10
shape <- 2
shape.mles <- replicate(50, {
    y.0 <- rweibull(m,shape=shape,scale=scale.0^(1/shape))
    y.1 <- rweibull(n,shape=shape,scale=scale.1^(1/shape))
    y <- c(y.0,y.1)
    objective <- function(shapes,y.0,y.1) {
        m <- length(y.0); n <- length(y.1)
        y <- c(y.0,y.1)
        sapply(shapes, function(shape) 
            ## sum(log(dweibull(y.0,shape=shape,scale=scale.0^(1/shape)))) + sum(log(dweibull(y.1,shape=shape,scale=scale.1^(1/shape))))) ## wikipedia parametrization of scale parameter
        (m+n)*log(shape)+(shape-1)*sum(log(y))-1/scale.0*sum(y.0^shape)-1/scale.1*sum(y.1^shape)-m*log(scale.0)-n*log(scale.1) ## musie parametrization
        )
    }
    optimize(objective,c(.1,5),y.0=y.0,y.1=y.1,maximum=TRUE)$maximum
})
hist(shape.mles)
abline(v=shape,col=2)


gen.data <- function(m,n,auc,shape=1,scale.0=1) {
    scale.1 <- scale.0*(auc/(1-auc))^(1/shape)
    y.0 <- rweibull(m,shape=shape,scale=scale.0^(1/shape))
    y.1 <- rweibull(n,shape=shape,scale=scale.1^(1/shape))
    return(list(y.0=y.0,y.1=y.1))
}


auc.mle <- function(y.0,y.1) {
    objective <- function(shapes,y.0,y.1) {
        m <- length(y.0); n <- length(y.1)
        y <- c(y.0,y.1)
        sapply(shapes, function(shape) 
        (m+n)/shape + sum(log(y)) - m*sum(y.0^shape*log(y.0))/sum(y.0^shape) - n*sum(y.1^shape*log(y.1))/sum(y.1^shape)
        )
    }
    shape.mle <- uniroot(objective,c(0.1,2),extendInt='yes',y.0=y.0,y.1=y.1)$root
    mean(y.1^shape.mle) / (mean(y.0^shape.mle) + mean(y.1^shape.mle))
}

## shape MLE--fixing Musie's formula
m <- n <- 50
shape <- 1
shape.mles <- replicate(50, {
    y <- gen.data(n,n,auc)
    y.0 <- y$y.0; y.1 <- y$y.1
    y <- unlist(y)
    objective <- function(shapes,y.0,y.1) {
        m <- length(y.0); n <- length(y.1)
        y <- c(y.0,y.1)
        sapply(shapes, function(shape) 
            ## (length(y)+sum(y^shape)) / (m*sum(y.0^shape*log(y.0))/sum(y.0^shape) +  n*sum(y.1^shape*log(y.1))/sum(y.1^shape)  )
        (m+n)/shape + sum(log(y)) - m*sum(y.0^shape*log(y.0))/sum(y.0^shape) - n*sum(y.1^shape*log(y.1))/sum(y.1^shape)
        )
    }
    ## objective <- Vectorize(objective,'shape')
    shape.mle <- uniroot(objective,c(0.1,2),extendInt='yes',y.0=y.0,y.1=y.1)$root
    ## optimize(objective,c(.1,2),y.0=y.0,y.1=y.1,maximum=TRUE)$maximum
})
hist(shape.mles)
abline(v=shape,col=2)

## auc MLE
m <- n <- 50
shape <- 1
auc.mles <- replicate(50, {
    y <- gen.data(n,n,auc)
    y.0 <- y$y.0; y.1 <- y$y.1
    y <- unlist(y)
    ## objective <- function(shapes,y.0,y.1) {
    ##     m <- length(y.0); n <- length(y.1)
    ##     y <- c(y.0,y.1)
    ##     sapply(shapes, function(shape) 
    ##         ## (length(y)+sum(y^shape)) / (m*sum(y.0^shape*log(y.0))/sum(y.0^shape) +  n*sum(y.1^shape*log(y.1))/sum(y.1^shape)  )
    ##     (m+n)/shape + sum(log(y)) - m*sum(y.0^shape*log(y.0))/sum(y.0^shape) - n*sum(y.1^shape*log(y.1))/sum(y.1^shape)
    ##     )
    ## }
    ## ## objective <- Vectorize(objective,'shape')
    ## shape.mle <- uniroot(objective,c(0.1,2),extendInt='yes',y.0=y.0,y.1=y.1)$root
    ## auc.mle <- mean(y.1^shape.mle) / (mean(y.0^shape.mle) + mean(y.1^shape.mle))
    ## ## optimize(objective,c(.1,2),y.0=y.0,y.1=y.1,maximum=TRUE)$maximum
    auc.mle(y.0,y.1)
})
hist(auc.mles)
abline(v=auc,col=2)


