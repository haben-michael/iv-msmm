## TODO:
## DONE--diff rather than just auc.
## full oracle (f-test) rather than just derivative oracle.
## maybe bad power in addition to bad fpr.
## coefs.lda should take x,g not x.0,x.1
## different strength of derivative
## maybe bootstrap estimator

## refactor--encapsulate routine for estimating the proposed estimator

auc.index.linearize <- function(x,g,beta,infl.fn,deriv.fn,params) {
    x.0 <- x[g==0,]; x.1 <- x[g==1,] # clean up
    approx.1 <- auc.hajek(x.0%*%beta,x.1%*%beta,terms.only=TRUE,IID=TRUE)
    infl.hat <- infl.fn(x,g,params)
    deriv.hat <- deriv.fn(x,g,beta)#grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat)
    ## sapply(list(deriv.hat.1=deriv.hat.1,deriv.hat.2=deriv.hat.2,deriv.hat.3=deriv.hat.3), function(deriv.hat) {
    approx.2 <- deriv.hat%*%infl.hat
    ## obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc
    approx <- as.numeric(approx.1+approx.2)
    ## var(approx) / length(approx)
}


require(mvtnorm)
require(parallel)
require(numDeriv)
source('misc.R')
sim <- function(n,p,taylor.part,alpha=.05,reps=3e2) {
    ## browser()
    uu <- 4*p*pi*exp(1)*taylor.part^2
    pi.0.min <- 1/2*(1-uu+sqrt(uu*(uu+2)))
    stopifnot(pi.0.min>1/2)
    pi.0 <- runif(1,pi.0.min,1)
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
    ## n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
    stopifnot(n.0<n)
    ## z.stats <- replicate(5e2, {
    z.stats <- mclapply(1:reps, mc.cores=detectCores()-3, FUN= function(dd){
        x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
        x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
        x <- rbind(x.0,x.1)
        g <- c(rep(0,n.0),rep(1,n.1))
        beta.hat <- coefs.lda(x.0,x.1,params=list())
        params.hat <- list(mu.0=colMeans(x.0),mu.1=colMeans(x.1),pi.0=nrow(x.0)/(nrow(x.0)+nrow(x.1)),Sigma.0=var(x.0),Sigma.1=var(x.1))
        ## approx.1 <- auc.hajek.lda.gaussian(x,g,beta.hat,params.hat,terms.only=TRUE,IID=TRUE)
        ## approx.1 <- auc.hajek(x.0%*%beta.hat,x.1%*%beta.hat,terms.only=TRUE,IID=TRUE)
        ## infl <- infl.lda(x,g,params.hat,terms.only=TRUE)
        ## deriv.star.hat.1 <- rep(0,p)
        ## deriv.star.hat.2 <- deriv.star
        ## deriv.star.hat.3 <- grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat)
        ## sapply(list(deriv.star.hat.1=deriv.star.hat.1,deriv.star.hat.2=deriv.star.hat.2,deriv.star.hat.3=deriv.star.hat.3), function(deriv.star.hat) {
        ##     obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc
        ##     approx.2 <- deriv.star.hat%*%infl
        ##     approx <- as.numeric(approx.1+approx.2)
        ##     ## var.hat <- delong.var.adjusted(x,g,beta.hat,infl.lda,deriv.fn,params.hat)
        ##     obs / sqrt(var(approx) / length(approx))        
        ## })    
        deriv.fn.1 <- function(x,g,beta)rep(0,ncol(x))
        deriv.fn.2 <- function(x,g,beta)deriv.star
        deriv.fn.3 <- function(x,g,beta)grad(function(u)auc.hat(x[g==0,]%*%u,x[g==1,]%*%u),beta)
        sapply(list(deriv.fn.1=deriv.fn.1,deriv.fn.2=deriv.fn.2,deriv.fn.3=deriv.fn.3), function(deriv.fn) {
            obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc
            ## approx.2 <- deriv.fn%*%infl
            ## approx <- as.numeric(approx.1+approx.2)
            linearized <- auc.index.linearize(x,g,beta.hat,infl.lda,deriv.fn,params.hat)
            var.hat <- var(linearized) / length(linearized)
            obs / sqrt(var.hat)
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
legend('bottomleft',legend=c('delong','oracle','delong adjusted'),col=1:3,lty=1)
abline(h=1-alpha,lty=2)
## matplot(max.dist,pch=1,type='l',col=1:3,lty=1)
## legend('bottomleft',legend=c('ignored','oracle','estimated'),col=1:3,lty=1)
## abline(h=1-alpha,lty=2)
## save.image('sessions/15e.RData')








## refactor--encapsulate data generation for this adversarial
## model

# sampler for adversarial data. conditioning on n.0,n.1, since
# otherwise couldn't avoid sampling n.0==0.
sampler.init <- function(n,p,taylor.part) {
    uu <- 4*p*pi*exp(1)*taylor.part^2
    pi.0.min <- 1/2*(1-uu+sqrt(uu*(uu+2)))
    stopifnot(pi.0.min>1/2)
    pi.0 <- runif(1,pi.0.min,1)
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
    ## n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
    stopifnot(n.0<n)
    sample <- function() {
        x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
        x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
        x <- rbind(x.0,x.1)
        d <- c(rep(0,n.0),rep(1,n.1))
        return(list(x=x,d=d))
    }
    return(list(n.0=n.0,n.1=n.1,deriv.star=deriv.star,auc=auc,params=params,sample=sample))
}
auc.index.linearize <- function(x,d,beta,infl.fn,deriv.fn=NULL) {
    x.0 <- x[d==0,]; x.1 <- x[d==1,] # clean up
    if(is.null(deriv.fn)) deriv.fn <- function(x,d,beta)grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta)
    approx.1 <- auc.hajek(x.0%*%beta,x.1%*%beta,terms.only=TRUE,IID=TRUE)
    infl.hat <- infl.fn(x,d)
    deriv.hat <- deriv.fn(x,d,beta)#grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat)
    ## sapply(list(deriv.hat.1=deriv.hat.1,deriv.hat.2=deriv.hat.2,deriv.hat.3=deriv.hat.3), function(deriv.hat) {
    approx.2 <- deriv.hat%*%infl.hat
    ## obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc
    approx <- as.numeric(approx.1+approx.2)
    ## var(approx) / length(approx)
}
require(mvtnorm)
require(parallel)
require(numDeriv)
source('../../misc.R')
sim <- function(n,p,taylor.part,alpha=.05,reps=3e2) {
    ## browser()
    cat('.')
    sampler <- sampler.init(n,p,taylor.part)
    ## n.0 <- sampler$n.0; n.1 <- sampler$n.1
    ## deriv.star <- sampler$deriv.star
    ## auc <- sampler$auc
    ## params <- sampler$params
    ## sample <- sampler$sample
    ## z.stats <- replicate(5e2, {
    z.stats <- mclapply(1:reps, mc.cores=detectCores()-3, FUN= function(dd){
        data <- sampler$sample()
        x <- data$x; d <- data$d
        x.0 <- x[d==0,]; x.1 <- x[d==1,]
        ## print(dim(x.1))
        ## params.hat <- data$params.hat
        ## beta.hat <- data$beta.hat
        beta.hat <- coefs.lda(x.0,x.1)
        ## params.hat <- list(mu.0=colMeans(x.0),mu.1=colMeans(x.1),pi.0=nrow(x.0)/(nrow(x.0)+nrow(x.1)),Sigma.0=var(x.0),Sigma.1=var(x.1))
        deriv.fn.1 <- function(x,d,beta)rep(0,ncol(x))
        deriv.fn.2 <- function(x,d,beta)sampler$deriv.star
        deriv.fn.3 <- NULL#function(x,d,beta)grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta)
        sapply(list(delong=deriv.fn.1,oracle=deriv.fn.2,proposed=deriv.fn.3), function(deriv.fn) {
            obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - sampler$auc
            ## approx.2 <- deriv.fn%*%infl
            ## approx <- as.numeric(approx.1+approx.2)
            linearized <- auc.index.linearize(x,d,beta.hat,infl.fn=function(x,d)infl.lda(x,d,params=NULL,var.equal=FALSE),deriv.fn=deriv.fn)
            var.hat <- var(linearized) / length(linearized)
            obs / sqrt(var.hat)
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
## options(error=browser)
start <- Sys.time()
set.seed(1)
p <- 8
## n <- 5e3
taylor.part <- 1/10
ns <- round(seq(3e2,5e3,len=30))
## ns <- round(seq(3e2,3e3,len=10))
reps <- 1e2#1e3
alpha <- .05
sim0 <- lapply(ns, function(n) sim(n,p,taylor.part,alpha=alpha,reps=reps))
fpr <- do.call(rbind,lapply(sim0,function(m)m[,'fpr']))
max.dist <- do.call(rbind,lapply(sim0,function(m)m[,'max.dist']))
Sys.time() - start
matplot(x=ns,fpr,pch=1,type='l',lty=1:3,col=1)
## legend('bottomleft',legend=c('delong','oracle','delong adjusted'),col=1:3,lty=1)
legend('bottomleft',legend=colnames(fpr),lty=1:3)
abline(h=1-alpha,lty=2)










## delta auc

sampler.init <- function(n,p,adjust.size) {
    u <- 4*p*pi*exp(1)*adjust.size^2
    pi.0.min <- 1/2*(1-u+sqrt(u*(u+2)))
    stopifnot(pi.0.min>1/2)
    pi.0 <- runif(1,pi.0.min,1)
    pi.1 <- 1-pi.0
    d <- sqrt(2/p)
    a <- (1-2*pi.0)^2/16/p/pi/exp(1)
    b <- -2*a-pi.0*adjust.size^2
    c <- a+(pi.0-1)*adjust.size^2
    eps <- min(1/2/a*(-b+c(-1,1)*sqrt(b^2-4*a*c)))
    stopifnot(eps>0)
    beta.star <- rep(1,p )*d
    Sigma.1 <- diag(p) 
    Sigma.0 <- diag(c(rep(eps,p/2),rep(2-eps,p/2)) )
    Sigma.pi <- pi.0*Sigma.0+pi.1*Sigma.1
    Sigma <- Sigma.0+Sigma.1
    mu.0 <- rep(0,p)
    mu.1 <- (pi.0*Sigma.0+pi.1*Sigma.1)%*%beta.star
    params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0,beta=beta.star,pi.0=pi.0)
    params$deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
    params$auc <- auc.lda.gaussian(beta.star,params)
    ## n <- 5e3
    n.0 <- round(n*pi.0); n.1 <- n-n.0
    ## n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
    stopifnot(n.0<n)
    sample <- function() {
        x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
        x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
        x <- rbind(x.0,x.1)
        d <- c(rep(0,n.0),rep(1,n.1))
        return(list(x=x,d=d))
    }
    return(list(n.0=n.0,n.1=n.1,params=params,sample=sample))
}
auc.index.linearize <- function(x,d,beta,infl.fn,deriv.fn=NULL) {
    n <- nrow(x)
    x.0 <- x[d==0,]; x.1 <- x[d==1,] # clean up
    if(is.null(deriv.fn)) deriv.fn <- function(x,d,beta)numDeriv::grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta,method='simple',method.args=list(eps=1/n^(.5)))
    approx.1 <- auc.hajek(x.0%*%beta,x.1%*%beta,terms.only=TRUE,IID=TRUE)
    infl.hat <- infl.fn(x,d)
    deriv.hat <- deriv.fn(x,d,beta)#grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat)
    ## sapply(list(deriv.hat.1=deriv.hat.1,deriv.hat.2=deriv.hat.2,deriv.hat.3=deriv.hat.3), function(deriv.hat) {
    approx.2 <- deriv.hat%*%infl.hat
    ## obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc
    approx <- as.numeric(approx.1+approx.2)
    ## var(approx) / length(approx)
}
require(mvtnorm)
require(parallel)
require(numDeriv)
source('../../misc.R')
sim <- function(n,p.full,p.red,adjust.size,alpha=.05,reps=3e2) {
    ## browser()
    cat('.')
    sampler <- sampler.init(n,p.full,adjust.size)
    params.full <- sampler$params
    params.red <- with(params.full, list(mu.0=mu.0[1:p.red],mu.1=mu.1[1:p.red],Sigma.0=Sigma.0[1:p.red,1:p.red],Sigma.1=Sigma.1[1:p.red,1:p.red],pi.0=pi.0))
    params.red$beta <- with(params.red,solve(pi.0*Sigma.0+(1-pi.0)*Sigma.1)%*%(mu.1-mu.0)) 
    params.red$auc <- auc.lda.gaussian(params.red$beta,params.red)
    params.red$deriv.star <- auc.deriv.lda.gaussian(params.red$beta,params.red) # ==0 for p.red=p.full=2 
    ## z.stats <- replicate(5e2, {
    z.stats <- mclapply(1:reps, mc.cores=detectCores()-3, FUN= function(dd){
        xd <- sampler$sample()
        ## x <- xd$x; d <- xd$d
        params.full$x <- xd$x; params.full$d <- xd$d # maybe rename to "dataset"?
        params.red$x <- xd$x[,1:p.red]; params.red$d <- xd$d
        ## datasets <- list(full=list(x=x,d=d,auc=params.full$auc,deriv.star=params.full$deriv.star),
                         ## reduced=list(x=x[,1:p.red],d=d,auc=params.red$auc,deriv.star=params.red$deriv.star))
        out <- lapply(list(full=params.full,red=params.red), function(params) {
            with(params, {
                x.0 <- x[d==0,]; x.1 <- x[d==1,]
                beta.hat <- coefs.lda(x.0,x.1)
                obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc
                linearized <- lapply(list(delong= function(x,d,beta)rep(0,ncol(x)),oracle=function(x,d,beta)deriv.star, proposed=NULL), function(deriv.fn)
                    auc.index.linearize(x,d,beta.hat,infl.fn=function(x,d)infl.lda(x,d,params=NULL,var.equal=FALSE),deriv.fn=deriv.fn))
                linearized <- simplify2array(linearized)
                list(obs=obs,linearized=linearized)
            })
        })
        diff.linearized <- out$full$linearized - out$red$linearized
        diff.obs <- out$full$obs - out$red$obs
        ## diff.linearized <-  out$red$linearized
        ## diff.obs <-  out$red$obs
        diff.var.hat <- apply(diff.linearized,2,function(iid)var(iid) / length(iid))
        c(diff.obs) / sqrt(diff.var.hat)
    })
    z.stats <- simplify2array(z.stats)
    fpr <- rowMeans(abs(z.stats)<qnorm(1-alpha/2))
}
## options(error=browser)
start <- Sys.time()
set.seed(2)
p.full <- 8
p.red <- 6
## n <- 5e3
adjust.size <- 1/8
ns <- round(seq(3e2,5e3,len=30))
## ns <- round(seq(3e2,5e3,len=10))
reps <- 3e2#1e3
alpha <- .05
fpr <- sapply(ns, function(n) sim(n,p.full,p.red,adjust.size,alpha=alpha,reps=reps))
Sys.time() - start
matplot(x=ns,t(fpr),pch=1,type='l',lty=1:3,col=1)
legend('bottomleft',legend=rownames(fpr),lty=1:3)
abline(h=1-alpha,lty=2)









## 23-3-3 adding parameter to control AUC

sampler.init <- function(n,p,auc,epsilon=NULL,pi.0=NULL,adjust.size=NULL) {
    stopifnot((!is.null(epsilon) && !is.null(pi.0)) || !is.null(adjust.size))
    if(is.null(epsilon) || is.null(pi.0)) { ## this branch may need updating
        u <- 4*p*pi*exp(1)*adjust.size^2
        pi.0.min <- 1/2*(1-u+sqrt(u*(u+2)))
        stopifnot(pi.0.min>1/2)
        ## pi.0 <- runif(1,pi.0.min,1)
        pi.0 <- pi.0.min*.9  + .1*1
        stopifnot(auc>.5 && auc < 1)
        a <- (1-2*pi.0)^2/16/p/pi/exp(1)
        b <- -2*a-pi.0*adjust.size^2
        c <- a+(pi.0-1)*adjust.size^2
        epsilon <- min(1/2/a*(-b+c(-1,1)*sqrt(b^2-4*a*c)))
        stopifnot(epsilon>0)
    } else {
        if(is.null(adjust.size)) {
            adjust.size <- (1-2*pi.0)*(1-epsilon)/4/sqrt(p.full*pi*exp(1))/sqrt(1+(epsilon-1)*pi.0)
        } else stop() ## can remove initial stopifnot
    }
    pi.1 <- 1-pi.0
    d <- uniroot(function(x)pnorm(x*sqrt(p/2)) - auc,interval=c(0,1))$root
    ## d <- sqrt(2/p)
    beta.star <- rep(1,p )*d
    Sigma.1 <- diag(p) 
    Sigma.0 <- diag(c(rep(epsilon,p/2),rep(2-epsilon,p/2)) )
    Sigma.pi <- pi.0*Sigma.0+pi.1*Sigma.1
    Sigma <- Sigma.0+Sigma.1
    mu.0 <- rep(0,p)
    mu.1 <- (pi.0*Sigma.0+pi.1*Sigma.1)%*%beta.star
    params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0,beta=beta.star,pi.0=pi.0,adjust.size=adjust.size,epsilon=epsilon)
    params$deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
    params$auc <- auc.lda.gaussian(beta.star,params)
    ## n <- 5e3
    n.0 <- round(n*pi.0); n.1 <- n-n.0
    ## n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
    stopifnot(n.0<n)
    sample <- function() {
        x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
        x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
        x <- rbind(x.0,x.1)
        d <- c(rep(0,n.0),rep(1,n.1))
        return(list(x=x,d=d))
    }
    return(list(n.0=n.0,n.1=n.1,params=params,sample=sample))
}
auc.index.linearize <- function(x,d,beta,infl.fn,deriv.fn=NULL) {
    n <- nrow(x)
    x.0 <- x[d==0,]; x.1 <- x[d==1,] # clean up
    if(is.null(deriv.fn)) deriv.fn <- function(x,d,beta)numDeriv::grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta,method='simple',method.args=list(eps=1/n^(.5)))
    approx.1 <- auc.hajek(x.0%*%beta,x.1%*%beta,terms.only=TRUE,IID=TRUE)
    infl.hat <- infl.fn(x,d)
    deriv.hat <- deriv.fn(x,d,beta)#grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat)
    ## sapply(list(deriv.hat.1=deriv.hat.1,deriv.hat.2=deriv.hat.2,deriv.hat.3=deriv.hat.3), function(deriv.hat) {
    approx.2 <- deriv.hat%*%infl.hat
    ## obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc
    approx <- as.numeric(approx.1+approx.2)
    ## var(approx) / length(approx)
}
suppressPackageStartupMessages({
    require(mvtnorm)
    require(parallel)
    require(numDeriv)
})
source('../../misc.R')
sim <- function(n,p.full,p.red,auc=NULL,epsilon=NULL,pi.0=NULL,adjust.size=NULL,alpha=.05,reps=3e2) {
    ## browser()
    ## cat('.')
    sampler <- sampler.init(n,p.full,auc=auc,epsilon=epsilon,pi.0=pi.0,adjust.size=adjust.size)
    params.full <- sampler$params
    print(params.full$adjust.size)
    ## params.red <- with(params.full, list(mu.0=mu.0[1:p.red],mu.1=mu.1[1:p.red],Sigma.0=Sigma.0[1:p.red,1:p.red],Sigma.1=Sigma.1[1:p.red,1:p.red],pi.0=pi.0))
    ## params.red$beta <- with(params.red,solve(pi.0*Sigma.0+(1-pi.0)*Sigma.1)%*%(mu.1-mu.0)) 
    ## params.red$auc <- auc.lda.gaussian(params.red$beta,params.red)
    ## params.red$deriv.star <- auc.deriv.lda.gaussian(params.red$beta,params.red) # ==0 for p.red=p.full=2 
    terms <- replicate(reps, {
    ## z.stats <- mclapply(1:reps, mc.cores=detectCores()-3, FUN= function(dd){
        xd <- sampler$sample()
        params.full$x <- xd$x; params.full$d <- xd$d # maybe rename to "dataset"? these arent just params anymore
        ## params.red$x <- xd$x[,1:p.red]; params.red$d <- xd$d
        ## out <- lapply(list(full=params.full,red=params.red), function(params) {
        out <- lapply(list(full=params.full), function(params) {
            with(params, {
                x.0 <- x[d==0,]; x.1 <- x[d==1,]
                beta.hat <- coefs.lda(x.0,x.1)
                obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) #- auc
                linearized <- lapply(list(delong= function(x,d,beta)rep(0,ncol(x)),oracle=function(x,d,beta)deriv.star, proposed=NULL), function(deriv.fn)
                    auc.index.linearize(x,d,beta.hat,infl.fn=function(x,d)infl.lda(x,d,params=NULL,var.equal=FALSE),deriv.fn=deriv.fn))
                linearized <- simplify2array(linearized)
                list(obs=obs,true=auc,linearized=linearized)
            })
        })
        linearized <- out$full$linearized# - out$red$linearized
        obs <- out$full$obs #- out$red$obs
        true <- out$full$true #- out$red$true
        var.hat <- apply(linearized,2,function(iid)var(iid) / length(iid))
        ## c(obs) / sqrt(var.hat)
        c(obs=obs,true=true,var.hat)
    })
}
## options(error=browser)
start <- Sys.time()
## set.seed(1)
p.full <- 8
p.red <- 6
n <- 5e2
## adjust.size <- 1/12
pi.0 <- .97
epsilon <- .0
ns <- round(seq(3e2,5e3,len=30))
ns <- round(seq(3e2,5e3,len=10))
reps <- 1e3
alpha <- .05
auc <- .8
 terms <- ## tryCatch(
    sim(n,p.full,p.red,auc=auc,epsilon=epsilon,pi.0=pi.0,adjust.size=NULL,alpha=alpha,reps=reps)
    ## error=function(e)NA)
Sys.time() - start
z.stats <- sapply(c('delong','oracle','proposed'), function(var.hat)(terms['obs',] - terms['true',]) / sqrt(terms[var.hat,]))
coverage <- colMeans(abs(z.stats) < qnorm(1-alpha/2))
std.error <- apply(z.stats,2,sd)
coverage







## delong performance is too good--debugging--looking at parts of
## linearization. [result: seemed to be a function of sample
## size. increasing sims from around 2k to 10k.]

suppressPackageStartupMessages({
    require(mvtnorm)
    require(parallel)
    require(numDeriv)
})
source('../../misc.R')
source('misc.R')
## sim <- function(n,p.full,p.red,auc=NULL,epsilon=NULL,pi.0=NULL,adjust.size=NULL,alpha=.05,reps=3e2) {
    ## browser()
    ## cat('.')
start <- Sys.time()
## set.seed(1)
p.full <- 8
p.red <- 6
n <- 5e2
## adjust.size <- 1/12
pi.0 <- .5
epsilon <- 1
ns <- round(seq(3e2,5e3,len=30))
ns <- round(seq(3e2,5e3,len=10))
reps <- 1e3
alpha <- .05
auc <- .9
sampler <- sampler.init(n,p.full,auc=auc,epsilon=epsilon,pi.0=pi.0,adjust.size=adjust.size)
params <- params.full <- sampler$params
## try <- with(sampler$params, {
##     Sigma.pi <- pi.0*Sigma.0+(1-pi.0)*Sigma.1
##     deriv.star%*%solve(expm::sqrtm(Sigma.pi))
## })
## try
## (1-2*pi.0)*(1-epsilon)/4/sqrt(p.full*pi*exp(1))/sqrt(1+(epsilon-1)*pi.0)
xd <- sampler$sample()
x <- xd$x; d <- xd$d
x.0 <- x[d==0,]; x.1 <- x[d==1,]
beta.hat <- coefs.lda(x.0,x.1)
beta <- params.full$beta
## obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) #- auc
## linearized <- lapply(list(delong= function(x,d,beta)rep(0,ncol(x)),oracle=function(x,d,beta)deriv.star, proposed=NULL), function(deriv.fn)
## auc.index.linearize(x,d,beta.hat,infl.fn=function(x,d)infl.lda(x,d,params=NULL,var.equal=FALSE),deriv.fn=deriv.fn))
approx.1 <- auc.hajek(x.0%*%beta.hat,x.1%*%beta.hat,terms.only=TRUE,IID=TRUE)
infl.hat <- infl.lda(x,d,params=NULL,var.equal=FALSE)
deriv.hat <- params.full$deriv.star#deriv.fn(x,d,beta)#grad(function(u)auc.hat(x.0%*%u,x.1%*%u),beta.hat)
## sapply(list(deriv.hat.1=deriv.hat.1,deriv.hat.2=deriv.hat.2,deriv.hat.3=deriv.hat.3), function(deriv.hat) {
approx.2 <- as.numeric(deriv.hat%*%infl.hat)
## obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc
approx <- as.numeric(approx.1+approx.2)


suppressPackageStartupMessages({
    require(mvtnorm)
    require(parallel)
    require(numDeriv)
})
set.seed(1)
source('../../misc.R')
source('misc.R')
p <- 8
## p.red <- 6
n <- 5e2
## adjust.size <- 1/12
pi.0 <- .9
epsilon <- .0
## ns <- round(seq(3e2,5e3,len=30))
ns <- round(seq(3e2,5e3,len=20))
reps <- 1e3
alpha <- .05
auc <- .9
by.n <- sapply(ns, simplify=FALSE, FUN=function(n) {
    sampler <- sampler.init(n,p,auc=auc,epsilon=epsilon,pi.0=pi.0)
    params <- sampler$params
    beta.star <- params$beta
    deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
    parts <- replicate(1e2, {
        xd <- sampler$sample()
        x <- xd$x; d <- xd$d
        x.0 <- x[d==0,]; x.1 <- x[d==1,]
        beta.hat <- coefs.lda(x.0,x.1,params=params)
        ## obs.1 <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc.lda.gaussian(beta.hat,params=params)
        ## approx.1 <- auc.hajek.lda.gaussian(x,d,beta.hat,params,terms.only=FALSE)
        ## ## approx.1 <- auc.hajek(x.0%*%beta.hat,x.1%*%beta.hat,terms.only=FALSE)          obs.1 - approx.1
        obs.2 <- auc.lda.gaussian(beta.hat,params) - auc.lda.gaussian(beta.star,params)
        infl <- infl.lda(x,d,params,var.equal=FALSE,terms.only=FALSE)
        ## infl <- infl.lda(x,d,params,var.equal=FALSE,terms.only=TRUE)
        ## approx.2 <- deriv.star%*%(beta.hat - beta.star)
        approx.2 <- deriv.star%*%infl
        obs.2 - approx.2
        ## obs.1+obs.2 - (approx.1+approx.2)
        ## ## beta.hat-beta.star
        ## approx.1 <- auc.hajek.lda.gaussian(x,d,beta.hat,params,terms.only=TRUE,IID=TRUE)
        ## infl <- infl.lda(x,d,params,var.equal=FALSE,terms.only=TRUE)
        ## approx.2 <- deriv.star%*%infl
        ## c(approx.1.sd=sd(approx.1),approx.2.sd=sd(approx.2))
    })
})
by.n <- simplify2array(by.n)
mad <- apply(abs(by.n),2,mean)
plot(mad)
## sds <- apply(by.n-true,2,sd)
## plot(sds)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)


## op <- par(mfrow=c(1,dim(by.n)[1]))
## for(i in 1:dim(by.n)[1]) {
##     stat <- colMeans(by.n[i,,])
##     lm0 <- lm(log(stat)~log(ns))
##     print(coef(lm0))
##     plot(ns,stat,main=rownames(by.n)[i])
##     curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)
## }
## par(op)


suppressPackageStartupMessages({
    require(mvtnorm)
    require(parallel)
    require(numDeriv)
})
set.seed(1)
start <- Sys.time()
source('../../misc.R')
source('misc.R')
p <- 8
## p.red <- 6
n <- 1e3
## adjust.size <- 1/12
pi.0 <- .8
epsilon <- .01
## ns <- round(seq(3e2,5e3,len=30))
## ns <- round(seq(3e2,3e3,len=10))
reps <- 1e3
## alpha <- .05
auc <- .9
## by.n <- sapply(ns, function(n) {
sampler <- sampler.init(n,p,auc=auc,epsilon=epsilon,pi.0=pi.0)
params <- sampler$params
beta.star <- params$beta
deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
z.stats <- replicate(3e2, {
    xd <- sampler$sample()
    x <- xd$x; d <- xd$d
    x.0 <- x[d==0,]; x.1 <- x[d==1,]
    beta.hat <- coefs.lda(x.0,x.1,params=params)
    ## obs.1 <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc.lda.gaussian(beta.hat,params=params)
    ## approx.1 <- auc.hajek.lda.gaussian(x,d,beta.hat,params,terms.only=FALSE)
    ## ## approx.1 <- auc.hajek(x.0%*%beta.hat,x.1%*%beta.hat,terms.only=FALSE)          obs.1 - approx.1
    ## obs.2 <- auc.lda.gaussian(beta.hat,params) - auc.lda.gaussian(beta.star,params)
    ## infl <- infl.lda(x,d,params,var.equal=FALSE,terms.only=FALSE)
    ## ## infl <- infl.lda(x,d,params,var.equal=FALSE,terms.only=TRUE)
    ## ## approx.2 <- deriv.star%*%(beta.hat - beta.star)
    ## approx.2 <- deriv.star%*%infl
    ## obs.2 - approx.2
    ## obs.1+obs.2 - (approx.1+approx.2)
    obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - params$auc
    approx.1 <- auc.hajek.lda.gaussian(x,d,beta.hat,params,terms.only=TRUE,IID=TRUE)
    infl <- infl.lda(x,d,params,var.equal=FALSE,terms.only=TRUE)
    approx.2 <- as.numeric(deriv.star%*%infl)
    approx <- approx.1+approx.2
    obs / sqrt(var(approx) / length(approx))
## mean(abs(beta.hat-beta.star))
    ## obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc.lda.gaussian(beta.star,params)
## })
})
Sys.time() - start
op <- par(mfrow=c(1,2))
qqnorm(z.stats); abline(0,1)
hist(z.stats)
par(op)
sd(z.stats)
median(z.stats)



##  behavior on reduced data
suppressPackageStartupMessages({
    require(mvtnorm)
    require(parallel)
    require(numDeriv)
})
set.seed(1)
source('../../misc.R')
source('misc.R')
p <- 8
p.red <- 6
## p.red <- 6
n <- 5e2
## adjust.size <- 1/12
pi.0 <- .9
epsilon <- .0
## ns <- round(seq(3e2,5e3,len=30))
ns <- round(seq(3e2,3e3,len=20))
reps <- 1e3
alpha <- .05
auc <- .9
by.n <- sapply(ns, simplify=FALSE, FUN=function(n) {
    sampler <- sampler.init(n,p,auc=auc,epsilon=epsilon,pi.0=pi.0)
    ## params <- sampler$params
    ## params$deriv <- auc.deriv.lda.gaussian(params$beta,params)
    params <- with(sampler$params, list(mu.0=mu.0[1:p.red],mu.1=mu.1[1:p.red],Sigma.0=Sigma.0[1:p.red,1:p.red],Sigma.1=Sigma.1[1:p.red,1:p.red],pi.0=pi.0))
    params$beta <- with(params,solve(pi.0*Sigma.0+(1-pi.0)*Sigma.1)%*%(mu.1-mu.0))
    params$auc <- auc.lda.gaussian(params$beta,params)
    params$deriv <- auc.deriv.lda.gaussian(params$beta,params)
    parts <- replicate(1e2, {
        xd <- sampler$sample()
        x <- xd$x[,1:p.red]; d <- xd$d
        x.0 <- x[d==0,]; x.1 <- x[d==1,]
        beta.hat <- coefs.lda(x.0,x.1,params=params)
        ## obs.1 <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) - auc.lda.gaussian(beta.hat,params=params)
        ## approx.1 <- auc.hajek.lda.gaussian(x,d,beta.hat,params,terms.only=FALSE)
        ## ## ## approx.1 <- auc.hajek(x.0%*%beta.hat,x.1%*%beta.hat,terms.only=FALSE)
        ## obs.1 - approx.1
        obs.2 <- auc.lda.gaussian(beta.hat,params) - auc.lda.gaussian(params$beta,params)
        infl <- infl.lda(x,d,params,var.equal=FALSE,terms.only=FALSE)
        ## infl <- infl.lda(x,d,params,var.equal=FALSE,terms.only=TRUE)
        ## approx.2 <- deriv.star%*%(beta.hat - beta.star)
        approx.2 <- params$deriv%*%infl
        obs.2 - approx.2
        ## obs.1+obs.2 - (approx.1+approx.2)
        ## ## beta.hat-beta.star
        ## approx.1 <- auc.hajek.lda.gaussian(x,d,beta.hat,params,terms.only=TRUE,IID=TRUE)
        ## infl <- infl.lda(x,d,params,var.equal=FALSE,terms.only=TRUE)
        ## approx.2 <- deriv.star%*%infl
        ## c(approx.1.sd=sd(approx.1),approx.2.sd=sd(approx.2))
    })
})
by.n <- simplify2array(by.n)
mad <- apply(abs(by.n),2,mean)
plot(mad)
## sds <- apply(by.n-true,2,sd)
## plot(sds)
lm0 <- lm(log(mad)~log(ns))
coef(lm0)



op <- par(mfrow=c(1,dim(by.n)[1]))
for(i in 1:dim(by.n)[1]) {
    stat <- colMeans(by.n[i,,])
    lm0 <- lm(log(stat)~log(ns))
    print(coef(lm0))
    plot(ns,stat,main=rownames(by.n)[i])
    curve(exp(coef(lm0)[1]+coef(lm0)[2]*log(x)),add=TRUE)
}
par(op)



## bootstrap


suppressPackageStartupMessages({
    require(mvtnorm)
    require(parallel)
    require(numDeriv)
})
set.seed(1)
source('../../misc.R')
source('misc.R')
set.seed(1)
start <- Sys.time()
source('../../misc.R')
source('misc.R')
p <- 8
## p.red <- 6
n <- 1e3
## adjust.size <- 1/12
pi.0 <- .8
epsilon <- .01
## ns <- round(seq(3e2,5e3,len=30))
## ns <- round(seq(3e2,3e3,len=10))
reps <- 1e3
alpha <- .05
auc <- .9
## by.n <- sapply(ns, function(n) {
sampler <- sampler.init(n,p,auc=auc,epsilon=epsilon,pi.0=pi.0)
params <- sampler$params
beta.star <- params$beta
deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
CIs <- replicate(3e2, {
    xd <- sampler$sample()
    x <- xd$x; d <- xd$d
    x.0 <- x[d==0,]; x.1 <- x[d==1,]
    beta.hat <- coefs.lda(x.0,x.1)
    obs <-  auc.hat(x.0%*%beta.hat,x.1%*%beta.hat)
    auc.bs <- replicate(1e2, {
        idx <- sample(1:nrow(x),replace=TRUE)
        x.bs <- x[idx,]; d.bs <- d[idx]
        x.0.bs <- x[d.bs==0,]; x.1.bs <- x.bs[d.bs==1,]
        beta.hat.bs <- coefs.lda(x.0.bs,x.1.bs)
        auc.hat(x.0.bs%*%beta.hat.bs,x.1.bs%*%beta.hat.bs)
    })
    obs + c(-1,1)*qnorm(1-alpha/2)*sd(auc.bs)
})
Sys.time() - start
mean(apply(CIs-auc,2,prod)<0)
