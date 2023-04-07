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
        stopifnot(epsilon<=1 && epsilon >=0)
        if(is.null(adjust.size)) {
            adjust.size <- (1-2*pi.0)*(1-epsilon)/4/sqrt(p*pi*exp(1))/sqrt(1+(epsilon-1)*pi.0)
        } else stop() ## can remove initial stopifnot
    }
    pi.1 <- 1-pi.0
    d <- qnorm(auc)/sqrt(p/2)
    ## d <- sqrt(2/p)
    beta.star <- rep(1,p )*d
    Sigma.1 <- diag(p) 
    Sigma.0 <- diag(c(rep(epsilon,p/2),rep(2-epsilon,p/2)) )
    Sigma.pi <- pi.0*Sigma.0+pi.1*Sigma.1
    Sigma <- Sigma.0+Sigma.1
    mu.0 <- rep(0,p)
    mu.1 <- as.numeric((pi.0*Sigma.0+pi.1*Sigma.1)%*%beta.star)
    params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0,beta=beta.star,auc=auc,pi.0=pi.0,adjust.size=adjust.size,epsilon=epsilon)
    params$deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
    ## params$auc <- auc.lda.gaussian(beta.star,params)
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
