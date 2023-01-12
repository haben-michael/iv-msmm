egger.test <- function(y,v,robust=FALSE) {
    lm0 <- if(robust) {
               estimatr::lm_robust(I(y/sqrt(v)) ~ I(1/sqrt(v)))
           } else {
               lm(I(y/sqrt(v)) ~ I(1/sqrt(v)))
           }
    structure(unname(coef(summary(lm0))[1,c(1,4)]), names=c('stat','pval'))
}
## checked against metafor in #13a
begg.test <- function(y,v,method='kendall',...) {
    theta.fe <- sum(y/v)/sum(1/v)
    with(cor.test((y-theta.fe)/sqrt(v-1/sum(1/v)),v,method=method,...),
         structure(unname(c(stat=estimate,pval=p.value)),names=c('stat','pval')))    
}
lin.test <- function(y,v) {
    with(altmeta::metapb(y=y,s2=v,model='FE'),
         structure(unname(c(stat=skewness,pval=skewness.pval)),names=c('stat','pval')))             
}
f.test <- function(y,v) {
    n <- length(y)
    theta.fe <- sum(y/v)/sum(1/v)
    regressand <- (y-theta.fe)/sqrt(v)
    lm0 <- lm(regressand ~ I(1/sqrt(v)))
    f.stat <- ((sum(regressand^2)-sum(resid(lm0)^2)) / 2) / (sum(resid(lm0)^2) / (n-2))
    pval <- 1-pf(f.stat,2,n-2)
    structure(c(unname(coef(lm0)),pval), names=c('stat.egger','stat.begg','pval'))
}
rlin <- function(n,mean=0,keep.prob=1,rsigma=function(n)runif(n,1,2),var.between=0,method,alpha.cutoff=.05) {
    y  <- sigma <- numeric()
    tries <- 0
    while(length(y)<n) {
        tries <- tries+1
        sigma.try <- rsigma(1)
        ## sigma.try <- rchisq(1,df=1)
        ## sigma.try <- rbeta(1,1/2,1/2)
        y.try <- rnorm(1,mean,sd=sqrt(var.between+sigma.try^2))
        cutoff <- sigma.try*qnorm(1-alpha.cutoff)
        if(switch(method,
                  '1' = y.try>cutoff || rbinom(1,1,keep.prob),
                  '2' = y.try>cutoff || sigma.try<=1.5 || rbinom(1,1,keep.prob),
                  '3a' = y.try>cutoff/sigma.try || rbinom(1,1,keep.prob)
                  ))
        {
            y <- c(y,y.try); sigma <- c(sigma,sigma.try)
        }
    }
    ## print(tries)
    return(list(y=y,sigma=sigma))
}





rnorm.trunc <- function(n,mean,sd,lower) {u <- runif(n); mean+sd*qnorm(u+pnorm((lower-mean)/sd)*(1-u))}
dnorm.trunc <- function(x,mean,sd,lower) (x>=lower) * 1/sd * dnorm((x-mean)/sd) / (1-pnorm((lower-mean)/sd))
pnorm.trunc <- function(q,mean,sd,lower) (q>=lower)*(pnorm((q-mean)/sd) - pnorm((lower-mean)/sd)) / (1-pnorm((lower-mean)/sd))

rnorm.effect.cutoff <- function(n,grand.mean=0,rvar.within=runif,var.between=0,effect.cutoff=-Inf) {
        var.within <- rvar.within(n)
        y <- rnorm.trunc(n,mean=grand.mean,sd=sqrt(var.within+var.between),lower=effect.cutoff) 
        ## est.v.bw(y,v.wi)
        return(list(y=y,v=var.within))
}


## usual method of moments 'tau^2' estimator
est.var.between <- function(y,var.within) {
    n <- length(y)
    theta.fe <- sum(y/var.within/sum(1/var.within))
    Q <- sum((y-theta.fe)^2/var.within)
    v.bw.hat <- max(0, (Q-n+1)/(sum(1/var.within)-sum(1/var.within^2)/sum(1/var.within)))    
}







## formulas for uniform model at null (pubbias #5)
E.cond <- function(q,theta2,u3) {
    ## pos.idx <- theta2>0
    ## with(list(q=q[pos.idx],theta2=theta2[pos.idx],u3=u3[pos.idx]), {
    lim1 <- pmin(u3,pmax(0,-q/theta2))#pmax(0,pmin(u3,-q/theta2))
    lim2 <- pmax(lim1,pmin(u3,(1-q)/theta2))#pmin(u3,pmax(lim1,(1-q)/theta2))
    lim3 <- pmax(0,pmin(u3,-q/theta2))#pmin(u3,pmax(0,-q/theta2))
    lim4 <- pmin(lim3,pmax(0,(1-q)/theta2))#pmax(0,pmin(lim3,(1-q)/theta2))
    (lim2 - q*(lim2-lim1)-theta2/2*(lim2^2-lim1^2))   *  (theta2>0) + 
        (u3-lim4-q*(lim3-lim4)-theta2/2*(lim3^2-lim4^2))  *  (theta2<0) +
        (1-pmax(0,pmin(1,q)))*u3   *   (theta2==0)
}
E.diff.cond <- E.cond ## deprecate E.cond
E.diff <- function(q,theta.2)
    if (theta.2<1) {
        (1-1/2/theta.2*(1-q)^2)*(1>=q & q>=1-theta.2) + (q+theta.2/2)*(q<=1-theta.2 & q>=0) + (theta.2/2*(1+q/theta.2)^2)*(q<=0 & q>=-theta.2) + (q>=1)
    } else {
        (1-1/2/theta.2*(1-q)^2)*(1>=q & q>=0) + (1-1/2/theta.2+q/theta.2)*(q<=0 & q>=1-theta.2) + (theta.2/2*(1+q/theta.2)^2)*(q<=1-theta.2 & q>=-theta.2) + (q>=1)
    }
E.uncond  <- function(theta2) { # P((T1-theta2*T2)*T2<0)
    L0 <- pmin(1/theta2,1)
    L2 <- pmin(1,-1/theta2)
    ## 2 * ((L2/2+1/2*(theta2-1/2)*L2^2+1/3*(theta2^2/2-theta2)*L2^3-theta2^2/8*L2^4) * (theta2<0) +   (1/2-L0/2+(theta2/2+1/4)*L0^2-(theta2^2/6+theta2/3)*L0^3+1/8*theta2^2*L0^4)*(theta2>0) + (1/4)*(theta2==0))
    2*ifelse(theta2<0, (L2/2+1/2*(theta2-1/2)*L2^2+1/3*(theta2^2/2-theta2)*L2^3-theta2^2/8*L2^4), ifelse(theta2>0,(1/2-L0/2+(theta2/2+1/4)*L0^2-(theta2^2/6+theta2/3)*L0^3+1/8*theta2^2*L0^4),1/4))
    ## 2/theta^2 * ( (theta>0)*(theta^2/2-1/2*pmin(theta,theta^2)+(theta/2+1/4)*pmin(1,theta^2)-(theta/6+1/3)*pmin(1,theta^3)+1/8*pmin(1,theta^4)) + (theta<0)*( 1/2*pmin(theta^2,abs(theta)) + 1/2*(theta-1/2)*pmin(theta^2,1) + 1/3*(1+abs(theta)/2)*pmin(abs(theta)^3,theta^2)-1/8*pmin(1,theta^4)) + (theta==
}
mu0.uncond <- function(theta2) 2*E.uncond(theta2) - 1
## mu0 <- function(z1,s1,theta.2)4*E.diff.cond(z1+1/2-theta.2*s1, theta.2, s1) - 2*E.diff.cond(z1+1/2-theta.2*s1, theta.2,1)-2*s1+1
## 2*P((z1-z2-theta.2*(s1-s2))*(s1-s2)<0 | z1,s1) - 1
mu0 <- function(z1,s1,theta.2,mean.s=1/2)4*E.diff.cond(z1+1/2-theta.2*(s1-(mean.s-1/2)), theta.2, s1-(mean.s-1/2)) - 2*E.diff.cond(z1+1/2-theta.2*(s1-mean.s+1/2), theta.2,1)-2*(s1-mean.s+1/2)+1
hajek.kernel.unif <- function(z,s,theta,mean.s=mean.s) {
        2*mu0(z,s,theta,mean.s=mean.s) - mu0.uncond(theta)#2*theta/3*(1-theta/4)
    }
hajek.unif <- function(z,s,theta,mean.s=1/2){
    ## theta.fe <- sum(z*s)/sum(s^2)
    2*mean(mu0(z,s,theta,mean.s=mean.s)) - mu0.uncond(theta)#2*theta/3*(1-theta/4)#- mu0(z,s,0,mean.s=mean.s))
}
## ## renaming mu0 to mu.cond, using mu to refer to the true mean
## mu.cond <- function(z1,s1,theta.2,mean.s=1/2)4*E.diff.cond(z1+1/2-theta.2*(s1-(mean.s-1/2)), theta.2, s1-(mean.s-1/2)) - 2*E.diff.cond(z1+1/2-theta.2*(s1-mean.s+1/2), theta.2,1)-2*(s1-mean.s+1/2)+1
## mu <- function(theta)
## hajek.kernel.unif <- function(z,s,theta,mean.s=mean.s) {
##         2*mu.cond(z,s,theta,mean.s=mean.s) - 2*theta/3*(1-theta/4)
##     }
## hajek.unif <- function(z,s,theta,mean.s=1/2){
##     ## theta.fe <- sum(z*s)/sum(s^2)
##     2*mean(mu.cond(z,s,theta,mean.s=mean.s)) - 2*theta/3*(1-theta/4)#- mu.cond(z,s,0,mean.s=mean.s))
## }
## E.uncond <- function(theta2) {
##     L0 <- pmin(1/theta2,1)
##     L2 <- pmin(1,-1/theta2)
##     (L2/2+1/2*(theta2-1/2)*L2^2+1/3*(theta2^2/2-theta2)*L2^3-theta2^2/8*L2^4) * (theta2<0) +   (1/2-L0/2+(theta2/2+1/4)*L0^2-(theta2^2/6+theta2/3)*L0^3+1/8*theta2^2*L0^4)*(theta2>0) + (1/4)*(theta2==0)
## }



## formulas for uniform model at alternative (pubbias #5)
## this unif.Z does not have var(Z)=1
unif.Z <- local({
    rZ <- function(n)runif(n,-1/2,1/2)
    dZ <- function(z)dunif(z,-1/2,1/2)
    pZ <- function(q)punif(q,-1/2,1/2)
    qZ <- function(u) (u-1/2)*(u>0)*(u<1)
    ## dZ.c <- function(z,cutoff) dunif(z,cutoff,1/2)
    rZ.c <- function(n,cutoff) runif(n,cutoff,1/2)#qZ(runif(n)*(1-pZ(cutoff))+pZ(cutoff))
    dZ.c <- function(z,cutoff) dunif(z,cutoff,1/2)
    ## pZ.c <- function(q,cutoff)(q>cutoff)*(pZ(q)-pZ(cutoff))/(1-pZ(cutoff))
    theta.to.cutoff <- function(theta)(2*theta-1/2)
    ## kernel.mean <- 1/6
    F.delta.c <- function(q,cutoff) (q>=0 & q<=1/2-cutoff)*(-q^2/2/(1/2-cutoff)^2+1/2+q/(1/2-cutoff)) + (q<0 & q>=-1/2+cutoff)*(1/2*(q/(1/2-cutoff)+1)^2) + 1*(q>1/2-cutoff)
    c.prime <- function(theta)2
    theta.max <- 1/3
    supp.Z <- c(-1/2,1/2)
    E.f.Z <- function(cutoff=-1/2)1/(1/2-cutoff) #deprecate
    E.f <- function(cutoff=-1/2)1/(1/2-cutoff)
    E.F.Z <- function(cutoff=-1/2)(1/3*(1/8-cutoff^3)-cutoff/2*(1/4-cutoff^2))/(1/2-cutoff)^2 #deprecate
    E.ZF <- function(cutoff=-1/2)(1/3*(1/8-cutoff^3)-cutoff/2*(1/4-cutoff^2))/(1/2-cutoff)^2
    var.Z <- function(cutoff=-1/2)1/12 / (1/2-cutoff)^2
    return(mget(ls()))
})
## updated version with var(Z)=1--need to finish
unif.Z <- local({
    s <- 1/sqrt(1/12) # Z=s*Unif(-1/2,1/2)
    rZ <- function(n)runif(n,-1/2,1/2)*s
    dZ <- function(z)dunif(z/s,-1/2,1/2)/s
    pZ <- function(q)punif(q/s,-1/2,1/2)
#    qZ <- function(p)qunif(p,-1/2,1/2)# (u-1/2)*(u>0)*(u<1)
    ## ## dZ.c <- function(z,cutoff) dunif(z,cutoff,1/2)
    ## rZ.c <- function(n,cutoff) runif(n,cutoff,1/2)#qZ(runif(n)*(1-pZ(cutoff))+pZ(cutoff))
    ## dZ.c <- function(z,cutoff) dunif(z,cutoff,1/2)
    ## ## pZ.c <- function(q,cutoff)(q>cutoff)*(pZ(q)-pZ(cutoff))/(1-pZ(cutoff))
    ## theta.to.cutoff <- function(theta)(2*theta-1/2)
    ## ## kernel.mean <- 1/6
    ## F.delta.c <- function(q,cutoff) (q>=0 & q<=1/2-cutoff)*(-q^2/2/(1/2-cutoff)^2+1/2+q/(1/2-cutoff)) + (q<0 & q>=-1/2+cutoff)*(1/2*(q/(1/2-cutoff)+1)^2) + 1*(q>1/2-cutoff)
    ## c.prime <- function(theta)2
    ## theta.max <- 1/3
    ## supp.Z <- c(-1/2,1/2)
    ## E.f.Z <- function(cutoff=-1/2)1/(1/2-cutoff) #deprecate
    E.f.c <- function(cutoff=-1/2)1/s/(1/2-cutoff)
    ## E.F.Z <- function(cutoff=-1/2)(1/3*(1/8-cutoff^3)-cutoff/2*(1/4-cutoff^2))/(1/2-cutoff)^2 #deprecate
    E.ZF.c <- function(cutoff=-1/2) s/(1/2-cutoff)*(1/3*(1/4+cutoff/2+cutoff^2)-cutoff/2*(1/2+cutoff))#(1/3*(1/8-cutoff^3)-cutoff/2*(1/4-cutoff^2))/(1/2-cutoff)^2
    E.f <- E.ZF <- 1/s
    var.Z <- function(cutoff=-1/2)1/12 / (1/2-cutoff)^2
    return(mget(ls()))
})
normal.Z <- local({
    rZ <- rnorm
    dZ <- dnorm
    pZ <- pnorm
    qZ <- qnorm
    rZ.c <- function(n,c) rnorm.trunc(n,0,1,c)
    mu <- function(x)exp(log(dnorm(x))-log(1-pnorm(x)))
    theta.to.cutoff <- function(theta)uniroot(function(x)mu(x)-theta,c(-1,1),extendInt='yes')$root  #function(theta)2*theta-1/2
    kernel.mean <- integrate(function(x)dnorm(x)^2,-Inf,Inf)$val
    theta.max <- 1
    supp.Z <- c(-Inf,Inf)
    E.f <- E.ZF <- 1/(2*sqrt(pi))
    return(mget(ls()))
})

power.Z <- function(p) {
    stopifnot(p>=-1)
    const <- ((p+1)/2)^(1/(p+1))
    sigma <- sqrt(2/(p+3)*const^(p+3))
    dZ <- function(z)abs(z)^p*sigma^(p+1)*(abs(z)<=const/sigma)
    pZ <- function(q) 1/2+abs(q)^p*q*sigma^(p+1)/(p+1)*(abs(q)<=const/sigma) + (q>=const/sigma)
    pZ <- function(q) (1/2+abs(q)^p*q*sigma^(p+1)/(p+1))*(abs(q)<=const/sigma) + (q>=const/sigma)
    qZ <- function(y) (abs(y-1/2)*(p+1)/sigma^(p+1))^(1/(p+1))*(y>=0)*(y<=1)*sign(y-1/2)
    rZ <- function(n) qZ(runif(n))
    E.f <- 2/(2*p+1)*sqrt(2/(p+3))*((p+1)/2)^(5/2)
    E.ZF <- const^(p+2)/sigma/(p+1)*(  1/2*(p+1)/(p+2) + const^(p+1)/(p+2)/(2*p+3) )
## rhs <- function(p) {
## c <- ((p+1)/2)^(1/(p+1))
## sigma <- sqrt(2/(p+3)*c^(p+3))
##     2*c^(p+2)/sigma/(p+1)*(  1/2*(p+1)/(p+2) + c^(p+1)/(p+2)/(2*p+3) )
## }
##     E.f <- 1/sigma*(1/2-const)
##     E.ZF <- sigma/(1/2-const)^2*(1/3*(1/4+const/2+const^2)-const/2*(1/2+const))
    return(mget(ls()))
    ## return(list(        
    ##     dZ <- function(z)abs(z)^p*sigma^(p+1)*(abs(z)<=const/sigma),
    ##     pZ <- function(q) 1/2+abs(q)^p*q*sigma^(p+1)/(p+1)*(abs(q)<=const/sigma) + (q>=const/sigma),
    ##     qZ <- function(y) (abs(y-1/2)*(p+1)/sigma^(p+1))^(1/(p+1))*(y>=0)*(y<=1)*sign(y-1/2),
    ##     rZ <- function(n) qZ(runif(n)),
    ## ))
}

beta.Z <- function(ab=c(1,1)) {
    a <- ab[1]; b <- ab[2]
    var.beta <- a*b/(a+b)^2/(a+b+1)
    sd.beta <- sqrt(var.beta)
    rZ <- function(n)(rbeta(n,a,b)-1/2)/sd.beta
    dZ <- function(x)dbeta(1/2+x*sd.beta,a,b)*sd.beta
    pZ <- function(q)pbeta(1/2+q*sd.beta,a,b)
    return(mget(ls()))
}


unif.S <- function(supp.S=c(0,1)) {
    return(list(
    dS=function(s)dunif(s,supp.S[1],supp.S[2]),
    pS=function(q)punif(q,supp.S[1],supp.S[2]),
    rS=function(n)runif(n,supp.S[1],supp.S[2]),
    E.S2=diff(supp.S)^2/12 + mean(supp.S)^2,
    E.S1=mean(supp.S),
    m=function(k)integrate(function(s)s^k*dS(s),supp.S[1],supp.S[2]),
    mean.S.pair=with(list(b=supp.S[2],a=supp.S[1]), 1/(b-a)*((a^2+b^2)/3-2/3*a*b))
    ## return(mget(ls()))
    ))
}
unif.S <- function(supp.S=c(0,1)) {
    ## return(list(
    dS <- function(s)dunif(s,supp.S[1],supp.S[2])
    pS <- function(q)punif(q,supp.S[1],supp.S[2])
    rS <- function(n)runif(n,supp.S[1],supp.S[2])
    E.S2 <- diff(supp.S)^2/12 + mean(supp.S)^2
    E.S1 <- mean(supp.S)
    mu.S <- function(k)integrate(function(s)s^k*dS(s),supp.S[1],supp.S[2])$value
    mean.S.pair <- with(list(b=supp.S[2],a=supp.S[1]), 1/(b-a)*((a^2+b^2)/3-2/3*a*b))
    return(mget(ls()))    
}

tau <- function(theta,cutoff=-1/2) { 
    l0 <- pmin(1,pmax(0,(1/2-cutoff)/abs(theta)))
    triangle <- -theta^2/2/(1/2-cutoff)^2*(l0^3/3-l0^4/4)+abs(theta)/(1/2-cutoff)*(l0^2/2-l0^3/3)+l0^2/4-l0/2+1/2    -1/4
    E.uncond <- 2*(1/4+sign(theta)*triangle)
    2*E.uncond - 1
}
kernel.cond <- function(z,s,theta,cutoff) {
    ## stopifnot(theta>0)
    z1 <- z; s1 <- s
    l1 <- 1/theta*(cutoff-z1+theta*s1)
    l2 <- 1/theta*(1/2-z1+theta*s1)
    1-s1 + 2*
        (if(theta>0) with(list(thresh=function(x)pmin(s1,pmax(0,x))), ((1/2-cutoff)*thresh(l1)+(1/2-z1+theta*s1)*(thresh(l2)-thresh(l1))-theta/2*(thresh(l2)^2-thresh(l1)^2))/(1/2-cutoff)) else   with(list(thresh=function(x)pmin(s1,pmax(0,x))), ((1/2-cutoff)*(s1-thresh(l1))+(1/2-z1+theta*s1)*(thresh(l1)-thresh(l2))-theta/2*(thresh(l1)^2-thresh(l2)^2))/(1/2-cutoff)))    -               with(list(thresh=function(x)pmin(1,pmax(0,x))),(theta<0)+sign(theta)*(((1/2-cutoff)*thresh(l1)+(1/2-z1+theta*s1)*(thresh(l2)-thresh(l1))-theta/2*(thresh(l2)^2-thresh(l1)^2))/(1/2-cutoff)))
    }       
tau.hat <- function(z,s,theta) {
    n <- length(z)
    2*mean(apply(combn(n,2),2,function(idx)(z[idx[1]]-z[idx[2]])/(s[idx[1]]-s[idx[2]])<theta)) - 1
    }
## same as cor.test((y-theta)/sqrt(v),v,method='kendall')$estimate
tau.hat.pi <- function(z,s,theta,cutoff=-1/2) 4*mean(kernel.cond(z,s,theta,cutoff))-2-tau(theta,cutoff)
