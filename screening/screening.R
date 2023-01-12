## 18 effect of screening on meta-analysis -- begg, asymptotics

## 18a screened vs. non-screened theta.hats. normal z (no effect) and
## power law (noticeable effect). With power law z, effect is stronger
## for more conservative p-vals.
source('../begg/2/misc.R')
n <- 4
B <- 1e4
alpha <- .1
## theta.fes <- c()
pairs <- replicate(B, {
    ## z <- with(power.Z(-.7),rZ(n))
    ## z <- rnorm(n)
    ## z <- rt(n,df=2.4) #standardize!
    z <- with(beta.Z(c(.3,.3)),rZ(n))
    s <- runif(n)+1
    y <- z/s
    v <- 1/s^2
    ## v <- (runif(n,1,2)*8)^2
    ## y <- rnorm(n,mean=2,sd=sqrt(v))
    ## pval <- egger.test(y,v)['pval']
    pval <-  begg.test(y,v)['pval']    
    ## pval <- lin.test(y,v)['pval']
    c(pval=unname(pval),theta.fe=sum(y/v)/sum(1/v) )
    ## if (pval < alpha) {
        ## theta.fes <<- c(theta.fes, sum(y/v)/sum(1/v) * sqrt(sum(1/v)) )
    ## }
})
null.idx <- pairs['pval',] > alpha
null.theta <- pairs['theta.fe',null.idx]
alt.theta <- pairs['theta.fe',!null.idx]
qqplot(null.theta,alt.theta,cex=.5); abline(0,1,col=2)

## qqnorm(theta.fes);abline(0,1)
## mean(theta.fes)

source('../begg/2/misc.R')
n <- 15
B <- 1e4
pairs <- with(power.Z(1), {
    replicate(B, {
        z <- rZ(n)
        s <- runif(n)
        y <- z/s
        v <- 1/s^2
        ## theta.fe <- sum(z*s)/sum(s^2)
        ## sqrt(n)*tau.hat(z,s,theta.fe) #/ sqrt(asy.var) #1
        pval <-  begg.test(y,v)['pval']    
        c(pval=unname(pval),theta.fe=sum(y/v)/sum(1/v) )
        ## if (pval < alpha) {
        ## theta.fes <<- c(theta.fes, sum(y/v)/sum(1/v) * sqrt(sum(1/v)) )
        ## }
    })
})    
alpha <- .05
null.idx <- pairs['pval',] > alpha
null.theta <- pairs['theta.fe',null.idx]
alt.theta <- pairs['theta.fe',!null.idx]
qqplot(null.theta,alt.theta,cex=.5); abline(0,1,col=2)




## 18aa p-value of screened theta.hats. For z ~ student's t, the high
## fpr isn't from screening but from student's t slow clt convergence
## (see following sim). #2 bimodal beta -- underpowered. #3 power Z
## with p close to -1 -- underpowered.
source('../begg/2/misc.R')
n <- 2e2
B <- 1e3
## a <- b <- 1; var.beta <- a*b/(a+b)^2/(a+b+1) #2
pairs <- with(power.Z(-.5), {
    replicate(B, {
        z <- rZ(n) #3
        ## z <- rnorm(n)
        ## z <- rt(n,df=3) #1 standardize!
        ## z <- (rbeta(n,a,b)-1/2) / sqrt(var.beta) #2
        s <- runif(n)+1
        y <- z/s
        v <- 1/s^2
        ## theta.fe <- sum(z*s)/sum(s^2)
        ## sqrt(n)*tau.hat(z,s,theta.fe) #/ sqrt(asy.var) #1
        pb0 <-  begg.test(y,v)
        ## pb0 <-  egger.test(y,v)        
        theta.fe <- sum(y/v)/sum(1/v)
        c(pb.pval=unname(pb0['pval']),pb.stat=unname(pb0['stat']),theta.fe=theta.fe,ma.stat=theta.fe/sqrt(1/sum(1/v)))
        ## if (pval < alpha) {
        ## theta.fes <<- c(theta.fes, sum(y/v)/sum(1/v) * sqrt(sum(1/v)) )
        ## }
    })
})    
alpha <- .05
null.idx <- pairs['pb.pval',] > alpha
theta.fe <- pairs['theta.fe',]
test.stat <- pairs['ma.stat',]
null.theta <- theta.fe[null.idx]
alt.theta <- theta.fe[!null.idx]
qqplot(null.theta,alt.theta,cex=.5); abline(0,1,col=2)
qqplot(null.theta,theta.fe,cex=.5); abline(0,1,col=2)
mean(1-pnorm(abs(test.stat[null.idx])) < .05)
## mean(1-pnorm(abs(test.stat[!null.idx])) < .05)
mean(1-pnorm(abs(test.stat)) < .05)
hist(1-pnorm(abs(test.stat[null.idx])))

## 18ab for t with df<4, third moment is infinite, berry esseen doesnt
## apply
n <- 1e4
pvals <- replicate(1e3, {
    z <- rt(n,df=3.5) 
    ## z <- with(power.Z(-.7),rZ(n))
    ## ## z <- rnorm(n)
    ## s <- runif(n)+1
    ## y <- z/s
    ## v <- 1/s^2
    ## ## lm0 <- lm(I(y/sqrt(v)) ~ I(1/sqrt(v)))
    ## ## coef(summary(lm0))[1,4]
    ## theta.fe <- sum(z*s)/sum(s^2)
    ## ma.stat <- theta.fe / sqrt(1/sum(s^2))
    ## 1-pnorm(ma.stat)
    1-pnorm(sqrt(n)*mean(z))
})
hist(pvals)
mean(pvals<.05)

dd



## 18ac plotting monte carlo power for both egger and begg. for power z,
## E(f(Z)) is infinite for p<=-1/2.
source('../begg/2/misc.R')
n <- 1e1
B <- 1e3
sim <- with(power.Z(-.85), {
    replicate(B, {
        z <- rZ(n)
        ## z <- rnorm(n)
        ## z <- rt(n,df=3) #standardize!
        s <- runif(n)+1
        y <- z/s
        v <- 1/s^2
        theta.fe <- sum(y/v)/sum(1/v)
        ## pb0 <-  begg.test(y,v)
        ## egger0 <-  egger.test(y,v)
        ## begg0 <- begg.test(y,v)
        ## begg0 <- with(unif.S(),  {
        ##         mean.S.pair <- with(unif.S(),mean.S.pair)
        ##         tau.hat.pi <- function(z,s,theta) {
        ##             stopifnot(theta==0) # finish
        ##             4*mean(pS(s)+pZ(z)-2*pS(s)*pZ(z)) - 2
        ##         }
        ##         D <- 2*E.f*mean.S.pair
        ##         c(stat=tau.hat.pi(z,s,0) + theta.fe*D)
        ## })
        ## begg0 <- c(stat=tau.hat(z,s,theta.fe))
        ## c(egger=egger.test(y,v),begg=begg.test(y,v,exact=TRUE),ma.stat=theta.fe/sqrt(1/sum(1/v)))
            begg.stat <- tau.hat(z,s,theta.fe)
            begg.pval <- (1-pnorm(abs(begg.stat)*sqrt(9*n/4)))*2
            c(egger=egger.test(y,v),begg.stat=begg.stat,begg.pval=begg.pval,ma.stat=theta.fe/sqrt(1/sum(1/v)))
        ## c(pb.pval=unname(pb0['pval']),pb.stat=unname(pb0['stat']),theta.fe=theta.fe,ma.stat=theta.fe/sqrt(1/sum(1/v)))
        ## if (pval < alpha) {
        ## theta.fes <<- c(theta.fes, sum(y/v)/sum(1/v) * sqrt(sum(1/v)) )
        ## }
    })
})    
pb.alphas <- seq(0.05,.3,len=50)
ma.power <- sapply(pb.alphas, function(alpha) {
    null.idx <- sim['begg.pval',] > alpha
    begg.power <- mean(1-pnorm(abs(sim['ma.stat',null.idx])) < .05)
    null.idx <- sim['egger.pval',] > alpha
    egger.power <- mean(1-pnorm(abs(sim['ma.stat',null.idx])) < .05)
    c(begg=begg.power,egger=egger.power)
})
plot(pb.alphas,ma.power['begg',],type='l',ylim=range(ma.power),col=2)
lines(pb.alphas,ma.power['egger',],type='l',ylim=range(ma.power))
legend('bottomleft',lty=1,col=1:2,legend=c('egger','begg'))
n*cov(sim['begg.stat',],sim['ma.stat',])

## 18ad fpr versus sample size with z~power(p). doesn't seem v large effect
## for p<.8 esp when pb.alpha=ma.alpha. when ma.alpha<pb.alpha,
## noticeable effect. one suggestion is to match the alpha
## levels. even if correlation is perfect between ma and pb test stat,
## might not be able to detect difference if alpha levels are the
## same. Either way, with p<.5, both egger and begg fpr's are almost
## the same, and p<.5 is the only time the asy normality formula is
## valid.
start <- Sys.time()
source('../begg/2/misc.R')
require(parallel)
B <- 1e2
pb.alpha <- .05
ma.alpha <- .05
ns <- round(seq(1e1,3e2,len=10))
## by.n <- mclapply(ns, mc.cores=detectCores()-2, FUN=function(n) {
    by.n <- lapply(ns, FUN=function(n) {
    cat('.')
    stats <-      replicate(B, {
        ## z <- with(power.Z(1), rZ(n))
        ## z <- rnorm(n)
        ## z <- with(list(df=3),rt(n,df=df)/sqrt(df/(df-2)))
        z <- with(beta.Z(c(.2,.2)),rZ(n))
        s <- runif(n)
        y <- z/s
        v <- 1/s^2
        theta.fe <- sum(y/v)/sum(1/v)
        ## c(egger=egger.test(y,v),begg=begg.test(y,v,exact=TRUE),ma.stat=theta.fe/sqrt(1/sum(1/v)))
        ## browser()
        ## begg.stat <- tau.hat.pi(z,s,0) + theta.fe*D
        begg.stat <- tau.hat(z,s,theta.fe)
        begg.pval <- (1-pnorm(abs(begg.stat)*sqrt(9*n/4)))*2
        c(egger=egger.test(y,v),begg.stat=begg.stat,begg.pval=begg.pval,ma.stat=theta.fe/sqrt(1/sum(1/v)))        
    })
    null.idx <- stats['begg.pval',] > pb.alpha
    begg.power <- mean(1-pnorm(abs(stats['ma.stat',null.idx])) < ma.alpha/2)
    null.idx <- stats['egger.pval',] > pb.alpha
    egger.power <- mean(1-pnorm(abs(stats['ma.stat',null.idx])) < ma.alpha/2)
    unconditional.power <- mean(1-pnorm(abs(stats['ma.stat',])) < ma.alpha/2)
    c(begg=begg.power,egger=egger.power,unconditional=unconditional.power) 
    })
by.n <- simplify2array(by.n)
plot(ns,by.n['egger',],type='l',ylim=range(by.n))
lines(ns,by.n['begg',],col=2)
lines(ns,by.n['unconditional',],col=3)
legend('bottomleft',lty=1,col=1:3,legend=c('egger','begg','unconditional'))
Sys.time()-start


pb.alphas <- seq(.01,.3,len=60)
means <- sapply(pb.alphas, function(pb.alpha) mean(abs(stats['ma.stat',stats['begg.pval',]>pb.alpha])))
## means <- sapply(pb.alphas, function(pb.alpha) mean(abs(stats['egger.stat',stats['egger.pval',]>pb.alpha])))
plot(pb.alphas,means)


## 18ae fpr versus sample size using hajek projection of tau. is formula right?
source('../begg/2/misc.R')
require(parallel)
B <- 5e1
pb.alpha <- .1
ma.alpha <- .05
ns <- round(seq(1e1,2e4,len=10))
by.n <- mclapply(ns, mc.cores=detectCores()-2, FUN=function(n) {
    cat('.')
    stats <- with(unif.S(), with(power.Z(0), {
        mean.S.pair <- with(unif.S(),mean.S.pair)
        tau.hat.pi <- function(z,s,theta) {
            stopifnot(theta==0) # finish
            4*mean(pS(s)+pZ(z)-2*pS(s)*pZ(z)) - 2
        }
        D <- 2*E.f*mean.S.pair
        replicate(B, {
            z <- rZ(n)
            ## z <- rnorm(n) 
            ## z <- rt(n,df=3)
            s <- runif(n)+1
            y <- z/s
            v <- 1/s^2
            theta.fe <- sum(y/v)/sum(1/v)
            begg.stat <- tau.hat.pi(z,s,0) + theta.fe*D
            begg.pval <- (1-pnorm(abs(begg.stat)*sqrt(9*n/4)))/2
            c(egger=egger.test(y,v),begg.pval=begg.pval,ma.stat=theta.fe/sqrt(1/sum(1/v)))
        })
    }))
    null.idx <- stats['begg.pval',] > pb.alpha
    begg.power <- mean(1-pnorm(abs(stats['ma.stat',null.idx])) < ma.alpha)
    null.idx <- stats['egger.pval',] > pb.alpha
    egger.power <- mean(1-pnorm(abs(stats['ma.stat',null.idx])) < ma.alpha)
    c(begg=begg.power,egger=egger.power)    
})
by.n <- simplify2array(by.n)
plot(ns,by.n['egger',],type='l',ylim=range(by.n))
lines(ns,by.n['begg',],col=2)
legend('bottomleft',lty=1,col=1:2,legend=c('egger','begg'))


dd

## 18af looking at E.f-E.ZF over a range of values

## power Z. note E.f=infinity for p<=-1/2
ps <- seq(-.9,0,len=1e3)
plot(ps,sapply(ps, function(p)with(power.Z(p),E.f-E.ZF)),type='l')
abline(h=0,,v=-.5,lty=2)

df <- 2.001
sd.z <- sqrt(df/(df-2))
rZ <- function(n)rt(n,df=df)/sd.z
pZ <- function(q)pt(q*sd.z,df=df)
dZ <- function(z)dt(z*sd.z,df=df)*sd.z
integrate(function(z)dZ(z)^2,-50,50)$val - integrate(function(z)z*pZ(z)*dZ(z),-50,50)$val

     
eps <- .001
ps <- seq(-.5-eps,-.5+eps,len=1e1)
asy.covs <- sapply(ps, function(p){
    with(unif.S(), with(power.Z(p), {
        D <- 2*E.f*mean.S.pair
        cov.theta.pi <- -2/E.S2*mean.S.pair*E.ZF
        var.Z <- 1
        asy.cov <- D/E.S2 + cov.theta.pi
        ## print(asy.cov)
        ## print(2*mean.S.pair/E.S2*(E.f-E.ZF))
    }))
})
plot(ps,asy.covs,type='l')
abline(h=0,,v=-.5,lty=2)

p <- runif(1,-1,6)
B <- 1e5
with(power.Z(p), abs(diff(c(E.f,mean(dZ(rZ(B)))))))
with(power.Z(p), {
    z <- rZ(B)
    abs(diff(c(E.ZF,mean(z*pZ(z)))))
})

## beta Z

source('../begg/2/misc.R')
ps <- seq(.3,1,len=1e1)
by.p <- sapply(ps, function(p) {
    print(p)
    rZ <- with(beta.Z(c(p,p)), rZ)
    dZ <- with(beta.Z(c(p,p)), dZ)
    pZ <- with(beta.Z(c(p,p)), pZ)
    ## integrate(function(z)dZ(z)^2,-5,5)$val - integrate(function(z)z*pZ(z)*dZ(z),-5,5)$val
    beta(2*p-1,2*p-1)/beta(p,p)^2 / (2*sqrt(2*p+1))#- integrate(function(z)z*pZ(z)*dZ(z),-5,5)$val
})
plot(ps,by.p)

dd

## plot(2*(1-pnorm(abs(pairs['pb.stat',])*sqrt(9*n/4))), pairs['pb.pval',]);abline(0,1) # for begg



conv <- function(z)integrate(function(x)pnorm(x)*pnorm(z-x),-10,10)$val
conv <- Vectorize(conv)
curve(conv,-2,2)

test.stats <- replicate(B, {
    z <- rnorm(n)
    s <- runif(n)
    y <- z/s
    v <- 1/s^2
    theta.fe <- sum(z*s)/sum(s^2)
    var.theta.fe <- 1/sum(1/v)
    theta.fe / sqrt(var.theta.fe)
})
qqnorm(test.stats); abline(0,1,col=2)

dd



## 18b #1 asy variance of tau.hat. #2 asy covariance of tau.hat and
## theta.fe

require(parallel)
## source('misc.R')
source('../begg/2/misc.R')
B <- 5e3
ns <- round(seq(1e1,1e3,len=40))
## n <- 2e2
supp.S <- c(0,1)+1
mean.S.pair <- with(unif.S(supp.S),mean.S.pair)
E.S2 <- with(unif.S(supp.S), E.S2)
rS <- with(unif.S(supp.S),rS)
by.n <- with(power.Z(-.49), {
    ## E.f <- with(unif.Z,E.f())#with(normal.Z,E.f)
    ## E.ZF <- with(unif.Z,E.ZF())#with(normal.Z,E.ZF)
    ## rZ <- with(unif.Z,rZ)#with(normal.Z,E.ZF)
    D <- 2*E.f*mean.S.pair
    cov.theta.pi <- -2/E.S2*mean.S.pair*E.ZF
    var.Z <- 1
    asy.var <- 4/9 + var.Z*D^2/E.S2 + 2*D*cov.theta.pi
    ## asy.var <- with(unif.S(), with(unif.Z, 4/9 + 4*mean.S.pair^2/E.S2*E.f.Z()*(E.f.Z()*var.Z()-2*E.F.Z())  ))
    ## asy.var <- 4/9 + var.Z*4*E.f*mean.S.pair^2/E.S2*(E.f-2*E.ZF)
    asy.cov <- 2*mean.S.pair/sqrt(E.S2)*(E.f-E.ZF)
    print(asy.cov)
    print(asy.var)
    mclapply(ns, mc.cores=detectCores()-3, FUN=function(n) {
        cat('.')
        ## by.n <- sapply(ns, function(n)
        stats <- replicate(B, {
            z <- rZ(n)
            s <- rS(n)
            ## y <- z/s
            ## v <- 1/s^2
            theta.fe <- sum(z*s)/sum(s^2)
            ## sqrt(n)*tau.hat(z,s,theta.fe) #/ sqrt(asy.var) #1
            ## c(theta.fe*sqrt(sum(s^2)),sqrt(n)*tau.hat(z,s,theta.fe)) #2
            tau.hat.pi <- function(z,s,theta) {
                stopifnot(theta==0) # finish
                with(unif.S(supp.S),4*mean(pS(s)+pZ(z)-2*pS(s)*pZ(z)) - 2)
            }
            sqrt(n)*(tau.hat.pi(z,s,0) + theta.fe*D)
        })
        ## print(         cov(t(stats))[1,2] )
        print(var(stats))
        var(stats) - asy.var #1
        ## cov(t(stats))[1,2] - asy.cov #2
    })
})
by.n <- simplify2array(by.n)
plot(ns,by.n)
abline(h=0,lty=2)
sd.by.n <- print(sd(by.n))
abline(h=mean(by.n),col=2)



bdd

## op <- par(mfrow=c(1,2))
## hist(stats,freq=FALSE)
## curve(dnorm,add=TRUE)
## qqnorm(stats);abline(0,1,col=2)
## par(op)


## 18c L2 equivalence of taylor/hajek simplification and tau.hat, for
## power-law distributed Z. Seems equivalent except for p close to
## -1/2 (perhaps as should have been expected).
## source('misc.R')
start <- Sys.time()
source('../begg/2/misc.R')
require(parallel)
B <- 1e2
supp.S <- c(1,2)+2
## n <- 10
ns <- round(seq(1e2,1e3,len=10))
## by.n <- lapply(ns,  FUN=function(n) {
by.n <- mclapply(ns, mc.cores=detectCores()-3, FUN=function(n) {
            cat('.')
            with(unif.S(supp.S), {
                ## mean.S.pair <- mean.S.pair
                with(power.Z(-.45), {
                    tau.hat.pi <- function(z,s,theta) {
                        stopifnot(theta==0) # finish
                        4*mean(pS(s)+pZ(z)-2*pS(s)*pZ(z)) - 2
                    }
                    ## asy.cov <- D/E.S2 + cov.theta.pi
                    ## E.f <- with(unif.Z,E.f())#with(normal.Z,E.f)
                    ## E.ZF <- with(unif.Z,E.ZF())#with(normal.Z,E.ZF)
                    ## rZ <- with(unif.Z,rZ)#with(normal.Z,E.ZF)
                    D <- 2*E.f*mean.S.pair
                    ## by.n <- lapply(ns,  FUN=function(n) {
                    ## cutoff <- unif.Z$theta.to.cutoff(1/sqrt(n))
                    ## print(cutoff)
                    ## by.n <- sapply(ns, function(n)
                    stats <- replicate(B, {
                        z <- rZ(n)
                        s <- rS(n)
                        y <- z/s
                        v <- 1/s^2
                        theta.fe <- sum(z*s)/sum(s^2)
                        ## c(tau.hat.pi(z,s,theta.fe,-1/2),tau.hat(z,s,theta.fe))
                        ## c(tau.hat.pi(z,s,0) + theta.fe*D,tau.hat(z,s,theta.fe))
                        ## c(tau.hat.pi(z,s,0) + theta.fe*D,begg.test(y,v,exact=NULL)['stat'])
                        c(tau.hat.pi(z,s,0) + theta.fe*D,tau.hat(z,s,theta.fe))
                    })
                    ## n*mean(stats^2)
                    n*mean((stats[1,]-stats[2,])^2)
                })
            })
})
by.n <- simplify2array(by.n)
plot(ns,by.n); abline(h=0,lty=2)
Sys.time() - start

## 18d formula for terms in asy covariance, and also asy covariance itself, for power law distribution.
## source('misc.R')
source('../begg/2/misc.R')
require(parallel)
B <- 3e1
## n <- 10
ns <- round(seq(1e2,3e3,len=50))
ns <- round(seq(1e1,1e3,len=50))
## by.n <- lapply(ns,  FUN=function(n) {
with(unif.S(), {
    with(power.Z(-.49), {
        tau.hat.pi <- function(z,s,theta) {
            stopifnot(theta==0) # finish
            4*mean(pS(s)+pZ(z)-2*pS(s)*pZ(z)) - 2
        }
        D <- 2*E.f*mean.S.pair
        cov.theta.pi <- -2/E.S2*mean.S.pair*E.ZF
        var.Z <- 1
        asy.cov <- D/E.S2 + cov.theta.pi
        ## by.n <- lapply(ns,  FUN=function(n) {
        by.n <- mclapply(ns, mc.cores=detectCores()-3, FUN=function(n) {
            cat('.')
            stats <- replicate(B, {
                z <- rZ(n)
                s <- rS(n)
                y <- z/s
                v <- 1/s^2
                theta.fe <- sum(z*s)/sum(s^2)
                ## c(tau.hat.pi(z,s,0) + theta.fe*D,theta.fe)
                c(begg.test(y,v)['stat'],theta.fe)
                ## c(tau.hat(z,s,theta.fe),theta.fe)
                ## c(tau.hat(z,s,theta.fe),tau.hat.pi(z,s,0) + theta.fe*D)
                ## c(D*theta.fe,tau.hat.pi(z,s,0))
            })
            ## n*mean((stats[1,]-stats[2,])^2)
            ## n*var(stats)
            n*cov(t(stats))[1,2]
        })
        by.n <- simplify2array(by.n)
        ## plot(ns,by.n); abline(h=1/E.S2,lty=2)
        ## plot(ns,by.n); abline(h=-2*mean.S.pair*E.ZF/E.S2,lty=2)
        plot(ns,by.n); abline(h=asy.cov)
    })
})



## 18e effect of choice of publication bias test alpha on MA stat
## for captions with tikzDevice got idea from https://stackoverflow.com/questions/6445439/tikzdevice-add-caption-and-label-to-tikz-diagram-using-r
## start <- Sys.time()
source('../begg/2/misc.R')
B <- 3e3
n <- 50
rZs <- list(power=with(power.Z(-.85), rZ), beta=with(beta.Z(c(.15,.15)), rZ))
means <- lapply(rZs, function(rZ) {
    stats <-      replicate(B, {
        ## z <- with(power.Z(-.8), rZ(n))
        ## z <- rnorm(n)
        ## z <- with(list(df=3),rt(n,df=df)/sqrt(df/(df-2)))
        ## z <- with(beta.Z(c(.2,.2)),rZ(n))
        z <- rZ(n)
        s <- runif(n)
        y <- z/s
        v <- 1/s^2
        theta.fe <- sum(y/v)/sum(1/v)
        begg.stat <- tau.hat(z,s,theta.fe)
        begg.pval <- (1-pnorm(abs(begg.stat)*sqrt(9*n/4)))*2
        c(egger=egger.test(y,v),begg.stat=begg.stat,begg.pval=begg.pval,ma.stat=theta.fe/sqrt(1/sum(1/v)))        
    })
    pb.alphas <- seq(.01,.3,len=60)
    means.by.begg <- sapply(pb.alphas, function(pb.alpha) mean(abs(stats['ma.stat',stats['begg.pval',]>pb.alpha])))
    means.by.egger <- sapply(pb.alphas, function(pb.alpha) mean(abs(stats['ma.stat',stats['egger.pval',]>pb.alpha])))
    cbind(begg=means.by.begg, egger=means.by.egger)
})

for(i in 1:2) {
    pdf(paste0('fig_alpha_0_',names(rZs)[i],'.pdf'))
    plot(pb.alphas,means[[i]][,'begg'],type='l',ylim=range(means[[i]]),xlab=expression(paste(alpha[0])),ylab='mean meta-analysis test statistic',lty=1)
    lines(pb.alphas,means[[i]][,'egger'],lty=2)
    legend('bottomleft',legend=c('Begg','Egger'),lty=1:2)
    dev.off()
}

dd






## version using tikzDevice--decided against

## for captions with tikzDevice got idea from https://stackoverflow.com/questions/6445439/tikzdevice-add-caption-and-label-to-tikz-diagram-using-r
## start <- Sys.time()
source('../begg/2/misc.R')
library(tikzDevice)
## options(tikzLatexPackages
##         =c(getOption( "tikzLatexPackages" ),"\\usepackage{amssymb}","\\usepackage{amsthm}","\\usepackage{bbm}"))
## require(parallel)
B <- 1e2
n <- 50
rZs <- list(power=with(power.Z(-.8), rZ), beta=with(beta.Z(c(.2,.2)), rZ))
for (i in 1:2) {
    stats <-      replicate(B, {
        z <- with(power.Z(-.8), rZ(n))
        ## z <- rnorm(n)
        ## z <- with(list(df=3),rt(n,df=df)/sqrt(df/(df-2)))
        ## z <- with(beta.Z(c(.2,.2)),rZ(n))
        z <- rZs[[i]](n)
        s <- runif(n)
        y <- z/s
        v <- 1/s^2
        theta.fe <- sum(y/v)/sum(1/v)
        begg.stat <- tau.hat(z,s,theta.fe)
        begg.pval <- (1-pnorm(abs(begg.stat)*sqrt(9*n/4)))*2
        c(egger=egger.test(y,v),begg.stat=begg.stat,begg.pval=begg.pval,ma.stat=theta.fe/sqrt(1/sum(1/v)))        
    })
    pb.alphas <- seq(.01,.3,len=60)
    means.by.begg <- sapply(pb.alphas, function(pb.alpha) mean(abs(stats['ma.stat',stats['begg.pval',]>pb.alpha])))
    means.by.egger <- sapply(pb.alphas, function(pb.alpha) mean(abs(stats['ma.stat',stats['egger.pval',]>pb.alpha])))
    ## means <- sapply(pb.alphas, function(pb.alpha) mean(abs(stats['egger.stat',stats['egger.pval',]>pb.alpha])))
    ## plot(pb.alphas,means.by.begg,type='l',ylim=range(c(means.by.begg,means.by.egger)))
    ## lines(pb.alphas,means.by.egger,col=2)
    ## sink(paste0('fig_alpha_0',names(rZs)[i],'.tex'))
    ## cat("\\begin{figure}\n")
    ## tikz('simpleEx.tex',width=3.5,height=3.5)
    ## tikz(console=TRUE,width=3.5,height=3.5)
    tikz(paste0('fig_alpha_0_',names(rZs)[i],'.tex'),width=3.5,height=3.5)
    plot(pb.alphas,means.by.begg,type='l',ylim=range(c(means.by.begg,means.by.egger)),xlab='$\\alpha_0$')
    lines(pb.alphas,means.by.egger,col=2)
    dev.off()
    ## invisible(dev.off())
    ## cat(paste("\\caption{figure}\\label{fig:alpha_0}\n",sep=""))
    ## cat("\\end{figure}\n")
    ## sink()
}

dd


## 19 egger
source('misc.R')
require(parallel)
B <- 4e4
## n <- 10
ns <- round(seq(1e2,5e2,len=5))
## with(power.Z(1), {
with(unif.S(c(2,3)), {
    with(normal.Z, {
        E.f <- E.f#with(normal.Z,E.f)#with(unif.Z,E.f())
        E.ZF <- E.ZF#with(normal.Z,E.ZF)#with(unif.Z,E.ZF())
        D <- 2*E.f*mean.S.pair
        by.n <- lapply(ns, FUN=function(n) {
            ## by.n <- mclapply(ns, mc.cores=detectCores()-3, FUN=function(n) {
            cat('.')
            ## by.n <- sapply(ns, function(n)
            stats <- replicate(B, {
                browser()
                z <- rZ(n)#runif(n)-1/2
                s <- rS(n)
                y <- z/s
                v <- 1/s^2
                theta.fe <- sum(z*s)/sum(s^2)
                ## hajek <- 2*mean(mu0(z,s,theta.fe) - mu0(z,s,0))
                ## begg <- cor.test((y-theta.fe)/sqrt(v-1/sum(1/v)),v,method='kendall')$estimate - cor.test(y/sqrt(v-1/sum(1/v)),v,method='kendall')$estimate
                ## c(hajek=hajek,begg=unname(begg))
                egger.lm <- lm(I(y/sqrt(v)) ~ I(1/sqrt(v)))
                RSS <- sum(resid(egger.lm)^2)
                tstat <- coef(summary(egger.lm))['(Intercept)','t value']
                m <- sapply(1:4,function(k)mean(s^k))
                numer <- sqrt(n)*(mean(z)*m[2] - mean(z*s)*m[1])
                denom <- sqrt(m[2]*(m[2]-m[1]^2))
                ## tstat <- numer / denom
                sqrt(n)*c(tstat,theta.fe)
            })
            cov(t(stats))[1,2]
        })
        by.n <- simplify2array(by.n)
        plot(ns,by.n)
        asy.cov <- (mu.S(3)*mu.S(2)-mu.S(1)*mu.S(4))/mu.S(2)^(3/2)/sqrt(mu.S(2)-mu.S(1)^2)
        print(asy.cov)
        abline(h=asy.cov)
        abline(h=mean(by.n),col=2)
    })
})

rZ <- rnorm
rS <- runif
n <- 5e1
B <- 2e3
with(unif.S(c(3,4)), {
    with(normal.Z, {
        ## asy.cov <- mu.S(3)-mu.S(1)*mu.S(4)/mu.S(2)
        asy.cov <- (mu.S(3)*mu.S(2)-mu.S(1)*mu.S(4))/mu.S(2)^(3/2)/sqrt(mu.S(2)-mu.S(1)^2)
        print(asy.cov)
        stats <- replicate(1e2, {
            z <- matrix(rZ(n*B),nrow=B)
            s <- matrix(rS(n*B),nrow=B)
            theta.fe <- rowSums(z*s)/rowSums(s^2)
            m <- sapply(1:4,function(k)rowMeans(s^k))
            numer <- rowMeans(z*s^2)*m[,2] - rowMeans(z*s^3)*m[,1]
            denom <- sqrt(m[,2]*(m[,2]-m[,1]^2))
            tstat <- numer / denom
            ## n*cov(theta.fe,numer) - asy.cov
            n*cov(theta.fe,tstat) - asy.cov
        })
        hist(stats)
        abline(v=0); abline(v=mean(stats),col=2)
    })
})



## 20 checking formulas for finite-sample covariance of begg test

## 20a bimodal but seems formula is right, based on hist/abline
n <- 1e1
B1 <- 3e3
## B2 <- 3e2
## pairs <- replicate(B2, {
    s <- runif(n)+3
    m <- sapply(1:4,function(k)mean(s^k))
    obs <- replicate(B1, {
        z <- rnorm(n)
        ## z[1]* (z[1]*(1/(s[1]-s[2])-s[1]/sum(s^2)) < z[2]/(s[1]-s[2])+sum(z[-1]*s[-1])/sum(s^2))
        ## z[1]*((z[1]-z[2])/(s[1]-s[2])<sum(z*s)/sum(s^2))
                z[1]*((z[2]-z[1])/(s[2]-s[1])<sum(z*s)/sum(s^2))
    })
    try <- integrate(function(z)-z*pnorm(z*(1/(s[1]-s[2])-s[1]/n/m[2])/sqrt(  (1/(s[1]-s[2]))^2+ (s[1]+s[2])/n/m[2]/(s[1]-s[2]) - (s[1]/n/m[2])^2 ))*dnorm(z),-5,5)$val
    hist(obs)
    abline(v=c(try,mean(obs)),col=2:3)
    c(mean(obs),try)
## })
## plot(pairs[1,],pairs[2,]);abline(0,1)


n <- 1e3
B1 <- 3e2
B2 <- 3e2
## pairs <- replicate(B2, {
    s <- runif(n)*10
    m <- sapply(1:4,function(k)mean(s^k))
    obs <- replicate(B1, {
        z <- rnorm(n)
        ## z[1]* (z[1]*(1/(s[1]-s[2])-s[1]/sum(s^2)) < z[2]/(s[1]-s[2])+sum(z[-1]*s[-1])/sum(s^2))
        z[1]*( (z[2]-z[3])/(s[2]-s[3]) < sum(z*s)/sum(s^2)  )
    })
    try <- integrate(function(z)-z*pnorm(z*s[1]/n/m[2]/sqrt(2/(s[2]-s[3])^2-1/n/m[2]-(s[1]/n/m[2])^2))*dnorm(z),-5,5)$val
    hist(obs)
    abline(v=c(try,mean(obs)),col=2:3)
    c(mean(obs),try)
## })
## plot(pairs[1,],pairs[2,]);abline(0,1)


c <- runif(1)+1
integrate(function(z)z^2*dnorm(z*c)*dnorm(z),-10,10)$val
(c^2+1)^(-3/2)/sqrt(2*pi)


c <- runif(1)+1
integrate(function(z)z*pnorm(z*c)*dnorm(z),-10,10)$val
1/sqrt(2*pi)*c/sqrt(c^2+1)

## 20b formulas for finite-sample cov(begg stat, theta.fe) conditional
## on sigma

## 20ba
n <- 1e1
B <- 1e4
theta <- runif(1,-1,1)
sigma <- runif(n,1,2)
sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
m <- sapply(1:2,function(k)mean((1/sigma)^k))
## v <- function(sigma)( (-sigma[2]/sigma.bar[2]+1/sigma[2]/n/m[2]*(1/sigma.bar[2]-1/sigma.bar[1]))^2 + (n-2)*(1/n/m[2]*(1/sigma.bar[1]-1/sigma.bar[2]))^2 )*(sigma[1]-sigma[2])^2 # wrong formula
v <- function(sigma) (sigma[2]^2*(-1/sigma.bar[2]+1/sigma[2]^2/n/m[2]*(1/sigma.bar[2]-1/sigma.bar[1]))^2 + (1/n/m[2]*(1/sigma.bar[1]-1/sigma.bar[2]))^2*(m[2]-1/sigma[1]^2-1/sigma[2]^2))*(sigma[1]-sigma[2])^2
v <- function(sigma)(sigma[2]^2/sigma.bar[2]^2-2/sigma.bar[2]/n/m[2]*(1/sigma.bar[2]-1/sigma.bar[1])+1/n^2/m[2]*(1/sigma.bar[1]-1/sigma.bar[2])^2*(1-1/m[2]/sigma[1]^2))*(sigma[1]-sigma[2])^2
obs <- replicate(B, {
    y <- rnorm(n,mean=theta,sd=sigma)
    theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
    ((y[2]-theta)*(-1/sigma.bar[2]+1/sigma[2]^2/n/m[2]/sigma.bar[2]-1/sigma[2]^2/n/m[2]/sigma.bar[1]) - sum((y[3:n]-theta)/sigma[3:n]^2)/n/m[2]*(1/sigma.bar[1]-1/sigma.bar[2]) ) * (sigma[1]-sigma[2])            
})
op <- par(mfrow=c(1,2))
hist(obs/sqrt(v(sigma)),prob=TRUE)
curve(dnorm,add=TRUE,col=2)
qqnorm(obs/sqrt(v(sigma))); abline(0,1,col=2)
par(op)



## ## 20bb 
## n <- 1e1
## B <- 1e4
## theta <- runif(1,-1,1)
## sigma <- runif(n,1,2)
## sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
## m <- sapply(1:2,function(k)mean((1/sigma)^k))
## v <- function(sigma)( (-sigma[2]/sigma.bar[2]+1/sigma[2]/n/m[2]*(1/sigma.bar[2]-1/sigma.bar[1]))^2 + (n-2)*(1/n/m[2]*(1/sigma.bar[1]-1/sigma.bar[2]))^2 )*(sigma[1]-sigma[2])^2
## obs <- replicate(B, {
##     y <- rnorm(n,mean=theta,sd=sigma)
##     theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
##     y[1]*( ((y[1]-theta.fe)/sigma.bar[1] - (y[2]-theta.fe)/sigma.bar[2]) * (sigma[1]-sigma[2]) > 0)
##     ## ((y[1]-theta.fe)/sigma.bar[1] - (y[2]-theta.fe)/sigma.bar[2]) * (sigma[1]-sigma[2])
##     ## ((y[1]-theta)*(1/sigma.bar[1]-1/sigma[1]^2/n/m[2]/sigma.bar[1]+1/sigma[1]^2/n/m[2]/sigma.bar[2]) + (y[2]-theta)*(-1/sigma.bar[2]+1/sigma[2]^2/n/m[2]/sigma.bar[2]-1/sigma[2]^2/n/m[2]/sigma.bar[1]) - sum((y[3:n]-theta)/sigma[3:n]^2)/n/m[2]*(1/sigma.bar[1]-1/sigma.bar[2]) ) * (sigma[1]-sigma[2])
##     ## ((y[2]-theta)*(-1/sigma.bar[2]+1/sigma[2]^2/n/m[2]/sigma.bar[2]-1/sigma[2]^2/n/m[2]/sigma.bar[1]) - sum((y[3:n]-theta)/sigma[3:n]^2)/n/m[2]*(1/sigma.bar[1]-1/sigma.bar[2]) ) * (sigma[1]-sigma[2])            
## })
## ## try <- integrate(function(y) (y)*(1-pnorm(-(y-theta)*(1/sigma.bar[1]-1/sigma[1]^2/n/m[2]/sigma.bar[1]+1/sigma[1]^2/n/m[2]/sigma.bar[2])*(sigma[1]-sigma[2])/sqrt(v(sigma))))*dnorm(y,mean=theta,sd=sigma[1]),theta-5,theta+5)$va 
## ## try <- integrate(function(y) (y-theta)*(1-pnorm(-(y-theta)*(1/sigma.bar[1]-1/sigma[1]^2/n/m[2]/sigma.bar[1]+1/sigma[1]^2/n/m[2]/sigma.bar[2])*(sigma[1]-sigma[2])/sqrt(v(sigma))))*dnorm(y,mean=theta,sd=sigma[1]),theta-5,theta+5)$val + theta/2
## ## try <- -integrate(function(y) (y-theta)*(pnorm(-(y-theta)*(1/sigma.bar[1]-1/sigma[1]^2/n/m[2]/sigma.bar[1]+1/sigma[1]^2/n/m[2]/sigma.bar[2])*(sigma[1]-sigma[2])/sqrt(v(sigma))))*dnorm(y,mean=theta,sd=sigma[1]),theta-5,theta+5)$val + theta/2
## try <- integrate(function(y) (y-theta)*(pnorm((y-theta)*(1/sigma.bar[1]-1/sigma[1]^2/n/m[2]/sigma.bar[1]+1/sigma[1]^2/n/m[2]/sigma.bar[2])*(sigma[1]-sigma[2])/sqrt(v(sigma))))*dnorm(y,mean=theta,sd=sigma[1]),theta-5,theta+5)$val + theta/2
## hist(obs)
## abline(v=c(try,mean(obs)),col=2:3)
## c(mean(obs),try)
##  20bb (fixed)
require(parallel)
n <- 4
B1 <- 1e4
B2 <- 1e2
pairs <- mclapply(1:B2, mc.cores=detectCores()-3,FUN=function(yy) {
## pairs <- replicate(B2,{
    theta <- 8#runif(1,-1,1)
    sigma <- runif(n,1,2)*5
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    jk <- sample(n,2,replace=FALSE); j <- jk[1]; k <- jk[2]
    obs <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        y[j]/sigma[j]^2*( ((y[j]-theta.fe)/sigma.bar[j] - (y[k]-theta.fe)/sigma.bar[k]) * (sigma[j]-sigma[k]) > 0)
    })
    ## try <- integrate(function(y) (y-theta) * pnorm((y-theta)*alpha(sigma,j,k))*dnorm(y,mean=theta,sd=sigma[j]),theta-5,theta+5)$val +theta/2
    numer <- (1/sigma.bar[j]+1/sigma[j]^2/n/m[2]*(1/sigma.bar[k] -1/sigma.bar[j]))
    denom <- sigma[j]^2/sigma.bar[j]^2+sigma[k]^2/sigma.bar[k]^2+(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[j]))^2*(1/sigma[j]^2+1/sigma[k]^2+n-2)-2/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[j])^2 #wrong formula
    denom <- (sigma[j]/sigma.bar[j])^2+(sigma[k]/sigma.bar[k])^2-1/n/m[2]*(1/sigma.bar[j]-1/sigma.bar[k])^2
    ## try <-     sigma[j]^2/sqrt(2*pi)*sign(sigma[j]-sigma[k])*numer/sqrt(denom) +theta/2
    try <-     1/sqrt(2*pi)*sign(sigma[j]-sigma[k])*numer/sqrt(denom) +theta/2/sigma[j]^2
    ## try <- integrate(function(y) (y-theta) * pnorm((y-theta)*beta(sigma,j,k,l))*dnorm(y,mean=theta,sd=sigma[j]),theta-5,theta+5)$val +theta/2
    ## try <-     sigma[j]^2/sqrt(2*pi)*((1/sigma.bar[l]-1/sigma.bar[k])/n/m[2]/sigma[j]^2*sign(sigma[k]-sigma[l])) / sqrt(   (sigma[k]/sigma.bar[k])^2+(sigma[l]/sigma.bar[l])^2 - 1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])^2  ) + theta/2
## ##     ## c(mean(obs),try(y,sigma,theta,j,k,l))
    ## c(mean(obs),mean(try))
    c(mean(obs),try)
})
pairs <- simplify2array(pairs)
plot(pairs[1,],pairs[2,]); abline(0,1)



## double sum version
require(parallel)
n <- 4
B1 <- 1e2
B2 <- 1e2
pairs <- mclapply(1:B2, mc.cores=detectCores()-3,FUN=function(yy) {
## pairs <- replicate(B2,{
    theta <- 2#runif(1,-1,1)
    sigma <- runif(n,1,2)
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    ## jk <- sample(n,2,replace=FALSE); j <- jk[1]; k <- jk[2]
    ## obs <- replicate(B1, {
    ## y <- rnorm(n,mean=theta,sd=sigma)
    ## theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
    ##     y[j]*( ((y[j]-theta.fe)/sigma.bar[j] - (y[k]-theta.fe)/sigma.bar[k]) * (sigma[j]-sigma[k]) > 0)
    ## })
    obs.double.sum <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        sum(sapply(1:n, function(j)
            sapply((1:n)[-j], function(k) 
                y[j]/sigma[j]^2*( ((y[j]-theta.fe)/sigma.bar[j] - (y[k]-theta.fe)/sigma.bar[k]) * (sigma[j]-sigma[k]) > 0)            
                )))
    })    
    try.double.sum <- sum(sapply(1:n, function(j)
        sapply((1:n)[-j], function(k) {
            numer <- (1/sigma.bar[j]+1/sigma[j]^2/n/m[2]*(1/sigma.bar[k] -1/sigma.bar[j]))
            ## denom <- sigma[j]^2/sigma.bar[j]^2+sigma[k]^2/sigma.bar[k]^2+(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[j]))^2*(1/sigma[j]^2+1/sigma[k]^2+n-2)-2/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[j])^2 #wrong formula
            denom <- (sigma[j]/sigma.bar[j])^2+(sigma[k]/sigma.bar[k])^2-1/n/m[2]*(1/sigma.bar[j]-1/sigma.bar[k])^2
            1/sqrt(2*pi)*sign(sigma[j]-sigma[k])*numer/sqrt(denom) +theta/2/sigma[j]^2
        })))
    c(mean(obs.double.sum),try.double.sum)
})
pairs <- simplify2array(pairs)
plot(pairs[1,],pairs[2,]); abline(0,1)


dd


double.sum <- sum(sapply(1:n, function(j)
    sapply((1:n)[-j], function(k)
        1/sigma[j]^2*y[j] * (((y[j]-theta.fe)/sigma.bar[j]-(y[k]-theta.fe)/sigma.bar[k])*(sigma[j]-sigma[k])>0)
        )))




## 20bc
n <- 1e1
B <- 1e4
theta <- runif(1,-1,1)
sigma <- runif(n,1,2)
sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
m <- sapply(1:2,function(k)mean((1/sigma)^k))
## j <- 1;k <- 2
alpha <- function(sigma,j,k) {
    numer <- (1/sigma.bar[j]-1/sigma[j]^2/n/m[2]*(1/sigma.bar[j] -1/sigma.bar[k]))*(sigma[j]-sigma[k])
    denom <- ( (-sigma[k]/sigma.bar[k]+1/n/m[2]/sigma[k]*(1/sigma.bar[k]-1/sigma.bar[j]))^2 + (n-2)*(1/n/m[2]*(1/sigma.bar[j]-1/sigma.bar[k]))^2 ) *(sigma[j]-sigma[k])^2
    numer/sqrt(denom)
}
try <- function(y,sigma,theta,j,k)integrate(function(y) (y-theta)* pnorm((y-theta)*alpha(sigma,j,k))*dnorm(y,mean=theta,sd=sigma[j]),theta-5,theta+5)$val + theta/2
jk <- sample(n,2,replace=FALSE); j <- jk[1]; k <- jk[2]
obs <- replicate(B, {
    y <- rnorm(n,mean=theta,sd=sigma)
    theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
    y[j]*( ((y[j]-theta.fe)/sigma.bar[j] - (y[k]-theta.fe)/sigma.bar[k]) * (sigma[j]-sigma[k]) > 0)
})
hist(obs)
abline(v=c(try(y,sigma,theta,j,k),mean(obs)),col=2:3)
c(mean(obs),try(y,sigma,theta,j,k))



## 20bd
n <- 1e1
B <- 1e4
theta <- runif(1,-1,1)
sigma <- runif(n,1,2)
sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
m <- sapply(1:2,function(k)mean((1/sigma)^k))
beta <- function(sigma,k,l) {
    denom <- (sigma[k]-sigma[l])^2*(1/sigma.bar[k]^2+1/sigma.bar[l]^2+(n-1)*(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[l]))^2)
    numer <- 1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])*(sigma[k]-sigma[l])
    numer/sqrt(denom)
}
## try <- function(y,sigma,theta,j,k,l)integrate(function(y) (y)* pnorm((y/sigma[j]^2)*beta(sigma,k,l))*dnorm(y,mean=theta,sd=sigma[j]),theta-5,theta+5)$val
try <- function(y,sigma,theta,j,k,l)sigma[j]*integrate(function(y) (y-theta)/sigma[j]* pnorm((y/sigma[j]^2)*beta(sigma,k,l))*dnorm(y,mean=theta,sd=sigma[j]),theta-5,theta+5)$val +theta/2
jkl <- sample(n,3,replace=FALSE); j <- jkl[1]; k <- jkl[2]; l <- jkl[3]
obs <- replicate(B, {
    y <- rnorm(n,mean=theta,sd=sigma)
    theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
    y[j]*( ((y[k]-theta.fe)/sigma.bar[k] - (y[l]-theta.fe)/sigma.bar[l]) * (sigma[k]-sigma[l]) > 0)
})
hist(obs)
abline(v=c(try(y,sigma,theta,j,k,l),mean(obs)),col=2:3)
c(mean(obs),try(y,sigma,theta,j,k,l))
try(y,sigma,theta,j,k,l)
b <- beta(sigma,k,l)
theta/2+sigma[j]/sqrt(2*pi)*b/sqrt(b^2+sigma[j]^2)*exp(-1/2*theta^2*b^2/(b^2+sigma[j]^2)/sigma[j]^2)


## ## 20bda
## n <- 4
## B1 <- 4e4
## B2 <- 1e2
## theta <- 1#runif(1,-1,1)
## sigma <- runif(n,1,2)
## sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
## m <- sapply(1:2,function(k)mean((1/sigma)^k))
## jkl <- sample(n,3,replace=FALSE); j <- jkl[1]; k <- jkl[2]; l <- jkl[3]
## obs <- replicate(B1, {
##     y <- rnorm(n,mean=theta,sd=sigma)
##     theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
##     ((y[k]-theta.fe)/sigma.bar[k] - (y[l]-theta.fe)/sigma.bar[l]) * (sigma[k]-sigma[l])
## })
## ## mean1 <- -theta*(sigma[k]-sigma[l])/sigma[j]^2*(1/sigma.bar[l]-1/sigma.bar[k])
## mean1 <- theta*(sigma[k]-sigma[l])*(1/sigma.bar[k]-1/sigma.bar[l])*(1+1/n/m[2]/sigma[j]^2)
## var1 <-  (sigma[k]-sigma[l])^2*(sigma[k]^2/sigma.bar[k]^2+sigma[l]^2/sigma.bar[l]^2+(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[l]))^2*(n*m[2]-1/sigma[j]^2))
## mean2 <- theta/n/m[2]/sigma[j]^2*(sigma[k]-sigma[l])*(1/sigma.bar[l]-1/sigma.bar[k])
## var2 <-(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[l]))^2*(sigma[k]-sigma[l])^2
## try <- rnorm(B1,mean=mean1,sd=sqrt(var1)) + rnorm(B1, mean=mean2,sd=sqrt(var2))    
## ## c(mean(obs),try(y,sigma,theta,j,k,l))
## qqplot(obs,try);abline(0,1,col=2)


## n <- 4
## B1 <- 4e4
## B2 <- 1e2
## theta <- 1#runif(1,-1,1)
## sigma <- runif(n,1,2)
## sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
## m <- sapply(1:2,function(k)mean((1/sigma)^k))
## jkl <- sample(n,3,replace=FALSE); j <- jkl[1]; k <- jkl[2]; l <- jkl[3]
## obs <- replicate(B1, {
##     y <- rnorm(n,mean=theta,sd=sigma)
##     theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
##     ((y[k]-theta.fe)/sigma.bar[k] - (y[l]-theta.fe)/sigma.bar[l]) * (sigma[k]-sigma[l])
## })
## ## mean1 <- -theta*(sigma[k]-sigma[l])/sigma[j]^2*(1/sigma.bar[l]-1/sigma.bar[k])
## mean1 <- theta*(sigma[k]-sigma[l])*(1/sigma.bar[k]-1/sigma.bar[l])*(1+1/n/m[2]/sigma[j]^2)
## var1 <-  (sigma[k]-sigma[l])^2*(sigma[k]^2/sigma.bar[k]^2+sigma[l]^2/sigma.bar[l]^2+(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[l]))^2*(n*m[2]-1/sigma[j]^2))
## mean2 <- theta/n/m[2]/sigma[j]^2*(sigma[k]-sigma[l])*(1/sigma.bar[l]-1/sigma.bar[k])
## var2 <-(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[l]))^2*(sigma[k]-sigma[l])^2
## try <- rnorm(B1,mean=mean1,sd=sqrt(var1)) + rnorm(B1, mean=mean2,sd=sqrt(var2))    
## ## c(mean(obs),try(y,sigma,theta,j,k,l))
## qqplot(obs,try);abline(0,1,col=2)



## 20bda
require(parallel)
n <- 4
B1 <- 5e3
B2 <- 1e2
pairs <- mclapply(1:B2, mc.cores=detectCores()-3,FUN=function(yy) {
## pairs <- replicate(B2,{
    theta <- 10#runif(1,-1,1)
    sigma <- runif(n,1,2)
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    jkl <- sample(n,3,replace=FALSE); j <- jkl[1]; k <- jkl[2]; l <- jkl[3]
    ## mean1 <- theta*(sigma[k]-sigma[l])*(1/sigma.bar[k]-1/sigma.bar[l])*(1+1/n/m[2]/sigma[j]^2)
    mean1 <- theta*(sigma[k]-sigma[l])*(1/sigma.bar[k]-1/sigma.bar[l] +1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])*(n*m[2]-1/sigma[j]^2))
    var1 <-  (sigma[k]-sigma[l])^2*( sigma[k]^2*(1/sigma.bar[k]+1/sigma[k]^2/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k]))^2 + sigma[l]^2*(-1/sigma.bar[l]+1/sigma[l]^2/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k]))^2 + (1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[l]))^2*(n*m[2]-1/sigma[j]^2-1/sigma[k]^2-1/sigma[l]^2))
    var1 <- (sigma[k]-sigma[l])^2*( (sigma[k]/sigma.bar[k])^2+(sigma[l]/sigma.bar[l])^2 - 1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])^2*(1+1/n/m[2]/sigma[j]^2))
    obs <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        ## y[j]* ( y[j]/n/m[2]/sigma[j]^2*(1/sigma.bar[l]-1/sigma.bar[k])*(sigma[k]-sigma[l])   +   (y[k]/sigma.bar[k]-y[l]/sigma.bar[l]+1/n/m[2]*sum(y[-j]/sigma[-j]^2)*(1/sigma.bar[l]-1/sigma.bar[k]))*(sigma[k]-sigma[l])    >0)
        ## y[j]/n/m[2]/sigma[j]^2*(1/sigma.bar[l]-1/sigma.bar[k])*(sigma[k]-sigma[l])   +
        (y[k]/sigma.bar[k]-y[l]/sigma.bar[l]+1/n/m[2]*sum(y[-j]/sigma[-j]^2)*(1/sigma.bar[l]-1/sigma.bar[k]))*(sigma[k]-sigma[l])
        ## (sigma[k]-sigma[l])*(y[k]*(1/sigma.bar[k]+1/sigma[k]^2/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])) + y[l]*(-1/sigma.bar[l]+1/sigma[l]^2/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])) + 1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])*sum((y/sigma^2)[-c(j,k,l)]))
        ## (sigma[k]-sigma[l])*(y[k]*(1/sigma.bar[k]+1/sigma[k]^2/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])))
    })
    ## yj <- rnorm(B1,mean=theta,sd=sigma[j])
    ## try <- yj* ( yj/n/m[2]/sigma[j]^2*(1/sigma.bar[l]-1/sigma.bar[k])*(sigma[k]-sigma[l])+   rnorm(B1,mean1,sqrt(var1))   >0 )
    ## try <- yj/n/m[2]/sigma[j]^2*(1/sigma.bar[l]-1/sigma.bar[k])*(sigma[k]-sigma[l])+
    try <- rnorm(B1,mean1,sqrt(var1))   
    c(mean(obs),mean1)#mean(try))
    ##     try <- integrate(function(y) y * pnorm((y-theta*n*m[2])*beta(sigma,j,k,l))*dnorm(y,mean=theta,sd=sigma[j]),theta-5,theta+5)$val
##     ## c(mean(obs),try(y,sigma,theta,j,k,l))
##     c(mean(obs),try)
})
pairs <- simplify2array(pairs)
plot(pairs[1,],pairs[2,]); abline(0,1)

## 20be
require(parallel)
n <- 4
B1 <- 1e4
B2 <- 1e2
pairs <- mclapply(1:B2, mc.cores=detectCores()-3,FUN=function(yy) {
## pairs <- replicate(B2,{
    theta <- 8#runif(1,-1,1)
    sigma <- runif(n,1,2)*5
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    ## beta <- function(sigma,k,l) {
    ##     denom <- (sigma[k]-sigma[l])^2*(1/sigma.bar[k]^2+1/sigma.bar[l]^2+(n-1)*(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[l]))^2)
    ##     numer <- 1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])*(sigma[k]-sigma[l])
    ##     numer/sqrt(denom)
    ## }
    jkl <- sample(n,3,replace=FALSE); j <- jkl[1]; k <- jkl[2]; l <- jkl[3]
    ## mean1 <- theta*(sigma[k]-sigma[l])*(1/sigma.bar[k]-1/sigma.bar[l] +1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])*(n*m[2]-1/sigma[j]^2))
    ## mean1 <-  theta*(sigma[k]-sigma[l])*(1/sigma.bar[k]-1/sigma.bar[l])*(1/n/m[2]/sigma[j]^2)
    ## var1 <-   (sigma[k]-sigma[l])^2*( (sigma[k]/sigma.bar[k])^2+(sigma[l]/sigma.bar[l])^2 - 1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])^2*(1+1/n/m[2]/sigma[j]^2))
    ## beta <- function(sigma,j,k,l) {
    ##     numer <- (1/sigma.bar[l]-1/sigma.bar[k])*(sigma[k]-sigma[l])/n/m[2]/sigma[j]^2
    ##     denom <- var1#(sigma[k]-sigma[l])^2*(sigma[k]^2/sigma.bar[k]^2+sigma[l]^2/sigma.bar[l]^2+(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[l]))^2*(n*m[2]-1/sigma[j]^2))
    ##     numer/sqrt(denom)
    ## }
    ## try <- function(y,sigma,theta,j,k,l)
        ## sigma[j]*integrate(function(y) (y-theta)/sigma[j]* pnorm((y/sigma[j]^2)*beta(sigma,k,l))*dnorm(y,mean=theta,sd=sigma[j]),theta-5,theta+5)$val +theta/2
    obs <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        y[j]*( ((y[k]-theta.fe)/sigma.bar[k] - (y[l]-theta.fe)/sigma.bar[l]) * (sigma[k]-sigma[l]) > 0)
        ## (y[k]/sigma.bar[k]-y[l]/sigma.bar[l]+1/n/m[2]*sum(y[-j]/sigma[-j]^2)*(1/sigma.bar[l]-1/sigma.bar[k]))*(sigma[k]-sigma[l])
    })
    ## mean2 <- theta/n/m[2]/sigma[j]^2*(sigma[k]-sigma[l])*(1/sigma.bar[l]-1/sigma.bar[k])
    ## var2 <-(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[l]))^2*(sigma[k]-sigma[l])^2
    ## try <- rnorm(B1,mean=mean1,sd=sqrt(var1)) + rnorm(B1, mean=mean2,sd=sqrt(var2))
    ## try <- replicate(B1, {
    ##     y <- rnorm(n,mean=theta,sd=sigma)
    ##     ## y[j]* ( y[j]/n/m[2]/sigma[j]^2*(1/sigma.bar[l]-1/sigma.bar[k])*(sigma[k]-sigma[l])   +   (y[k]/sigma.bar[k]-y[l]/sigma.bar[l]+1/n/m[2]*sum(y[-j]/sigma[-j]^2)*(1/sigma.bar[l]-1/sigma.bar[k]))*(sigma[k]-sigma[l])    >0)
    ##     ## y[j]* ( y[j]/n/m[2]/sigma[j]^2*(1/sigma.bar[l]-1/sigma.bar[k])*(sigma[k]-sigma[l])+   rnorm(1,mean1,sqrt(var1))   >0)
    ##     j.term <- y[j]*1/n/m[2]/sigma[j]^2*(1/sigma.bar[l]-1/sigma.bar[k])*(sigma[k]-sigma[l])
    ##     y[j]* pnorm( (1/sigma.bar[l]-1/sigma.bar[k])*(sigma[k]-sigma[l])/n/m[2]/sigma[j]^2/sqrt(var1)*(y[j]-theta))
    ## })
    ## try <- yj* (yj/n/m[2]/sigma[j]^2*(sigma[k]-sigma[l])*(1/sigma.bar[l]-1/sigma.bar[k])  + rnorm(1, mean=mean2,sd=sqrt(var2)) > 0)
    ## try <- integrate(function(y) y * pnorm((y-theta)*beta(sigma,j,k,l))*dnorm(y,mean=theta,sd=sigma[j]),theta-5,theta+5)$val
    ## try <- integrate(function(y) (y-theta) * pnorm((y-theta)*beta(sigma,j,k,l))*dnorm(y,mean=theta,sd=sigma[j]),theta-5,theta+5)$val +theta/2
    try <-     sigma[j]^2/sqrt(2*pi)*((1/sigma.bar[l]-1/sigma.bar[k])/n/m[2]/sigma[j]^2*sign(sigma[k]-sigma[l])) / sqrt(   (sigma[k]/sigma.bar[k])^2+(sigma[l]/sigma.bar[l])^2 - 1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])^2  ) + theta/2
## ##     ## c(mean(obs),try(y,sigma,theta,j,k,l))
    c(mean(obs),mean(try))
##     c(mean(obs),try)
})
pairs <- simplify2array(pairs)
plot(pairs[1,],pairs[2,]); abline(0,1)


## triple sum version
require(parallel)
n <- 4
B1 <- 3e2
B2 <- 1e2
pairs <- mclapply(1:B2, mc.cores=detectCores()-3,FUN=function(yy) {
## pairs <- replicate(B2,{
    theta <- 2#runif(1,-1,1)
    sigma <- runif(n,1,2)
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    jkl <- sample(n,3,replace=FALSE); j <- jkl[1]; k <- jkl[2]; l <- jkl[3]
    ## obs <- replicate(B1, {
    ##     y[j]*( ((y[j]-theta.fe)/sigma.bar[j] - (y[k]-theta.fe)/sigma.bar[k]) * (sigma[j]-sigma[k]) > 0)
    ## })
    obs.triple.sum <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        triple.sum <- sum(sapply(1:n, function(j) {
            double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
                k <- kl[1]; l <- kl[2]
                y[j]/sigma[j]^2*( ((y[k]-theta.fe)/sigma.bar[k] - (y[l]-theta.fe)/sigma.bar[l]) * (sigma[k]-sigma[l]) > 0)
            }))
        }))
    })
    ## obs.double.sum <- sum(sapply(1:n, function(j)
    ##     sapply((1:n)[-j], function(k)
    ##         1/sigma[j]^2*y[j] * (((y[j]-theta.fe)/sigma.bar[j]-(y[k]-theta.fe)/sigma.bar[k])*(sigma[j]-sigma[k])>0)
    ##         )))
    try.triple.sum <-
        sum(sapply(1:n, function(j) {
            double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
                k <- kl[1]; l <- kl[2]
                ( sigma[j]^2/sqrt(2*pi)*((1/sigma.bar[l]-1/sigma.bar[k])/n/m[2]/sigma[j]^2*sign(sigma[k]-sigma[l])) / sqrt(   (sigma[k]/sigma.bar[k])^2+(sigma[l]/sigma.bar[l])^2 - 1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])^2  ) + theta/2 ) / sigma[j]^2
            }))
        }))
    ## c(obs.triple.sum,try.triple.sum)
    c(mean(obs.triple.sum),try.triple.sum)
})
pairs <- simplify2array(pairs)
plot(pairs[1,],pairs[2,]); abline(0,1)



dd



## 20bf full final formmula for
## cov(tau.hat,theta.hat|sigma). convergence much clearer for #1 than
## #2, though they're equivalent
source('../begg/2/misc.R')
require(parallel)
n <- 4
B1 <- 3e2
B2 <- 3e2
pairs <- mclapply(1:B2, mc.cores=detectCores()-3,FUN=function(yy) {
    theta <- 2#runif(1,-1,1)
    sigma <- runif(n,1,2)*1
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    ## y <- rnorm(n,mean=theta,sd=sigma)
    ## theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
    obs <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        ## ## theta.fe*begg.test(y,sigma^2)['stat']
        ## double.sum <- sum(sapply(1:n, function(j)
        ##     sapply((1:n)[-j], function(k)
        ##         1/sigma[j]^2*y[j] * (((y[j]-theta.fe)/sigma.bar[j]-(y[k]-theta.fe)/sigma.bar[k])*(sigma[j]-sigma[k])>0)
        ##         )))
        ## triple.sum <- sum(sapply(1:n, function(j) {
        ##     double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
        ##         k <- kl[1]; l <- kl[2]
        ##         1/sigma[j]^2*y[j]*(((y[k]-theta.fe)/sigma.bar[k]-(y[l]-theta.fe)/sigma.bar[l])*(sigma[k]-sigma[l])>0)
        ##     }))
        ## }))
        ## double.sum+triple.sum
        ## ( theta.fe*begg.test(y,sigma^2)['stat']+theta ) / (2/choose(n,2)/n/m[2]) #1
        theta.fe*begg.test(y,sigma^2)['stat'] #2
    })
    try.triple.sum <-
        sum(sapply(1:n, function(j) {
            double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
                k <- kl[1]; l <- kl[2]
                ( sigma[j]^2/sqrt(2*pi)*((1/sigma.bar[l]-1/sigma.bar[k])/n/m[2]/sigma[j]^2*sign(sigma[k]-sigma[l])) / sqrt(   (sigma[k]/sigma.bar[k])^2+(sigma[l]/sigma.bar[l])^2 - 1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])^2  ) + theta/2 ) / sigma[j]^2
            }))
        }))
    try.double.sum <- sum(sapply(1:n, function(j)
        sapply((1:n)[-j], function(k) {
            numer <- (1/sigma.bar[j]+1/sigma[j]^2/n/m[2]*(1/sigma.bar[k] -1/sigma.bar[j]))
            ## denom <- sigma[j]^2/sigma.bar[j]^2+sigma[k]^2/sigma.bar[k]^2+(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[j]))^2*(1/sigma[j]^2+1/sigma[k]^2+n-2)-2/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[j])^2 #wrong formula
            denom <- (sigma[j]/sigma.bar[j])^2+(sigma[k]/sigma.bar[k])^2-1/n/m[2]*(1/sigma.bar[j]-1/sigma.bar[k])^2
            1/sqrt(2*pi)*sign(sigma[j]-sigma[k])*numer/sqrt(denom) +theta/2/sigma[j]^2
        })))
    try <- try.double.sum+try.triple.sum
    ## try <- 2/choose(n,2)/n/m[2]*try - theta #2
    c(try,mean(obs))
})
pairs <- simplify2array(pairs)
plot(pairs[1,],pairs[2,]); abline(0,1)


## 20bg same as 20bf, using d and e notation (notes B-III)
source('../begg/2/misc.R')
require(parallel)
n <- 4
B1 <- 3e2
B2 <- 3e2
pairs <- mclapply(1:B2, mc.cores=detectCores()-3,FUN=function(yy) {
    theta <- 2#runif(1,-1,1)
    sigma <- runif(n,1,2)*8
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    ## y <- rnorm(n,mean=theta,sd=sigma)
    ## theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
    obs <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        theta.fe*begg.test(y,sigma^2)['stat'] 
    })
    triple.sum <-
        sum(sapply(1:n, function(j) {
            double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
                k <- kl[1]; l <- kl[2]
                numer <- 1/sqrt(2*pi)*((1/sigma.bar[l]-1/sigma.bar[k])/n/m[2]/sigma[j]^2*sign(sigma[k]-sigma[l]))
                denom <-   (sigma[k]/sigma.bar[k])^2+(sigma[l]/sigma.bar[l])^2 - 1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])^2   # e_{jk}
                2/choose(n,2)/n/m[2]* numer / sqrt(denom)
            }))
        }))
    ## double.sum <- sum(sapply(1:n, function(j)
    ##     sapply((1:n)[-j], function(k) {
    ##         numer <-     1/sqrt(2*pi)*sign(sigma[j]-sigma[k])*(1/sigma.bar[j]+1/sigma[j]^2/n/m[2]*(1/sigma.bar[k] -1/sigma.bar[j]))
    ##         denom <- sigma[j]^2/sigma.bar[j]^2+sigma[k]^2/sigma.bar[k]^2+(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[j]))^2*(1/sigma[j]^2+1/sigma[k]^2+n-2)-2/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[j])^2 # d_{jk}
    ##     numer / sqrt(denom) 
    ##     })))
    double.sum <- sum(apply(combn((1:n),2),2, function(jk) {
        j <- jk[1]; k <- jk[2]
        ## numer <-     1/sqrt(2*pi)*sign(sigma[j]-sigma[k])*(1/sigma.bar[j]+1/sigma[j]^2/n/m[2]*(1/sigma.bar[k] -1/sigma.bar[j]))
        numer <- 1/sqrt(2*pi)*sign(sigma[j]-sigma[k])*(1/sigma.bar[j]-1/sigma.bar[k]+1/n/m[2]*(1/sigma.bar[k] -1/sigma.bar[j])*(1/sigma[j]^2+1/sigma[k]^2))
        ## denom <- sigma[j]^2/sigma.bar[j]^2+sigma[k]^2/sigma.bar[k]^2+(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[j]))^2*(1/sigma[j]^2+1/sigma[k]^2+n-2)-2/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[j])^2 # d_{jk}
        denom <- (sigma[j]/sigma.bar[j])^2+(sigma[k]/sigma.bar[k])^2-1/n/m[2]*(1/sigma.bar[j]-1/sigma.bar[k])^2
        2/choose(n,2)/n/m[2]* numer / sqrt(denom) 
        }))
    try <- double.sum+triple.sum
    c(try,mean(obs))
})
pairs <- simplify2array(pairs)
plot(pairs[1,],pairs[2,]); abline(0,1)


## 20bh growth rates. cancellation between double mean and triple mean
## terms. difference is converging to 0.
require(parallel)
## n <- 4
B <- 5e3
ns <- round(seq(5,1e2,len=50))
by.n <- mclapply(ns, mc.cores=detectCores()-3, FUN=function(n) {
## by.n <- sapply(ns, function(n) {
    theta <- 2#runif(1,-1,1)
    ## mean(replicate(B, {
    sigma <- runif(n,1,2)*8
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    ## mean(replicate(B, {
    ##     y <- rnorm(n,mean=theta,sd=sigma)
    ##     theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
    ##     theta.fe*begg.test(y,sigma^2)['stat'] 
    ## }))
    triple.sum <-
        sum(sapply(1:n, function(j) {
            double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
                k <- kl[1]; l <- kl[2]
                numer <- 1/sqrt(2*pi)*((1/sigma.bar[l]-1/sigma.bar[k])/n/m[2]/sigma[j]^2*sign(sigma[k]-sigma[l]))
                denom <-   (sigma[k]/sigma.bar[k])^2+(sigma[l]/sigma.bar[l])^2 - 1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])^2   # e_{jk}
                2/choose(n,2)/n/m[2]* numer / sqrt(denom)
            }))
        }))
    double.sum <- sum(apply(combn((1:n),2),2, function(jk) {
        j <- jk[1]; k <- jk[2]
        ## numer <-     1/sqrt(2*pi)*sign(sigma[j]-sigma[k])*(1/sigma.bar[j]+1/sigma[j]^2/n/m[2]*(1/sigma.bar[k] -1/sigma.bar[j]))
        numer <- 1/sqrt(2*pi)*sign(sigma[j]-sigma[k])*(1/sigma.bar[j]-1/sigma.bar[k]+1/n/m[2]*(1/sigma.bar[k] -1/sigma.bar[j])*(1/sigma[j]^2+1/sigma[k]^2))
        ## denom <- sigma[j]^2/sigma.bar[j]^2+sigma[k]^2/sigma.bar[k]^2+(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[j]))^2*(1/sigma[j]^2+1/sigma[k]^2+n-2)-2/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[j])^2 # d_{jk}
        denom <- (sigma[j]/sigma.bar[j])^2+(sigma[k]/sigma.bar[k])^2-1/n/m[2]*(1/sigma.bar[j]-1/sigma.bar[k])^2
        2/choose(n,2)/n/m[2]* numer / sqrt(denom) 
    }))
    c(double.sum=double.sum,triple.sum=triple.sum)
})
by.n <- simplify2array(by.n)
plot(ns,ns*by.n['double.sum',],ylim=range(ns*t(by.n)),type='l')
lines(ns,ns*by.n['triple.sum',])
lines(ns,ns*(by.n['double.sum',]+by.n['triple.sum',]))
abline(h=0)

dd

## 20ca
require(parallel)
source('misc.R')
n <- 4
B1 <- 5e2
B2 <- 2e2
pairs <- replicate(B2, {
    theta <- runif(1,-1,1)
    sigma <- runif(n,1,2)
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    alpha <- function(sigma,j,k) {
        numer <- (1/sigma.bar[j]-1/sigma[j]^2/n/m[2]*(1/sigma.bar[j] -1/sigma.bar[k]))*(sigma[j]-sigma[k])
        denom <- ( (-sigma[k]/sigma.bar[k]+1/n/m[2]/sigma[k]*(1/sigma.bar[k]-1/sigma.bar[j]))^2 + (n-2)*(1/n/m[2]*(1/sigma.bar[j]-1/sigma.bar[k]))^2 ) *(sigma[j]-sigma[k])^2
        numer/sqrt(denom)
    }
    obs <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        double.sum <- sum(apply(combn(n,2),2, function(kl) {
            k <- kl[1]; l <- kl[2]
            ((y[k]-theta.fe)/sigma.bar[k]-(y[l]-theta.fe)/sigma.bar[l])*(sigma[k]-sigma[l])>0
        }))
        sum(y/sigma^2)*double.sum                   
    })
    try <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        double.sum <- sapply(1:n, function(j)
            sapply((1:n)[-j], function(k)
                -1/sigma[j]^2*((y[j]-theta)*pnorm((y[j]-theta)*alpha(sigma,j,k))+theta/2)
                ))
        double.sum + m[2]*theta/4*n*(n-1)*(n-2)
    })
    c(mean(obs),mean(try))
})
plot(pairs[1,],pairs[2,]); abline(0,1)







## 20caa
require(parallel)
source('misc.R')
n <- 4
B1 <- 5e2
B2 <- 5e2
pairs <- replicate(B2, {
    theta <- runif(1,-1,1)
    sigma <- runif(n,1,2)
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    alpha <- function(sigma,j,k) {
        numer <- (1/sigma.bar[j]-1/sigma[j]^2/n/m[2]*(1/sigma.bar[j] -1/sigma.bar[k]))*(sigma[j]-sigma[k])
        denom <- ( (-sigma[k]/sigma.bar[k]+1/n/m[2]/sigma[k]*(1/sigma.bar[k]-1/sigma.bar[j]))^2 + (n-2)*(1/n/m[2]*(1/sigma.bar[j]-1/sigma.bar[k]))^2 ) *(sigma[j]-sigma[k])^2
        numer/sqrt(denom)
    }
    ## try <- replicate(B1, {
    ##     y <- rnorm(n,mean=theta,sd=sigma)
    ##     theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
    ##     double.sum <- sum(apply(combn(n,2),2, function(kl) {
    ##         k <- kl[1]; l <- kl[2]
    ##         ((y[k]-theta.fe)/sigma.bar[k]-(y[l]-theta.fe)/sigma.bar[l])*(sigma[k]-sigma[l])>0
    ##     }))
    ##     2/choose(n,2)/m[2]*sum(y/sigma^2)*double.sum-theta
    ## })
    obs <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        c(theta.fe=theta.fe,tau=unname(begg.test(y,sigma^2)['stat']))
    })
    ## c(cov(obs['theta.fe',],obs['tau',]),mean(try))
    c(cov(obs['theta.fe',],obs['tau',]),mean(obs['theta.fe',]*obs['tau',]))
           ## c(mean(obs),mean(try))
})
plot(pairs[1,],pairs[2,]); abline(0,1)



















## ## 20cb
## require(parallel)
## source('misc.R')
## n <- 4
## B1 <- 1e3
## B2 <- 1e2
## pairs <- mclapply(1:B2, mc.cores=detectCores()-3,FUN=function(yy) {
## ## pairs <- replicate(B2, {
##     theta <- runif(1,-1,1)
##     sigma <- runif(n,1,2)/10
##     sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
##     m <- sapply(1:2,function(k)mean((1/sigma)^k))
##     alpha <- function(sigma,j,k) {
##         numer <- (1/sigma.bar[j]-1/sigma[j]^2/n/m[2]*(1/sigma.bar[j] -1/sigma.bar[k]))*(sigma[j]-sigma[k])
##         denom <- ( (-sigma[k]/sigma.bar[k]+1/n/m[2]/sigma[k]*(1/sigma.bar[k]-1/sigma.bar[j]))^2 + (n-2)*(1/n/m[2]*(1/sigma.bar[j]-1/sigma.bar[k]))^2 ) *(sigma[j]-sigma[k])^2
##         numer/sqrt(denom)
##     }
##     obs <- replicate(B1, {
##         y <- rnorm(n,mean=theta,sd=sigma)
##         theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
##         theta.fe*begg.test(y,sigma^2)['stat']
##     })
##     try <- replicate(B1, {
##         y <- rnorm(n,mean=theta,sd=sigma)
##         theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
##         theta.fe*begg.test(y,sigma^2)['stat']
##     })
##     ## try <- replicate(B1, {
##     ##     y <- rnorm(n,mean=theta,sd=sigma)
##     ##     theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
##     ##     double.sum <- sum(apply(combn(n,2),2, function(kl) {
##     ##         k <- kl[1]; l <- kl[2]
##     ##         ((y[k]-theta.fe)/sigma.bar[k]-(y[l]-theta.fe)/sigma.bar[l])*(sigma[k]-sigma[l])>0
##     ##     }))
##     ##     2/choose(n,2)*double.sum*theta.fe-theta
##     ## })
##     c(mean(obs), mean(try))
##            ## c(mean(obs),mean(try))
## })
## pairs <- simplify2array(pairs)
## plot(pairs[1,],pairs[2,]); abline(0,1)















## 20cc
require(parallel)
source('../begg/2/misc.R')
n <- 4
B1 <- 1e2
B2 <- 1e2
pairs <- mclapply(1:B2, mc.cores=detectCores()-3,FUN=function(yy) {
    ## pairs <- replicate(B2, {
    theta <- 1#runif(1,-1,1)
    sigma <- runif(n,1,2)
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    alpha <- function(sigma,j,k) {
        numer <- (1/sigma.bar[j]-1/sigma[j]^2/n/m[2]*(1/sigma.bar[j] -1/sigma.bar[k]))*(sigma[j]-sigma[k])
        denom <- ( (-sigma[k]/sigma.bar[k]+1/n/m[2]/sigma[k]*(1/sigma.bar[k]-1/sigma.bar[j]))^2 + (n-2)*(1/n/m[2]*(1/sigma.bar[j]-1/sigma.bar[k]))^2 ) *(sigma[j]-sigma[k])^2
        numer/sqrt(denom)
    }
    obs <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        double.sum <- sum(sapply(1:n, function(j)
            sapply((1:n)[-j], function(k)
                1/sigma[j]^2*y[j] * (((y[j]-theta.fe)/sigma.bar[j]-(y[k]-theta.fe)/sigma.bar[k])*(sigma[j]-sigma[k])>0)
                )))
        double.sum #+ m[2]*theta/4*n*(n-1)*(n-2)
    })
    try <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        double.sum <- sum(sapply(1:n, function(j)
            sapply((1:n)[-j], function(k)
                1/sigma[j]^2*((y[j]-theta)*pnorm((y[j]-theta)*alpha(sigma,j,k))+theta/2)
                )))
        double.sum #+ m[2]*theta/4*n*(n-1)*(n-2)
    })
    c(mean(obs),mean(try))
})
pairs <- simplify2array(pairs)
plot(pairs[1,],pairs[2,]); abline(0,1)



## 20cd corrected in 20cdd
require(parallel)
source('../begg/2/misc.R')
n <- 4
B1 <- 1e3
B2 <- 1e2
pairs <- mclapply(1:B2, mc.cores=detectCores()-3,FUN=function(yy) {
    ## pairs <- replicate(B2, {
    theta <- 1#runif(1,-1,1)/1e2 # theta variation was dominating
    sigma <- runif(n,1,2)*5
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    obs <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        triple.sum <- sum(sapply(1:n, function(j) {
            double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
                k <- kl[1]; l <- kl[2]
                ((y[k]-theta.fe)/sigma.bar[k]-(y[l]-theta.fe)/sigma.bar[l])*(sigma[k]-sigma[l])>0
            }))
            1/sigma[j]^2*y[j]*double.sum
        }))
        triple.sum 
    })
    c(mean(obs),try=theta*n*m[2]*(n-1)*(n-2)/4)
})
pairs <- simplify2array(pairs)
plot(pairs[1,],pairs[2,]); abline(0,1)




## 20cdd after realizing y_j dependence hidden in theta.hat
require(parallel)
source('../begg/2/misc.R')
n <- 4
B1 <- 3e3
B2 <- 1e2
pairs <- mclapply(1:B2, mc.cores=detectCores()-3,FUN=function(yy) {
    ## pairs <- replicate(B2, {
    theta <- 1#runif(1,-1,1)/1e2 # theta variation was dominating
    sigma <- runif(n,1,2)*5
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    beta <- function(sigma,k,l) {
        denom <- (sigma[k]-sigma[l])^2*(1/sigma.bar[k]^2+1/sigma.bar[l]^2+(n-1)*(1/n/m[2]*(1/sigma.bar[k]-1/sigma.bar[l]))^2)
        numer <- 1/n/m[2]*(1/sigma.bar[l]-1/sigma.bar[k])*(sigma[k]-sigma[l])
        numer/sqrt(denom)
    }
    obs <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        triple.sum <- sum(sapply(1:n, function(j) {
            double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
                k <- kl[1]; l <- kl[2]
                ((y[k]-theta.fe)/sigma.bar[k]-(y[l]-theta.fe)/sigma.bar[l])*(sigma[k]-sigma[l])>0
            }))
            1/sigma[j]^2*y[j]*double.sum
        }))
        triple.sum 
    })
    try <- sum(sapply(1:n, function(j) {
        double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
            k <- kl[1]; l <- kl[2]
            b <- beta(sigma,k,l)           
            ## 1/sigma[j]^2*  dnorm(theta,sd=sigma[j]*sqrt(b^2+sigma[j]^2)/b) + theta/2
            1/sigma[j]^2*(theta/2 + 1/sqrt(2*pi)*sigma[j]*b/sqrt(b^2+sigma[j]^2)*exp(-1/2*theta^2*b^2/(b^2+sigma[j]^2)/sigma[j]^2))
        }))
    }))
    c(mean(obs),try=try)
})
pairs <- simplify2array(pairs)
plot(pairs[1,],pairs[2,]); abline(0,1)

b <- beta(sigma,k,l)           
dnorm(theta,sd=sigma[j]*sqrt(b^2+sigma[j]^2)/b)
1/sqrt(2*pi*sigma[j]^2)*b/sqrt(b^2+sigma[j]^2)*exp(-1/2*theta^2*b^2/(b^2+sigma[j]^2)/sigma[j]^2)

dd

    ##  sum(apply(combn((1:n),3),2, function(jkl) {
## +          j <- jkl[1]; k <- jkl[2]; l <- jkl[3]
## +          1/sigma[j]^2*y[j] * (((y[k]-theta.fe)/sigma.bar[k]-(y[l]-theta.fe)/sigma.bar[l])*(sigma[k]-sigma[l])>0)
## +      }))
## + 
##     +
    
## 20ce
source('../begg/2/misc.R')
theta <- runif(1,-1,1)
sigma <- runif(n,1,2)*1
sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
m <- sapply(1:2,function(k)mean((1/sigma)^k))
y <- rnorm(n,mean=theta,sd=sigma)
theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
## double.sum <- sum(sapply(1:n, function(j)
##     sapply((1:n)[-j], function(k)
##         1/sigma[j]^2*y[j] * (((y[j]-theta.fe)/sigma.bar[j]-(y[k]-theta.fe)/sigma.bar[k])*(sigma[j]-sigma[k])>0)
##         )))
## triple.sum <- sum(sapply(1:n, function(j) {
##     double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
##         k <- kl[1]; l <- kl[2]
##         ((y[k]-theta.fe)/sigma.bar[k]-(y[l]-theta.fe)/sigma.bar[l])*(sigma[k]-sigma[l])>0
##     }))
##     1/sigma[j]^2*y[j]*double.sum
## }))
double.sum <- sum(sapply(1:n, function(j)
    sapply((1:n)[-j], function(k)
        1/sigma[j]^2*y[j] * (((y[j]-theta.fe)/sigma.bar[j]-(y[k]-theta.fe)/sigma.bar[k])*(sigma[j]-sigma[k])>0)
        )))
triple.sum <- sum(sapply(1:n, function(j) {
    double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
        k <- kl[1]; l <- kl[2]
        ((y[k]-theta.fe)/sigma.bar[k]-(y[l]-theta.fe)/sigma.bar[l])*(sigma[k]-sigma[l])>0
    }))
    1/sigma[j]^2*y[j]*double.sum
}))
## double.sum+triple.sum
2/choose(n,2)/n/m[2]*(double.sum+triple.sum)
theta.fe*begg.test(y,sigma^2)['stat']+theta.fe

## double.sum <- sum(apply(combn(n,2),2, function(kl) {
##     k <- kl[1]; l <- kl[2]
##     ((y[k]-theta.fe)/sigma.bar[k]-(y[l]-theta.fe)/sigma.bar[l])*(sigma[k]-sigma[l])>0
## }))
## ## sum(y/sigma^2)*double.sum
## 2/choose(n,2)/n/m[2]*sum(y/sigma^2)*double.sum




## 20cf
require(parallel)
source('../begg/2/misc.R')
n <- 4
B1 <- 1e2
B2 <- 1e2
pairs <- replicate(B2, {
    theta <- 2#runif(1,-1,1)
    sigma <- runif(n,1,2)
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    alpha <- function(sigma,j,k) {
        numer <- (1/sigma.bar[j]-1/sigma[j]^2/n/m[2]*(1/sigma.bar[j] -1/sigma.bar[k]))*(sigma[j]-sigma[k])
        denom <- ( (-sigma[k]/sigma.bar[k]+1/n/m[2]/sigma[k]*(1/sigma.bar[k]-1/sigma.bar[j]))^2 + (n-2)*(1/n/m[2]*(1/sigma.bar[j]-1/sigma.bar[k]))^2 ) *(sigma[j]-sigma[k])^2
        numer/sqrt(denom)
    }
    obs <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        theta.fe*begg.test(y,sigma^2)['stat']+theta
        ## double.sum <- sum(sapply(1:n, function(j)
        ##     sapply((1:n)[-j], function(k)
        ##         1/sigma[j]^2*y[j] * (((y[j]-theta.fe)/sigma.bar[j]-(y[k]-theta.fe)/sigma.bar[k])*(sigma[j]-sigma[k])>0)
        ##         )))
        ## double.sum #+ m[2]*theta/4*n*(n-1)*(n-2)
    })
    try <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        ## double.sum <- sum(sapply(1:n, function(j)
        ##     sapply((1:n)[-j], function(k)
        ##         1/sigma[j]^2*y[j] * (((y[j]-theta.fe)/sigma.bar[j]-(y[k]-theta.fe)/sigma.bar[k])*(sigma[j]-sigma[k])>0)
        ##         )))
        ## triple.sum <- sum(sapply(1:n, function(j) {
        ##     double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
        ##         k <- kl[1]; l <- kl[2]
        ##         ((y[k]-theta.fe)/sigma.bar[k]-(y[l]-theta.fe)/sigma.bar[l])*(sigma[k]-sigma[l])>0
        ##     }))
        ##     1/sigma[j]^2*y[j]*double.sum
        ## }))
        ## double.sum <- sum(sapply(1:n, function(j)
        ##     sapply((1:n)[-j], function(k)
        ##         1/sigma[j]^2*((y[j]-theta)*pnorm((y[j]-theta)*alpha(sigma,j,k))+theta/2)
        ##         )))
        double.sum <- sum(sapply(1:n, function(j)
            sapply((1:n)[-j], function(k)
                1/sigma[j]^2*((y[j]-theta)*pnorm((y[j]-theta)*alpha(sigma,j,k)))
                )))
        triple.sum <- theta*n*m[2]*(n-1)*(n-2)/4
        ## double.sum+triple.sum
        ## 2/choose(n,2)/n/m[2]*(double.sum+triple.sum)
        ## 2/choose(n,2)/n/m[2]*double.sum + (n-2)/n*theta
        2/choose(n,2)/n/m[2]*double.sum + theta
    })
    c(mean(obs),mean(try))
})
plot(pairs[1,],pairs[2,]); abline(0,1)


dd

## ## 20cg
## require(parallel)
## source('../begg/2/misc.R')
## n <- 4
## B1 <- 1e2
## B2 <- 2e2
## ## pairs <- replicate(B2, {
## pairs <- mclapply(1:B2, mc.cores=detectCores()-3,FUN=function(yy) {
##     theta <- runif(1,-1,1)
##     sigma <- runif(n,1,2)*10
##     sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
##     m <- sapply(1:2,function(k)mean((1/sigma)^k))
##     alpha <- function(sigma,j,k) {
##         numer <- (1/sigma.bar[j]-1/sigma[j]^2/n/m[2]*(1/sigma.bar[j] -1/sigma.bar[k]))*(sigma[j]-sigma[k])
##         denom <- ( (-sigma[k]/sigma.bar[k]+1/n/m[2]/sigma[k]*(1/sigma.bar[k]-1/sigma.bar[j]))^2 + (n-2)*(1/n/m[2]*(1/sigma.bar[j]-1/sigma.bar[k]))^2 ) *(sigma[j]-sigma[k])^2
##         numer/sqrt(denom)
##     }
##     obs <- replicate(B1, {
##         y <- rnorm(n,mean=theta,sd=sigma)
##         theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
##         theta.fe*begg.test(y,sigma^2)['stat']+theta
##         ## double.sum <- sum(sapply(1:n, function(j)
##         ##     sapply((1:n)[-j], function(k)
##         ##         1/sigma[j]^2*y[j] * (((y[j]-theta.fe)/sigma.bar[j]-(y[k]-theta.fe)/sigma.bar[k])*(sigma[j]-sigma[k])>0)
##         ##         )))
##         ## double.sum #+ m[2]*theta/4*n*(n-1)*(n-2)
##     })
##     try <- replicate(B1, {
##         y <- rnorm(n,mean=theta,sd=sigma)
##         theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
##         ## double.sum <- sum(sapply(1:n, function(j)
##         ##     sapply((1:n)[-j], function(k)
##         ##         1/sigma[j]^2*y[j] * (((y[j]-theta.fe)/sigma.bar[j]-(y[k]-theta.fe)/sigma.bar[k])*(sigma[j]-sigma[k])>0)
##         ##         )))
##         ## triple.sum <- sum(sapply(1:n, function(j) {
##         ##     double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
##         ##         k <- kl[1]; l <- kl[2]
##         ##         ((y[k]-theta.fe)/sigma.bar[k]-(y[l]-theta.fe)/sigma.bar[l])*(sigma[k]-sigma[l])>0
##         ##     }))
##         ##     1/sigma[j]^2*y[j]*double.sum
##         ## }))
##         double.sum <- sum(sapply(1:n, function(j)
##             sapply((1:n)[-j], function(k)
##                 1/sigma[j]^2*((y[j]-theta)*pnorm((y[j]-theta)*alpha(sigma,j,k))+theta/2)
##                 )))
##         ## double.sum <- sum(sapply(1:n, function(j)
##             ## sapply((1:n)[-j], function(k) {
##                 ## 1/sigma[j]^2*((y[j]-theta)*pnorm((y[j]-theta)*alpha(sigma,j,k)))
##                 ## alpha <- alpha(sigma,j,k)
##                 ## alpha/sqrt(1+sigma[j]^2*alpha^2)
##             ## }
##             ## )))
##         ## triple.sum <- theta*n*m[2]*(n-1)*(n-2)/4
##         ## double.sum+triple.sum
##         ## 2/choose(n,2)/n/m[2]*(double.sum+triple.sum)
##         2/choose(n,2)/n/m[2]*double.sum + (n-2)/n*theta
##         ## 2/choose(n,2)/n/m[2]*double.sum# + theta
##     })
##     c(mean(obs),mean(try))
## })
## pairs <- simplify2array(pairs)
## plot(pairs[1,],pairs[2,]); abline(0,1)









require(parallel)
source('../begg/2/misc.R')
n <- 4
B1 <- 5e2
B2 <- 1e2
pairs <- mclapply(1:B2, mc.cores=detectCores()-3,FUN=function(yy) {
    ## pairs <- replicate(B2, {
    theta <- runif(1,-1,1)
    sigma <- runif(n,1,2)
    sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
    m <- sapply(1:2,function(k)mean((1/sigma)^k))
    alpha <- function(sigma,j,k) {
        numer <- (1/sigma.bar[j]-1/sigma[j]^2/n/m[2]*(1/sigma.bar[j] -1/sigma.bar[k]))*(sigma[j]-sigma[k])
        denom <- ( (-sigma[k]/sigma.bar[k]+1/n/m[2]/sigma[k]*(1/sigma.bar[k]-1/sigma.bar[j]))^2 + (n-2)*(1/n/m[2]*(1/sigma.bar[j]-1/sigma.bar[k]))^2 ) *(sigma[j]-sigma[k])^2
        numer/sqrt(denom)
    }
    obs <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        theta.fe*begg.test(y,sigma^2)['stat']
        ## double.sum <- sum(sapply(1:n, function(j)
        ##     sapply((1:n)[-j], function(k)
        ##         1/sigma[j]^2*y[j] * (((y[j]-theta.fe)/sigma.bar[j]-(y[k]-theta.fe)/sigma.bar[k])*(sigma[j]-sigma[k])>0)
        ##         )))
        ## double.sum #+ m[2]*theta/4*n*(n-1)*(n-2)
    })
    try <- replicate(B1, {
        y <- rnorm(n,mean=theta,sd=sigma)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        ## double.sum <- sum(sapply(1:n, function(j)
        ##     sapply((1:n)[-j], function(k)
        ##         1/sigma[j]^2*y[j] * (((y[j]-theta.fe)/sigma.bar[j]-(y[k]-theta.fe)/sigma.bar[k])*(sigma[j]-sigma[k])>0)
        ##         )))
        ## triple.sum <- sum(sapply(1:n, function(j) {
        ##     double.sum <- sum(apply(combn((1:n)[-j],2),2, function(kl) {
        ##         k <- kl[1]; l <- kl[2]
        ##         ((y[k]-theta.fe)/sigma.bar[k]-(y[l]-theta.fe)/sigma.bar[l])*(sigma[k]-sigma[l])>0
        ##     }))
        ##     1/sigma[j]^2*y[j]*double.sum
        ## }))
        ## double.sum <- sum(sapply(1:n, function(j)
        ##     sapply((1:n)[-j], function(k)
        ##         1/sigma[j]^2*((y[j]-theta)*pnorm((y[j]-theta)*alpha(sigma,j,k))+theta/2)
        ##         )))
        double.sum <- sum(sapply(1:n, function(j)
            sapply((1:n)[-j], function(k)
                1/sigma[j]^2*((y[j]-theta)*pnorm((y[j]-theta)*alpha(sigma,j,k)))
                )))
        triple.sum <- theta*n*m[2]*(n-1)*(n-2)/4
        ## double.sum+triple.sum
        ## 2/choose(n,2)/n/m[2]*(double.sum+triple.sum)
        ## 2/choose(n,2)/n/m[2]*double.sum + (n-2)/n*theta
        2/choose(n,2)/n/m[2]*double.sum 
    })
    c(mean(obs),mean(try))
})
pairs <- simplify2array(pairs)
plot(pairs[1,],pairs[2,]); abline(0,1)





alpha <- runif(1); beta <- runif(1)
## integrate(function(z)z^2*dnorm(alpha*z+beta)*dnorm(z),-5,5)$val
## 1/(2*pi)*integrate(function(z)z^2*exp(-1/2*((alpha^2+1)*(z+alpha*beta/(alpha^2+1))^2+beta^2-alpha^2*beta^2/(alpha^2+1))),-5,5)$val
## 1/(2*pi)*integrate(function(z)z^2*exp(-1/2*((alpha^2+1)*z^2+2*alpha*beta*z+beta^2)),-5,5)$val
## 1/sqrt(2*pi)*exp(-beta^2/2/(alpha^2+1))*(alpha^2+1)^(-1/2)*(1/(alpha^2+1)+(alpha*beta/(alpha^2+1))^2)
## 1/sqrt(2*pi)*exp(-beta^2/2/(alpha^2+1))*(alpha^2*(beta^2+1)+1)*(alpha^2+1)^(-5/2)
integrate(function(z)z*dnorm(alpha*z+beta)*dnorm(z),-5,5)$val
1/sqrt(2*pi)*exp(-beta^2/2/(alpha^2+1))*(alpha^2+1)^(-1/2)*(-alpha*beta)/(alpha^2+1)
integrate(function(z)z*pnorm(alpha*z+beta)*dnorm(z),-5,5)$val
1/sqrt(2*pi)*alpha/sqrt(alpha^2+1)*exp(-beta^2/2/(alpha^2+1))




n <- 3
replicate(1e2, {
sigma <- abs(rcauchy(n))#runif(n)
u <- rnorm(n)
sigma.bar <- sqrt(sigma^2 - 1/sum(1/sigma^2))
diffs <- apply(combn((1:n),2),2, function(jk) {
    j <- jk[1]; k <- jk[2]
    (u[j]*sigma.bar[k]-u[k]*sigma.bar[j])*(u[j]*sigma[k]-u[k]*sigma[j])
})
## if(sum(diffs<0)>0)browser()
sum(diffs<0)
})


z <- rnorm(n)
sigma <- runif(n)
theta.hat <- runif(1)
diffs1 <- apply(combn((1:n),2),2, function(jk) {
    j <- jk[1]; k <- jk[2]
    ((y[j]-theta.hat)/sigma[j] - (y[k]-theta.hat)/sigma[k])*(sigma[j]-sigma[k]) > 0
})
s <- 1/sigma
z <- y/sigma
diffs2 <- apply(combn((1:n),2),2, function(jk) {
    j <- jk[1]; k <- jk[2]
   (z[j]-z[k])/(s[j]-s[k]) < theta.hat
})
mean(diffs1);mean(diffs2)



## n <- 3
ns <- round(seq(3,5e1,len=20))
a <- -.6
by.n <- lapply(ns, function(n) {
    ustats <- replicate(1e2, {
        sigma <- 1/(runif(n)/10)
        ## sigma <- abs(rcauchy(n))
        sigma.theta <- sqrt(1/sum(1/sigma^2))
        sigma.bar <- sqrt(sigma^2 - sigma.theta^2)
        ## y <- rnorm(n,mean=theta,sd=sigma)
        y <- rgamma(n,.01,10)#rbeta(n,.1,4)
        theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        u <- y-theta.fe
        ## diffs <- apply(combn((1:n),2),2, function(jk) {
        ##     j <- jk[1]; k <- jk[2]
        ##     ## (u[j]/sigma.bar[j]-u[k]/sigma.bar[k])*(u[j]/sigma[j]-u[k]/sigma[k]) < 0
        ##     abs(u[j]/sigma[j]-u[k]/sigma[k])/(abs(u[j]/sigma[j])+abs(u[k]/sigma[k]))  < n^a
        ##     ## 1-sqrt(1-max(sigma.theta/sigma[j],sigma.theta/sigma[k])^2) > n^a
        ##     ## 1-max(sigma.theta/sigma[j],sigma.theta/sigma[k])^2 < (1- n^a)^2  
        ##     ## abs(u[j]/sigma[j]-u[k]/sigma[k])/(abs(u[j]/sigma[j])+abs(u[k]/sigma[k]))  < 1-sqrt(1-max(sigma.theta/sigma[j],sigma.theta/sigma[k])^2) 
        ##     ## sum(1/sigma^2) < 1/(2*n^a - n^(2*a))
        ## })
        ## if(sum(diffs<0)>0)browser()
        ## diffs <-  1-(sigma.theta/sigma)^2 < (1-n^a)^2
        ## diffs <-  (sigma.theta/sigma)^2 > 1-(1-n^a)^2#-n^(2*a)+2^n^a
        ## diffs <- sigma[1]*sum(1/sigma[-1]^2) < 1/(2*n^a-n^(2*a))
        ## diffs <- sum(1/sigma[-1]^2) < 1/(2*n^a-n^(2*a))
        ## diffs <-  1-sqrt(1-(sigma.theta/sigma)^2) > n^a
        ## mean(diffs)
        mean(1/sum(1/sigma^2))*n
    })
    mean(ustats^2)
})
by.n <- simplify2array(by.n)
plot(ns,ns^1*by.n)



## 1/n convergence of sigma.1^2/sigma.theta^2
ns <- round(seq(1e3,5e3,len=50))
a <- -.9
by.n <- lapply(ns, function(n) {
    stats <- replicate(1e3, {
        sigma <- 1/(runif(n)/100)
        s <- 1/sigma
        sigma.theta <- sqrt(1/sum(1/sigma^2))
        ## n^(abs(a)) * 1/sum(1/sigma^2)
        ## 1-sqrt(1-(sigma.theta/sigma)^2) > n^a
        ## sigma[1]^2*sum(1/sigma^2) < n^(-a)
        ## s[1]^2/sum(s[-1]^2) > n^a
        ## abs(1 - sqrt(1-(sigma.theta/sigma[1])^2) / sqrt(1-(sigma.theta/sigma[2])^2))
        abs(1 - sqrt(1-(sigma.theta/sigma[1])^2) / sqrt(1-(sigma.theta/sigma[2])^2)) > 1/n^(.9)
    })
    mean(stats)
})
by.n <- simplify2array(by.n)
plot(ns,ns^(.5)*by.n)





ns <- round(seq(1e3,5e3,len=50))
## a <- -.9
theta <- 1
by.n <- lapply(ns, function(n) {
    stats <- replicate(5e2, {
        ## sigma <- 1/(runif(n)/100)
        ## y <- rnorm(n,mean=theta,sd=sigma)
        ## ## s <- 1/sigma
        ## sigma.theta <- sqrt(1/sum(1/sigma^2))
        ## theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
        ## ((y[1]-theta.fe)/sigma[1])/((y[2]-theta.fe)/sigma[2])
        
    })
    abs(mean(stats<1)-pcauchy(1))
})
by.n <- simplify2array(by.n)
plot(ns,ns^(1)*by.n)


B <- 1e5
rZ <- rcauchy
rZ <- function(n)rgamma(n,.01,10)
rZ <- runif
z1 <- rZ(B); z2 <- rZ(B)
hist(abs(z1-z2)/(abs(z1)+abs(z2)),prob=TRUE)
lines(density(abs(z1-z2)/(abs(z1)+abs(z2))),col=2)

B <- 1e5
n <- 1e2
a <- 1
## rZ <- rcauchy
## rZ <- function(n)rgamma(n,.01,10)
rZ <- runif
rZ <- function(n)rbeta(n,1,4)
z1 <- rZ(B); z2 <- rZ(B)
## hist(abs(z1-z2)/(abs(z1)+abs(z2)),prob=TRUE)
dens1 <- density(abs(z1-z2)/(abs(z1)+abs(z2)))
dens2 <- density((abs(z1-z2)-1/n^a)/(abs(z1)+abs(z2)+1/n^a))
plot(dens1,col=1,ylim=range(c(dens1$y,dens2$y)),xlim=range(c(dens1$x,dens2$x)))
lines(dens2,col=2)


## z.diff
B <- 1e4
a <- b <- .3
var.beta <- a*b/(a+b)^2/(a+b+1)
sd.beta <- sqrt(var.beta)
rZ <- function(n)(rbeta(n,a,b)-1/2)/sd.beta
dZ <- function(x)dbeta(1/2+x*sd.beta,a,b)*sd.beta
pZ <- function(q)pbeta(1/2+q*sd.beta,a,b)
## hist(rZ(1e4),prob=TRUE)
## curve(dZ(x),add=TRUE,col=2)
## plot(ecdf(rZ(1e4)))
## curve(pZ(x),add=TRUE,col=2)
z <- rZ(B)
mean(dZ(z)) - mean(z*pZ(z))

B <- 4e3
abs <- matrix(runif(2*B)*10+1/2,ncol=2)
z.diffs <- apply(abs,1, function(ab) {
    a <- ab[1]; b <- ab[2]
    a <- b <- ab[1]
    var.beta <- a*b/(a+b)^2/(a+b+1)
    sd.beta <- sqrt(var.beta)
    rZ <- function(n)(rbeta(n,a,b)-1/2)/sd.beta
    dZ <- function(x)dbeta(1/2+x*sd.beta,a,b)*sd.beta
    pZ <- function(q)pbeta(1/2+q*sd.beta,a,b)
    ## hist(rZ(1e4),prob=TRUE)
    ## curve(dZ(x),add=TRUE,col=2)
    ## plot(ecdf(rZ(1e4)))
    ## curve(pZ(x),add=TRUE,col=2)
    ## mean(dZ(z)) - mean(z*pZ(z))
    obj.1 <- sd.beta*beta(2*a-1,2*b-1)/beta(a,b)^2
    z <- rZ(B)
    obj.2 <- mean(z*pZ(z))
    obj.1-obj.2
})
min(z.diffs)

ab.min <- abs[which.min(z.diffs),]
var.beta <- a*b/(a+b)^2/(a+b+1)
dZ <- function(x)dbeta(1/2+x*sd.beta,ab.min[2],ab.min[1])*sqrt(var.beta)
curve(dZ,-5,5)


as <- seq(.6,20,len=2000)
z.diffs <- sapply(as, function(a) {
    var.beta <- a^2/(2*a)^2/(2*a+1)
    sd.beta <- sqrt(var.beta)
obj.1 <- sd.beta*beta(2*a-1,2*a-1)/beta(a,a)^2
obj.2 <- integrate(function(x)pbeta(x,a,a)*(1-pbeta(x,a,a)),0,1)$val / sd.beta / 2
    obj.1-obj.2
})
plot(as,z.diffs,type='l'); abline(h=0)
abline(v=1,col=2)

dd


n <- 5
y <- rnorm(n)
sigma <- rnorm(n)
theta.fe <- sum(y/sigma^2)/sum(1/sigma^2)
coef(lm(I(y/sigma) ~ I(1/sigma)-1))
coef(lm(I((y-theta.fe)/sigma) ~ I(1/sigma)))
coef(lm(I(y/sigma) ~ I(1/sigma)))

x <- rnorm(n)
theta.fe <- sum(y*x^2)/sum(x^2)
coef(lm(I(y*x)~x))
proj <- coef(lm(I(y*x)~x-1))
coef(lm(I(y*x-proj*x)~x))



x <- rnorm(n)
## theta.fe <- sum(y*x^2)/sum(x^2)
coef(lm(y~x))
proj <- coef(lm(y~x-1))
coef(lm(I(y-proj*x)~x))


a <- runif(1)*5;b <- runif(1)*5
a <- 1; b <- 5
var.gamma <- a/b^2
integrate(function(x)x*pgamma(x,a,b)*dgamma(x,a,b),0,10)$val
1/2*integrate(function(x)1-pgamma(x,a,b)^2,0,10)$val
B <- 1e4
mean(abs(rgamma(B,a,b)-rgamma(B,a,b)))/2 
integrate(function(x)pgamma(x,a,b)*(1-pgamma(x,a,b)),0,100)$val

a <- 1
bs <- seq(1,5,length=10)
obs <- sapply(bs, function(b)mean(abs(rgamma(B,a,b)-rgamma(B,a,b)))/2 )
fla <- sapply(bs, function(b)integrate(function(x)pgamma(x,a,b)*(1-pgamma(x,a,b)),0,100)$val)
plot(bs,obs,type='l')
lines(bs,fla,col=2)

F0 <- function(x)pgamma(x,a,.1)
F1 <- function(x)pgamma(x,a,10)
lambdas <- seq(0,1,len=20)
## curve(1/2*(pgamma(x,a,1)+pgamma(x,a,5)),0,5,add=TRUE,col=2)
obj <- function(F)integrate(function(x)F(x)*(1-F(x)),0,100)$val
obs <- sapply(lambdas, function
fla <- sapply(lambdas, function(lambda)obj(function(x)lambda*F0(x)+(1-lambda)*F1(x)))
plot(lambdas,fla)


obj(function(x)1/2*(pgamma(x,a,1)+pgamma(x,a,5)))

a <- .15;b <- .39
s2.beta <- a*b/(a+b)^2/(a +b+1) + (a/(a+b))^2
2*integrate(function(x)pbeta(x,a,b)*(1-pbeta(x,a,b)),0,1)$val/sqrt(s2.beta)
mean(abs(rbeta(B,a,b)-rbeta(B,a,b)))/sqrt(mean(rbeta(B,a,b)^2))

while(1) {
abc <- runif(3)
stopifnot(sum(abc^2)-3*a*b>0)
}


2*b*acosh(-b/a)-1






start <- Sys.time()
source('../begg/2/misc.R')
require(parallel)

sim <- function(B,n,rZ,rS=runif) {
    ## B <- 1e2
    ## pb.alpha <- .05
    ## ma.alpha <- .05
    ## ns <- round(seq(1e1,3e2,len=10))
    ## ## by.n <- mclapply(ns, mc.cores=detectCores()-2, FUN=function(n) {
    ##     by.n <- lapply(ns, FUN=function(n) {
    ##     cat('.')
    replicate(B, {
        ## z <- with(power.Z(1), rZ(n))
        ## z <- rnorm(n)
        ## z <- rt(n,df=10)
        ## z <- with(beta.Z(c(.25,.25)),rZ(n))
        z <- rZ(n)
        s <- rS(n)
        y <- z/s
        v <- 1/s^2
        theta.fe <- sum(y/v)/sum(1/v)
        var.theta.fe <- 1/sum(1/v)
        ## c(egger=egger.test(y,v),begg=begg.test(y,v,exact=TRUE),ma.stat=theta.fe/sqrt(1/sum(1/v)))
        ## begg.stat <- tau.hat.pi(z,s,0) + theta.fe*D
        begg.stat <- tau.hat(z,s,theta.fe)
        begg.pval <- (1-pnorm(abs(begg.stat)*sqrt(9*n/4)))*2
        c(egger=egger.test(y,v),begg.stat=begg.stat,begg.pval=begg.pval,ma.stat=theta.fe/sqrt(var.theta.fe))        
    })
}


process <- function(stats,pb.alpha=.05,ma.alpha=.05) {    
    null.idx <- stats['begg.pval',] > pb.alpha
    begg.power <- mean(1-pnorm(abs(stats['ma.stat',null.idx])) < ma.alpha/2)
    null.idx <- stats['egger.pval',] > pb.alpha
    egger.power <- mean(1-pnorm(abs(stats['ma.stat',null.idx])) < ma.alpha/2)
    unconditional.power <- mean(1-pnorm(abs(stats['ma.stat',])) < ma.alpha/2)
    c(begg=begg.power,egger=egger.power,unconditional=unconditional.power) 
}


dfs <- seq(2.5,7,len=10)
rZs <- lapply(dfs,function(df)
    function(n)rt(n,df=df)/sqrt(df/(df-2)) )
n <- 20
B <- 1e3
stats.out <- lapply(rZs, function(rZ)sim(B=B,n=n,rZ,rZ))
process.out <- lapply(stats.out, function(stats) process(stats))
process.out <- simplify2array(process.out)
plot(dfs,process.out['unconditional',],ylim=range(process.out),type='l')
lines(dfs,process.out['egger',],col=2)
lines(dfs,process.out['begg',],col=3)
legend('bottomleft',lty=1,col=1:3,legend=c('unconditional','egger','begg'))




## 21 zdiff optimization

obj <- function(par) {
    a <- par[1]; b <- par[2]
    if(a>0 | b< -a)return(Inf)
    u <- (1+2*(a+b)/sqrt(1+2*a/(b-a)))/(2*b)
    eqn1 <- acosh(-b/a) - u
    eqn2 <- 4*b*u+2/3*b*u^3+2*(a-b)*sqrt((a+b)/(b-a))*(2+u^2)
    eqn1^2+eqn2^2
    }

optim(c(-1,1.0011),obj)

diff(abs(optim(c(-103,104),obj)$par))


z.diff <- function(par){
    a <- par[1]; b <- par[2]; d <- par[3]
    -a^2*d+2*b^2*d-b*d^2+(4*b^2*d^3)/3+4*a*b*d*cosh(d/sqrt(2))+a^2*d*cosh(sqrt(2)*d)-sqrt(2)*a*d*sinh(d/sqrt(2))+2*sqrt(2)*a*b*d^2*sinh(d/sqrt(2))+sqrt(2)*a^2*sinh(sqrt(2)*d)
}
moment0 <- function(par) {
    a <- par[1]; b <- par[2]; d <- par[3]
    2*b*d+a*2*sqrt(2)*sinh(d/sqrt(2))
    }
moment2 <- function(par) {
    a <- par[1]; b <- par[2]; d <- par[3]
    2*b*d^3/3-8*a*d*cosh(d/sqrt(2))+2*sqrt(2)*a*(4+d^2)*sinh(d/sqrt(2))
    }
obj <- function(par) {
    a <- par[1]; b <- par[2]; d <- par[3]
    if(a<=0  && (-b/a<=1 || d<=0 || d>sqrt(2)*acosh(-b/a)))return(Inf)
    if(a>=0 && (a+b<=0 || d<=0))return(Inf)
    ## moment0 <- 2*b*d+a*2*sqrt(2)*sinh(d/sqrt(2))
    ## moment2 <- 2*b*d^2/3-8*a*d*cosh(d/sqrt(2))+2*sqrt(2)*a*(4+d^2)*sinh(d/sqrt(2))
    (moment0(par)-1)^2 + (moment2(par)-1)^2
}
## optim0 <- optim(c(-1,2,1.2),obj)
optim0 <- optim(c(.2,.04,.2),obj)
par0
par0 <- optim0$par
c(moment0(par0),moment2(par0))
z.diff(par0)

init.pars <- matrix(runif(1e4*3),ncol=3)
z.diffs <- apply(init.pars,1, function(x)z.diff(optim(x,obj,control=list(maxit=5e4,abstol=1e-12,reltol=1e-12))$par))
par0 <- optim(init.pars[which.min(z.diffs),],obj)$par
par0
min(z.diffs)


a <- par0[1]; b <- par0[2]; d <- par0[3]
pdf <- function(x)(a*cosh(x/sqrt(2))+b)*(abs(x)<=d)
cdf <- function(x)b*d+b*x+sqrt(2)*a*(sinh(d/sqrt(2))+sinh(x/sqrt(2)))
## curve(a*cosh(x/sqrt(2))+b,-d,d)
curve(pdf,-d,d,ylim=c(0,1))
curve(cdf,add=TRUE,col=2)

integrate(function(x)pdf(x)^2,-d,d)
integrate(function(x)pdf(x),-d,d)
integrate(function(x)pdf(x)*x^2,-d,d)




## repeating constraining endpoints to be intersection with x-axis, since that is where min appears to be anyway 


z.diff <- function(par){
    a <- par[1]; b <- par[2]; d <- sqrt(2)*acosh(-b/a)
    -a^2*d+2*b^2*d-b*d^2+(4*b^2*d^3)/3+4*a*b*d*cosh(d/sqrt(2))+a^2*d*cosh(sqrt(2)*d)-sqrt(2)*a*d*sinh(d/sqrt(2))+2*sqrt(2)*a*b*d^2*sinh(d/sqrt(2))+sqrt(2)*a^2*sinh(sqrt(2)*d)
}
moment0 <- function(par) {
    a <- par[1]; b <- par[2]; d <- sqrt(2)*acosh(-b/a)
    2*b*d+a*2*sqrt(2)*sinh(d/sqrt(2))
    }
moment2 <- function(par) {
    a <- par[1]; b <- par[2]; d <- sqrt(2)*acosh(-b/a)
    2*b*d^3/3-8*a*d*cosh(d/sqrt(2))+2*sqrt(2)*a*(4+d^2)*sinh(d/sqrt(2))
    }
obj <- function(par) {
    a <- par[1]; b <- par[2]; d <- sqrt(2)*acosh(-b/a)
    if(a<=0  && (-b/a<=1 || d<=0 ))return(Inf)
    ## if(a>=0 && (a+b<=0 || d<=0))return(Inf)
    ## moment0 <- 2*b*d+a*2*sqrt(2)*sinh(d/sqrt(2))
    ## moment2 <- 2*b*d^2/3-8*a*d*cosh(d/sqrt(2))+2*sqrt(2)*a*(4+d^2)*sinh(d/sqrt(2))
    (moment0(par)-1)^2 + (moment2(par)-1)^2
}
optim0 <- optim(c(-1,2,1.2),obj,control=list(maxit=5e4,abstol=1e-16,reltol=1e-16))
## optim0 <- optim(c(.2,.04,.2),obj)
par0
par0 <- optim0$par
c(moment0(par0),moment2(par0))
z.diff(par0)



dd

## egger orthogonality figure

s <- runif(2)
v <- 1/s^2
y <- rnorm(2,sd=sqrt(v))
theta.fe <- sum(y/v)/sum(1/v)
plot(s[1],s[2],xlim=range(c(s,1)),ylim=range(c(s,1)))
points(1,1)


## zdiff for student's t

B <- 1e5
df <- 2.01
dfs <- seq(2.000001,2.005,len=60)
pairs <- sapply(dfs, function(df) {
sd <- sqrt(df/(df-2))
dZ <- function(z)dt(z*sd,df=df)*sd
rZ <- function(n)rt(n,df=df)/sd
pZ <- function(q)pt(q*sd,df=df)
E.f <- mean(dZ(rZ(B)))
z <- rZ(B)
E.ZF <- mean(z*pZ(z))
c(E.f,E.ZF)
})
plot(dfs,pairs[1,]-pairs[2,])





