## 1. Q can be written as a quad form in the blocks, what does this
## mean?  2. discuss m=n=1 case of phi form. sum is auc, etc. obj
## worth studying on its own, to motivate description of general
## phi. 3. check harmonic functions, others with extrema at boundary

## 1. comparing jk and obu variance estimates
source('misc.R')
require(mvtnorm)
I <- 1e1
m <- 10
n <- 10
rho <- .3
Sigma.XX <- matrix(0,nrow=m,ncol=m)
Sigma.XX <- rho^(abs(row(Sigma.XX)-col(Sigma.XX)))
Sigma.YY <- matrix(0,nrow=n,ncol=n)
Sigma.YY <- rho^(abs(row(Sigma.YY)-col(Sigma.YY)))
Sigma.XY <- matrix(0,nrow=m,ncol=n)
Sigma.XY <- rho^(abs(row(Sigma.XY)-col(Sigma.XY)))/2
Sigma <- rbind(cbind(Sigma.XX,Sigma.XY),cbind(t(Sigma.XY),Sigma.YY))
mu.X <- 0; mu.Y <- .5
auc.hats <- replicate(5e2,{
    xy <- rmvnorm(I,c(rep(mu.X,m),rep(mu.Y,n)),sigma=Sigma)
    x <- t(xy[,1:m]); y <- t(xy[,(m+1):(m+n)])
    x <- unname(as.list(as.data.frame(x)))
    y <- unname(as.list(as.data.frame(y)))
    rbind(obu=auc.obu(x,y), jk=auc.jk(x,y))
})
diff(c(var(auc.hats['obu','theta.hat',]),mean(auc.hats['obu','var.hat',])))
diff(c(var(auc.hats['jk','theta.hat',]),mean(auc.hats['jk','var.hat',])))

summary(abs(auc.hats['obu','var.hat',] - auc.hats['jk','var.hat',]))

## varying I, seems max abs diff is O(1/I^2.5). not very useful
## measure, since both estimators are consistent for the same target,
## use relative error instead.
require(parallel)
require(mvtnorm)
Is <- round(seq(8,1e2,len=30))
by.I <- mclapply(Is,mc.cores=detectCores()-3,FUN=function(I){
    print(I)
    auc.hats <- replicate(1e2,{
        xy <- rmvnorm(I,c(rep(mu.X,m),rep(mu.Y,n)),sigma=Sigma)
        x <- t(xy[,1:m]); y <- t(xy[,(m+1):(m+n)])
        x <- unname(as.list(as.data.frame(x)))
        y <- unname(as.list(as.data.frame(y)))
        rbind(obu=auc.obu(x,y), jk=auc.jk(x,y))
    })
    c(max=max(abs(auc.hats['obu','var.hat',] - auc.hats['jk','var.hat',])),
      mad=mean(abs(auc.hats['obu','var.hat',] - auc.hats['jk','var.hat',])),
      mean=mean(auc.hats['obu','var.hat',] - auc.hats['jk','var.hat',])
      )
})
by.I <- simplify2array(by.I)
## save.image('210809.RData')
plot(Is,by.I['mad',],type='o')

plot(Is,by.I['max',]*Is^2.5,type='o')
plot(Is,by.I['mad',]*Is^2.3,type='o')
plot(Is,by.I['mean',]*Is^2.5,type='o')



## relative difference, appears O(I^(3/2))
source('misc.R')
require(parallel)
require(mvtnorm)
m <- 5
n <- 5
rho <- .3
Sigma.XX <- matrix(0,nrow=m,ncol=m)
Sigma.XX <- rho^(abs(row(Sigma.XX)-col(Sigma.XX)))
Sigma.YY <- matrix(0,nrow=n,ncol=n)
Sigma.YY <- rho^(abs(row(Sigma.YY)-col(Sigma.YY)))
Sigma.XY <- matrix(0,nrow=m,ncol=n)
Sigma.XY <- rho^(abs(row(Sigma.XY)-col(Sigma.XY)))/2
Sigma <- rbind(cbind(Sigma.XX,Sigma.XY),cbind(t(Sigma.XY),Sigma.YY))
mu.X <- 0; mu.Y <- .5
Is <- round(seq(8,100,len=100))
by.I <- mclapply(Is,mc.cores=detectCores()-3,FUN=function(I){
    print(I)
    auc.hats <- replicate(1e2,{
        xy <- rmvnorm(I,c(rep(mu.X,m),rep(mu.Y,n)),sigma=Sigma)
        x <- t(xy[,1:m]); y <- t(xy[,(m+1):(m+n)])
        x <- unname(as.list(as.data.frame(x)))
        y <- unname(as.list(as.data.frame(y)))
        rbind(obu=auc.obu(x,y), jk=auc.jk(x,y))
    })
    max(abs(auc.hats['obu','var.hat',] - auc.hats['jk','var.hat',])/auc.hats['obu','var.hat',])
})
by.I <- simplify2array(by.I)
plot(Is,by.I,type='o')
plot(Is,by.I*Is^(1.5),type='o')

## relative difference, variable m,n. same outcome
source('misc.R')
require(parallel)
mu.X <- 0; mu.Y <- .5
Is <- round(seq(8,50,len=50))
by.I <- mclapply(Is,mc.cores=detectCores()-3,FUN=function(I){
    print(I)
    auc.hats <- replicate(1e2,{
        m <- sample(1:10,I,replace=TRUE); n <- sample(1:10,I,replace=TRUE)
        x <- lapply(m,function(m.i)rnorm(m.i))
        y <- lapply(n,function(n.i)rnorm(n.i))
        rbind(obu=auc.obu(x,y), jk=auc.jk(x,y))
    })
    max(abs(auc.hats['obu','var.hat',] - auc.hats['jk','var.hat',])/auc.hats['obu','var.hat',])
})
by.I <- simplify2array(by.I)
plot(Is,by.I,type='o')
plot(Is,by.I*Is^(1.5),type='o')

dd
## under delong setup
M <- 7
N <- 10
x <- c(rep(list(list()),N),as.list(rnorm(M)))
y <- c(as.list(rnorm(N)),rep(list(list()),M))
auc.obu(x,y)
auc.jk(x,y)
auc.jk2(x,y)


## |mean-median| for comparison is O(1/sqrt(I). not very relevant
## |since in our case the target is going to 0
Is <- round(seq(1,500,len=20))
pairs <- lapply(Is, function(I) {
    replicate(1e2, {
    x <- rnorm(I)
    c(mean=mean(x),median=median(x))
    })
})
pairs <- simplify2array(pairs)
diffs <- abs(pairs['mean',,] - pairs['median',,])
plot(Is,colMeans(diffs)*sqrt(Is))


## exact variance formulas for m=n=1 case
source('misc.R')
require(mvtnorm)
I <- 5e1
rho <- .3
Sigma <- matrix(c(1,rho,rho,1),2)
mu.X <- 0.5; mu.Y <- .7
## auc.hats <- replicate(5e2,{
## xy <- rmvnorm(2,c(mu.X,mu.Y),sigma=Sigma)
## x <- xy[,1]; y <- xy[,2]
## ## x <- t(xy[,1:m]); y <- t(xy[,(m+1):(m+n)])
## mean(sapply(1:I, function(i) mean(apply(combn((1:I)[-i],2),2,function(idx)(x[i]<y[idx[1]])*(x[1]<y[idx[2]])))))
## B <- 1e4
## mean(rnorm(B,mu.X)<pmin(rnorm(B,mu.Y),rnorm(B,mu.Y)))
## mean((1-pnorm(rnorm(1e4,mu.X),mean=mu.Y))^2)
ijik <- E.10 <- integrate(function(x)(1-pnorm(x,mean=mu.Y))^2*dnorm(x,mu.X),mu.X-7,mu.X+7)$val
## mean(rnorm(B,mu.Y)>pmax(rnorm(B,mu.X),rnorm(B,mu.X)))
## mean(pnorm(rnorm(1e4,mu.Y),mean=mu.X)^2)
ijkj <- E.01 <- integrate(function(x)pnorm(x,mean=mu.X)^2*dnorm(x,mu.Y),mu.Y-7,mu.Y+7)$val
## xy <- rmvnorm(B,c(mu.X,mu.Y),sigma=Sigma)
## x <- xy[,1]; y <- xy[,2]
## mean((rnorm(B,mu.X)<y)*(x<rnorm(B,mu.Y)))
## x <- .4
ijki <- E.11 <- integrate(  Vectorize(function(x)  integrate(Vectorize(function(y)pnorm(y,mu.X)*(1-pnorm(x,mu.Y))*dmvnorm(c(x,y),c(mu.X,mu.Y),Sigma)),mu.Y-7,mu.Y+7)$val), mu.X-7,mu.X+7)$val
theta <- integrate(function(y)pnorm(y,mu.X)*dnorm(y,mu.Y),mu.Y-7,mu.Y+7)$val
pop.var <- E.10 - theta^2 + E.01 - theta^2 + 2*(E.11-theta^2)
pop.var/I

xy <- rmvnorm(I,c(mu.X,mu.Y),sigma=Sigma)
x <- t(xy[,1]); y <- t(xy[,2])
x <- unname(as.list(as.data.frame(x)))
y <- unname(as.list(as.data.frame(y)))
auc.obu(x,y)['var.hat']
## theta.hat <- concord.xy(unlist(x),unlist(y)) / I^2

## these variance estimators are sqrt(I) rate, but of the population variance, so after normalizing by I
source('misc.R')
require(mvtnorm)
rho <- .3
Sigma <- matrix(c(1,rho,rho,1),2)
mu.X <- 0.5; mu.Y <- .7
E.10 <- integrate(function(x)(1-pnorm(x,mean=mu.Y))^2*dnorm(x,mu.X),mu.X-7,mu.X+7)$val
E.01 <- integrate(function(x)pnorm(x,mean=mu.X)^2*dnorm(x,mu.Y),mu.Y-7,mu.Y+7)$val
E.11 <- integrate(  Vectorize(function(x)  integrate(Vectorize(function(y)pnorm(y,mu.X)*(1-pnorm(x,mu.Y))*dmvnorm(c(x,y),c(mu.X,mu.Y),Sigma)),mu.Y-7,mu.Y+7)$val), mu.X-7,mu.X+7)$val
theta <- integrate(function(y)pnorm(y,mu.X)*dnorm(y,mu.Y),mu.Y-7,mu.Y+7)$val
pop.var <- E.10 - theta^2 + E.01 - theta^2 + 2*(E.11-theta^2)
Is <- round(seq(5,20,len=20))
var.hats <- lapply(Is, function(I) {
    print(I)
    var.hats <- replicate(1e2,{
        xy <- rmvnorm(I,c(mu.X,mu.Y),sigma=Sigma)
        x <- t(xy[,1]); y <- t(xy[,2])
        x <- unname(as.list(as.data.frame(x)))
        y <- unname(as.list(as.data.frame(y)))
        c(jk=unname(auc.jk(x,y)['var.hat']),jk2=unname(auc.jk2(x,y)['var.hat']))
    })
})
var.hats <- simplify2array(var.hats)

matplot(Is,Is*t(var.hats['jk',,]),col=1,pch=1)
abline(h=pop.var,col=2)

matplot(Is,(Is*t(var.hats['jk',,])-pop.var)*sqrt(Is),col=1,pch=1)

matplot(Is,Is*t(var.hats['jk2',,]),col=1,pch=1)
abline(h=pop.var,col=2)
matplot(Is,(Is*t(var.hats['jk2',,])-pop.var)*sqrt(Is),col=1,pch=1)

hist(var.hats)
abline(v=pop.var/I,col='red')
abline(v=mean(var.hats),col='blue')

## m=n=1 case as before, but now also no intracluster correlation (for
## fair comparison with 2-sample jk)
m <- 10
n <- 10
mu.x <- 1; sd.x <- 2
mu.y <- -1; sd.y <- 1
auc.hats <- replicate(1e3, {
x <- rnorm(m,mu.x,sd.x)
y <- rnorm(n,mu.y,sd.y)
auc.hat <- mean(outer(x,y,`<`))
})
var(auc.hats)*(m+n)
E.yxx <- integrate(function(x)pnorm(x,mu.y,sd.y)^2*dnorm(x,mu.x,sd.x),mu.x-7,mu.x+7)$val
E.xyy <- integrate(function(y)pnorm(y,mu.x,sd.x)^2*dnorm(y,mu.y,sd.y),mu.y-7,mu.y+7)$val
theta <- integrate(function(x)pnorm(x,mu.y,sd.y)*dnorm(x,mu.x,sd.x),mu.x-7,mu.x+7)$val
(m+n)*((E.yxx-theta^2)/m + (E.xyy-(1-theta)^2)/n)

dd

## 2. two-sample jk


auc.jk2 <- function(x,y,alpha=.05) {
    concord <- function(x,y) mean(outer(x,y,`<`))
    concord.sum <- function(x,y)sum(outer(x,y,`<`))
    m <- sapply(x,length); n <- sapply(y,length)
    M <- sum(m); N <- sum(n)
    I <- length(x); stopifnot(length(x)==length(y))

    theta.hat <- concord(unlist(x),unlist(y))
    ## theta.del <- sapply(1:I,function(i)concord(unlist(x[-i]),unlist(y[-i])))
    theta.del.x <- sapply(1:I,function(i)concord(unlist(x[-i]),unlist(y)))
    theta.del.y <- sapply(1:I,function(j)concord(unlist(x),unlist(y[-j])))
    pseudo.x <- I*theta.hat - (I-1)*theta.del.x
    pseudo.y <- I*theta.hat - (I-1)*theta.del.y
    ## pseudo <- I*theta.hat - (I-1)*theta.del
    theta.jk <- mean(c(pseudo.x,pseudo.y))
    var.hat <- (var(pseudo.x)+var(pseudo.y))/I
    q <- qnorm(1-alpha/2)

    return(c(theta.hat=theta.jk,var.hat=var.hat,CI.lower=theta.jk-q*sqrt(var.hat),CI.upper=theta.jk+q*sqrt(var.hat)))
}

require(mvtnorm)
I <- 1e1
m <- 5
n <- 5
rho <- .3
Sigma.XX <- matrix(0,nrow=m,ncol=m)
Sigma.XX <- rho^(abs(row(Sigma.XX)-col(Sigma.XX)))
Sigma.YY <- matrix(0,nrow=n,ncol=n)
Sigma.YY <- rho^(abs(row(Sigma.YY)-col(Sigma.YY)))
Sigma.XY <- matrix(0,nrow=m,ncol=n)
Sigma.XY <- rho^(abs(row(Sigma.XY)-col(Sigma.XY)))/2
Sigma <- rbind(cbind(Sigma.XX,Sigma.XY),cbind(t(Sigma.XY),Sigma.YY))
mu.X <- 0; mu.Y <- .5

auc.hats <- replicate(5e2,{
    xy <- rmvnorm(I,c(rep(mu.X,m),rep(mu.Y,n)),sigma=Sigma)
    x <- t(xy[,1:m]); y <- t(xy[,(m+1):(m+n)])
    x <- unname(as.list(as.data.frame(x)))
    y <- unname(as.list(as.data.frame(y)))
    rbind(obu=auc.obu(x,y), jk=auc.jk(x,y), jk2=auc.jk2(x,y))
})

diff(c(var(auc.hats['obu','theta.hat',]),mean(auc.hats['obu','var.hat',])))
diff(c(var(auc.hats['jk','theta.hat',]),mean(auc.hats['jk','var.hat',])))

summary(abs(auc.hats['obu','var.hat',] - auc.hats['jk','var.hat',]))
summary(abs(auc.hats['obu','var.hat',] - auc.hats['jk2','var.hat',]))
summary(abs(auc.hats['jk','var.hat',] - auc.hats['jk2','var.hat',]))
## not as close to obu estimator as one-sample jk

## 2a varying I, relative difference
require(mvtnorm)
require(parallel)
I <- 1e1
m <- 5
n <- 5
rho <- .3
Sigma.XX <- matrix(0,nrow=m,ncol=m)
Sigma.XX <- rho^(abs(row(Sigma.XX)-col(Sigma.XX)))
Sigma.YY <- matrix(0,nrow=n,ncol=n)
Sigma.YY <- rho^(abs(row(Sigma.YY)-col(Sigma.YY)))
Sigma.XY <- matrix(0,nrow=m,ncol=n)
Sigma.XY <- rho^(abs(row(Sigma.XY)-col(Sigma.XY)))/2
Sigma <- rbind(cbind(Sigma.XX,Sigma.XY),cbind(t(Sigma.XY),Sigma.YY))
mu.X <- 0; mu.Y <- .5

Is <- round(seq(8,500,len=50))
by.I <- mclapply(Is,mc.cores=detectCores()-3,FUN=function(I){
## by.I <- lapply(Is,FUN=function(I){
    print(I)
    auc.hats <- replicate(50,{
        xy <- rmvnorm(I,c(rep(mu.X,m),rep(mu.Y,n)),sigma=Sigma)
        x <- t(xy[,1:m]); y <- t(xy[,(m+1):(m+n)])
        x <- unname(as.list(as.data.frame(x)))
        y <- unname(as.list(as.data.frame(y)))
        rbind(obu=auc.obu(x,y), jk=auc.jk(x,y), jk2=auc.jk2(x,y))
    })
    c(jk.diff=max(abs(auc.hats['obu','var.hat',] - auc.hats['jk','var.hat',])),jk2.diff=max(abs(auc.hats['obu','var.hat',] - auc.hats['jk2','var.hat',])),jk.reldiff=max(abs(auc.hats['obu','var.hat',] - auc.hats['jk','var.hat',])/auc.hats['obu','var.hat',]),jk2.reldiff=max(abs(auc.hats['obu','var.hat',] - auc.hats['jk2','var.hat',])/auc.hats['obu','var.hat',]))
})
by.I <- simplify2array(by.I)

reldiff.idx <- 1:nrow(by.I) %in% grep('reldiff',rownames(by.I))
op <- par(mfrow=c(1,2))
plot(Is,by.I['jk.diff',],ylim=range(by.I[!reldiff.idx,]),type='o')
lines(Is,by.I['jk2.diff',],type='o',lty=2)
plot(Is,by.I['jk.reldiff',],ylim=range(by.I[reldiff.idx,]),type='o')
lines(Is,by.I['jk2.reldiff',],type='o',lty=2)
par(op)
## similar outcome, single sample jk much closer to obu estimator, relative diff ot two sample jk to obu seems to stabilize around 2 (but see below)

plot(Is,by.I['jk2.reldiff',]*Is^(.2),type='o')

## ## 2b repeat 2a with uncorrelated data (faster)
## require(mvtnorm)
## require(parallel)
## I <- 1e1
## m <- 5
## n <- 5
## Is <- round(seq(8,100,len=50))
## by.I <- mclapply(Is,mc.cores=detectCores()-3,FUN=function(I){
## ## by.I <- lapply(Is,FUN=function(I){
##     print(I)
##     auc.hats <- replicate(1e2,{
##         x <- lapply(1:I,function(m)rnorm(m))
##         y <- lapply(1:I,function(n)rnorm(n))
##         rbind(obu=auc.obu(x,y), jk=auc.jk(x,y), jk2=auc.jk2(x,y))
##     })
##     c(jk.diff=max(abs(auc.hats['obu','var.hat',] - auc.hats['jk','var.hat',])),jk2.diff=max(abs(auc.hats['obu','var.hat',] - auc.hats['jk2','var.hat',])),jk.reldiff=max(abs(auc.hats['obu','var.hat',] - auc.hats['jk','var.hat',])/auc.hats['obu','var.hat',]),jk2.reldiff=max(abs(auc.hats['obu','var.hat',] - auc.hats['jk2','var.hat',])/auc.hats['obu','var.hat',]))
## })
## by.I <- simplify2array(by.I)

## reldiff.idx <- 1:nrow(by.I) %in% grep('reldiff',rownames(by.I))
## op <- par(mfrow=c(1,2))
## plot(Is,by.I['jk.diff',],ylim=range(by.I[!reldiff.idx,]),type='o')
## lines(Is,by.I['jk2.diff',],type='o',lty=2)
## plot(Is,by.I['jk.reldiff',],ylim=range(by.I[reldiff.idx,]),type='o')
## lines(Is,by.I['jk2.reldiff',],type='o',lty=2)
## par(op)

## plot(Is,by.I['jk2.reldiff',]*Is^(.2),type='o')




## 3. checking formulas
source('misc.R')
concord.mean <- function(x,y)mean(outer(x,y,`<`))
concord.sum <- function(x,y)sum(outer(x,y,`<`))
m <- sapply(x,length); n <- sapply(y,length)
M <- sum(m); N <- sum(n)
I <- length(x)
theta.hat <- concord(unlist(x),unlist(y))
## theta.del <- sapply(1:I,function(i)concord(unlist(x[-i]),unlist(y[-i])))
theta.del.x <- sapply(1:I,function(i)concord(unlist(x[-i]),unlist(y)))
theta.del.y <- sapply(1:I,function(j)concord(unlist(x),unlist(y[-j])))
pseudo.x <- I*theta.hat - (I-1)*theta.del.x
pseudo.y <- I*theta.hat - (I-1)*theta.del.y

I <- 5
m <- 10; n <- 5
M <- I*m; N <- I*n
x <- replicate(I,rnorm(m),simplify=FALSE)
y <- replicate(I,rnorm(n),simplify=FALSE)
phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
psi <- phi*m*n
theta.hat <- concord.mean(unlist(x),unlist(y))
theta.del <- sapply(1:I,function(i)concord.mean(unlist(x[-i]),unlist(y[-i])))
V.x <- sapply(x,function(x.i)concord.sum(x.i,unlist(y)))/N
V.y <- sapply(y,function(y.i)concord.sum(unlist(x),y.i))/M
V.bar.x <- V.x/m; V.bar.y <- V.y/n
## theta.del - mean(theta.del)
## (I-1)/I*theta.del - sapply(1:I,function(i)sum(theta.del[-i]))/I
theta.del
1/((I-1)^2*m*n) * (sum(psi) - colSums(psi) - rowSums(psi) + diag(psi))
(I/(I-1))^2*theta.hat - I/(I-1)^2*(V.y/n + V.x/m) + diag(phi)/(I-1)^2
mean(theta.del)
(I/(I-1))^2*theta.hat - 2*I/(I-1)^2*theta.hat + sum(diag(phi))/I/(I-1)^2
theta.del - mean(theta.del)
-I/(I-1)^2*(V.bar.x + V.bar.y - 2*theta.hat) + 1/(I-1)^2*(diag(phi)-mean(diag(phi)))
var(theta.del)*(I-1)^2/I
auc.jk(x,y)['var.hat']
I/(I-1)^3*sum( (V.bar.x + V.bar.y - 2*theta.hat - 1/I*(diag(phi)-mean(diag(phi))))^2 )

auc.obu(x,y)['var.hat']
1/I/(I-1) * sum((V.bar.x + V.bar.y - 2*theta.hat)^2)

a <- (2*I-1)/I/(I-1)^3 * sum( (V.bar.x + V.bar.y - 2*theta.hat)^2 )
b <- 2/(I-1)^3 * sum( (V.bar.x+V.bar.y-2*theta.hat)*(diag(phi)-mean(diag(phi))) )
c <- 1/I/(I-1)^3 * sum( (diag(phi)-mean(diag(phi)))^2 )
auc.jk(x,y)['var.hat'] - auc.obu(x,y)['var.hat']
a - b + c


## 3a relationships among parts of difference
source('misc.R')
I <- 10
m <- n <- 5
M <- I*m; N <- I*n
parts <- replicate(1e2,{
    x <- replicate(I,runif(m),simplify=FALSE)
    y <- replicate(I,rnorm(n),simplify=FALSE)
    phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
    ## phi <- matrix(0,I,I) #1
    ## phi <- ((row(phi)>=col(phi))+0)[,I:1]
    ## psi <- phi*m*n
    ## theta.hat <- concord.mean(unlist(x),unlist(y))
    ## ## theta.del <- sapply(1:I,function(i)concord.mean(unlist(x[-i]),unlist(y[-i])))
    ## V.x <- sapply(x,function(x.i)concord.sum(x.i,unlist(y)))/N
    ## V.y <- sapply(y,function(y.i)concord.sum(unlist(x),y.i))/M
    ## V.bar.x <- V.x/m; V.bar.y <- V.y/n
    V.bar.x <- rowMeans(phi); V.bar.y <- colMeans(phi)
    theta.hat <- mean(phi)
    ## theta.del - mean(theta.del)
    ## (I-1)/I*theta.del - sapply(1:I,function(i)sum(theta.del[-i]))/I
    a <- (2*I-1)/I/(I-1)^3 * sum( (V.bar.x + V.bar.y - 2*theta.hat)^2 )
    a2 <- (2*I)/I/(I-1)^3 * sum( (V.bar.x + V.bar.y - 2*theta.hat)^2 )
    b <- 2/(I-1)^3 * sum( (V.bar.x+V.bar.y-2*theta.hat)*(diag(phi)-mean(diag(phi))) )
    c <- 1/I/(I-1)^3 * sum( (diag(phi)-mean(diag(phi)))^2 )
    c2 <- 1/I/(I-1)^3 * sum( (diag(phi)-mean(diag(phi)))^2 ) - 1/I/(I-1)^3 * sum( (V.bar.x + V.bar.y - 2*theta.hat)^2 )
    a.tilde <- sum( (V.bar.x + V.bar.y - 2*theta.hat)^2 )
    b.tilde <- sum((V.bar.x+V.bar.y)*(diag(phi)-mean(diag(phi))))
    c(a=a,a2=a2,b=b,c=c,c2=c2,a.tilde=a.tilde,b.tilde=b.tilde)
})
summary(parts['a',] - parts['b',]+parts['c',])

plot(parts['a.tilde',],parts['b.tilde',]); abline(0,1)

plot(parts['a',],type='l',ylim=range(parts[c('a','b','c'),]))
lines(parts['b',],col=2)
lines(parts['c',],col=3)
lines(parts['b',]-parts['a',],col=4)
lines(-(parts['a',]-parts['b',]+parts['c',]),col=5,lty=2)

plot(parts['a',]-parts['b',],-parts['c',])

plot(parts['a',]-parts['b',],col=4,type='l')
lines(-parts['c',],col=3)

cor(parts['a',]-parts['b',],parts['c',])
## part c anti-correlated with part a - part b, but much smaller
## magnitude as I grows

plot(parts['a',],type='l')
lines(parts['b',],col=2)
lines(parts['a2',],col=3)

plot(parts['a',]-parts['b',],type='l')
lines(parts['a',]-parts['a2',],col=2)

## OK to ignore O(1/I) terms? seems so.
plot(parts['a2',]-parts['b',],type='l')
lines(parts['c2',],col=2)


## 3b Adversarial phi--lower right triangular. Quad form is O(1/I^2). I had previously suspected the quad form maxed at O(1/I^(2.5)).
source('misc.R')
Is <- floor(seq(3,100,len=20))
by.I <- sapply(Is, function(I) {
    phi <- matrix(0,I,I) #1
    phi <- ((row(phi)>=col(phi))+0)[,I:1]
    V.bar.x <- rowMeans(phi); V.bar.y <- colMeans(phi)
    theta.hat <- mean(phi)
    a <- (2*I-1)/I/(I-1)^3 * sum( (V.bar.x + V.bar.y - 2*theta.hat)^2 )
    b <- 2/(I-1)^3 * sum( (V.bar.x+V.bar.y-2*theta.hat)*(diag(phi)-mean(diag(phi))) )
    c <- 1/I/(I-1)^3 * sum( (diag(phi)-mean(diag(phi)))^2 )
    c(obj=a-b+c)*(I-1)^3
})
plot(Is,by.I/Is)

phi <- matrix(0,I,I)
phi <- ((row(phi)>=col(phi))+0)[,I:1]
V.bar.x <- rowMeans(phi); V.bar.y <- colMeans(phi)
theta.hat <- mean(phi)
var(V.bar.x+V.bar.y)/I
I/(I-1)^2*sum( (V.bar.x+V.bar.y-2*theta.hat - diag(phi)/I+mean(diag(phi))/I)^2)
var.obu <- (I+1)/3/I^2
var.jk <- 1/(I-1)^2*(I^2/3-I/2-1/12)

## I/(I-1)^2*sum( (V.bar.x+V.bar.y-2*theta.hat - diag(phi)/I+mean(diag(phi))/I)^2)
## ks <- 1:I
## sum( (2*ks/I-(I+1)/I-(ks>I/2)/I+1/2/I)^2)
## sum( (2*ks-I-1/2-(ks>I/2) )^2)
## sum(c(sum((2*ks)^2),I*(-I-1/2)^2,sum((ks>I/2)^2), 2*sum(2*ks*(-I-1/2)), -4*sum(ks*(ks>I/2)), 2*sum((-I-1/2)*(-(ks>I/2)))))
## sum(c(2/3*I*(I+1)*(2*I+1), I*(I+1/2)^2,I/2,-(2*I+1)*I*(I+1),-2*I*(3/4*I+1/2),I*(I+1/2)))
## -1/3*I*(I+1)*(2*I+1)+I*(I+1/2)^2+I*(I+1)-2*I*(3/4*I+1/2)
## 1/I/(I-1)^2* (-1/3*I*(I+1)*(2*I+1)+I*(I+1/2)^2+I/2*(I+3/2)-2*I*(3/4*I+1/2))
## (I^2-3*I-1)/3/(I-1)^2
## 1/12/(I-1)^2*(4*I^2-6*I-1)

xy <- sort(runif(2*I))
y <- xy[(1:I)*2]; x <- rev(xy[2*(1:I)-1])
phi <- outer(x,y,'<')+0
auc.jk(x,y)

concord.mean <- function(x,y)mean(outer(x,y,`<`))
concord.sum <- function(x,y)sum(outer(x,y,`<`))
I <- 5
m <- n <- 3
M <- I*m; N <- I*n
parts <- replicate(1e3,{
    x <- replicate(I,runif(m,0,10),simplify=FALSE)
    y <- replicate(I,rcauchy(n),simplify=FALSE)
    phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
    psi <- phi*m*n
    theta.hat <- concord.mean(unlist(x),unlist(y))
    ## theta.del <- sapply(1:I,function(i)concord.mean(unlist(x[-i]),unlist(y[-i])))
    V.x <- sapply(x,function(x.i)concord.sum(x.i,unlist(y)))/N
    V.y <- sapply(y,function(y.i)concord.sum(unlist(x),y.i))/M
    V.bar.x <- V.x/m; V.bar.y <- V.y/n
    ## theta.del - mean(theta.del)
    ## (I-1)/I*theta.del - sapply(1:I,function(i)sum(theta.del[-i]))/I
    ## c(
    ##     sum((V.bar.x+V.bar.y)*(V.bar.x+V.bar.y-diag(phi))),
    ##     2*I*theta.hat*(2*theta.hat-mean(diag(phi)))
    ## )
    c(
        cov((V.bar.x+V.bar.y),diag(phi)),
        var(V.bar.x+V.bar.y)
    )
    ## c(
    ##     sum((V.bar.x+V.bar.y)^2),
    ##     4*I*theta.hat^2
    ## )
    ## c(
    ##     sum((V.bar.x+V.bar.y)^2),
    ##     sum(V.bar.x+V.bar.y)^2/I
    ## )
})
plot(parts[1,],parts[2,]); abline(0,1)

concord.mean <- function(x,y)mean(outer(x,y,`<`))
concord.sum <- function(x,y)sum(outer(x,y,`<`))
I <- 10
m <- n <- 5
M <- I*m; N <- I*n
parts <- replicate(1e3,{
    x <- replicate(I,runif(m,0,10),simplify=FALSE)
    y <- replicate(I,rcauchy(n),simplify=FALSE)
    phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
    psi <- phi*m*n
    theta.hat <- concord.mean(unlist(x),unlist(y))
    ## theta.del <- sapply(1:I,function(i)concord.mean(unlist(x[-i]),unlist(y[-i])))
    V.x <- sapply(x,function(x.i)concord.sum(x.i,unlist(y)))/N
    V.y <- sapply(y,function(y.i)concord.sum(unlist(x),y.i))/M
    V.bar.x <- V.x/m; V.bar.y <- V.y/n
## (    sum((V.bar.x+V.bar.y)*(V.bar.x+V.bar.y-diag(phi))) -     2*I*theta.hat*(2*theta.hat-mean(diag(phi))) ) / I
        cov((V.bar.x+V.bar.y),diag(phi)) -        var(V.bar.x+V.bar.y)
})
plot(parts)



## 4. m=n=1


I <- 1000
parts <- replicate(1e3,{

    probs <- (1:(I+1))^3; probs <- probs/sum(probs)
    phi <- t(sapply(sample(0:I,I,replace=TRUE,prob=probs ), function(i)c(rep(0,I-i),rep(1,i))))
    ## phi <- matrix(sample(0:1,I^2,replace=TRUE,prob=c(.1,.9)),nrow=I)
    theta.hat <- sum(phi)/I^2
    V.x <- rowSums(phi)/I
    V.y <- colSums(phi)/I
    V.bar.x <- V.x; V.bar.y <- V.y
    sum((V.x + V.y - 2*theta.hat)*(V.x + V.y - diag(phi)))
    sum(((V.x + V.y - 2*theta.hat)*(V.x + V.y - diag(phi)))[diag(phi)==0]);     sum(((V.x + V.y - 2*theta.hat)*(V.x + V.y - diag(phi)))[diag(phi)==1])
    
    ## ## plot(V.x + V.y - 2*theta.hat,V.x + V.y - diag(phi),col=diag(phi)+1,xlim=c(-1,1)*2,ylim=c(-1,1)*2)
    ## cuts <- cut(abs((V.x + V.y - 2*theta.hat)*(V.x + V.y - diag(phi))),200)
    ## plot(V.x + V.y - 2*theta.hat  ,V.x + V.y - diag(phi),col=terrain.colors(length(levels(cuts)))[cuts])
    ## abline(0,1)
    ## cor((V.x + V.y - 2*theta.hat),(V.x + V.y - diag(phi)))
    ## prods <- (V.x + V.y - 2*theta.hat)*(V.x + V.y - diag(phi))
    ## hist(abs(prods)[prods>0])
    ## hist(abs(prods)[prods<=0],col=2,add=TRUE)
    ## plot(V.x + V.y - 2*theta.hat,V.x + V.y - diag(phi))
    ## plot((V.x + V.y - 2*theta.hat)*(V.x + V.y - diag(phi)))
    ## plot(V.x + V.y,V.x + V.y - diag(phi))
})
summary(parts)


plot(V.x + V.y - 2*theta.hat,V.x + V.y - diag(phi),col=diag(phi)+1,xlim=c(-1,1)*2,ylim=c(-1,1)*2)
cuts <- cut(abs((V.x + V.y - 2*theta.hat)*(V.x + V.y - diag(phi))),200)
plot(V.x + V.y - 2*theta.hat  ,V.x + V.y - diag(phi),col=rev(gray.colors(length(levels(cuts))))[cuts],ylim=c(-2,2),xlim=c(-2,2))
abline(0,1)
abline(2*theta.hat,1,lty=2)
abline(2*theta.hat-1,1,lty=3)
points(V.x + V.y - 2*theta.hat,V.x + V.y - 2*theta.hat)
abline(h=0);abline(v=0)


## V.x + V.y | diag(phi)==0 (stochastically increasing in row/col #)
I <- 10
parts <- replicate(1e3,{
    probs <- (1:(I+1))^0; probs <- probs/sum(probs)
    phi <- t(sapply(sample(0:I,I,replace=TRUE,prob=probs ), function(i)c(rep(0,I-i),rep(1,i))))
    ## phi <- matrix(sample(0:1,I^2,replace=TRUE,prob=c(.1,.9)),nrow=I)
    theta.hat <- sum(phi)/I^2
    V.x <- rowSums(phi)/I
    V.y <- colSums(phi)/I
    V.bar.x <- V.x; V.bar.y <- V.y
    out <- (V.x + V.y)
    out[diag(phi)==0] <- NA
    out
})
boxplot(t(parts))


I <- 10
parts <- replicate(1e3,{
    phi <- t(sapply(sample(0:I,I,replace=TRUE), function(i)c(rep(0,I-i),rep(1,i))))
    theta.hat <- sum(phi)/I^2
    V.x <- rowSums(phi)/I
    V.y <- colSums(phi)/I
    V.bar.x <- V.x; V.bar.y <- V.y
    a <- (2*I-1)/I/(I-1)^3 * sum( (V.bar.x + V.bar.y - 2*theta.hat)^2 )
    b <- 2/(I-1)^3 * sum( (V.bar.x+V.bar.y-2*theta.hat)*diag(phi) )
    c <- 1/I/(I-1)^3 * sum( (diag(phi)-mean(diag(phi)))^2 )
    c(a=a,b=b,c=c)
})
summary(abs(parts['a',] - parts['b',] + parts['c',]))
summary(abs(parts['a',] - parts['b',]))



Is <- seq(5,50,len=10)
by.I <- sapply(Is, function(I) {
    parts <- replicate(4e3,{
        phi <- t(sapply(sample(0:I,I,replace=TRUE), function(i)c(rep(0,I-i),rep(1,i))))
        theta.hat <- sum(phi)/I^2
        V.x <- rowSums(phi)/I
        V.y <- colSums(phi)/I
        V.bar.x <- V.x; V.bar.y <- V.y
        a <- (2*I-1)/I/(I-1)^3 * sum( (V.bar.x + V.bar.y - 2*theta.hat)^2 )
        a2 <- (2*I)/I/(I-1)^3 * sum( (V.bar.x + V.bar.y - 2*theta.hat)^2 )
        b <- 2/(I-1)^3 * sum( (V.bar.x+V.bar.y-2*theta.hat)*(diag(phi)-mean(diag(phi))) )
        c <- 1/I/(I-1)^3 * sum( (diag(phi)-mean(diag(phi)))^2 )
        a.tilde <- sum( (V.bar.x + V.bar.y - 2*theta.hat)^2 )
        b.tilde <- sum((V.bar.x+V.bar.y)*(diag(phi)-mean(diag(phi))))
        c(a=a,a2=a2,b=b,c=c,a.tilde=a.tilde,b.tilde=b.tilde)
    })
    c(max(abs(parts['a',] - parts['b',]+parts['c',])),
      max(abs(parts['a',] - parts['b',])))
})
plot(Is,by.I[1,],ylim=range(by.I),type='l')
lines(Is,by.I[2,],col=2)

plot(Is,by.I[1,]-by.I[2,])


## 4a. checking formulas in m=n=1 case, prob interpretation


I <- 10
b <- sort(sample(0:I,I,replace=TRUE))
b <- c(0,b,I)
f.B <- function(j)diff(b)[j]/I * (j>=1 & j<=I+1)
## f.B <- function(j) ( (b[j]-b[j-1])/I  )*(j>=1 & j<=I+1)
## f.B <- function(j) ( (b[j]-b[j-1])/I  )*(j>=1 & j<=I+1)
F.B <- function(j) b[j+1]/I
## sum((1:(I+1))*f.B(1:(I+1)))
## I+2-sum(b)/I
sum(F.B(1:(I+1))*f.B(1:(I+1)))
(sum(b^2) - sum(b[-1] * b[-(I+2)]))/I^2
sum(b[-1]*diff(b))/I^2


I <- 10
a <- sample(0:I,I,replace=TRUE)
b <- numeric(I+1)
b <- with(environment(ecdf(I+1-a)), {b[x] <- y; b})
for(j in seq_len(length(b)-1)+1) if(b[j]==0) b[j] <- b[j-1]
js <- 0:I; c(sum((I-js)^2*(b[js+1]-c(0,b[js])))*I,
             sum(a^2))




I <- 10
probs <- (1:(I+1))^-1; probs <- probs/sum(probs)
parts <- replicate(1e3, {
    a <- sample(0:I,I,replace=TRUE,prob=probs )
    phi <- t(sapply(a, function(i)c(rep(0,I-i),rep(1,i))))
    ## phi <- matrix(sample(0:1,I^2,replace=TRUE,prob=c(.1,.9)),nrow=I)
    theta.hat <- sum(phi)/I^2
    V.x <- rowSums(phi)/I
    V.y <- colSums(phi)/I
    V.bar.x <- V.x; V.bar.y <- V.y
    sum((V.x + V.y - 2*theta.hat)*(V.x + V.y - diag(phi)))
    sum(V.x^2) + 2*sum(V.x*V.y) + sum(V.y^2) - sum((V.x+V.y)*diag(phi)) - 4*I*theta.hat^2 + 2*theta.hat*sum(diag(phi))
    Is <- 1:I
    F <- ecdf(a)
    sum(a^2)/I^2 + 2/I*sum(a*(1-F(I-(Is)))) + sum((1-F(I-(Is)))^2) - 1/I*sum(a*(a>I-(Is))) - sum((1-F(I-(Is)))*(a>I-(Is))) - 4*sum(a)^2/I^3 + 2/I^2*sum(a)*sum(a>I-Is)
    ## c( 2/I*sum(a*(1-F(I-(Is)))) ,  2/I^2*sum(a)*sum(a>I-Is), sum(a^2)/I^2 , sum((1-F(I-(Is)))^2),- 1/I*sum(a*(a>I-(Is))), - sum((1-F(I-(Is)))*(a>I-(Is))), - 4*sum(a)^2/I^3 )
    c( 2/I*sum(a*(1-F(I-(Is)))) ,  2/I^2*sum(a)*sum(a>I-Is), sum(a^2)/I^2 , sum((1-F(I-(Is)))^2),- 1/I*sum(a*(a>I-(Is))), - sum((1-F(I-(Is)))*(a>I-(Is))), - 4*sum(a)^2/I^3 )[]
    ## |#3 + #5| <= 5*I/12 - 3/8 + 1/12/I
})
boxplot(t(parts))
summary(colSums(parts))
## rowMeans(parts)


I <- 10
probs <- (1:(I+1))^2; probs <- probs/sum(probs)
parts <- replicate(1e4, {
    a <- sample(0:I,I,replace=TRUE,prob=probs )
    phi <- t(sapply(a, function(i)c(rep(0,I-i),rep(1,i))))
    ## phi <- matrix(sample(0:1,I^2,replace=TRUE,prob=c(.1,.9)),nrow=I)
    theta.hat <- sum(phi)/I^2
    V.x <- rowSums(phi)/I
    V.y <- colSums(phi)/I
    V.bar.x <- V.x; V.bar.y <- V.y
    Is <- 1:I
    F <- ecdf(a)
    ## parts 3 and 5
    c( sum(a^2)/I^2 ,sum(a*(a>I-(Is)))/I )
    ## c( sum(a^1)/I^1 ,sum(a*(a>I-(Is)))/I )
    c( sum(V.x^2), sum(V.x*diag(phi)) )
    ## parts 4 and 6
    ## c( sum((1-F(I-(Is)))^2), sum((1-F(I-(Is)))*(a>I-(Is))))
    ## c( sum(V.y^2), sum(V.y*diag(phi)))
})
plot(parts[1,],parts[2,]); abline(0,1)
summary(abs(parts[1,]-parts[2,]))


## parts 3 and 5
phi <- rbind(t(sapply(0:(I/2), function(i)c(rep(0,i),rep(1,I-i)))),t(replicate(I/2-1,c(rep(0,I/2),rep(1,I/2)))))
diag(phi)[1:5] <- 0
## phi <- t(sapply(a, function(i)c(rep(0,I-i),rep(1,i))))
V.x <- rowSums(phi)/I
V.y <- colSums(phi)/I
## abs(sum(a^2)/I^2 -sum(a*(a>I-(Is)))/I )
## abs(sum(V.x^2) - sum(V.x*diag(phi)))
## sum((abs(V.x^2 - V.x*diag(phi)))[1:5])
## 7*I/24 - 3/8 + 1/12/I
sum((abs(V.x^2 - V.x*diag(phi)))[])
1/8*I + 7*I/24 - 3/8 + 1/12/I



## I=2, m=n=1
I <- 2
probs <- (1:(I+1))^-1; probs <- probs/sum(probs)
a <- sample(0:I,I,replace=TRUE,prob=probs )
phi <- t(sapply(a, function(i)c(rep(0,I-i),rep(1,i))))
## phi <- matrix(sample(0:1,I^2,replace=TRUE,prob=c(.1,.9)),nrow=I)
theta.hat <- sum(phi)/I^2
V.x <- rowSums(phi)/I
V.y <- colSums(phi)/I
V.bar.x <- V.x; V.bar.y <- V.y
## sum((V.x + V.y - 2*theta.hat)*(V.x + V.y - diag(phi)))
## sum(V.x^2) + 2*sum(V.x*V.y) + sum(V.y^2) - sum((V.x+V.y)*diag(phi)) - 4*I*theta.hat^2 + 2*theta.hat*sum(diag(phi))
Is <- 1:I
F <- ecdf(a)
x <- a[1]; y <- a[2]
    c( 2/I*sum(a*(1-F(I-(Is)))) ,  2/I^2*sum(a)*sum(a>I-Is), sum(a^2)/I^2 , sum((1-F(I-(Is)))^2),- 1/I*sum(a*(a>I-(Is))), - sum((1-F(I-(Is)))*(a>I-(Is))), - 4*sum(a)^2/I^3 )[]
-1/4*(x==1) + 3/4*(y==1) + (x==2) + (y==2) + 1/4*((min(x,y)==0) + (max(x,y)==0) + (min(x,y)<=1) + (max(x,y)<=1))
1/4*(x^2+y^2) + x + y - x/2*((x<=1) + (y<=1)) - y/2*((x<=0) + (y<=0)) + 1 - (x<=0) - (y<=0) + 1/4*((x<=0)+(y<=0))^2 + 1 - (x<=1) - (y<=1) +1/4*((x<=1)+(y<=1))^2 +    1/2*(x+y)*((x>1) + (y>0))- (x>1) + (x>1)*((x<=1) + (y<=1))/2 - (y>0) + (y>0)*((x<=0)+(y<=0))/2 - 1/2*(x^2+y^2+2*x*y) - 1/2*(x*(x>1)+y*(y>0))
-1/4*(x^2+y^2) - x*y +x  + y + 2 - 3/4*(x<=0)-3/4*(y<=0) + 1/2*(x<=0)*(y<=0) - 3/4*(x<=1)-3/4*(y<=1)+1/2*(x<=1)*(y<=1) + 1/2*((y>0)*(x<=0) + (x>1)*(y<=1)) - (x>1) - (y>0) + x/2*((y>0) - (x<=1) - (y<=1)) + y/2*((x>1) -(x<=0) - (y<=0))
-1/4*(x^2+y^2) +(x==2)+y/2*(y>0)+1/4*((x<=1)+(y==0)) - 3/4*((x<=0)+(y<=1)) + 1/2*((x==0)*(y==0)+(x<=1)*(y<=1)+(y!=0)*(x==0)+(x>1)*(y!=2))
-1/2*((x==0) + (y==0)+(y==1)) + 1/2*((x==0)*(y==0)+(x<=1)*(y<=1)+(y!=0)*(x==0)+(x>1)*(y!=2))




## 5. I=3, m=n=1
I <- 3
probs <- (1:(I+1))^-1; probs <- probs/sum(probs)
a <- sample(0:I,I,replace=TRUE,prob=probs )
phi <- t(sapply(a, function(i)c(rep(0,I-i),rep(1,i))))
## phi <- matrix(sample(0:1,I^2,replace=TRUE,prob=c(.1,.9)),nrow=I)
theta.hat <- sum(phi)/I^2
V.x <- rowSums(phi)/I
V.y <- colSums(phi)/I
V.bar.x <- V.x; V.bar.y <- V.y
## sum((V.x + V.y - 2*theta.hat)*(V.x + V.y - diag(phi)))
## sum(V.x^2) + 2*sum(V.x*V.y) + sum(V.y^2) - sum((V.x+V.y)*diag(phi)) - 4*I*theta.hat^2 + 2*theta.hat*sum(diag(phi))
Is <- 1:I
F <- ecdf(a)
sum(a^2)/I^2 + 2/I*sum(a*(1-F(I-(Is)))) + sum((1-F(I-(Is)))^2) + 2/I^2*sum(a)*sum(a>I-Is)- sum((1-F(I-(Is)))*(a>I-(Is))) - 4*sum(a)^2/I^3 - 1/I*sum(a*(a>I-(Is))) 
x <- a[1]; y <- a[2]; z <- a[3]
-1/27*(x^2+y^2+z^2) - 8/27*(x*y+x*z+y*z)+2/3*(x+y+z) + 3 - 2*sum(F(0:2)) + sum(F(0:2)^2)- (1-F(2))*(x>2)-(1-F(1))*(y>1)-(1-F(0))*(z>0) - 2/3*(x*F(2)+y*F(1)+z*F(0))  - 1/9*(x*(x>2)+y*(y>1)+z*(z>0)) + 2/9*(x*((y>1)+(z>0)) + y*((x>2)+(z>0)) + z*((x>2)+(y>1)))
## z*(-2/3 + 2/9*(x+y - 1 + (z>0)/2 + (x<=1) - (y==3))) + y*(-2/3 + 2/9*(x+z - 1 + (y>1)/2 + (x==0) - (z==3)) ) + x*(-2/3 + 2/9*(z+y - 2 + (x>2)/2 + (y==0) + (z<=1)))
## -1/27*(x^2+y^2+z^2) + 4/27*(x*y+x*z+y*z) + 3 - 2*sum(F(0:2)) + sum(F(0:2)^2)- (1-F(2))*(x>2)-(1-F(1))*(y>1)-(1-F(0))*(z>0) + z*2/9*( (z>0)/2 - (x>1) - (y==3)) + y*2/9*(  (y>1)/2 - (x>0) - (z==3))  + x*2/9*( (x>2)/2 - (y>0) - (z>1))
## js <- 0:2
## sum((1-F(js))^2) - (1-F(2))*(x>2)-(1-F(1))*(y>1)-(1-F(0))*(z>0) 
## sum(1-F(js)) - sum(a>I-js) + (y>2)*(z>2)+(x>1)*(z>1)+(x>0)*(y>0)
## sum((1-F(js))^2)
## 1/9*sum(((x>js)+(y>js)+(z>js))^2)
## 1/3*sum(1-F(js)) + 2/9*sum((x>js)*(y>js)+(x>js)*(z>js)+(y>js)*(z>js))
## sum((1-F(js))^2) - (1-F(2))*(x>2)-(1-F(1))*(y>1)-(1-F(0))*(z>0) 
## 1/3*sum(1-F(js)) - 1/3*sum(a>I-(1:3))+ 2/9*sum((x>js)*(y>js)+(x>js)*(z>js)+(y>js)*(z>js))  -2/9*(((y>2)+(z>2))*(x>2)+((x>1)+(z>1))*(y>1)+((x>0)+(y>0))*(z>0)) -1/9*(((y>2)+(z>2))*(x>2)+((x>1)+(z>1))*(y>1)+((x>0)+(y>0))*(z>0))
## 1/3*sum(1-F(js)) - 1/3*sum(a>I-(1:3))+ 2/9*((y>2)*(z>2) + (x>1)*(z>1) + (x>0)*(y>0))  -1/9*(((y>2)+(z>2))*(x>2)+((x>1)+(z>1))*(y>1)+((x>0)+(y>0))*(z>0))
## -1/27*(x^2+y^2+z^2) + 4/27*(x*y+x*z+y*z) +2/9*z*((z>0)/2-(x>1)-(y>2)) + 2/9*y*((y>1)/2-(x>0)-(z>2)) + 2/9*x*((x>2)/2-(y>0)-(z>1)) + 1/3*sum(1-F(0:2)) - 1/3*sum(a>I-(1:3)) +  2/9*((y>2)*(z>2) + (x>1)*(z>1) + (x>0)*(y>0))  -1/9*(((y>2)+(z>2))*(x>2)+((x>1)+(z>1))*(y>1)+((x>0)+(y>0))*(z>0))
-1/27*(x^2+y^2+z^2) + 4/27*(x*y+x*z+y*z) +2/9*z*((z>0)/2-(x>1)-(y>2)) + 2/9*y*((y>1)/2-(x>0)-(z>2)) + 2/9*x*((x>2)/2-(y>0)-(z>1)) + sum(a)/9 - 1/3*sum(a>I-(1:3)) +  2/9*((y>2)*(z>2) + (x>1)*(z>1) + (x>0)*(y>0))  -1/9*(((y>2)+(z>2))*(x>2)+((x>1)+(z>1))*(y>1)+((x>0)+(y>0))*(z>0))
## -1/27*(x^2+y^2+z^2) +1/9*(z*(z>0)+y*(y>1)+x*(x>2)) - 1/3*sum(a>I-(1:3)) + 1/9*(x+y+z)
## 1/27*(-4*(z==1)-(z==2)+2*(x==1)+2*(x==2)+2*(y==1)-(y==2))
 4/27*(x*y+x*z+y*z) +2/9*z*(-(x>1)-(y>2)) + 2/9*y*(-(x>0)-(z>2)) + 2/9*x*(-(y>0)-(z>1))  +  2/9*((y>2)*(z>2) + (x>1)*(z>1) + (x>0)*(y>0))  -1/9*(((y>2)+(z>2))*(x>2)+((x>1)+(z>1))*(y>1)+((x>0)+(y>0))*(z>0)) + 1/27*(-4*(z==1)-(z==2)+2*(x==1)+2*(x==2)+2*(y==1)-(y==2))
4/27*(x*y+x*z+y*z) +2/9*z*(-(x>1)-(y>2)) + 2/9*y*(-(x>0)-(z>2)) + 2/9*x*(-(y>0)-(z>1))  +  2/9*((y>2)*(z>2) + (x>1)*(z>1) + (x>0)*(y>0))  -1/9*(((y>2)+(z>2))*(x>2)+((x>1)+(z>1))*(y>1)+((x>0)+(y>0))*(z>0)) + 1/27*(-4*(z==1)-(z==2)+2*(x==1)+2*(x==2)+2*(y==1)-(y==2))
-2/27*(x*y+x*z+y*z) + 6/27*((x>1)*(y>2)+(x>2)*(y>1)+(x>0)*(z>2)+(x>2)*(z>0)+(y>1)*(z>0)+(y>0)*(z>1)) + 3/27*((x>1)*(y>1)+(x>2)*(y>2)+(x>0)*(z>0)+(x>2)*(z>2)+(y>0)*(z>0)+(y>1)*(z>1))+ 1/27*(-4*(z==1)-(z==2)+2*(x==1)+2*(x==2)+2*(y==1)-(y==2))


## 5a. matrix part
require(abind)
I <- 3
probs <- (1:(I+1))^-1; probs <- probs/sum(probs)
a <- sample(0:I,I,replace=TRUE,prob=probs )
x <- a[1]; y <- a[2]; z <- a[3]
-2/27*(x*y+x*z+y*z) + 6/27*((x>1)*(y>2)+(x>2)*(y>1)+(x>0)*(z>2)+(x>2)*(z>0)+(y>1)*(z>0)+(y>0)*(z>1)) + 3/27*((x>1)*(y>1)+(x>2)*(y>2)+(x>0)*(z>0)+(x>2)*(z>2)+(y>0)*(z>0)+(y>1)*(z>1))#+ 1/27*(-4*(z==1)-(z==2)+2*(x==1)+2*(x==2)+2*(y==1)-(y==2))
m.xy <- matrix(c(-2,-2,-2,-2,1,4,-2,4,1),nrow=3,byrow=TRUE)
m.xz <- matrix(c(1,-2,4,-2,-2,-2,4,-2,1),nrow=3,byrow=TRUE)
m.yz <- matrix(c(1,4,-2,4,1,-2,-2,-2,-2),nrow=3,byrow=TRUE)
v.x <- x > (0:2); v.y <- y>(0:2); v.z <- z>(0:2)
(v.x%*%m.xy%*%v.y + v.x%*%m.xz%*%v.z + v.y%*%m.yz%*%v.z)/27
## xy <- matrix(c(-2,-4,-6,-4,-5,-3,-6,-3,0),nrow=3,byrow=TRUE)
## m <- abind(xy+matrix(c(2,6,4,0,4,2,4,8,6),nrow=3,byrow=TRUE),
## xy+matrix(c(4,9,5,0,5,-1,2,7,3),nrow=3,byrow=TRUE),
## xy+matrix(c(6,9,3,0,3,-3,3,6,0),nrow=3,byrow=TRUE),along=3)
## m[a[1],a[2],a[3]]/27
xy <- matrix(c(-2,-4,-6,-4,-5,-3,-6,-3,0),nrow=3,byrow=TRUE)
xy <- rbind(rep(0,4),cbind(rep(0,3),xy))
m <- abind(xy,
           xy+matrix(c(0,1,5,3,1,2,6,4,-1,0,4,2,-3,4,8,6),nrow=4,byrow=TRUE),
           xy+matrix(c(0,5,10,6,-1,4,9,5,-5,0,5,-1,-3,2,7,3),nrow=4,byrow=TRUE),
           xy+matrix(c(0,3,6,0,3,6,9,3,-3,0,3,-3,0,3,6,0),nrow=4,byrow=TRUE),along=3)
m[a[1]+1,a[2]+1,a[3]+1]/27



## 5b. augmented matrix 
require(abind)
I <- 3
probs <- (1:(I+1))^-1; probs <- probs/sum(probs)
a <- sample(0:I,I,replace=TRUE,prob=probs )
x <- a[1]; y <- a[2]; z <- a[3]
-2/27*(x*y+x*z+y*z) + 6/27*((x>1)*(y>2)+(x>2)*(y>1)+(x>0)*(z>2)+(x>2)*(z>0)+(y>1)*(z>0)+(y>0)*(z>1)) + 3/27*((x>1)*(y>1)+(x>2)*(y>2)+(x>0)*(z>0)+(x>2)*(z>2)+(y>0)*(z>0)+(y>1)*(z>1))+ 1/27*(-4*(z==1)-(z==2)+2*(x==1)+2*(x==2)+2*(y==1)-(y==2))
xy <- matrix(c(-2,-4,-6,-4,-5,-3,-6,-3,0),nrow=3,byrow=TRUE)
xy <- rbind(rep(0,4),cbind(rep(0,3),xy))
m <- abind(xy,
           xy+matrix(c(0,1,5,3,1,2,6,4,-1,0,4,2,-3,4,8,6),nrow=4,byrow=TRUE),
           xy+matrix(c(0,5,10,6,-1,4,9,5,-5,0,5,-1,-3,2,7,3),nrow=4,byrow=TRUE),
           xy+matrix(c(0,3,6,0,3,6,9,3,-3,0,3,-3,0,3,6,0),nrow=4,byrow=TRUE),along=3)
m[,,2] <- m[,,2]-4; m[,,3] <- m[,,3]-1
m[2,,] <- m[2,,]+2; m[3,,] <- m[3,,]+2
m[,2,] <- m[,2,]+2; m[,3,] <- m[,3,]-1
m[a[1]+1,a[2]+1,a[3]+1]/27
(x*y+x*z+y*z)/27


## I=4
I <- 3
probs <- (1:(I+1))^0; probs <- probs/sum(probs)
as <- replicate(1e4,sample(0:I,I,replace=TRUE,prob=probs ))
diffs <- sapply(1:ncol(as), function(jj){
    a <- as[,jj]
    ## phi <- t(sapply(a, function(i)c(rep(0,I-i),rep(1,i))))
    ## ## phi <- matrix(sample(0:1,I^2,replace=TRUE,prob=c(.1,.9)),nrow=I)
    ## theta.hat <- sum(phi)/I^2
    ## V.x <- rowSums(phi)/I
    ## V.y <- colSums(phi)/I
    ## V.bar.x <- V.x; V.bar.y <- V.y
    ## ## sum((V.x + V.y - 2*theta.hat)*(V.x + V.y - diag(phi)))
    ## ## sum(V.x^2) + 2*sum(V.x*V.y) + sum(V.y^2) - sum((V.x+V.y)*diag(phi)) - 4*I*theta.hat^2 + 2*theta.hat*sum(diag(phi))
    Is <- 1:I
    F <- ecdf(a)
    I^2*(sum(a^2)/I^2 + 2/I*sum(a*(1-F(I-(Is)))) + sum((1-F(I-(Is)))^2) + 2/I^2*sum(a)*sum(a>I-Is)- sum((1-F(I-(Is)))*(a>I-(Is))) - 4*sum(a)^2/I^3 - 1/I*sum(a*(a>I-(Is))) )
})
## 3:8/27=(2/3)^3, 4:9/16=(3/4)^2, 5:96/125
## for(d in 1:2000) if(sum(round(unique(diffs)*d,6) %% 1)==0)print(d)
## idx <- which(abs(diffs)==max(abs(diffs)))
idx <- which(abs(abs(diffs)-max(abs(diffs)))<1e-6)
unique(apply(as.matrix(as[,idx]),2,paste,collapse=''))

## 5c. full I^2 x I^2 matrix
require(Matrix)
I <- 3
E <- diag(I); ones <- matrix(1,nrow=I,ncol=1)
e <- lapply(1:I,function(i)matrix(E[i,],ncol=1))
Exch <- E[I:1,]
m.a <-  kronecker(E,ones%*%t(ones))
m.b <- kronecker(t(ones),kronecker(Exch,ones))
m.c <- kronecker(ones%*%t(ones),E)
m.d <- kronecker(matrix(as.numeric(Exch),nrow=1),kronecker(ones,ones))
m.e <- do.call(rbind,lapply(I:1,function(i)kronecker(t(ones),e[[i]]%*%t(e[[i]]))))
m.f <- kronecker(ones%*%t(ones),ones%*%t(ones))
m.g <- as.matrix(bdiag(apply(Exch,1,function(r)kronecker(matrix(r,nrow=1),ones),simplify=FALSE)))
probs <- (1:(I+1))^0; probs <- probs/sum(probs)
a <- sample(0:I,I,replace=TRUE,prob=probs )
Is <- 1:I
F <- ecdf(a)
## I^2*c(   sum(a^2)/I^2 ,+ 2/I*sum(a*(1-F(I-(Is)))), + sum((1-F(I-(Is)))^2), + 2/I^2*sum(a)*sum(a>I-Is),- sum((1-F(I-(Is)))*(a>I-(Is))), - 4*sum(a)^2/I^3 ,- 1/I*sum(a*(a>I-(Is)))   )
v.lst <- lapply(1:I,function(i)as.numeric(a[i]>=1:I))
v <- matrix(unlist(v.lst),ncol=1)
m.lst <- list(a=m.a,b=2*m.b,c=m.c,d=2*m.d,e=-I*m.e,f=-4/I*m.f,g=-I*m.g)
## sapply(m.lst,function(mm)t(v)%*%mm%*%v)
m.sum <- Reduce(`+`,m.lst)



## 5d. wrapper
M <- function(I,partition=FALSE) {
    E <- diag(I); ones <- matrix(1,nrow=I,ncol=1)
    e <- lapply(1:I,function(i)matrix(E[i,],ncol=1))
    Exch <- E[I:1,]
    m.a <-  kronecker(E,ones%*%t(ones))
    m.b <- kronecker(t(ones),kronecker(Exch,ones))
    m.c <- kronecker(ones%*%t(ones),E)
    m.d <- kronecker(matrix(as.numeric(Exch),nrow=1),kronecker(ones,ones))
    m.e <- do.call(rbind,lapply(I:1,function(i)kronecker(t(ones),e[[i]]%*%t(e[[i]]))))
    m.f <- kronecker(ones%*%t(ones),ones%*%t(ones))
    m.g <- as.matrix(bdiag(apply(Exch,1,function(r)kronecker(matrix(r,nrow=1),ones),simplify=FALSE)))
    ans <- Reduce(`+`,list(a=m.a,b=2*m.b,c=m.c,d=2*m.d,e=-I*m.e,f=-4/I*m.f,g=-I*m.g))
    if(partition) {
        idx.lst <- split(1:(I^2),rep(1:I,each=I))
        ans <- apply(expand.grid(1:I,1:I),1,function(x)ans[idx.lst[[x[2]]],idx.lst[[x[1]]]],simplify=FALSE)
    }
    return(ans)
}


I <- 4
probs <- (1:(I+1))^0; probs <- probs/sum(probs)
compare <- replicate(1e3, {
    a <- sample(0:I,I,replace=TRUE,prob=probs )
    Is <- 1:I
    F <- ecdf(a)
    true <- sum(I^2*c(   sum(a^2)/I^2 ,+ 2/I*sum(a*(1-F(I-(Is)))), + sum((1-F(I-(Is)))^2), + 2/I^2*sum(a)*sum(a>I-Is),- sum((1-F(I-(Is)))*(a>I-(Is))), - 4*sum(a)^2/I^3 ,- 1/I*sum(a*(a>I-(Is)))   ))
    v.lst <- lapply(1:I,function(i)as.numeric(a[i]>=1:I))
    v <- matrix(unlist(v.lst),ncol=1)
    m.lst <- list(a=m.a,b=2*m.b,c=m.c,d=2*m.d,e=-I*m.e,f=-4/I*m.f,g=-I*m.g)
    ## sapply(m.lst,function(mm)t(v)%*%mm%*%v)
    ## sum(sapply(m.lst,function(mm)t(v)%*%mm%*%v))
    ## Reduce(`+`,m.lst)
    test <- t(v)%*%M(I)%*%v
    c(true,test)
})
max(abs(compare[1,]-compare[2,]))

I <- 7
m <- M(I)
## as <- replicate(1e5,sample(0:I,I,replace=TRUE ))
as <- t(expand.grid(rep(list(0:I),I)))
diffs <- sapply(1:ncol(as), function(jj) {
    a <- as[,jj]
    ## v.lst <- lapply(1:I,function(i)as.numeric(a[i]>=1:I))
    ## v <- matrix(unlist(v.lst),ncol=1)
    v <- as.matrix(as.numeric(sapply(a,function(a.i)a.i>=(1:I))),ncol=1)
    t(v)%*%m%*%v 
    ## v <- as.numeric(sapply(a,function(a.i)a.i>=(1:I)))
    ## sum(m[,v][v,])
    ## diff <- 0
    ## for(i in 1:I)
    ##     for(j in 1:I) {
    ##         idx <- I*(i-1)+j
    ##         diff <- diff + sum(m[[idx]][1:a[i],1:a[j]])
    ##     }
    ## diff
})
idx <- which(abs(abs(diffs)-max(abs(diffs)))<1e-6)
unique(apply(as.matrix(as[,idx]),2,paste,collapse=''))


## phi totally positive? no.
## det(phi) gets very small as I grows
concord.mean <- function(x,y)mean(outer(x,y,`<`))
I <- 5
m <- n <- 10
M <- I*m; N <- I*n
dets <- replicate(1e3,{
x <- replicate(I,runif(m),simplify=FALSE)
y <- replicate(I,runif(n),simplify=FALSE)
phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
psi <- phi*m*n
det(phi[sample(1:I,2),sample(1:I,2)])
det(phi)
})
summary(dets)

I <- 10
m <- n <- 10
M <- I*m; N <- I*n
x <- replicate(I,runif(m),simplify=FALSE)
y <- replicate(I,runif(n),simplify=FALSE)
phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
sapply(2:I,function(i)0+(phi[1,]<phi[i,]))


## diff in I>1 case <= I=1 case?
source('misc.R')
I <- 5
m <- n <- 10
M <- I*m; N <- I*n
obj <- replicate(1e4, {
x <- replicate(I,runif(m),simplify=FALSE)
y <- replicate(I,runif(n),simplify=FALSE)
phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
## phi <- t(sapply(a, function(i)c(rep(0,I-i),rep(1,i))))
## phi <- matrix(sample(0:1,I^2,replace=TRUE,prob=c(.1,.9)),nrow=I)
theta.hat <- sum(phi)/I^2
V.x <- rowSums(phi)/I
V.y <- colSums(phi)/I
V.bar.x <- V.x; V.bar.y <- V.y
sum((V.x + V.y - 2*theta.hat)*(V.x + V.y - diag(phi)))
})
max(abs(obj))


## 7 [out of order] quad form with lower order term


I <- 3
E <- diag(I); ones <- matrix(1,nrow=I,ncol=1)
m.a <- kronecker(E,ones%*%t(ones))
m.b <- kronecker(ones%*%t(ones),E)
m.c <- kronecker(ones,kronecker(E,t(ones)))
m.d <- Reduce(`+`,lapply(1:I,function(j)kronecker(matrix(E[j,],nrow=1),kronecker(kronecker(E[,j],E[,j]),t(ones)))))
m.e <- Reduce(`+`,lapply(1:I,function(j)kronecker(E[,j],kronecker(t(ones),E[,j]%*%t(E[,j])))))
m.f <- kronecker(ones%*%t(ones),ones%*%t(ones))
m.g <- kronecker(matrix(E,nrow=1),kronecker(ones,ones))
m.lst <- list(a=m.a/I^2,b=m.b/I^2,c=m.c*2/I^2,d=-1/I*m.d,e=-1/I*m.e,f=-4/I^3*m.f,g=2/I^2*m.g)
ans.1 <- Reduce(`+`,m.lst)
## if(symmetric) ans <- (ans+t(ans))/2
m.h <- Reduce(`+`,lapply(1:I,function(j)kronecker(E[,j],E[,j])%*%t(kronecker(E[,j],E[,j]))))
m.i <- as.numeric(E)%*%t(as.numeric(E))
m2.lst <- list(a=-1/I^3*m.a,b=-1/I^3*m.b,c=-2/I^3*m.c,f=4/I^4*m.f,h=1/I*m.h,i=-1/I^2*m.i)
ans.2 <- Reduce(`+`,m2.lst)
m.lst <- list(a=m.a*(1/I^2-1/I^3),b=m.b*(1/I^2-1/I^3),c=m.c*(2/I^2-2/I^3),d=-1/I*m.d,e=-1/I*m.e,f=(-4/I^3+4/I^4)*m.f,g=2/I^2*m.g,h=1/I*m.h,i=-1/I^2*m.i)
ans <- Reduce(`+`,m.lst)

I <- 3
m <- n <- 5
M <- I*m; N <- I*n
## parts <- replicate(1e2,{
x <- replicate(I,runif(m),simplify=FALSE)
y <- replicate(I,rnorm(n),simplify=FALSE)
phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
psi <- phi*m*n
phi.vec <- matrix(t(phi),ncol=1)
theta.hat <- concord.mean(unlist(x),unlist(y))
V.x <- sapply(x,function(x.i)concord.sum(x.i,unlist(y)))/N
V.y <- sapply(y,function(y.i)concord.sum(unlist(x),y.i))/M
V.bar.x <- V.x/m; V.bar.y <- V.y/n
c(a=sum(V.bar.x^2),b=sum(V.bar.y^2),c=2*sum(V.bar.x*V.bar.y),d=-sum(diag(phi)*V.bar.x),e=-sum(diag(phi)*V.bar.y),f=-2*theta.hat*sum(V.bar.x+V.bar.y),g=2*theta.hat*sum(diag(phi)),h=sum(diag(phi)^2)/I,i=(sum(diag(phi))/I)^2)
## })

sum((V.bar.x+V.bar.y-2*theta.hat)*(V.bar.x+V.bar.y-diag(phi)))
t(phi.vec)%*%ans.1%*%phi.vec

1/I*sum((diag(phi)-mean(diag(phi)))^2 - (V.bar.x+V.bar.y-2*theta.hat)^2)
## 1/I*sum(diag(phi)^2)-mean(diag(phi))^2 - mean((V.bar.x+V.bar.y)^2)+4*theta.hat^2
t(phi.vec)%*%ans.2%*%phi.vec

sum((V.bar.x+V.bar.y-2*theta.hat)*(V.bar.x+V.bar.y-diag(phi)))+1/I*sum((diag(phi)-mean(diag(phi)))^2 - (V.bar.x+V.bar.y-2*theta.hat)^2)
t(phi.vec)%*%(ans.1+ans.2)%*%phi.vec
t(phi.vec)%*%(ans)%*%phi.vec


## 6 m,n fixed but >= 1

## 6a quad form 
source('misc.R')

partition <- function(m) {
    I <- as.integer(sqrt(nrow(m)))
    idx.lst <- split(1:(I^2),rep(1:I,each=I))
    apply(expand.grid(1:I,1:I),1,function(x)m[idx.lst[[x[2]]],idx.lst[[x[1]]]],simplify=FALSE)
}
M <- function(I,partition=FALSE,symmetric=FALSE) {
    E <- diag(I); ones <- matrix(1,nrow=I,ncol=1)
    m.a <- kronecker(E,ones%*%t(ones))
    m.b <- kronecker(ones%*%t(ones),E)
    m.c <- kronecker(ones,kronecker(E,t(ones)))
    m.d <- Reduce(`+`,lapply(1:I,function(j)kronecker(matrix(E[j,],nrow=1),kronecker(kronecker(E[,j],E[,j]),t(ones)))))
    m.e <- Reduce(`+`,lapply(1:I,function(j)kronecker(E[,j],kronecker(t(ones),E[,j]%*%t(E[,j])))))
    m.f <- kronecker(ones%*%t(ones),ones%*%t(ones))
    m.g <- kronecker(matrix(E,nrow=1),kronecker(ones,ones))
    m.lst <- list(a=m.a/I^2,b=m.b/I^2,c=m.c*2/I^2,d=-1/I*m.d,e=-1/I*m.e,f=-4/I^3*m.f,g=2/I^2*m.g)
    ans <- Reduce(`+`,m.lst)
    if(symmetric) ans <- (ans+t(ans))/2
    ## sapply(m.lst,function(m)t(v)%*%m%*%v)
    if(partition) {
        ans <- partition(ans)
    }
    ans
}

I <- 5
m <- n <- 5
M <- I*m; N <- I*n
## parts <- replicate(1e2,{
x <- replicate(I,runif(m),simplify=FALSE)
y <- replicate(I,rnorm(n),simplify=FALSE)
phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
psi <- phi*m*n
phi.vec <- matrix(t(phi),ncol=1)
theta.hat <- concord.mean(unlist(x),unlist(y))
V.x <- sapply(x,function(x.i)concord.sum(x.i,unlist(y)))/N
V.y <- sapply(y,function(y.i)concord.sum(unlist(x),y.i))/M
V.bar.x <- V.x/m; V.bar.y <- V.y/n
c(a=sum(V.bar.x^2),b=sum(V.bar.y^2),c=2*sum(V.bar.x*V.bar.y),d=-sum(diag(phi)*V.bar.x),e=-sum(diag(phi)*V.bar.y),f=-2*theta.hat*sum(V.bar.x+V.bar.y),g=2*theta.hat*sum(diag(phi)))
## })

sum((V.bar.x+V.bar.y-2*theta.hat)*(V.bar.x+V.bar.y-diag(phi)))
t(phi.vec)%*%M(I)%*%phi.vec

## M symmetric on vertical + horizontal flips
sum(abs(M(I)[I^2:1,][,I^2:1] - M(I)))
## row sums = - col sums
sum(abs(rowSums(M(I))+colSums(M(I))))



## symmetric form Q
source('misc.R')
I <- 3
eigen((M(I)+t(M(I))))$val
for(I in 2:40)print(sprintf('%d: %f',I,max(eigen((M(I)+t(M(I))))$val)))
## ...M(2)+t(M(2))==0. 2*(I-1) nonzero eigenvalues, all equal in
## magnitude, half positive and half neg
## evectors for nonzero evals sum to 0 (since rowsums==0):
colSums(eigen(M.obj(I,symm=TRUE))$vec)[c(1:(I-1),(I^2-I+2):(I^2))]
## each IxI block sums to zero, whether symm=TRUE or not:
for(I in 2:20)print(max(abs(sapply(M(I,partition=TRUE),sum))))
for(I in 2:20)print(max(abs(sapply(M(I,partition=TRUE,symmetric=TRUE),sum))))
## not block row/col sums though (but row/col sums of entire matrices do):
sapply(M(I,partition=TRUE,symm=TRUE),rowSums)*I^2
## IxI blocks all symmetric (untrue for symm=FALSE):
for(I in 2:20)print(sum(sapply(M(I,partition=TRUE,symmetric=TRUE),function(m)sum(abs(m-t(m))))))
for(I in 2:20)print(sum(sapply(M(I,partition=TRUE,symmetric=FALSE),function(m)sum(abs(m-t(m))))))
## symmtric about anti-diagonal also
m <- M.obj(I,sym=TRUE)
sum(abs(do.call(cbind,apply(m,2,rev,simplify=FALSE))-do.call(rbind,apply(m,2,rev,simplify=FALSE))))
## t(kronecker(P,P))%*%Q%*%kronecker(P,P) == Q, for permutation mtx
## P. Same as permuting clusters.
I <- 5
Q <- M.obj(I,sym=TRUE)
P <- do.call(cbind,unname((sapply(sample(I),function(i)within(list(v=rep(0,I)),v[i] <- 1))))) # IxI perm matrix
sum(abs(t(kronecker(P,P))%*%Q%*%kronecker(P,P) - Q))




I <- 5
M0 <- M(I,symm=TRUE)*I^3*2
P <- matrix(runif(I^2),nrow=I)
P <- kronecker(diag(I),P/rowSums(P))
M1 <- t(P)%*%M0%*%P
max(sapply(partition(M1),function(m)sum(abs(m-t(m)))))
max(sapply(partition(M1),function(m)sum(m)))
## uniform P gives 0 matrix:
P <- kronecker(diag(I),matrix(1/I,nrow=I,ncol=I))
t(P)%*%M0%*%P
## identity P
P <- kronecker(diag(I),diag(I)[I:1,])
t(P)%*%M0%*%P
## blocks are indefinite
sapply(M(I,symm=TRUE,part=TRUE),function(m)eigen(m)$val)
## blocks have zero trace
sapply(M(I,symm=TRUE,part=TRUE),function(m)sum(eigen(m)$val))

## shift symmetry -- just shifting indices of clusers
source('misc.R')
Q <- M.obj(I,sym=TRUE)
I <- 5
diffs <- replicate(1e2,{
s <- sample(1:(I-1),1)
phi <- matrix(runif(I^2),nrow=I)
phi.vec <- as.numeric(t(phi))
phi.vec2 <- as.numeric(t(phi[c((I-s+1):I,1:(I-s)),c((I-s+1):I,1:(I-s))]))
t(phi.vec)%*%Q%*%phi.vec - t(phi.vec2)%*%Q%*%phi.vec2
}) 
sum(abs(diffs))

dd

## 2I phi.vecs mutually orthogonal in Q basis. Correspond to constant
## column and constant row phis. follows from blocks summing to 0 and
## Q(p,q,r,s) symmetry in arguments. shows row and col sums of Q are
## 0. Also includes constant diagonal. [update: actually in nullspace]
source('misc.R')
I <- 5
E <- diag(I); ones <- matrix(1,nrow=I,ncol=1)
phi.vecs <- cbind(kronecker(E,ones),kronecker(ones,E))
sum(abs(t(phi.vecs)%*%M(I,sym=TRUE)%*%phi.vecs))
m <- Q<- M.obj(I,sym=TRUE)
a <- kronecker(E[,3],ones)
b <- kronecker(E[,2],ones)
t(a)%*%m%*%b ## blocks sum to 0
c <- kronecker(ones,E[,3]) ## r,s/p,q symmetry
d <- kronecker(ones,E[,2]) ## r,p/p,q symmetry
t(c)%*%Q%*%d
t(a)%*%Q%*%c
qq <- Q%*%phi.vecs
qq[,1:I]-qq[,(I+1):(2*I)]
## I vectors in nullspace of Q
sum(abs(Q%*%(phi.vecs[,1:I]-phi.vecs[,(I+1):(2*I)])))
vv <- (eigen(M.obj(I,symm=TRUE))$vec)[,c(1:(I-1),(I^2-I+2):(I^2))]
sapply(1:(2*I-2),function(k)sapply(1:I,function(i)sum(vv[(I*(i-1)):(I*i),k])))

phi <- matrix(runif(I^2),nrow=I)
t(as.numeric(phi))%*%Q%*%as.numeric(phi)
t(as.numeric(t(phi)))%*%Q%*%as.numeric(t(phi))
rbind(as.numeric(phi),as.numeric(t(phi)))

## remains Q-orthogonal on adding constant diagonal matrices 
phi.vec <- c(diag(I))
phi.vec%*%Q%*%phi.vec
phi.vecs <- cbind(kronecker(E,ones),kronecker(ones,E))
sum(abs(t(phi.vecs)%*%Q%*%phi.vecs))
phi.vecs <- cbind(phi.vecs,c(diag(I)))
sum(abs(t(phi.vecs)%*%Q%*%phi.vecs))


dd


## WW is an identity for Q, so ww-E in the nullspace. ww-E has I(I-1)
## nonzero columns, and for each column its negation is also present,
## so at most I(I-1)/2 dimensional subspace of the nullspace of Q in
## ww-E. [show actually equal to I(I-1)/2]. The fact that ww is a
## right identity for Q is just the symmetry of Q[p,q,r,s] in q,s, fact
## that it is a left identity is the p,r, symmetry, combined with ww
## being symmetric. As above, these symmetries related to exchanging
## x,y roles in the data.
source('misc.R')
I <- 5
Q <- M.obj(I,sym=TRUE)
pos.vecs <- with(eigen(Q),vectors[,values>1e-7])
neg.vecs <- with(eigen(Q),vectors[,values< -1e-7])
N <- with(eigen(Q),vectors[,abs(values)<1e-7])
e <- Vectorize(function(j)unlist(within(list(v=rep(0,I^2)),v[j] <- 1)))
flip <- function(v)as.numeric(matrix(v,nrow=I,byrow=TRUE))
j <- sample(1:I,1); k <- sample(1:I,1)
t(e(j))%*%Q%*%e(j)
t(flip(e(j)))%*%Q%*%flip(e(j))
## t(flip(e(1)))%*%Q%*%flip(e(5))
## t(e(1))%*%Q%*%e(5)
t(flip(e(j))+flip(e(k)))%*%Q%*%(flip(e(j))+flip(e(k)))
t(e(j)+e(k))%*%Q%*%(e(j)+e(k))
c(t(flip(e(j)))%*%Q%*%flip(e(j)),t(flip(e(k)))%*%Q%*%flip(e(k)),t(flip(e(j)))%*%Q%*%flip(e(k)))
c(t(e(j))%*%Q%*%e(j),t(e(k))%*%Q%*%e(k),t(e(j))%*%Q%*%e(k))
## t(e(flip(j)))%*%Q%*%e(j) + t(e(flip(k)))%*%Q%*%e(k)
ww <- t(sapply(1:(I^2),function(i)flip(e(i))))
vv <- diag(I^2)
sum(abs(Q%*%ww-Q)); sum(abs(ww%*%Q-Q))
## t(ww)%*%Q%*%ww - t(vv)%*%Q%*%vv
## t(ww)%*%Q%*%vv - Q
## Q%*%(ww-vv)
## vv <- cbind(diag(1:(I^2)),sapply(1:(I^2),function(i)flip(e(i))))
## t(vv)%*%Q%*%vv
## 0,1,3,6,10,15,21,28,36 ie j(j-1)/2

## phi.vecs defined above equally inclinded to pos and neg eigenvectors
## phi.vecs corresponding to all 1s on row or col of phi
source('misc.R')
I <- 3
Q <- M.obj(I,sym=TRUE)
E <- diag(I); ones <- matrix(1,nrow=I,ncol=1)
pos.vecs <- with(eigen(Q),vectors[,values>1e-7])
neg.vecs <- with(eigen(Q),vectors[,values< -1e-7])
N <- with(eigen(Q),vectors[,abs(values)<1e-7])
e <- Vectorize(function(j)unlist(within(list(v=rep(0,I^2)),v[j] <- 1)))
phi.vecs <- cbind(kronecker(E,ones),kronecker(ones,E))
colSums((t(N)%*%phi.vecs)^2)
colSums((t(pos.vecs)%*%phi.vecs)^2)
colSums((t(neg.vecs)%*%phi.vecs)^2)
## phi.vecs as initial row segments 
M.cumsum <- matrix(0,nrow=I,ncol=I); M.cumsum <- (row(M.cumsum)<=col(M.cumsum))+0
phi.vecs <- kronecker(E, M.cumsum)
## phi.vecs <- sapply(1:I^2,function(i)rowSums(e(1:i)))
projs <-     cbind(N=colSums((t(N)%*%phi.vecs)^2),
pos=colSums((t(pos.vecs)%*%phi.vecs)^2),
neg=colSums((t(neg.vecs)%*%phi.vecs)^2))
projs
apply(projs,2,diff) #not all increasing...
op <- par(mfrow=c(1,3))
plot(projs[,'neg']);abline(0,1)
plot(projs[,'pos']);abline(0,1)
plot(projs[,'N']);abline(0,1)
par(op)
## phi.vecs as initial col segments. same as projections of initial
## row segments.
projs.row <- projs
phi.vecs <- do.call(cbind,lapply(I:1,function(r)sapply(I:1,function(c) {
    phi <- matrix(0,nrow=I,ncol=I)
    phi[r,I:c] <- 1
    as.numeric(phi)
})))
projs <-     cbind(N=colSums((t(N)%*%phi.vecs)^2),
pos=colSums((t(pos.vecs)%*%phi.vecs)^2),
neg=colSums((t(neg.vecs)%*%phi.vecs)^2))
projs
apply(projs,2,diff) # increasing within each block
op <- par(mfrow=c(1,3))
plot(projs[,'neg']);abline(0,1)
plot(projs[,'pos']);abline(0,1)
plot(projs[,'N']);abline(0,1)
par(op)
sum(abs(projs.row - projs))
## sum of all projections onto pos eigenspace equal to sum onto neg
## eigenspace, but different spreads
boxplot(projs.row)



## Q in general invariant on applying the same permutation to both
## rows and cols of phi--this is just permuting the indices of he
## clusters.
sum(abs(replicate(1e2,{
    I <- sample(3:10,1)
    Q <- M.obj(I,symm=TRUE)
    phi <- matrix(runif(I^2),I)
    v <- as.numeric(phi)
    perm <- sample(1:I)
    v.perm <- as.numeric(phi[perm,perm])
    t(v)%*%Q%*%v - t(v.perm)%*%Q%*%v.perm
})))

dd


## looking further at pos and neg eigenspaces. off diagonal blocks of
## pos.qf and neg.qf all have sum -1/4, for all I, diagonal blocks sum
## to I/4-1/4
Is <- 3:10; names(Is) <- Is
sapply(Is, function(I) {
    Q <- M.obj(I,sym=TRUE)
    pos.vecs <- with(eigen(Q),vectors[,values>1e-7])
    neg.vecs <- with(eigen(Q),vectors[,values< -1e-7])
    neg.qf <- neg.vecs%*%t(neg.vecs)
    pos.qf <- pos.vecs%*%t(pos.vecs)
    c(sum(abs(unlist(lapply(partition(neg.qf),sum)[(1:(I^2)-1)%%(I+1)!=0]) - -.25)),
      sum(abs(unlist(lapply(partition(pos.qf),sum)[(1:(I^2)-1)%%(I+1)!=0]) - -.25)),
      sum(abs(unlist(lapply(partition(pos.qf),sum)[(1:(I^2)-1)%%(I+1)==0]) - (I/4-1/4))))
})

## ratio of pos.qf to neg.qf has simple structure
source('misc.R')
I <- 3
Q <- M.obj(I,sym=TRUE)
E <- diag(I); ones <- matrix(1,nrow=I,ncol=1)
pos.vecs <- with(eigen(Q),vectors[,values>1e-7])
neg.vecs <- with(eigen(Q),vectors[,values< -1e-7])
## a <- kronecker(E[,3],ones)
## b <- kronecker(E[,2],ones)
## c <- kronecker(ones,E[,3]) 
## d <- kronecker(ones,E[,2]) 
## sum((t(pos.vecs)%*%a)^2)
## sum((t(neg.vecs)%*%a)^2)
## sum((t(pos.vecs)%*%b)^2)
## sum((t(neg.vecs)%*%b)^2)
neg.qf <- neg.vecs%*%t(neg.vecs)
pos.qf <- pos.vecs%*%t(pos.vecs)
op <- par(mfrow=c(1,2))
image(pos.qf);image(neg.qf)
par(op)
## round(pos.qf/neg.qf,2)
phi.vecs1 <- sapply(1:I,function(i)e(1+(i-1)*(I+1)))
phi.vecs2 <- e(which(rowSums(phi.vecs1)==0))
## (t(phi.vecs1)%*%pos.qf%*%phi.vecs1) / (t(phi.vecs1)%*%neg.qf%*%phi.vecs1)
## (t(phi.vecs2)%*%pos.qf%*%phi.vecs2) / (t(phi.vecs2)%*%neg.qf%*%phi.vecs2)
## (t(phi.vecs1)%*%pos.qf%*%phi.vecs2) / (t(phi.vecs1)%*%neg.qf%*%phi.vecs2)
a <- phi.vecs2[,3]
b <- phi.vecs1[,1]
## phi.vecs1 and phi.vecs2 orthogonal wrt pos.qf+neg.qf (proportional
## to Q^2)
## t(a+b)%*%(pos.qf + neg.qf)%*%(a+b)
## t(a)%*%(pos.qf + neg.qf)%*%b
## t(a)%*%(pos.qf + neg.qf)%*%a + t(b)%*%(pos.qf + neg.qf)%*%b
## t(a+b)%*%(pos.qf - neg.qf)%*%(a+b)
## t(a)%*%(pos.qf - neg.qf)%*%a + t(b)%*%(pos.qf - neg.qf)%*%b
## (t(b)%*%(pos.qf - neg.qf)%*%b) / t(a)%*%(pos.qf - neg.qf)%*%a
## t(phi.vecs1)%*%Q%*%phi.vecs1
## t(phi.vecs2)%*%Q%*%phi.vecs2
unique(round(as.numeric(t(phi.vecs2)%*%Q%*%phi.vecs2) /as.numeric(t(phi.vecs1)%*%Q%*%phi.vecs1),4))

op <- par(mfrow=c(1,2))
image(t(phi.vecs1)%*%Q%*%phi.vecs1)
image(t(phi.vecs2)%*%Q%*%phi.vecs2)
par(op)

I <- 10
Is <- 3:20; names(Is) <- Is
e <- Vectorize(function(j,I)unlist(within(list(v=rep(0,I^2)),v[j] <- 1)))
ratios <- lapply(Is, function(I) {
    print(I)
    Q <- M.obj(I,sym=TRUE)   
    phi.vecs1 <- sapply(1:I,function(i)e(1+(i-1)*(I+1),I))
    phi.vecs2 <- e(which(rowSums(phi.vecs1)==0),I)
    sort(unique(round(as.numeric(t(phi.vecs2)%*%Q%*%phi.vecs2) /as.numeric(t(phi.vecs1)%*%Q%*%phi.vecs1),4)))
})


## 6aa
Is <- 3:20; names(Is) <- Is
ratios <- sapply(Is, function(I) {
    Q <- M.obj(I,sym=TRUE)
    pos.vecs <- with(eigen(Q),vectors[,values>1e-7])
    neg.vecs <- with(eigen(Q),vectors[,values< -1e-7])
    neg.qf <- neg.vecs%*%t(neg.vecs)
    pos.qf <- pos.vecs%*%t(pos.vecs)
    (neg.qf/pos.qf)[1]
})
plot(Is,ratios)
ratios
## can't identify sequence...  3: 17-12*sqrt(2), 4: 0, 5: 49-20*sqrt(6), 10: 1/9, 20: 1/4. 2*I^2-1 - (4I)*sqrt(...). ... = ((2*I^2-1 - x)/(4*I))^2
with(list(x=1/4,I=20),((2*I^2-1 - x)/(4*I))^2)
## seems conjugate is inverse. f below refers to fundamental unit of
## the field a+b\sqrt(D).
3: 17+12*sqrt(2) = f(2)^4
4: Inf
5: 49+20*sqrt(6) = f(6)^2
6: 17+12*sqrt(2) = 17+6*sqrt(8) = f(2)^4 (= case 3)
8: 7+4*sqrt(3) = 7+2*sqrt(12) = f(3)^2
10: 9
12: 7/2+3/2*sqrt(5) = f(5)^4
20: 4

x <- ratios['13']
for(x in ratios[11:length(ratios)]) {
    print(x); print('---')
d <- 0
while(sqrt(d)<x) {
    a <- 0
    while(a<x) {
        b <- 0
        while(b<x) {
            if(abs(a+b*sqrt(d)-x)<1e-4)cat(a,' ',b,' ',d,'\n')
            b <- b+1
        }
        a <- a+1
    }
    d <- d+1
}
}

dd

## 6b. represent phi as PF (arbitrary m,n)
## data/phi => (P,F)
source('misc.R')
I <- 3
m <- sample(1:10,I,replace=TRUE)
n <- sample(1:10,I,replace=TRUE)
M <- sum(m); N <- sum(n)
cumsum.mat <- matrix(0,nrow=M+N,ncol=M+N)
cumsum.mat[row(cumsum.mat)>=col(cumsum.mat)] <- 1
x <- lapply(m,function(m.i)runif(m.i))
y <- lapply(n,function(n.i)runif(n.i))
x.rank <- split(rank(c(unlist(x),unlist(y)))[1:M],rep(1:I,m))
y.rank <- split(rank(c(unlist(y),unlist(x)))[1:N],rep(1:I,n))
F <- simplify2array(sapply(x.rank,function(idx)within(list(v=rep(0,M+N)),v[idx] <- 1)))
## F <- cumsum.mat%*%F%*%(diag(1/m))
## sum(abs(cumsum.mat%*%F - apply(F,2,cumsum)))
P <- t(simplify2array(sapply(y.rank,function(idx)within(list(v=rep(0,M+N)),v[idx] <- 1))))
## P <- diag(1/n)%*%t(simplify2array(P))
phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord(x.i,y.i))))
phi.try <- t(diag(1/n)%*%P%*%cumsum.mat%*%F%*%diag(1/m))
sum(abs(phi-phi.try))

## (P,F) => data
source('misc.R')
I <- 3
m <- sample(1:10,I,replace=TRUE)
n <- sample(1:10,I,replace=TRUE)
M <- sum(m); N <- sum(n)
cumsum.mat <- matrix(0,nrow=M+N,ncol=M+N)
cumsum.mat[row(cumsum.mat)>=col(cumsum.mat)] <- 1
F <- simplify2array(sapply(1:I,function(i)within(list(v=rep(0,M+N)),v[sample(1:(M+N),m[i])] <- 1)))
P <- simplify2array(sapply(1:I,function(i)within(list(v=rep(0,M+N)),v[sample(1:(M+N),n[i])] <- 1)))
phi.PF <- t(  diag(1/n)%*%t(P)%*%cumsum.mat%*%F%*%diag(1/m)  )
phi.PF <- diag(1/m)%*%t(F)%*%t(cumsum.mat)%*%P%*%diag(1/n)
x <- apply(F,2,function(v)which(v==1))
y <- apply(P,2,function(v)which(v==1))
phi.data <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord(x.i,y.i))))
sum(abs(phi.data-phi.PF))

## 6c. vectorized phi using P,F decomposition
I <- 3
m <- sample(1:10,I,replace=TRUE)
n <- sample(1:10,I,replace=TRUE)
M <- sum(m); N <- sum(n)
cumsum.mat <- matrix(0,nrow=M+N,ncol=M+N)
cumsum.mat[row(cumsum.mat)>=col(cumsum.mat)] <- 1
F <- simplify2array(sapply(1:I,function(i)within(list(v=rep(0,M+N)),v[sample(1:(M+N),m[i])] <- 1)))
P <- simplify2array(sapply(1:I,function(i)within(list(v=rep(0,M+N)),v[sample(1:(M+N),n[i])] <- 1)))
phi <- diag(1/m)%*%t(F)%*%t(cumsum.mat)%*%P%*%diag(1/n)
## kronecker(diag(I),diag(1/m)%*%t(F)) %*% matrix((t(cumsum.mat)%*%P%*%diag(1/n)),ncol=1) # vectorizes phi but column major
F.diag <- do.call(rbind,apply(diag(1/m)%*%t(F),1,function(v)kronecker(diag(I),matrix(v,nrow=1)),simplify=FALSE))
phi.vec <- F.diag %*% matrix((t(cumsum.mat)%*%P%*%diag(1/n)),ncol=1)

## 6d. compute main obj using (P,F), assuming m,n fixed again
source('misc.R')
I <- 3
m <- rep(5,I); n <- rep(3,I)
M <- sum(m); N <- sum(n)
## parts <- replicate(1e2,{
## x <- replicate(I,runif(m),simplify=FALSE)
## y <- replicate(I,rnorm(n),simplify=FALSE)
## phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
## psi <- phi*m*n
## phi.vec <- matrix(t(phi),ncol=1)
cumsum.mat <- matrix(0,nrow=M+N,ncol=M+N)
cumsum.mat[row(cumsum.mat)>=col(cumsum.mat)] <- 1
F <- simplify2array(sapply(1:I,function(i)within(list(v=rep(0,M+N)),v[sample(1:(M+N),m[i])] <- 1)))
P <- simplify2array(sapply(1:I,function(i)within(list(v=rep(0,M+N)),v[sample(1:(M+N),n[i])] <- 1)))
phi <- diag(1/m)%*%t(F)%*%t(cumsum.mat)%*%P%*%diag(1/n)
## kronecker(diag(I),diag(1/m)%*%t(F)) %*% matrix((t(cumsum.mat)%*%P%*%diag(1/n)),ncol=1) # vectorizes phi but column major
F.diag <- do.call(rbind,apply(diag(1/m)%*%t(F),1,function(v)kronecker(diag(I),matrix(v,nrow=1)),simplify=FALSE))
phi.vec <- F.diag %*% matrix((t(cumsum.mat)%*%P%*%diag(1/n)),ncol=1)
theta.hat <- sum(phi)/I^2
V.bar.x <- rowSums(phi)/I#sapply(x,function(x.i)concord(x.i,unlist(y),mean=FALSE))/N/m
V.bar.y <- colSums(phi)/I#sapply(y,function(y.i)concord(unlist(x),y.i,mean=FALSE))/M
## x <- apply(F,2,function(v)which(v==1))
## y <- apply(P,2,function(v)which(v==1))
## x <- split(x,rep(1:I,m)); y <- split(y,rep(1:I,n))
## rowSums(phi)/I
## sapply(x,function(x.i)concord(x.i,unlist(y),mean=FALSE))/N/m
sum((V.bar.x+V.bar.y-2*theta.hat)*(V.bar.x+V.bar.y-diag(phi)))
t(phi.vec)%*%M.obj(I)%*%phi.vec
t(phi.vec)%*%M.obj(I,symm=TRUE)%*%phi.vec
matrix(t(diag(1/n)%*%t(P)%*%cumsum.mat),nrow=1)%*%t(F.diag)%*%M.obj(I,symm=TRUE)%*%F.diag %*% matrix((t(cumsum.mat)%*%P%*%diag(1/n)),ncol=1)
rowMeans(phi)
phi
phi-rowMeans(phi) # M*phi == M*(phi-rowmeans)
t(matrix(phi*0+rowMeans(phi),ncol=1))%*%M.obj(I,symm=TRUE)%*%matrix(phi*0+rowMeans(phi),ncol=1)

## 6e nullspace of quad form
source('misc.R')
I <- 3
pos.idx <- c(1:(I-1),(I^2-I+2):(I^2))
## eigen(M.obj(I,symm=TRUE))$vec[,pos.idx]
M0 <- M.obj(I,sym=TRUE)
v <- c(rep(1,I),rep(0,I^2-I))
v <- rep(runif(I),each=I)
t(v)%*%M0%*%v
## root.M <- with(eigen(M0), vectors %*% diag(sqrt(as.complex(values))) %*% t(vectors))
## sum(abs(root.M%*%root.M-M0))
## t(matrix(phi*0+rowMeans(phi),ncol=1))%*%root.M%*%root.M%*%matrix(phi*0+rowMeans(phi),ncol=1)
## (root.M%*%matrix(phi*0+rowMeans(phi),ncol=1))^2
## not useful, need hermitian sqrt
eigen(M0)$val[pos.idx]
eigen(M0)$vec[,pos.idx]
## for wolfram:
paste0("{{",paste(apply(round(M0*I^3*2,7),1,paste,collapse=','),collapse="},{"),"}}")


## M2 <- function(I,partition=FALSE,symmetric=FALSE) {
##     E <- diag(I); ones <- matrix(1,nrow=I,ncol=1)
##     m.a <- kronecker(E,ones%*%t(ones))
##     m.b <- kronecker(ones%*%t(ones),E)
##     m.c <- kronecker(ones,kronecker(E,t(ones)))
##     m.d <- Reduce(`+`,lapply(1:I,function(j)kronecker(matrix(E[j,],nrow=1),kronecker(kronecker(E[,j],E[,j]),t(ones)))))
##     m.e <- Reduce(`+`,lapply(1:I,function(j)kronecker(E[,j],kronecker(t(ones),E[,j]%*%t(E[,j])))))
##     m.f <- kronecker(ones%*%t(ones),ones%*%t(ones))
##     m.g <- kronecker(matrix(E,nrow=1),kronecker(ones,ones))
##     m.lst <- list(a=m.a/I^2,b=m.b/I^2,c=m.c*2/I^2,d=-1/I*m.d,e=-1/I*m.e,f=-4/I^3*m.f,g=2/I^2*m.g)
##     ans <- Reduce(`+`,m.lst)
##     if(symmetric) ans <- (ans+t(ans))/2
##     ## sapply(m.lst,function(m)t(v)%*%m%*%v)
##     if(partition) {
##         ans <- partition(ans)
##     }
##     ans
## }


## 6f powers of quad form
source('misc.R')
I <- 5
mm0 <- with(list(m0=M.obj(I,sym=TRUE)*I^2),m0%*%m0)
mm0 <- round(mm0,7)
ans <- list()
rs <- which(mm0[1,]!=0)
cs <- which(mm0[,1]!=0)
for(r in 2:length(rs))
    for(c in 2:length(cs))
        ans <- c(ans,list(mm0[(rs[r-1]+1):(rs[r]-1),(cs[c-1]+1):(cs[c]-1)]))
lapply(ans[1:(I-1)],function(m)(m==-I)+0)

## normalized even powers are orthog projections
## normalization is 2*I^2/(I-2)
I <- 3
m0 <- M.obj(I,sym=TRUE)*I^2
mm0 <- 2*m0%*%m0 / I^2/(I-2)
norm(mm0%*%mm0-mm0)

m1 <- M.obj(I,sym=TRUE)*I*sqrt(2/(I-2))

m0 <- M.obj(I,sym=TRUE)
mm0 <- m0%*%m0
k <- 5; norm(mat.pow(mm0,k) - (2*I^2/(I-2))^(1-k)*mm0)

## M is eigenvectors of M%*%M with eigenvalue 
I <- 3
m0 <- M.obj(I,sym=TRUE)
mm0 <- m0%*%m0 
(mm0%*%m0) / m0
## 9/2,16,75/2,72 = 1/2*I^2*(I-2)
## (m0%*%m0 * I^4) %*% m0*I^2 = 1/2*I^2*(I-2) * m0*I^2
norm((mm0 * I^4) %*% m0*I^2  - 1/2*I^2*(I-2) * m0*I^2)

## trace(M%*%M)=(I-2)(I-1)/I^2
sapply(2:9,function(I)sum(diag(mat.pow(M.obj(I,sym=TRUE),2)))- (I-2)*(I-1)/I^2)
## trace((M%*%M)^(2k))=2(I-2)(2I^2/(I-1))^-k
I <- 5
sapply(2:9,function(k)sum(diag(mat.pow(M.obj(I,sym=TRUE),2*k))) - 2*(I-1)*(2*I^2/(I-2))^(-k))


tk <- function(I,k) if(k%%2==1) 0 else 2*(I-1)*(2*I^2/(I-2))^(-k)



## 6g delta description of quad form
source('misc.R')
Q <- function(p,q,r,s,I=3)
    1/I^2*sum(apply(combn(1:4,2),2,function(idx)abs(diff(c(p,q,r,s)[idx]))<1e-7)) -1/(2*I)*sum(apply(combn(1:4,3),2,function(idx)sum(abs(diff(c(p,q,r,s)[idx])))<1e-7)) - 4/I^3

I <- 3
m0 <- M.obj(I,sym=TRUE)
p <- 1;q <- 1;r <- 1;s <- 1
m0[(p-1)*I+r,(q-1)*I+s]
Q(p,q,r,s,I=I)
sum(apply(expand.grid(rep(list(1:I),4)),1,function(pqrs)m0[(pqrs[1]-1)*I+pqrs[3],(pqrs[2]-1)*I+pqrs[4]]!=Q(pqrs[1],pqrs[2],pqrs[3],pqrs[4],I=I)))




Q1 <- function(p,q,r,s,I=3) sum(apply(combn(1:4,2),2,function(idx)abs(diff(c(p,q,r,s)[idx]))<1e-7))
Q1.try <- function(p,q,r,s,I=3) 2*I*((p==q)+(q==r)+(r==s)+(p==s))+5*I*((p==r)+(q==s))+I^2*(q==s)*(p==r)+I+16
Q1.try <- function(p,q,r,s,I=3) 2*I*Q1(p,q,r,s)+3*I*((p==r)+(q==s))+I^2*(q==s)*(p==r)+I+16
p <- sample(1:I,1);q <- sample(1:I,1);r <- sample(1:I,1);s <- sample(1:I,1)
sum(apply(expand.grid(u=1:3,v=1:3),1,function(uv)Q1(uv['u'],q,uv['v'],s)*Q1(p,uv['u'],r,uv['v'])))
Q1.try(p,q,r,s)


Q2 <- function(p,q,r,s,I=3)sum(apply(combn(1:4,3),2,function(idx)sum(abs(diff(c(p,q,r,s)[idx])))<1e-7)) 
Q2.try <- function(p,q,r,s,I=3) 2*Q2(p,q,r,s,I=I)+(p==q)+(q==r)+(r==s)+(p==s)+2*(q==s)*(p==r)+2*I*(p==q)*(q==r)*(r==s)
Q2.try <- function(p,q,r,s,I=3) 2*Q2(p,q,r,s,I)+Q1(p,q,r,s,I)-((q==s)-(p==r))^2+2*I*(p==q)*(q==r)*(r==s)
p <- sample(1:I,1);q <- sample(1:I,1);r <- sample(1:I,1);s <- sample(1:I,1)
sum(apply(expand.grid(u=1:3,v=1:3),1,function(uv)Q2(uv['u'],q,uv['v'],s)*Q2(p,uv['u'],r,uv['v'])))
Q2.try(p,q,r,s)

Q12.try <- function(p,q,r,s,I=3) 2*Q1(p,q,r,s)+4*(p==r)+2*I*((p==q)*(q==r)+(p==r)*(r==s)+(p==r)*(q==s))+2
p <- sample(1:I,1);q <- sample(1:I,1);r <- sample(1:I,1);s <- sample(1:I,1)
sum(apply(expand.grid(u=1:3,v=1:3),1,function(uv)Q1(uv['u'],q,uv['v'],s)*Q2(p,uv['u'],r,uv['v'])))
Q12.try(p,q,r,s)


Q21.try <- function(p,q,r,s,I=3) 2*Q1(p,q,r,s)+4*(q==s)+2*I*((p==q)*(q==s)+(q==r)*(q==s)+(p==r)*(q==s))+2
p <- sample(1:I,1);q <- sample(1:I,1);r <- sample(1:I,1);s <- sample(1:I,1)
sum(apply(expand.grid(u=1:3,v=1:3),1,function(uv)Q2(uv['u'],q,uv['v'],s)*Q1(p,uv['u'],r,uv['v'])))
Q21.try(p,q,r,s)


Qa <- function(p,q,r,s,I=3) 1/I^2*Q1(p,q,r,s,I)-1/2/I*Q2(p,q,r,s,I)
p <- sample(1:I,1);q <- sample(1:I,1);r <- sample(1:I,1);s <- sample(1:I,1)
sum(apply(expand.grid(u=1:3,v=1:3),1,function(uv)Qa(uv['u'],q,uv['v'],s)*Qa(p,uv['u'],r,uv['v'])))
1/I^4*Q1.try(p,q,r,s)+1/4/I^2*Q2.try(p,q,r,s)-1/2/I^3*(Q12.try(p,q,r,s)+Q21.try(p,q,r,s))
1/I^4*Q1.try(p,q,r,s)+1/4/I^2*Q2.try(p,q,r,s)-2/I^3*(Q1(p,q,r,s,I)+1/2*I*Q2(p,q,r,s,I)+(p==r)+(q==s)+I*(p==r)*(q==s)+1)
with(list(x1=p==q,x2=p==s,x3=r==q,x4=r==s,x5=q==s,x6=p==r),1/4/I^4*(64-4*I+I^2*(x2+x3+x4 +x1*(1+2*x3*(-1+I*x4)-2*x5)) +4*I*x5+I^2*x5-2*I^2*x3*x5-I^2*x5^2+4*I*x6+I^2*x6-2*I^2*x4*x6-2*I^2*x5*x6-I^2*x6^2))
with(list(x1=p==q,x2=p==s,x3=r==q,x4=r==s,x5=q==s,x6=p==r),1/4/I^4*(64-4*I+I^2*Q1(p,q,r,s,I) +I^2*x1*(2*x3*(-1+I*x4)-2*x5) +4*I*x5-2*I^2*x3*x5-I^2*x5^2+4*I*x6-2*I^2*(x4*x6+x5*x6)-I^2*x6^2))
1/4/I^4*(64-4*I+I^2*Q1(p,q,r,s,I)-I^2*(2*Q2(p,q,r,s,I)+((p==r)+(q==s))^2-2*I*(p==q)*(q==r)*(r==s))+4*I*((q==s)+(p==r)))

p <- sample(1:I,1);q <- sample(1:I,1);r <- sample(1:I,1);s <- sample(1:I,1)
sum(apply(expand.grid(u=1:3,v=1:3),1,function(uv)Q(uv['u'],q,uv['v'],s))+4/I^3)
4/I


## Q <- function(p,q,r,s,I=3)
##     1/I^2*sum(apply(combn(1:4,2),2,function(idx)abs(diff(c(p,q,r,s)[idx]))<1e-7)) -1/(2*I)*sum(apply(combn(1:4,3),2,function(idx)sum(abs(diff(c(p,q,r,s)[idx])))<1e-7)) - 4/I^3
## Q <- function(p,q,r,s,I=3)Qa(p,q,r,s,I)-4/I^3
p <- sample(1:I,1);q <- sample(1:I,1);r <- sample(1:I,1);s <- sample(1:I,1)
sum(apply(expand.grid(u=1:3,v=1:3),1,function(uv)Q(uv['u'],q,uv['v'],s)*Q(p,uv['u'],r,uv['v'])))
1/4/I^4*(64-4*I+I^2*Q1(p,q,r,s,I)-I^2*(2*Q2(p,q,r,s,I)+((p==r)+(q==s))^2-2*I*(p==q)*(q==r)*(r==s))+4*I*((q==s)+(p==r))) - 32/I^4 + 16/I^4
1/4/I^4*(-4*I+I^2*Q1(p,q,r,s,I)-I^2*(2*Q2(p,q,r,s,I)+((p==r)+(q==s))^2-2*I*(p==q)*(q==r)*(r==s))+4*I*((q==s)+(p==r)))


I <- 3
m0 <- M.obj(I,sym=TRUE)
mm0 <- m0%*%m0
Q.sqr <- function(p,q,r,s,I=3)  1/4/I^4*(-4*I+I^2*Q1(p,q,r,s,I)-I^2*(2*Q2(p,q,r,s,I)+((p==r)+(q==s))^2-2*I*(p==q)*(q==r)*(r==s))+4*I*((q==s)+(p==r)))
p <- sample(1:I,1);q <- sample(1:I,1);r <- sample(1:I,1);s <- sample(1:I,1)
mm0[(p-1)*I+r,(q-1)*I+s]
Q.sqr(p,q,r,s)
sum(apply(expand.grid(rep(list(1:I),4)),1,function(pqrs)abs(mm0[(pqrs[1]-1)*I+pqrs[3],(pqrs[2]-1)*I+pqrs[4]]-Q.sqr(pqrs[1],pqrs[2],pqrs[3],pqrs[4],I=I))))

sapply(2:10,function(I) {
m0 <- M.obj(I,sym=TRUE)
mm0 <- m0%*%m0
mmm0 <- mm0%*%m0
alpha <- (I-2)/2/I^2
norm(mmm0-alpha*m0)
})


## 6g nullspace of Q^2

source('misc.R')
I <- 5
m0 <- M.obj(I,sym=TRUE)
mm0 <- m0%*%m0 * (2*I^2)/(I-2)
evs <- sapply(0:(I-2), function(k)within(list(v=rep(0,I^2)),{
    v[k+2] <- -1
    v[I+1+k*I] <- 1
    }))
evs <- simplify2array(evs)
norm(mm0%*%evs)
norm(m0%*%evs)
ev <- within(list(v=rep(0,I^2)),v[1+(I+1)*(0:(I-1))] <- 1)
ev <- simplify2array(ev)
norm(m0%*%ev)

evs <- lapply(0:(I-2),function(k)t(simplify2array(sapply(1:(I-2),function(j) {
    v <- rep(0,I^2)
    r <- 1+(I-1)*k+j
    ## if(r>=5)browser()
    v[I+r+floor((r-2)/I)+1] <- 1
    v
    }))))
evs <- do.call(rbind,evs)

## I <- 7
## ones <- matrix(1,nrow=I-2)
## a <- 2/(I-2)*kronecker(ones,t(ones))
## a <- a+diag(-1,I-2)
## ## a <- cbind(a,0,diag(I-2))
## dup.j <- as.matrix(bdiag(diag(2),matrix(c(-3/5,1),ncol=2),diag(2)))

evs <- lapply(0:(I-2),function(k)t(simplify2array(sapply(1:(I-2),function(j) {
    v <- rep(0,I^2)
    p <- k+1
    r <- 1+(I-1)*k+j
    ## if(r==8)browser()
    v[I+r+floor((r-2)/I)+1] <- 1
    idx <- 1:length(v)
    v[(2<=idx)&(idx<=min(p,j))] <- 2/(I-2)
    v[1+min(p,j)] <- v[1+max(p,j)+(p<=j)] <- 2/(I-2)-1
    v[1+min(p,j)<idx & idx<1+max(p,j)+(p<=j)] <- 2/(I-2)
    v[1+max(p,j)+(p<=j) < idx & idx<=I] <- 2/(I-2)
    v
    }))))
evs <- do.call(rbind,evs)
norm(m0%*%t(evs))

##6e permutation symmetry of Q
idx <- function(p,r,q=NULL,s=NULL) {
    if(is.null(q)&is.null(s)) return(as.numeric(I*(p-1)+r))
    return(as.numeric(c(I*(p-1)+r,I*(q-1)+s)))
}
## idx <- function(pqrs) {
##     if(length(pqrs)==2) return(as.numeric(I*(pqrs['p']-1)+pqrs['r']))
##     return(as.numeric(c(I*(pqrs['p']-1)+pqrs['r'],I*(pqrs['q']-1)+pqrs['s'])))
## }
pi <- function(j)3-j
pqrs <- expand.grid(structure(rep(list(1:2),4),names=c('p','q','r','s')))
vv <- cbind(pqrs,a=(apply(pqrs,1,function(x)with(as.list(x),Q[matrix(idx(p,q,r,s),1)]*phi.vec[matrix(idx(p,r),1)]*phi.vec[matrix(idx(q,s),1)]))),
b=(apply(pqrs,1,function(x)with(as.list(x),Q[matrix(idx(p,q,r,s),1)]*phi.vec[matrix(idx(pi(p),pi(r)),1)]*phi.vec[matrix(idx(pi(q),pi(s)),1)])))
)
vv
vv[,'a']-rev(vv[,'b'])
vv2 <- cbind(pqrs,a=(apply(pqrs,1,function(x)with(as.list(x),Q[matrix(idx(p,q,r,s),1)]*phi.vec[matrix(idx(p,r),1)]*phi.vec[matrix(idx(q,s),1)]))),
b=(apply(pqrs,1,function(x)with(as.list(x),Q[matrix(idx(pi(p),pi(q),pi(r),pi(s)),1)]*phi.vec[matrix(idx(pi(p),pi(r)),1)]*phi.vec[matrix(idx(pi(q),pi(s)),1)])))
)
vv2

with(as.list(pqrs[1,]),
c(Q[matrix(idx(p,q,r,s),1)]*phi.vec[matrix(idx(p,r),1)]*phi.vec[matrix(idx(q,s),1)], Q[matrix(idx(pi(p),pi(q),pi(r),pi(s)),1)]*phi.vec[matrix(idx(pi(p),pi(r)),1)]*phi.vec[matrix(idx(pi(q),pi(s)),1)])
)

i1 <- sample(1:I,1);i2 <- sample((1:I)[-i1],1)
pqrs <- expand.grid(structure(rep(list(c(i1,i2)),4),names=c('p','q','r','s')))
pi <- function(j)c(i2,i1)[which(c(i1,i2)==j)]
sapply(1:I,function(i)with(as.list(pqrs[i,]),
c(Q[matrix(idx(p,q,r,s),1)], Q[matrix(idx(pi(p),pi(q),pi(r),pi(s)),1)])
))

sapply(1:I,function(j)Q[matrix(idx(j,j,j,j),1)])
pqrs <- expand.grid(structure(rep(list(1:2),4),names=c('p','q','r','s')))

sapply(1:nrow(pqrs),function(i)sd(replicate(1e3,{
    pi <- sample(1:I)
    with(as.list(pqrs[i,]),Q[matrix(idx(pi[p],pi[q],pi[r],pi[s]),1)])
})))

sum(abs(replicate(1e3,{
pqrs <- structure(sample(1:I,4,replace=TRUE),names=c('p','q','r','s'))
sd(replicate(1e2,{
    pi <- sample(1:I)
    with(as.list(pqrs),Q[matrix(idx(pi[p],pi[q],pi[r],pi[s]),1)])
}))
})))



## 8 general m,n
## checking formulas
I <- 3
m <- sample(1:10,I,replace=TRUE); n <- sample(1:10,I,replace=TRUE)
## m  <- n <- m*0+5
x <- lapply(m,function(m.i)rnorm(m.i))
y <- lapply(n,function(n.i)rnorm(n.i))
M <- sum(m); N <- sum(n)
theta.hat <- concord(unlist(x),unlist(y))
theta.del <- sapply(1:I,function(i)concord(unlist(x[-i]),unlist(y[-i])))
pseudo <- I*theta.hat - (I-1)*theta.del
theta.del.s <- sapply(1:I,function(i)concord.sum(unlist(x[-i]),unlist(y[-i])))
pseudo <- I*theta.hat - (I-1)*theta.del
pseudo.s <- 
theta.jk <- mean(pseudo)
psi <- outer(x,y,Vectorize(concord.sum))
## ## 8a
## (I-2)*(sum(psi)-sum(diag(psi))) + (I-1)*sum(diag(psi))
## sum(theta.del.s)
## pseudo
## I*theta.hat - (I-1)*theta.del
## I/M/N*sum(psi) - (I-1)/(M-m)/(N-n)*(sapply(1:I,function(i)sum(psi[-i,-i])))
## mean(pseudo)
## I/M/N*sum(psi) - (I-1)/I/(M-m)/(N-n)*((I-2)*(sum(psi)-sum(diag(psi)))+(I-1)*sum(diag(psi)))
## theta.jk.try <- sum(psi)*I/M/N
## for(j in 1:I) {
##     theta.jk.try <- theta.jk.try-(I-1)/I*psi[j,j]*sum((1/(M-m)/(N-n))[-j])
##     for(k in 1:I)
##         if(j!=k) {
##             theta.jk.try <- theta.jk.try-(I-1)/I*psi[j,k]*sum((1/(M-m)/(N-n))[-c(j,k)])
##         }
## }
## theta.jk.try
## theta.jk
## ## 8b
## theta.jk.try <- 0
## for(j in 1:I) {
##     theta.jk.try <- theta.jk.try+psi[j,j]*(I/M/N-(I-1)/I*sum((1/(M-m)/(N-n))[-j]))
##     for(k in 1:I)
##         if(j!=k) {
##             theta.jk.try <- theta.jk.try+psi[j,k]*(I/M/N - (I-1)/I*sum((1/(M-m)/(N-n))[-c(j,k)]))
##         }
## }
## theta.jk.try
## theta.jk
theta.jk.try <- 0
for(j in 1:I) {
    ## browser()
    theta.jk.try <- theta.jk.try+1/M/N*psi[j,j]*(I-(I-1)*M*N/I*sum((1/(M-m)/(N-n))[-j]))
    for(k in 1:I)
        if(j!=k) {
            theta.jk.try <- theta.jk.try+1/M/N*psi[j,k]*(I - (I-1)*M*N/I*sum((1/(M-m)/(N-n))[-c(j,k)]))
        }
}
theta.jk.try
theta.jk




## 10 [out of order] block structure of quad form
source('misc.R')
I <- 5#sample(3:10,1)
diag.idx <- ((1:I^2)-1)%%(I+1)==0
parts <- M.obj(I,symm=TRUE,part=TRUE)
## lapply(parts,function(m)round(eigen(m)$vec,3))
norm(apply(t(simplify2array(lapply(parts,function(m)eigen(m)$val)[!diag.idx])),2,diff))##offdiag blocks all have the same eigenvals
norm(apply(t(simplify2array(lapply(parts,function(m)eigen(m)$val)[diag.idx])),2,diff)) ##as do blocks on the diagonal
## t(simplify2array(lapply(parts,function(m)round(eigen(m)$val,3))[diag.idx]))
## off diagonal blocks have I-2 evals all 1/I^2, and 2 evals -(I-2)/2/I^2:
norm(simplify2array(lapply(parts,function(m)eigen(m)$val)[!diag.idx]) - c(rep(1/I^2,I-2),rep(-(I-2)/2/I^2,2)))
## diagonal blocks also have I-2 evals all 1/I^2, and then 2 evals (I-2)/2/I^2*(-1\pm  I^(1/2)):
norm(simplify2array(lapply(parts,function(m)sort(eigen(m)$val))[diag.idx]) -
     sort(c(rep(1/I^2,I-2),(I-2)/2/I^2*(-1-sqrt(I)),(I-2)/2/I^2*(-1+sqrt(I)))))


## 9 phi^TQphi growth rate

## ## 9a
## ## M() too slow, quad form computing directly
## source('misc.R')
## require(parallel)
## m <- n <- 5
## Is <- round(seq(2,70,len=50))
## ## by.I <- sapply(Is,function(I) {
## by.I <- mclapply(Is,mc.cores=detectCores()-4,FUN=function(I) {
##     print(I)
##     qfs <- replicate(1e2, {
##         M <- I*m; N <- I*n
##         x <- replicate(I,runif(m),simplify=FALSE)
##         y <- replicate(I,rnorm(n),simplify=FALSE)
##         phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
##         ## psi <- phi*m*n
##         phi.vec <- matrix(t(phi),ncol=1)
##         theta.hat <- concord.mean(unlist(x),unlist(y))
##         V.x <- sapply(x,function(x.i)concord.sum(x.i,unlist(y)))/N
##         V.y <- sapply(y,function(y.i)concord.sum(unlist(x),y.i))/M
##         V.bar.x <- V.x/m; V.bar.y <- V.y/n
##         ## c(a=sum(V.bar.x^2),b=sum(V.bar.y^2),c=2*sum(V.bar.x*V.bar.y),d=-sum(diag(phi)*V.bar.x),e=-sum(diag(phi)*V.bar.y),f=-2*theta.hat*sum(V.bar.x+V.bar.y),g=2*theta.hat*sum(diag(phi)))
##         ## t(phi.vec)%*%M.obj(I)%*%phi.vec
##         sum((V.bar.x+V.bar.y-2*theta.hat)*(V.bar.x+V.bar.y-diag(phi)))
##     })
##     ## c(max=max(abs(qfs)),median(qfs)
##     qfs
## })
## by.I <- simplify2array(by.I)

## maxes <- apply(by.I,2,function(x)max(abs(x)))
## plot(Is,maxes)

## plot(Is,maxes/Is^(1/2))
## ## save.image('081221.RData')
## ## load('081221.RData')
## ## maxes are O(I^(1/2)), as expected


## 9b growth rate of convex part of Q
source('misc.R')
require(parallel)
I <- 2
Is <- round(seq(2,30,len=30))
## by.I <- sapply(Is,function(I) {
by.I <- mclapply(Is,mc.cores=detectCores()-4,FUN=function(I) {
    print(I)
    Q <- M.obj(I,symm=TRUE) # slow here
    vals <- eigen(Q)$val; vecs <- eigen(Q)$vec
    neg.idx <- vals< -1e-8;     pos.idx <- vals>1e-8
    ## pos.vecs <- vecs[,pos.idx]*vals[pos.idx]
    ## neg.vecs <- vecs[,neg.idx]*vals[neg.idx]
    m <- 5; n <- 5
    qfs <- replicate(20, {
        x <- replicate(I,runif(m),simplify=FALSE)
        y <- replicate(I,rnorm(n),simplify=FALSE)
        phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
        phi.vec <- matrix(t(phi),ncol=1)
        ## ## phi.vec <- matrix(rbinom(I^2,1,.5),ncol=1)
        ## ## theta.hat <- concord.mean(unlist(x),unlist(y))
        ## phi.vec <- matrix(runif(I^2),ncol=1)
        ## F <- t(replicate(I,sort(runif(I*n))))
        ## wts <- sample(rep(1:I,n))
        ## phi <- sapply(1:I,function(i)F%*%(wts==i)/n)
        ## ## F <- t(replicate(I,sample(0:1,I*n,replace=TRUE)))
        ## ## F <- t(apply(F,1,cumsum))
        ## ## F <- matrix(rnorm(I^2),nrow=I,ncol=I)
        ## M.cumsum <- with(list(m=matrix(nrow=I,ncol=I)),(row(m)<=col(m))+0)
        ## M.cumsum <- diag(I)
        ## ## F <- F%*%M.cumsum
        ## phi <- matrix(rnorm(I^2),nrow=I,ncol=I)%*%M.cumsum%*%matrix(rexp(I^2,4),nrow=I,ncol=I)
        ## phi <- matrix(rnorm(I^2,0),nrow=I,ncol=I)%*%matrix(runif(I^2),nrow=I,ncol=I)
        ## phi.vec <- matrix(runif(I^2,-5,5),ncol=1)
        ## phi.vec <- matrix(phi,ncol=1)
        c(pos=sum(vals[pos.idx]*(t(vecs[,pos.idx])%*%phi.vec)^2),neg=sum(vals[neg.idx]*(t(vecs[,neg.idx])%*%phi.vec)^2),total=t(phi.vec)%*%Q%*%phi.vec)
    })
    ## c(max=max(abs(qfs)),median(qfs)
    qfs
})
by.I <- simplify2array(by.I)
maxes.pos <- apply(by.I['pos',,],2,function(x)max(abs(x)))
maxes.neg <- apply(by.I['neg',,],2,function(x)max(abs(x)))
maxes.total <- apply(by.I['neg',,]+by.I['pos',,],2,function(x)max(abs(x)))
plot(Is,maxes.pos,ylim=range(c(maxes.pos,maxes.neg,maxes.total)))
points(Is,maxes.neg,col=2)
points(Is,maxes.total,col=3)

plot(Is,maxes.total/sqrt(Is))
plot(Is,maxes.pos/Is)
## pos and neg part growing at a linear or superlinear rate [but above
## doesn't include 1/sqrt(I) eigenvalues].  difference seems to be
## sqrt rate [would be constant rate if including eigenval].

plot(Is,maxes.total/sqrt(Is))



## 9c phi^TQphi growth rate, tracking down source of correlation between neg and pos parts
source('misc.R')
require(parallel)
I <- 5
Q <- M.obj(I,symm=TRUE) # slow here
vals <- eigen(Q)$val; vecs <- eigen(Q)$vec
neg.idx <- vals< -1e-8;     pos.idx <- vals>1e-8
pos.vecs <- vecs[,pos.idx]#*vals[pos.idx]
neg.vecs <- vecs[,neg.idx]#*vals[neg.idx]
method <- 'a'
m <- 5; n <- 5
qfs <- replicate(1e3, {
    phi.vec <- switch(method,
                      'a'={
                          x <- replicate(I,runif(m),simplify=FALSE)
                          y <- replicate(I,runif(n),simplify=FALSE)
                          phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
                          phi.vec <- matrix(t(phi),ncol=1)
                      },
                      'b'={    
                          ## F <- t(replicate(I,sort(runif(I^2),decr=FALSE)))
                          ## phi <- t(replicate(I,as.numeric(with(list(wt=runif(I^2)),F%*%wt/sum(wt)))))
                          ## phi.vec <- matrix(t(phi),ncol=1)
                          F <- t(replicate(I,sort(runif(I*n),decr=FALSE)))
                          wts <- sample(rep(1:I,n))
                          phi <- sapply(1:I,function(i)F%*%(wts==i)/n)
                          phi.vec <- matrix(t(phi),ncol=1)
                          ## tmp <- phi.vec[1];phi.vec[1] <- phi.vec[I^2];phi.vec[I^2] <- tmp
                          ## jj <- 1;phi.vec[sample(1:I^2,jj)] <- runif(jj)
                          ## phi.vec[1:7] <- .7#runif(2)
                          ## phi.vec[c(1,7)] <- runif(1)
                          idx0 <- ((1:I^2)-1)%%(I+1)==0
                          phi.vec[idx0] <- (runif(sum(diag.idx)))#1#seq(0,1,len=sum(diag.idx))#
                          phi.vec[!diag.idx] <- 1#(runif(sum(!diag.idx)))
                          phi.vec
                          ## FF <- kronecker(diag(I),F)
                      },
                      'c'={
                          ## phi.vec <- matrix(rbinom(I^2,1,.5),ncol=1)
                          phi.vec <- matrix(runif(I^2),ncol=1)
                          ## phi.vec <- sort(runif(1:I^2))
                          })
    ## c(pos=sum((t(vecs[,pos.idx])%*%phi.vec)^2),neg=sum((t(vecs[,neg.idx])%*%phi.vec)^2),total=t(phi.vec)%*%Q%*%phi.vec)
        c(pos=sum((t(pos.vecs)%*%phi.vec)^2),neg=sum((t(neg.vecs)%*%phi.vec)^2),total=t(phi.vec)%*%Q%*%phi.vec/vals[1])
})
plot(qfs['pos',],abs(qfs['neg',]));abline(0,1)
cor(qfs['pos',],qfs['neg',])
## cor(qfs['neg2',],abs(qfs['neg',]));abline(0,1)
phi.vec=sort(runif(I^2))
projs <- replicate(5e3,with(list(phi.vec=sort(runif(I^2))),cbind(t(pos.vecs)%*%phi.vec,t(neg.vecs)%*%phi.vec)))
## cors <- cor(t(projs[,1,]),t(projs[,2,]))
## max(cors)
## with(list(idx=which(cors==max(cors),arr.ind=TRUE)),{a<<-pos.vecs[,idx[1]];b<<-neg.vecs[,idx[2]]})
head(t(qfs))
colMeans(t(qfs))
head(qfs['pos',]/qfs['neg',])
## correlation is 1 when diagonal entries constant or offdiagonal
## entries constant. correlation 1 but not equal.

plot(neg.vecs[,3],ylim=c(-1,1))
points(neg.vecs[,4],col=2)
points(pos.vecs[,1],col=3)

a <- pos.vecs[,1];b <- neg.vecs[,2];c <- neg.vecs[,4]
projs <- replicate(1e2,with(list(phi.vec=sort(runif(I^2))),c(a%*%phi.vec,b%*%phi.vec,c%*%phi.vec)))
cor(t(projs))


## 9d looking at phi.vec as inc vector, except for some indices
## i%(I+1) indices special
source('misc.R')
require(parallel)
I <- 4
Q <- M.obj(I,symm=TRUE) # slow here
vals <- eigen(Q)$val; vecs <- eigen(Q)$vec
neg.idx <- vals< -1e-8;     pos.idx <- vals>1e-8
pos.vecs <- vecs[,pos.idx]#*vals[pos.idx]
neg.vecs <- vecs[,neg.idx]#*vals[neg.idx]
method <- 'b'
m <- 5; n <- 5
idx <- 1:I^2
by.idx <- sapply(idx,function(idx) {
    qfs <- replicate(1e3, {
        F <- t(replicate(I,sort(runif(I*n),decr=FALSE)))
        wts <- sample(rep(1:I,n))
        phi <- sapply(1:I,function(i)F%*%(wts==i)/n)
        phi.vec <- matrix(phi,ncol=1)
        ## tmp <- phi.vec[1];phi.vec[1] <- phi.vec[I^2];phi.vec[I^2] <- tmp
        ## jj <- 1;phi.vec[sample(1:I^2,jj)] <- runif(jj)
        phi.vec[idx] <- .5#runif(1)
        phi.vec
        ## FF <- kronecker(diag(I),F)
        c(pos=sum((t(vecs[,pos.idx])%*%phi.vec)^2),neg=sum((t(vecs[,neg.idx])%*%phi.vec)^2),total=t(phi.vec)%*%Q%*%phi.vec)
    })
    cor(qfs['pos',],qfs['neg',])
})
plot(idx,by.idx)
which(by.idx<.5)
diag.idx <- ((1:I^2)-1)%%(I+1)==0
points(idx[diag.idx],by.idx[diag.idx],col=2)
mean(by.idx[!diag.idx])-mean(by.idx[diag.idx])

dd
## 9d varying substitute quantity
source('misc.R')
require(parallel)
I <- 4
Q <- M.obj(I,symm=TRUE) # slow here
vals <- eigen(Q)$val; vecs <- eigen(Q)$vec
neg.idx <- vals< -1e-8;     pos.idx <- vals>1e-8
pos.vecs <- vecs[,pos.idx]#*vals[pos.idx]
neg.vecs <- vecs[,neg.idx]#*vals[neg.idx]
method <- 'b'
m <- 5; n <- 5
idx <- 1:I^2
repls <- seq(0,1,len=40)
by.repl <- sapply(repls,function(repl) {
    print(repl)
    by.idx <- sapply(idx,function(idx) {
        qfs <- replicate(2e2, {
            F <- t(replicate(I,sort(runif(I*n),decr=FALSE)))
            wts <- sample(rep(1:I,n))
            phi <- sapply(1:I,function(i)F%*%(wts==i)/n)
            phi.vec <- matrix(phi,ncol=1)
            phi.vec[idx] <- repl
            phi.vec
            c(pos=sum((t(vecs[,pos.idx])%*%phi.vec)^2),neg=sum((t(vecs[,neg.idx])%*%phi.vec)^2),total=t(phi.vec)%*%Q%*%phi.vec)
        })
        cor(qfs['pos',],qfs['neg',])
    })
    ## plot(idx,by.idx)
    ## which(by.idx<.5)
    diag.idx <- ((1:I^2)-1)%%(I+1)==0
    ## points(idx[diag.idx],by.idx[diag.idx],col=2)
    mean(by.idx[!diag.idx])-mean(by.idx[diag.idx])
})
plot(repls,by.repl)



## 9d look of phi.vec
I <- 5
n <- 5
Q <- M.obj(I,symm=TRUE) 
phi.vecs <- replicate(1e4,{
    F <- t(replicate(I,sort(runif(I*n),decr=FALSE)))
    wts <- sample(rep(1:I,n))
    phi <- sapply(1:I,function(i)F%*%(wts==i)/n)
    phi.vec <- as.numeric(t(phi))
})
plot(rowMeans(phi.vecs)) #not much pattern, maybe increasing

I <- 5
n <- 5
F <- t(replicate(I,sort(runif(I*n),decr=FALSE)))
wts <- sample(rep(1:I,n))
phi <- sapply(1:I,function(i)F%*%(wts==i)/n)
phi.vec <- as.numeric(t(phi))
t(apply(phi,1,rank)) #ranks usually presrved after mixing

plot(phi.vec)

dd

##9e arg-min phi.vec
source('misc.R')
require(parallel)
I <- 5
m <- n <- 1
diag.idx <- ((1:I^2)-1)%%(I+1)==0
Q <- M.obj(I,symm=TRUE) 
vals <- eigen(Q)$val; vecs <- eigen(Q)$vec
lambda <- max(abs(vals))
neg.idx <- vals< -1e-8;     pos.idx <- vals>1e-8
pos.vecs <- vecs[,pos.idx]#*vals[pos.idx]
neg.vecs <- vecs[,neg.idx]#*vals[neg.idx]
phi.vecs <- replicate(1e4,{
    ## F <- t(replicate(I,sort(runif(I*n),decr=FALSE)))
    ## wts <- sample(rep(1:I,n))
    ## phi <- sapply(1:I,function(i)F%*%(wts==i)/n)
    ## phi.vec <- as.numeric(phi)
    ## phi.vec
    x <- replicate(I,runif(m),simplify=FALSE)
    y <- replicate(I,runif(n),simplify=FALSE)
    phi <- t(sapply(x,function(x.i)sapply(y,function(y.i)concord.mean(x.i,y.i))))
    phi.vec <- as.numeric(t(phi))
} )
qfs <- sapply(1:ncol(phi.vecs),function(i) {
    phi.vec <- phi.vecs[,i]
    c(pos=sum((t(pos.vecs)%*%phi.vec)^2),neg=sum((t(neg.vecs)%*%phi.vec)^2),total=t(phi.vec)%*%Q%*%phi.vec)
})
min.idx <- which.min(abs(qfs['total',]))
min.vec <- phi.vecs[,min.idx]
plot(min.vec,col=diag.idx+1)
min(abs(qfs['total',]))
qfs[c('pos','neg'),min.idx]
image2 <- function(m)image(t(m)[,nrow(m):1])
phi <- t(matrix(min.vec,I))
image2(phi)

## symmetry groups
pos.vecs <- with(eigen(Q),vectors[,values>1e-7])
neg.vecs <- with(eigen(Q),vectors[,values< -1e-7])
neg.qf <- neg.vecs%*%t(neg.vecs)
pos.qf <- pos.vecs%*%t(pos.vecs)
min.vecs <- phi.vecs[,(abs(qfs['total',])<1e-10)]
min.vecs <- min.vecs[,colSums(min.vecs^2)>0]
min.vecs.str <- apply(min.vecs,2,function(v)paste(v,collapse=''));
min.vecs.str <- min.vecs.str[!duplicated(min.vecs.str)] #remove dups
signatures <- sapply(min.vecs.str,function(str) {
    phi <- t(matrix(as.numeric(strsplit(str,'')[[1]]),I))
    ## rotflips <- with(list(rot = function(m) t(m[nrow(m):1,,drop=FALSE]),
    ##     flip=function(m)m[,I:1]),
    ##     lapply(c(list(rot(rot(rot(rot(phi))))),list(phi),list(rot(phi)),list(rot(rot(phi))),list(rot(rot(rot(phi))))),flip))
    ## symmetries form 4-member subgroup of S4
    s4sub <- function(m)list(id=m,rot180=m[I:1,I:1],diag=t(m),antidiag=t(m[I:1,I:1]))
    shifts <- function(m)lapply(1:(I-1),function(s)(m[c((I-s+1):I,1:(I-s)),c((I-s+1):I,1:(I-s))]))
    symmetries <- c(s4sub(phi),s4sub(1-phi))#,shifts(phi),shifts(1-phi))
    symmetries <- lapply(symmetries,shifts)
    symmetries <- do.call(c,symmetries)
    stopifnot(length(unique(round(sapply(symmetries,function(phi)as.numeric(phi)%*%Q%*%as.numeric(phi)),7)))==1)
    ## ## transposed <- t(phi)
    ## ## antitransposed <- t(phi[I:1,I:1])
    ## ## flipped <- phi[I:1,I:1]
    ## negated <- list()#1-phi
    ## shifts <- list()#lapply(1:(I-1),function(s)(phi[c((I-s+1):I,1:(I-s)),c((I-s+1):I,1:(I-s))]))
    ## symmetries <- c(list(phi),list(transposed),list(antitransposed),list(flipped),list(negated),shifts)
    signature <- paste(sort(sapply(symmetries, function(m)paste(m,collapse=''))),collapse='')
})
duplicated(signatures)
sum(duplicated(signatures))
min.vecs.str <- min.vecs.str[!duplicated(signatures)]
min.vecs <- unname(sapply(min.vecs.str,function(x)as.numeric(strsplit(x,'')[[1]])))
min.phis <- apply(min.vecs,2,function(phi.vec)t(matrix(phi.vec,I)),simplify=FALSE)


op <- par(mfrow=c(1,length(min.phis)))
lapply(min.phis,image2)
par(op)
plot(sort(abs(t(pos.vecs)%*%phi.vecs[,min.idx])),
sort(abs(t(neg.vecs)%*%phi.vecs[,min.idx])));abline(0,1)
order(abs(t(pos.vecs)%*%phi.vecs[,min.idx]));order(abs(t(neg.vecs)%*%phi.vecs[,min.idx]))
table(sapply(min.phis,sum))


min.vec <- phi.vecs[,min.idx]
min.vec[!diag.idx ] <- runif(1)
min.vec[diag.idx] <- (0:4)
t(min.vec)%*%Q%*%min.vec

F <- t(replicate(I,sort(runif(I*n),decr=FALSE)))
wts <- sample(rep(1:I,n))
phi <- sapply(1:I,function(i)F%*%(wts==i)/n)
## diag(phi) <- runif(1)

phi[4,] <- runif(1)
vec0 <- as.numeric(phi)
cbind(pos=t(pos.vecs)%*%vec0,neg=t(neg.vecs)%*%vec0)^2*lambda
colSums(cbind(pos=t(pos.vecs)%*%vec0,neg=t(neg.vecs)%*%vec0)^2*lambda)
t(vec0)%*%Q%*%vec0


batch.idx <- sapply(2:sum(diag.idx),function(i)which(diag.idx)[i-1]<(1:I^2) & (1:I^2)<which(diag.idx)[i])
vec0 <- phi.vecs[,sample(1:100,1)]

sub <- c(0,1,0,0,0)
## sub <- sample(sub)
## print(sub)
vec0 <- as.numeric(I^2)
## vec0[!diag.idx ] <- runif(1)
## vec0[diag.idx] <- 1:5#runif(I)#0:4+1
## vec0[!diag.idx ] <- 1:4
for(i in 1:ncol(batch.idx))vec0[batch.idx[,i]] <- sub[i]
vec0[diag.idx] <- c(0,0,0,0,0,1)
cbind(pos=t(pos.vecs)%*%vec0,neg=t(neg.vecs)%*%vec0)^2*lambda
colSums(cbind(pos=t(pos.vecs)%*%vec0,neg=t(neg.vecs)%*%vec0)^2)*lambda
sum((t(pos.vecs)%*%vec0)^2) / sum((t(neg.vecs)%*%vec0)^2)
t(vec0)%*%Q%*%vec0


batch.idx <- sapply(2:sum(diag.idx),function(i)which(diag.idx)[i-1]<(1:I^2) & (1:I^2)<which(diag.idx)[i])
vec0 <- rnorm(I^2)
## vec0 <- phi.vecs[,sample(1:100,1)]
vec0[diag.idx] <- runif(1)
vec0[!diag.idx] <- 1:15
t(vec0)%*%Q%*%vec0
for(j in 1:ncol(batch.idx)) vec0[batch.idx[,j]] <- mean(vec0[batch.idx[,j]])
t(vec0)%*%Q%*%vec0
sapply( 1:ncol(batch.idx),function(j)unique(vec0[batch.idx[,j]]))


I <- 3
Q <- M.obj(I,symm=TRUE) 
vals <- eigen(Q)$val; vecs <- eigen(Q)$vec
lambda <- max(abs(vals))
neg.idx <- vals< -1e-8;     pos.idx <- vals>1e-8
pos.vecs <- vecs[,pos.idx]#*vals[pos.idx]
neg.vecs <- vecs[,neg.idx]#*vals[neg.idx]
diag.idx <- ((1:I^2)-1)%%(I+1)==0
batch.idx <- sapply(2:sum(diag.idx),function(i)which(diag.idx)[i-1]<(1:I^2) & (1:I^2)<which(diag.idx)[i])
vec0 <- rnorm(I^2)
## vec0 <- phi.vecs[,sample(1:100,1)]
vec0[diag.idx] <- runif(1)
diffs <- sapply(1:sum(!diag.idx),function(jj) {
    vec0[!diag.idx] <- 1:jj
    a <- t(vec0)%*%Q%*%vec0
    for(j in 1:ncol(batch.idx)) vec0[batch.idx[,j]] <- mean(vec0[batch.idx[,j]])
    b <- t(vec0)%*%Q%*%vec0
    b-a
    })
which(abs(diffs)<1e-8)



I <- 5
n <- 1
Q <- M.obj(I,symm=TRUE)
vals <- eigen(Q)$val; vecs <- eigen(Q)$vec
lambda <- max(abs(vals))
neg.idx <- vals< -1e-8;     pos.idx <- vals>1e-8
pos.vecs <- vecs[,pos.idx]#*vals[pos.idx]
neg.vecs <- vecs[,neg.idx]#*vals[neg.idx]
phi <- matrix(runif(I^2),nrow=I)
phi.vec <- as.numeric(phi)
t(phi.vec)%*%Q%*%phi.vec
idx1 <- sample(1:I,1);idx2 <- sample(1:I,1)
tmp <- phi[idx1,]
phi[idx1,] <- phi[idx2,]
phi[idx2,] <- tmp
tmp <- phi[,idx1]
phi[,idx1] <- phi[,idx2]
phi[,idx2] <- tmp
phi.vec <- as.numeric(phi)
t(phi.vec)%*%Q%*%phi.vec
## Q invariant to swapping row i,j and column i,j, corresponding to
## swapping indices on x,y; related to p,q,r,s symmetry? No an
## additional symmetry--see 6e

source('misc.R')
require(parallel)
I <- 5
n <- 1
Is <- round(seq(3,30,len=30))
by.I <- sapply(Is,function(I) {
## by.I <- mclapply(Is,mc.cores=detectCores()-4,FUN=function(I) {
    print(I)
    diag.idx <- ((1:I^2)-1)%%(I+1)==0
    Q <- M.obj(I,symm=TRUE)
    ## vecs <- matrix(runif(I^2*1e3),ncol=I^2)
    ## projs <- apply(vecs,1,function(r) t(r)%*%Q%*%r)
    ## print(vecs[which.max(abs(projs)),])
    ## max(abs(projs))
    vals <- eigen(Q)$val; vecs <- eigen(Q)$vec
    lambda <- max(abs(vals))
    neg.idx <- vals< -1e-8;     pos.idx <- vals>1e-8
    pos.vecs <- vecs[,pos.idx]#*vals[pos.idx]
    neg.vecs <- vecs[,neg.idx]#*vals[neg.idx]
    vec0 <- (pos.vecs[,2]>0)+0
    t(vec0)%*%Q%*%vec0
    ## sum(apply(pos.vecs,2,function(c)sum(c[c>0])^2))
})
## by.I <- simplify2array(by.I)
plot(Is,by.I)
plot(Is,by.I/Is^(1))
## for this phi.vec quad form grows faster than I^(1/2), maybe even I
## this phi.vec is in [0,1] but maybe not a valid phi.vec

## 9f mixing Q with phi parts--Q with wts 1. wts%*%Q%*%wts also has
## projection-like eigenvalues. eigenvals of Qwtd same as Q/n. 2. for
## n=1 (wts.mat=wts.mat.big a perm matrix) evecs of Qwtd are
## permutations of evecs of Q: Q.wtd.P^-1.v = P^-1.Q.P.P^-1.v = lambda
## P^-1.v, where P is kronecker(E,P'), P' is a permutation matrix, E
## is the IxI identity matrix, and v' = P^-1v. So in fact the evecs of
## Qwtd aren't just any permutations of the evecs of Q but each
## I-segment is permuted within itself. 3. For general n: wts.mat=P'
## is rectangular I x (In) so wts.mat.big=P is I^2 x (I^2n) but P
## still has t(P) for a right inverse, so t(P).v is an eigenvector of
## Qwtd whenever v is an eigenvec of Q, with the same eigenval. Ie
## since t(P) is its inverse, Q and Q.wtd are similar so have the same
## char poly, same evals. 5. for 2--why not symmetric in m,n? m comes
## from F.vec
source('misc.R')
I <- 3
n <- 3
Q <- M.obj(I,symm=TRUE)
F <- t(replicate(I,sort(runif(I*n),decr=FALSE)))
wts <- sample(rep(1:I,n))
## phi <- sapply(1:I,function(i)F%*%(wts==i)/n)
## phi.vec <- as.numeric(t(phi))
## t(phi.vec)%*%Q%*%phi.vec
wts.mat <- t(sapply(1:I,function(i)(wts==i)))/n
F.vec <- as.numeric(t(F))
## phi.vec <- kronecker(diag(I),wts.mat)%*%F.vec
## t(phi.vec)%*%Q%*%phi.vec
wts.mat.big <- kronecker(diag(I),wts.mat)
Qwtd <- t(wts.mat.big)%*%Q%*%wts.mat.big
t(F.vec)%*%Qwtd%*%F.vec 
unique(round(eigen(Q)$val,5)) 
unique(round(eigen(Qwtd)$val,5))*n
t(wts.mat.big)%*%eigen(Q)$vec/ Qwtd%*%t(wts.mat.big)%*%eigen(Q)$vec # 3

dd
## Qwtd2 <- t(kronecker(diag(I),wts.mat)[,1:(I*n)])%*%Q%*%(kronecker(diag(I),wts.mat)[,1:(I*n)])
## unique(round(eigen(Qwtd2)$val,5))*n

## 9g mixing Q with phi parts--Q with wts, and F.vec
source('misc.R')
I <- 3
n <- 1
Q <- M.obj(I,symm=TRUE)
wts <- sample(rep(1:I,n))
## phi <- sapply(1:I,function(i)F%*%(wts==i)/n)
## phi.vec <- as.numeric(t(phi))
## t(phi.vec)%*%Q%*%phi.vec
wts.mat <- t(sapply(1:I,function(i)(wts==i)))/n
## F <- t(replicate(I,sort(runif(I*n),decr=FALSE)))
F <- sample(rep(1:I,n))
F.mat <- t(sapply(1:I,function(i)(F==i)))/n
F.vec <- as.numeric(t(F.mat))
## phi.vec <- kronecker(diag(I),wts.mat)%*%F.vec
## t(phi.vec)%*%Q%*%phi.vec
wts.mat.big <- kronecker(diag(I),wts.mat)
Qwtd <- t(wts.mat.big)%*%Q%*%wts.mat.big
## t(F.vec)%*%Qwtd%*%F.vec 
unique(round(eigen(Q)$val,5)) 
## unique(round(eigen(Qwtd)$val,5))*n
M.cumsum <- with(list(m=matrix(nrow=I*n,ncol=I*n)),(row(m)>=col(m))+0)
M.cumsum <- kronecker(diag(I),M.cumsum)
Qwtdcs <- t(M.cumsum)%*%Qwtd%*%M.cumsum
eigen(Qwtdcs)$vec
t(M.cumsum%*%F.vec)%*%Qwtd%*%M.cumsum%*%F.vec

## 9gg m=n=1 case, checking formulas [ipad auc#4]
source('misc.R')
Qi <- function(p,q,r,s,I=3)
    1/I^2*sum(apply(combn(1:4,2),2,function(idx)abs(diff(c(p,q,r,s)[idx]))<1e-7)) -1/(2*I)*sum(apply(combn(1:4,3),2,function(idx)sum(abs(diff(c(p,q,r,s)[idx])))<1e-7)) - 4/I^3
I <- 5
n <- 1
Q <- M.obj(I,symm=TRUE)
M.cumsum <- with(list(m=matrix(nrow=I*n,ncol=I*n)),(row(m)>=col(m))+0)
M.cumsum <- kronecker(diag(I),M.cumsum)
## diffs <- replicate(1e3,{
F <- sample(rep(1:I,n))
F.mat <- t(sapply(1:I,function(i)(F==i)))/n
F.vec <- as.numeric(t(F.mat))
f <- apply(F.mat,1,function(r)which(r==1))
phi.vec <- M.cumsum%*%F.vec
pqrs <- expand.grid(structure(rep(list(1:I),4),names=c('p','q','r','s')))
t(phi.vec)%*%Q%*%phi.vec

phi <- matrix(phi.vec,nrow=I,byrow=TRUE)
phi2 <- phi[,c(2,1,3:I)]
phi2 <- phi2[c(2,1,3:I),]
phi.vec2 <- as.numeric(t(phi2))
t(phi.vec2)%*%Q%*%phi.vec2

dd
## } )
## sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])*(s>=f[q])*Qi(p,q,r,s,I=I))))
## cbind(pqrs,a=apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(p==q))))
## sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(p==q))))
## sum((I+1-f)^2)
## I*(I+1)*(2*I+1)/6
## sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(r==s))))
## I^2*(I+1) - sum(apply(expand.grid(1:I,1:I),1,function(x)max(x[1],x[2])))
## sum(apply(expand.grid(1:I,1:I),1,function(x)max(x[1],x[2])))
## sum((1:I)^2) + sum(sapply(1:(I-1),function(p)sum((p+1):I)))
## sum((1:I)^2) + 1/2*I*(I^2-1) - 1/2*sum((2:I)*(1:(I-1)))
## sum((1:I)^2) + 1/2*I*(I^2-1) -1/2*sum((1:(I-1))^2)-1/4*I*(I-1)
## I^3/2+3/4*I^2+I/4-1/12*I*(I+1)*(2*I+1)
## I^3/3+I^2/2+I/6
## sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(p==r))))
## sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(q==s))))
## A <- sum(I +1-f)
## A*sum(f<=(1:I))
## sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(p==s))))
## sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(q==r))))
## sum((1:I)*(I+1-f))
## sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(p==q)&(q==r))))
## sum((f<=(1:I))*(I+1-f))
## sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(p==r)&(r==s))))
## sum((f<=(1:I))*(1:I))
A <- sum(I +1-f)
1/I^2*(I*(I+1)*(2*I+1)/3 + I*(I+1)*sum(f<=(1:I)) + I*(I+1)^2 - 2*sum((1:I)*f)) - 1/I*(sum((f<=(1:I))*(I+1-f)) + sum((f<=(1:I))*(1:I))) - 4/I^3*A^2
(I+1)*(2*I+1)/3/I - 2/I^2*sum((1:I)*f) + 1/I*sum((f-(1:I))*(f<=(1:I)))
c((I+1)*(2*I+1)/3/I, - 2/I^2*sum((1:I)*f), + 1/I*sum((f-(1:I))*(f<=(1:I))))

## 9ggg m=n=1 case, checking formulas [ipad auc#4] -- redoing with general
## F. Previous version was not general, in assuming F is an orthogonal
## matrix. row of F can have a 1 in any position (even the same
## between rows) or be all 0s.
source('misc.R')
Qi <- function(p,q,r,s,I=3)
    1/I^2*sum(apply(combn(1:4,2),2,function(idx)abs(diff(c(p,q,r,s)[idx]))<1e-7)) -1/(2*I)*sum(apply(combn(1:4,3),2,function(idx)sum(abs(diff(c(p,q,r,s)[idx])))<1e-7)) - 4/I^3
I <- 5
n <- 1
Q <- M.obj(I,symm=TRUE)
M.cumsum <- with(list(m=matrix(nrow=I*n,ncol=I*n)),(row(m)>=col(m))+0)
M.cumsum <- kronecker(diag(I),M.cumsum)
## F <- sample(rep(1:I,n))
## F.mat <- t(sapply(1:I,function(i)(F==i)))/n
F.mat <- matrix(0,nrow=I,ncol=I); for(i in 1:I) F.mat[i,sample(0:I,1)] <- 1
F.vec <- as.numeric(t(F.mat))
f <- apply(F.mat,1,function(r)which(r==1))
f[sapply(f,length)==0] <- I+1
f <- as.numeric(f)
phi.vec <- M.cumsum%*%F.vec
pqrs <- expand.grid(structure(rep(list(1:I),4),names=c('p','q','r','s')))
t(phi.vec)%*%Q%*%phi.vec
sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])*(s>=f[q])*Qi(p,q,r,s,I=I))))
## cbind(pqrs,a=apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(p==q))))
sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(p==q))))
## sum((I+1-f)^2)
## I*(I+1)*(2*I+1)/6
I*(I+1)^2-2*(I+1)*sum(f)+sum(f^2) #1
sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(r==s))))
## I^2*(I+1) - sum(apply(expand.grid(1:I,1:I),1,function(x)max(f[x[1]],f[x[2]])))
I^2*(I+1) - sum(sort(f)*(2*(1:I)-1)) #2
sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(p==r))))
sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(q==s))))
A <- sum(I +1-f)
A*sum(f<=(1:I)) #3,4
sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(p==s))))
sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(q==r))))
## sum((1:I)*(I+1-f))
sum(ecdf(f)(1:I)*I*(I+1-f)) #5,6
sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(p==q)&(q==r))))
sum((f<=(1:I))*(I+1-f)) #a,b
sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])&(s>=f[q])&(p==r)&(r==s))))
## sum((f<=(1:I))*(1:I))
## A <- sum(I +1-f)
sum((f<=(1:I))*ecdf(f)(1:I)*I) #c,d
sum(apply(pqrs,1,function(x)with(as.list(x),(r>=f[p])*(s>=f[q]))))
(I*(I+1) - sum(f))^2

1/I^2*(I*(I+1)^2-2*(I+1)*sum(f)+sum(f^2) + I^2*(I+1) - sum(sort(f)*(2*(1:I)-1)) + 2*(A*sum(f<=(1:I))) + 2*sum(ecdf(f)(1:I)*I*(I+1-f))) - 1/I*(sum((f<=(1:I))*(I+1-f)) + sum((f<=(1:I))*ecdf(f)(1:I)*I)) - 4/I^3*((I*(I+1) - sum(f))^2)



## maximizing permutations
require(gtools)
I <- 9
n <- 1
## fs <- replicate(1e5,sample(1:I))
Is <- 3:9
obj.f <- function(f,I)(I+1)*(2*I+1)/3/I - 2/I^2*sum((1:I)*f) + 1/I*sum((f-(1:I))*(f<=(1:I)))
obj.delta <- function(delta,I)sum(((delta<0)-2*(1:I)/I)*delta)/I
perms <- lapply(Is,function(I) {
    ## fs <- t(permutations(I,I))
    ## diffs <- apply(fs,2,function(f) 2/I^2*sum((1:I)*f) - 1/I*sum((f-(1:I))*(f<=(1:I))))
    fs <- t(permutations(I,I))
    deltas <- fs-(1:I)
    diffs <- apply(deltas,2,function(delta) obj.delta(delta,I))
    list(range=range(diffs),
         max=deltas[,which.max(diffs)],
         min=deltas[,which.min(diffs)],
         maxabs=deltas[,which.max(abs(diffs))])
})
names(perms) <- Is
sapply(perms,function(perm)perm$min)
sapply(perms,function(perm)perm$max)
sapply(perms,function(perm)perm$maxabs)
## for f (not delta version): can kind of see a pattern in argmins,
## not for argmaxes

## delta representation
sapply(perms,function(perm)perm$min - (1:length(perm$min)))

cycle.decomposition <- function(perm) {
    names(perm) <- 1:length(perm)
    cycles <- list()
    while(length(perm)>0) {
        cycle <- as.integer(names(perm[1]))
        repeat {
            last.entry <- cycle[length(cycle)]
            cycle <- c(cycle,perm[as.character(last.entry)])
            if(cycle[1]==cycle[length(cycle)])break
        }
        cycle <- cycle[-length(cycle)]
        perm <- perm[!(names(perm)%in%as.character(cycle))]
        cycles <- c(cycles,list(unname(cycle)))
    }
    cycles
}
perm <- structure(sample(1:7),names=1:7)
cycle.decomposition(perm)

sapply(perms,function(perm)cycle.decomposition(perm$min))
sapply(perms,function(perm)cycle.decomposition(perm$max))

f <- perms[[3]]$min
I <- length(f)
f0 <- 1:I
2/I^2*sum((1:I)*f0) + 1/I*sum(((1:I)-f0)*(f0<=(1:I)))
2/I^2*sum((1:I)*f) + 1/I*sum(((1:I)-f)*(f<=(1:I)))
-1+3/I-2/I^2

## [DETOUR] delta as a distance on permutations isomorphic to dot prod with
## identity perm? 
require(gtools)
I <- 5
## fs <- replicate(1e5,sample(1:I))
fs <- t(permutations(I,I))
delta <- colSums(abs(fs-(1:I)))
delta <- colSums((fs-(1:I))^2)
dotprods <- apply(fs,2,function(perm)sum(perm*(1:I)))
plot(delta,dotprods)
## for equally spaced numbers in f, f%*%[I] is strictly monotonically
## decreasing in |f-[I]|_{L2}

## comparing f to other permutations than [I]
I <- 10
x <- rnorm(I); y <- rnorm(I)
x <- 1:I; y <- (1:I)/20 # equally vs inequally spaced
## x <- 1:I; y <- (1:I)^3 
## x <- (1:I)^3; y <- (1:I)^3 
pairs <- replicate(1e4, {
x.pi <- sample(x); y.pi <- sample(y)
c(L2=sum(abs(rank(x.pi)-rank(y.pi))^2),dot=x.pi%*%y.pi) #L1 vs L2
})
plot(pairs['L2',],pairs['dot',])
maxmin <- simplify2array(by(data.frame(t(pairs)),pairs['L2',],function(df)c(L2=unique(df$L2),max=max(df$dot),min=min(df$dot))))
lines(maxmin['L2',],maxmin['max',],col=2)
lines(maxmin['L2',],maxmin['min',],col=2)
## again for equally spaced x and equally spaced y (spacings can
## differ) pi(x)%*%pi(y) is monotonically decreasing in L2
## distance. [ie residual inverse to projection.] when not equally
## spaced, seems that the max and the min pi(x)%*%pi(y) for fixed L2
## distance are almost but not quite decreasing? They do seem to
## decrease for L1 distance.  [end detour]

dd
    
## 9h m=n=1 case still, assuming F orthogonal, delta formula for
## objective [AUC #5]
I <- 7
f0 <- 1:I
f <- c(I,2:(I-1),1)
obj <- function(f,I)(I+1)*(2*I+1)/3/I - 2/I^2*sum((1:I)*f) + 1/I*sum((f-(1:I))*(f<=(1:I)))
obj(f,I);obj(f0,I)
delta <- f-f0
sum((-(delta<0)*delta/I*(1-2*(1:I)/I)) + (delta>0)*2*(1:I)*delta/I^2)

source('misc.R')
obj.delta <- function(delta,I)sum(((delta<0)-2*(1:I)/I)*delta)/I
I <- 6
n <- 1
Q <- M.obj(I,symm=TRUE)
M.cumsum <- with(list(m=matrix(nrow=I*n,ncol=I*n)),(row(m)>=col(m))+0)
M.cumsum <- kronecker(diag(I),M.cumsum)
F <- sample(rep(1:I,n))
F.mat <- t(sapply(1:I,function(i)(F==i)))/n
F.vec <- as.numeric(t(F.mat))
f <- apply(F.mat,1,function(r)which(r==1))
(I+1)*(2*I+1)/3/I - 2/I^2*sum((1:I)*f) + 1/I*sum((f-(1:I))*(f<=(1:I)))
## c((I+1)*(2*I+1)/3/I, - 2/I^2*sum((1:I)*f), + 1/I*sum((f-(1:I))*(f<=(1:I))))
delta <- f-(1:I)
sum((delta<0)*delta/I*(1-2*(1:I)/I) - (delta>0)*2*(1:I)*delta/I^2)
phi.vec <- M.cumsum%*%F.vec
t(phi.vec)%*%Q%*%phi.vec
(sum((delta<0)*delta) - sum(2*(1:I)*delta)/I)/I
obj.delta(delta,I)
## delta.star <- c(-(I-1),-(I-2),rep(0,I-4),I-2,I-1)
obj.delta(c(-(I-1),-(I-2),rep(0,I-4),I-2,I-1),I)
obj.delta(c(-(I-1),rep(0,I-2),I-1),I)

jstar <- floor(I/4+1/2)
2/I^2 * (1:jstar)*(I-(2*(1:jstar)-1))
(delta>0)*2*(1:I)*delta/I^2
2/I^2*js*(I-(2*js-1))
(delta<0)*delta/I*(1-2/I*(1:I)) 
-(I-(2*(1:jstar)-1))*(1-2/I*(I+1-(1:jstar)))/I
js <- 1:jstar
-1/I*(I+1-2*js)*(-1-2/I+2/I*js)
sum(-1/I*((I+1)*(-1-2/I) +  (2/I*(I+1)-2*(-1-2/I))*js - 4/I*js^2))
-1/I*(-(I+1)*(1+2/I)*jstar + sum((4+6/I)*js) - sum(4/I*js^2))
sum((delta<0)*delta/I*(1-2/I*(1:I))  - (delta>0)*2*(1:I)*delta/I^2)
1/I*(I+1)*(1+2/I)*jstar + sum(js)*(-6/I-8/I^2) + sum(js^2)*8/I^2
1/I*(I+1)*(1+2/I)*jstar  - (6/I+8/I^2)*1/2*jstar*(jstar+1) + 8/I^2/6*jstar*(jstar+1)*(2*jstar+1)


## 9i relating delta form for objective [AUC#5] to old form [AUC#1,4a above]

source('misc.R')
obj.delta <- function(delta,I)sum(((delta<0)-2*(1:I)/I)*delta)/I
I <- 6
n <- 1
Q <- M.obj(I,symm=TRUE)
M.cumsum <- with(list(m=matrix(nrow=I*n,ncol=I*n)),(row(m)>=col(m))+0)
M.cumsum <- kronecker(diag(I),M.cumsum)
## F <- sample(rep(1:I,n))
F.mat <- matrix(0,nrow=I,ncol=I); for(i in 1:I) F.mat[i,sample(1:I,1)] <- 1
F.vec <- as.numeric(t(F.mat))
f <- apply(F.mat,1,function(r)which(r==1))
## (I+1)*(2*I+1)/3/I - 2/I^2*sum((1:I)*f) + 1/I*sum((f-(1:I))*(f<=(1:I)))
delta <- f-(1:I)
## sum((delta<0)*delta/I*(1-2*(1:I)/I) - (delta>0)*2*(1:I)*delta/I^2)
phi.vec <- M.cumsum%*%F.vec
phi <- matrix(phi.vec,nrow=I,byrow=TRUE)
## t(phi.vec)%*%Q%*%phi.vec
(sum((delta<0)*delta) - sum(2*(1:I)*delta)/I)/I

a <- rowSums(phi)
theta.hat <- sum(phi)/I^2
V.x <- rowSums(phi)/I
V.y <- colSums(phi)/I
V.bar.x <- V.x; V.bar.y <- V.y
Is <- 1:I
F <- ecdf(a)
## sum((V.x + V.y - 2*theta.hat)*(V.x + V.y - diag(phi)))
## sum(V.x^2) + 2*sum(V.x*V.y) + sum(V.y^2) - sum((V.x+V.y)*diag(phi)) - 4*I*theta.hat^2 + 2*theta.hat*sum(diag(phi))
sum(a^2)/I^2 + 2/I*sum(a*(1-F(I-(Is)))) + sum((1-F(I-(Is)))^2) - 1/I*sum(a*(a>I-(Is))) - sum((1-F(I-(Is)))*(a>I-(Is))) - 4*sum(a)^2/I^3 + 2/I^2*sum(a)*sum(a>I-Is)




## 11 arithmetic progression as isotropic vectors, new dispersion measure

## 11a Figure out why matrix(1:9,3) is isotropic under Q.  Subspace mapped
## by Q form to 0. 1. Any permutation of phi=matrix(1:I^2,I) is mapped
## to 0 by quad form Q. But not in nullspace, ie, the pos and neg
## parts of Q are equal. 2. More generally, also holds for phi where
## rows are arithmetic progressions, with the same step for all
## rows. step=0 is the constant phi case observed earlier. But,
## constant case is in the nullsapce of Q, in general only the quad
## form mapping it to 0. 3. Still holds when rows and cols are
## permuted independently: The result is still a matrix with arith
## progs for rows or cols with common step. In general the same
## permutation must be applied to row and cols for invariance of Q
## (ie, shuffling the cluster indices). 4. Still holds when phi is
## formed by mixing columns of an I x In matrix the rows of which are
## all arithmetic progressions with the same step: Again, cols are
## arithmetic progressions with the same step. Independently permuting
## rows/cols also holds.  e <-
## Vectorize(function(j,I)unlist(within(list(v=rep(0,I^2)),v[j] <-
## 1))) 5. Proof idea: Fixing r,s, for all p,q, the Q(p,q,r,s) block
## element shows up the same number of times in each block of Qwtdcs,
## so ones.Qwtdcs.ones being 0 reduces to sum_{p,q}Q(p,q,r,s)=0, which
## by symmetry is sum_{r,s}Q(p,q,r,s)=0, and we know each block sums
## to 0.
sum(abs(replicate(1e2,{   
    I <- sample(3:10,1)
    n <- 3
    Q <- M.obj(I,symm=TRUE)
    ## phi <- t(matrix(1:I^2,I)) #1
    step <- runif(1)
    phi <- t(replicate(I,runif(1)+step*(1:I))) #2
    ## orthog <- diag(I)[,sample(1:I)]
    ## phi <- orthog%*%phi
    ## phi <- matrix(1:(I^2*n),nrow=I,byrow=TRUE) 
    phi <- t(replicate(I,runif(1)+step*(1:(I*n)))) #4
    idx <- sample(rep(1:I,n))
    P <- sapply(1:I,function(i)idx==i)+0
    phi <- phi%*%P
    phi <- phi[sample(1:I),sample(1:I)] #3
    v <- as.numeric(phi)
    t(v)%*%Q%*%v
    ## wts <- sample(rep(1:I,n)) #5
    ## wts.mat <- t(sapply(1:I,function(i)(wts==i)))/n
    ## wts.mat.big <- kronecker(diag(I),wts.mat)
    ## Qwtd <- t(wts.mat.big)%*%Q%*%wts.mat.big
    ## M.cumsum <- with(list(m=matrix(nrow=I*n,ncol=I*n)),(row(m)>=col(m))+0)
    ## M.cumsum <- kronecker(diag(I),M.cumsum)
    ## Qwtdcs <- t(M.cumsum)%*%Qwtd%*%M.cumsum
    ## F.vec <- matrix(1,nrow=I^2*n)
    ## ## t(M.cumsum%*%F.vec)%*%Qwtd%*%M.cumsum%*%F.vec
    ## t(F.vec)%*%Qwtdcs%*%F.vec
    ## eigen(Qwtd)$val
    ## t(F.vec)%*%Qwtd%*%F.vec
    ## t(v)%*%t(orthog)%*%Q%*%orthog%*%v
    ## pos.vecs <- with(eigen(Q),vectors[,values>1e-7])
    ## neg.vecs <- with(eigen(Q),vectors[,values< -1e-7])
    ## neg.qf <- neg.vecs%*%t(neg.vecs)
    ## pos.qf <- pos.vecs%*%t(pos.vecs)
    ## t(v)%*%pos.qf%*%v;     t(v)%*%neg.qf%*%v
})))


## Tried to extend Q-basis. Arith progressions with steps alternating
## 1,-1,1,... is, for even I (and variants below), Q-orthogonal to
## itself and to C. But not Q-orthog to B basis vectors.
I <- 6
Q <- M.obj(I,symm=TRUE)
B <- lapply(1:I,function(i)within(list(m=matrix(0,I,I)),m[i,] <- 1)[[1]] / sqrt(I))
C <- with(list(m=matrix(1,nrow=I)%*%matrix(1:I,1)),m-rowMeans(m))
phi.vec <- (I^2:1)-4
phi.vec%*%Q%*%phi.vec

phi <- t(sapply(1:I,function(i)(1:I)))
phi2 <- t(sapply(1:I,function(i)(-1)^(i<=(I/2))*(1:I)))
phi2 <- t(sapply(1:I,function(i)(-1)^(i)*(1:I)))
sum(phi*phi2)
phi2.vec <- c(t(phi2))
phi2.vec%*%Q%*%phi2.vec
sapply(1:I,function(i)c(t(B[[i]]))%*%Q%*%phi2.vec)

phi2.vec%*%Q%*%phi2.vec
step <- runif(1)
phi2 <- t(sapply(1:I,function(i)(-1)^i*step*(1:I))) 
phi2.vec <- c(t(phi2))
phi2.vec%*%Q%*%phi2.vec
c(t(B[[2]]))%*%Q%*%phi2.vec

dd


## 11b Basis for V=subspace consisting of matrices with rows that are
## arithmetic progressions with the same step. Ignoring mixing step
## first. [update: not just a orthogonal wrt dot product but also wrt
## Q i.e. (B[[i]], Q.B[[j]])=0 for i\neq j, and (C, Q.B[[j]])=0]
I <- sample(3:10,1)
Q <- M.obj(I,symm=TRUE)
B <- lapply(1:I,function(i)within(list(m=matrix(0,I,I)),m[i,] <- 1)[[1]] / sqrt(I))
C <- with(list(m=matrix(1,nrow=I)%*%matrix(1:I,1)),m-rowMeans(m))
## sum(C^2) - I^2*(I^2-1)/12 # fla check
C <- C/sqrt(sum(C^2))
basis <- c(B,list(C))
sum(abs(t(basis.mat)%*%Q%*%basis.mat))

phi <- replicate(I,sort(runif(I)))
proj <- Reduce('+',lapply(basis,function(m)sum(phi*m)*m))
jklm <- expand.grid(j=1:I,k=1:I,l=1:I,m=1:I)
sum(apply(jklm,1,function(jklm)with(as.list(jklm),  phi[j,k]*phi[j,l]/I^2 + 12/I^2/(I^2-1)*(k-(I+1)/2)*phi[j,k]*(m-(I+1)/2)*phi[l,m]    )))
sum(proj^2)
dd
## sum(apply(jkl,1,function(jkl)with(as.list(jkl),  phi[j,k]*phi[j,l]*(I^3+3*I^2+(5-6*(k+l))*I+3*(2*k-1)*(2*l-1))/I^2/(I^2-1)    )))
## sum(apply(jkl,1,function(jkl)with(as.list(jkl),  phi[j,k]*phi[j,l]*(1/I+12*(k-(I+1)/2)*(l-(I+1)/2)/I^2/(I^2-1))    )))
## sum(proj^2)
## sum(Reduce('+',lapply(B,function(m)sum(phi*m)*m))^2)
## sum(apply(jkl,1,function(jkl)with(as.list(jkl),  phi[j,k]*phi[j,l]/I   )))
## sum((C*phi)^2)
## sum(apply(jkl,1,function(jkl)with(as.list(jkl),  phi[j,k]*phi[j,l]*(k-(I+1)/2)*(l-(I+1)/2)    ))) * 12/I^2/(I^2-1)
resid <- phi-proj
## proj belongs to V
stopifnot(length(unique(round(as.numeric(apply(proj,1,diff)),8))) == 1)
## resid orthogonal to V
stopifnot(sum(abs(sapply(basis,function(m)sum(m*resid)))) <1e-10)
phi.vec <- as.numeric(phi)
proj.vec <- as.numeric(proj); resid.vec <- as.numeric(resid)
phi.vec%*%Q%*%phi.vec
proj.vec%*%Q%*%proj.vec + resid.vec%*%Q%*%resid.vec + 2*resid.vec%*%Q%*%proj.vec
c(proj.vec%*%Q%*%proj.vec, resid.vec%*%Q%*%resid.vec, 2*resid.vec%*%Q%*%proj.vec)

## 11c growth rate of orthogonal components of phiQphi
source('misc.R')
require(parallel)
Is <- 3:30
## parts <- lapply(Is, function(I) {
parts <- mclapply(Is, mc.cores=detectCores()-3,FUN=function(I) {
    print(I)
    Q <- M.obj(I,symm=TRUE)
    B <- lapply(1:I,function(i)within(list(m=matrix(0,I,I)),m[i,] <- 1)[[1]] / sqrt(I))
    C <- with(list(m=matrix(1,nrow=I)%*%matrix(1:I,1)),m-rowMeans(m))
    C <- C/sqrt(sum(C^2))
    basis <- c(B,list(C))
    ## pos.vecs <- with(eigen(Q),vectors[,values>1e-7])
    ## neg.vecs <- with(eigen(Q),vectors[,values< -1e-7])
    ## phi.vec <- pos.vecs[,1]
    ## phi <- matrix(phi.vec,I,byrow=TRUE)
    ## proj <- Reduce('+',lapply(basis,function(m)sum(phi*m)*m))
    ## resid <- phi-proj
    ## phi.vec <- as.numeric(phi)
    ## proj.vec <- as.numeric(proj); resid.vec <- as.numeric(resid)
    ## c(proj=proj.vec%*%Q%*%proj.vec, resid=resid.vec%*%Q%*%resid.vec, cross=2*resid.vec%*%Q%*%proj.vec)
    ## phi.vec%*%Q%*%phi.vec
    replicate(1e2, {
        phi <- t(replicate(I,sort(runif(I))))
        ## phi <- matrix(c(rep(0,I-1),1),nrow=I,ncol=I,byrow=TRUE) #2
        ## phi <- matrix(runif(I^2),I)
        ## phi <- matrix(rbinom(I^2,1,.5),I)
        proj <- Reduce('+',lapply(basis,function(m)sum(phi*m)*m))
        resid <- phi-proj
        phi.vec <- as.numeric(phi)
        proj.vec <- as.numeric(proj); resid.vec <- as.numeric(resid)
        c(proj=proj.vec%*%Q%*%proj.vec, resid=resid.vec%*%Q%*%resid.vec, cross=2*resid.vec%*%Q%*%proj.vec,resid.norm2=sum(resid.vec^2),proj.norm2=sum(proj.vec^2))
    })
})

resid <- sapply(parts,function(m)m['resid',])
cross <- sapply(parts,function(m)m['cross',])
matplot(t(cross),pch='.',col=rgb(1,0,0,.2),cex=5,ylim=range(c(resid,cross)),x=Is)
matplot(t(resid),pch='.',cex=5,col=rgb(0,1,0,.2),add=TRUE)

resid.max <- apply(resid,2,function(x)max(abs(x)))
cross.max <- apply(cross,2,function(x)max(abs(x)))
plot(Is,resid.max,ylim=range(c(resid.max,cross.max)))
points(Is,cross.max,col=2)

## with(apply(parts,c(1,3)
with(list(parts=simplify2array(parts)),matplot(t(parts['resid',,]/parts['resid.norm2',,]),pch='.'))

with(list(parts=simplify2array(parts)),matplot(t(parts['proj.norm2',,]/parts['resid.norm2',,]),pch='.')) # ratio |proj|^2/|resid|^2 doesnt go to 0

total <- sapply(parts,colSums)
total.max <- apply(total,2,function(x)max(abs(x)))
plot(Is,total.max/sqrt(Is))
## total seems to grow at a less than sqrt(I) rate when phi has random
## entries, but taking the conterexample of phi.vec=pos.vecs[,1]>0
## shows the max actually grows faster
parts <- simplify2array(parts)
plot(Is,colSums(parts)/sqrt(Is))


## 11d growth rate of parts of projection
source('misc.R')
require(parallel)
Is <- 3:30
## parts <- lapply(Is, function(I) {
parts <- mclapply(Is, mc.cores=detectCores()-4,FUN=function(I) {
    print(I)
    Q <- M.obj(I,symm=TRUE)
    B <- lapply(1:I,function(i)within(list(m=matrix(0,I,I)),m[i,] <- 1)[[1]] / sqrt(I))
    C <- with(list(m=matrix(1,nrow=I)%*%matrix(1:I,1)),m-rowMeans(m))
    C <- C/sqrt(sum(C^2))
    basis <- c(B,list(C))
    ## jklm <- expand.grid(j=1:I,k=1:I,l=1:I,m=1:I)
    ## pos.vecs <- with(eigen(Q),vectors[,values>1e-7])
    ## neg.vecs <- with(eigen(Q),vectors[,values< -1e-7])
    ## phi.vec <- pos.vecs[,1]
    ## phi <- matrix(phi.vec,I,byrow=TRUE)
    ## proj <- Reduce('+',lapply(basis,function(m)sum(phi*m)*m))
    ## resid <- phi-proj
    ## phi.vec <- as.numeric(phi)
    ## proj.vec <- as.numeric(proj); resid.vec <- as.numeric(resid)
    ## c(proj=proj.vec%*%Q%*%proj.vec, resid=resid.vec%*%Q%*%resid.vec, cross=2*resid.vec%*%Q%*%proj.vec)
    ## phi.vec%*%Q%*%phi.vec
    replicate(5e2, {
        ## phi <- t(replicate(I,sort(runif(I))))
        phi <- t(replicate(I,sort(sample(1:I,I,replace=TRUE))))/I
        ## phi <- matrix(runif(I^2),I)
        ## phi <- matrix(rbinom(I^2,1,.5),I)
        ## rowSums(apply(jklm,1,function(jklm)with(as.list(jklm),  c(part.b=phi[j,k]*phi[j,l]/I^2,part.c= 12/I^2/(I^2-1)*(k-(I+1)/2)*phi[j,k]*(m-(I+1)/2)*phi[l,m])    )))
        proj.B.sqr <- sum(sapply(B,function(m)sum(phi*m)^2))
        proj.C.sqr <- sum(C*phi)^2
        ## resid <- phi-proj
        ## phi.vec <- as.numeric(phi)
        ## proj.vec <- as.numeric(proj); resid.vec <- as.numeric(resid)
        ## c(proj=proj.vec%*%Q%*%proj.vec, resid=resid.vec%*%Q%*%resid.vec, cross=2*resid.vec%*%Q%*%proj.vec)
        c(proj.B.sqr=proj.B.sqr,proj.C.sqr=proj.C.sqr,demeaned=sum(phi^2)-proj.B.sqr ,resid.sqr=sum(phi^2)-proj.B.sqr-proj.C.sqr)
    })
})
parts <- simplify2array(parts)
with(as.data.frame(t(apply(parts,c(1,3),mean))),{
    op <- par(mfrow=c(1,2))
    plot(Is,proj.B.sqr,ylim=range(c(proj.B.sqr,proj.C.sqr)))
    points(Is,proj.C.sqr,col=2)
    plot(Is,proj.B.sqr/proj.C.sqr)
    par(op)
    ## plot(Is,resid.sqr);abline(lm(resid.sqr ~ Is)) #1
    ## plot(Is,demeaned,col=2)
})
## For general unstructured phi proj.B.sqr grows at a rate much faster
## than proj.C.sqr, but for phi with monotonically increasing rows,
## seems closer rates. proj.B.sqr is the same under both types of phi,
## being invariant to permutations of the rows, but proj.C.sqr is
## larger when phi is structured.  Appears proj.B.sqr always >
## proj.C.sqr:
summary(as.numeric(parts['proj.B.sqr',,] - parts['proj.C.sqr',,]))
## Line #1 shows that norm squared of residuals are growing linearly
## in I, can't hope to get O(sqrt(I)) bound for proj.Q.resid by
## bounding resid norm

## 11e Alternative dispersion measure. D(ax+b)=aD(x).
## ((I/2):1)%*%apply(cbind(x,rev(x))[1:(I/2),],1,diff)
D <- function(x){
    I <- length(x);
    x <- sort(x);
    wts <- -(1:I - (I+1)/2)[1:(I/2)]
    sqrt(12/I/(I^2-1)) *
        wts%*%apply(cbind(x,rev(x))[1:(I/2),,drop=FALSE],1,diff)
}
vars <- replicate(1e4, {
    x <- sort(runif(I))
    D.apprx <- 1/2*(I+2)/(I-1) * sum(apply(cbind(x,rev(x))[1:(I/2),,drop=FALSE],1,diff)^2)
c(D=D(x)^2,var=sum((x-mean(x))^2),D.apprx=D.apprx)
})
plot(vars['D',],vars['var',]); abline(0,1)
plot(vars['D.apprx',],vars['D',]); abline(0,1)
plot(vars['D.apprx',],vars['var',]); abline(0,1)

## growth rate of difference of dispersion measures
Is <- 3:50
spreads <- sapply(Is, function(I) {
    replicate(1e3, {
        x <- sort(runif(I))
        ## sd(x)*sqrt(I-1) - D(x)
        c(sd=sd(x)*sqrt(I-1),D=D(x))
    })
},simplify=FALSE)
spreads <- simplify2array(spreads)
boxplot(spreads['sd',,]-spreads['D',,])

plot(Is,apply(spreads['sd',,]-spreads['D',,],2,max)*sqrt(Is))
plot(Is,apply(spreads['sd',,]^2-spreads['D',,]^2,2,max))

plot(colMeans(spreads['sd',,]))
points(colMeans(spreads['D',,]),col=2)

plot(Is,apply(spreads['sd',,],2,max),ylim=c(0,2))
points(Is,apply(spreads['D',,],2,min),col=2)

## growth rate of difference of dispersion measures
Is <- 3:40
spreads <- sapply(Is, function(I) {
    replicate(1e2, {
        ## x <- sort(runif(I))
        phi <- t(replicate(I,sort(runif(I))))
        ## sd(x)*sqrt(I-1) - D(x)
        ## c(sd=sd(x)*sqrt(I-1),D=D(x))
        c(var=sum(apply(phi,1,var))*(I-1), D.sqr=sum(apply(phi,1,D))^2/I)
    })
},simplify=FALSE)
spreads <- simplify2array(spreads)
boxplot(spreads['var',,]-spreads['D.sqr',,])

plot(Is,apply(spreads['sd',,]-spreads['D',,],2,max)*sqrt(Is))
plot(Is,apply(spreads['sd',,]^2-spreads['D',,]^2,2,max))

## 3a_2 + a_1 = s fixed, a_1<=a_2
s <- 4
a1s <- seq(0,s/4,len=100)
plot(a1s,sapply(a1s,function(a1)var(c(-(s-a1)/6,-a1/2,a1/2,(s-a1)/6))))

a1 <- .3
x4 <- a2 <- 1/3*(s-a1)
x2s <- seq(0,a2-a1,len=100)
## plot(a1s,sapply(x2s,function(x2)D(c(0,x2,x2+a1,a2))))
plot(a1s,sapply(x2s,function(x2)var(c(0,x2,x2+a1,a2))))

I <- 8
wts <- -(1:I - (I+1)/2)[1:(I/2)]
D0 <- 1
as <- replicate(1e5, {
    a <- sort(runif(floor(I/2)))
    a <- a * D0 / c(rev(a)%*%wts)
})
xs <- apply(as,2,function(a) {
    ## a <- sort(runif(floor(I/2)))
    x <- cumsum(c(0,runif(floor(I/2)-1,0,rev(diff(a)))))
    x <- c(x,rev(x+rev(a)))
})
vars <- apply(xs,2,function(x)    var(x)*(I-1))
## hist(vars)
## a.star <- D0/sum(wts)
## var.star <- I*a.star^2/4
## x.star <- c(rep(0,floor(I/2)),rep(a.star,floor(I/2)))
## D(c(rep(0,I/2),rep(a.star,I/2)))
as[,which.max(vars)]
xs[,which.max(vars)]
max(vars)
## a.star <- D0/max(wts)
## var(c(0,a.star/2,a.star/2,a.star))*(I-1)
var(c(0,rep(D0/max(wts),I-1)))*(I-1)

I <- 4
x <- sort(runif(I))
diffs <- outer(x,x,'-')
diffs[row(diffs)<col(diffs)] <- 0
sum(diffs^2)/I #=var
## sum((x-mean(x))^2)
antidiag <- diffs[row(diffs)+col(diffs)==I+1] 
wts <- -(1:I - (I+1)/2)[1:(I/2)]
12/I/(I^2-1)*(wts%*%antidiag[1:length(wts)])^2 #D(x)^2

Is <- 3:30
## parts <- lapply(Is, function(I) {
replicate(5e2, {
    I <- sample(2:20,1)
    Q <- M.obj(I,symm=TRUE)
    B <- lapply(1:I,function(i)within(list(m=matrix(0,I,I)),m[i,] <- 1)[[1]] / sqrt(I))
    C <- with(list(m=matrix(1,nrow=I)%*%matrix(1:I,1)),m-rowMeans(m))
    C <- C/sqrt(sum(C^2))
    basis <- c(B,list(C))
    phi <- t(replicate(I,sort(runif(I))))
    proj.B.sqr <- sum(sapply(B,function(m)sum(phi*m)^2))
    proj.C.sqr <- sum(C*phi)^2
    resid.sqr <- sum(phi^2) - proj.B.sqr - proj.C.sqr
    ## phi.vec <- as.numeric(phi)
    ## proj.vec <- as.numeric(proj); resid.vec <- as.numeric(resid)
    ## c(proj=proj.vec%*%Q%*%proj.vec, resid=resid.vec%*%Q%*%resid.vec, cross=2*resid.vec%*%Q%*%proj.vec)
    ## sum(apply(phi,1,function(r)sum((r-mean(r))^2)))   -   12/I^2/(I^2-1) * sum(apply(phi,1,function(r)r%*%(1:I - (I+1)/2)))^2
    sum(apply(phi,1,function(r)sum((r-mean(r))^2)))   -   12/I^1/(I^2-1) * sum(apply(phi,1,function(r)r%*%(1:I - (I+1)/2))^2)
    
    c(proj.B=proj.B,proj.C.sqr=proj.C.sqr,sum((phi-rowMeans(phi))^2) - proj.C.sqr)#,sum(phi^2)-proj.B-proj.C.sqr)
})

,dd

## 11ee alternative dispersion measure as squared covariance of
## ordered x with (1,...,I) scaled
D <- function(x) {I <- length(x); cov(x,(1:I)/sd(1:I))^2 }


dd

## 11f 1. adding symmetries to the orthogonal Q-basis. Didn't work
## since no longer a basis. #2 After removing I(I-1)/2 dim subspace of
## nullspace of Q, |resid|^2 still growing with I.
phi.symmetries <- function(phi) {
    I <- nrow(phi)
    s4sub <- function(m)list(id=m,rot180=m[I:1,I:1],diag=t(m),antidiag=t(m[I:1,I:1]))
    shifts <- function(m)lapply(1:(I-1),function(s)(m[c((I-s+1):I,1:(I-s)),c((I-s+1):I,1:(I-s))]))
    symmetries <- c(s4sub(phi))#,s4sub(1-phi))
    symmetries <- lapply(symmetries,shifts)
    symmetries <- do.call(c,symmetries)
    symmetries
}


I <- sample(3:10,1)
Is <- 3:30
resids <- lapply(Is, function(I) {
    print(I)
    Q <- M.obj(I,symm=TRUE)
    B <- lapply(1:I,function(i)within(list(m=matrix(0,I,I)),m[i,] <- 1)[[1]] / sqrt(I))
    C <- with(list(m=matrix(1,nrow=I)%*%matrix(1:I,1)),m-rowMeans(m))
    ## sum(C^2) - I^2*(I^2-1)/12 # fla check
    C <- C/sqrt(sum(C^2))
    basis <- c(B,list(C))
    ## basis2 <- c(B,phi.symmetries(C),phi.symmetries(B[[2]])) #1
    ## basis2.mat <- sapply(basis2,function(m)as.numeric(t(m)))
    ## sum(abs(t(basis2.mat)%*%Q%*%basis2.mat))
    e <- Vectorize(function(j)unlist(within(list(v=rep(0,I^2)),v[j] <- 1))) #2
    flip <- function(v)as.numeric(matrix(v,nrow=I,byrow=TRUE))
    N0 <- t(sapply(1:(I^2),function(i)flip(e(i)))) - diag(I^2)
    N0 <- with(eigen(N0), vectors[,abs(values)>1e-7])
    replicate(1e2, {
        phi <- t(replicate(I,sort(runif(I))))
        proj <- Reduce('+',lapply(basis,function(m)sum(phi*m)*m))
        ## proj2 <- Reduce('+',lapply(basis2,function(m)sum(phi*m)*m))
        resid <- phi-proj
        ## resid <- phi-proj
        phi.vec <- as.numeric(phi)
        proj.vec <- as.numeric(proj); resid.vec <- as.numeric(resid)
        resid2 <- resid.vec - N0%*%(t(N0)%*%resid.vec)
        ## c(residQresid=resid.vec%*%Q%*%resid.vec)
        c(resid.sqr=sum(resid^2),resid2.sqr=sum(resid2^2))
    })
})
resids <- simplify2array(resids)

matplot(t(resids['resid.sqr',,]),pch='.',cex=5,col=1)
matplot(t(resids['resid2.sqr',,]),pch='.',cex=5,col=2,add=TRUE) #2


## 11g residQresid small due to null space of Q or null vectors?
source('misc.R')
require(parallel)
Is <- 3:30
## parts <- lapply(Is, function(I) {
parts <- mclapply(Is, mc.cores=detectCores()-3,FUN=function(I) {
    print(I)
    Q <- M.obj(I,symm=TRUE)
    Q2 <- Q%*%Q
    lambda <- sqrt((I-2)/2/I^2)
    B <- lapply(1:I,function(i)within(list(m=matrix(0,I,I)),m[i,] <- 1)[[1]] / sqrt(I))
    C <- with(list(m=matrix(1,nrow=I)%*%matrix(1:I,1)),m-rowMeans(m))
    C <- C/sqrt(sum(C^2))
    basis <- c(B,list(C))
    ## ww <- t(sapply(1:(I^2),function(i)flip(e(i))))
    ## N0 <- ww - diag(I^2)
    pos.vecs <- with(eigen(Q),vectors[,values>1e-7])
    neg.vecs <- with(eigen(Q),vectors[,values< -1e-7])
    ## phi.vec <- pos.vecs[,1]
    ## phi <- matrix(phi.vec,I,byrow=TRUE)
    ## proj <- Reduce('+',lapply(basis,function(m)sum(phi*m)*m))
    ## resid <- phi-proj
    ## phi.vec <- as.numeric(phi)
    ## proj.vec <- as.numeric(proj); resid.vec <- as.numeric(resid)
    ## c(proj=proj.vec%*%Q%*%proj.vec, resid=resid.vec%*%Q%*%resid.vec, cross=2*resid.vec%*%Q%*%proj.vec)
    ## phi.vec%*%Q%*%phi.vec
    replicate(1e3, {
        phi <- t(replicate(I,sort(runif(I))))
        ## phi <- matrix(runif(I^2),I)
        ## phi <- matrix(rbinom(I^2,1,.5),I)
        phi <- matrix(c(rep(0,I-1),1),nrow=I,ncol=I,byrow=TRUE) ##a
        proj <- Reduce('+',lapply(basis,function(m)sum(phi*m)*m))
        resid <- as.numeric(t(phi-proj))
        ## resid.approx <- resid - 
        ## phi.vec <- as.numeric(phi)
        ## proj.vec <- as.numeric(proj); resid.vec <- as.numeric(resid)
        ## c(proj=proj.vec%*%Q%*%proj.vec, resid=resid.vec%*%Q%*%resid.vec, cross=2*resid.vec%*%Q%*%proj.vec,resid.norm2=sum(resid.vec^2),proj.norm2=sum(proj.vec^2))
        c(resid.norm2=sum(resid^2),Q.resid.norm2=sum(((Q/lambda)%*%resid)^2),resid.Q.resid.norm2=resid%*%(Q/lambda)%*%resid,proj.pos2=sum((t(pos.vecs)%*%resid)^2),proj.neg2=sum((t(neg.vecs)%*%resid)^2),Q2.resid.norm2=sum(((Q2/lambda^2)%*%resid)^2))#,proj.demeaned.norm2=sum((proj-mean(proj))^2))
    })
})
parts <- simplify2array(parts)

## matplot(Is,t(parts['resid.norm2',,])/Is^2,pch='.',col=rgb(1,0,0,.2),cex=5) #5
matplot(t(parts['resid.norm2',,]),pch='.',col=rgb(1,0,0,.2),cex=5,ylim=range(parts),x=Is)
## matplot(t(parts['Q.resid.norm2',,]),pch='.',col=rgb(0,1,0,.2),cex=5,add=TRUE) #1
## matplot(t(parts['resid.Q.resid.norm2',,]),pch='.',col=rgb(0,0,1,.2),cex=5,add=TRUE)
## matplot(t(parts['proj.pos2',,]),pch='.',col=rgb(0,1,0,.2),cex=5,add=TRUE) #2
## matplot(t(parts['proj.neg2',,]),pch='.',col=rgb(0,0,1,.2),cex=5,add=TRUE)
matplot(t(parts['Q2.resid.norm2',,]),pch='.',col=rgb(0,0,1,.2),cex=5) #3
## matplot(t(parts['proj.demeaned.norm2',,]),pch='.',col=rgb(0,0,1,.2),cex=5) #4

## 1. projection of resid onto col space of Q is O(1) so proj.Q.resid
## is O(I^(1/2)). 2. Projection of resid onto pos or neg parts is also
## O(1). So don't need to look at the difference. Can use convexity of
## positive and negative parts of Q. 3. Equivalently, projection of
## resid onto col space of Q^2 is O(1).  4. O(I) growth for proj norm
## even after demeaning (checking since ones vector is orthogonal to
## Q.rho) [update: all wrong. counterex phi in line ##a. further
## update: actually this looks like it asymptotes, so O(1).]. 


matplot(t(parts['Q.resid.norm2',,]),pch='.',col=rgb(0,1,0,.2),cex=5,add=TRUE)


with(as.data.frame(t(apply(parts,c(1,3),mean))), {
    matplot(t(resid.norm2),pch='.',col=rgb(1,0,0,.2),cex=5,ylim=range(c(resid,cross)),x=Is)
    ## matplot(t(resid),pch='.',cex=5,col=rgb(0,1,0,.2),add=TRUE)
})


## 11h checking formulas
source('misc.R')
I <- 5
Q <- M.obj(I,symm=TRUE)
Q2 <- Q%*%Q
lambda <- sqrt((I-2)/2/I^2)
B <- lapply(1:I,function(i)within(list(m=matrix(0,I,I)),m[i,] <- 1)[[1]] / sqrt(I))
C <- with(list(m=matrix(1,nrow=I)%*%matrix(1:I,1)),m-rowMeans(m))
C <- C/sqrt(sum(C^2))
basis <- c(B,list(C))
phi <- t(replicate(I,sort(runif(I))))
proj <- Reduce('+',lapply(basis,function(m)sum(phi*m)*m))
rho <- resid <- as.numeric(t(phi-proj))
phi.vec <- as.numeric(t(phi))
proj.vec <- as.numeric(t(proj))
alpha <- sum(phi*C)*sqrt(12/I^2/(I^2-1))
ps <- 1:I
idx <- function(v,p,r,I)sapply(p,function(p0)sapply(r,function(r0)v[I*(p0-1)+r0]))
q <- sample(1:I,1);s <- sample(1:I,1)
## ## ## 11h-0 formula for generic elt of residual
## sum(abs(sapply(1:I,function(q)sapply(1:I,function(s) {
##     idx(rho,q,s,I)  -  (phi[q,s] - mean(phi[q,]) - alpha*(s-(I+1)/2))
## }))))
## sum(rho^2)
## (I-1)*sum(apply(phi,1,var)) - I*(I-1)*cov(colMeans(phi),ps/sd(ps))^2
## 11h-1 formula for generic elt of Q.v
## v <- matrix(runif(I^2),I)
## v.vec <- as.numeric(t(v))
## p <- sample(1:I,1);r <- sample(1:I,1)
## sum(abs(sapply(1:I,function(p)sapply(1:I,function(r) {
##     idx(Q%*%v.vec,p,r,I) - (1/I^2*(sum(v[p,])+sum(v[r,])+sum(v[,r])+sum(v[,p]) + sum(diag(v)) + (p==r)*sum(v)) - 1/2/I*((p==r)*(sum(v[p,])+sum(v[,p])) + v[p,p] + v[r,r]) - 4/I^3*sum(v))
##     }))))
## 11h-2 formula for generic elt of Q.resid
## p <- sample(1:I,1);r <- sample(1:I,1)
## sum(abs(sapply(1:I,function(p)sapply(1:I,function(r) {
## idx(Q%*%rho,p,r,I) - (
## 1/I^2*(sum(phi[,r])+sum(phi[,p])+sum(diag(phi))-3/I*sum(phi)-I*alpha*(p+r-(I+1))) - 1/2/I*( phi[p,p]+phi[r,r]-mean(phi[p,])-mean(phi[r,])-alpha*(p+r-(I+1))+(p==r)*(sum(phi[,p])-1/I*sum(phi)-I*alpha*(p-(I+1)/2)))
## )}))))
## ## 11h-3 formula for generic elt of proj
## p <- sample(1:I,1);r <- sample(1:I,1)
## sum(abs(sapply(1:I,function(p)sapply(1:I,function(r) {
## idx(proj.vec,p,r,I)-
## (rowMeans(phi)[p] + cov(colMeans(phi),ps/sd(ps))*(r-(I+1)/2)/sd(ps))
## }))))
## ## 11h-4 formula for generic elt of proj * Q.resid
## p <- sample(1:I,1);r <- sample(1:I,1)
## ps <- 1:I
## sum(abs(sapply(1:I,function(p)sapply(1:I,function(r) {
## idx(proj.vec*Q%*%rho,p,r,I) - (
## (rowMeans(phi)[p]+cov(colMeans(phi),ps/sd(ps))*(r-(I+1)/2)/sd(ps)) * 
## (1/I^2*(sum(phi[,r])+sum(phi[,p])+sum(diag(phi))-3/I*sum(phi)-I*alpha*(p+r-(I+1))) - 1/2/I*( phi[p,p]+phi[r,r]-mean(phi[p,])-mean(phi[r,])-alpha*(p+r-(I+1))+(p==r)*(sum(phi[,p])-1/I*sum(phi)-I*alpha*(p-(I+1)/2))))
## )}))))
##  ## 11h-5 {p==r} entries of proj * Q.resid
## sum(sapply(1:I, function(p)(rowMeans(phi)[p]+cov(colMeans(phi),ps/sd(ps))*(p-(I+1)/2)/sd(ps)) * (colMeans(phi)[p]-mean(phi)-(p-(I+1)/2)*cov(colMeans(phi),ps)/var(ps))))
## sum((rowMeans(phi) - mean(phi) + cov(colMeans(phi),ps/sd(ps))*(ps-(I+1)/2)/sd(ps))*(colMeans(phi) - mean(phi) - cov(colMeans(phi),ps/sd(ps))*(ps-(I+1)/2)/sd(ps)))
## sum((rowMeans(phi) - mean(phi)) * (colMeans(phi) - mean(phi))) - sum(cov(colMeans(phi),ps/sd(ps))^2*(I-1)) + sum(cov(colMeans(phi),ps/sd(ps))*(ps-(I+1)/2)/sd(ps)*(colMeans(phi) - mean(phi))) - sum(cov(rowMeans(phi),ps/sd(ps))*(ps-(I+1)/2)/sd(ps)*(colMeans(phi) - mean(phi)))
## sum((rowMeans(phi)-mean(phi))*(colMeans(phi)-mean(phi))) - sum(cov(colMeans(phi),ps/sd(ps))^2*(I-1)) + sum(((colMeans(phi)-mean(phi))*(ps-(I+1)/2)/sd(ps) - (rowMeans(phi)-mean(phi))*(ps-(I+1)/2)/sd(ps))*cov(colMeans(phi),ps/sd(ps)))
## (I-1)*cov(rowMeans(phi),colMeans(phi)) - cov(colMeans(phi),ps/sd(ps))^2*(I-1)+ sum(((colMeans(phi)-mean(phi))*(ps-(I+1)/2)/sd(ps) - (rowMeans(phi)-mean(phi))*(ps-(I+1)/2)/sd(ps))*cov(colMeans(phi),ps/sd(ps)))
## (I-1)*(cov(rowMeans(phi),colMeans(phi)) - cov(colMeans(phi),ps/sd(ps))^2) + sum((colMeans(phi)-rowMeans(phi))*(ps-(I+1)/2)/sd(ps))*cov(colMeans(phi),ps/sd(ps))
## (I-1)*(cov(rowMeans(phi),colMeans(phi)) - cov(colMeans(phi),ps/sd(ps))^2) +(I-1)*cov(colMeans(phi)-rowMeans(phi),ps/sd(ps))*cov(colMeans(phi),ps/sd(ps))
## (I-1)*(cov(rowMeans(phi),colMeans(phi)) - cov(rowMeans(phi),ps/sd(ps))* cov(colMeans(phi),ps/sd(ps)))
## ## 11h-6 formula for generic elt of resid * Q.resid
## p <- sample(1:I,1);r <- sample(1:I,1)
## sum(abs(sapply(1:I,function(p)sapply(1:I,function(r) {
## idx(rho * Q%*%rho,p,r,I) - (
## 1/I^2*(sum(phi[,r])+sum(phi[,p])+sum(diag(phi))-3/I*sum(phi)-I*alpha*(p+r-(I+1))) - 1/2/I*( phi[p,p]+phi[r,r]-mean(phi[p,])-mean(phi[r,])-alpha*(p+r-(I+1))+(p==r)*(sum(phi[,p])-1/I*sum(phi)-I*alpha*(p-(I+1)/2)))
## ) * (phi[p,r]-rowMeans(phi)[p] - alpha*(r-(I+1)/2)     )   }))))
## ## separating out O(1) terms from {p==r} terms
## sum(abs(sapply(1:I,function(p)sapply(1:I,function(r) {
##     idx(rho * Q%*%rho,p,r,I) - (
##         (
##     1/I^2*(sum(phi[,r])+sum(phi[,p])+sum(diag(phi))-3/I*sum(phi)-I*alpha*(p+r-(I+1))) -
##     1/2/I*( phi[p,p]+phi[r,r]-mean(phi[p,])-mean(phi[r,])-alpha*(p+r-(I+1)))
##             )* (phi[p,r]-rowMeans(phi)[p] - alpha*(r-(I+1)/2)) -
##     1/2*(p==r)*(mean(phi[,p])-mean(phi)-alpha*(p-(I+1)/2)) * (phi[p,r]-rowMeans(phi)[p] - alpha*(r-(I+1)/2))
## )    }))))
## ## C.pre.Q.resid
## ps <- 1:I
## C.pre <- as.numeric(replicate(I,ps))
## C.pre%*%Q%*%resid
## sum(sapply(1:I,function(p)sapply(1:I,function(r) {
##     r*(
##         1/I*(colMeans(phi)[p]+colMeans(phi)[r]+1/2*(rowMeans(phi)[p]+rowMeans(phi)[r])) + sum(diag(phi))/I^2 - 3/I^3*sum(phi) - alpha/2/I*(p+r-I-1) - 1/2/I*(diag(phi)[p]+diag(phi)[r] + (p==r)*(colSums(phi)[p] - sum(phi)/I - I*alpha*(p-(I+1)/2)))
##         )})))
## 3/2*mean(phi)*I*(I+1)/2 + sum(ps*colMeans(phi)) + 1/2*sum(ps*rowMeans(phi)) + mean(diag(phi))*1/2*I*(I+1)-3*mean(phi)/2*I*(I+1)-alpha/2/I*(1/4*I^2*(I+1)^2+1/6*I^2*(I+1)*(2*I+1)-I^2*(I+1)^2/2) - 1/2/I*(I*sum(ps*diag(phi))+sum(diag(phi))*I*(I+1)/2 + sum(ps*colSums(phi)) - sum(phi)*(I+1)/2 - I*alpha*(1/6*I*(I+1)*(2*I+1)-1/4*I*(I+1)^2))
## 1/2*sum(ps*(colMeans(phi)+rowMeans(phi))) + (I+1)/4*sum(diag(phi)) - alpha/2*I*(I^2-1)/12 - 1/2*sum(ps*diag(phi)) -1/2*(I+1)/I*sum(phi) + alpha/4*I*(I^2-1)/6
## 1/2*sum(ps*(colMeans(phi)+rowMeans(phi))) -1/2*(I+1)/I*sum(phi)+ (I+1)/4*sum(diag(phi))  - 1/2*sum(ps*diag(phi)) 
## (I-1)/2*cov(ps,colMeans(phi)+rowMeans(phi)-diag(phi))
## cor(ps,colMeans(phi)+rowMeans(phi)-diag(phi))
## sd(colMeans(phi)+rowMeans(phi)-diag(phi))
## c(t(C))%*%Q%*%resid*sqrt(sum((C.pre-mean(C.pre[1:I]))^2))
## c(t(C))%*%Q%*%resid*sqrt(I*(I-1)*var(ps))
## sum((C.pre-mean(C.pre[1:I]))^2)
## I*sum((C.pre[1:I]-mean(C.pre))^2)
## I*(I-1)*var(ps)
## I^2*(I^2-1)/12
## (C.pre - (I+1)/2)/sqrt(I*(I-1)*var(ps)) - t(C)
## c(t(C))%*%Q%*%rho
## cov(ps,colMeans(phi)+rowMeans(phi)-diag(phi)) * sqrt(3)/I*sqrt((I-1)/(I+1))
## 1/2*sqrt((I-1)/I)*cov(ps/sd(ps),colMeans(phi)+rowMeans(phi)-diag(phi))
## sum(phi*C)
## (I-1)*cov(colSums(phi),ps) / sqrt(I*(I-1)*var(ps))
## sum(phi*C) * c(t(C))%*%Q%*%rho
## (I-1)/2*cov(colMeans(phi),ps/sd(ps))*cov(ps/sd(ps),colMeans(phi)+rowMeans(phi)-diag(phi))
## ## (B,Q.rho)
## j <- 2
## sum(sapply(1:I,function(p)sapply(1:I,function(r) {
##     (p==j)*(
##         1/I*(colMeans(phi)[p]+colMeans(phi)[r]+1/2*(rowMeans(phi)[p]+rowMeans(phi)[r])) + sum(diag(phi))/I^2 - 3/I^3*sum(phi) - alpha/2/I*(p+r-I-1) - 1/2/I*(diag(phi)[p]+diag(phi)[r] + (p==r)*(colSums(phi)[p] - sum(phi)/I - I*alpha*(p-(I+1)/2)))
##     )
## })))
## 1/2*(rowMeans(phi)[j] + colMeans(phi)[j]) - mean(phi) + 1/2*mean(diag(phi)) - phi[j,j]/2
## sum(phi*B[[j]])*sum(c(t(B[[j]]))%*%Q%*%rho)
## rowMeans(phi)[j] * ((rowMeans(phi)[j]+colMeans(phi)[j])/2 - mean(phi) + mean(diag(phi))/2 - phi[j,j]/2)
## sum(sapply(1:I, function(j)sum(phi*B[[j]])*sum(c(t(B[[j]]))%*%Q%*%rho)))
## (I-1)/2*cov(rowMeans(phi),rowMeans(phi)+colMeans(phi)-diag(phi)) 

ps <- 1:I
phi0 <- matrix(0,I,I)
for(j in (I/2+1):I)for(k in j:I) phi0[j,k] <- phi0[k,j] <- 1
diag(phi0 ) <- 0
cov(colMeans(phi0)+rowMeans(phi0)-diag(phi0),ps/sd(ps))

dd


## growth rate of (phi,C)(C,Q.rho) on adversarial phi
source('misc.R')
require(parallel)
Is <- 3:40
parts <- mclapply(Is, mc.cores=detectCores()-3,FUN=function(I) {
## parts <- lapply(Is,FUN=function(I) {
    print(I)
    Q <- M.obj(I,symm=TRUE)
    B <- lapply(1:I,function(i)within(list(m=matrix(0,I,I)),m[i,] <- 1)[[1]] / sqrt(I))
    C <- with(list(m=matrix(1,nrow=I)%*%matrix(1:I,1)),m-rowMeans(m))
    C <- C/sqrt(sum(C^2))
    basis <- c(B,list(C))
    ## phi <- matrix(0,I,I)
    ## for(j in ceiling(I/2):I)for(k in j:I) phi[j,k] <- phi[k,j] <- 1
    phi <- matrix(1,I,I)
    for(j in 1:(I/2))for(k in 1:(I/2)) phi[j,k] <- phi[k,j] <- 0
    ## diag(phi) <- 0
    proj <- Reduce('+',lapply(basis,function(m)sum(phi*m)*m))
    rho.vec <- resid.vec <- as.numeric(t(phi-proj))
    c(factor1=sum(phi*C), factor2= c(t(C))%*%Q%*%rho.vec, full=c(t(proj))%*%Q%*%rho.vec)
})
parts <- simplify2array(parts)
plot(Is,parts['factor1',],ylim=range(parts))
points(Is,parts['factor2',],col=2)
points(Is,parts['factor1',]*parts['factor2',],col=3) # 1
plot(Is,parts['full',]) # 2
plot(Is,parts['factor2',]) 
## can't find a valid phi with linear rate. always seems O(1).

I <- 7
while(TRUE) {
    as <- sample(0:I,I,replace=TRUE)
    phi <- t(sapply(as,function(a)c(rep(0,a),rep(1,I-a))))
    diffs <- diff(colMeans(phi)+rowMeans(phi)-diag(phi))
    if(sum(diffs<0)==0) {
        print(as)
        ## print(diffs)
        print(max(diffs))
        print(phi)
        ## break
    }
}



I <- 7
Is <- 3:30
max.sds <- sapply(Is, function(I) {
    print(I)
max(replicate(2e4, {
    as <- sample(0:I,I,replace=TRUE)
    phi <- t(sapply(as,function(a)c(rep(0,a),rep(1,I-a))))
    sd(colMeans(phi)+rowMeans(phi)-diag(phi))
    }))
} )
plot(Is,max.sds)

## I-as == rowSums(phi), I*ecdf(as)[0:(I-1)] == colSums(phi), as<(1:I)
## == diag(phi)
sum(I-as != rowSums(phi))
sum(I*ecdf(as)(0:(I-1)) != colSums(phi))
sum((as<(1:I)) != diag(phi))
colMeans(phi)+rowMeans(phi)-diag(phi)
ecdf(as)(0:(I-1)) + 1-as/I - (as<(1:I))
var(colMeans(phi)+rowMeans(phi)-diag(phi))
var(ecdf(as)(0:(I-1))  - as/I - (as<(1:I)))
var(ecdf(as)(0:(I-1)))+ var(as/I + (as<(1:I))) - 2*cov(ecdf(as)(0:(I-1)), as/I + (as<(1:I)))

I  <- 13
as <- sample(0:I,I,replace=TRUE)
parts <- replicate(1e3, {
idx <- sample(1:I)
with(list(as=as[idx]),c(var(ecdf(as)(0:(I-1))), var(as/I + (as<(1:I))), - 2*cov(ecdf(as)(0:(I-1)), as/I + (as<(1:I)))))
})

boxplot(t(parts))
plot(parts[2,],parts[3,])
summary(parts[2,]+parts[3,])
lm(parts[2,]~parts[3,])

dd


## ## 11i offdiag elts of rho * Q.rho
## source('misc.R')
## I <- 5
## Q <- M.obj(I,symm=TRUE)
## ## Q2 <- Q%*%Q
## lambda <- sqrt((I-2)/2/I^2)
## B <- lapply(1:I,function(i)within(list(m=matrix(0,I,I)),m[i,] <- 1)[[1]] / sqrt(I))
## C <- with(list(m=matrix(1,nrow=I)%*%matrix(1:I,1)),m-rowMeans(m))
## C <- C/sqrt(sum(C^2))
## basis <- c(B,list(C))
## phi <- t(replicate(I,sort(runif(I))))
## proj <- Reduce('+',lapply(basis,function(m)sum(phi*m)*m))
## rho <- resid <- as.numeric(t(phi-proj))
## phi.vec <- as.numeric(t(phi))
## proj.vec <- as.numeric(t(proj))
## ## with(list(m=rho*(Q%*%rho)),sum(m[row(m)!=col(m)]))


## 11i growth rate of diag/offdiag parts of cross term proj*Q.rho, rho*Q.rho
source('misc.R')
require(parallel)
Is <- 3:40
parts <- mclapply(Is, mc.cores=detectCores()-3,FUN=function(I) {
## parts <- lapply(Is,FUN=function(I) {
    print(I)
    Q <- M.obj(I,symm=TRUE)
    B <- lapply(1:I,function(i)within(list(m=matrix(0,I,I)),m[i,] <- 1)[[1]] / sqrt(I))
    C <- with(list(m=matrix(1,nrow=I)%*%matrix(1:I,1)),m-rowMeans(m))
    C <- C/sqrt(sum(C^2))
    basis <- c(B,list(C))
    replicate(5e1, {
        phi <- t(replicate(I,sort(runif(I))))
        proj <- Reduce('+',lapply(basis,function(m)sum(phi*m)*m))
        rho.vec <- resid.vec <- as.numeric(t(phi-proj))
        phi.vec <- as.numeric(t(phi))
        proj.vec <- as.numeric(t(proj))
        chi.vec <- proj.vec * (Q%*%rho.vec)
        chi.mat <- matrix(t(chi.vec),nrow=I,byrow=TRUE)
        rho.mat <- matrix(t(rho.vec),nrow=I,byrow=TRUE)
        c(cross.diag=sum(diag(chi.mat)),cross.offdiag=sum( chi.mat[row(chi.mat)!=col(chi.mat)]),rhoQrho.diag=sum(diag(rho.mat)),rhoQrho.offdiag=sum(rho.mat[row(rho.mat)!=col(rho.mat)]))
        ## sum(chi.vec)
    })
})
parts <- simplify2array(parts)

matplot(Is,t(parts['cross.offdiag',,]),pch='.',cex=4,col=1)
matplot(Is,t(parts['cross.diag',,]),pch='.',cex=4,col=2,add=TRUE)
## matplot(Is,t(parts['offdiag',,]+parts['diag',,]),pch='.',cex=4,col=1)

matplot(Is,t(parts['rhoQrho.offdiag',,]),pch='.',cex=4,col=1)
matplot(Is,t(parts['rhoQrho.diag',,]),pch='.',cex=4,col=2,add=TRUE)

## 1. seems sum of offdiag terms of crossterm is much larger than diag
## terms. 2. the two sums are linearly related.
plot(parts['cross.diag',,20],parts['cross.offdiag',,20])
plot(parts['cross.diag',,3],parts['cross.offdiag',,3])
## resid(lm(parts[1,,20]~parts[2,,20]))


## Q.rho and proj almost orthogonal, like rho and proj. why?
source('misc.R')
require(parallel)
Is <- 3:40
parts <- mclapply(Is, mc.cores=detectCores()-3,FUN=function(I) {
## parts <- lapply(Is,FUN=function(I) {
    print(I)
    Q <- M.obj(I,symm=TRUE)
    B <- lapply(1:I,function(i)within(list(m=matrix(0,I,I)),m[i,] <- 1)[[1]] / sqrt(I))
    C <- with(list(m=matrix(1,nrow=I)%*%matrix(1:I,1)),m-rowMeans(m))
    C <- C/sqrt(sum(C^2))
    basis <- c(B,list(C))
    replicate(5e1, {
        phi <- t(replicate(I,sort(runif(I))))
        proj <- Reduce('+',lapply(basis,function(m)sum(phi*m)*m))
        rho.vec <- resid.vec <- as.numeric(t(phi-proj))
        phi.vec <- as.numeric(t(phi))
        proj.vec <- as.numeric(t(proj))
        chi.vec <- proj.vec * (Q%*%rho.vec)
        chi.mat <- matrix(t(chi.vec),nrow=I,byrow=TRUE)
        rho.mat <- matrix(t(rho.vec),nrow=I,byrow=TRUE)
        Q.rho.mat <- matrix(t(Q%*%rho.vec),nrow=I,byrow=TRUE)
        Q.rho.mat.nodiag <- Q.rho.mat;  diag(Q.rho.mat.nodiag) <- 0
        Q.rho.mat.onlydiag <- diag(diag(Q.rho.mat))
        c(proj.norm2=sum(proj^2),Q.rho.norm2=sum((Q%*%rho.vec)^2),proj.Q.norm2=sum((Q%*%proj.vec)^2),cross=proj.vec%*%c(Q.rho.mat),cross.nodiag=proj.vec%*%c(Q.rho.mat.nodiag),cross.onlydiag=proj.vec%*%c(Q.rho.mat.onlydiag),C.Q.rho=c(C)%*%c(Q.rho.mat/sqrt(sum(Q.rho.mat^2))))
        ## sum(chi.vec)
    })
})
parts <- simplify2array(parts)

matplot(Is,t(parts['C.Q.rho',,])*sqrt(Is),pch='.',cex=4,col=2) #1

matplot(Is,t(parts['proj.norm2',,]),pch='.',cex=4,col=1)
matplot(Is,t(parts['Q.rho.norm2',,]),pch='.',cex=4,col=1,add=TRUE)
matplot(Is,t(parts['cross',,]),pch='.',cex=4,col=2,add=TRUE)

matplot(Is,t(parts['cross.nodiag',,]),pch='.',cex=4,col=1)
matplot(Is,t(parts['cross.onlydiag',,]),pch='.',cex=4,col=2,add=TRUE)

matplot(Is,t(parts['proj.Q.norm2',,]),pch='.',cex=4,col=1)
## line #1 suggests angle between C and Q.rho is O(1/sqrt(I))


## growth rate of |Q.rho|^2, seems O(1/I) for uniform dist
## phi. Counterexample below more like O(sqrt(I)) rate [update: no,
## seems to asymptote, so O(1).
source('misc.R')
require(parallel)
Is <- 3:60
idx <- function(v,p,r,I)sapply(p,function(p0)sapply(r,function(r0)v[I*(p0-1)+r0]))
## parts <- lapply(Is, function(I) {
parts <- mclapply(Is, mc.cores=detectCores()-3,FUN=function(I) {
    print(I)
    C <- with(list(m=matrix(1,nrow=I)%*%matrix(1:I,1)),m-rowMeans(m))
    C <- C/sqrt(sum(C^2))
    replicate(5e1, {
        phi <- t(replicate(I,sort(runif(I))))
        phi <- matrix(c(rep(0,I-1),1),nrow=I,ncol=I,byrow=TRUE)
        ## alpha <- sum(phi*C)*sqrt(12/I^2/(I^2-1))
        ## Q.rho <- sapply(1:I,function(p)sapply(1:I,function(r) 1/I^2*(sum(phi[,r])+sum(phi[,p])+sum(diag(phi))-3/I*sum(phi)-I*alpha*(p+r-(I+1))) - 1/2/I*( phi[p,p]+phi[r,r]-mean(phi[p,])-mean(phi[r,])-alpha*(p+r-(I+1))+(p==r)*(sum(phi[,p])-1/I*sum(phi)-I*alpha*(p-(I+1)/2)))))
        ## sum(Q.rho^2)
        ## sum(diag(Q.rho)^2)
        ## Q.rho.diag <- sapply(1:I,function(p)1/2/I*(sum(phi[,p])-1/I*sum(phi)-I*alpha*(p-(I+1)/2)))
        ## sum(Q.rho.diag^2)  *4
        (I-1)*var(colMeans(phi)) -12/I*(I-1)/(I+1)*cov(colMeans(phi),1:I)^2
    })
})
parts <- simplify2array(parts)
matplot(Is,t(parts),pch='.',cex=4,col=1)

## cross term is the difference of innerproduct (cov) and a semidefinite form
I <- 5
Is <- 1:30
by.I <- sapply(Is, function(I) {
    ps <- 1:I
    diffs <- replicate(1e3, {
        x <- runif(I); y <- runif(I)
        cov(x,y) - cov(x,ps/sd(ps)) * cov(y,ps/sd(ps))
    })
})
matplot(t(by.I)*sqrt(Is),pch='.',cex=5,col=1)
## For uniform x,y, seems of order sqrt(I). Maybe not for adversarial
## x,y. For adversarial x,y may need to use structure of
## x=colMeans(phi),y=rowMeans(phi)

hist(diffs)

## checking bounds from wolfram for cov(x,y) subject to cov(x,p)=a,cov(y,p)=b,|x|^2=|y^2|=1
## I=3 case
I <- 3
ps <- 1:I / sd(1:I)
diffs <- replicate(1e3, {
    x <- rnorm(I); x <- (x-mean(x))/sqrt(sum((x-mean(x))^2))
    y <- rnorm(I); y <- (y-mean(y))/sqrt(sum((y-mean(y))^2))
    a <- x%*%ps; b <- y%*%ps
    bounds <- c(lower=1/2*(a*b-sqrt(2-a^2)*sqrt(2-b^2)),upper=1/2*(a*b+sqrt(2-a^2)*sqrt(2-b^2)))
    bounds - c(x%*%y)
})
max(diffs[1,])
min(diffs[2,])
summary(diffs[1,])
summary(diffs[2,])


## general I case
I <- 6
ps <- 1:I / sd(1:I)
diffs <- replicate(1e3, {
    x <- rnorm(I); x <- (x-mean(x))/sd(x)/sqrt(I-1)#sqrt(sum((x-mean(x))^2))
    y <- rnorm(I); y <- (y-mean(y))/sqrt(sum((y-mean(y))^2))
    a <- x%*%ps; b <- y%*%ps
    bounds <- c(lower=1/(I-1)*(a*b-sqrt((I-1)-a^2)*sqrt((I-1)-b^2)),upper=1/(I-1)*(a*b+sqrt((I-1)-a^2)*sqrt((I-1)-b^2)))
    bounds - c(x%*%y)
})
max(diffs[1,])
min(diffs[2,])
summary(diffs[1,])
summary(diffs[2,])

## checking jensens gap, but needed it to be O(1/I). looks O(1).
I <- 5
Is <- 2:30
diffs <- sapply(Is, function(i) {
    replicate(1e2, {
        phi <- t(replicate(I,sort(runif(I))))
        mean(apply(phi,1,var)) - var(colMeans(phi))
    })
})
matplot(Is,t(diffs),pch='.',col=1,cex=5)


phi <- matrix(runif(9),nrow=3)




## 12. fuller Q-orthogonal basis. Includes constant row, constant
## column, and constant diagonal vectors.
source('misc.R')
I <- 6
Q <- M.obj(I,symm=TRUE)
phi <- matrix(runif(I^2),I)
ones <- matrix(1,nrow=I)
proj <- kronecker(t(ones),rowMeans(phi)) + t(kronecker(t(ones),colMeans(phi))) - kronecker(ones,t(ones))*mean(phi) + (sum(diag(phi))-I*mean(phi))/(I-1)*(diag(I)-1/I)
proj.vec <- c(t(proj))
proj.vec%*%Q%*%proj.vec


## 12a checking formulas
source('misc.R')
I <- 6
Q <- M.obj(I,symm=TRUE)
phi <- matrix(runif(I^2),I)
ones <- matrix(1,nrow=I)
proj <- kronecker(t(ones),rowMeans(phi)) + t(kronecker(t(ones),colMeans(phi))) - kronecker(ones,t(ones))*mean(phi) + (sum(diag(phi))-I*mean(phi))/(I-1)*(diag(I)-1/I)
rho <- resid <- as.numeric(t(phi-proj))
phi.vec <- c(t(phi))
proj.vec <- c(t(proj))
idx <- function(v,p,r,I)sapply(p,function(p0)sapply(r,function(r0)v[I*(p0-1)+r0]))
## ## 12a-1 generic elt of projection
## i <- sample(1:I,1);j <- sample(1:I,1)
## sum(abs(sapply(1:I,function(i)sapply(1:I,function(j) {
##     idx(proj.vec,i,j,I)  -  (rowMeans(phi)[i]+colMeans(phi)[j]-mean(phi)+(sum(diag(phi))-I*mean(phi))/(I-1)*((i==j)-1/I))
## }))))
## 12a-2 generic elt of residual
sum(abs(sapply(1:I,function(i)sapply(1:I,function(j) {
    idx(rho,i,j,I)  -  (phi[i,j] - rowMeans(phi)[i]-colMeans(phi)[j]+mean(phi)-(sum(diag(phi))-I*mean(phi))/(I-1)*((i==j)-1/I))
}))))
##  12a-3 norm^2 of projection
## sum(proj.vec^2)
## I^2*mean(phi)^2 + (I-1)/I*(var(colSums(phi))+var(rowSums(phi))) + ((sum(diag(phi))-I*mean(phi))/(I-1))^2*(I-1)
## sum((kronecker(ones,t(ones))*mean(phi))^2)  + sum((kronecker(t(ones),rowMeans(phi)) - mean(phi))^2) + sum((t(kronecker(t(ones),colMeans(phi))) - mean(phi))^2) + sum(((sum(diag(phi))-I*mean(phi))/(I-1)*(diag(I)-1/I))^2)
## ## kronecker(ones,t(ones))*mean(phi) + (kronecker(t(ones),rowMeans(phi)) - mean(phi)) + (t(kronecker(t(ones),colMeans(phi))) - mean(phi)) + ((sum(diag(phi))-I*mean(phi))/(I-1)*(diag(I)-1/I))
## I^2*mean(phi)^2 + I*(I-1)*(var(colMeans(phi))+var(rowMeans(phi))) + (sum(diag(phi))-I*mean(phi))^2/(I-1)
## I^2*(2-I)/(I-1)*mean(phi)^2 + I*(sum(colMeans(phi)^2) + sum(rowMeans(phi)^2) )+ (sum(diag(phi))^2 -2*sum(diag(phi))*I*mean(phi) )/(I-1)
## ## 12a-4 generic elt of Q.rho
## p <- sample(1:I,1);r <- sample(1:I,1)
## sum(abs(sapply(1:I,function(p)sapply(1:I,function(r) {
##    idx(Q%*%rho,p,r,I) - (
##     -1/2/I*(phi[p,p]+phi[r,r]-rowMeans(phi)[p]-colMeans(phi)[p]-rowMeans(phi)[r]-colMeans(phi)[r]+4*mean(phi)-2/I*sum(diag(phi)))
##    )
##    }))))



## 12b growth rate of |Q.rho|^2. Seems to be O(1) for unstructured
## phi, at least O(1/sqrt(I)) for phi with increasing rows. But for
## lower right triangular phi (a valid phi ie increasing rows) it is
## O(1).
source('misc.R')
require(parallel)
Is <- 3:30
## parts <- lapply(Is, function(I) {
parts <- mclapply(Is, mc.cores=detectCores()-3,FUN=function(I) {
    print(I)
    Q <- M.obj(I,symm=TRUE)
    replicate(1e1, {
        ## phi <- t(replicate(I,sort(runif(I))))
        ## phi <- phi[,sample(1:I)]
        phi <- matrix(0,I,I)
        phi <- ((row(phi)>=col(phi))+0)[,I:1]
        ## ## phi <- matrix(runif(I^2),I)
        ones <- matrix(1,nrow=I)
        proj <- kronecker(t(ones),rowMeans(phi)) + t(kronecker(t(ones),colMeans(phi))) - kronecker(ones,t(ones))*mean(phi) + (sum(diag(phi))-I*mean(phi))/(I-1)*(diag(I)-1/I)
        rho <- resid <- as.numeric(t(phi-proj))
        phi.vec <- c(t(phi))
        proj.vec <- c(t(proj))
        sum((Q%*%rho)^2)
        c(Q.rho.norm2=var(colMeans(phi)+rowMeans(phi)-diag(phi))*(I-1)/2/I,proj.norm2=sum(proj.vec^2))
        c(Q.rho.norm2=sum((Q%*%rho)^2),Q.proj.norm2=sum((Q%*%proj.vec)^2),proj.norm2=sum(proj.vec^2),cross=(proj.vec%*%Q%*%rho),obj=phi.vec%*%Q%*%phi.vec)
    })
})
parts <- simplify2array(parts)

matplot(t(parts)*Is^(1/2),pch='.',col=rgb(1,0,0,.2),cex=5,x=Is)
plot(Is,apply(parts['Q.rho.norm2',,],2,max))
plot(Is,apply(parts['Q.proj.norm2',,],2,max))
plot(Is,apply(parts['proj.norm2',,],2,max))
plot(Is,apply(parts['cross',,],2,max))
plot(Is,apply(parts['obj',,],2,max))


## 12c. upper right triangular phi also Q-isotropic. 1. Q-orthogonal
## to constant row,col and diag matrices. [update: upper right trianular
## matrix is actually in the nullspace of Q, and so are constant diag matrices]
source('misc.R')
I <- 6
Q <- M.obj(I,sym=TRUE)
T <- with(list(m=matrix(0,I,I)),row(m)<=col(m))+0
T.vec <- c(t(T))
T.vec%*%Q%*%T.vec
ones <- matrix(1,nrow=I)
Id <- diag(I)
b <- kronecker(Id[,1],ones)
t(b)%*%Q%*%T.vec
basis <- cbind(kronecker(ones,Id),kronecker(Id,ones),c(t(diag(I))),T.vec)
sum(abs(t(basis)%*%Q%*%basis)) #1

proj.MD <- function(phi){# project mat onto constant row,col,diag matrices
    I <- nrow(phi)
    ones <- matrix(1,nrow=I)
    kronecker(t(ones),rowMeans(phi)) + t(kronecker(t(ones),colMeans(phi))) - kronecker(ones,t(ones))*mean(phi) + (sum(diag(phi))-I*mean(phi))/(I-1)*(diag(I)-1/I)
}

## checking formula for proj.MD(T)
idx <- function(v,p,r,I)sapply(p,function(p0)sapply(r,function(r0)v[I*(p0-1)+r0]))
proj.MD(T)
j <- sample(1:I,1);i <- sample(1:I,1)
colMeans(T)[j]+rowMeans(T)[i]-mean(T)+(sum(diag(T))-I*mean(T))/(I-1)*((i==j)-1/I)
j/I+(I-i+1)/I-(I+1)/2/I+1/2*((i==j)-1/I)
(j-i)/I+(i==j)/2+1/2
## formula for |T-proj.MD(T)|^2
rho <- resid <- T-proj.MD(T)
sum(rho^2)
js <- 1:(I-1)
sum(-js^3/I^2+2/I*js^2-5/4*js+I/4)
sum((I-js)*(js/I-1/2)^2)
1/12*(I-1)*(I-2)

i <- sample(1:I,1);j <- sample(1:I,1)
sum(abs(sapply(1:I,function(i)sapply(1:I,function(j) {
    ## proj.MD(T)[i,j]-(    (j-i)/I+(i==j)/2+1/2)
    (T-proj.MD(T))[i,j]-    ((i-j)/I-(i==j)/2+(-1)^(j<i)/2 )
}))))




I <- 5
Q <- M.obj(I,sym=TRUE)
T <- with(list(m=matrix(0,I,I)),row(m)<=col(m))+0
T.vec <- c(t(T))
P <- do.call(cbind,unname((sapply(sample(I),function(i)within(list(v=rep(0,I)),v[i] <- 1))))) # IxI perm matrix
Q%*%kronecker(P,P)%*%T.vec


dd



## 13. distributional properties

## 13a. expectation matrix of phi.Q.phi. Seems nonzero, due to
## |{i,j,k,l}|=3 terms. Checked with formulas and simulation below.
require(mvtnorm)
source('misc.R')
I <- 5
Q <- M.obj(I,symm=TRUE)
## phi.vec <- c(t(matrix(runif(1),I,I) + diag(I)*runif(1))) # mu'Qmu term of E
## phi.vec%*%Q%*%phi.vec
ijkl <- expand.grid(rep(list(1:I),4))
ijkl <- ijkl[,4:1]
counts <- apply(ijkl,1,function(r)length(unique(c(r))))
## table(counts)==rev(c(I*(I-1)*(I-2)*(I-3),6*I*(I-1)*(I-2),7*I*(I-1),I))
## u <- unlist(ijkl[493,])
## outer(c('i','j','k','l'),c('i','j','k','l'))
sigs <- apply(ijkl,1,function(r) {
    ## eq.mat <- outer(r,r,`==`)
    ## paste(apply(eq.mat,2,function(c)which.max(c)),collapse='')
    paste(factor(r,labels=c('i','j','k','l')[1:length(unique(r))]),collapse='')
    paste(factor(r,labels=c('i','j','k','l')[order(unique(r))]),collapse='')
})
for(var in c('ijki','ijik','ijkj','ijjj','ijii','ijij','ijji','iiii'))assign(var,runif(1))
rho <- .3
Sigma <- matrix(c(1,rho,rho,1),2)
mu.X <- 0.5; mu.Y <- .7
ijik <- E.10 <- integrate(function(x)(1-pnorm(x,mean=mu.Y))^2*dnorm(x,mu.X),mu.X-7,mu.X+7)$val
ijkj <- E.01 <- integrate(function(x)pnorm(x,mean=mu.X)^2*dnorm(x,mu.Y),mu.Y-7,mu.Y+7)$val
ijki <- E.11 <- integrate(  Vectorize(function(x)  integrate(Vectorize(function(y)pnorm(y,mu.X)*(1-pnorm(x,mu.Y))*dmvnorm(c(x,y),c(mu.X,mu.Y),Sigma)),mu.Y-7,mu.Y+7)$val), mu.X-7,mu.X+7)$val
theta <- integrate(function(y)pnorm(y,mu.X)*dnorm(y,mu.Y),mu.Y-7,mu.Y+7)$val
## use mc for other probabilities
integrate(  Vectorize(function(x)  integrate(Vectorize(function(y)pnorm(y,mu.X)*(1-pnorm(x,mu.Y))*dmvnorm(c(x,y),c(mu.X,mu.Y),Sigma)),mu.Y-7,mu.Y+7)$val), mu.X-7,mu.X+7)$val
## ## ijjj and ijii seem equal, why?
## ijjj <- integrate( Vectorize(function(x)integrate(Vectorize(function(y)pnorm(y,mu.X)*dmvnorm(c(x,y),c(mu.X,mu.Y),Sigma)),x,x+7)$val), mu.X-7,mu.X+7)$val
## ijii <- integrate( Vectorize(function(x)integrate(Vectorize(function(y)(1-pnorm(x,mu.Y))*dmvnorm(c(x,y),c(mu.X,mu.Y),Sigma)),x,x+7)$val), mu.X-7,mu.X+7)$val
probs <- rowMeans(replicate(1e5, {
    xy <- rmvnorm(2,mean=c(mu.X,mu.Y),sigma=Sigma)
    x <- xy[,1]; y <- xy[,2]
    c(ijjj=(x[1]<y[2])*(x[2]<y[2]), ijii=(x[1]<y[2])*(x[1]<y[1]), ijij=x[1]<y[2], ijji=(x[1]<y[2])*(x[2]<y[1]), theta.11=x[1]<y[1])
}))
theta.11 <- probs['theta.11']
covs <- rep(0,I^4)
covs[sigs=='ijki' | sigs=='ijjk'] <- ijki-theta^2
covs[sigs=='ijik' | sigs=='ijkj'] <- ijik-theta^2
## covs[sigs=='ijkj'] <- ijkj
covs[sigs=='ijjj' | sigs=='iiji'] <- probs['ijjj']-theta*theta.11
covs[sigs=='ijii' | sigs=='iiij'] <- probs['ijii']-theta*theta.11
covs[sigs=='ijij'] <- probs['ijij']-theta^2
covs[sigs=='ijji'] <- probs['ijji']-theta^2
covs[sigs=='iiii'] <- theta.11*(1-theta.11)
cov.mat <- matrix(covs,I^2,byrow=TRUE)
sum(diag(Q%*%cov.mat))

quad.forms <- replicate(5e3,{
    xy <- rmvnorm(I,c(mu.X,mu.Y),sigma=Sigma)
    x <- xy[,1]; y <- xy[,2]
    phi <- outer(x,y,'<=')+0
    ## phi <- matrix(0,I,I) #1
    ## phi <- ((row(phi)>=col(phi))+0)[,I:1]
    phi.vec <- c(t(phi))
    t(phi.vec)%*%Q%*%phi.vec
})
mean(quad.forms)

