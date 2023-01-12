## require(pracma)
require(mvtnorm)
require(parallel)
source('misc.R')
I <-3e1
mean.X <- .0; mean.Y <- .5
cor.XX <- .0; cor.YY <- .3; cor.XY <- .2
n.range <- m.range <- 1:10
ests <- replicate(1.5e2, {
    m <- sample(m.range,I,replace=TRUE)
    n <- sample(n.range,I,replace=TRUE)
    ## data <- lapply(1:I, function(i) {
    data <- mclapply(1:I, mc.cores=detectCores()-3,FUN=function(i) {
        sigma <- with(list(m.x=matrix(cor.XX,m[i],m[i]), m.xy=matrix(cor.XY,m[i],n[i]),m.y=matrix(cor.YY,n[i],n[i])), rbind(cbind(m.x,m.xy),cbind(t(m.xy),m.y)))
        diag(sigma) <- 1
        xy <- rmvnorm(1,mean=c(rep(mean.X,m[i]),rep(mean.Y,n[i])),sigma=sigma)
        x <- xy[1:m[i]]; y <- xy[(m[i]+1):(m[i]+n[i])]
        list(x=x,y=y)
    })
    xs <- lapply(data,function(xy)xy$x)
    ys <- lapply(data,function(xy)xy$y)
    ## mapply(auc,xs,ys)
    phi <- outer(xs,ys,Vectorize(auc))
    ## phi.cols <- split(t(phi),1:I)
    ## ## f <- function(x,y)sum(x*y)
    ## m <- outer(phi.cols,phi.cols,Vectorize(function(x,y)sum(x*y)))
    ## diag(m) <- 0
    ## ijik <- sum(m)/I^2/(I-1)
    ijik <- (sum(rowSums(phi)^2)-sum(phi*phi))/I^2/(I-1)    
    ## phi.rows <- split(phi,1:I)
    ## m <- outer(phi.rows,phi.rows,Vectorize(function(x,y)sum(x*y)))
    ## diag(m) <- 0
    ## ijkj <- sum(m)/I^2/(I-1)
    ijkj <- (sum(colSums(phi)^2)-sum(phi*phi))/I^2/(I-1)    
    ## mean(apply(combn(I,2),2,function(idx)phi[idx[1],]*phi[idx[2],]))
    ijki <- (sum(rowSums(phi)*colSums(phi)) - sum(phi*t(phi)))/I^2/(I-1)
    ijij <- (sum(phi^2)-sum(diag(phi)^2))/I/(I-1)
    c(ijik=ijik,ijkj=ijkj,ijki=ijki,ijij=ijij)
})
## integrate(function(x)(1-pnorm(x,mean.Y))^2*dnorm(x,mean.X),mean.X-5,mean.X+5)$val
## x2 <- .3
sigma.XX <- matrix(c(1,cor.XX,cor.XX,1),2)
sigma.YY <- matrix(c(1,cor.YY,cor.YY,1),2)
sigma.XY <- matrix(c(1,cor.XY,cor.XY,1),2)
## integrate(Vectorize(function(x2)integrate(Vectorize(function(x1)(1-pnorm(x1,mean.Y))*(1-pnorm(x2,mean.Y))/(2*pi*det(sigma.XX))*exp(-1/2*(t(c(x1,x2)-rep(mean.X,2))%*%solve(sigma.XX)%*%(c(x1,x2)-rep(mean.X,2))   ))), mean.X-7,mean.X+7)$val),mean.X-5,mean.X+5)$val
## integral2(Vectorize(function(x1,x2)(1-pnorm(x1,mean.Y))*(1-pnorm(x2,mean.Y))*dmvnorm(x1,x2,mean=rep(mean.X,2),sigma=matrix(c(1,corr.XY,corr.XY,1),2))),mean.X-5,mean.X+5,mean.X-5,mean.X+5)
## true <- mean(1/(1:I)) * integrate(function(x)(1-pnorm(x,mean.Y))^2*dnorm(x,mean.X),mean.X-5,mean.X+5)$val  +  mean(((1:I)-1)/(1:I)) * integrate(Vectorize(function(x2)integrate(Vectorize(function(x1)(1-pnorm(x1,mean.Y))*(1-pnorm(x2,mean.Y))/(2*pi*det(sigma.XX))*exp(-1/2*(t(c(x1,x2)-rep(mean.X,2))%*%solve(sigma.XX)%*%(c(x1,x2)-rep(mean.X,2))   ))), mean.X-7,mean.X+7)$val),mean.X-5,mean.X+5)$val
op <- par(mfrow=c(1,4))
true.ijik <-  mean(1/m.range) * integrate(function(x)(1-pnorm(x,mean.Y))^2*dnorm(x,mean.X),mean.X-5,mean.X+5)$val  +  mean((m.range-1)/m.range) * integral2.gaussian(function(x1,x2)(1-pnorm(x1,mean.Y))*(1-pnorm(x2,mean.Y)),rep(mean.X,2),sigma.XX)
hist(ests['ijik',])
abline(v=true.ijik,col=2)
abline(v=mean(ests['ijik',]),col=4)
true.ijkj <-  mean(1/n.range) * integrate(function(y)(pnorm(y,mean.X))^2*dnorm(y,mean.Y),mean.Y-5,mean.Y+5)$val  +  mean((n.range-1)/n.range) * integral2.gaussian(function(y1,y2)pnorm(y1,mean.X)*pnorm(y2,mean.X),mean=rep(mean.Y,2),sigma=sigma.YY)
hist(ests['ijkj',])
abline(v=true.ijkj,col=2)
abline(v=mean(ests['ijkj',]),col=4)
true.ijki <- integral2.gaussian(function(x,y)(1-pnorm(x,mean.Y))*pnorm(y,mean.X),c(mean.X,mean.Y),sigma.XY)
hist(ests['ijki',])
abline(v=true.ijki,col=2)
abline(v=mean(ests['ijki',]),col=4)
f1 <- Vectorize(function(y1,y2)integral2.gaussian(function(x1,x2)1,rep(mean.X,2),sigma.XX,xmax=y1,ymax=y2))
f2 <- Vectorize(function(y)integral2.gaussian(function(x1,x2)1,rep(mean.X,2),sigma.XX,xmax=y1,ymax=y1))
true.ijij <- mean((m.range-1)/m.range)*mean((n.range-1)/n.range)*integral2.gaussian(fun=f1,mean=rep(mean.Y,2),sigma=sigma.YY,vectorized=FALSE) + mean((m.range-1)/m.range)*mean(1/n.range)*integrate(function(y)f2(y)*dnorm(y,mean.Y),mean.Y-5,mean.Y+5)$val + mean((n.range-1)/n.range)*mean(1/m.range^2)*integral2.gaussian(function(y1,y2)pnorm(y1,mean.X)*pnorm(y2,mean.X),rep(mean.Y,2),sigma.YY) + mean(1/m.range)*mean(1/n.range)*integrate(function(y)pnorm(y,mean.X)*dnorm(y,mean.Y),mean.Y-5,mean.Y+5)$val
hist(ests['ijij',])
abline(v=true.ijij,col=2)
abline(v=mean(ests['ijij',]),col=4)
par(op)




## 2. dependence between m,n and cluster auc
source('misc.R')
I <- 3e2
max.cluster.size <- 10
ests <- replicate(1e2, {
    m <- sample(1:max.cluster.size,I,replace=TRUE)
    n <- sample(1:max.cluster.size,I,replace=TRUE)
    ## m <- n <- sample(1:max.cluster.size,I,replace=TRUE) #1
    x <- lapply(m,function(m.i)rnorm(m.i,0))
    y <- lapply(n,function(n.i)rnorm(n.i,n.i/10))
    auc.obu(x,y)
})
hist(ests['var.hat',])
abline(v=var(ests['theta.hat',]),col=2)
abline(v=mean(ests['var.hat',]),col=3)


## 2a. toy example
require(parallel)
source('../auc/misc.R')
I <- 5e2
Is <- round(seq(10,500,len=10))
## by.I <- lapply(Is, function(I) {
by.I <- mclapply(Is, mc.cores=detectCores()-3, FUN=function(I) {
    print(I)
    ests <- replicate(1e2, {
        ## m <- sample(1:2,I,replace=TRUE)
        ## n <- sample(1:2,I,replace=TRUE)
        ## x <- lapply(m,function(m.i)rep(m.i,m.i))
        ## y <- lapply(n,function(n.i)rep(n.i,n.i))
        ## m <- 1+rbinom(I,1,1/2)
        ## n <- 1+rbinom(I,1,1/2)
        ## M <- sum(m); N <- sum(n)
        ## x <- lapply(m,function(m.i) if(m.i==1) 0 else c(0,1))
        ## y <- lapply(n,function(n.i) if(n.i==1) 1 else c(0,1))
    m <- sample(1:(max.cluster.size/2),I,replace=TRUE)
    n <- sample(1:max.cluster.size,I,replace=TRUE)
    ## m <- n <- sample(1:max.cluster.size,I,replace=TRUE) #1
    x <- lapply(m,function(m.i)rnorm(m.i,0))
    y <- lapply(n,function(n.i)rnorm(n.i,n.i/10))
        auc.obu(x,y)
    })
    list(var.mc=I*var(ests['theta.hat',]), var.obu=I*ests['var.hat',])
})
var.obu <- sapply(by.I,function(lst)lst$var.obu)
matplot(Is,t(var.obu),col=1,pch=1)
var.mc <- sapply(by.I,function(lst)lst$var.mc) 
lines(Is,var.mc,col=2)


## 2aa. formula for variance 

reps <- 1e2
Is <- 3:300
by.I <- sapply(Is, function(I) c(mc=var(1/(I+rbinom(reps,I,1/2))),apprx=4/81/I^3))
plot(by.I[1,],by.I[2,]); abline(0,1)


m## 2ab. fla for obu variance
source('misc.R')
I.01 <- I.10 <- I <- 10
m <- 1+rbinom(I,1,1/2)
n <- 1+rbinom(I,1,1/2)
M <- sum(m); N <- sum(n)
x <- lapply(m,function(m.i) if(m.i==1) 0 else c(0,1))
y <- lapply(n,function(n.i) if(n.i==1) 1 else c(0,1))
## ## auc.obu(x,y)
## (theta.hat <- psi(unlist(x),unlist(y)) / (M*N))
## I/2 * (1/M + 1/N)
## all.equal((V.10 <- sapply(x,function(x.i)psi(x.i,unlist(y)))/N), (I*m/N+1)/2)
## all.equal((V.01 <- sapply(y,function(y.i)psi(unlist(x),y.i))/M), (M+I*n)/2/M)
## (S.10 <- sum((V.10 - m*theta.hat)^2) * I.10 / ((I.10-1)*M))
## I/M/(I-1)*sum((1/2-I*m/2/M)^2)
## I^2/M/(I-1)/4*(I*sum(m^2)/M^2-1)
## (S.01 <- sum((V.01 - n*theta.hat)^2) * I.01 / ((I.01-1)*N))
## I^2/4/N/(I-1)*(I*sum(n^2)/N^2-1)
## (S.11 <- sum((V.10 - m*theta.hat)*(V.01 - n*theta.hat)) * I / (I-1))
## I/4/(I-1)*(I*sum(m*n)/M/N-1)
## I^3/4/(I-1)*sum((m/M^2+n/N^2)^2) - I^2/4/(I-1)*(1/M+1/N)^2
## auc.obu(x,y)['var.hat']

source('misc.R')
I.01 <- I.10 <- I <- 10
Is <- seq(5,100,len=5)
by.I <- sapply(Is, function(I) {
    replicate(1e2, {
        m <- 1+rbinom(I,1,1/2)
        n <- 1+rbinom(I,1,1/2)
        M <- sum(m); N <- sum(n)
        x <- lapply(m,function(m.i) if(m.i==1) 0 else c(0,1))
        y <- lapply(n,function(n.i) if(n.i==1) 1 else c(0,1))
        I*auc.obu(x,y)['var.hat']
    })
})
matplot(t(by.I),col=1,pch=1)
abline(h=/81,col=2)



m <- 1+rbinom(4,1,1/2)
n <- 1+rbinom(4,1,1/2)
Is <- round(seq(5,300,len=20))
by.I <- sapply(Is, function(I) {
    pairs <- replicate(1e2, {
        m <- c(m,1+rbinom(I-4,1,1/2))
        n <- c(n,1+rbinom(I-4,1,1/2))
        M <- sum(m); N <- sum(n)
        x <- lapply(m,function(m.i) if(m.i==1) 0 else c(0,1))
        y <- lapply(n,function(n.i) if(n.i==1) 1 else c(0,1))
        psi.12 <- psi(x[[1]],y[[2]])
        psi.34 <- psi(x[[3]],y[[4]])
        c(psi.12,psi.34)*I^2/M/N
    })
    I*cov(t(pairs))[1,2]
})
plot(Is,by.I)

source('misc.R')
max.cluster.size <- 10
Is <- round(seq(5,500,len=50))
by.I <- sapply(Is, function(I) {
    pairs <- replicate(1e4, {
        ## m <- 1+rbinom(I,1,1/2)
        ## n <- 1+rbinom(I,1,1/2)
        ## M <- sum(m); N <- sum(n)
        ## x <- lapply(m[1:4],function(m.i) if(m.i==1) 0 else c(0,1))
        ## y <- lapply(n[1:4],function(n.i) if(n.i==1) 1 else c(0,1))
        m <- sample(1:(max.cluster.size/2),I,replace=TRUE)
        n <- sample(1:max.cluster.size,I,replace=TRUE)
        M <- sum(m); N <- sum(n)
        x <- lapply(m[1:4],function(m.i)rnorm(m.i,0))
        y <- lapply(n[1:4],function(n.i)rnorm(n.i,n.i/10))
        psi.12 <- psi(x[[1]],y[[2]])
        psi.34 <- psi(x[[3]],y[[4]])
        c(psi.12,psi.34)*I^2/M/N
    })
    cov(t(pairs))[1,2]
})
plot(Is,Is*by.I)
plot(Is,by.I)

n <- 50
ns <- seq(10,5e2,len=30)
by.n <- sapply(ns, function(n) {
    pairs <- replicate(5e3, {
        z <- rnorm(n)
        w <- rnorm(n)
        ## c((1+z/sqrt(n))*rnorm(1),(1+z/sqrt(n))*rnorm(1))
        c((1+mean(z)*mean(w))*(rnorm(1)+(30+z[1])),(1+mean(z)*mean(w))*(rnorm(1)+(30+w[1])))
    })
    n*cov(t(pairs))[1,2]
})
plot(ns,by.n)

dd

## m <- 1+rbinom(4,1,1/2)
## n <- 1+rbinom(4,1,1/2)
Is <- round(seq(5,300,len=50))
by.I <- sapply(Is, function(I) {
    pairs <- replicate(1e4, {
        m <- 1+rbinom(I,1,1/2)
        n <- 1+rbinom(I,1,1/2)
        M <- sum(m); N <- sum(n)
        ## x <- lapply(m,function(m.i) if(m.i==1) 0 else c(0,1))
        ## y <- lapply(n,function(n.i) if(n.i==1) 1 else c(0,1))
        psi.12.cond <- (m[1]+n[2])/2
        psi.34.cond <- (m[3]+n[4])/2
        ## c(psi.12*I^2/M/N,psi.34*I^2/M/N)
        c( (M*N/I^2 - (2/3)^2)*psi.12.cond, psi.34.cond)
        ## c( (M*N/I^2 - (2/3)^2), psi.34.cond)
    })
    cov(t(pairs))[1,2]
})
plot(Is,Is*by.I)
abline(h=mean(Is*by.I),col=2)
abline(h=0,col=1)
## matplot(t(by.I),col=1,pch=1)



require(parallel)
source('misc.R')
I <- 5e2
max.cluster.size <- 10
Is <- round(seq(10,500,len=50))
## by.I <- lapply(Is, function(I) {
by.I <- mclapply(Is, mc.cores=detectCores()-3, FUN=function(I) {
    ## print(I)
    ests <- replicate(1e2, {
        ## m <- sample(1:2,I,replace=TRUE)
        ## n <- sample(1:2,I,replace=TRUE)
        ## x <- lapply(m,function(m.i)rep(m.i,m.i))
        ## y <- lapply(n,function(n.i)rep(n.i,n.i))
        ## m <- 1+rbinom(I,1,1/2)
        ## n <- 1+rbinom(I,1,1/2)
        ## M <- sum(m); N <- sum(n)
        ## x <- lapply(m,function(m.i) if(m.i==1) 0 else c(0,1))
        ## y <- lapply(n,function(n.i) if(n.i==1) 1 else c(0,1))
        m <- sample(1:(max.cluster.size/2),I,replace=TRUE)
        n <- sample(1:max.cluster.size,I,replace=TRUE)
        M <- sum(m); N <- sum(n)
        ## m <- n <- sample(1:max.cluster.size,I,replace=TRUE) #1
        x <- lapply(m,function(m.i)rnorm(m.i,0))
        y <- lapply(n,function(n.i)rnorm(n.i,n.i/10))
        V.10 <- sapply(x,function(x.i)psi(x.i,unlist(y)))/N
        V.01 <- sapply(y,function(y.i)psi(unlist(x),y.i))/M
        c(I*(I-1)*(I-2)*(I-3)/M^2/N^2*auc(x[[1]],y[[2]])*auc(x[[3]],y[[4]]),I^2/M/N*auc(x[[1]],y[[2]]))
    })
    mean(ests[1,]) - mean(ests[2,])^2
})
by.I <- simplify2array(by.I)
plot(Is,Is*by.I)






## 3. variance formula for obu estimator


## ## 3b
## I <- 5e1
## rho <- .3
## mu.m <- 3; mu.n <- 5; mu.psi <- 2
## var.mnpsi <- with(list(m=matrix(runif(9),3)),m%*%t(m))
## sd.mnpsi <- with(eigen(var.mnpsi),vectors%*%sqrt(diag(values))%*%t(vectors))
## ## Is <- round(seq(5,500,length.out=40))
## ## by.I <- sapply(Is, function(I) {
## data <- replicate(3, {
##     mnpsi <- with(list(xyz=matrix(rnorm(3*I),nrow=3)),sd.mnpsi%*%xyz)
##     ## m <- mu.m+mnpsi[1,]; n <- mu.n+mnpsi[2,]; psi <- mu.psi+mnpsi[3,]
##     c(m=mu.m+mnpsi[1,], n=mu.n+mnpsi[2,], psi=mu.psi+mnpsi[3,]   )
## },simplify=FALSE)
## names(data) <- c('i','j','k')
## R <- data$i$m*data$j$n*data

## I*var(Ss)/mean(Ss)^2
## ## })
## ## plot(Is,by.I)
## ## abline(h=-4 + 2*(rho+mu.m*mu.n)/(mu.m*mu.n) + (1+mu.m^2)/mu.m^2 + (1+mu.n^2)/mu.n^2,col=2)
## ## abline(h=mean(by.I),col=3)



## 3a second order taylor expansion. error seems about 1/sqrt(I).
I <- 5e1
rho <- .3
mu.m <- 3; mu.n <- 5
var.mn <- matrix(c(1,rho,rho,1),2)
cov.mn <- var.mn[1,2]; var.m <- var.mn[1,1]; var.n <- var.mn[2,2]
sd.mn <- with(eigen(var.mn),vectors%*%sqrt(diag(values))%*%t(vectors))
Is <- round(seq(5,500,length.out=40))
by.I <- sapply(Is, function(I) {
    data <- replicate(1e2, {
        mn <- with(list(xy=matrix(rnorm(2*I),nrow=2)),sd.mn%*%xy)
        m <- mu.m+mn[1,]; n <- mu.n+mn[2,]
        psi <- outer(m,n,'+') + rnorm(I^2)*3
        R <- (sum(psi)-sum(diag(psi)))/I/(I-1)
        S <- mean(m)*mean(n)
        c(R=R,S=S)
    })
    data <- as.data.frame(t(data))
    c(mc=I*var(data$R/data$S),
      taylor=with(data, I*(mean(R)/mean(S))^2 * (var(R)/mean(R)^2 + var(S)/mean(S)^2 - 2*cov(R,S)/mean(R)/mean(S))))
})

plot(Is,by.I['mc',])
points(Is,by.I['taylor',],col=2)
plot(Is,(by.I['mc',]-by.I['taylor',])*sqrt(Is))


## 3b checking formulas and estimators for asy var

## 3ba S (denominator) term
I <- 5e2
rho <- .3
mu.m <- 3; mu.n <- 5
var.mn <- matrix(c(1,rho,rho,1),2)
sd.mn <- with(eigen(var.mn),vectors%*%sqrt(diag(values))%*%t(vectors))
Is <- round(seq(5,500,length.out=40))
by.I <- sapply(Is, function(I) {
    Ss <- replicate(1e3, {
        mn <- with(list(xy=matrix(rnorm(2*I),nrow=2)),sd.mn%*%xy)
        m <- mu.m+mn[1,]; n <- mu.n+mn[2,]
        S <- mean(m)*mean(n)
    })
    I*var(Ss)/mean(Ss)^2
})
plot(Is,by.I)
abline(h=-4 + 2*(rho+mu.m*mu.n)/(mu.m*mu.n) + (1+mu.m^2)/mu.m^2 + (1+mu.n^2)/mu.n^2,col=2)
abline(h=mean(by.I),col=3)


I <- 2e2
rho <- .3
mu.m <- 3; mu.n <- 5
var.mn <- matrix(c(1,rho,rho,1),2)
cov.mn <- var.mn[1,2]; var.m <- var.mn[1,1]; var.n <- var.mn[2,2]
sd.mn <- with(eigen(var.mn),vectors%*%sqrt(diag(values))%*%t(vectors))
Is <- round(seq(1e2,5e2,length.out=40))
by.I <- sapply(Is, function(I) {
    data <- replicate(1e2, {
        mn <- with(list(xy=matrix(rnorm(2*I),nrow=2)),sd.mn%*%xy)
        m <- mu.m+mn[1,]; n <- mu.n+mn[2,]
        ## psi <- outer(m,n,'+') + rnorm(I^2)
        ## R <- with(list(mat=t(psi*m)*n), (sum(mat)-sum(diag(mat)))/I/(I-1))
        S <- mean(m)*mean(n)
       ## c(a= with(list(mat=t(psi*m^2)*n), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  b=with(list(mat=t(psi*m)*n^2), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  c=with(list(mat=t(psi*m*n)*n), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  d=with(list(mat=t(psi*m)*m*n), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  e=with(list(mat=t(psi*m)*n), (sum(mat)-sum(diag(mat)))/I/(I-1)),
        ##  R=R,S=S)
        S.ratio.est <- -4 + 2*mean(m*n)/mean(m)/mean(n) + mean(m^2)/mean(m)^2 + mean(n^2)/mean(n)^2
        c(S=S,S.ratio.est=S.ratio.est)
    })    
    data <- as.data.frame(t(data))
    c(S.ratio.mc=with(data, I*var(S)/mean(S)^2),
      S.ratio.est=mean(data$S.ratio.est))
})
## plot(Is, by.I['cross.mc',]-by.I['cross.formula',])
## plot(Is, by.I['E.RS.mc',]-by.I['E.RS.formula',])
## plot(Is, by.I['cov.RS.mc',]-by.I['cov.RS.formula',])
plot(Is, by.I['S.ratio.mc',],ylim=c(0,.4))
points(Is, by.I['S.ratio.est',],col=2)


## 3bb cov term
I <- 2e2
rho <- .3
mu.m <- 3; mu.n <- 5
var.mn <- matrix(c(1,rho,rho,1),2)
cov.mn <- var.mn[1,2]; var.m <- var.mn[1,1]; var.n <- var.mn[2,2]
sd.mn <- with(eigen(var.mn),vectors%*%sqrt(diag(values))%*%t(vectors))
Is <- round(seq(1e2,5e2,length.out=40))
by.I <- sapply(Is, function(I) {
    data <- replicate(1e2, {
        mn <- with(list(xy=matrix(rnorm(2*I),nrow=2)),sd.mn%*%xy)
        m <- mu.m+mn[1,]; n <- mu.n+mn[2,]
        psi <- outer(m,n,'+') + rnorm(I^2)
        R <- (sum(psi)-sum(diag(psi)))/I/(I-1)#with(list(mat=t(psi*m)*n), (sum(mat)-sum(diag(mat)))/I/(I-1))
        S <- mean(m)*mean(n)
       ## c(a= with(list(mat=t(psi*m^2)*n), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  b=with(list(mat=t(psi*m)*n^2), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  c=with(list(mat=t(psi*m*n)*n), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  d=with(list(mat=t(psi*m)*m*n), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  e=with(list(mat=t(psi*m)*n), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  R=R,S=S)
        ## c(c(a= m[1]^2*n[2],b=m[1]*n[2]^2,c=m[1]*n[1]*n[2],d=m[1]*m[2]*n[2],e=m[1]*n[2])*psi[1,2],
        ## R=R,S=S)
        cross.est <- 1/R/I^2 * sum((colSums(psi)+rowSums(psi))*(m/mean(m)+n/mean(n))) - 4
        c(R=R,S=S,cross.est=cross.est)
    })    
    data <- as.data.frame(t(data))
    c(cross.mc=with(data, I*cov(R,S)/mean(R)/mean(S)),
      cross.est=mean(data$cross.est))
    ## c(cross.mc=with(data, I*cov(R,S)/mean(R)/mean(S)),
    ##   cross.formula=with(as.list(colMeans(data)),-4+(a/mu.m+b/mu.n+c/mu.n+d/mu.m)/e),
    ##   E.RS.mc=with(data, mean(R*S)),
    ##   E.RS.formula=with(as.list(colMeans(data)),e + 1/I*(a*mu.n+b*mu.m+c*mu.m+d*mu.n)),
    ##   cov.RS.mc=with(data, I*cov(R,S)),
    ##   cov.RS.formula=with(as.list(colMeans(data)), -4*e*mu.m*mu.n),
    ##   E.R.E.S.mc=with(data, mean(R)*mean(S)),
    ##   E.R.E.S.formula=mean(data$e*mu.m*mu.n))
})
## plot(Is, by.I['cross.mc',]-by.I['cross.formula',])
## plot(Is, by.I['E.RS.mc',]-by.I['E.RS.formula',])
## plot(Is, by.I['cov.RS.mc',]-by.I['cov.RS.formula',])
plot(Is, by.I['cross.mc',],ylim=c(0,.4))
points(Is, by.I['cross.est',],col=2)



## 3bc R (numerator) term
I <- 2e2
rho <- .3
mu.m <- 3; mu.n <- 5
var.mn <- matrix(c(1,rho,rho,1),2)
cov.mn <- var.mn[1,2]; var.m <- var.mn[1,1]; var.n <- var.mn[2,2]
sd.mn <- with(eigen(var.mn),vectors%*%sqrt(diag(values))%*%t(vectors))
Is <- round(seq(1e2,5e2,length.out=40))
by.I <- sapply(Is, function(I) {
    data <- replicate(1e2, {
        mn <- with(list(xy=matrix(rnorm(2*I),nrow=2)),sd.mn%*%xy)
        m <- mu.m+mn[1,]; n <- mu.n+mn[2,]
        psi <- outer(m,n,'+') + rnorm(I^2)
        ## R <- with(list(mat=t(psi*m)*n), (sum(mat)-sum(diag(mat)))/I/(I-1))
        R <- (sum(psi)-sum(diag(psi)))/I/(I-1)
        S <- mean(m)*mean(n)
       ## c(a= with(list(mat=t(psi*m^2)*n), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  b=with(list(mat=t(psi*m)*n^2), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  c=with(list(mat=t(psi*m*n)*n), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  d=with(list(mat=t(psi*m)*m*n), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  e=with(list(mat=t(psi*m)*n), (sum(mat)-sum(diag(mat)))/I/(I-1)),
       ##  R=R,S=S)
        ## c(a=m[1]^2*n[2]*n[3]*psi[1,2]*psi[1,3],b=m[1]*m[3]*n[2]^2*psi[1,2]*psi[3,2],c=m[1]*m[3]*n[2]*n[1]*psi[1,2]*psi[3,1],d=m[1]*n[2]*psi[1,2],
        ## R=R,S=S)
        c(R=R,S=S,
          R.ratio.est=sum((rowSums(psi)+colSums(psi))^2)/I^3/R^2 - 4
          )
    })    
    data <- as.data.frame(t(data))
    c(R.ratio.mc=with(data,I*var(R)/mean(R)^2),
      R.ratio.est=mean(data$R.ratio.est))
})
## plot(Is, by.I['cross.mc',]-by.I['cross.formula',])
## plot(Is, by.I['E.RS.mc',]-by.I['E.RS.formula',])
## plot(Is, by.I['cov.RS.mc',]-by.I['cov.RS.formula',])
plot(Is, by.I['R.ratio.mc',],ylim=c(0,.4))
points(Is, by.I['R.ratio.est',],col=2)

plot(Is, by.I['R.ratio.mc',]-by.I['R.ratio.est',])
abline(h=mean(by.I['R.ratio.mc',]-by.I['R.ratio.est',]),col=2)
abline(h=0,col=3)



## 3bd obu var estimator in terms of psi
source('misc.R')
I <- 10
max.cluster.size <- 10
m <- sample(max.cluster.size,I,repl=TRUE)
n <- sample(max.cluster.size,I,repl=TRUE)
M <- sum(m); N <- sum(n)
x <- sapply(m, function(m.i)rnorm(m.i))
y <- sapply(n, function(n.i)rnorm(n.i))
psi <- outer(x,y,Vectorize(auc),normalized=FALSE)
theta.hat <- sum(psi)/M/N
ijk <- expand.grid(i=1:I,j=1:I,k=1:I)
sum(apply(ijk, 1, function(r)
    with(as.list(r),
         psi[i,j]*(psi[i,k]/N^2/mean(m)^2+psi[k,j]/M^2/mean(n)^2+2*psi[k,i]/I^2/(mean(m)*mean(n))^2) -2*theta.hat/I*psi[i,j]*(m[i]/N/mean(m)^2 + n[j]/mean(n)^2/M + m[j]/M/mean(m)/mean(n) + n[i]/N/mean(m)/mean(n))+ theta.hat^2/I^2*(m[i]^2/mean(m)^2 + n[i]^2/mean(n)^2 + 2*m[i]*n[i]/mean(m)/mean(n))
         ))) / (I-1)
I*auc.obu(x,y)['var.hat']



## 3be full taylor apprx estimator
source('misc.R')
I <- 2e2
rho <- .3
mu.m <- 3; mu.n <- 5
var.mn <- matrix(c(1,rho,rho,1),2)
cov.mn <- var.mn[1,2]; var.m <- var.mn[1,1]; var.n <- var.mn[2,2]
sd.mn <- with(eigen(var.mn),vectors%*%sqrt(diag(values))%*%t(vectors))
Is <- round(seq(1e1,1e2,length.out=10))
by.I <- sapply(Is, function(I) {
    print(I)
    data <- replicate(1e2, {
        ## mn <- with(list(xy=matrix(rnorm(2*I),nrow=2)),sd.mn%*%xy)
        ## m <- mu.m+mn[1,]; n <- mu.n+mn[2,]
        ## M <- sum(m); N <- sum(n)
        ## psi <- outer(m,n,'+') + rnorm(I^2)
        m <- sample(max.cluster.size,I,repl=TRUE)
        n <- sample(max.cluster.size,I,repl=TRUE)
        M <- sum(m); N <- sum(n)
        x <- sapply(m, function(m.i)rnorm(m.i))
        y <- sapply(n, function(n.i)rnorm(n.i))
        psi <- outer(x,y,Vectorize(auc),normalized=FALSE)
        theta.hat <- sum(psi)/M/N
        R <- (sum(psi)-sum(diag(psi)))/I/(I-1)
        S <- mean(m)*mean(n)
        ## R.ratio.est <- sum((rowSums(psi)+colSums(psi))^2)/I^3/R^2 - 4
        ## ## S.ratio.est <- -4 + 2*mean(m*n)/mean(m)/mean(n) + mean(m^2)/mean(m)^2 + mean(n^2)/mean(n)^2
        ## S.ratio.est <- -4 + I*sum((m/M+n/N)^2)
        ## cross.est <- 1/R/I^2 * sum((colSums(psi)+rowSums(psi))*(m/mean(m)+n/mean(n))) - 4
        ## var.taylor.est <- (R/mean(m)/mean(n))^2*(R.ratio.est+S.ratio.est-2*cross.est)
        var.taylor.est <- (R/mean(m)/mean(n))^2*I*sum((1/R/I^2*(colSums(psi)+rowSums(psi)) - (m/M+n/N))^2)
        c(R=R,S=S,
          var.taylor.est=var.taylor.est,
          var.obu=unname(auc.obu(x,y)['var.hat']))
    })    
    data <- as.data.frame(t(data))
    c(var.taylor.mc=with(data,I*(mean(R)/mean(S))^2 * (var(R)/mean(R)^2 + var(S)/mean(S)^2 - 2*cov(R,S)/mean(R)/mean(S))),
      var.taylor.est=mean(data$var.taylor.est),
      var.obu=I*mean(data$var.obu),
      var.mc=with(data, I*var(R/S)))
})
## plot(Is, by.I['cross.mc',]-by.I['cross.formula',])
## plot(Is, by.I['E.RS.mc',]-by.I['E.RS.formula',])
## plot(Is, by.I['cov.RS.mc',]-by.I['cov.RS.formula',])
plot(Is, by.I['var.mc',])
points(Is, by.I['var.taylor.est',],col=2)
points(Is, by.I['var.obu',],col=3)

plot(Is, by.I['var.taylor.mc',]-by.I['var.taylor.est',])
abline(h=mean(by.I['var.taylor.mc',]-by.I['var.taylor.est',]),col=2)
abline(h=0,col=3)


## rest of #3 below mostly mistakes etc probably not needed

## I*(I/(I-1)/M^2*sum((rowSums(psi)/N-m*theta.hat)^2) + I/(I-1)/N^2*sum((colSums(psi)/M-n*(theta.hat))^2) + 2*I/(I-1)/M/N*sum((rowSums(psi)/N-m*theta.hat)*(colSums(psi)/M-n*theta.hat)))
## sum((rowSums(psi)/N-m*theta.hat)^2)/(I-1)/mean(m)^2 + sum((colSums(psi)/M-n*theta.hat)^2)/(I-1)/mean(n)^2 + sum((rowSums(psi)/N-m*theta.hat)*(colSums(psi)/M-n*theta.hat))*2/(I-1)/mean(m)/mean(n)
## sum(apply(ijk, 1, function(r)
##     with(as.list(r),
##          ## 1/N^2/mean(m)^2*psi[i,j]*psi[i,k] + theta.hat^2/mean(m)^2*m[i]^2/I^2 - 2*theta.hat/N/mean(m)^2/I*m[i]*psi[i,j]
##          ## 1/M^2/mean(n)^2*psi[i,j]*psi[k,j] + theta.hat^2/mean(n)^2*n[i]^2/I^2 - 2*theta.hat/M/mean(n)^2/I*n[j]*psi[i,j]
##          2/mean(m)^2/mean(n)^2/I^2*psi[i,j]*psi[k,i] - 2*theta.hat/I*psi[i,j]*(m[j]/M/mean(m)/mean(n) + n[i]/N/mean(m)/mean(n)) +2*theta.hat^2/mean(m)/mean(n)/I^2* m[i]*n[i]
##     )))/(I-1)
    
## ## 3c estimators of parts of asy variance

## ## 3ca S (denominator) term

## I <- 2e2
## rho <- .3
## mu.m <- 3; mu.n <- 5
## var.mn <- matrix(c(1,rho,rho,1),2)
## cov.mn <- var.mn[1,2]; var.m <- var.mn[1,1]; var.n <- var.mn[2,2]
## sd.mn <- with(eigen(var.mn),vectors%*%sqrt(diag(values))%*%t(vectors))
## Is <- round(seq(1e2,5e2,length.out=40))
## by.I <- sapply(Is, function(I) {
##     data <- replicate(1e2, {
##         mn <- with(list(xy=matrix(rnorm(2*I),nrow=2)),sd.mn%*%xy)
##         m <- mu.m+mn[1,]; n <- mu.n+mn[2,]
##         M <- sum(m); N <- sum(n)
##         ## psi <- outer(m,n,'+') + rnorm(I^2)
##         ## R <- with(list(mat=t(psi*m)*n), (sum(mat)-sum(diag(mat)))/I/(I-1))
##         S <- mean(m)*mean(n)
##         ## S.ratio.est <- -4 + 2*mean(m*n)/mean(m)/mean(n) + mean(m^2)/mean(m)^2 + mean(n^2)/mean(n)^2
##         S.ratio.est <- -4 + I*sum((m/M + n/N)^2)
##         c(S=S,S.ratio.est=S.ratio.est)
##     })    
##     data <- as.data.frame(t(data))
##     c(S.ratio.mc=with(data, I*var(S)/mean(S)^2),
##       S.ratio.est=mean(data$S.ratio.est))
## })
## plot(Is, by.I['S.ratio.mc',],ylim=c(0,.4))
## points(Is, by.I['S.ratio.est',],col=2)


## ## 3cb cov term
## I <- 2e2
## rho <- .3
## mu.m <- 3; mu.n <- 5
## var.mn <- matrix(c(1,rho,rho,1),2)
## cov.mn <- var.mn[1,2]; var.m <- var.mn[1,1]; var.n <- var.mn[2,2]
## sd.mn <- with(eigen(var.mn),vectors%*%sqrt(diag(values))%*%t(vectors))
## Is <- round(seq(1e1,1e2,length.out=20))
## by.I <- sapply(Is, function(I) {
##     print(I)
##     data <- replicate(1e2, {
##         mn <- with(list(xy=matrix(rnorm(2*I),nrow=2)),sd.mn%*%xy)
##         m <- mu.m+mn[1,]; n <- mu.n+mn[2,]
##         M <- sum(m); N <- sum(n)
##         psi <- outer(m,n,'+') + rnorm(I^2)
##         theta.hat <- sum(psi)/M/N
##         R <- with(list(mat=t(psi*m)*n), (sum(mat)-sum(diag(mat)))/I/(I-1))
##         S <- mean(m)*mean(n)
##         mat.small <- psi * outer(m,n)
##         mat.big <- mat.small * (outer(m,m,'+')/mean(m) + outer(n,n,'+')/mean(n))
##         c(R=R,S=S,cross.est= (sum(mat.big)-sum(diag(mat.big)))/(sum(mat.small)-sum(diag(mat.small))) -4
##           )
##     })    
##     data <- as.data.frame(t(data))
##     c(cross.mc=with(data, I*cov(R,S)/mean(R)/mean(S)),
##       cross.est=mean(data$cross.est))
## })
## plot(Is, by.I['cross.mc',],ylim=c(0,1))
## points(Is, by.I['cross.est',],col=2)


## ## 3cc R (numerator) term
## I <- 5
## ijk <- expand.grid(i=1:I,j=1:I,k=1:I)
## ijk.unique <- ijk[apply(ijk,1,function(r)sum(duplicated(r))==0),]
## g <- function(x,y)x*y
## h <- function(x,y)x^2/y
## sum(apply(ijk.unique,1,function(row)with(as.list(row),g(i,j)*h(i,k))))
## sum(apply(ijk,1,function(row)with(as.list(row),g(i,j)*h(i,k)-1/I*(g(i,i)*h(i,j)+h(i,i)*g(i,j)+g(i,j)*h(i,j)) + 2/I^2*g(i,i)*h(i,i))))



## I <- 2e2
## rho <- .3
## mu.m <- 3; mu.n <- 5
## var.mn <- matrix(c(1,rho,rho,1),2)
## cov.mn <- var.mn[1,2]; var.m <- var.mn[1,1]; var.n <- var.mn[2,2]
## sd.mn <- with(eigen(var.mn),vectors%*%sqrt(diag(values))%*%t(vectors))
## Is <- round(seq(1e2,5e2,length.out=40))
## by.I <- sapply(Is, function(I) {
##     print(I)
##     data <- replicate(1e2, {
##         mn <- with(list(xy=matrix(rnorm(2*I),nrow=2)),sd.mn%*%xy)
##         m <- mu.m+mn[1,]; n <- mu.n+mn[2,]
##         M <- sum(m); N <- sum(n)
##         psi <- outer(m,n,'+') + rnorm(I^2)
##         R <- (sum(psi)-sum(diag(psi)))/I/(I-1)
##         ## R <- with(list(mat=t(psi*m)*n), (sum(mat)-sum(diag(mat)))/I/(I-1))
##         S <- mean(m)*mean(n)
##         mat.small <- psi * outer(m,n)
##         R.ratio.est <- sum((rowSums(mat.small)+colSums(mat.small))^2) * (I-1) / sum(mat.small)^2 - 4
##         c(R=R,R.ratio.est=R.ratio.est)
##     })    
##     data <- as.data.frame(t(data))
##     c(R.ratio.mc=with(data,I*var(R)/mean(R)^2),
##       R.ratio.est=mean(data$R.ratio.est))
## })
## ## plot(Is, by.I['cross.mc',]-by.I['cross.formula',])
## ## plot(Is, by.I['E.RS.mc',]-by.I['E.RS.formula',])
## ## plot(Is, by.I['cov.RS.mc',]-by.I['cov.RS.formula',])
## plot(Is, by.I['R.ratio.mc',],ylim=c(0,.4))
## points(Is, by.I['R.ratio.est',],col=2)




## ## 3cd full  estimator of first order variance approx
## source('misc.R')
## I <- 2e2
## rho <- .3
## mu.m <- 3; mu.n <- 5
## var.mn <- matrix(c(1,rho,rho,1),2)
## cov.mn <- var.mn[1,2]; var.m <- var.mn[1,1]; var.n <- var.mn[2,2]
## sd.mn <- with(eigen(var.mn),vectors%*%sqrt(diag(values))%*%t(vectors))
## Is <- round(seq(1e2,2e2,length.out=15))
## by.I <- sapply(Is, function(I) {
##     print(I)
##     data <- replicate(1e1, {
##         ## mn <- with(list(xy=matrix(rnorm(2*I),nrow=2)),sd.mn%*%xy)
##         ## m <- mu.m+mn[1,]; n <- mu.n+mn[2,]
##         ## M <- sum(m); N <- sum(n)
##         ## psi <- outer(m,n,'+') + rnorm(I^2)
##         max.cluster.size <- 10
## m <- sample(max.cluster.size,I,repl=TRUE)
## n <- sample(max.cluster.size,I,repl=TRUE)
## M <- sum(m); N <- sum(n)
## x <- sapply(m, function(m.i)rnorm(m.i))
## y <- sapply(n, function(n.i)rnorm(n.i))
## psi <- outer(x,y,Vectorize(auc),normalized=FALSE)
## theta.hat <- sum(psi)/M/N
##         R <- with(list(mat=t(psi*m)*n), (sum(mat)-sum(diag(mat)))/I/(I-1))
##         S <- mean(m)*mean(n)
##         mat.small <- psi * outer(m,n)
##         mat.big <- mat.small * (outer(m,m,'+')/mean(m) + outer(n,n,'+')/mean(n))
##         ## cross.est <- (sum(mat.big)-sum(diag(mat.big)))/(sum(mat.small)-sum(diag(mat.small))) -4
##         cross.est <- sum(mat.big)/sum(mat.small) #- 4      
##         R.ratio.est <- sum((rowSums(mat.small)+colSums(mat.small))^2) * (I-1) / sum(mat.small)^2 #- 4
##         S.ratio.est <- I*sum((m/M + n/N)^2) #- 4
##         c(R=R,S=S,
##           ## R.ratio.est=R.ratio.est,
##           ## cross.est= cross.est,
##           ## cross.est= sum(mat.big)/sum(mat.small) -4,
##           ## S.ratio.est=S.ratio.est,
##           var.est=(R/S)^2*(R.ratio.est + S.ratio.est - 2*cross.est),
##           var.obu=I*auc.obu(x,y)['var.hat']
##           )
##     })    
##     data <- as.data.frame(t(data))
##     ## with(data,c(I*var(R)/mean(R)^2, I*var(S)/mean(S)^2, - 2*I*cov(R,S)/mean(R)/mean(S)))
##     ## with(data, c(mean(R.ratio.est), + mean(S.ratio.est), - 2*mean(cross.est)))
##        c(var.mc=with(data, I*(mean(R)/mean(S))^2 * (var(R)/mean(R)^2 + var(S)/mean(S)^2 - 2*cov(R,S)/mean(R)/mean(S))),
##          var.est=mean(data$var.est),
##          var.obu=mean(data$var.obu)
##                   )
## })
## ## plot(Is, by.I['cross.mc',]-by.I['cross.formula',])
## ## plot(Is, by.I['E.RS.mc',]-by.I['E.RS.formula',])
## ## plot(Is, by.I['cov.RS.mc',]-by.I['cov.RS.formula',])
## plot(Is, by.I['var.mc',],ylim=c(0,1))
## points(Is, by.I['var.est',],col=2)
## points(Is, by.I['var.obu',],col=3)



## 4. formula for covariance. Note R,S have different defns from #3.

## 4a var(S)
I <- 100
Is <- round(seq(10,400,len=30))
max.cluster.size <- 10
by.I <- sapply(Is, function(I) {
    sim <- replicate(1e2, {
        m <- sample(max.cluster.size, I, repl=TRUE)
        n <- sample(max.cluster.size, I, repl=TRUE)
        S <- mean(m)*mean(n)
        var.S.hat <- (I-1)*(I-1)/I^3*(mean(m^2)*mean(n)^2+mean(n^2)*mean(m)^2+2*mean(m*n)*mean(m)*mean(n)-4*(mean(m)*mean(n))^2)
        c(S=S,var.S.hat=var.S.hat)
    })
    sim <- as.data.frame(t(sim))
    I*c(var.S.mc=var(sim$S),var.S.hat=mean(sim$var.S.hat))
})
plot(Is,by.I['var.S.mc',])
points(Is,by.I['var.S.hat',],col=2)


## 4b
require(parallel)
I <- 2e2
rho <- .3
mu.m <- 3; mu.n <- 5
psi.12 <- mu.m+mu.n
phi.11 <- log((mu.m+.5)/(mu.m-.5)) + log((mu.n+.5)/(mu.n-.5))
var.mn <- matrix(c(1,rho,rho,1),2)
cov.mn <- var.mn[1,2]; var.m <- var.mn[1,1]; var.n <- var.mn[2,2]
sd.mn <- with(eigen(var.mn),vectors%*%sqrt(diag(values))%*%t(vectors))
Is <- round(seq(1e1,7e2,length.out=5e1))
## by.I <- sapply(Is, function(I) {
by.I <- mclapply(Is, mc.cores=detectCores()-3, FUN=function(I) {
    data <- replicate(1e1, {
        mn <- with(list(xy=matrix(runif(2*I,-.5,.5),nrow=2)),sd.mn%*%xy)
        m <- mu.m+mn[1,]; n <- mu.n+mn[2,]
        M <- sum(m); N <- sum(n)
        psi <- outer(m,n,'+') + rnorm(I^2)
        phi <- diag(psi)/m/n
        theta.hat.11 <- mean(diag(psi)/m/n)
        theta.hat.12 <- (sum(psi)-sum(diag(psi)))/M/N * I/(I-1)
        S <- mean(m)*mean(n)
        S.hat <- cov(m,n)/I+mean(m)*mean(n)
        R <- theta.hat.11*theta.hat.12*S
        ## R <- (sum(psi)-sum(diag(psi)))*sum(diag(psi)/m/n)/I^2/(I-1)
        psi.12.hat <- (sum(psi)-sum(diag(psi)))/I/(I-1)
        ## hat <- 1/I/mean(m)/mean(n) * sum( diag(psi)/m/n*(colMeans(psi)+rowMeans(psi)) - (m/mean(m)+n/mean(n))^2 + theta.hat.11*(colMeans(psi)+rowMeans(psi))*(n/mean(n)+m/mean(m)))  - 2*theta.hat.11*theta.hat.12 + 4/mean(m)/mean(n)
        hat.a.1 <- 1/S*mean((colMeans(psi)+rowMeans(psi))*phi)  - 2*theta.hat.12*theta.hat.11
        hat.a.2 <- ( mean((colMeans(psi)+rowMeans(psi))*(n*mean(m)+m*mean(n))) - 4*mean(m)*mean(n)*psi.12.hat) * theta.hat.11 / S^2
        hat.a.3 <- (mean((m/mean(m)+n/mean(n))^2)-4) * theta.hat.11 / S
        hat.a <- 1/S * mean( phi*(colMeans(psi)+rowMeans(psi)) - psi.12.hat*theta.hat.11*(m/mean(m)+n/mean(n))^2 + theta.hat.11*(colMeans(psi)+rowMeans(psi))*(m/mean(m)+n/mean(n)) ) - 6*theta.hat.11*theta.hat.12 + 4*theta.hat.11/S*psi.12.hat
        hat.b <- 1/S*theta.hat.11*psi.12.hat*mean( (m/mean(m)+n/mean(n))*(rowMeans(psi)/psi.12.hat+phi/mean(phi)) + colMeans(psi)/psi.12.hat*(m/mean(m)+n/mean(n)) - 6)
        hat.c <- theta.hat.11*theta.hat.12*(mean((m/mean(m)+n/mean(n))^2)-4)
        hat <- hat.a-hat.b+hat.c
        hat <- mean(phi/mean(phi)*((colMeans(psi)+rowMeans(psi))/mean(psi) - m/mean(m)-n/mean(n))) * theta.hat.11*theta.hat.12
        c(R=R,S=S,theta.hat.11=theta.hat.11,theta.hat.12=theta.hat.12,psi.12.hat=psi.12.hat,hat.a.1=hat.a.1,hat.a.2=hat.a.2,hat.a.3=hat.a.3,hat.a=hat.a,hat.b=hat.b,hat.c=hat.c,hat=hat)
    })
    data <- as.data.frame(t(data))
    mc.a.1 <- with(data, I/mean(S)*(mean(R) - psi.12*phi.11))
    mc.a.2 <- with(data, I*phi.11/mean(S)^2*cov(psi.12.hat,S))
    mc.a.3 <- with(data, I*phi.11/mean(S)^3*var(S))
    mc.a <- mc.a.1-mc.a.2+mc.a.3
    mc.a <- with(data, I*(mean(R)/mean(S) - phi.11*mean(theta.hat.12)))
    mc.b <- with(data, I*cov(R,S)/mean(S)^2)
    mc.c <- with(data, I*var(S)*mean(R)/mean(S)^3)
    mc <- with(data, I*cov(theta.hat.11,theta.hat.12))
    ## mc <- with(data, mean(R/S) - mean(theta.hat.11)*mean(theta.hat.12))
    taylor <- with(data, I*(mean(R)/mean(S) -mean(theta.hat.11)*mean(theta.hat.12)- cov(R,S)/mean(S)^2 + var(S)*mean(R)/mean(S)^3))
    taylor <- with(data, I*(mean(R)/mean(S) - phi.11*mean(theta.hat.12)) - I*cov(R,S)/mean(S)^2+I*var(S)*mean(R)/mean(S)^3)
    ## mc <- with(data, I*(mean(R)/mean(S)-phi.11*mean(theta.hat.12)-cov(R,S)/mean(S)^2+var(S)*mean(R)/mean(S)^3))
    ## mc <- with(data, I*(mean(R/S) - mean(theta.hat.11)*mean(theta.hat.12)))
    c(mc.a.1=mc.a.1,mc.a.2=mc.a.2,mc.a.3=mc.a.3,
      mc.a=mc.a,mc.b=mc.b,mc.c=mc.c,mc=mc,taylor=taylor,
      hat.a.1=mean(data$hat.a.1),hat.a.2=mean(data$hat.a.2),hat.a.3=mean(data$hat.a.3),
      hat.a=mean(data$hat.a),
      hat.b=mean(data$hat.b),
      hat.c=mean(data$hat.c),
      hat=mean(data$hat)
      )
})
by.I <- simplify2array(by.I)
op <- par(mfrow=c(1,3))
plot(Is,by.I['mc.a',]-by.I['hat.a',]); abline(h=0)
plot(Is,by.I['mc.b',]-by.I['hat.b',]); abline(h=0)
plot(Is,by.I['mc.c',]-by.I['hat.c',]); abline(h=0)
par(op)

op <- par(mfrow=c(1,3))
plot(Is,(by.I['mc',]-by.I['taylor',])); abline(h=0); abline(h=mean(by.I['mc',]-by.I['taylor',]),col=2)
plot(Is,by.I['taylor',]-by.I['hat',]); abline(h=0); abline(h=mean(by.I['taylor',]-by.I['hat',]),col=2)
plot(Is,by.I['mc',]-by.I['hat',]); abline(h=0);  abline(h=mean(by.I['mc',]-by.I['hat',]),col=2)
par(op)

op <- par(mfrow=c(1,3))
plot(Is,by.I['mc.a.1',]-by.I['hat.a.1',]); abline(h=0); abline(h=mean(by.I['mc.a.1',]-by.I['hat.a.1',]),col=2)
plot(Is,by.I['mc.a.2',]-by.I['hat.a.2',]); abline(h=0); abline(h=mean(by.I['mc.a.2',]-by.I['hat.a.2',]),col=2)
plot(Is,by.I['mc.a.3',]-by.I['hat.a.3',]); abline(h=0); abline(h=mean(by.I['mc.a.3',]-by.I['hat.a.3',]),col=2)
par(op)


op <- par(mfrow=c(1,3))
plot(Is,by.I['mc.a.1.prime',]-by.I['hat.a.1',]);abline(h=0)
plot(Is,by.I['taylor2',]-by.I['hat2',]);abline(h=0)
plot(Is,by.I['taylor3',]-by.I['hat3',]);abline(h=0)
par(op)

rowMeans(by.I)
plot(Is, by.I['mc.cov',])
points(Is, by.I['hat.cov',],col=2)
plot(Is, by.I['mc.cov',]-by.I['hat.cov',]);abline(h=0)
summary(by.I['mc.cov',]-by.I['hat.cov',])

plot(Is,by.I['mc',])
points(Is,by.I['taylor',],col=2)
points(Is,by.I['hat',],col=3)


plot(Is,by.I['mc',]-by.I['hat',])
points(Is,by.I['taylor',]-by.I['hat',],col=2)
abline(h=0)
abline(h=mean(by.I['mc',]-by.I['hat',]),col=2)

## save.image('4b.RData')
## 297F.02
## s <- 'abcdefghijklmn'
## ## first split on 'f'
## first.split <- strsplit(s, 'f')[[1]]
## first.split

## ## then split on 'k' (two ways given below: first with sapply, then with a loop)
## second.split <- sapply(first.split, function(smaller.string) strsplit(smaller.string, 'k'))
## second.split
## for (j in 1:length(first.split)) print(strsplit(first.split[j], 'k'))



## 4.5 check hajek projection formula using #1 hajek projection itself, #2 defective kernel, #3 M=N=1 case.
source('misc.R')
## I <- 10
Delta <- .5
E.M <- E.N <- 5
Is <- round(seq(5,1e2,len=30))
by.I <- sapply(Is, function(I) {
    residues <- replicate(1e1, {
        M <- rpois(I,E.M-1)+1
        N <- rpois(I,E.N-1)+1
        ## M <- N <- rep(E.M,I)
        x <- lapply(M,function(m)rnorm(m))
        y <- lapply(N,function(n)rnorm(n)+Delta)
        ## U <- mean(outer(x,y,Vectorize(auc),norm=FALSE))
        ## U.hat <- E.N*(E.M - sum(pnorm(unlist(x)-Delta))/I) + E.M*sum(pnorm(unlist(y)))/I - 2*E.M*E.N*pnorm(Delta/sqrt(2))
        ## E.U <- E.M*E.N*pnorm(Delta/sqrt(2))
        ## I*(U-E.U-U.hat)
        table <- outer(x,y,Vectorize(auc),norm=FALSE) ## #1 hajek projection
        diag(table) <- 0
        W <- sum(table)/I/(I-1)
        U.hat <- E.N*mean(M) - E.N*sum(pnorm(unlist(x)-Delta))/I + E.M*sum(pnorm(unlist(y)))/I - 2*E.M*E.N*pnorm(Delta/sqrt(2))
        E.U <- E.M*E.N*pnorm(Delta/sqrt(2))
        W-E.U-U.hat
        ## table <- outer(x,y,Vectorize(auc),norm=FALSE) ## #2 using defective kernel
        ## diag(table) <- 0
        ## E.U <- E.M*E.N*pnorm(Delta/sqrt(2))
        ## given.x <- E.N*(M-sapply(split(pnorm(unlist(x)-Delta),rep(1:I,M)),sum))
        ## given.y <- E.M*sapply(split(pnorm(unlist(y)),rep(1:I,N)),sum)
        ## psi.bar.table <- table - outer(given.x,given.y,`+`) + E.U
        ## diag(psi.bar.table) <- 0
        ## sum(psi.bar.table)/I/(I-1)
        ## x <- rnorm(I);y <- rnorm(I) ## #3 M=N=1 case
        ## table <- outer(unlist(x),unlist(y),'<')
        ## diag(table ) <- 0
        ## U <- sum(table)/I/(I-1)
        ## E.U <- pnorm(Delta/sqrt(2))
        ## U.hat <- mean(1-pnorm(unlist(x)-Delta) - E.U) + mean(pnorm(unlist(y))-E.U)
 ## bb <- U-E.U-U.hat
        ## stopifnot(abs(aa-bb)<1e5)
        ## bb
    })
    mean(residues^2)
})
plot(Is,Is^2*by.I)

source('misc.R')
## I <- 10
E.M <- E.N <- 5
Is <- round(seq(5,1e2,len=30))
by.I <- sapply(Is, function(I) {
    residues <- replicate(1e1, {
        M <- rpois(I,E.M-1)+1
        N <- rpois(I,E.N-1)+1
        W <- sum(outer(M,N))/I/(I-1) - E.M*E.N
        U.hat <- mean(E.N*M + N*E.M) - 2*E.M*E.N
        W-U.hat
        ## ## M <- N <- rep(E.M,I)
        ## x <- lapply(M,function(m)rnorm(m))
        ## y <- lapply(N,function(n)rnorm(n)+Delta)
        ## ## U <- mean(outer(x,y,Vectorize(auc),norm=FALSE))
        ## ## U.hat <- E.N*(E.M - sum(pnorm(unlist(x)-Delta))/I) + E.M*sum(pnorm(unlist(y)))/I - 2*E.M*E.N*pnorm(Delta/sqrt(2))
        ## ## E.U <- E.M*E.N*pnorm(Delta/sqrt(2))
        ## ## I*(U-E.U-U.hat)
        ## table <- outer(x,y,Vectorize(auc),norm=FALSE) ## #1 hajek projection
        ## diag(table) <- 0
        ## W <- sum(table)/I/(I-1)
        ## U.hat <- E.N*mean(M) - E.N*sum(pnorm(unlist(x)-Delta))/I + E.M*sum(pnorm(unlist(y)))/I - 2*E.M*E.N*pnorm(Delta/sqrt(2))
        ## E.U <- E.M*E.N*pnorm(Delta/sqrt(2))
        ## W-E.U-U.hat
    })
    mean(residues^2)
    ## var(residues)
    ## abs(mean(residues))
})
plot(Is,Is^2*by.I)






aa <- integrate(function(x)(1-pnorm(x-Delta))^2*dnorm(x),-5,5)$val + integrate(function(x)pnorm(x+Delta)^2*dnorm(x),-5,5)$val + 2*integrate(function(x)(1-pnorm(x-Delta))*dnorm(x),-5,5)$val*integrate(function(x)pnorm(x+Delta)*dnorm(x),-5,5)$val
theta.12 <- pnorm(Delta/sqrt(2))
var.formula <- aa - 4*theta.12^2

I <- 10
Delta <- .3
Is <- round(seq(10,1e3,len=50))
by.I <- sapply(Is, function(I) {
    theta.12 <- pnorm(Delta/sqrt(2))
    theta.12.hats <- replicate(1e2, {
        x <- rnorm(I)
        y <- rnorm(I)+Delta
        theta.12.hat <- mean(outer(x,y,'<'))
    })
    var.obs <- var(sqrt(I)*(theta.12.hats - theta.12))
})
plot(Is,by.I)
abline(h=var.formula)

source('misc.R')
I <- 1e2
Delta <- .2
## Is <- round(seq(10,1e3,len=30))
## by.I <- sapply(Is, function(I) {
theta.12 <- pnorm(Delta/sqrt(2))
hats <- replicate(3e2, {
    ## cat('.')
    x <- rnorm(I)
    y <- rnorm(I)+Delta
    psi <- outer(x,y,'<')
    theta.12.hat <- mean(psi)
    var.hat <- mean(((rowSums(psi)+colSums(psi))/I-theta.12.hat*2)^2  )
    c(theta.12.hat=sqrt(I)*(theta.12.hat-theta.12),var.hat=var.hat)
    ## with(auc.cluster(x,y),c(var.hat=I*vcov.hat[2,2],theta.12.hat=theta.12.hat))
})
hist(hats['var.hat',])
abline(v=var(hats['theta.12.hat',]),col=2)
abline(v=mean(hats['var.hat',]),col=3)







n <- 1e4
var1 <- mu1 <- 5
var2 <- mu2 <- 7
samples <- replicate(1e3, {
x <- rpois(n,mu1)
y <- rpois(n,mu2)
## sqrt(n)*(mean(x)-mu1)/sqrt(var1)
var.ratio <- 1/mu2^2*(var1+(mu1/mu2)^2*var2)
sqrt(n)*(mean(x)/mean(y) - mu1/mu2)/sqrt(var.ratio)
})
qqnorm(samples); abline(0,1)

dd


## 5 relationship between 2 auc's

I <- 3
reps <- 1e2
gs <- replicate(reps, {
    g <- t(replicate(I,sort(runif(I))))
    g <- g - apply(g,1,min)
    g <- g / apply(g,1,max)
})
fs <-  replicate(reps, rep(1/I,I))#with(list(u=runif(I)),u/sum(u)))
diffs <- sapply(1:reps, function(j) {
    g <- gs[,,j]; f <- fs[,j]
    sum(f*diag(g)) - f%*%g%*%f
})
summary(diffs)
orders <- apply(apply(gs[,2,],2,order),2,paste,collapse='')
order.diffs <- data.frame(orders,diffs)[order(orders),]


I <- 2
    theta.hat.11 <- mean(x<y)
    phi <- outer(x,y,'<')
    theta.hat.12 <- (sum(phi)-sum(diag(phi)))/I/(I-1)
    diff(c(theta.hat.11,theta.hat.12))
})
summary(diffs)
hist(diffs)


I <- 7
diffs <- replicate(1e3, {
g <- replicate(I, ecdf(runif(I)))
y.atoms <- runif(I)
theta.11 <- mean(sapply(1:I,function(i)g[[i]](y.atoms[i])))
theta.12 <- mean(sapply(g, function(cdf)mean(cdf(y.atoms))))
theta.11-theta.12
})
summary(diffs)



I <- 3
diffs <- replicate(1e2,{
    x <- sort(runif(I)); y <- sort(runif(I))
    f.xy <- with(list(m=matrix(runif(I^2),I)), m/sum(m))
    phi <- outer(x,y,'<')
    theta.11 <- sum(f.xy*phi)
    f.x <- rowSums(f.xy); f.y <- colSums(f.xy)
    f.xy.ind <- outer(f.x,f.y,'*')
    theta.12 <- sum(f.xy.ind*phi)
    abs(theta.11-1/2)-abs(theta.12-1/2)
})



I <- 200
pairs <- replicate(5e3, {
    ## phi <- matrix(sample(0:1,4,replace=TRUE),2)
    phi <- t(sapply(sample(0:I,I,repl=TRUE), function(a)c(rep(0,I-a),rep(1,a))))
    p <- with(list(u=runif(I)),u/sum(u))
    theta.11 <- diag(phi)%*%p
    theta.12 <- p%*%phi%*%p
    ## c(theta.11=(theta.11-1/2),theta.12=(theta.12-1/2)) # without abs values, symmetric about origin.
    c(theta.11=abs(theta.11-1/2),theta.12=abs(theta.12-1/2))
})
## plot(pairs[1,],pairs[2,],asp=1);abline(0,1)
mean(pairs[1,]>pairs[2,])
mean(pairs[1,]==pairs[2,])





## 6. simulations

## 6a my own version of mixtools::ellipse
require(mvtnorm)
require(mixtools)
draw.ellipse <- function(mean.xy=c(0,0),cov.xy=diag(2),alpha=.05,resolution=1e2,...) {
    q <- qchisq(1-alpha,df=2)
    cov.xy.sqrt <- with(eigen(cov.xy), vectors%*%diag(sqrt(values))%*%t(vectors))
    theta <- seq(0,2*pi,len=resolution)
    circle <- rbind(cos(theta),sin(theta))
    ellipse <- mean.xy + cov.xy.sqrt%*%(sqrt(q)*circle)
    ## browser()
    lines(t(ellipse),...)
}


resolution <- 1e2
cov.xy <- matrix(runif(4),nrow=2); cov.xy <- cov.xy%*%t(cov.xy)
mean.xy <- runif(2)
alpha <- .05
xy <- mvtnorm::rmvnorm(1e3,mean=mean.xy,sigma=cov.xy)
plot(xy)
draw.ellipse(mean.xy,cov.xy,lty=1)
mixtools::ellipse(mu=mean.xy,sigma=cov.xy,col=2)


## 6b simulation to test conf ellipse
source('misc.R')
I <- 2e2
rho <- .7
mu.m <- 3; mu.n <- 5
psi.12 <- mu.m+mu.n
phi.11 <- log((mu.m+.5)/(mu.m-.5)) + log((mu.n+.5)/(mu.n-.5))
var.mn <- matrix(c(1,rho,rho,1),2)
cov.mn <- var.mn[1,2]; var.m <- var.mn[1,1]; var.n <- var.mn[2,2]
sd.mn <- with(eigen(var.mn),vectors%*%sqrt(diag(values))%*%t(vectors))
Is <- round(seq(1e1,7e2,length.out=5e1))
## by.I <- sapply(Is, function(I) {
## by.I <- mclapply(Is, mc.cores=detectCores()-3, FUN=function(I) {
data <- replicate(1e3, {
    mn <- with(list(xy=matrix(runif(2*I,-.5,.5),nrow=2)),sd.mn%*%xy)
    m <- mu.m+mn[1,]; n <- mu.n+mn[2,]
    M <- sum(m); N <- sum(n)
    psi <- outer(m,n,'+') + rnorm(I^2)
    phi <- diag(psi)/m/n
    theta.11.hat <- mean(diag(psi)/m/n)
    theta.12.hat <- (sum(psi)-sum(diag(psi)))/M/N * I/(I-1)
    cov.hat <- mean(phi/mean(phi)*((colMeans(psi)+rowMeans(psi))/mean(psi) - m/mean(m)-n/mean(n))) * theta.11.hat*theta.12.hat
    var.11.hat <- var(phi)
    ## R <- (sum(psi)-sum(diag(psi)))/I/(I-1)
    var.12.hat <- theta.12.hat^2*I*sum((1/theta.12.hat/M/N*(colSums(psi)+rowSums(psi)) - (m/M+n/N))^2)
    c(theta.11.hat=theta.11.hat,theta.12.hat=theta.12.hat,cov.hat=cov.hat,var.11.hat=var.11.hat,var.12.hat=var.12.hat)
})
data <- data.frame(t(data))

sigma.mc <- with(data, matrix(c(var(theta.11.hat),rep(cov(theta.11.hat,theta.12.hat),2),var(theta.12.hat)),2))
mean.hat <- c(mean(data$theta.11.ha),mean(data$theta.12.hat))
plot(theta.12.hat ~ theta.11.hat, data=data)
draw.ellipse(mean.hat,sigma.mc,col=2)
## with(data,c(cov(theta.11.hat,theta.12.hat),mean(cov.hat)/I))
## with(data,c(var(theta.11.hat),mean(var.11.hat)/I))
## with(data,c(var(theta.12.hat),mean(var.12.hat)/I))
sigma.hat.mc <- with(as.list(colMeans(data)),matrix(c(var.11.hat,cov.hat,cov.hat,var.12.hat),2))/I
draw.ellipse(mean.hat,sigma.hat.mc,col=3)
sigma.hat.mc <- with(data[1,],matrix(c(var.11.hat,cov.hat,cov.hat,var.12.hat),2))/I
draw.ellipse(mean.hat,sigma.hat.mc,col=4)


## 6c routine to get asy vcov matrix theta_11, theta_12. 
source('misc.R')
auc.cluster <- function(x,y,alpha=.05) {
    stopifnot(length(x)==length(y))
    I <- length(x)
    ## if(is.null(psi))  {
        psi <- outer(x,y,Vectorize(auc),normalized=FALSE)
        m <- sapply(x,length); n <- sapply(y,length)
        M <- sum(m); N <- sum(n)
    ## }
    phi <- diag(psi)/m/n
    theta.11.hat <- mean(diag(psi)/m/n)
    theta.12.hat <- (sum(psi)-sum(diag(psi)))/M/N * I/(I-1)
    ## theta.12.hat <- sum(psi)/M/N
    cov.hat <- mean(phi/mean(phi)*((colMeans(psi)+rowMeans(psi))/mean(psi) - m/mean(m)-n/mean(n))) * theta.11.hat*theta.12.hat
    var.11.hat <- var(phi)
    ## R <- (sum(psi)-sum(diag(psi)))/I/(I-1)
    ## var.12.hat <- (R/mean(m)/mean(n))^2*I*sum((1/R/I^2*(colSums(psi)+rowSums(psi)) - (m/M+n/N))^2) * (I/(I-1))^2
    var.12.hat <- theta.12.hat^2*I*sum((1/theta.12.hat/M/N*(colSums(psi)+rowSums(psi)) - (m/M+n/N))^2) * (I/(I-1))^2
    list(theta.11.hat=theta.11.hat,theta.12.hat=theta.12.hat,vcov.hat=matrix(c(var.11.hat,cov.hat,cov.hat,var.12.hat),2))
}


I <- 10
cov.xx <- .3; cov.xy <- .4; cov.yy <- .5
var.x <- 1; var.y <- 2; var.z <- 3
mean.x <- 0; mean.y <- 1
theta.11 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy))
theta.12 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z))
hats <- replicate(1e3, {
    m <- 1+rpois(I,5); n <- 1+rpois(I,5)
    data <- rbinormal(I=I,m=m,n=n,mean.x=mean.x,mean.y=mean.y,cov.xx=cov.xx,cov.xy=cov.xy,cov.yy=cov.yy,var.x=var.x,var.y=var.y,var.z=var.z)
    x <- data$x; y <- data$y
    ## psi <- outer(x,y,Vectorize(auc),normalized=FALSE)
    ## theta.12.hat <- (sum(psi)-sum(diag(psi))) / sum(m)/sum(n) * I/(I-1)
    ## theta.11.hat <- mean(diag(psi)/m/n)
    ## c(theta.11.hat=theta.11.hat,theta.12.hat=theta.12.hat)
    auc.cluster(x,y)
},simplify=FALSE)
theta.11.hats <- sapply(hats,function(lst)lst$theta.11.hat)
theta.12.hats <- sapply(hats,function(lst)lst$theta.12.hat)
vcov.hats <- simplify2array(lapply(hats,function(lst)lst$vcov.hat))
vcov.hat.mean <- apply(vcov.hats,1:2,mean)
c(var(theta.11.hats)*I,vcov.hat.mean[1,1])
c(var(theta.12.hats)*I,vcov.hat.mean[2,2])
c(cov(theta.11.hats,theta.12.hats)*I,vcov.hat.mean[1,2])

## 6d formulas for binormal model

require(mvtnorm)
I <- 5
cov.xx <- .3; cov.xy <- .4; cov.yy <- .5
var.x <- 1; var.y <- 2; var.z <- 3
mean.x <- 0; mean.y <- 1
theta.11 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy))
theta.12 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z))
hats <- replicate(1e3, {
    m <- 1+rpois(I,5); n <- 1+rpois(I,5)
    data <- lapply(1:I, function(i) {
        ## m <- ms[i]; n <- ns[i]
        mean.xy <- c(rep(mean.x,m[i]),rep(mean.y,n[i]))
        vcov.xy <- rbind(cbind(matrix(cov.xx,m[i],m[i]),matrix(cov.xy,m[i],n[i])),
                         cbind(matrix(cov.xy,n[i],m[i]),matrix(cov.yy,n[i],n[i]))) +
            diag(c(rep(var.x-cov.xx,m[i]),rep(var.y-cov.yy,n[i])))
        xy <- rmvnorm(1,mean=mean.xy,sigma=vcov.xy) + rnorm(1,sd=sqrt(var.z))
        x <- xy[1:(m[i])]; y <- xy[(m[i]+1):(m[i]+n[i])]
        list(x=x,y=y)
    })
    x <- lapply(data,function(xy)xy$x)
    y <- lapply(data,function(xy)xy$y)
    psi <- outer(x,y,Vectorize(auc),normalized=FALSE)
    theta.12.hat <- (sum(psi)-sum(diag(psi))) / sum(m)/sum(n) * I/(I-1)
    theta.11.hat <- mean(diag(psi)/m/n)
    c(theta.11.hat=theta.11.hat,theta.12.hat=theta.12.hat)
})
op <- par(mfrow=c(1,2))
hist(hats['theta.11.hat',]); abline(v=theta.11,col=2); abline(v=mean(hats['theta.11.hat',]),col=3)
hist(hats['theta.12.hat',]); abline(v=theta.12,col=2); abline(v=mean(hats['theta.12.hat',]),col=3)
par(op)
mean(hats['theta.11.hat',]) - theta.11
mean(hats['theta.12.hat',]) - theta.12

## 6d-1 routine to sample from binormal model with exchangeable vcov matrix with
## optional random effect
require(mvtnorm)
source('misc.R')
rbinormal <- function(I,m,n,mean.x,mean.y,cov.xx,cov.xy,cov.yy,var.x,var.y,var.z) {
    data <- lapply(1:I, function(i) {
        ## m <- ms[i]; n <- ns[i]
        mean.xy <- c(rep(mean.x,m[i]),rep(mean.y,n[i]))
        vcov.xy <- rbind(cbind(matrix(cov.xx,m[i],m[i]),matrix(cov.xy,m[i],n[i])),
                         cbind(matrix(cov.xy,n[i],m[i]),matrix(cov.yy,n[i],n[i]))) +
            diag(c(rep(var.x-cov.xx,m[i]),rep(var.y-cov.yy,n[i])))
        xy <- rmvnorm(1,mean=mean.xy,sigma=vcov.xy) + rnorm(1,sd=sqrt(var.z))
        x <- xy[1:(m[i])]; y <- xy[(m[i]+1):(m[i]+n[i])]
        list(x=x,y=y)
    })
    x <- lapply(data,function(xy)xy$x)
    y <- lapply(data,function(xy)xy$y)
    return(list(x=x,y=y))
}

I <- 5
cov.xx <- .3; cov.xy <- .4; cov.yy <- .5
var.x <- 1; var.y <- 2; var.z <- 3
mean.x <- 0; mean.y <- 1
theta.11 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy))
theta.12 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z))
hats <- replicate(5e2, {
    m <- 1+rpois(I,5); n <- 1+rpois(I,5)
    data <- rbinormal(I=I,m=m,n=n,mean.x=mean.x,mean.y=mean.y,cov.xx=cov.xx,cov.xy=cov.xy,cov.yy=cov.yy,var.x=var.x,var.y=var.y,var.z=var.z)
    x <- data$x; y <- data$y
    psi <- outer(x,y,Vectorize(auc),normalized=FALSE)
    theta.12.hat <- (sum(psi)-sum(diag(psi))) / sum(m)/sum(n) * I/(I-1)
    theta.11.hat <- mean(diag(psi)/m/n)
    c(theta.11.hat=theta.11.hat,theta.12.hat=theta.12.hat)
})
op <- par(mfrow=c(1,2))
hist(hats['theta.11.hat',]); abline(v=theta.11,col=2); abline(v=mean(hats['theta.11.hat',]),col=3)
hist(hats['theta.12.hat',]); abline(v=theta.12,col=2); abline(v=mean(hats['theta.12.hat',]),col=3)
par(op)
mean(hats['theta.11.hat',]) - theta.11
mean(hats['theta.12.hat',]) - theta.12


I <- 5
cov.xx <- .3; cov.xy <- .4; cov.yy <- .5
var.x <- 1; var.y <- 2; var.z <- 100
mean.x <- 0; mean.y <- 3
theta.11 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy))
theta.12 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z))
## hats <- replicate(5e2, {
m <- 1+rpois(I,5); n <- 1+rpois(I,5)
data <- rbinormal(I=I,m=m,n=n,mean.x=mean.x,mean.y=mean.y,cov.xx=cov.xx,cov.xy=cov.xy,cov.yy=cov.yy,var.x=var.x,var.y=var.y,var.z=var.z)
x <- data$x; y <- data$y
plot(0,type='n',xlim=range(unlist(data)),ylim=c(0,I))
sapply(1:I, function(i) points(c(data$x[[i]],data$y[[i]]),y=rep(i,m[i]+n[i]),col=rep(c(1,2),c(m[i],n[i]))))
points(c(unlist(data$x),unlist(data$y)),y=rep(0,sum(m)+sum(n)),col=rep(c(1,2),c(sum(m),sum(n))))

## 6d-2 formula to check if point lies in ellipse given parametrically
## as center + sigma * circle(theta)
source('misc.R')
alpha <- .05
q <- qchisq(1-alpha,df=2)
mean.xy <- c(.4,.6)
cov.xy <- matrix(runif(4),2); cov.xy <- cov.xy%*%t(cov.xy)
cov.xy.inv.sqrt <- with(eigen(solve(cov.xy)), vectors%*%diag(sqrt(values))%*%t(vectors))
axis <- 3
plot(0,type='n',xlim=c(-1,1)*axis,ylim=c(-1,1)*axis)
draw.ellipse(mean.xy,cov.xy,alpha=alpha)
reps <- 1e4
pt <- matrix(runif(2*reps,-axis,axis),nrow=2)
## idx <- colSums((cov.xy.inv.sqrt%*%(pt-mean.xy))^2) < q
## points(t(pt[,idx]))
inside.ellipse <- function(pt,center,qform,radius) {
    pt <- as.matrix(pt,nrow=2)
    qform.inv.sqrt <- with(eigen(solve(qform)), vectors%*%diag(sqrt(values))%*%t(vectors))
    colSums((qform.inv.sqrt%*%(pt-center))^2) < radius^2
}
idx <- inside.ellipse(pt, mean.xy, cov.xy, sqrt(q))
points(t(pt[,idx]))

## 6e type 1 error rate for H0: theta.11==theta.12, tested at null of (1/2+eps,1/2+eps)

## 6e-1 just fpr
require(mvtnorm)
source('misc.R')
alpha <- .05
I <- 10
cov.xx <- .3; cov.xy <- .2; cov.yy <- .3
var.x <- 1; var.y <- 2; var.z <- 3
mean.x <- 0; mean.y <- 0
theta.11 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy))
theta.12 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z))
reject <- replicate(1e3, {
    m <- 1+rpois(I,5); n <- 1+rpois(I,5)
    data <- rbinormal(I=I,m=m,n=n,mean.x=mean.x,mean.y=mean.y,cov.xx=cov.xx,cov.xy=cov.xy,cov.yy=cov.yy,var.x=var.x,var.y=var.y,var.z=var.z)
    x <- data$x; y <- data$y
    plot(0,type='n',xlim=c(0,1),ylim=c(0,1))
    contrast <- matrix(c(1,-1),nrow=2)
    z.stat <- with(auc.cluster(x,y), sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
    abs(z.stat) > qnorm(1-alpha/2)
})
mean(reject)

## 6e-2 fpr, bias, etc by I
require(mvtnorm)
require(parallel)
source('misc.R')
alpha <- .05
## I <- 10
Is <- c(10,30,100)
cov.xx <- .3; cov.xy <- .2; cov.yy <- .3
var.x <- 1; var.y <- 2; var.z <- 3
mean.x <- 0; mean.y <- 0
theta.11 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy))
theta.12 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z))
by.I <- mclapply(Is, mc.cores=detectCores()-3, FUN=function(I) {
    hats <- replicate(1e3, {
        m <- 1+rpois(I,5); n <- 1+rpois(I,5)
        data <- rbinormal(I=I,m=m,n=n,mean.x=mean.x,mean.y=mean.y,cov.xx=cov.xx,cov.xy=cov.xy,cov.yy=cov.yy,var.x=var.x,var.y=var.y,var.z=var.z)
        x <- data$x; y <- data$y
        ## z.stat <- with(auc.cluster(x,y), sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
        ## abs(z.stat) > qnorm(1-alpha/2)
        auc.cluster(x,y)
    },simplify=FALSE)
    contrast <- matrix(c(1,-1),nrow=2)
    var.hats <- sapply(hats, function(hat) with(hat, t(contrast)%*%vcov.hat%*%contrast))
    theta.11.hats <- sapply(hats, function(hat) hat$theta.11.hat)
    theta.12.hats <- sapply(hats, function(hat) hat$theta.12.hat)
    z.stats <- sapply(hats, function(hat) with(hat, sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast)))
    c(pt.11.bias=mean(theta.11.hats)-theta.11,pt.12.bias=mean(theta.12.hats)-theta.12,var.hat.norm=mean(var.hats),var.mc.norm=I*var(theta.11.hats-theta.12.hats),fpr=mean(abs(z.stats) > qnorm(1-alpha/2)))
})

by.I <- t(simplify2array(by.I))
## save.image('6e-2.RData')


dd

## 6f type 2 error rate


## 6f-1 spacing deltas to get more evenly spaced theta.11-theta.12
epsilon <- .05
var.x <- 1; var.y <- 1
cov.xx <- .3; cov.xy <- -.2; cov.yy <- .5
delta.null <- qnorm(1/2+epsilon)*sqrt(var.x+var.y-2*cov.xy)
deltas <- delta.null + seq(0,1,len=20)
## deltas <- deltas^(1/10)
## deltas <- qnorm(.5+epsilon+seq(0,.2,len=6))*sqrt(var.x+var.y-2*cov.xy)
plot(sapply(deltas,function(delta) {
    mean.y <- mean.x+delta
    var.z <- (1/qnorm(1/2+epsilon)^2*delta^2-var.x-var.y)/2
    theta.11 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy))
    theta.12 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z))
    theta.11-theta.12
}))

## 6f-2 null is theta11=theta12=1/2+epsilon
require(mvtnorm)
require(parallel)
source('misc.R')
alpha <- .05
epsilon <- .05
mean.x <- 0
var.x <- 1; var.y <- 1
cov.xx <- .3; cov.xy <- -.2; cov.yy <- .5
delta.null <- qnorm(1/2+epsilon)*sqrt(var.x+var.y-2*cov.xy)
deltas <- qnorm(.5+epsilon+seq(0,.2,len=6))*sqrt(var.x+var.y-2*cov.xy)
Is <- c(10,30,100)
by.I <- mclapply(Is, mc.cores=detectCores(), FUN=function(I) {
## by.I <- lapply(Is,function(I) {
    print(I)
    by.delta <- sapply(deltas, function(delta) {
        ## delta <- mean.y <- delta.null+.2
        mean.y <- mean.x+delta
        var.z <- (1/qnorm(1/2+epsilon)^2*delta^2-var.x-var.y)/2
        theta.11 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy))
        theta.12 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z))
        ## hats <- replicate(5e2, {
        ## m <- 1+rpois(I,5); n <- 1+rpois(I,5)
        ## data <- rbinormal(I=I,m=m,n=n,mean.x=mean.x,mean.y=mean.y,cov.xx=cov.xx,cov.xy=cov.xy,cov.yy=cov.yy,var.x=var.x,var.y=var.y,var.z=var.z)
        c(theta.11,theta.12)
        reject <- replicate(1e1, {
            m <- 1+rpois(I,5); n <- 1+rpois(I,5)
            data <- rbinormal(I=I,m=m,n=n,mean.x=mean.x,mean.y=mean.y,cov.xx=cov.xx,cov.xy=cov.xy,cov.yy=cov.yy,var.x=var.x,var.y=var.y,var.z=var.z)
            x <- data$x; y <- data$y
            ## plot(0,type='n',xlim=c(0,1),ylim=c(0,1))
            contrast <- matrix(c(1,-1),nrow=2)
            z.stat <- with(auc.cluster(x,y), sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
            abs(z.stat) > qnorm(1-alpha/2)
            ## auc.cluster(x,y)
        })
        ## vcov <- apply(simplify2array(reject[3,]),1:2,mean)
        ## c(var(unlist(reject['theta.11.hat',])),vcov[1,1]/I)
        c(theta.11=theta.11,theta.12=theta.12,power=mean(reject))
    })
})
op <- par(mfrow=c(1,3))
lapply(by.I, function(by.delta) {
    with(data.frame(t(by.delta)), plot(theta.11-theta.12, power,xlab=expression(paste(theta[11]-theta[12])), type='l',ylim=c(0,1)) )
})
par(op)
## save.image('6f-2.RData')

plot(0,type='n',xlim=range(simplify2array(by.I)['theta.11',,] - simplify2array(by.I)['theta.12',,]),ylim=c(0,1),xlab=expression(paste(theta[11]-theta[12])),ylab='power')
sapply(1:length(by.I), function(j) 
    with(data.frame(t(by.I[[j]])), lines(theta.11-theta.12, power, lty=j)) )
legend('topleft',lty=1:length(by.I),legend=paste0('I=',Is))

## plot sequence of alternatives

plot(0,type='n',xlim=c(0,1),ylim=c(0,1),xlab=expression(theta[11]),ylab=expression(theta[12]))
abline(0,1)
points(by.I[[1]]['theta.11',],by.I[[1]]['theta.12',])
text(by.I[[1]][1,1]-.03,by.I[[1]][2,1]+.03,expression(H[0]))

## region |y-1/2|>|x-1/2|
n <- 1e4
y <- runif(n)
x <- runif(n)
idx <- abs(y-1/2)>abs(x-1/2)
plot(0,type='n',xlim=c(0,1),ylim=c(0,1))
points(x[idx],y[idx])



## 6g power surface

## 6g-1 check grid of theta.11,theta.12 have only valid points under
## binormal model, vcov matrix for data is psd

## curve(pnorm(sqrt(2)*qnorm(x)),0,1) ## valid (theta_11,theta_12) for binormal

require(mvtnorm)
source('misc.R')
I <- 10
alpha <- .05
B <- 1e2
theta.11 <- seq(1/2,1,length.out=B)
theta.12 <- seq(0.05,.95,length.out=B)
## theta.12 <- runif(B,1-pnorm(sqrt(2)*qnorm(theta.11)),pnorm(sqrt(2)*qnorm(theta.11)))
mean.x <- 0
var.x <- var.y <- 5; var.z <- 0
## cov.xys <- -(var.x+var.y)/2*( (qnorm(theta.12)/qnorm(theta.11))^2*(1+2*var.z/(var.x+var.y)) - 1)
reject.rates <- matrix(NA,B,B)
for(i in 1:B)
    for(j in 1:B) {
        if(qnorm(theta.11[i])*qnorm(theta.12[j])<0) next
        if(theta.12[j] < 1-pnorm(sqrt(2)*qnorm(theta.11[i])) || theta.12[j] > pnorm(sqrt(2)*qnorm(theta.11[i]))) next
        c12 <- 2*var.z/(var.x+var.y)
        ## cov.xy <- -(var.x+var.y)/2*( (qnorm(theta.11[i])/qnorm(theta.12[j]))^2*(1+c12) - 1)
        c11 <- (1+c12)*(qnorm(theta.12[j])/qnorm(theta.11[i]))^2-1
        cov.xy <- -(var.x+var.y)/2*c11
        cov.xx <- cov.yy <- abs(cov.xy) ## ancillary but must be >= cov.xy for psd vcov
        ## print(cov.xy)
        ## if(is.nan(cov.xy)) browser()
        ## if(abs(cov.xy)>cov.xx) next
        if(abs(cov.xy)>mean(c(var.x,var.y))-.01) next
        ## c11 <- -2*cov.xy/(var.x+var.y)
        mean.y <- qnorm(theta.11[i])*sqrt(1+c11)*sqrt(var.x+var.y)
        delta <- (mean.y-mean.x)/sqrt(var.x+var.y)
        if(abs(theta.11[i]-pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy)))>1e-5)browser()
        if(abs(theta.12[j]-pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z)))>1e-5)
            browser()
        ## replicate(1e1, {
        ##     m <- 1+rpois(I,5); n <- 1+rpois(I,5)
        ##     sapply(1:I, function(k) {
        ##         vcov.xy <- rbind(cbind(matrix(cov.xx,m[k],m[k]),matrix(cov.xy,m[k],n[k])),
        ##                          cbind(matrix(cov.xy,n[k],m[k]),matrix(cov.yy,n[k],n[k]))) +
        ##             diag(c(rep(var.x-cov.xx,m[k]),rep(var.y-cov.yy,n[k])))
        ##         if(sum(eigen(vcov.xy)$val<0)>0) browser()
        ##         ## print(eigen(vcov.xy)$val)
        ##         ## rmvnorm(1,sigma=vcov.xy)
        ##     })
        ## })
        reject.rates[i,j] <- 1
    }
reject.rates[is.na(reject.rates)] <- 0
image(reject.rates)


theta.11 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy))
theta.12 <- pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z))

## 6g-2 (moved to 211120 sim)
require(mvtnorm)
source('misc.R')
I <- 10
alpha <- .05
B <- 1e1
theta.11 <- seq(1/2,1,length.out=B)
theta.12 <- seq(0.05,.95,length.out=B)
## theta.11 <- .7
## theta.12 <- .5
mean.x <- 0
var.x <- var.y <- 5; var.z <- 0
## reject.rates <- matrix(NA,B,B)
reject <- as.list(numeric(B^2))
dim(reject) <- c(B,B)
for(i in 1:B)
    for(j in 1:B) {
        reject[i,j] <- NA
        if(qnorm(theta.11[i])*qnorm(theta.12[j])<0) next
        if(theta.12[j] < 1-pnorm(sqrt(2)*qnorm(theta.11[i])) || theta.12[j] > pnorm(sqrt(2)*qnorm(theta.11[i]))) next
        c12 <- 2*var.z/(var.x+var.y)
        ## cov.xy <- -(var.x+var.y)/2*( (qnorm(theta.11[i])/qnorm(theta.12[j]))^2*(1+c12) - 1)
        c11 <- (1+c12)*(qnorm(theta.12[j])/qnorm(theta.11[i]))^2-1
        cov.xy <- -(var.x+var.y)/2*c11
        cov.xx <- cov.yy <- abs(cov.xy) ## ancillary but must be >= cov.xy for psd vcov
        ## print(cov.xy)
        ## if(is.nan(cov.xy)) browser()
        ## if(abs(cov.xy)>cov.xx) next
        if(abs(cov.xy)>mean(c(var.x,var.y))-.01) next
        ## c11 <- -2*cov.xy/(var.x+var.y)
        mean.y <- qnorm(theta.11[i])*sqrt(1+c11)*sqrt(var.x+var.y)
        delta <- (mean.y-mean.x)/sqrt(var.x+var.y)
        if(abs(theta.11[i]-pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy)))>1e-5)browser()
        if(abs(theta.12[j]-pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z)))>1e-5)
            browser()
        if(abs(cov.xy)>mean(c(var.x,var.y))-.01) next
        reject[[i,j]] <- replicate(2e0, {
            ## cat('.')
            m <- 1+rpois(I,5); n <- 1+rpois(I,5)
            data <- rbinormal(I=I,m=m,n=n,mean.x=mean.x,mean.y=mean.y,cov.xx=cov.xx,cov.xy=cov.xy,cov.yy=cov.yy,var.x=var.x,var.y=var.y,var.z=var.z,plot=TRUE)
            x <- data$x; y <- data$y
            ## plot(0,type='n',xlim=c(0,1),ylim=c(0,1))
            contrast <- matrix(c(1,-1),nrow=2)
            z.stat <- with(auc.cluster(x,y), sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
            abs(z.stat) > qnorm(1-alpha/2)
            ## auc.cluster(x,y)
        })
        ## reject.rates[i,j] <- mean(reject)
    }

power <- matrix(NA,B,B)
for(i in 1:B) for(j in 1:B) power[i,j] <- mean(reject[[i,j]])
image(theta.11,theta.12,power)
abline(a=0,b=1);abline(a=1,b=-1)

## image(x=theta.11,y=theta.12,z=reject.rates)




## 7 binormal model--checking formula for pop auc

I <- 5
max.cluster.size <- 10
m <- sample(1:max.cluster.size,I,replace=TRUE); n <- sample(1:max.cluster.size,I,replace=TRUE)
Delta <- 0
var.z <- 1
z <- rnorm(I)
x <- lapply(1:I, function(i) z[i]+  rnorm(m[i]))
y <- lapply(1:I, function(i) z[i]+  rnorm(n[i]))


## 7a binormal model -- following Obuchowski 97 sim
require(mvtnorm)
require(parallel)
source('misc.R')
I <- 1e2
alpha <- .05
k <- 10
## D.corr <- .4
## intra.corr <- .4
## theta.12 <- auc.pop <- .7
params <- expand.grid(list(D.corr=c(0,.4,.8,.95),intra.corr=c(0,.1,.4,.8),theta.12=c(.7,.8)))
by.param <- apply(params,1,function(r) {
    print(r)
    D.corr <- unname(r['D.corr']); intra.corr <- unname(r['intra.corr']); theta.12 <- unname(r['theta.12'])
    Delta <- sqrt(2)*qnorm(theta.12)
    theta.11 <- unname(pnorm(Delta / sqrt(2-2*intra.corr)))
    Sigma <- matrix(D.corr,k,k) + diag(k)*(1-D.corr)
    ## coverage <- replicate(3e2, {
    ## sim <- mclapply(1:3e2, mc.cores=detectCores()-3, FUN=function(jj){
        ## cat('.')
    sim <- lapply(1:3e1, FUN=function(jj){
        ## cat('.')
        D.latent <- rmvnorm(I,mean=rep(0,k),Sigma)
        D.idx <- D.latent>0
        m <- pmax(pmin(rowSums(D.idx),k-1),1) # enforce m>0,n>0
        n <- k-m
        data <- rbinormal(I=I,m=m,n=n,mean.x=0,mean.y=Delta,cov.xx=intra.corr,cov.xy=intra.corr,cov.yy=intra.corr)
        hat <- auc.cluster(data$x,data$y)
        q <- qchisq(1-alpha,df=2)
        ## if(sum(is.na(hat$vcov.hat))>0)
        ## browser()
        ## plot(0,type='n',xlim=c(0,1),ylim=c(0,1))
        ## draw.ellipse(   c(hat$theta.11.hat,hat$theta.12.hat), hat$vcov.hat/I,alpha=alpha)
        ## points(theta.11,theta.12)
        ## if(prod(eigen(hat$vcov.hat)$values)<.000000001)
            ## browser()
        cover <-  #tryCatch(
            inside.ellipse(pt=c(theta.11,theta.12),  center=c(hat$theta.11.hat,hat$theta.12.hat), qform=hat$vcov.hat/I, radius=sqrt(q))
           #,error=function(e){browser();print(1)})  ## small I--degenerate AUC
        list(theta.11.hat=hat$theta.11.hat,theta.12.hat=hat$theta.12.hat,cover=cover,vcov.hat=hat$vcov.hat)
        ## CI <- hat$theta.12.hat + c(-1,1)*sqrt(hat$vcov[2,2]/I)*qnorm(1-alpha/2)
        ## prod(CI-theta.12)<0
        ## CI <- hat$theta.11.hat + c(-1,1)*sqrt(hat$vcov[1,1]/I)*qnorm(1-alpha/2)
        ## prod(CI-theta.11)<0
        ## with(as.list(auc.obu(data$x,data$y,alpha=alpha)),prod(c(CI.lower,CI.upper)-theta.12)<0)
    })
    theta.11.hats <- sapply(sim,function(ll)ll$theta.11.hat)
    theta.12.hats <- sapply(sim,function(ll)ll$theta.12.hat)
    coverage <- mean(sapply(sim,function(ll)ll$cover))
    vcov.hats <- simplify2array(lapply(sim,function(ll)ll$vcov.hat))
    list(theta.11=theta.11,coverage=coverage,bias.theta.11=mean(theta.11.hats)-theta.11, bias.theta.12=mean(theta.12.hats)-theta.12,vcov.hat=apply(vcov.hats,1:2,mean),vcov.mc=matrix(c(var(theta.11.hats),rep(cov(theta.11.hats,theta.12.hats),2),var(theta.12.hats))*I,2))
})


out <- sapply(by.param, function(lst)unlist(lst)[c('theta.11','coverage','bias.theta.11','bias.theta.12','vcov.hat1','vcov.hat2','vcov.hat4','vcov.mc1','vcov.mc2','vcov.mc4')])
out <- cbind(params,t(out))





## 7b binormal model range of theta.11 values for each theta.12 value
## moved to manuscript.R


require(mvtnorm)
require(parallel)
source('misc.R')
I <- 1e2
alpha <- .05
k <- 10
## D.corr <- .4
## intra.corr <- .4
## theta.12 <- auc.pop <- .7
theta.12s <- c(.6,.8)
theta.11s <- sapply(theta.12s, function(theta.12)seq(theta.12,.95,len=3))
## theta.11s <- split(t(theta.11s),rep(1:ncol(theta.11s),nrow(theta.11s)))
aucs <- do.call(rbind,lapply(1:length(theta.12s),function(i)cbind(theta.12=theta.12s[i],theta.11=theta.11s[,i])))
params <- expand.grid(list(D.corr=c(0,.4,.8,.95)))
params <- as.data.frame(cbind(aucs[rep(1:nrow(aucs),each=nrow(params)),], params[rep(1:nrow(params),nrow(aucs)),,drop=FALSE]))
by.param <- apply(params,1,function(r) {
    print(r)
    D.corr <- unname(r['D.corr']);  theta.12 <- unname(r['theta.12']);  theta.11 <- unname(r['theta.11'])
    Delta <- sqrt(2)*qnorm(theta.12)
    rho <- 1 - 1/2*Delta^2 / qnorm(theta.11)^2
    rho <- round(rho,9)
    Sigma <- matrix(D.corr,k,k) + diag(k)*(1-D.corr)
    sim <- lapply(1:3e1, FUN=function(jj){
        ## cat('.')
        ## print(rho)
        D.latent <- rmvnorm(I,mean=rep(0,k),Sigma)
        D.idx <- D.latent>0
        m <- pmax(pmin(rowSums(D.idx),k-1),1) # enforce m>0,n>0
        n <- k-m
        data <- rbinormal(I=I,m=m,n=n,mean.x=0,mean.y=Delta,cov.xy=rho)
        hat <- auc.cluster(data$x,data$y)
        q <- qchisq(1-alpha,df=2)
        ## if(sum(is.na(hat$vcov.hat))>0)
        ## browser()
        ## plot(0,type='n',xlim=c(0,1),ylim=c(0,1))
        ## draw.ellipse(   c(hat$theta.11.hat,hat$theta.12.hat), hat$vcov.hat/I,alpha=alpha)
        ## points(theta.11,theta.12)
        ## if(prod(eigen(hat$vcov.hat)$values)<.000000001)
            ## browser()
        cover <-  #tryCatch(
            inside.ellipse(pt=c(theta.11,theta.12),  center=c(hat$theta.11.hat,hat$theta.12.hat), qform=hat$vcov.hat/I, radius=sqrt(q))
           #,error=function(e){browser();print(1)})  ## small I--degenerate AUC
        list(theta.11.hat=hat$theta.11.hat,theta.12.hat=hat$theta.12.hat,cover=cover,vcov.hat=hat$vcov.hat)
        ## CI <- hat$theta.12.hat + c(-1,1)*sqrt(hat$vcov[2,2]/I)*qnorm(1-alpha/2)
        ## prod(CI-theta.12)<0
        ## CI <- hat$theta.11.hat + c(-1,1)*sqrt(hat$vcov[1,1]/I)*qnorm(1-alpha/2)
        ## prod(CI-theta.11)<0
        ## with(as.list(auc.obu(data$x,data$y,alpha=alpha)),prod(c(CI.lower,CI.upper)-theta.12)<0)
    })
    theta.11.hats <- sapply(sim,function(ll)ll$theta.11.hat)
    theta.12.hats <- sapply(sim,function(ll)ll$theta.12.hat)
    coverage <- mean(sapply(sim,function(ll)ll$cover))
    vcov.hats <- simplify2array(lapply(sim,function(ll)ll$vcov.hat))
    list(coverage=coverage,bias.theta.11=mean(theta.11.hats)-theta.11, bias.theta.12=mean(theta.12.hats)-theta.12,vcov.hat=apply(vcov.hats,1:2,mean),vcov.mc=matrix(c(var(theta.11.hats),rep(cov(theta.11.hats,theta.12.hats),2),var(theta.12.hats))*I,2))
})


require(xtable)
## out <- sapply(by.param, function(lst)unlist(lst)[c('coverage','bias.theta.11','bias.theta.12','vcov.hat1','vcov.hat2','vcov.hat4','vcov.mc1','vcov.mc2','vcov.mc4')])
## out <- cbind(params,t(out))
out <- sapply(by.param, function(lst)with(lst,c(coverage=coverage,bias.theta.11=bias.theta.11,bias.theta.12=bias.theta.12,structure((vcov.hat-vcov.mc)[c(1,2,4)],names=c('vcov.11','vcov.12','vcov.22')))))
out <- cbind(params,t(out))









############## summer 2022


## 8 checking formulas for z+\xi model
n <- 1e4
beta <- .3
Z <- rnorm(n)
xi <- rnorm(n)
B <- Z+xi
D <- rbinom(n,1,prob=pnorm(beta*Z))
plot(ecdf(Z[D==0]))
lines(ecdf(Z[D==1]),col=2)

z0 <- 1
mean(Z[D==1]<=z0)
zz <- rnorm(1e4)
pnorm(z0)*mean(pnorm(beta*zz[zz<z0]))/mean(pnorm(zz))
mean(Z[D==0]<=z0)
zz <- rnorm(1e4)
pnorm(z0)*mean(1-pnorm(beta*zz[zz<z0]))/(1-mean(pnorm(zz)))





## 9. censored normal model

## 9a. individual auc formula
a  <- bnd <- abs(qnorm(.05))
n <- 1e6
rho <- .4
delta <- .7
rnorm.censored <- function(n,bnd,delta,rho) {
    x <- rnorm(n)
    y <- rnorm(n)
    y <- rho*x+sqrt(1-rho^2)*y+delta
    ## var(y)
    ## cor(x,y)
    x[x < -bnd] <- -bnd
    x[x > bnd] <- bnd
    y[y < -bnd] <- -bnd
    y[y > bnd] <- bnd
    return(data.frame(x=x,y=y))
}
xy <- rnorm.censored(n,a,delta,rho)
x <- xy$x; y <- xy$y
## mean(x<y)
f <- Vectorize(function(x,y)mvtnorm::dmvnorm(c(x,y),mean=c(0,delta),sigma=matrix(c(1,rho,rho,1),2)))
F <- function(lower,upper)mvtnorm::pmvnorm(lower,upper,mean=c(0,delta),sigma=matrix(c(1,rho,rho,1),2))
interior <- integrate(Vectorize(function(x)integrate(function(y)f(x,y),x,a)$val),-a,a)$val
bnd  <- pnorm(-a) + 1-pnorm(a-delta) - F(c(-Inf,-Inf),c(-a,-a)) - F(c(a,a),c(Inf,Inf)) - F(c(-Inf,a),c(-a,Inf))
## mean(x<y)
## interior + bnd
diag <- F(c(-Inf,-Inf),c(-a,-a)) + F(c(a,a),c(Inf,Inf))
mean(x<y) + 1/2*mean(x==y)
interior + bnd + 1/2*diag

auc.norm.censored <- function(bnd,delta,rho) {
    f <- Vectorize(function(x,y)mvtnorm::dmvnorm(c(x,y),mean=c(0,delta),sigma=matrix(c(1,rho,rho,1),2)))
    F <- function(lower,upper)mvtnorm::pmvnorm(lower,upper,mean=c(0,delta),sigma=matrix(c(1,rho,rho,1),2))
    interior <- integrate(Vectorize(function(x)integrate(function(y)f(x,y),x,bnd)$val),-bnd,bnd)$val
    boundary  <- pnorm(-bnd) + 1-pnorm(bnd-delta) - F(c(-Inf,-Inf),c(-bnd,-bnd)) - F(c(bnd,bnd),c(Inf,Inf)) - F(c(-Inf,bnd),c(-bnd,Inf))
    diag <- F(c(-Inf,-Inf),c(-bnd,-bnd)) + F(c(bnd,bnd),c(Inf,Inf))
    interior + boundary + 1/2*diag
    ## integrate(Vectorize(function(x)integrate(function(y)f(x,y),x,bnd)$val),-bnd,bnd)$val + pnorm(-bnd) + 1-pnorm(bnd-delta) - 1/2*F(c(-Inf,-Inf),c(-bnd,-bnd)) - 1/2*F(c(bnd,bnd),c(Inf,Inf)) - F(c(-Inf,bnd),c(-bnd,Inf))
}

## check auc.norm.censored
n <- 1e4
pairs <- replicate(1e2, {
    rho <- runif(1,-1,1)
    bnd <- abs(qnorm(runif(1)/2))
    delta <- runif(1)*2
    formula <- auc.norm.censored(bnd,delta,rho)
    xy <- rnorm.censored(n,bnd,delta,rho)
    x <- xy$x; y <- xy$y
    observed <- mean(x<y) + 1/2*mean(x==y)
    c(formula=formula,observed=observed)
})
plot(pairs['formula',],pairs['observed',]);abline(0,1)



auc.norm.censored(qnorm(.0),.7,.4)
pnorm(.7/sqrt(2*(1-.4)))

rho <- -.4
bnd <- qnorm(.0)
deltas <- seq(0,5,length=30)
aucs <- sapply(deltas, function(delta)auc.norm.censored(bnd,delta,rho=rho))
plot(deltas,aucs)





rho <- .4
auc <- .7
uniroot(Vectorize(function(delta)auc.norm.censored(bnd=qnorm(.2),delta=delta,rho=rho)-.7),interval=c(0,.5),extendInt='yes')



## 9b population AUC formula

## rnorm.censored <- function(n,bnd,delta,rho) {
rnorm.censored <- function(n,bnd,delta,rho) {
    stopifnot(bnd>=0)
    x <- rnorm(n)
    y <- rnorm(n)
    y <- rho*x+sqrt(1-rho^2)*y+delta
    ## var(y)
    ## cor(x,y)
    x[x < -bnd] <- -bnd
    x[x > bnd] <- bnd
    y[y < -bnd] <- -bnd
    y[y > bnd] <- bnd
    return(data.frame(x=x,y=y))
}

auc.norm.censored <- function(bnd,delta,rho) {
    if(rho==0)return (
                  -integrate(function(x)dnorm(x)*pnorm(x-delta),-bnd,bnd)$val + 1/2*(pnorm(bnd)-pnorm(bnd-delta)-pnorm(-bnd-delta)+pnorm(bnd)*(pnorm(-bnd-delta)+pnorm(bnd-delta))+1)
              )
    f <- Vectorize(function(x,y)mvtnorm::dmvnorm(c(x,y),mean=c(0,delta),sigma=matrix(c(1,rho,rho,1),2)))
    F <- function(lower,upper)mvtnorm::pmvnorm(lower,upper,mean=c(0,delta),sigma=matrix(c(1,rho,rho,1),2))
    ## integrate(Vectorize(function(x)integrate(function(y)f(x,y),x,bnd)$val),-bnd,bnd)$val + pnorm(-bnd) + 1-pnorm(bnd-delta) - 1/2*F(c(-Inf,-Inf),c(-bnd,-bnd)) - 1/2*F(c(bnd,bnd),c(Inf,Inf)) - F(c(-Inf,bnd),c(-bnd,Inf))
    interior <- integrate(Vectorize(function(x)integrate(function(y)f(x,y),x,bnd)$val),-bnd,bnd)$val
    boundary  <- pnorm(-bnd) + 1-pnorm(bnd-delta) - F(c(-Inf,-Inf),c(-bnd,-bnd)) - F(c(bnd,bnd),c(Inf,Inf)) - F(c(-Inf,bnd),c(-bnd,Inf))
    diag <- F(c(-Inf,-Inf),c(-bnd,-bnd)) + F(c(bnd,bnd),c(Inf,Inf))
    interior + boundary + 1/2*diag
}


n <- 1e3
bnd <- abs(pnorm(.2))
rho <- 0
delta <- .5
xy <- rnorm.censored(n,a,delta,rho)
x <- xy$x; y <- xy$y
auc.norm.censored(bnd,delta,rho)
integrate(function(x)dnorm(x)*(pnorm(bnd-delta)-pnorm(x-delta)),-bnd,bnd)$value + pnorm(-bnd)+1-pnorm(bnd-delta)-pnorm(-bnd)*(1-pnorm(bnd-delta))-1/2*pnorm(-bnd)*pnorm(-bnd-delta)-1/2*(1-pnorm(bnd))*(1-pnorm(bnd-delta))
-integrate(function(x)dnorm(x)*pnorm(x-delta),-bnd,bnd)$val + 1/2*(pnorm(bnd)-pnorm(bnd-delta)-pnorm(-bnd-delta)+pnorm(bnd)*(pnorm(-bnd-delta)+pnorm(bnd-delta))+1)

a <- runif(1)
delta <- runif(1)
-integrate(function(x)dnorm(x-delta)*dnorm(x),-a,a)$val
-1/2/sqrt(pi)*exp(-delta^2/4)*(pnorm(-delta/sqrt(2)+a*sqrt(2))-pnorm(-delta/sqrt(2)-a*sqrt(2)))

bnd <- abs(qnorm(.1))
deltas <- seq(0,4,len=10)
aucs <- sapply(deltas, function(delta)
    c(pop=auc.norm.censored(bnd,delta,0),indiv=auc.norm.censored(bnd,delta,.9)))
plot(deltas,aucs['pop',],type='l',ylim=c(0,1))
lines(deltas,aucs['indiv',],col=2)

auc.indiv <- .6
bnd=qnorm(.1)
auc.pop=.4
obj <- function(par,...) {
    delta <- par[1]; rho <- par[2]
    (auc.norm.censored(bnd,delta,rho) - auc.indiv)^2 + (auc.norm.censored(bnd,delta,0) - auc.pop)^2
}

optim0 <- optim(par=c(.1,.2),fn=obj)

auc.pop <- .6
delta.star <- uniroot(function(delta)auc.norm.censored(bnd,delta,rho=0) - auc.pop,c(0,2),extendInt='yes')$root

auc.norm.censored(bnd,delta=delta.star,rho=.95) - auc.indiv

auc.indiv <- .8
obj <- Vectorize(function(rho)abs(auc.norm.censored(bnd,delta.star,rho)-auc.indiv))
## rho <- uniroot(f=f,interval=c(.9,.95))$root
## uniroot(function(rho)auc.norm.censored(bnd,delta.star,rho=rho) - auc.indiv,c(.7,.8),extendInt='upX')$root
## rho.star <- optim(par=.5,fn=obj)$par
rho.star <- optimize(f=obj,interval=c(0,1))$min

aucs.to.params <- function(theta.12,theta.11) {
    auc.pop <- theta.12; auc.indiv <- theta.11
    delta.star <- uniroot(function(delta)auc.norm.censored(bnd,delta,rho=0) - auc.pop,c(0,2),extendInt='yes')$root
rho.star <- optimize(f=obj,interval=c(0,1))$min
return(c(delta=delta.star,rho=rho.star))
}

aucs.to.params(auc.pop,auc.indiv)


I <- 1e2
alpha <- .05
k <- 10
## D.corr <- .4
## intra.corr <- .4
## theta.12 <- auc.pop <- .7
theta.12s <- c(.6,.8)
theta.11s <- sapply(theta.12s, function(theta.12)seq(theta.12,.95,len=3))
## theta.11s <- split(t(theta.11s),rep(1:ncol(theta.11s),nrow(theta.11s)))
aucs <- do.call(rbind,lapply(1:length(theta.12s),function(i)cbind(theta.12=theta.12s[i],theta.11=theta.11s[,i])))
params <- expand.grid(list(D.corr=c(0,.4,.8,.95)))
params <- as.data.frame(cbind(aucs[rep(1:nrow(aucs),each=nrow(params)),], params[rep(1:nrow(params),nrow(aucs)),,drop=FALSE]))


apply(params, 1, function(r) aucs.to.params(r['theta.12'],r['theta.11']))

## 9c check data generated under rbinormal.censored auc against theoretical
## population auc
source('misc.R')
I <- 1e2
m <- n <- rep(10,I)
pairs <- replicate(1e2, {
    rho <- runif(1)
    ## m <- sample(1:10,I,replace=TRUE)
    ## n <- sample(1:10,I,replace=TRUE)
    mean.y <- runif(1)
    bound <- abs(qnorm(.2))
    xy <- rbinormal.censored(I=I,m=m,n=n,bound=bound,mean.y=mean.y,cov.xy=rho, cov.xx=rho,cov.yy=rho)
    x <- xy$x; y <- xy$y
    obs <- auc.cluster(x,y,get.vcov=FALSE)['theta.12.hat']
    ## x <- unlist(x); unlist(y) <- unlist(unlist(y))
    theoretical <- auc.norm.censored(bound,mean.y,rho=0)
    c(obs=unname(obs),theoretical=theoretical)
})
plot(pairs['obs',],pairs['theoretical',])
abline(0,1)

## mean(outer(unlist(x),unlist(y),`<`)) + 1/2*mean(outer(unlist(x),unlist(y),`==`))

## xy <- rnorm.censored(1e4,bound,mean.y,rho)
## x <- xy$x; y <- xy$y
## mean(x<y)+1/2*mean(x==y)
## auc.norm.censored(bound,mean.y,rho=0)


## individual--setting rho=rho rather than rho=0 in call to
## auc.norm.censored for theoretical theta.11
source('misc.R')
I <- 1e2
## m <- n <- rep(10,I)
pairs <- replicate(1e2, {
    m <- sample(1:10,I,replace=TRUE)
    n <- sample(1:10,I,replace=TRUE)
    rho <- runif(1)
    mean.y <- runif(1)
    bound <- abs(qnorm(.2))
    xy <- rbinormal.censored(I=I,m=m,n=n,bound=bound,mean.y=mean.y,cov.xy=rho, cov.xx=rho,cov.yy=rho)
    x <- xy$x; y <- xy$y
    obs <- auc.cluster(x,y,get.vcov=FALSE)['theta.11.hat']
    ## x <- unlist(x); unlist(y) <- unlist(unlist(y))
    theoretical <- auc.norm.censored(bound,mean.y,rho=rho)
    c(obs=unname(obs),theoretical=theoretical)    
})
plot(pairs['obs',],pairs['theoretical',])
abline(0,1)



## 9d coverage is poor for theta11=.95. check how obu estimator
## fares. result: obu is fine, but the max theta.12 or theta.11 is <1 for the censored binormal model.
source('misc.R')
I <- 1e1
m  <- n <- replicate(10,I)
bound <- abs(qnorm(.2))
theta.12 <- .93; theta.11 <- .95
f <- Vectorize(function(x)auc.binormal.censored(bound,x,rho=.5) )
curve(f, 0,8)
abline(h=1/2*(pnorm(bound)+1))


Delta.star <- uniroot(function(Delta)auc.binormal.censored(bound,Delta,rho=0) - theta.12,c(0,5),extendInt='yes')$root

Delta.rho <- auc.to.params.binormal.censored(theta.12,theta.11,bound)
Delta <- Delta.rho['Delta']; rho <- Delta.rho['rho']

coverage <- replicate(5e2, {
xy <- rbinormal.censored(I=I,m=m,n=n,bound=bound,mean.y=Delta,cov.xy=rho)
x <- xy$x; y <- xy$y
with(as.list(auc.obu(x,y)), theta.12 > CI.lower & theta.12 < CI.upper)
})
mean(coverage)



n <- 1e4
pairs <- replicate(1e3, {
    delta <- runif(1)
    rho <- runif(1,-1,1)
    x <- rnorm(n)
    y <- rnorm(n)
    y <- rho*x+sqrt(1-rho^2)*y + delta
    c(mean(x<y),pnorm(delta/sqrt(2*(1-rho))))
})
plot(pairs[1,],pairs[2,]);abline(0,1)


## 10 data analysis

## 10a obuchowski 97 data. personalized auc only defined for the "left"
## diseased and "right" diseased arteries data. outcome: theta.11 and
## theta.12 not significantly different for either reader.

source('misc.R')
mra <- as.numeric(scan('obuchowski1997_table3.txt',what=''))
top <- mra[1:(16*6)]
top <- t(array(top,c(16,6)))
bot <- mra[(16*6+1):(length(mra)-8)]
bot <- t(array(bot,c(8,5)))
both <- top[,1:4]
left <- top[,5:8]
right <- rbind(top[,9:12],bot[,1:4])
neither <- rbind(top[,13:16],bot[,5:8],matrix(mra[(length(mra)-7):length(mra)],nrow=2,byrow=TRUE))
colnames(both) <- colnames(left) <- colnames(right) <- colnames(neither) <- c('left.1','right.1','left.2','right.2')
mra.1 <- rbind(
    rbind(cbind(both[,'left.1'],1),cbind(both[,'right.1'],1)),
    rbind(cbind(left[,'left.1'],1),cbind(left[,'right.1'],0)),
    rbind(cbind(right[,'left.1'],0),cbind(right[,'right.1'],1)),
    rbind(cbind(neither[,'left.1'],0),cbind(neither[,'right.1'],0)))
mra.2 <- rbind(
    rbind(cbind(both[,'left.2'],1),cbind(both[,'right.2'],1)),
    rbind(cbind(left[,'left.2'],1),cbind(left[,'right.2'],0)),
    rbind(cbind(right[,'left.2'],0),cbind(right[,'right.2'],1)),
    rbind(cbind(neither[,'left.2'],0),cbind(neither[,'right.2'],0)))

## put in format suitable for auc.cluster
mra.1 <- rbind(left[,c('right.1','left.1')], right[,c('left.1','right.1')])
mra.2 <- rbind(left[,c('right.2','left.2')], right[,c('left.2','right.2')])
mra.1 <- na.omit(mra.1)
mra.2 <- na.omit(mra.2)
colnames(mra.1) <- colnames(mra.2) <- c('control','case')
auc.1 <- auc.cluster(x=mra.1[,'control'],y=mra.1[,'case'])
auc.2 <- auc.cluster(x=mra.2[,'control'],y=mra.2[,'case'])

contrast <- matrix(c(1,-1),nrow=2)
alpha <- .05
I <- nrow(mra.1)
z.stat <- with(auc.1, sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
abs(z.stat) > qnorm(1-alpha/2)

contrast <- matrix(c(1,-1),nrow=2)
alpha <- .05
I <- nrow(mra.2)
z.stat <- with(auc.2, sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
abs(z.stat) > qnorm(1-alpha/2)



## progabide data. outcome: no good--clusters are entirely case or control.
load('progabide.RData')


## 10b blip data. outcome: test that two aucs differ has p-val a little over 5% for biomarker CD4. around 6--7% for CD4P.
source('misc.R')
blip <- read.csv('data/final.csv')
blip$blip <- as.numeric(blip$VL>1e3)
blip <- subset(blip,select=c('ID','Blip_YN50','CD4','CD4P','Visits','blip','Sex','Transmission','trtm','AgeDiag','TimeDiag','AGE','DOV','TIME'))
blip <- subset(blip,select=c('ID','CD4','CD4P','Visits','blip','Blip_YN50','TIME'))
blip <- dplyr::rename(blip,Visit=Visits)
blip$CD4[blip$CD4>6000] <- NA
blip$CD4P[blip$CD4P>=99] <- NA
blip$Visit <- unlist(sapply(rle(blip$ID)$lengths,function(x)seq_len(x)))
## blip$blip <- factor(blip$blip)
blip$ID <- factor(blip$ID)
blip <- droplevels(blip)

clusters <- split(blip,blip$ID)
M <- sapply(clusters,function(cluster)sum(cluster$blip==0,na.rm=TRUE))
N <- sapply(clusters,function(cluster)sum(cluster$blip==1,na.rm=TRUE))
clusters <- clusters[M>0 & N>0]
x <- lapply(clusters,function(cluster)na.omit(cluster$CD4[cluster$blip==0]))
y <- lapply(clusters,function(cluster)na.omit(cluster$CD4[cluster$blip==1]))
(( auc0 <- auc.cluster(x=x,y=y) ))

contrast <- matrix(c(1,-1),nrow=2)
alpha <- .05
I <- length(x)
z.stat <- with(auc0, sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
1-pnorm(z.stat)
abs(z.stat) > qnorm(1-alpha/2)



## 10c nypd stop and frisk data 2021. theta11 and theta12 significantly
## different when predicting black/non-black by length of stop. need to find meaning of officer_command_code (used for clustering). direction of auc is negative: black status associated with shorter durations. 
source('misc.R')
snf <- read.csv('data/sqf-2021.csv')
snf$SUSPECT_REPORTED_AGE_FLAG <- as.numeric(snf$SUSPECT_REPORTED_AGE)>20
snf$race <- snf$SUSPECT_RACE_DESCRIPTION=='BLACK'#'BLACK HISPANIC'##BLACK, WHITE HISPANIC
clusters <- split(snf,snf$ISSUING_OFFICER_COMMAND_CODE)
biomarker <- 'SUSPECT_WEIGHT'
biomarker <- 'SUSPECT_REPORTED_AGE_FLAG'
biomarker <- 'STOP_DURATION_MINUTES'
## biomarker <- 'OBSERVED_DURATION_MINUTES'
## status <- 'SUSPECT_ARRESTED_FLAG'
## x <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=SUSPECT_ARRESTED_FLAG=='N',select=biomarker)))))
## y <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=SUSPECT_ARRESTED_FLAG=='Y',select=biomarker)))))
x <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=race==FALSE,select=biomarker)))))
y <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=race==TRUE,select=biomarker)))))
## x <- lapply(clusters, function(cluster)na.omit(cluster$SUSPECT_HEIGHT[cluster$SUSPECT_ARRESTED_FLAG=='N']))
## y <- lapply(clusters, function(cluster)na.omit(cluster$SUSPECT_HEIGHT[cluster$SUSPECT_ARRESTED_FLAG=='Y']))
omit.idx <- sapply(x,length)==0 | sapply(y,length)==0
x <- x[!omit.idx]; y <- y[!omit.idx]
(( auc0 <- auc.cluster(x=x,y=y) ))
contrast <- matrix(c(1,-1),nrow=2)
alpha <- .05
I <- length(x)
z.stat <- with(auc0, sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
(1-pnorm(abs(z.stat)))*2
abs(z.stat) > qnorm(1-alpha/2)
plot(auc0,resol=1e2)



aggregate(STOP_DURATION_MINUTES ~ SUSPECT_RACE_DESCRIPTION, mean, data=snf)

snf$STOP_DURATION_MINUTES[snf$STOP_DURATION_MINUTES>60] <- NA
hist(snf$STOP_DURATION_MINUTES[snf$race==TRUE])
hist(snf$STOP_DURATION_MINUTES[snf$race==FALSE],add=TRUE,col=2)


## 10cc nypd data--looking at several years data
source('misc.R')
filelist <- dir('data')
filelist <- filelist[grep('^sqf[-0-9]+\\.csv',filelist)]
filelist <- paste0('data/',filelist)
## snf <- lapply(filelist[3:4], read.csv)
snf <- lapply(filelist[], function(file) {
    with(read.csv(file),
         data.frame(STOP_DURATION_MINUTES,SUSPECT_RACE_DESCRIPTION,ISSUING_OFFICER_COMMAND_CODE))
})
snf <- do.call(rbind,snf)
snf$race <- snf$SUSPECT_RACE_DESCRIPTION=='WHITE' #| snf$SUSPECT_RACE_DESCRIPTION=='WHITE HISPANIC'#BLACK HISPANIC'##BLACK, WHITE HISPANIC
## snf$race <- snf$SUSPECT_RACE_DESCRIPTION=='BLACK HISPANIC' | snf$SUSPECT_RACE_DESCRIPTION=='BLACK'
## snf$race <- snf$SUSPECT_RACE_DESCRIPTION=='BLACK' | snf$SUSPECT_RACE_DESCRIPTION=='BLACK HISPANIC'
snf$race <- snf$SUSPECT_RACE_DESCRIPTION=='BLACK HISPANIC' | snf$SUSPECT_RACE_DESCRIPTION=='WHITE HISPANIC'
clusters <- split(snf,snf$ISSUING_OFFICER_COMMAND_CODE)
biomarker <- 'SUSPECT_WEIGHT'
biomarker <- 'SUSPECT_REPORTED_AGE_FLAG'
biomarker <- 'STOP_DURATION_MINUTES'
## biomarker <- 'OBSERVED_DURATION_MINUTES'
## status <- 'SUSPECT_ARRESTED_FLAG'
## x <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=SUSPECT_ARRESTED_FLAG=='N',select=biomarker)))))
## y <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=SUSPECT_ARRESTED_FLAG=='Y',select=biomarker)))))
x <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=race==FALSE,select=biomarker)))))
y <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=race==TRUE,select=biomarker)))))
## x <- lapply(clusters, function(cluster)na.omit(cluster$SUSPECT_HEIGHT[cluster$SUSPECT_ARRESTED_FLAG=='N']))
## y <- lapply(clusters, function(cluster)na.omit(cluster$SUSPECT_HEIGHT[cluster$SUSPECT_ARRESTED_FLAG=='Y']))
omit.idx <- sapply(x,length)==0 | sapply(y,length)==0
x <- x[!omit.idx]; y <- y[!omit.idx]
(( auc0 <- auc.cluster(x=x,y=y) ))
contrast <- matrix(c(1,-1),nrow=2)
alpha <- .05
I <- length(x)
z.stat <- with(auc0, sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
(1-pnorm(abs(z.stat)))*2
abs(z.stat) > qnorm(1-alpha/2)
plot(auc0,resol=1e2)



## 10d boston stop n frisk data.
source('misc.R')
years <- 2019:2021 # stop_duration coded differently for year<=2018
snfs <- lapply(years, function(year) {
    ## if(year==2017)browser()
    snf <- read.csv(paste0('data/fio_',year,'.csv'))
    snf.names <- read.csv(paste0('data/fio_names_',year,'.csv'))
    snf <- merge(snf,snf.names,by='fc_num')
    ## browser()
    ## print(year)
    ## print(table(snf$stop_duration))
    snf$biomarker <- snf$stop_duration
    snf$biomarker[snf$biomarker=='NULL'] <- NA
    snf$biomarker <- as.numeric(snf$biomarker)
    snf$status <- snf$race=='White'  & snf$ethnicity=='Not of Hispanic Origin'#& snf$ethnicity!='Hispanic Origin'
    snf$status <- snf$race=='White'  & snf$ethnicity!='Hispanic Origin'
    snf$status <- snf$ethnicity=='Hispanic Origin'
    subset(snf,select=c(fc_num,status,biomarker,race,contact_officer))
})
snf <- do.call(rbind,snfs)
snf$biomarker[snf$biomarker>200] <- NA
snf <- snf[order(snf$fc_num),]
single.encounters <- with(rle(snf$fc_num), values[lengths==1])
snf <- snf[snf$fc_num %in% single.encounters,]
clusters <- split(snf,snf$contact_officer)
x <- lapply(clusters, function(cluster)na.omit(cluster$biomarker[cluster$status==FALSE]))
y <- lapply(clusters, function(cluster)na.omit(cluster$biomarker[cluster$status==TRUE]))
omit.idx <- sapply(x,length)==0 | sapply(y,length)==0
x <- x[!omit.idx]; y <- y[!omit.idx]
(( auc0 <- auc.cluster(x=x,y=y) ))
contrast <- matrix(c(1,-1),nrow=2)
alpha <- .05
stopifnot(length(x)==length(y))
I <- length(x)
z.stat <- with(auc0, sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
(1-pnorm(abs(z.stat)))*2
abs(z.stat) > qnorm(1-alpha/2)
plot(auc0)
plot(auc0,alpha=.01,add=TRUE)


aggregate(biomarker ~ race, data=snf, FUN=mean)
aggregate(biomarker ~ race, data=snf, FUN=sd)

## little correlation bw means and cluster sizes
lengths <- sapply(x,length)+sapply(y,length)
means <- (sapply(x,sum)+sapply(y,sum)) / lengths
