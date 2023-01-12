## 1. figures showing large indiv auc vs small pop auc and vice versa


## 1a. large indiv auc, small pop auc

require(mvtnorm)
source('../misc.R')
I <- 10
max.cluster.size <- 20
m <- sample(1:max.cluster.size,I,replace=TRUE); n <- sample(1:max.cluster.size,I,replace=TRUE)
Delta <- 2
var.z <- 1e2
data <- rbinormal(I=I,m=m,n=n,mean.x=0,mean.y=Delta,cov.xx=0,cov.xy=0,cov.yy=0,var.x=1,var.y=1,mean.z=rep(0,I),var.z=var.z,plot=FALSE)
x <- data$x; y <- data$y
auc.cluster(x,y,get.vcov=FALSE)
## save.image('211223a.RData')

pch.control <- '|'
pch.case <- '-'
plot(0,type='n',xlim=range(unlist(c(x,y))),ylim=c(0,I))
for(i in 1:I) points(x[[i]],y=rep(i,length(x[[i]])),pch=pch.control)
for(i in 1:I) points(y[[i]],y=rep(i,length(y[[i]])),pch=pch.case)
points(unlist(c(x,y)),y=rep(0,sum(m)+sum(n)),pch=rep(c(pch.control,pch.case),c(sum(m),sum(n))))


## 1b. large pop auc, small indiv auc

require(mvtnorm)
source('../misc.R')
I <- 15
max.cluster.size <- 10
## m <- sample(1:max.cluster.size,I,replace=TRUE); n <- sample(1:max.cluster.size,I,replace=TRUE)
m <- rpois(I,1)+1
m <- m + rbinom(I,1,1/2)*(20-2*m)
n <- 20-m
Delta <- 0
var.z <- 1
data <- rbinormal(I=I,m=m,n=n,mean.x=0,mean.y=Delta,cov.xx=0,cov.xy=0,cov.yy=0,var.x=1,var.y=1,mean.z=(n*(n>m)-m*(n<=m))/5,var.z=var.z,plot=FALSE)
x <- data$x; y <- data$y
auc.cluster(x,y,get.vcov=FALSE)
## save.image('211223b.RData')


## 1c. plotting routine after 1a and 1b have been run

## op <- par(mfrow=c(1,2))
rm(list=ls())
for(save.file in c('211223a','211223b')) {
    load(paste0(save.file,'.RData'))
    pch.control <- '|'
    pch.case <- '-'
    png(paste0(save.file,'.png'))
    plot(0,type='n',xlim=range(unlist(c(x,y))),ylim=c(0,I),xlab='predictor',ylab='clusters',xaxt='n',yaxt='n',mgp=c(1,1,0))   
    for(i in 1:I) points(x[[i]],y=rep(i,length(x[[i]])),pch=pch.control)
    for(i in 1:I) points(y[[i]],y=rep(i,length(y[[i]])),pch=pch.case)
    abline(h=1/2)
    points(unlist(c(x,y)),y=rep(0,sum(m)+sum(n)),pch=rep(c(pch.control,pch.case),c(sum(m),sum(n))))
    dev.off()
}
## par(op)

## 1d. check of formula for theta_12 under location shift, RE model
## I <- 1e4
## m1 <- sample(1:max.cluster.size,I,replace=TRUE); n1 <- sample(1:max.cluster.size,I,replace=TRUE)
## z1 <- rnorm(I,n1-m1)
## m2 <- sample(1:max.cluster.size,I,replace=TRUE); n2 <- sample(1:max.cluster.size,I,replace=TRUE)
## z2 <- rnorm(I,n2-m2)
## mean(m1*n2*(z1-z2 +2*rnorm(I)*0< Delta)) / (mean(m)*mean(n))








## 2. figures showing large indiv auc vs small pop auc and vice versa
## re-doing #1 without random effect. also using Obuchowski '97 method
## of generating correlated m,n.


## 2a. large indiv auc, small pop auc

require(mvtnorm)
source('../misc.R')
I <- 10
max.cluster.size <- 20
m <- sample(1:max.cluster.size,I,replace=TRUE); n <- sample(1:max.cluster.size,I,replace=TRUE)
Delta <- .3
cc <- .99
## data <- rbinormal(I=I,m=m,n=n,mean.x=0,mean.y=Delta,cov.xx=0,cov.xy=0,cov.yy=0,var.x=1,var.y=1,mean.z=rep(0,I),var.z=var.z,plot=FALSE)
data <- rbinormal(I=I,m=m,n=n,mean.x=0,mean.y=Delta,cov.xx=cc,cov.xy=cc,cov.yy=cc,var.x=1,var.y=1,plot=FALSE)
x <- data$x; y <- data$y
auc.cluster(x,y,get.vcov=FALSE)
## save.image('211231a.RData')

pch.control <- '|'
pch.case <- '-'
plot(0,type='n',xlim=range(unlist(c(x,y))),ylim=c(0,I))
for(i in 1:I) points(x[[i]],y=rep(i,length(x[[i]])),pch=pch.control)
for(i in 1:I) points(y[[i]],y=rep(i,length(y[[i]])),pch=pch.case)
points(unlist(c(x,y)),y=rep(0,sum(m)+sum(n)),pch=rep(c(pch.control,pch.case),c(sum(m),sum(n))))


## 2b. large pop auc, small indiv auc

require(mvtnorm)
source('../misc.R')
I <- 10
max.cluster.size <- k <- 10
## m <- sample(1:max.cluster.size,I,replace=TRUE); n <- sample(1:max.cluster.size,I,replace=TRUE)
m <- rpois(I,1)+1
m <- m + rbinom(I,1,1/2)*(20-2*m)
n <- 20-m
D.corr <- .99 # controls how u-shaped m,n are
Sigma <- matrix(D.corr,k,k) + diag(k)*(1-D.corr)
D.latent <- rmvnorm(I,mean=rep(0,k),Sigma)
D.idx <- D.latent>0
m <- pmax(pmin(rowSums(D.idx),k-1),1) # enforce m>0,n>0
n <- k-m
Delta <- 0
## data <- rbinormal(I=I,m=m,n=n,mean.x=0,mean.y=Delta,cov.xx=0,cov.xy=0,cov.yy=0,var.x=1,var.y=1,mean.z=(n*(n>m)-m*(n<=m))/5,var.z=var.z,plot=FALSE)
mm <- (n*(n>m)-m*(n<=m))/5
data <- rbinormal(I=I,m=m,n=n,mean.x=mm,mean.y=mm+Delta,cov.xx=0,cov.xy=0,cov.yy=0,var.x=1/2,var.y=1/2,plot=FALSE)
x <- data$x; y <- data$y
auc.cluster(x,y,get.vcov=FALSE)
## save.image('211231b.RData')

pch.control <- '|'
pch.case <- '-'
plot(0,type='n',xlim=range(unlist(c(x,y))),ylim=c(0,I))
for(i in 1:I) points(x[[i]],y=rep(i,length(x[[i]])),pch=pch.control)
for(i in 1:I) points(y[[i]],y=rep(i,length(y[[i]])),pch=pch.case)
points(unlist(c(x,y)),y=rep(0,sum(m)+sum(n)),pch=rep(c(pch.control,pch.case),c(sum(m),sum(n))))


## 2c. plotting routine after 1a and 1b have been run

## op <- par(mfrow=c(1,2))
rm(list=ls())
for(save.file in c('211231a','211231b')) {
    load(paste0(save.file,'.RData'))
    pch.control <- '|'
    pch.case <- '-'
    png(paste0(save.file,'.png'))
    plot(0,type='n',xlim=range(unlist(c(x,y))),ylim=c(0,I),xlab='predictor',ylab='clusters',xaxt='n',yaxt='n',mgp=c(1,1,0))   
    for(i in 1:I) points(x[[i]],y=rep(i,length(x[[i]])),pch=pch.control)
    for(i in 1:I) points(y[[i]],y=rep(i,length(y[[i]])),pch=pch.case)
    abline(h=1/2)
    points(unlist(c(x,y)),y=rep(0,sum(m)+sum(n)),pch=rep(c(pch.control,pch.case),c(sum(m),sum(n))))
    dev.off()
}
## par(op)















## 2. coverage simulation
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
    ## sim <- lapply(1:3e0, FUN=function(jj){
    sim <- mclapply(1:3e1, mc.cores=detectCores()-3,FUN=function(jj){
        D.latent <- rmvnorm(I,mean=rep(0,k),Sigma)
        D.idx <- D.latent>0
        m <- pmax(pmin(rowSums(D.idx),k-1),1) # enforce m>0,n>0
        n <- k-m
        data <- rbinormal(I=I,m=m,n=n,mean.x=0,mean.y=Delta,cov.xy=rho)
        hat <- auc.cluster(data$x,data$y)
        q <- qchisq(1-alpha,df=2)
        cover <-  #tryCatch(
            inside.ellipse(pt=c(theta.11,theta.12),  center=c(hat$theta.11.hat,hat$theta.12.hat), qform=hat$vcov.hat/I, radius=sqrt(q))
        ## cover.11 <- prod(hat$theta.11.hat + c(-1,1)*sqrt(hat$vcov.hat[1,1]/I)*qnorm(1-alpha/2) - theta.11)<0
        ## cover.12 <- prod(hat$theta.12.hat + c(-1,1)*sqrt(hat$vcov.hat[2,2]/I)*qnorm(1-alpha/2) - theta.12)<0
        ## cover.12a <- with(as.list(auc.obu(data$x,data$y)),CI.lower<theta.12 && CI.upper>theta.12)
        cover.11 <- cover.12 <- cover.12a <- NA
        ## if(cover.12!=cover.12a) browser()
           #,error=function(e){browser();print(1)})  ## small I--degenerate AUC
        list(theta.11.hat=hat$theta.11.hat,theta.12.hat=hat$theta.12.hat,cover=cover,cover.11=cover.11,cover.12=cover.12,cover.12a=cover.12a,vcov.hat=hat$vcov.hat)
    })
    theta.11.hats <- sapply(sim,function(ll)ll$theta.11.hat)
    theta.12.hats <- sapply(sim,function(ll)ll$theta.12.hat)
    coverage <- mean(sapply(sim,function(ll)ll$cover))
    coverage.11 <- mean(sapply(sim,function(ll)ll$cover.11))
    coverage.12 <- mean(sapply(sim,function(ll)ll$cover.12))
    coverage.12a <- mean(sapply(sim,function(ll)ll$cover.12a))
    vcov.hats <- simplify2array(lapply(sim,function(ll)ll$vcov.hat))
    list(coverage=coverage,coverage.11=coverage.11,coverage.12=coverage.12,coverage.12a=coverage.12a,bias.theta.11=mean(theta.11.hats)-theta.11, bias.theta.12=mean(theta.12.hats)-theta.12,vcov.hat=apply(vcov.hats,1:2,mean),vcov.mc=matrix(c(var(theta.11.hats),rep(cov(theta.11.hats,theta.12.hats),2),var(theta.12.hats))*I,2))
})


require(xtable)
## out <- sapply(by.param, function(lst)unlist(lst)[c('coverage','bias.theta.11','bias.theta.12','vcov.hat1','vcov.hat2','vcov.hat4','vcov.mc1','vcov.mc2','vcov.mc4')])
## out <- cbind(params,t(out))
## out <- sapply(by.param, function(lst)with(lst,c(coverage=coverage,coverage.11=coverage.11,coverage.12=coverage.12,coverage.12a=coverage.12a,bias.theta.11=bias.theta.11,bias.theta.12=bias.theta.12,structure((vcov.hat-vcov.mc)[c(1,2,4)],names=c('vcov.11','vcov.12','vcov.22')))))
out <- sapply(by.param, function(lst)with(lst,c(coverage=coverage,bias.theta.11=bias.theta.11,bias.theta.12=bias.theta.12,structure((vcov.hat-vcov.mc)[c(1,2,4)],names=c('vcov.11','vcov.12','vcov.22')))))
out <- cbind(params,t(out))
## sink('220105.tex')
print.xtable(xtable(out))
## sink()
## lines <- scan('220105.tex',what='',sep='\n')
## start.idx <- grep('begin\\{tabular\\}',lines)
## end.idx <- grep('end\\{tabular\\}',lines)
## writeLines(lines[start.idx:end.idx],'220105.tex')


require(mvtnorm)
I <- 5
rho <- -1/(I-1) #- .00000003
sigma <- matrix(rho,I,I) + (1-rho)*diag(I)
rmvnorm(1,sigma=sigma)

