## 1. figures showing large indiv auc vs small pop auc and vice versa


## 1a. large indiv auc, small pop auc

require(mvtnorm)
source('../misc.R')
I <- 15
max.cluster.size <- 20
m <- sample(1:max.cluster.size,I,replace=TRUE); n <- sample(1:max.cluster.size,I,replace=TRUE)
Delta <- 2
var.z <- 1e2
data <- rbinormal_old(I=I,m=m,n=n,mean.x=0,mean.y=Delta,cov.xx=0,cov.xy=0,cov.yy=0,var.x=1,var.y=1,mean.z=rep(0,I),var.z=var.z,plot=FALSE)
x <- data$x; y <- data$y
auc.cluster(x,y,get.vcov=FALSE)
## save.image('211223a.RData')
## save.image('220816a.RData')

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
for(save.file in c('220816a','211223b')) {
    load(paste0(save.file,'.RData'))
    pch.control <- '|'
    pch.case <- '-'
    png(paste0(save.file,'.pdf'))
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















## ## 2. coverage simulation



## require(mvtnorm)
## require(parallel)
## source('../misc.R')
## I <- 1e2
## alpha <- .05
## q <- qchisq(1-alpha,df=2)
## k <- 5
## ## D.corr <- .4
## ## intra.corr <- .4
## ## theta.12 <- auc.pop <- .7
## theta.12s <- c(.6,.8)
## theta.11s <- sapply(theta.12s, function(theta.12)seq(theta.12,.95,len=3))
## aucs <- do.call(rbind,lapply(1:length(theta.12s),function(i)cbind(theta.12=theta.12s[i],theta.11=theta.11s[,i])))
## params <- expand.grid(list(D.corr=c(0,.4,.8,.95)))
## params <- as.data.frame(cbind(aucs[rep(1:nrow(aucs),each=nrow(params)),], params[rep(1:nrow(params),nrow(aucs)),,drop=FALSE]))
## by.param <- apply(params,1,function(r) {
##     print(r)
##     D.corr <- unname(r['D.corr']);  theta.12 <- unname(r['theta.12']);  theta.11 <- unname(r['theta.11'])
##     Delta <- sqrt(2)*qnorm(theta.12)
##     rho <- 1 - 1/2*Delta^2 / qnorm(theta.11)^2
##     rho <- round(rho,9)
##     Sigma <- matrix(D.corr,k,k) + diag(k)*(1-D.corr)
##     ## sim <- lapply(1:3e0, FUN=function(jj){
##     sim <- mclapply(1:1e2, mc.cores=detectCores()-3,FUN=function(jj){
##         D.latent <- rmvnorm(I,mean=rep(0,k),Sigma)
##         D.idx <- D.latent>0
##         m <- pmax(pmin(rowSums(D.idx),k-1),1) # enforce m>0,n>0
##         n <- k-m
##         data <- rbinormal(I=I,m=m,n=n,mean.x=0,mean.y=Delta,cov.xy=rho)
##         hat <- auc.cluster(data$x,data$y)
##         cover <-  #tryCatch(
##             inside.ellipse(pt=c(theta.11,theta.12),  center=c(hat$theta.11.hat,hat$theta.12.hat), qform=hat$vcov.hat/I, radius=sqrt(q))
##         ## cover.11 <- prod(hat$theta.11.hat + c(-1,1)*sqrt(hat$vcov.hat[1,1]/I)*qnorm(1-alpha/2) - theta.11)<0
##         ## cover.12 <- prod(hat$theta.12.hat + c(-1,1)*sqrt(hat$vcov.hat[2,2]/I)*qnorm(1-alpha/2) - theta.12)<0
##         ## cover.12a <- with(as.list(auc.obu(data$x,data$y)),CI.lower<theta.12 && CI.upper>theta.12)
##         cover.11 <- cover.12 <- cover.12a <- NA
##         ## if(cover.12!=cover.12a) browser()
##            #,error=function(e){browser();print(1)})  ## small I--degenerate AUC
##         list(theta.11.hat=hat$theta.11.hat,theta.12.hat=hat$theta.12.hat,cover=cover,cover.11=cover.11,cover.12=cover.12,cover.12a=cover.12a,vcov.hat=hat$vcov.hat)
##     })
##     theta.11.hats <- sapply(sim,function(ll)ll$theta.11.hat)
##     theta.12.hats <- sapply(sim,function(ll)ll$theta.12.hat)
##     coverage <- mean(sapply(sim,function(ll)ll$cover))
##     coverage.11 <- mean(sapply(sim,function(ll)ll$cover.11))
##     coverage.12 <- mean(sapply(sim,function(ll)ll$cover.12))
##     coverage.12a <- mean(sapply(sim,function(ll)ll$cover.12a))
##     vcov.hats <- simplify2array(lapply(sim,function(ll)ll$vcov.hat))
##     list(coverage=coverage,coverage.11=coverage.11,coverage.12=coverage.12,coverage.12a=coverage.12a,bias.theta.11=mean(theta.11.hats)-theta.11, bias.theta.12=mean(theta.12.hats)-theta.12,vcov.hat=apply(vcov.hats,1:2,mean),vcov.mc=matrix(c(var(theta.11.hats),rep(cov(theta.11.hats,theta.12.hats),2),var(theta.12.hats))*I,2))
## })


## require(xtable)
## ## out <- sapply(by.param, function(lst)unlist(lst)[c('coverage','bias.theta.11','bias.theta.12','vcov.hat1','vcov.hat2','vcov.hat4','vcov.mc1','vcov.mc2','vcov.mc4')])
## ## out <- cbind(params,t(out))
## ## out <- sapply(by.param, function(lst)with(lst,c(coverage=coverage,coverage.11=coverage.11,coverage.12=coverage.12,coverage.12a=coverage.12a,bias.theta.11=bias.theta.11,bias.theta.12=bias.theta.12,structure((vcov.hat-vcov.mc)[c(1,2,4)],names=c('vcov.11','vcov.12','vcov.22')))))
## out <- sapply(by.param, function(lst)with(lst,c(coverage=coverage,bias.theta.11=bias.theta.11,bias.theta.12=bias.theta.12,structure((vcov.hat-vcov.mc)[c(1,2,4)],names=c('vcov.11','vcov.12','vcov.22')))))
## out <- cbind(params,t(out))
## rownames(out) <- NULL
## ## sink('220105.tex')
## colnames(out) <- c('$\\theta_{12}$','$\\theta_{11}$','$\\rho_{MN}$','coverage','$\\theta_{11}$','$\\theta_{12}$','$\\Sigma_{11}$','$\\Sigma_{12}$','$\\Sigma_{22}$')
## addtorow <- list()
## addtorow$pos <- list(0,0,0,x0)
## addtorow$command <- #c("& \\multicolumn{8}{c}{Grade 6} \\\\\n",
##     c('\\multicolumn{3}{|c||}{parameters} & coverage &\\multicolumn{5}{c|}{bias} \\\\\n',
##       '\\hline\n',
##       '\\hline\n',
##                       '$\\theta_{12}$ & $\\theta_{11}$ & $\\rho_{MN}$ &  & $\\theta_{11}$ & $\\theta_{12}$ & $\\Sigma_{11}$ & $\\Sigma_{12}$ & $\\Sigma_{22}$ \\\\\n')
## ## "Grade 3 & A & B & C & D \\\\\n")
## col.wid <- 'p{1.4cm}'
## c('c',rep(c('|',rep(col.wid,3),'|'),3))
## align.str <- paste(c('c',rep('p{1.4cm}|',9)),collapse='')
## align.str <-  "|c|C{1.4cm}|C{1.4cm}|C{1.4cm}||C{1.4cm}||C{1.4cm}|C{1.4cm}|C{1.4cm}|C{1.4cm}|C{1.4cm}|"
## print.xtable(xtable(out, align=align.str),include.rownames=FALSE,include.colnames=FALSE,sanitize.text.function=function(x)x, add.to.row=addtorow)
## ## sink()
## lines <- scan('220105.tex',what='',sep='\n')
## start.idx <- grep('begin\\{tabular\\}',lines)
## end.idx <- grep('end\\{tabular\\}',lines)
## writeLines(lines[start.idx:end.idx],'220105.tex')



## require(mvtnorm)
## I <- 5
## rho <- -1/(I-1) #- .00000003
## sigma <- matrix(rho,I,I) + (1-rho)*diag(I)
## rmvnorm(1,sigma=sigma)



## 2. coverage simulation
## (8/10/22)
## slightly re-factored so the same code also works with both #1 binormal and #2 censored binormal models

require(mvtnorm)
require(parallel)
source('../misc.R')

bound <- abs(qnorm(.2)) # 2
auc.to.params <- function(theta.12,theta.11)auc.to.params.binormal.censored(theta.12,theta.11,bound) #2
sampler <- rbinormal.censored #2
## auc.to.params <- auc.to.params.binormal # binormal #1
## sampler <- rbinormal #1

B <- 1e3
I <- 60
alpha <- .05
q <- qchisq(1-alpha,df=2)
ks <- 2:5
## D.corr <- .4
## intra.corr <- .4
## theta.12 <- auc.pop <- .7
theta.12s <- c(.7,.8)
## theta.11s <-sapply(theta.12s, function(theta.12)seq(theta.12,.95,len=3))
theta.11s <-sapply(theta.12s, function(theta.12)c(theta.12,theta.12+.1))
aucs <- do.call(rbind,lapply(1:length(theta.12s),function(i)cbind(theta.12=theta.12s[i],theta.11=theta.11s[,i])))
## aucs <- expand.grid(theta.12=theta.12s,theta.11=theta.11s)
params <- expand.grid(list(D.corr=c(0,.1,.4,.8)))
## params <- expand.grid(list(D.corr=c(.4,.8)))#!!!!!!
params <- as.data.frame(cbind(aucs[rep(1:nrow(aucs),each=nrow(params)),], params[rep(1:nrow(params),nrow(aucs)),,drop=FALSE]))
params <- round(params,1)
rownames(params) <- NULL
by.param <- apply(params,1,function(r) {
    ## browser()
    print(r)
    k <- sample(ks,1)
    D.corr <- unname(r['D.corr']);  theta.12 <- unname(r['theta.12']);  theta.11 <- unname(r['theta.11'])
    Delta.rho <- auc.to.params(theta.12,theta.11)
    rho <- round(Delta.rho['rho'],9)
    Delta <- Delta.rho['Delta']
    ## browser()
    Sigma <- matrix(D.corr,k,k) + diag(k)*(1-D.corr)
    sim <- lapply(1:B, FUN=function(jj){
    ## sim <- mclapply(1:B, mc.cores=detectCores()-3,FUN=function(jj){
        D.latent <- rmvnorm(I,mean=rep(0,k),Sigma)
        D.idx <- D.latent>0
        m <- pmax(pmin(rowSums(D.idx),k-1),1) # enforce m>0,n>0
        n <- k-m
        data <- sampler(I=I,m=m,n=n,mean.x=0,mean.y=Delta,cov.xy=rho)
        hat <- auc.cluster(data$x,data$y)
        cover <-  #tryCatch(
            inside.ellipse(pt=c(theta.11,theta.12),  center=c(hat$theta.11.hat,hat$theta.12.hat), qform=hat$vcov.hat/I, radius=sqrt(q))#,error=function(e){cat('!!!');browser()})
        ## cover.11 <- prod(hat$theta.11.hat + c(-1,1)*sqrt(hat$vcov.hat[1,1]/I)*qnorm(1-alpha/2) - theta.11)<0
        ## cover.12 <- prod(hat$theta.12.hat + c(-1,1)*sqrt(hat$vcov.hat[2,2]/I)*qnorm(1-alpha/2) - theta.12)<0
        ## cover.12a <- with(as.list(auc.obu(data$x,data$y)),CI.lower<theta.12 && CI.upper>theta.12)
        cover.11 <- cover.12 <- cover.12a <- NA
        ## if(cover.12!=cover.12a) browser()
           #,error=function(e){browser();print(1)})  ## small I--degenerate AUC
        list(theta.11.hat=hat$theta.11.hat,theta.12.hat=hat$theta.12.hat,cover=cover,cover.11=cover.11,cover.12=cover.12,cover.12a=cover.12a,vcov.hat=hat$vcov.hat)
    })
    ## browser()
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
out <- sapply(by.param, function(lst)with(lst,c(coverage=coverage,bias.theta.12=bias.theta.12,bias.theta.11=bias.theta.11,structure((vcov.hat-vcov.mc)[c(1,2,4)],names=c('vcov.11','vcov.12','vcov.22')))))
out <- cbind(params,t(out))
rownames(out) <- NULL
colnames(out) <- c('$\\theta_{12}$','$\\theta_{11}$','$\\rho_{MN}$','coverage','$\\theta_{12}$','$\\theta_{11}$','$\\Sigma_{11}$','$\\Sigma_{12}$','$\\Sigma_{22}$')
addtorow <- list()
addtorow$pos <- list(0,0,0,0)
addtorow$command <- #c("& \\multicolumn{8}{c}{Grade 6} \\\\\n",
    c('\\multicolumn{3}{|c||}{parameters} & coverage &\\multicolumn{5}{c|}{bias} \\\\\n',
      '\\hline\n',
      '\\hline\n',
                      '$\\theta_{12}$ & $\\theta_{11}$ & $\\rho_{MN}$ &  & $\\theta_{12}$ & $\\theta_{11}$  & $\\Sigma_{11}$ & $\\Sigma_{12}$ & $\\Sigma_{22}$ \\\\\n')
## "Grade 3 & A & B & C & D \\\\\n")
col.wid <- 'p{1.4cm}'
c('c',rep(c('|',rep(col.wid,3),'|'),3))
align.str <- paste(c('c',rep('p{1.4cm}|',9)),collapse='')
align.str <-  "|c|C{1.4cm}|C{1.4cm}|C{1.4cm}||C{1.4cm}||C{1.4cm}|C{1.4cm}|C{1.4cm}|C{1.4cm}|C{1.4cm}|"
## save.file <- 'figs/220105.tex' #1
save.file <- 'figs/220811.tex' #2
## sink(save.file)
print.xtable(xtable(out, align=align.str),include.rownames=FALSE,include.colnames=FALSE,sanitize.text.function=function(x)x, add.to.row=addtorow)
## sink() 
lines <- scan(save.file,what='',sep='\n')
start.idx <- grep('begin\\{tabular\\}',lines)
end.idx <- grep('end\\{tabular\\}',lines)
writeLines(lines[start.idx:end.idx],save.file)

## write out stats for \input
cat(I,'%',file='input/sim_coverage_I.txt')
cat(format(B,big.mark=','),'%',file='input/sim_coverage_reps.txt')
cat(paste(ks,collapse=', '),'%',file='input/sim_coverage_ks.txt')
for(col in colnames(params)) 
    cat(paste(unique(params[,col]),collapse=', '),'%',file=paste0('input/sim_coverage_',col,'s.txt'))
    
dd




## ## 3. data analysis
## source('../misc.R')
## data.dir <- '../data/'
## statuses.nyc <- c(quote(SUSPECT_RACE_DESCRIPTION=='BLACK' | SUSPECT_RACE_DESCRIPTION=='BLACK HISPANIC'),
##               quote(SUSPECT_RACE_DESCRIPTION=='BLACK'),
##               quote(SUSPECT_RACE_DESCRIPTION=='WHITE' | SUSPECT_RACE_DESCRIPTION=='WHITE HISPANIC'),
##               quote(SUSPECT_RACE_DESCRIPTION=='WHITE'),
##               quote(SUSPECT_RACE_DESCRIPTION=='WHITE HISPANIC' | SUSPECT_RACE_DESCRIPTION=='BLACK HISPANIC'))
## statuses <- c(quote(race=='black.nonhisp' | race=='black.hisp'),
##               quote(race=='black.nonhisp'),
##               quote(race=='white.nonhisp' | race=='white.hisp'),
##               quote(race=='white.nonhisp'),
##               quote(race=='white.hisp' | race=='black.hisp'))
## statuses.boston <- c(quote(race=='Black'),
##               quote(race=='Black' & ethnicity!='Hispanic Origin'),
##               quote(race=='White'),
##               quote(race=='White' & ethnicity!='Hispanic Origin'),
##               quote(ethnicity=='Hispanic Origin'))
## status.descriptions <- c('Black','Black, non-Hispanic','White','White, non-Hispanic','Hispanic')



## ## nyc data
## filelist <- dir(data.dir)
## filelist <- filelist[grep('^sqf[-0-9]+\\.csv',filelist)]
## filelist <- paste0(data.dir,filelist)
## ## snf <- lapply(filelist[3:4], read.csv)
## snf <- lapply(filelist[], function(file) {
##     with(read.csv(file),
##          data.frame(STOP_DURATION_MINUTES,SUSPECT_RACE_DESCRIPTION,ISSUING_OFFICER_COMMAND_CODE))
## })
## snf <- do.call(rbind,snf)

## ## clean-up
## snf <- within(snf, {
##     STOP_DURATION_MINUTES[STOP_DURATION_MINUTES > 300] <- NA
##     SUSPECT_RACE_DESCRIPTION[SUSPECT_RACE_DESCRIPTION=='MALE'] <- NA
##     SUSPECT_RACE_DESCRIPTION[SUSPECT_RACE_DESCRIPTION=='(null)'] <- NA
##     SUSPECT_RACE_DESCRIPTION[SUSPECT_RACE_DESCRIPTION %in% c('AMER IND','AMERICAN INDIAN/ALASKAN N','AMERICAN INDIAN/ALASKAN NATIVE','MIDDLE EASTERN/SOUTHWEST','MIDDLE EASTERN/SOUTHWEST ASIAN')] <- 'other'#'OTHER'
##     SUSPECT_RACE_DESCRIPTION[SUSPECT_RACE_DESCRIPTION %in% c('ASIAN/PAC.ISL','ASIAN / PACIFIC ISLANDER')] <- 'asian'
##     SUSPECT_RACE_DESCRIPTION[SUSPECT_RACE_DESCRIPTION=='BLACK'] <- 'black.nonhisp'
##     SUSPECT_RACE_DESCRIPTION[SUSPECT_RACE_DESCRIPTION=='BLACK HISPANIC'] <- 'black.hisp'
##     SUSPECT_RACE_DESCRIPTION[SUSPECT_RACE_DESCRIPTION=='WHITE'] <- 'white.nonhisp'
##     SUSPECT_RACE_DESCRIPTION[SUSPECT_RACE_DESCRIPTION=='WHITE HISPANIC'] <- 'white.hisp'
##     #'ASIAN / PACIFIC ISLANDER'
##     ## SUSPECT_RACE_DESCRIPTION[SUSPECT_RACE_DESCRIPTION=='ASIAN / PACIFIC ISLANDER'] <- 'ASIAN')
##     ## SUSPECT_RACE_DESCRIPTION[SUSPECT_RACE_DESCRIPTION=='MIDDLE EASTERN/SOUTHWEST'] <- 'MIDDLE EASTERN/SOUTHWEST ASIAN'
## })
## ## aggregate table
## mean.durations <- aggregate(STOP_DURATION_MINUTES ~ SUSPECT_RACE_DESCRIPTION, mean, data=snf)
## sd.durations <- aggregate(STOP_DURATION_MINUTES ~ SUSPECT_RACE_DESCRIPTION, sd, data=snf)
## colnames(mean.durations)  <- c('group','mean.duration')
## colnames(sd.durations) <- c('group','sd.duration')
## counts <- table(snf$SUSPECT_RACE_DESCRIPTION)
## counts <- data.frame(group=names(counts),count=unname(as.integer(counts)))
## summary.nyc <- merge(merge(mean.durations,sd.durations),counts)
## summary.nyc$freq <- summary.nyc$count / sum(summary.nyc$count)
## cat(format(nrow(snf),big.mark=','),'%',file='input/total_stops_nyc.txt')


## ## AUC analysis
## nyc <- sapply(statuses.nyc, function(status) {
##     browser()
##     snf$race <- with(snf, eval(status))
##     ## snf$race <- snf$SUSPECT_RACE_DESCRIPTION=='WHITE' | snf$SUSPECT_RACE_DESCRIPTION=='WHITE HISPANIC'#BLACK HISPANIC'##BLACK, WHITE HISPANIC
##     ## ## snf$race <- snf$SUSPECT_RACE_DESCRIPTION=='BLACK HISPANIC' | snf$SUSPECT_RACE_DESCRIPTION=='BLACK'
##     ## snf$race <- snf$SUSPECT_RACE_DESCRIPTION=='BLACK' | snf$SUSPECT_RACE_DESCRIPTION=='BLACK HISPANIC'
##     clusters <- split(snf,snf$ISSUING_OFFICER_COMMAND_CODE)
##     ## biomarker <- 'SUSPECT_WEIGHT'
##     ## biomarker <- 'SUSPECT_REPORTED_AGE_FLAG'
##     biomarker <- 'STOP_DURATION_MINUTES'
##     x <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=race==FALSE,select=biomarker)))))
##     y <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=race==TRUE,select=biomarker)))))
##     omit.idx <- sapply(x,length)==0 | sapply(y,length)==0
##     x <- x[!omit.idx]; y <- y[!omit.idx]
##     auc0 <- auc.cluster(x=x,y=y) 
##     contrast <- matrix(c(1,-1),nrow=2)
##     alpha <- .05
##     I <- length(x)
##     z.stat <- with(auc0, sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
##     ## (1-pnorm(abs(z.stat)))*2
##     ## abs(z.stat) > qnorm(1-alpha/2)
##     if(paste(status,collapse='')==paste(statuses.nyc[[1]],collapse='')) {
##         ## pdf('220823a.pdf')
##         plot(auc0,resol=1e2)
##         plot(auc0,alpha=.01,add=TRUE)
##         ## dev.off()
##     }
##     ## c(auc0,list(equality.pval=(1-pnorm(abs(z.stat)))*2))    
##     with(auc0,structure(c(theta.12.hat,theta.12.CI,theta.11.hat,theta.11.CI,(1-pnorm(abs(z.stat)))*2,I,length(unlist(x)),length(unlist(y))),
##                         names=c('theta.12.hat','theta.12.CI.lower','theta.12.CI.upper','theta.11.hat','theta.11.CI.lower','theta.11.CI.upper','equality.pval','I','M','N')))
## })
## nyc <- t(nyc)

## ## boston data
## years <- 2019:2021 # stop_duration coded differently for year<=2018
## snfs <- lapply(years, function(year) {
##     ## if(year==2017)browser()
##     snf <- read.csv(paste0(data.dir,'fio_',year,'.csv'))
##     snf.names <- read.csv(paste0(data.dir,'fio_names_',year,'.csv'))
##     snf <- merge(snf,snf.names,by='fc_num')
##     ## browser()
##     ## print(year)
##     ## print(table(snf$stop_duration))
##     ## snf$biomarker <- snf$stop_duration
##     ## snf$biomarker[snf$biomarker=='NULL'] <- NA
##     ## snf$biomarker <- as.numeric(snf$biomarker)
##     ## ## snf$status <- snf$race=='White'  & snf$ethnicity=='Not of Hispanic Origin'#& snf$ethnicity!='Hispanic Origin'
##     ## ## snf$status <- snf$race=='White'  & snf$ethnicity!='Hispanic Origin'
##     ## snf$status <- with(snf, eval(status))
##     subset(snf,select=c(fc_num,race,ethnicity,stop_duration,contact_officer))
## })
## snf <- do.call(rbind,snfs)
## ## clean up
## snf$biomarker <- snf$stop_duration
## snf$biomarker[snf$biomarker=='NULL'] <- NA
## snf$biomarker <- as.numeric(snf$biomarker)
## snf$biomarker[snf$biomarker>300] <- NA
## snf$race2 <- snf$race
## snf <- within(snf, {
##     race2[race2==''] <- NA
##     race2[race2 %in% c('Unknown','NULL')] <- NA
##     race2[race2 %in% c("American Indian or Alaskan Native","Native Hawaiian or Other Pacific Islander","Other")] <- 'other'
##     race2[race=='Black' & ethnicity=='Hispanic Origin'] <- 'Black Hispanic'
##     race2[race=='Black' & ethnicity=='Not of Hispanic Origin'] <- 'Black non-Hispanic'
##     race2[race=='White' & ethnicity=='Hispanic Origin'] <- 'White Hispanic'
##     race2[race=='White' & ethnicity=='Not of Hispanic Origin'] <- 'White non-Hispanic'
##     race2[!(ethnicity %in% c('Hispanic Origin','Not of Hispanic Origin'))] <- NA
## })
## ## aggregate data
## mean.durations <- aggregate(biomarker ~ race2, mean, data=snf)
## sd.durations <- aggregate(biomarker ~ race2, sd, data=snf)
## colnames(mean.durations)  <- c('group','mean.duration')
## colnames(sd.durations) <- c('group','sd.duration')
## counts <- table(snf$race2)
## counts <- data.frame(group=names(counts),count=unname(as.integer(counts)))
## summary.boston <- merge(merge(mean.durations,sd.durations),counts)
## summary.boston$freq <- summary.boston$count / sum(summary.boston$count)
## ## cat(format(nrow(snf),big.mark=','),'%',file='input/total_stops_boston.txt')
## ## AUC analysis
## boston <- sapply(statuses.boston, function(status) {
##     ## snfs <- lapply(years, function(year) {
##     ##     ## if(year==2017)browser()
##     ##     snf <- read.csv(paste0(data.dir,'fio_',year,'.csv'))
##     ##     snf.names <- read.csv(paste0(data.dir,'fio_names_',year,'.csv'))
##     ##     snf <- merge(snf,snf.names,by='fc_num')
##     ##     ## browser()
##     ##     ## print(year)
##     ##     ## print(table(snf$stop_duration))
##     ##     snf$biomarker <- snf$stop_duration
##     ##     snf$biomarker[snf$biomarker=='NULL'] <- NA
##     ##     snf$biomarker <- as.numeric(snf$biomarker)
##     ##     ## snf$status <- snf$race=='White'  & snf$ethnicity=='Not of Hispanic Origin'#& snf$ethnicity!='Hispanic Origin'
##     ##     ## snf$status <- snf$race=='White'  & snf$ethnicity!='Hispanic Origin'
##     ##     snf$status <- with(snf, eval(status))
##     ##     subset(snf,select=c(fc_num,status,biomarker,race,contact_officer))
##     ## })
##     ## snf <- do.call(rbind,snfs)
##     snf$status <- with(snf, eval(status))
##     snf <- snf[order(snf$fc_num),]
##     single.encounters <- with(rle(snf$fc_num), values[lengths==1])
##     snf <- snf[snf$fc_num %in% single.encounters,]
##     clusters <- split(snf,snf$contact_officer)
##     x <- lapply(clusters, function(cluster)na.omit(cluster$biomarker[cluster$status==FALSE]))
##     y <- lapply(clusters, function(cluster)na.omit(cluster$biomarker[cluster$status==TRUE]))
##     omit.idx <- sapply(x,length)==0 | sapply(y,length)==0
##     x <- x[!omit.idx]; y <- y[!omit.idx]
##     (( auc0 <- auc.cluster(x=x,y=y) ))
##     contrast <- matrix(c(1,-1),nrow=2)
##     alpha <- .05
##     stopifnot(length(x)==length(y))
##     I <- length(x)
##     z.stat <- with(auc0, sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
##     ## (1-pnorm(abs(z.stat)))*2
##     ## abs(z.stat) > qnorm(1-alpha/2)
##     ## plot(auc0,resol=1e2)
##     ## c(auc0,list(equality.pval=(1-pnorm(abs(z.stat)))*2))
##     if(paste(status,collapse='')==paste(statuses.boston[[1]],collapse='')) {
##         ## pdf('220823b.pdf')
##         plot(auc0,resol=1e2)
##         plot(auc0,alpha=.01,add=TRUE)
##         ## dev.off()
##     }
##     ## c(auc0,list(equality.pval=(1-pnorm(abs(z.stat)))*2))    
##     with(auc0,structure(c(theta.12.hat,theta.12.CI,theta.11.hat,theta.11.CI,(1-pnorm(abs(z.stat)))*2,I,length(unlist(x)),length(unlist(y))),
##                         names=c('theta.12.hat','theta.12.CI.lower','theta.12.CI.upper','theta.11.hat','theta.11.CI.lower','theta.11.CI.upper','equality.pval','I','M','N')))
##     })
## boston <- t(boston)
## cbind(status.descriptions,as.data.frame(boston))

## ## 3a table of AUC estimate--run after above to get nyc and boston data
## require(xtable)
## ## out <- as.data.frame(round(rbind(t(nyc),t(boston))[c(1,5,2,6,3,7,4,8),],2))
## idx <- as.numeric(t(matrix(1:(2*nrow(nyc)),ncol=2)))
## out <- as.data.frame(round(rbind(nyc,boston)[idx,],2))
## ## out <- cbind(c('nyc','boston'),out)
## ## out <- cbind(case.status=rep(status.descriptions,each=2),out)
## out$city <- c('NYC','Boston')
## out$case.status <- rep(status.descriptions,each=2)
## cat(out$theta.12.hat[out$case.status=='Black' & out$city=='NYC'],'%',file='input/da_black_nyc_theta12_mean.txt')
## cat(out$theta.12.CI.lower[out$case.status=='Black' & out$city=='NYC'],'---',out$theta.12.CI.upper[out$case.status=='Black' & out$city=='NYC'],'%',sep='',file='input/da_black_nyc_theta12_CI.txt')
## cat(out$theta.11.hat[out$case.status=='Black' & out$city=='NYC'],'%',file='input/da_black_nyc_theta11_mean.txt')
## cat(out$theta.11.CI.lower[out$case.status=='Black' & out$city=='NYC'],'---',out$theta.11.CI.upper[out$case.status=='Black' & out$city=='NYC'],'%',sep='',file='input/da_black_nyc_theta11_CI.txt')
## out <- within(out,theta.12 <- paste0(theta.12.hat,' [',theta.12.CI.lower,', ',theta.12.CI.upper,']'))
## out <- within(out,theta.11 <- paste0(theta.11.hat,' [',theta.11.CI.lower,', ',theta.11.CI.upper,']'))
## out <- subset(out,select=c(case.status,city,I,M,N,theta.12,theta.11,equality.pval))
## cat(out$theta.12[out$case.status=='Black' & out$city=='Boston'],'%',file='input/da_black_boston_theta12.txt')
## cat(out$theta.11[out$case.status=='Black' & out$city=='Boston'],'%',file='input/da_black_boston_theta11.txt')
## cat(out$theta.12[out$case.status=='White' & out$city=='NYC'],'%',file='input/da_white_nyc_theta12.txt')
## cat(out$theta.11[out$case.status=='White' & out$city=='NYC'],'%',file='input/da_white_nyc_theta11.txt')
## cat(out$theta.12[out$case.status=='White' & out$city=='Boston'],'%',file='input/da_white_boston_theta12.txt')
## cat(out$theta.11[out$case.status=='White' & out$city=='Boston'],'%',file='input/da_white_boston_theta11.txt')
## cat(out$equality.pval[out$case.status=='White' & out$city=='Boston'],'%',file='input/da_white_boston_equality_pval.txt')
## out <- within(out,case.status[duplicated(case.status)] <- '')
## colnames(out) <- c('case status','data set','I','$\\Sigma M_i$','$\\Sigma N_i$','$\\theta_{12}$','$\\theta_{11}$','$H_{0}:\\theta_{12}=\\theta_{11}$')
## ## addtorow <- list()
## ## addtorow$pos <- list(0,0,0,x0)
## ## addtorow$command <- #c("& \\multicolumn{8}{c}{Grade 6} \\\\\n",
## ##     c('\\multicolumn{3}{|c||}{parameters} & coverage &\\multicolumn{5}{c|}{bias} \\\\\n',
## ##       '\\hline\n',
## ##       '\\hline\n',
## ##                       '$\\theta_{12}$ & $\\theta_{11}$ & $\\rho_{MN}$ &  & $\\theta_{11}$ & $\\theta_{12}$ & $\\Sigma_{11}$ & $\\Sigma_{12}$ & $\\Sigma_{22}$ \\\\\n')
## ## "Grade 3 & A & B & C & D \\\\\n")
## ## col.wid <- 'p{1.4cm}'
## ## c('c',rep(c('|',rep(col.wid,3),'|'),3))
## ## align.str <- paste(c('c',rep('p{1.4cm}|',9)),collapse='')
## ## align.str <-  "|c|C{1.4cm}|C{1.4cm}|C{1.4cm}||C{1.4cm}||C{1.4cm}|C{1.4cm}|C{1.4cm}|C{1.4cm}|C{1.4cm}|"
## sink('220821.tex')
## print.xtable(xtable(out,digits=rep(0,ncol(out)+1)),include.rownames=FALSE,sanitize.text.function=function(x)x)
## sink()
## lines <- scan('220821.tex',what='',sep='\n')
## start.idx <- grep('begin\\{tabular\\}',lines)
## end.idx <- grep('end\\{tabular\\}',lines)
## writeLines(lines[start.idx:end.idx],'220821.tex')



## ## sink('test.txt')
## ## cat('12333')
## ## sink()


## ## 3b table of summary data--run after above to get nyc and boston data

## ## standardize race labels. should be done when reading in data. will
## ## require updating statuses.nyc/.boston. then only one set of
## ## statuses will be needed.

## summary.nyc$group <- c('Asian','Black non-Hispanic','Black Hispanic','other','White non-Hispanic','White Hispanic')
## summary.nyc <- summary.nyc[order(summary.nyc$group),]
## summary.boston <- summary.boston[order(summary.boston$group),]

## out <- cbind(nyc=summary.nyc,boston=subset(summary.boston,select=-group))
## out[,-1] <- round(out[-1],2)
## out$nyc.mean.duration <- paste0(out$nyc.mean.duration,' (',out$nyc.sd.duration,')')
## out$boston.mean.duration <- paste0(out$boston.mean.duration,' (',out$boston.sd.duration,')')
## out <- subset(out,select=-c(nyc.sd.duration,boston.sd.duration))
## colnames(out) <- c('group','mean duration (SD)','count','freq.','mean duration (SD)','count','freq.')
## filename <- '220823.tex'
## sink(filename)
## addtorow <- list()
## addtorow$pos <- list(0,0,0,0)
## addtorow$command <- c("& \\multicolumn{3}{c||}{NYC} & \\multicolumn{3}{c|}{Boston}\\\\\n",
##                       '\\hline\n',
##                       '\\hline\n',
##                       paste0(paste0(colnames(out),collapse=' & '), '\\\\\n'))
## align.str <-  "|c|C{3cm}||C{2cm}|C{2cm}|C{2cm}||C{2cm}|C{2cm}|C{2cm}|"#|C{2cm}|C{2cm}|"
## print.xtable(xtable(out,align=align.str), add.to.row = addtorow,include.colnames=FALSE,include.rownames=FALSE,sanitize.text.function=function(x)x)
## sink()
## lines <- scan(filename,what='',sep='\n')
## start.idx <- grep('begin\\{tabular\\}',lines)
## end.idx <- grep('end\\{tabular\\}',lines)
## writeLines(lines[start.idx:end.idx],filename)



























## 3. data analysis section

## 3a. read in and clean up nyc and boston data
data.dir <- '../data/'
## nyc data
filelist <- dir(data.dir)
filelist <- filelist[grep('^sqf[-0-9]+\\.csv',filelist)]
filelist <- paste0(data.dir,filelist)
## snf <- lapply(filelist[3:4], read.csv)
snf <- lapply(filelist[], function(file) {
    with(read.csv(file),
         data.frame(duration=STOP_DURATION_MINUTES,race=SUSPECT_RACE_DESCRIPTION,cluster=ISSUING_OFFICER_COMMAND_CODE))
})
snf <- do.call(rbind,snf)
snf <- within(snf, {
    duration[duration > 300] <- NA
    race[race=='MALE'] <- NA
    race[race=='(null)'] <- NA
    race[race %in% c('AMER IND','AMERICAN INDIAN/ALASKAN N','AMERICAN INDIAN/ALASKAN NATIVE','MIDDLE EASTERN/SOUTHWEST','MIDDLE EASTERN/SOUTHWEST ASIAN')] <- 'other'#'OTHER'
    race[race %in% c('ASIAN/PAC.ISL','ASIAN / PACIFIC ISLANDER')] <- 'asian'
    race[race=='BLACK'] <- 'black.nonhisp'
    race[race=='BLACK HISPANIC'] <- 'black.hisp'
    race[race=='WHITE'] <- 'white.nonhisp'
    race[race=='WHITE HISPANIC'] <- 'white.hisp'
})
nyc <- snf
## boston data
years <- 2019:2021 # stop_duration coded differently for year<=2018
snfs <- lapply(years, function(year) {
    ## if(year==2017)browser()
    snf <- read.csv(paste0(data.dir,'fio_',year,'.csv'))
    snf.names <- read.csv(paste0(data.dir,'fio_names_',year,'.csv'))
    snf <- merge(snf,snf.names,by='fc_num')
    snf <- snf[order(snf$fc_num),]
    single.encounters <- with(rle(snf$fc_num), values[lengths==1])
    snf <- snf[snf$fc_num %in% single.encounters,]
    ## subset(snf,select=c(fc_num,race,ethnicity,stop_duration,contact_officer))
    with(snf, data.frame(duration=stop_duration,race.old=race,ethnicity=ethnicity,cluster=contact_officer))
})
snf <- do.call(rbind,snfs)
snf$duration[snf$duration=='NULL'] <- NA
snf$duration <- as.numeric(snf$duration)
snf$duration[snf$duration>300] <- NA
snf$race <- snf$race.old
snf <- within(snf, {
    race[race==''] <- NA
    race[race %in% c('Unknown','NULL')] <- NA
    race[race %in% c("American Indian or Alaskan Native","Native Hawaiian or Other Pacific Islander","Other")] <- 'other'
    race[race=='Asian'] <- 'asian'
    race[race=='Black' & ethnicity=='Hispanic Origin'] <- 'black.hisp'
    ## race[race=='Black' & ethnicity=='Not of Hispanic Origin'] <- 'black.nonhisp'
    race[race=='Black' & ethnicity!='Hispanic Origin'] <- 'black.nonhisp'
    race[race=='White' & ethnicity=='Hispanic Origin'] <- 'white.hisp'
    ## race[race=='White' & ethnicity=='Not of Hispanic Origin'] <- 'white.nonhisp'
    race[race=='White' & ethnicity!='Hispanic Origin'] <- 'white.nonhisp'
    ## race[!(ethnicity %in% c('Hispanic Origin','Not of Hispanic Origin'))] <- NA
})
boston <- snf



## 3b aggregate tables
summary.tables <- lapply(list(nyc=nyc,boston=boston), function(snf) {
    mean.durations <- aggregate(duration ~ race, mean, data=snf)
    sd.durations <- aggregate(duration ~ race, sd, data=snf)
    colnames(mean.durations)  <- c('group','mean.duration')
    colnames(sd.durations) <- c('group','sd.duration')
    counts <- table(snf$race)
    counts <- data.frame(group=names(counts),count=unname(as.integer(counts)))
    out <- merge(merge(mean.durations,sd.durations),counts)
    out$freq <- out$count / sum(out$count)
    out
})




## 3c AUC analysis
source('../misc.R')
statuses <- c(quote(race=='black.nonhisp' | race=='black.hisp'),
              quote(race=='black.nonhisp'),
              quote(race=='black.hisp'),
              quote(race=='white.nonhisp' | race=='white.hisp'),
              quote(race=='white.nonhisp'),
              quote(race=='white.hisp'),
              quote(race=='white.hisp' | race=='black.hisp'))
status.descriptions <- c('Black','Black non-Hispanic','Black Hispanic','White','White non-Hispanic','White Hispanic','Hispanic')
## auc.estimates <- lapply(list(nyc=nyc,boston=boston), function(snf) {
auc.estimates <- Map(function(snf,city.name) {
    ## browser()
    sapply(statuses, function(status) {
        ## browser()
        snf$status <- with(snf, eval(status))
        clusters <- split(snf,snf$cluster)
        ## biomarker <- 'STOP_DURATION_MINUTES'
        x <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=status==FALSE,select=duration)))))
        y <- lapply(clusters, function(cluster)na.omit(as.numeric(unlist(subset(cluster,subset=status==TRUE,select=duration)))))
        omit.idx <- sapply(x,length)==0 | sapply(y,length)==0
        x <- x[!omit.idx]; y <- y[!omit.idx]
        auc0 <- auc.cluster(x=x,y=y) 
        contrast <- matrix(c(1,-1),nrow=2)
        alpha <- .05
        I <- length(x)
        z.stat <- with(auc0, sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
        if(paste(status,collapse='')==paste(statuses[[1]],collapse='')) {
            ## print(city.name)
            pdf(paste0('figs/220823_',city.name,'.pdf'))
            plot(auc0,resol=1e2)
            plot(auc0,alpha=.01,add=TRUE)
            dev.off()
        }
        with(auc0,structure(c(theta.12.hat,theta.12.CI,theta.11.hat,theta.11.CI,(1-pnorm(abs(z.stat)))*2,I,length(unlist(x)),length(unlist(y))),
                            names=c('theta.12.hat','theta.12.CI.lower','theta.12.CI.upper','theta.11.hat','theta.11.CI.lower','theta.11.CI.upper','pval','I','M','N')))
    })
}, list(nyc=nyc,boston=boston), c('nyc','boston'))
auc.estimates <- lapply(auc.estimates,t)
## auc.estimates <- lapply(auc.estimates,function(m)cbind(case=status.descriptions,as.data.frame(m)))
auc.estimates <- lapply(auc.estimates,function(m) {
    df <- cbind(case=status.descriptions,as.data.frame(m))
    rounded <- lapply(df[,grep('theta',colnames(df))],function(col)sprintf(col,fmt='%#.2f'))
    df$theta.12.str <- with(rounded, paste0(theta.12.hat,' [',theta.12.CI.lower,', ',theta.12.CI.upper,']'))
    df$theta.11.str <- with(rounded, paste0(theta.11.hat,' [',theta.11.CI.lower,', ',theta.11.CI.upper,']'))
    df$pval.str <- sprintf(df$pval,fmt='%#.2f')
    df
})


## 3d pretty print summaary data tables
require(xtable)
out <- lapply(summary.tables, function(x) {
    x <- x[order(x$group),]
    ## x[,-1] <- round(x[-1],2)
    ## x$mean.duration <- sprintf(x$mean.duration,fmt='%#.2f')
    ## x$sd.duration <- sprintf(x$sd.duration,fmt='%#.2f')
    ## x$freq <- sprintf(x$sd.duration,fmt='%#.2f')
    for(col in c('mean.duration','sd.duration','freq'))
        x[,col] <- sprintf(x[,col],fmt='%#.2f')
    x$mean.duration <- paste0(x$mean.duration,' (',x$sd.duration,')')
    subset(x,select=-sd.duration)
    })
out <- cbind(nyc=out[['nyc']],boston=subset(out[['boston']],select=-group))
colnames(out) <- c('group','mean duration (SD)','count','freq.','mean duration (SD)','count','freq.')
out$group <- c('Asian','Black Hispanic','Black non-Hispanic','other','White Hispanic','White non-Hispanic')
other.idx <- which(out$group=='other')
out <- out[c(1:(other.idx-1),(other.idx+1):nrow(out),other.idx),]
addtorow <- list()
addtorow$pos <- list(0,0,0,0)
addtorow$command <- c("& \\multicolumn{3}{c||}{NYC} & \\multicolumn{3}{c|}{Boston}\\\\\n",
                      '\\hline\n',
                      '\\hline\n',
                      paste0(paste0(colnames(out),collapse=' & '), '\\\\\n'))
align.str <-  "|c|l||C{3cm}|C{1cm}|C{1cm}||C{3cm}|C{1cm}|C{1cm}|"#|C{2cm}|C{2cm}|"
filename <- '220823.tex'
sink(filename)
print.xtable(xtable(out,align=align.str,digits=rep(0,ncol(out)+1)), add.to.row = addtorow,include.colnames=FALSE,include.rownames=FALSE,sanitize.text.function=function(x)x)
sink()
lines <- scan(filename,what='',sep='\n')
start.idx <- grep('begin\\{tabular\\}',lines)
end.idx <- grep('end\\{tabular\\}',lines)
writeLines(lines[start.idx:end.idx],filename)


## 3f data for \input
for(city in c('nyc','boston')) {
    cat(format(nrow(get(city)),big.mark=','),'%',file=paste0('input/da_total_stops_',city,'.txt'))
    apply(auc.estimates[[city]],1,function(row) {
        cat(row['theta.12.str'],'%',sep='',file=tolower(paste0('input/da_',row['case'],'_',city,'_theta12.txt')))
        cat(row['theta.11.str'],'%',sep='',file=tolower(paste0('input/da_',row['case'],'_',city,'_theta11.txt')))
        cat(row['pval.str'],'%',sep='',file=tolower(paste0('input/da_',row['case'],'_',city,'_pval.txt')))
    })
}
## rounded <- lapply(auc.estimates$nyc[,grep('theta',colnames(auc.estimates$nyc))],function(col)sprintf(col,fmt='%#.2f'))
with(auc.estimates$nyc, {
    fmt <- '%#.2f'
    cat(sprintf(fmt,theta.12.hat[case=='Black']),'%',file='input/da_black_nyc_theta12_mean.txt')
    cat(sprintf(fmt,theta.12.CI.lower[case=='Black']),'---',sprintf(fmt,theta.12.CI.upper[case=='Black']),'%',sep='',file='input/da_black_nyc_theta12_CI.txt')
    cat(sprintf(fmt,theta.11.hat[case=='Black']),'%',file='input/da_black_nyc_theta11_mean.txt')
    cat(sprintf(fmt,theta.11.CI.lower[case=='Black']),'---',sprintf(fmt,theta.11.CI.upper[case=='Black']),'%',sep='',file='input/da_black_nyc_theta11_CI.txt')
})



## 3e pretty print table of AUC estimate
require(xtable)
## out <- as.data.frame(round(rbind(t(nyc),t(boston))[c(1,5,2,6,3,7,4,8),],2))
idx <- as.numeric(t(matrix(1:(2*nrow(auc.estimates[[1]])),ncol=2)))
out <- as.data.frame(do.call(rbind,auc.estimates)[idx,])
out$city <- c('NYC','Boston')
out$case.status <- rep(status.descriptions,each=2)
out <- subset(out,select=c(case.status,city,I,M,N,theta.12.str,theta.11.str,pval.str))
out <- within(out,case.status[duplicated(case.status)] <- '')
colnames(out) <- c('case group','data set','I','$\\Sigma M_i$','$\\Sigma N_i$','$\\theta_{12}$','$\\theta_{11}$','$H_{0}:\\theta_{12}=\\theta_{11}$')
addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <- c(    '\\hline\n')
align.str <-  "|c|l||c|c|c|c|c|c|c|"#C{3cm}|C{1cm}|C{1cm}||C{3cm}|C{1cm}|C{1cm}|"
sink('220821.tex')
print.xtable(xtable(out,align=align.str,digits=rep(0,ncol(out)+1)),include.rownames=FALSE,add.to.row=addtorow,sanitize.text.function=function(x)x)
sink()
lines <- scan('220821.tex',what='',sep='\n')
start.idx <- grep('begin\\{tabular\\}',lines)
end.idx <- grep('end\\{tabular\\}',lines)
writeLines(lines[start.idx:end.idx],'220821.tex')


