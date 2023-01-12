## 1. replication of tables in binormal vs lehmann manuscript. There
## was a replication in 210923.R but it was producing strange results
## I think due to the weibull generation. This one gave auc.binormal
## the same interface as a auc.cox and auc.mle from 210924.R and used
## the weibull generation from there.

require(survival)
require(parallel)
## require(dplyr)
## require(ggplot2)
require(xtable)

auc.cox <- function(y.0,y.1,alpha=.05) {
    y <- c(y.0,y.1)
    d <- rep(c(0,1),c(length(y.0),length(y.1)))
    fit <-         coxph(Surv(y) ~ d)
    log.theta.hat <- unname(coef(fit))
    se.hat <- sqrt(vcov(fit))
    q <- qnorm(1-alpha/2)
    CI.lower <- log.theta.hat-q*se.hat
    CI.upper <- log.theta.hat+q*se.hat
    auc.hat <- 1/(exp(log.theta.hat)+1)
    ## est.cox <- c(auc.hat=auc.hat,se.hat=se.hat,CI.lower=1/(exp(CI.upper)+1),CI.upper=1/(exp(CI.lower)+1))
}

gen.data.weibull <- function(m,n,auc,shape=1,scale.0=1) {
    scale.1 <- scale.0*(auc/(1-auc))^(1/shape)
    y.0 <- rweibull(m,shape=shape,scale=scale.0)#^(1/shape))
    y.1 <- rweibull(n,shape=shape,scale=scale.1)#^(1/shape))
    return(list(y.0=y.0,y.1=y.1))
}


## auc.binormal <- function(y,d,alpha=.05) {
auc.binormal <- function(y.0,y.1,alpha=.05) {
    y <- c(y.0,y.1)
    d <- rep(c(0,1),c(length(y.0),length(y.1)))
    n <- sum(d) # assumes same # of diseased as non-diseased
    y.0 <- y[d==0]; y.1 <- y[d==1]
    s0 <- sd(y.0); s1 <- sd(y.1)
    auc.hat <- pnorm((mean(y.1)-mean(y.0))/sqrt(s0^2+s1^2))        
    ## var.hat <- 1/n + (mean(y.1)-mean(y.0))^2/2/(s0^2+s1^2)^3*(s0^4+s1^4)/(n-1)
    ## q <- qnorm(1-alpha/2)
    ## CI.lower <- qnorm(auc.hat)-q*sqrt(var.hat)
    ## CI.upper <- qnorm(auc.hat)+q*sqrt(var.hat)
    ## est.binormal <- c(auc.hat=auc.hat,se.hat=sqrt(var.hat),CI.lower=pnorm(CI.lower),CI.upper=pnorm(CI.upper))
}


reps <- 1e3
ns <- c(5,10,15,30,60,100,200)#seq(10,40,by=10)
aucs <- AUCs <- c(.55,.7,.9)#seq(.5,.9,by=.1)



by.aucs <- lapply(list(function(m,n,auc)gen.data.weibull(n,n,auc,shape=1,scale.0=1),
                       function(m,n,auc)gen.data.weibull(n,n,auc,shape=3.5,scale.0=1),
                       function(m,n,auc)gen.data.weibull(n,n,auc,shape=20,scale.0=1)), function(gen.data) {
                by.auc <- lapply(aucs, function(auc) {
                    ## by.n <- mclapply(ns, mc.cores=detectCores()-4, function(n) {
                    by.n <- lapply(ns, function(n) {
                        auc.hats <- replicate(reps, {
                            y <- gen.data(n,n,auc)
                            y.0 <- y$y.0; y.1 <- y$y.1
                            y <- unlist(y)
                            c(binormal=auc.binormal(y.0,y.1),cox=auc.cox(y.0,y.1))
                        })
                        bias <- rowMeans(auc.hats)-auc
                        rbind(bias=bias, mse=bias^2+apply(auc.hats,1,var))
                    })
                    by.n <- simplify2array(by.n)
                })
                by.auc <- simplify2array(by.auc)
                dimnames(by.auc)[[3]] <- ns
                dimnames(by.auc)[[4]] <- aucs
                by.auc
})






lapply(by.aucs, function(by.auc) {
    ftbl <- ftable(by.auc,row.vars=c(3,1),col.vars=c(4,2))
    ftbl <- round(ftbl,4)
    xftbl <- xtableFtable(ftbl,method='compact',digits=3)
    sink('211220_tables.tex',append=TRUE)
    print.xtableFtable(xftbl, booktabs = TRUE)    
    sink()
})





## 2. referee comments
## 21-12-21: referee: In Sec. 3.2, Tables 1-4, it may be helpful to
## add some measure of relative efficiency (RE) of the two approaches
## (ratio of MSEs?) Furthermore, RE could be visualized graphically as
## a function of the sample size to gain better insight into relative
## performance of the two methods.



png('211220.png')
op <- par(mfrow=c(1,3))
lapply(by.aucs, function(by.auc) {
ratios <- by.auc[2,1,,] / by.auc[2,2,,]
matplot(x=as.numeric(rownames(ratios)),y=ratios,type='l',lty=1,col=1:3,xlab='sample size',ylab='ratio of MSEs',main='MSE of binormal / MSE of Cox')
legend('bottomright',lty=1,col=1:3,legend=colnames(ratios))
})
par(op)
dev.off()
