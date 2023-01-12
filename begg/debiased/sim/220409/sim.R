start.time <- Sys.time()


args <- commandArgs(trailingOnly=TRUE)
B <- as.numeric(args[1])


## 2c-3. FPR with t distr. df=1, maybe df=1.25 cases seem to have
## unbounded bias. More pronounced with unif[0,1] than unif[1,2] for S
## distribution, seems to be some interaction. [update 22-2: shouldn't
## z be scaled to unit variance? moving to sim/220409.]
source('../../../2/misc.R')
n.rZ <- 5
## Z.params <- seq(2.5,1,len=n.rZ)
Z.params <- seq(2.5,1,len=n.rZ)
Z.params <- seq(2.1,2.01,len=n.rZ)
Z.params <- seq(1e-2,1e-6,len=n.rZ)
rZs <- lapply(Z.params,function(epsilon)function(n)rt(n,df=2+epsilon)*sqrt(epsilon/(2+epsilon)))
require(parallel)
ns <- round(seq(1e3,5e3,len=20))
## B <- 1e1
## by.n <- sapply(ns, function(n) {
## rZ <- rZs[[3]]
## by.n <- mclapply(ns, mc.cores=detectCores()-3,FUN=function(n) {
## by.Z <- mclapply(rZs, mc.cores=detectCores()-3,FUN=function(rZ) {
by.Z <- lapply(rZs, FUN=function(rZ) {
    ## cat('\n')
    ## by.n <- sapply(ns,FUN=function(n) {
    by.n <- mclapply(ns,mc.cores=detectCores()-4,FUN=function(n) {
        ## cat('.')
        mean(replicate(B, {
            z <- rZ(n) 
            s <- runif(n,1,2)
            ## theta.hat <- rowSums(z*s)/rowSums(s^2)
            ## x <- (z[,1]-z[,2])/(s[,1]-s[,2])
            ## y <- (z[,3]-z[,4])/(s[,3]-s[,4])
            ## cov(x > theta.hat, y > theta.hat)
            begg.test(y=z/s,v=1/s^2)['pval'] < .1
        }))
    })
    by.n <- simplify2array(by.n)
})
by.Z <- simplify2array(by.Z)

## png('neg_bias_students.png')
## colors <- rev(gray.colors(n.rZ))
## matplot(ns,by.Z,type='l',col=colors,lty=1,xlab='no. of studies',ylab='FPR')
## abline(h=.1,lty=2)
## legend('topleft',col=colors,lty=1,legend=Z.params,title='df')
## dev.off()


filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(B,ns,by.Z,n.rZ,Z.params,file=filename)#ns,B,rZs,by.rZ,file=filename)
print(Sys.time() - start.time)
