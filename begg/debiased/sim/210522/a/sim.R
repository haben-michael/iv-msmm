start.time <- Sys.time()


args <- commandArgs(trailingOnly=TRUE)
B <- as.numeric(args[1])

ns <- round(seq(6,3000,len=20))
## rZ <- rnorm
## rZ <- function(n)runif(n,-1/2,1/2)
## rZ <- function(n)rbeta(n,.2,.2) - 1/2
## rZ <-  function(n)rt(n,df=2.25)
rZs <- list(rnorm, function(n)runif(n,-1/2,1/2), function(n)rbeta(n,.2,.2)-1/2,function(n)rbeta(n,.3,.3)-1/2,function(n)rbeta(n,.4,.4)-1/2,function(n)rbeta(n,.5,.5)-1/2,function(n)rbeta(n,.6,.6)-1/2,function(n)rt(n,df=1),function(n)rt(n,df=1.5),function(n)rt(n,df=2),function(n)rt(n,df=2.25),function(n)rt(n,df=2.5),function(n)rt(n,df=3))

by.rZ <- lapply(rZs, function(rZ) {
    probs.full <- probs.2 <- probs.1 <- numeric()
    ## probs <- matrix(nrow=2,ncol=0)    ##medians <- mclapply(ns, mc.cores=detectCores()-2,FUN=function(n) {     
    ## medians <- sapply(ns, FUN=function(n) {
    for(n in ns) {
        ## print(n)
        z <- matrix(rZ(B*n),ncol=n)
        s <- matrix(runif(B*n),ncol=n)
        theta.hat <- rowSums(z*s)/rowSums(s^2)
        ## theta.del  <- rowSums((z*s)[,5:n])/rowSums(s[,1:n]^2)
        theta.ind <- rowSums((z*s)[,5:n])/rowSums(s[,5:n]^2)
        x <- (z[,1]-z[,2])/(s[,1]-s[,2])
        y <- (z[,3]-z[,4])/(s[,3]-s[,4])
        term1 <- x*y - theta.ind*(x+y)
        term2 <- (theta.ind-theta.hat)*(x+y)+theta.hat^2
        probs.full <- cbind(probs.full,mean(term1 + term2 < 0))
        probs.1 <- cbind(probs.1, mean(term1 < 0))
        probs.2 <- cbind(probs.2, mean(abs(term1)<abs(term2))) 
   }
    list(probs.full=probs.full,probs.1=probs.1,probs.2=probs.2)
})

## by.rZ.full <- sapply(by.rZ, function(lst)lst[['probs.full']])
## by.rZ.gt <- sapply(by.rZ, function(lst)lst[['probs.2']])

## op <- par(mfrow=c(4,3))
## for(j in 2:length(rZs)) {
##     plot(ns,(by.rZ.full[,j]-1/2)*ns,main=deparse(rZs[[j]])[2])
##     abline(h=0,lty=2)
##     }
## par(op)

## op <- par(mfrow=c(4,3))
## for(j in 1:(length(rZs)-1)) {
##     plot(ns,by.rZ.gt[,j]*ns^(1),main=deparse(rZs[[j]])[2])
##     abline(h=0,lty=2)
##     }
## par(op)


filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(ns,B,rZs,by.rZ,file=filename)
print(Sys.time() - start.time)
