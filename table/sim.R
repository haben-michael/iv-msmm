args <- commandArgs(trailingOnly=TRUE)
T <- as.numeric(args[1])
n <- as.numeric(args[2])
reps <- as.numeric(args[3])

require(parallel,quietly=TRUE)
source('iv.msmm.R')



delta <- function(A,Z,L,alpha)  pnorm(alpha['const']+alpha['L']*L) 
f.Z <- function(A,Z,L,gamma) plogis(gamma['const'] + gamma['L']*L)

sd.noise <- 1
A.space <- c(0,1)
beta <- rep(1,T)
sd.L <- 1
nu.bar <- structure(c(-0.16, 0.11, 0.26), names=c('const','L','U'))
alpha <- structure(c(0.96, 0.46), names=c('const','L'))
gamma <- structure(c( 0.53, -0.10    ), names=c('const','L'))
nu <- nu.bar[1:2]/sqrt(1+nu.bar['U']^2)
b.U <- rep(1,T) / sqrt(T)
b.L <- rep(1,T) / sqrt(T)

## fits <- mclapply(1:reps, mc.cores=detectCores()-2, FUN=function(jj) {
fits <- lapply(1:reps, FUN=function(jj) {
    U <- matrix(rnorm(n*T),nrow=n)
    A <- Z <- L <- A.mean <- L.mean <- matrix(nrow=n,ncol=T)
    L.mean <- cbind(1,L.mean)
    for(t in 1:T) {
        L[,t] <- L.mean[,t] + U[,t]*sd.L
        delta.cached <- delta(NA,NA,L[,t],alpha=alpha)
        Z[,t] <- rbinom(n,1,prob=f.Z(NA,Z=1,L=L[,t],gamma))
        A.mean[,t] <- pnorm(nu.bar['const'] + nu.bar['L']*L[,t] + nu.bar['U']*U[,t])*(1-delta.cached) + Z[,t]*delta.cached
        A[,t] <- rbinom(n,1,prob=A.mean[,t])
        L.mean[,t+1] <- A[,t]
    }
    L.mean <- L.mean[,-(T+1)]
    Y <- U%*%b.U + (L-L.mean)%*%b.L + A%*%beta + rnorm(n)

    nui.hat <- nuisance.est.1(A,Z,L)
    alpha.hat <- nui.hat$alpha; nu.hat <- nui.hat$nu; gamma.hat <- nui.hat$gamma
    delta.hat <- delta(A,Z,L,alpha.hat)

    causal=tryCatch(iv.msmm(Y,A,Z,L,design=A,f.Z=f.Z,
                                delta=delta,
                                nuisance.est=function(A,Z,L)list(alpha=alpha.hat,beta=NULL,nu=nu.hat,gamma=gamma.hat),
                                stabilization=(-1)^(1-A),get.sandwich=TRUE,get.bootstrap=TRUE), error=function(e)NA)
})
fits <- simplify2array(fits)

target.idx <- T
beta.causal <- fits[grep('beta.hat',rownames(fits))[target.idx],]
sd.sandwich <- fits[grep('sd.sandwich',rownames(fits))[target.idx],]
sd.bootstrap <- fits[grep('sd.bootstrap',rownames(fits))[target.idx],]
levels <- seq(0,1,length.out=21)
levels <- .05
coverages <- lapply(list(sandwich=sd.sandwich,bootstrap=sd.bootstrap),function(sd.hat) {
    coverage <- sapply(levels,function(level) {
        q <- qnorm(1-level/2)
        upper <- fits[target.idx,] + q*sd.hat
        lower <- fits[target.idx,] - q*sd.hat
        mean((upper > beta[target.idx]) & (lower < beta[target.idx]))
    })
    names(coverage) <- round(1-levels,2)
    coverage
})

coverages$bootstrap.empirical <- mean(apply(t(fits[grep('CI.bootstrap.empirical',rownames(fits)),][c(target.idx*2-1,target.idx*2),]) - fits[target.idx,],1,prod) < 0)
coverages$bootstrap.pct <- mean(apply(t(fits[grep('CI.bootstrap.pct',rownames(fits)),][c(target.idx*2-1,target.idx*2),]) - fits[target.idx,],1,prod) < 0)

filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(reps,n,T,target.idx,beta,beta.causal,coverages,sd.sandwich=sd.sandwich,sd.bootstrap=sd.bootstrap,file=filename)
