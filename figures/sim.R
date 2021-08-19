## (060721) fixing issues found when preparing jasa revision 
## (061921) adding correlation among the U, Z, normalizing b.U, b.L

args <- commandArgs(trailingOnly=TRUE)
T <- as.numeric(args[1])
n <- as.numeric(args[2])
reps <- as.numeric(args[3])
source('iv.msmm.R')
require(parallel,quietly=TRUE)
require(abind,quietly=TRUE)
delta <- function(A,Z,L,alpha)  pnorm(alpha['const']+alpha['L']*L) 
f.Z <- function(A,Z,L,gamma) .5^Z*(1-.5)^(1-Z)
f.Z <- function(A,Z,L,gamma) plogis(gamma['const'] + gamma['L']*L)

sd.noise <- 1
A.space <- c(0,1)
beta <- rep(1,T)
sd.L <- 1/3 ## change to lambda

nu.bar <- structure(c(0, .5 , .41  ), names=c('const','L','U'))  
alpha <- structure(c( 0.23, 0.18   ), names=c('const','L'))
gamma <- structure(c(0.97, 0.35   ), names=c('const','L'))
nu <- nu.bar[1:2]/sqrt(1+nu.bar['U']^2)
b.U <- rep(1,T) / sqrt(T)
b.L <- rep(1,T) / sqrt(T)

beta.hats <- lapply(1:reps, FUN=function(jj) {
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
    A.mean.obs <- pnorm(nu.hat['const'] + nu.hat['L']*L)*(1-delta.hat) + Z*delta.hat
    W.sra <- A.mean.obs^A*(1-A.mean.obs)^(1-A)
    W.sra <- apply(W.sra,1,prod)
    W.oracle <- A.mean^A*(1-A.mean)^(1-A)
    W.oracle <- apply(W.oracle,1,prod)

    rbind(
        assoc=unname(coef(lm(Y ~ A-1))),
        sra=lm.wt(Y,W.sra,design=A),
        causal=tryCatch(iv.msmm(Y,A,Z,L,design=A,f.Z=f.Z,
                                delta=delta,
                                nuisance.est=function(A,Z,L)list(alpha=alpha.hat,beta=NULL,nu=nu.hat,gamma=gamma.hat),
                                stabilization=(-1)^(1-A))[1:T], error=function(e)NA),
        oracle=lm.wt(Y,W.oracle,design=A)            
    )
})
beta.hats <- simplify2array(beta.hats)
filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(beta.hats,beta,T,n,reps,file=filename)
