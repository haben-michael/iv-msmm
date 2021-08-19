require(parallel,quietly=TRUE)
require(abind,quietly=TRUE)
source('iv.msmm.R')
start.time <- Sys.time()

args <- commandArgs(trailingOnly=TRUE)
J <- T <- as.numeric(args[1])
n <- as.numeric(args[2])
reps <- as.numeric(args[3])

lambda <- 1/sqrt(T)
L.coefs <- rep(lambda,T)
U.coefs <- rep(lambda,T)
sigma <- 1/4
grid <- expand.grid(0:1,0:1)
beta <- 1#runif(1,0,5)
Z.mean <- 1/2

p.L <- 0.22
p.U <- 0.07
q <- p.B <- 0.38 # clean up
alpha <- runif(2,-.5,.5)

f.LU.A <- function(data,t,n=nrow(data$A))if(t>1) with(data, ((1-q)*dbinom(L[,t]==A[,t-1],1,p.L) + q*dbinom(U[,t]==A[,t-1],1,p.U))/2 ) else matrix(1/4,nrow=n)
f.A.LU <- function(data,t)with(data, (1-q)*dbinom(L[,t]==A[,t],1,p.L)+q*dbinom(U[,t]==A[,t],1,p.U))
f.A.ZLU <- function(data,t)with(data, f.A.LU(data,t) + (2*Z[,t]-1)*(2*A[,t]-1)*alpha[L[,t]+1]/2 )
f.L.A <- function(data,t) if(t>1) (1-q)*dbinom(data$L[,t]==data$A[,t-1],1,p.L) + q/2 else 1/2
E.L.A <- function(data,t,n=nrow(data$A)) if(t>1) (1-q)*dbinom(data$A[,t-1],1,p.L) + q/2 else matrix(1/2,nrow=n)
E.U.A <- function(data,t,n=nrow(data$A)) if(t>1) q*dbinom(data$A[,t-1],1,p.U) + (1-q)/2 else matrix(1/2,nrow=n)
f.L.A <- function(data,t) { m <- E.L.A(data,t); m^data$L[,t] * (1-m)^(1-data$L[,t])}
f.U.AL <- function(data,t) f.LU.A(data,t)/f.L.A(data,t)
f.U.future <- function(data,t)     f.A.ZLU(data,t)*f.LU.A(data,t) / (f.A.AZL(data,t)*f.L.A(data,t)) # ie f(U_t | A_{t-1},AZL_t)
f.A.AZL <- function(data,t,n=nrow(data$A),T=ncol(data$A)) with(data,rowSums(do.call(cbind,sapply(0:1,function(u)f.A.ZLU(list(A=A,Z=Z,L=L,U=matrix(u,nrow=n,ncol=T)),t)*f.U.AL(list(A=A,Z=Z,L=L,U=matrix(u,nrow=n,ncol=T)),t),simplify=FALSE))))
design <- function(A)as.matrix(rowSums(A),ncol=1)
h <- function(A)as.matrix(rowSums(A),nrow=nrow(A))
delta <- function(A,Z,L,alpha.hat) alpha[L+1]*(-1)^(1-A)
f.Z <- function(A,Z,L,gamma) Z.mean^Z*(1-Z.mean)^(1-Z)
fits <- lapply(1:reps,FUN=function(jj) {
    A <- Z <- L <- U <- A.mean <- matrix(nrow=n,ncol=T+2)
    L[,1] <- rbinom(n,1,.5)
    U[,1] <- rbinom(n,1,.5)
    Z <- matrix(rbinom(n*(T+2),1,Z.mean),nrow=n,ncol=T+2)
    B <- matrix(rbinom(n*(T+2),1,q),nrow=n)
    for(t in 1:(T+1)) {
        A.mean[,t] <- f.A.ZLU(data=list(A=cbind(A[,seq_len(t-1)],1),Z=Z,L=L,U=U),t=t)
        stopifnot(all.equal((1-q)*p.L^L[,t]*(1-p.L)^(1-L[,t]) + q*p.U^U[,t]*(1-p.U)^(1-U[,t]) + (-1)^(1-Z[,t])*alpha[L[,t]+1]/2,A.mean[,t]))
        stopifnot(min(A.mean[,t])>0 & max(A.mean[,t])<1)
        A[,t] <- rbinom(n,1,prob=A.mean[,t])
        probs <- apply(grid, 1, function(r) ((1-q)*dbinom(r[1]==A[,t],1,p.L) + q*dbinom(r[2]==A[,t],1,p.U))/2)
        idx <- apply(probs, 1, function(prob)which(1==rmultinom(1,1,prob)))
        L[,t+1] <- grid[idx,1]
        U[,t+1] <- grid[idx,2]
    }
    A.lag <- matrix(A[,1:T],ncol=T)
    A <- matrix(A[,2:(T+1)],ncol=T); L <- matrix(L[,2:(T+1)],ncol=T); U <- U[,2:(T+1)]
    Z <- matrix(Z[,2:(T+1)],ncol=T); A.mean <- matrix(A.mean[,2:(T+1)],ncol=T)
    A.ancillary <- sapply(1:T,function(j) glm(A[,j] ~ cbind(1,A[,seq_len(j-1)]),family=binomial)$fitted)
    A.mean.obs <- sapply(1:T, function(t) {
        covariates <- list(A=A,Z=Z,L=L)
        covariates$A[,t] <- 1
        f.A.AZL(data=covariates,t)
    })
    eta <- (L - sapply(1:T,function(t) E.L.A(list(A=A),t=t))) %*% L.coefs + (U - sapply(1:T,function(t) E.U.A(list(A=A,NA,NA),t=t))) %*% U.coefs
    epsilon <- rnorm(n)*sigma
    Y <- eta + design(A)%*%beta + epsilon
    iv.msmm(Y=Y,A=A,L=L,Z=Z,design=design(A),f.Z=f.Z,delta=delta,stabilization=(-1)^(1-A),get.sandwich=TRUE,get.bootstrap=TRUE,nuisance.est=nuisance.est.2,nuisance.inf=nuisance.inf.2,bootstrap.reps=100)
})
fits <- simplify2array(fits)

beta.causal <- fits['beta.hat',]
sd.sandwich <- fits['sd.sandwich',]
sd.bootstrap <- fits['sd.bootstrap',]
levels <- seq(0,1,length.out=21)
levels <- .05
coverages <- lapply(list(sandwich=sd.sandwich,bootstrap=sd.bootstrap),function(sd.hat) {
    coverage <- sapply(levels,function(level) {
        q <- qnorm(1-level/2)
        upper <- fits['beta.hat',] + q*sd.hat
        lower <- fits['beta.hat',] - q*sd.hat
        mean((upper > beta) & (lower < beta))
    })
    names(coverage) <- round(1-levels,2)
    coverage
})


filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(reps,n,T,beta,beta.causal,coverages,sd.sandwich=sd.sandwich,sd.bootstrap=sd.bootstrap,file=filename)
