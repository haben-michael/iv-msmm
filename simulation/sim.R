require(parallel)
require(abind)
mc.cores <- detectCores()
print(mc.cores)

args <- commandArgs(trailingOnly=TRUE)
J <- as.numeric(args[1])
n <- as.numeric(args[2])
## reps <- as.numeric(args[3])
print(J)
print(n)

trim <- function(v,p){
    q <- quantile(v,c(p/2,1-p/2))
    v[v>q[1] & v<q[2]]
}
msm.design <- function(A) {
    design <- cbind(1,rowSums(A))
}
reps <- 5e2
bootstrap.reps <- 1e3
beta <- rep(1,J)
gamma <- c(-.3,.3)
alpha <- c(-.2,L=.2,U=.2)
beta <- rep(1,J)
par <- c(beta,alpha,gamma)
Y.coefs <- (1:J)#rep(1,J)
L.coefs <- c(0,lapply(seq_len(J-1),function(j)c(-.05/2,.05)))
Z.mean <- .5


ipw.est <- function(Y,W,design) {
    as.numeric(solve(t(design)%*%(design/W))%*%(t(design)%*%(Y/W)))
}

est.nuisance.par <- function(A,L,Z) { ## inherit--interface for getting estimates nuisance params
    A.pooled <- as.numeric(A)
    Z.pooled <- as.numeric(Z)
    L.pooled <- cbind(1,as.numeric(L))
    mle <- optim(c(0,0,0,0),fn=function(par){
        alpha <- par[1:2]; nu <- par[3:4]
        delta.par <- pnorm(L.pooled%*%alpha)
        A.mean.par <- as.numeric(pnorm(L.pooled%*%nu)*(1-delta.par) + Z.pooled*delta.par)
        loglik <- A.pooled*log(A.mean.par) + (1-A.pooled)*log(1-A.mean.par)
        -sum(loglik)
    })$par
    list(alpha=mle[1:2],beta=NULL,nu=mle[3:4])
}

## Estimate part of score and information matrix corresponding to the
## nuisance parameters \alpha, \beta, \nu, assuming linear model MSMM
## described in manuscript simulation section.
est.nuisance.inf <- function(A,Z,L,alpha.hat,beta.hat,nu.hat,score.beta) { ## write interface
    A.pooled <- as.numeric(A)
    Z.pooled <- as.numeric(Z)
    L.pooled <- cbind(1,as.numeric(L))
    delta.mle <- pnorm(L.pooled%*%alpha.hat)
    A.mean.mle <- pnorm(L.pooled%*%nu.hat)*(1-delta.mle) + Z.pooled*delta.mle
    dd <- as.numeric(A.pooled/A.mean.mle - (1-A.pooled)/(1-A.mean.mle)) 
    ddnn <- as.numeric(dnorm(L.pooled%*%nu.hat))
    ddaa <- as.numeric(dnorm(L.pooled%*%alpha.hat))
    ppnn <- as.numeric(pnorm(L.pooled%*%nu.hat))
    ppaa <- as.numeric(pnorm(L.pooled%*%alpha.hat))
    I.nui.11 <- t(dd*(Z.pooled-ppnn)*ddaa*as.numeric(-L.pooled%*%alpha.hat)*L.pooled)%*%L.pooled
    I.nui.12 <- -t(dd*ddnn*ddaa*L.pooled)%*%L.pooled
    I.nui.21 <- -t(dd*ddaa*ddnn*L.pooled)%*%L.pooled 
    I.nui.22 <- t(dd*as.numeric(1-pnorm(L.pooled%*%alpha.hat))*ddnn*as.numeric(-L.pooled%*%nu.hat)*L.pooled)%*%L.pooled
    A.mean.deriv <- cbind((Z.pooled-ppnn)*ddaa*L.pooled,(1-ppaa)*ddnn*L.pooled)
    I.nui.3 <- t(as.numeric(A.pooled/A.mean.mle^2 + (1-A.pooled)/(1-A.mean.mle)^2)*A.mean.deriv)%*%A.mean.deriv
    I.nui <- rbind(cbind(I.nui.11,I.nui.12),cbind(I.nui.21,I.nui.22)) - I.nui.3
    score.nui <- t(apply(array(dd*A.mean.deriv,dim=c(n,J,ncol(A.mean.deriv))), c(1,3), sum))
    I.beta.alpha <- with(list(ee=alpha.hat[1] + alpha.hat[2]*L), { d.scalar=dnorm(ee)/pnorm(ee); cbind(rowSums(d.scalar),rowSums(L*d.scalar))})
    I.beta.alpha <- score.beta%*%I.beta.alpha
    I.beta.nu <- matrix(0,nrow=length(beta),ncol=length(nu.hat))
    list(I.nuisance=I.nui,I.beta.alpha=I.beta.alpha,I.beta.gamma=NULL,I.beta.nu=I.beta.nu,score.nuisance=score.nui)
}

fits <- mclapply(1:reps,mc.cores=mc.cores-2,  FUN=function(i){
    L <- A <- delta <- A.mean <- L.mean <- delta.fitted <- delta.ancillary <- matrix(nrow=n,ncol=J)
    U <- matrix(rnorm(n*J),nrow=n)
    Z <- matrix(rbinom(n*J,1,.5),nrow=n)
    for(j in 1:J) {
        L.mean[,j] <- cbind(1,A[,j-1])%*%L.coefs[[j]] 
        L[,j] <- L.mean[,j] + rnorm(n)/5
        delta[,j] <- pnorm(cbind(1,L[,j])%*%gamma)
        delta.ancillary[,j] <- delta[,j]
        A.mean[,j] <- pnorm(cbind(1,L[,j],U[,j])%*%alpha)*(1-delta[,j]) + Z[,j]*delta[,j]#;hist(A.mean)
        A[,j] <- rbinom(n,1,prob=A.mean[,j])
    }
    A.ancillary <- sapply(1:J,function(j) glm(A[,j] ~ cbind(1,A[,seq_len(j-1)]),family=binomial)$fitted)
    design <- msm.design(A)
    Y <- 3*(L[,1]-L.mean[,1]) + U%*%Y.coefs + design%*%beta + rnorm(n)    
    nui.hats <- est.nuisance.par(A,L,Z)
    alpha.hat <- nui.hats[['alpha']];  beta.hat <- nui.hats['beta']; nu.hat <- nui.hats[['nu']]
    delta.fitted  <- pnorm(alpha.hat[1] + alpha.hat[2]*L)
    W.causal <- Z.mean^Z*(1-Z.mean)^(1-Z)*delta.fitted/(-1)^(2-Z-A)/(A.ancillary^A*(1-A.ancillary)^(1-A))/delta.ancillary
    W.causal <- apply(W.causal,1,prod)
    beta.ipw <- ipw.est(Y,W.causal,design)
      
    score.beta <- t(as.numeric(((Y-design%*%beta.ipw)/W.causal))*design)
    I.nui <- est.nuisance.inf(A=A,L=L,Z=Z,alpha.hat=alpha.hat,beta.hat=NULL,nu.hat=nu.hat,score.beta=score.beta)
    score <- rbind(score.beta,I.nui$score.nui)
    I.beta.beta <- t(-(1/W.causal)*design)%*%design # design=del(msm)/del(beta)
    I.top <- cbind(I.beta.beta,I.nui$I.beta.alpha,I.nui$I.beta.gamma,I.nui$I.beta.nu)
    I.bottom <- cbind(matrix(0,nrow=nrow(I.nui$I.nui),ncol=length(beta)),I.nui$I.nui)
    I <- rbind(I.top,I.bottom)
    I.inv <- solve(I)
    sd.sandwich <- sqrt(diag(I.inv%*%(score%*%t(score))%*%t(I.inv))[1:length(beta)])

    bootstrap.ests <- replicate(bootstrap.reps, {
        idx <- sample(n,replace=TRUE)
        ipw.est(Y=Y[idx],W=W.causal[idx],design=design[idx,])
    })
    sd.bootstrap <- apply(bootstrap.ests,1,function(v)sd(trim(v,p.trim)))
    c(beta.ipw=beta.ipw,sd.sandwich=sd.sandwich,sd.bootstrap=sd.bootstrap)
})
fits <- simplify2array(fits)


target.idx <- 2
beta.causal <- fits[target.idx,]
sd.sandwich <- fits[length(beta)+target.idx,]
sd.bootstrap <- fits[2*length(beta)+target.idx,]
levels <- seq(0,1,length.out=21)
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


filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(reps,n,J,target.idx,beta,beta.causal,coverages,sd.sandwich=sd.sandwich,sd.bootstrap=sd.bootstrap,file=filename)
