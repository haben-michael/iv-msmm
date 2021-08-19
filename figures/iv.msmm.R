lm.wt <- function(Y,W,design) {
    as.numeric(solve(t(design)%*%(design/W))%*%(t(design)%*%(Y/W)))
}
## probit model 
nuisance.est.1 <- function(A,Z,L,init=rep(0,4)) {
    A.pooled <- as.numeric(A)
    Z.pooled <- as.numeric(Z)
    L.pooled <- cbind(1,as.numeric(L))
    mle <- optim(init,fn=function(par){
        alpha <- par[1:2]; nu <- par[3:4]
        delta.par <- pnorm(L.pooled%*%alpha)
        A.mean.par <- pnorm(L.pooled%*%nu)*(1-delta.par) + Z.pooled*delta.par
        loglik <- A.pooled*log(A.mean.par) + (1-A.pooled)*log(1-A.mean.par)
        -sum(loglik)
    })$par
    gamma <- structure(coef(glm(as.numeric(Z) ~ as.numeric(L), family='binomial')), names=c('const','L'))
    list(alpha=structure(mle[1:2],names=c('const','L')),beta=NULL,gamma=gamma,nu=structure(mle[3:4],names=c('const','L')))
}


nuisance.inf.1 <- function(A,Z,L,alpha.hat,beta.hat,nu.hat,gamma.hat,score.beta) { ## write interface
    n <- nrow(A); T <- ncol(A); p <- length(beta.hat)
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
    I.nui.22 <- t(dd*(1-ppaa)*ddnn*as.numeric(-L.pooled%*%nu.hat)*L.pooled)%*%L.pooled
    A.mean.deriv <- cbind((Z.pooled-ppnn)*ddaa*L.pooled,(1-ppaa)*ddnn*L.pooled)
    I.nui.3 <- t(as.numeric(A.pooled/A.mean.mle^2 + (1-A.pooled)/(1-A.mean.mle)^2)*A.mean.deriv)%*%A.mean.deriv
    I.nui.alphanu <- rbind(cbind(I.nui.11,I.nui.12),cbind(I.nui.21,I.nui.22)) - I.nui.3
    I.nui.gamma <- with(list(sigma=as.numeric(plogis(L.pooled%*%gamma.hat))), t(L.pooled)%*%diag(sigma*(1-sigma))%*%L.pooled)
    I.nui <- rbind(cbind(I.nui.alphanu,matrix(0,nrow=nrow(I.nui.alphanu),ncol=ncol(I.nui.gamma))), cbind(matrix(0,ncol=ncol(I.nui.alphanu),nrow=nrow(I.nui.gamma)),I.nui.gamma))
    score.alphanu <- t(apply(array(dd*A.mean.deriv,dim=c(n,T,ncol(A.mean.deriv))), c(1,3), sum))
    score.gamma <- t(with(list(d.scalar=Z - plogis(gamma.hat['const'] + gamma.hat['L']*L)), cbind(rowSums(d.scalar),rowSums(L*d.scalar))))
    score.nui <- rbind(score.alphanu,score.gamma)
    I.beta.alpha <- with(list(ee=alpha.hat['const'] + alpha.hat['L']*L), { d.scalar=dnorm(ee)/pnorm(ee); cbind(rowSums(d.scalar),rowSums(L*d.scalar))})
    I.beta.alpha <- score.beta%*%I.beta.alpha
    I.beta.gamma <- score.beta%*%t(score.gamma)#I.beta.gamma
    I.beta.nu <- matrix(0,nrow=p,ncol=length(nu.hat))
    list(I.nuisance=I.nui,I.beta.alpha=I.beta.alpha,I.beta.gamma=I.beta.gamma,I.beta.nu=I.beta.nu,score.nuisance=score.nui)
}

## markov model
nuisance.est.2 <- function(A,Z,L) {
    A.pooled <- as.numeric(A)
    Z.pooled <- as.numeric(Z)
    L.pooled <- as.numeric(L)
    mle <- optim(c(rep(1/4,2),rep(.5,1)),control=list(trace=0,maxit=2e3),function(par){
        alpha <- par[1:2]; nu <- par[3]# alpha aka (\delta_0,\delta_1), nu aka p_L
        delta.par <- (-1)^(1-Z.pooled)*alpha[L.pooled+1]
        A.mean.par <- .5*(p.B + 2*(1-p.B)*nu[1]^(L.pooled==1)*(1-nu[1])^(L.pooled!=1)) + delta.par
        if(min(A.mean.par)<=0 | max(A.mean.par)>=1) {
            return(jitter(.Machine$double.xmax/1e4))
        }
        loglik <- A.pooled*log(A.mean.par) + (1-A.pooled)*log(1-A.mean.par)
        return(-sum(loglik))
    })$par
    list(alpha=mle[1:2],beta=NULL,gamma=NULL,nu=mle[3])
}

nuisance.inf.2 <- function(A,Z,L,alpha.hat,beta.hat,nu.hat,gamma.hat,score.beta) { ## write interface
    A.pooled <- as.numeric(A)
    Z.pooled <- as.numeric(Z)
    L.pooled <- as.numeric(L)
    A.mean.mle <- .5*(p.B + 2*(1-p.B)*nu.hat^(L.pooled==A.pooled)*(1-nu.hat)^(L.pooled!=A.pooled))
    A.mean.deriv <- cbind((-1)^Z.pooled*(1-L.pooled)*(alpha.hat[2]/alpha.hat[1])^L.pooled,(-1)^(1-Z.pooled)*L.pooled*(alpha.hat[1]/alpha.hat[2])^(1-L.pooled),(1-p.B)*L.pooled*(1/nu.hat-1)^(1-L.pooled) - (1-p.B)*(1-L.pooled)*(nu.hat/(1-nu.hat))^L.pooled)
    score.nui <- with(list(dd=A.pooled/A.mean.mle - (1-A.pooled)/(1-A.mean.mle)),t(apply(array(dd*A.mean.deriv,dim=c(n,T,ncol(A.mean.deriv))), c(1,3), sum)))
    I.nui <- t((- A.pooled/A.mean.mle^2 - (1-A.pooled)/(1-A.mean.mle)^2) * A.mean.deriv)%*%A.mean.deriv
    I.beta.alpha <- cbind(rowSums((1-L)*alpha.hat[1]),rowSums(L*alpha.hat[2]))
    I.beta.alpha <- score.beta%*%I.beta.alpha
    I.beta.nu <- matrix(0,nrow=length(beta),ncol=length(nu.hat))
    list(I.nuisance=I.nui,I.beta.alpha=I.beta.alpha,I.beta.gamma=NULL,I.beta.nu=I.beta.nu,score.nuisance=score.nui)
}



iv.msmm <- function(Y,A,Z,L,design,f.Z,delta,get.sandwich=FALSE,get.bootstrap=FALSE,nuisance.est=nuisance.est.1,nuisance.inf=nuisance.inf.1,stabilization=NULL,bootstrap.reps=1e2) {
    T <- ncol(A); n <- nrow(A); p <- ncol(design) # rename J -> T
    A.ancillary <- sapply(1:T,function(j) glm(A[,j] ~ cbind(1,A[,seq_len(j-1)]),family=binomial)$fitted)
    nui.hats <- nuisance.est(A,Z,L)
    alpha.hat <- nui.hats[['alpha']];  beta.hat <- nui.hats['beta']; nu.hat <- nui.hats[['nu']]; gamma.hat <- nui.hats[['gamma']]
    delta.cached <- delta(A,Z,L,alpha.hat)
    W.causal <- f.Z(A,Z,L,gamma.hat)*delta.cached*(-1)^(1-Z)*stabilization /(A.ancillary^A*(1-A.ancillary)^(1-A))
    W.causal <- apply(W.causal,1,prod)
    beta.hat <- lm.wt(Y,W.causal,design)
    sd.sandwich = NA; sd.bootstrap = NA; CI.bootstrap <- NA

    if(get.sandwich) {
        score.beta <- t(as.numeric(((Y-design%*%beta.hat)/W.causal))*design)
        I.nui <- nuisance.inf(A=A,L=L,Z=Z,alpha.hat=alpha.hat,beta.hat=beta.hat,gamma.hat=gamma.hat,nu.hat=nu.hat,score.beta=score.beta)
        score <- rbind(score.beta,I.nui$score.nui)
        I.beta.beta <- t(-(1/W.causal)*design)%*%design 
        I.top <- cbind(I.beta.beta,I.nui$I.beta.alpha,I.nui$I.beta.nu,I.nui$I.beta.gamma)
        I.bottom <- cbind(matrix(0,nrow=nrow(I.nui$I.nui),ncol=p),I.nui$I.nui)
        I <- rbind(I.top,I.bottom)
        extra.params <- rowSums(I)==0 & colSums(I)==0                                        
        if (sum(extra.params)>0) {
            warning('Dropping extra parameters from information matrix')
            I <- I[!extra.params,!extra.params]
            score <- score[!extra.params,]
        }
        I.inv <- solve(I)
        sd.sandwich <- sqrt(diag(I.inv%*%(score%*%t(score))%*%t(I.inv))[1:p])
    }

    if(get.bootstrap) {
        bootstrap.ests <- replicate(bootstrap.reps, {
            idx <- sample(n,replace=TRUE)
            lm.wt(Y=Y[idx],W=W.causal[idx],design=design[idx,])
        })
        if(p==1)bootstrap.ests <- matrix(bootstrap.ests,nrow=1) ## clean up
        sd.bootstrap <- apply(bootstrap.ests,1,sd)
    }
    return(c(beta.hat=beta.hat,sd.sandwich=sd.sandwich,sd.bootstrap=sd.bootstrap))
}



## iv.msmm(Y,A,Z,L,design=A,f.Z=f.Z,
##         delta=delta,
##         nuisance.est=function(A,Z,L)list(alpha=alpha,beta=NULL,nu=nu),
##         stabilization=(-1)^(1-A))
