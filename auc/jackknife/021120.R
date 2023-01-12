## n <- 50
## ns <- floor(seq(10,1000,length.out=20))
## by.n <- sapply(ns,function(n) {
##     rho.hats <- replicate(5e2, {
##         rho <- .4
##         sigma <- matrix(c(1,0,rho,sqrt(1-rho^2)),byrow=TRUE,nrow=2)
##         data <- sigma%*% rbind(rnorm(n),rnorm(n))
##         x <- data[1,]; y <- data[2,]
##         c(est=mean((x-mean(x))*(y-mean(y))),true=mean(x*y))
##     })
##     apply(rho.hats,1,sd)
## })
## plot(ns,by.n[1,])
## lines(ns,by.n[2,])


require(mvtnorm)
concord.xy <- function(x,y) {
    m <- length(x)
    n <- length(y)
    ## stopifnot(m>1,n>1)
    if(length(x)==0 | length(y)==0)return(0)

    sum(sapply(y,function(y.i)x<y.i))
    ## sum(outer(x,y,`<`))
}



auc.obu <- function(x,y,alpha=.05) {

    I <- length(x)
    stopifnot(length(x)==length(y))
    m <- sapply(x,length)
    n <- sapply(y,length)
    I.10 <- sum(m>0)
    I.01 <- sum(n>0)
    M <- sum(m)
    N <- sum(n)

    theta.hat <- concord.xy(unlist(x),unlist(y)) / (M*N)
    V.10 <- sapply(x,function(x.i)concord.xy(x.i,unlist(y)))/N
    V.01 <- sapply(y,function(y.i)concord.xy(unlist(x),y.i))/M

    S.10 <- sum((V.10 - m*theta.hat)^2) * I.10 / ((I.10-1)*M)
    S.01 <- sum((V.01 - n*theta.hat)^2) * I.01 / ((I.01-1)*N)
    S.11 <- sum((V.10 - m*theta.hat)*(V.01 - n*theta.hat)) * I / (I-1)

    var.hat <- S.10/M + S.01/N + 2*S.11/(M*N)
    q <- qnorm(1-alpha/2)
    return(c(theta.hat=theta.hat,var.hat=var.hat,CI.lower=theta.hat-q*sqrt(var.hat),CI.upper=theta.hat+q*sqrt(var.hat)))

}


auc.jk <- function(x,y,alpha=.05) {
    concord <- function(x,y) mean(outer(x,y,`<`))
    concord.sum <- function(x,y)sum(outer(x,y,`<`))
    m <- sapply(x,length); n <- sapply(y,length)
    M <- sum(m); N <- sum(n)
    I <- length(x); stopifnot(length(x)==length(y))

    theta.hat <- concord(unlist(x),unlist(y))
    theta.del <- sapply(1:I,function(i)concord(unlist(x[-i]),unlist(y[-i])))
    pseudo <- I*theta.hat - (I-1)*theta.del
    theta.jk <- mean(pseudo)
    var.hat <- var(pseudo)/I
    q <- qnorm(1-alpha/2)


    return(c(theta.hat=theta.jk,var.hat=var.hat,CI.lower=theta.jk-q*sqrt(var.hat),CI.upper=theta.jk+q*sqrt(var.hat)))
}


auc.true <- function(delta)
    integrate(function(u)pnorm(u)*dnorm(u,mean=delta),lower=-Inf,upper=Inf)$value

set.seed(1)
I <- 4e2
theta.true <- .7
delta <- uniroot(Vectorize(function(delta)auc.true(delta)-theta.true), interval=c(0,4))$root
cluster.size <- 1

rho.statuses <- rho.markers <- c(0,.4,.8)
rhos <- expand.grid(rho.status=rho.statuses,rho.marker=rho.markers)

var.hats.normal <- lapply(1:nrow(rhos),function(i) {
    rho.status <- rhos[i,'rho.status']
    Sigma.status <- matrix(rho.status,nrow=cluster.size,ncol=cluster.size)
    diag(Sigma.status) <- 1

    rho.marker <- rhos[i,'rho.marker']
    Sigma.marker <- matrix(rho.marker,nrow=cluster.size,ncol=cluster.size)
    diag(Sigma.marker) <- 1

    coverage <- replicate(2e1, {
        status <-  rmvnorm(I,sigma=Sigma.status)
        status <- status>0

        data <- rmvnorm(I,sigma=Sigma.marker)
        data[status] <- data[status]+delta

        x <- lapply(1:I,function(i)data[i,!status[i,]])
        y <- lapply(1:I,function(i)data[i, status[i,]])
        est0 <- c(obu=auc.obu(x,y, alpha=.05)['var.hat'],jk=auc.jk(x,y,alpha=.05)['var.hat'])
    })
    coverage
})

var.obu.normal <- unlist(lapply(var.hats.normal,function(m)m['obu.var.hat',]))
var.jk.normal <- unlist(lapply(var.hats.normal,function(m)m['jk.var.hat',]))
plot(var.jk.normal ~ var.obu.normal)
max(abs(var.jk.normal-var.obu.normal))
abline(0,1)


var.hats.nonnormal <- lapply(1:nrow(rhos),function(i) {
    rho.status <- rhos[i,'rho.status']
    Sigma.status <- matrix(rho.status,nrow=cluster.size,ncol=cluster.size)
    diag(Sigma.status) <- 1

    rho.marker <- rhos[i,'rho.marker']
    Sigma.marker <- matrix(rho.marker,nrow=cluster.size,ncol=cluster.size)
    diag(Sigma.marker) <- 1

    coverage <- replicate(2e3, {
        status <-  rmvnorm(I,sigma=Sigma.status)
        status <- status>0

        data <- rmvnorm(I,sigma=Sigma.marker)
        data[status] <- data[status]+delta

        x <- lapply(1:I,function(i)data[i,!status[i,]])
        y <- lapply(1:I,function(i)data[i, status[i,]])

        x <- lapply(x,function(x.i)x.i[x.i < -.84] <- -.84)
        y <- lapply(y,function(y.i)y.i[y.i > .84] <- .84)
        est0 <- c(obu=auc.obu(x,y, alpha=.05)['var.hat'],jk=auc.jk(x,y,alpha=.05)['var.hat'])
    })
    coverage
})

var.obu.nonnormal <- unlist(lapply(var.hats.nonnormal,function(m)m['obu.var.hat',]))
var.jk.nonnormal <- unlist(lapply(var.hats.nonnormal,function(m)m['jk.var.hat',]))

## save.image('fig9.RData')
## png('fig9.png')
plot(c(var.obu.normal,var.obu.nonnormal),c(var.jk.normal,var.jk.nonnormal),xlab='Obuchowski estimator',ylab='jackknife variance estimator')
abline(0,1)
## dev.off()


## max difference versus I
Is <- seq(1e2,1e3,length.out=10)
set.seed(1)
max.diffs <- sapply(Is,function(I) {
    print(I)
    ## I <- 4e2
    theta.true <- .7
    delta <- uniroot(Vectorize(function(delta)auc.true(delta)-theta.true), interval=c(0,4))$root
    cluster.size <- 2

    rho.statuses <- rho.markers <- c(0,.4,.8)
    rhos <- expand.grid(rho.status=rho.statuses,rho.marker=rho.markers)

    var.hats.normal <- lapply(1:nrow(rhos),function(i) {
        rho.status <- rhos[i,'rho.status']
        Sigma.status <- matrix(rho.status,nrow=cluster.size,ncol=cluster.size)
        diag(Sigma.status) <- 1

        rho.marker <- rhos[i,'rho.marker']
        Sigma.marker <- matrix(rho.marker,nrow=cluster.size,ncol=cluster.size)
        diag(Sigma.marker) <- 1

        coverage <- replicate(3e1, {
            status <-  rmvnorm(I,sigma=Sigma.status)
            status <- status>0

            data <- rmvnorm(I,sigma=Sigma.marker)
            data[status] <- data[status]+delta

            x <- lapply(1:I,function(i)data[i,!status[i,]])
            y <- lapply(1:I,function(i)data[i, status[i,]])
            est0 <- c(obu=auc.obu(x,y, alpha=.05)['var.hat'],jk=auc.jk(x,y,alpha=.05)['var.hat'])
        })
        coverage
    })
    var.obu.normal <- unlist(lapply(var.hats.normal,function(m)m['obu.var.hat',]))
    var.jk.normal <- unlist(lapply(var.hats.normal,function(m)m['jk.var.hat',]))
    max(abs(var.jk.normal-var.obu.normal))
})
plot(Is,max.diffs,xlab='# of clusters',ylab='maximum absolute error')
lm0 <- lm(max.diffs ~ I(1/Is^2))
curve(coef(lm0)[1] +  coef(lm0)[2]/x^2,add=TRUE,col=disease.color)





M <- 40
vars <- replicate(1e2, {
    phi <- matrix(rbinom(M^2,1,.5),nrow=M,ncol=M)
    theta.hat <- sum(phi)/M^2
    Vx <- rowSums(phi)/M
    Vy <- colSums(phi)/M
    Sx <- sum((Vx-theta.hat)^2)
    Sy <- sum((Vy-theta.hat)^2)
    Sxy <- sum((Vx-theta.hat)*(Vy-theta.hat))
    var.obu <- (Sx+Sy+2*Sxy)/M/(M-1)
    ## sum((rowSums(phi)/M+colSums(phi)/M - 2*theta.hat)^2)/M/(M-1)
    ## ( sum((rowSums(phi)+colSums(phi)-diag(phi))^2)/M^2 - sum(diag(phi))/M^2 - 4*M*theta.hat^2 + 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi))) ) /M/(M-1)
    theta.del <- sapply(1:M,function(j)sum(phi[-j,-j])/(M-1)^2)
    pseudo <- M*theta.hat - (M-1)*theta.del
    var.jk <- var(pseudo)/M
    ## var.jk2 <- sum((pseudo - (M*theta.hat - sum(diag(phi))/M)/(M-1))^2)/M/(M-1)
    ## pseudo - mean(pseudo)
    ## M*theta.hat - (M-1)*theta.del - M*theta.hat/(M-1) + sum(diag(phi))/M/(M-1)
    ## M*(M-2)*theta.hat/(M-1)+sum(diag(phi))/M/(M-1) - (M^2*theta.hat - rowSums(phi)-colSums(phi)+diag(phi))/(M-1)
    ## sum(((rowSums(phi) + colSums(phi) - diag(phi) - 2*M*theta.hat + sum(diag(phi))/M ) / (M-1) )^2) /M/(M-1)
    ## ( sum((rowSums(phi)+colSums(phi)-diag(phi))^2)/(M-1)^2 - M*((sum(diag(phi))-2*M^2*theta.hat)/M/(M-1))^2 ) /M/(M-1)
    M*(M-1)*c(var.obu,var.jk)
})
plot(vars[1,],vars[2,]); abline(a=0,b=1)




M <- 40
vars <- replicate(4e2, {
    phi <- matrix(rbinom(M^2,1,.5),nrow=M,ncol=M)
    theta.hat <- sum(phi)/M^2
    c(
        sum((rowSums(phi)+colSums(phi)-diag(phi))^2)/M^2 - sum(diag(phi))/M^2 - 4*M*theta.hat^2 + 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi))) ,
        sum((rowSums(phi)+colSums(phi)-diag(phi))^2)/(M-1)^2 - M*((sum(diag(phi))-2*M^2*theta.hat)/M/(M-1))^2
      )
    ## c(
    ##      -sum(diag(phi))/M^2 - 4*M*theta.hat^2 + 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi))) ,
    ##      -M*((sum(diag(phi))-2*M^2*theta.hat)/M/(M-1))^2
    ## )
    ## c(
    ##     sum((rowSums(phi)+colSums(phi)-diag(phi))^2)/(M-1)^2,
    ##     - M*((sum(diag(phi))-2*M^2*theta.hat)/M/(M-1))^2
    ##     )
})
plot(vars[1,],vars[2,]); abline(a=0,b=1)


sum((colSums(phi)+rowSums(phi)-diag(phi))^2)
ones <- matrix(1,nrow=M,ncol=1)
t(ones)%*%(t(phi)%*%phi+2*phi%*%phi+phi%*%t(phi))%*%ones + sum(diag(phi)) - 2*t(ones)%*%(phi+t(phi))%*%diag(phi)


( t(ones)%*%(t(phi)%*%phi + phi%*%t(phi) + 2*phi%*%phi - 4/M*phi%*%ones%*%t(ones)%*%phi)%*%ones - 2*t(ones)%*%( phi - diag(M)/2 + t(phi) - 2/M*phi%*%ones%*%t(ones) + diag(phi)%*%t(ones)/2/M )%*%diag(phi) ) / (M-1)^2

sum((colSums(phi)+rowSums(phi)-diag(phi))^2)/(M-1)^2 - ((sum(diag(phi)) - 2*M^2*theta.hat)/(M-1))^2/M
(sum((colSums(phi)+rowSums(phi)-diag(phi))^2) - t(ones)%*%diag(phi)%*%t(ones)%*%diag(phi)/M - 4/M*t(ones)%*%phi%*%ones%*%t(ones)%*%phi%*%ones + 4/M*t(ones)%*%phi%*%ones%*%t(ones)%*%diag(phi))/(M-1)^2

sum((rowSums(phi)/M+colSums(phi)/M-2*theta.hat)^2)
1/M^2*t(ones)%*%( t(phi)%*%phi + phi%*%t(phi) + 2*phi%*%phi )%*%ones - 4*theta.hat^2*M
1/M^2*t(ones)%*%( t(phi)%*%phi + phi%*%t(phi) + 2*phi%*%phi -4/M*phi%*%ones%*%t(ones)%*%phi )%*%ones 


(2*M-1)/(M*(M-1))^2* t(ones)%*%(t(phi)%*%phi + phi%*%t(phi) + 2*phi%*%phi - 4/M*phi%*%ones%*%t(ones)%*%phi)%*%ones   -    2*t(ones)%*%( phi - diag(M)/2 + t(phi) - 2/M*phi%*%ones%*%t(ones) + diag(phi)%*%t(ones)/2/M )%*%diag(phi)/(M-1)^2 


(   (M-1/2)/(M)^2* t(ones)%*%(t(phi)%*%phi + phi%*%t(phi) + 2*phi%*%phi - 4/M*phi%*%ones%*%t(ones)%*%phi)%*%ones   -    t(ones)%*%( phi - diag(M)/2 + t(phi) - 2/M*phi%*%ones%*%t(ones) + diag(phi)%*%t(ones)/2/M )%*%diag(phi)  ) 


M <- 4e1
ones <- matrix(1,nrow=M,ncol=1)
vars <- replicate(1e2, {
    phi <- matrix(rbinom(M^2,1,.1),nrow=M,ncol=M)
    theta.hat <- sum(phi)/M^2
    c(
    (M-1/2)/M^2* t(ones)%*%(t(phi)%*%phi + phi%*%t(phi) + 2*phi%*%phi - 4/M*phi%*%ones%*%t(ones)%*%phi)%*%ones,
    t(ones)%*%( phi - diag(M)/2 + t(phi) - 2/M*phi%*%ones%*%t(ones) + diag(phi)%*%t(ones)/2/M )%*%diag(phi) 
      ) 
})
plot(vars[1,],vars[2,],asp=1); abline(a=0,b=1)

## 1
M <- 4e1
ones <- matrix(1,nrow=M,ncol=1)
vars <- replicate(1e2, {
    phi <- matrix(rbinom(M^2,1,runif(1)),nrow=M,ncol=M)
    phi <- matrix(rnorm(M^2),nrow=M,ncol=M)
    theta.hat <- sum(phi)/M^2
    Vx <- rowSums(phi)/M
    Vy <- colSums(phi)/M
    Sx <- sum((Vx-theta.hat)^2)
    Sy <- sum((Vy-theta.hat)^2)
    Sxy <- sum((Vx-theta.hat)*(Vy-theta.hat))
    var.obu <- (Sx+Sy+2*Sxy)/M/(M-1)
    theta.del <- sapply(1:M,function(j)sum(phi[-j,-j])/(M-1)^2)
    pseudo <- M*theta.hat - (M-1)*theta.del
    var.jk <- var(pseudo)/M
    c(
    (M-1/2)/M^2* t(ones)%*%(t(phi)%*%phi + phi%*%t(phi) + 2*phi%*%phi - 4/M*phi%*%ones%*%t(ones)%*%phi)%*%ones -
    t(ones)%*%( phi - diag(M)/2 + t(phi) - 2/M*phi%*%ones%*%t(ones) + diag(phi)%*%t(ones)/2/M )%*%diag(phi) ,
    1/2*(M-1)^3*M*(var.jk-var.obu)
      ) #*2/(M-1)^2/M/(M-1)
})
plot(vars[1,],vars[2,],asp=1); abline(a=0,b=1)
max(abs(vars[1,]-vars[2,]))



M <- 1e1
ones <- matrix(1,nrow=M,ncol=1)
vars <- replicate(1e2, {
    phi <- matrix(rbinom(M^2,1,.5),nrow=M,ncol=M)
    phi <- matrix(0,nrow=M,ncol=M)
    phi[row(phi)<col(phi)] <- rbinom(M*(M-1)/2,1,.5)
    phi <- phi + t(phi)
    diag(phi) <- rbinom(M,1,.5)
    theta.hat <- sum(phi)/M^2
    c(
    (M-1/2)/M^2* t(ones)%*%(t(phi)%*%phi + phi%*%t(phi) + 2*phi%*%phi - 4/M*phi%*%ones%*%t(ones)%*%phi)%*%ones -
    t(ones)%*%( phi - diag(M)/2 + t(phi) - 2/M*phi%*%ones%*%t(ones) + diag(phi)%*%t(ones)/2/M )%*%diag(phi) 
      ) #*2/(M-1)^2/M/(M-1)
})
hist(vars);abline(v=0,col='red')


## 2 try all matrices for small M to find actual max error
M <- 5
combn <- as.matrix(expand.grid(rep(list(0:1),M^2)))
ones <- matrix(1,nrow=M,ncol=1)
errs <- sapply(1:nrow(combn),function(j) {
    phi <- matrix(combn[j,],nrow=M,ncol=M)
    theta.hat <- sum(phi)/M^2
    c(
    (M-1/2)/M^2* t(ones)%*%(t(phi)%*%phi + phi%*%t(phi) + 2*phi%*%phi - 4/M*phi%*%ones%*%t(ones)%*%phi)%*%ones,
    t(ones)%*%( phi - diag(M)/2 + t(phi) - 2/M*phi%*%ones%*%t(ones) + diag(phi)%*%t(ones)/2/M )%*%diag(phi) 
      ) #*2/(M-1)^2/M/(M-1)
})
op <- par(mfrow=c(1,2))
plot(errs[1,],errs[2,])
hist(diff(errs));abline(v=0,col='red')
par(op)
max(abs(diff(errs)))
## save.image('2.RData')

## max error by p
M <- 4e1
ps <- seq(0,1,length.out=20)
ones <- matrix(1,nrow=M,ncol=1)
by.p <- sapply(ps, function(p) {
    diffs <- replicate(1e2, {
        phi <- matrix(rbinom(M^2,1,p),nrow=M,ncol=M)
        theta.hat <- sum(phi)/M^2
        c(
        (M-1/2)/M^2* t(ones)%*%(t(phi)%*%phi + phi%*%t(phi) + 2*phi%*%phi - 4/M*phi%*%ones%*%t(ones)%*%phi)%*%ones -
        t(ones)%*%( phi - diag(M)/2 + t(phi) - 2/M*phi%*%ones%*%t(ones) + diag(phi)%*%t(ones)/2/M )%*%diag(phi) 
        ) #*2/(M-1)^2/M/(M-1)
    })
    diff(range(diffs))
})
plot(ps,by.p)



M <- 1e1
ones <- matrix(1,nrow=M,ncol=1)
vars <- replicate(1e2, {
    phi <- matrix(rbinom(M^2,1,.5),nrow=M,ncol=M)
    phi <- matrix(0,nrow=M,ncol=M)
    phi[row(phi)<col(phi)] <- rbinom(M*(M-1)/2,1,.5)
    phi <- phi + t(phi)
    diag(phi) <- rbinom(M,1,.5)
    theta.hat <- sum(phi)/M^2
    c(
    (M-1/2)/M^2* t(ones)%*%(t(phi)%*%phi + phi%*%t(phi) + 2*phi%*%phi - 4/M*phi%*%ones%*%t(ones)%*%phi)%*%ones -
    t(ones)%*%( phi - diag(M)/2 + t(phi) - 2/M*phi%*%ones%*%t(ones) + diag(phi)%*%t(ones)/2/M )%*%diag(phi) 
      ) #*2/(M-1)^2/M/(M-1)
})



ones <- matrix(1,nrow=M,ncol=1)
phi <- matrix(1,nrow=M,ncol=M)
diag(phi) <- 0
t(ones)%*%(t(phi)%*%phi + phi%*%t(phi) + 2*phi%*%phi - 4/M*phi%*%ones%*%t(ones)%*%phi)%*%ones 

hist(vars);abline(v=0,col='red')

M <- 4e1
ones <- matrix(1,nrow=M,ncol=1)
A <- diag(M) - ones%*%t(ones)/M
phi <- matrix(1,nrow=M,ncol=M)
t(ones)%*%phi%*%A%*%phi%*%ones


A <- diag(M) - ones%*%t(ones)/M
phi <- matrix(0,nrow=M,ncol=M)
phi[row(phi)<col(phi)] <- rbinom(M*(M-1)/2,1,.5)
phi <- phi + t(phi)
diag(phi) <- rbinom(M,1,.5)
## t(ones)%*%(t(phi)%*%phi + phi%*%t(phi) + 2*phi%*%phi - 4/M*phi%*%ones%*%t(ones)%*%phi)%*%ones 
phi <- matrix(rnorm(M^2),nrow=M,ncol=M); phi <- phi + t(phi)
s <- rowSums(phi)
t(ones)%*%phi%*%A%*%phi%*%ones
sum(s^2)-sum(s)^2/M
var(s)*(M-1)




M <- 1e2
ones <- matrix(1,nrow=M,ncol=1)
vars <- replicate(1e2, {
    phi <- matrix(rbinom(M^2,1,.5),nrow=M,ncol=M)
    phi <- matrix(runif(M^2,0,M),nrow=M,ncol=M)
    c(
        t(ones)%*%(phi%*%t(phi))%*%ones,
        t(ones)%*%( phi%*%ones%*%t(ones)%*%phi/M)%*%ones
      ) 
})
plot(vars[1,],vars[2,]);abline(a=0,b=1,col='red')

Ms <- floor(seq(1e1,5e2,length.out=1e1))
by.M <- sapply(Ms, function(M) {
## M <- 1e1
ones <- matrix(1,nrow=M,ncol=1)
vars <- replicate(2e3, {
    phi <- matrix(rbinom(M^2,1,.5),nrow=M,ncol=M)
    ## phi <- matrix(runif(M^2,0,M),nrow=M,ncol=M)
    ## colSums(phi)%*%(rowSums(phi) - mean(rowSums(phi)))
    ## colSums(phi)%*%rowSums(phi) - sum(phi)^2/M
    M*cov(colSums(phi),rowSums(phi)) 
})
max(abs(vars))
})
plot(Ms,by.M)
lines(Ms,Ms*sqrt(Ms))
## max of colSums %*% (rowSums - mean(rowSums)) is Msqrt(M) for binary
## phi when phi has bernoulli entries

by.M <- sapply(Ms,function(M){
phi <- matrix(0,nrow=M,ncol=M)
phi[row(phi)<=M/2] <- 1
phi[col(phi)<=M/2] <- 1
M*cov(rowSums(phi),colSums(phi))
})
plot(Ms,by.M)
lines(Ms,Ms^2.5,col='red')

Ms <- floor(seq(1e1,5e2,length.out=1e1))
by.M <- sapply(Ms, function(M) {
## M <- 1e1
ones <- matrix(1,nrow=M,ncol=1)
vars <- replicate(2e3, {
    ## phi <- matrix(rbinom(M^2,1,.5),nrow=M,ncol=M)
    ## phi <- matrix(runif(M^2,0,M),nrow=M,ncol=M)
    c <- rbinom(M,1,.5)
    ## rbinom(M,1,.5)%*%(c-mean(c))
    M*cov(rbinom(M,1,.5),c)
})
max(abs(vars))
})
plot(Ms,by.M)
lines(Ms,sqrt(Ms))
## max of phi[1,] %*% (phi[,1] - mean(phi[,1])) is sqrt(M) for binary
## phi




## 1 special upper triangular-like structure of phi
M <- sample(10:200,1)
## M <- 10
x <- rnorm(M); y <- rnorm(M)
x <- sort(x); y <- sort(y)
sum(outer(x,y,`<`))
phi <- t(sapply(x,function(x.i)as.numeric(x.i<y)))
sum(phi)
p <- sapply(x,function(x.i)which.max(x.i<c(y,Inf)))
p[x>max(y)] <- M+1
stopifnot(sum(p)==M*(M+1)-sum(phi))
sapply(1:M,function(j)sum(p<=j))
colSums(phi)
sum((colSums(phi)+rowSums(phi))*diag(phi))
sum(diag(phi))*M - sum(p*diag(phi))+ sum(diag(phi)) + sum(diag(phi)*(1:M-1)) + sum(diag(phi)*(1:M+1-p))
(M+1)*sum(diag(phi)) - sum(p*diag(phi)) + sum(diag(phi)*(2*(1:M)-p))
(M+1)*sum(diag(phi)) + sum(diag(phi)*(2*(1:M)-p) - p*diag(phi))
(M+1)*sum(diag(phi)) + 2*sum(diag(phi)*(1:M-p))
sum(diag(phi)*(1:M - p))
sum(phi[row(phi)>col(phi)])

M <- 10
x <- rnorm(M); y <- rnorm(M)
x <- sort(x); y <- sort(y)
phi <- t(sapply(x,function(x.i)as.numeric(x.i<y)))
phi


M <- sample(10:200,1)
p <- sort(sample(1:M,replace=TRUE))
phi <- t(sapply(p,function(j)c(rep(0,max(0,j-1)),rep(1,M-j+1))))
tr <- sum(diag(phi)); theta.hat <- sum(phi)/M^2
sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/(M-1)^2 - 1/M^2)
tr/M^2 + 4*M*theta.hat^2 - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))
1/M * ((tr-2*M^2*theta.hat)/(M-1))^2
4*M*theta.hat^2 + 1/M^2 * sum(diag(phi) - 4*diag(phi)*(1:M - p))
4*theta.hat^2*M*(1-(M/(M-1))^2)

## two parts inversely correlated though the first only depends on the sum of the p's
## all parts highly correlated
ans <- replicate(1e4,{
M <- sample(10:200,1)
p <- sort(sample(1:M,replace=TRUE))
phi <- t(sapply(p,function(j)c(rep(0,max(0,j-1)),rep(1,M-j+1))))
tr <- sum(diag(phi)); theta.hat <- sum(phi)/M^2
4*theta.hat^2*M*(1-(M/(M-1))^2) - 2/M^2*(2*M+1)*tr - 4/M^2*sum(diag(phi)*(1:M - p)) - 1/M*(tr/(M-1))^2 + 4*M*theta.hat*tr/(M-1)^2
## c(4*theta.hat^2*M*(1-(M/(M-1))^2),2/M^2*(2*M+1)*tr - 4/M^2*sum(diag(phi)*(1:M - p)) - 1/M*(tr/(M-1))^2 + 4*M*theta.hat*tr/(M-1)^2)
c(4*theta.hat^2*M*(1-(M/(M-1))^2),2/M^2*(2*M+1)*tr, - 4/M^2*sum(diag(phi)*(1:M - p)), - 1/M*(tr/(M-1))^2, + 4*M*theta.hat*tr/(M-1)^2)
})
cor(t(ans))
rowMeans(abs(ans))


## verify structures/formulas for largest and smallest upper sum for fixed trace
upper <- function(u) {
    d <-  length(u)
    ## print(d)
    if(d>=M) {
        ## browser()
        us <<- rbind(us,u)
        return()
        ## print('here')
    }
    for(m in max(d+1,u[d]):(M+1))
        upper(c(u,m))
}
u.to.m <- function(u) t(sapply(u,function(u.i)c(rep(0,u.i-1),rep(1,M-u.i+1))))
M <- 10
us <- numeric()
for(u1 in 1:(M+1))upper(u1)
tr <- colSums(t(us) == 1:M)
upper.sums <- rowSums(M+1 - us)
t = 1
u.to.m(rbind(us[tr==t,])[which.max(upper.sums[tr==t]),])
u.to.m(rbind(us[tr==t,])[which.min(upper.sums[tr==t]),])
## formulas for bounds on upper sum in terms of trace
max(upper.sums[tr==t])
1/2*M*(M+1) - (M-t)
min(upper.sums[tr==t])
1/2*M*(M+1) - 1/2*(M-t)*(M-t+1)

## 1a upper and lower
to.m <- function(pivots) t(sapply(pivots,function(pivot)c(rep(0,pivot-1),rep(1,M-pivot+1))))
M <- 5
pivots <- do.call(expand.grid,rep(list(1:(M+1)),M))
diffs <- t(apply(pivots,1,diff))
pivots <- pivots[rowSums(diffs<0)==0,]
tr <- colSums((t(pivots) - (1:M))<=0)
    
upper <- t(apply(pivots,1,function(r)pmax(r,1:M)))
upper.sums <- rowSums(M+1 - upper)
t <- sample(1:M,1)
## to.m(rbind(upper[tr==t,])[which.max(upper.sums[tr==t]),])
stopifnot(max(upper.sums[tr==t]) == 1/2*M*(M+1) - (M-t))
## to.m(rbind(upper[tr==t,])[which.min(upper.sums[tr==t]),])
stopifnot(min(upper.sums[tr==t]) == 1/2*M*(M+1) - 1/2*(M-t)*(M-t+1))
lower.sums <- colSums(pmax((1:M) - t(pivots),0))
stopifnot(max(lower.sums[tr==t]) == 1/2*t*(t-1))
stopifnot(min(lower.sums[tr==t]) == 0)

## 2  checking formulas
M <- sample(10:200,1)
M <- 7
p <- sort(sample(1:M,replace=TRUE))
## phi <- t(sapply(p,function(j)c(rep(0,max(0,j-1)),rep(1,M-j+1))))
## tr <- sum(diag(phi)); theta.hat <- sum(phi)/M^2
tr <- sample(1:M,1)
phi.max <- t(sapply(1:M,function(j)c(rep(1*(j<=tr),j-1 + (j>tr)),rep(1,M-j+1 - (j>tr)))))
theta.hat <- sum(phi.max)/M^2
theta.hat
1/2*(1-1/M+tr*(tr+1)/M^2) 

##2b
theta.hat^2*(1-(M/(M-1))^2)
1/4*(-2/M+1/M^2 - 2*tr*(tr+1)*(2*M^1-1)/M^3/(M-1) + (tr*(tr+1))^2*(1-2*M)/M^4/(M-1)^2)

##2c
4*M/(M-1)^2*theta.hat*tr
2*M/(M-1)^2*(tr-tr/M+tr^2*(tr+1)/M^2) 

##2d
4*M*theta.hat^2*(1-(M/(M-1))^2) + 4*M/(M-1)^2*theta.hat*tr - (2*M+1)*tr/M^2 - 1/M*(tr/(M-1))^2
-2 + 1/M - 2*tr*(tr+1)*(2*M-1)/M^2/(M-1) + (tr*(tr+1))^2/M^3/(M-1)^2*(1-2*M) + 2*M*tr/(M-1)^2 - 2*tr/(M-1)^2 + 2*tr^2*(tr+1)/M/(M-1)^2 - (2*M+1)/M^2*tr - tr^2/M/(M-1)^2
## (-1 + tr/M - tr^2/M^2 + tr^3/M^3 - tr^4/M^4)*2

M <- sample(10:200,1)
p <- sort(sample(1:M,replace=TRUE))
phi <- t(sapply(p,function(j)c(rep(0,max(0,j-1)),rep(1,M-j+1))))
tr <- sum(diag(phi)); theta.hat <- sum(phi)/M^2
tr/M^2+4*M*theta.hat^2-2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) - 1/M*((tr-2*M^2*theta.hat)/(M-1))^2 ##2e
1/M*((tr-2*M^2*theta.hat)/(M-1))^2
(tr/(M-1))^2/M - 4*M*theta.hat*tr/(M-1)^2 + 4*M^3*theta.hat^2/(M-1)^2 ##2f
tr/M^2+4*M*theta.hat^2-2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))
4*M*theta.hat^2-1/M^2*(sum(diag(phi)*(4*(1:M)-4*p+2*M+1))) ## 2g
4*M*theta.hat^2-1/M^2*(sum(diag(phi)*(4*(1:M)-4*p+2*M+1))) - ((tr/(M-1))^2/M - 4*M*theta.hat*tr/(M-1)^2 + 4*M^3*theta.hat^2/(M-1)^2) ## 2h==2e


## 2i first order approximtion to upper bound
M <- sample(10:200,1)
M <- 7
## p <- sort(sample(1:M,replace=TRUE))
## phi <- t(sapply(p,function(j)c(rep(0,max(0,j-1)),rep(1,M-j+1))))
## tr <- sum(diag(phi)); theta.hat <- sum(phi)/M^2
tr <- sample(1:M,1)
phi.max <- t(sapply(1:M,function(j)c(rep(1*(j<=tr),j-1 + (j>tr)),rep(1,M-j+1 - (j>tr)))))
theta.hat <- sum(phi.max)/M^2
p <- apply(cbind(phi.max,1),1,which.max)

tr/M^2 +4*M*theta.hat^2-2/M^2*sum(diag(phi.max)*(colSums(phi.max)+rowSums(phi.max))) - 1/M*((tr-2*M^2*theta.hat)/(M-1))^2
4*M*theta.hat^2-1/M^2*(sum(diag(phi.max)*(4*(1:M)-4*p+2*M+1)))   -    ((tr/(M-1))^2/M - 4*M*theta.hat*tr/(M-1)^2 + 4*M^3*theta.hat^2/(M-1)^2 )

4*M*theta.hat^2*(1-(M/(M-1))^2) + 4*M/(M-1)^2*theta.hat*tr - (2*M+1)*tr/M^2 - 1/M*(tr/(M-1))^2
-2 + 1/M - 2*tr*(tr+1)*(2*M-1)/M^2/(M-1) + (tr*(tr+1))^2/M^3/(M-1)^2*(1-2*M) + 2*M*tr/(M-1)^2 - 2*tr/(M-1)^2 + 2*tr^2*(tr+1)/M/(M-1)^2 - (2*M+1)/M^2*tr - tr^2/M/(M-1)^2

4/M^2*sum(diag(phi.max)*(1:M - p)) ## ignored term

4*M*theta.hat^2*(1-(M/(M-1))^2) + 4*M/(M-1)^2*theta.hat*tr - (2*M+1)*tr/M^2 - 1/M*(tr/(M-1))^2
-2 + 1/M - 2*tr*(tr+1)*(2*M-1)/M^2/(M-1) + (tr*(tr+1))^2/M^3/(M-1)^2*(1-2*M) + 2*M*tr/(M-1)^2 - 2*tr/(M-1)^2 + 2*tr^2*(tr+1)/M/(M-1)^2 - (2*M+1)/M^2*tr - tr^2/M/(M-1)^2
-2 - 4*tr^2/M^2-2*tr^4/M^4+2*tr^3/M^3
dd

## *M*theta.hat^2*(1-(M/(M-1))^2) + 4*M/(M-1)^2*theta.hat*tr - 2*(2*M+1)*tr/M^2 - 1/M*(tr/(M-1))^2
## ## (-1 + tr/M - tr^2/M^2 + tr^3/M^3 - tr^4/M^4)*2

## -2 - tr*(tr+1)*(2*M-1)/M^2/(M-1) + 1/M  + (tr*(tr+1))^2/M^3/(M-1)^2*(1-2*M) + 2*M*tr/(M-1)^2 - 2*tr/(M-1)^2 + 2*tr^2*(tr+1)/M/(M-1)^2 - 2*(2*M+1)/M^2*tr - tr^2/M/(M-1)^2
## -2 - tr^2*2/M/(M-1) + 2*tr^4/M^2/(M-1)^2 + 2*M*tr/(M-1)^2 + 2*tr^3/M/(M-1)^2 - 4*M*tr/M^2
## -1 + M*tr/(M-1)^2 - 2*M*tr/M^2 - tr^2/M/(M-1) + tr^4/M^2/(M-1)^2  + tr^3/M/(M-1)^2 
## -1 - tr/M + tr^2/M^2 - tr^3/M^3 + tr^4/M^4

## checking relative magnitude of the omitted term (lower bound sum) from upper bound  approximation
summands <- replicate(5e3,{
M <- sample(10:200,1)
M <- 10
p <- sort(sample(1:(M+1),size=M,replace=TRUE))
phi <- t(sapply(p,function(j)c(rep(0,max(0,j-1)),rep(1,M-j+1))))
tr <- sum(diag(phi)); theta.hat <- sum(phi)/M^2
## l <- sum(phi[row(phi)>col(phi)])
## u <- sum(phi[row(phi)<=col(phi)])
## 4*theta.hat^2*M*(1-(M/(M-1))^2) + 4*theta.hat*M*tr/(M-1)^2 - 4/M^2*sum(diag(phi)*(1:M - p))
## 4/M^3*(1-(M/(M-1))^2)*(u^2+2*u*l+l^2) + 4*tr/M/(M-1)^2*(u+l) - 4/M^2*l
## c(4*theta.hat^2*M*(1-(M/(M-1))^2) , 4*theta.hat*M*tr/(M-1)^2, - 4/M^2*sum(diag(phi)*(1:M - p)))
c(-2 + 1/M - 2*tr*(tr+1)*(2*M-1)/M^2/(M-1) + (tr*(tr+1))^2/M^3/(M-1)^2*(1-2*M) + 2*M*tr/(M-1)^2 - 2*tr/(M-1)^2 + 2*tr^2*(tr+1)/M/(M-1)^2 - (2*M+1)/M^2*tr - tr^2/M/(M-1)^2, 4/M^2*sum(diag(phi)*(1:M - p)))
})
rowMeans(abs(summands))




## preserving paired structure

## 3. sizes of summands
M <- 10
x <- rnorm(M); y <- rnorm(M)
x <- x[order(y)]; y <- y[order(y)]
phi <- t(sapply(x,function(x.i)x.i<y))+0
phi

M <- 5
ps <- sample(1:(M+1),M)
phi <- t(sapply(ps, function(p) c(rep(0,p-1),rep(1,M+1-p))))
t <- sum(diag(phi)); theta.hat <- sum(phi)/M^2
ps.to.m <- function(ps) t(sapply(ps, function(p) c(rep(0,p-1),rep(1,M+1-p))))
with(list(phi=phi[c(1,2,4,3,5),]),sum(diag(phi)*(colSums(phi)+rowSums(phi))))
with(list(phi=phi),sum(diag(phi)*(colSums(phi)+rowSums(phi))))

sum(diag(phi)*(colSums(phi)+rowSums(phi)))
phi

with(list(phi=phi[c(1,2,4,3,5),]),sum(diag(phi)*(colSums(phi))))
with(list(phi=phi),sum(diag(phi)*(colSums(phi))))


M <- 500
parts <- replicate(1e3,{
ps <- sample(1:(M+1),M,replace=TRUE)
phi <- t(sapply(ps, function(p) c(rep(0,p-1),rep(1,M+1-p))))
t <- sum(diag(phi)); theta.hat <- sum(phi)/M^2
## 4*M*theta.hat^2*(1-(M/(M-1))^2) + 4*M/(M-1)^2*theta.hat*t - 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi))) - (2*M+1)*t/M^2 - 1/M*(t/(M-1))^2
c(4*M*theta.hat^2*(1-(M/(M-1))^2),4*M/(M-1)^2*theta.hat*t, - 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi))) )
})
rowMeans(abs(parts))
cor(t(parts))


## 4. structure of minimizing and maximizing phis
ps.to.m <- function(ps) t(sapply(ps, function(p) c(rep(0,p-1),rep(1,M+1-p))))
M <- sample(10:200,1)
M <- 5
ps <- sample(1:(M+1),M,replace=TRUE)
## phi <- t(sapply(ps, function(p) c(rep(0,p-1),rep(1,M+1-p))))
## theta.hat <- sum(phi)/M^2
t <- sum(diag(phi))
S <- sum(M+1 - ps)
## stopifnot(sum(phi)==S)
## 4*M*theta.hat^2*(1-(M/(M-1))^2) + 4*M/(M-1)^2*theta.hat*t - 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi)))
4*S^2/M^3*(1-(M/(M-1))^2) + 4*S*t/M/(M-1)^2 - 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi)))

M <- 6
ps.list <- expand.grid(rep(list(1:(M+1)),M))
diffs <- apply(ps.list,1,function(ps) {
    t <- sum(ps<=(1:M))
    S <- sum(M+1 - ps)
    phi <- ps.to.m(ps)
    c(t,4*S^2/M^3*(1-(M/(M-1))^2) + 4*S*t/M/(M-1)^2 - 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi))))
})
by.t <- split(cbind(diffs=diffs[2,],ps.list),diffs[1,])
max.phis <- lapply(by.t,function(ps){
    ## ps.to.m(ps[which.max(ps$diffs),-1])
    lapply(split(ps[ps$diffs==max(ps$diffs),-1],1:sum(ps$diffs==max(ps$diffs))),ps.to.m)
})
min.phis <- lapply(by.t,function(ps){
    lapply(split(ps[ps$diffs==min(ps$diffs),-1],1:sum(ps$diffs==min(ps$diffs))),ps.to.m)
})
with(list(phi=min.phis[[7]][[1]]),{t <- sum(diag(phi)); S <- sum(phi); c(4*S^2/M^3*(1-(M/(M-1))^2), + 4*S*t/M/(M-1)^2, - 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi))))})
with(list(phi=ps.to.m(c(7,7,7,7,7,6))),{t <- sum(diag(phi)); S <- sum(phi); 4*S^2/M^3*(1-(M/(M-1))^2) + 4*S*t/M/(M-1)^2 - 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi)))})

## 4a
unlist(lapply(max.phis,function(l)lapply(l,sum)))
sapply(0:M,function(t)t*(t+1)/2)
unlist(lapply(max.phis,function(l)lapply(l,function(phi){sum(diag(phi)*(colSums(phi)+rowSums(phi)))})))
sapply(0:M,function(t)t*(t+1))
unlist(lapply(min.phis,function(l)lapply(l,sum)))
sapply(0:M,function(t)M^2-1/2*(M-t)*(M-t+1))
unlist(lapply(min.phis,function(l)lapply(l,function(phi){sum(diag(phi)*(colSums(phi)+rowSums(phi)))})))
sapply(0:M,function(t)2*M*t)

## 4b lower bound on first summands
with(list(phi=min.phis[[3]][[1]]),{
    t <- sum(diag(phi))
    S <- sum(phi)
    c(4*S^2/M^3*(1-(M/(M-1))^2)+ 4*S*t/M/(M-1)^2 - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))),
      4/M^3*(1-(M/(M-1))^2)*(M^2/2-t^2/2-M/2+t/2+M*t+t*M^2/2/(1-2*M))^2- t^2*M/(M-1)^2/(1-2*M)- 4*t/M
      )
    })

## 4c upper bound
with(list(phi=max.phis[[3]][[1]]),{
    t <- sum(diag(phi))
    S <- sum(phi)
    c(4*S^2/M^3*(1-(M/(M-1))^2)+ 4*S*t/M/(M-1)^2 - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))),
      4/M^3*(1-(M/(M-1))^2)*(t*(t+1)/2 + t*M^2/2/(1-2*M))^2 - t^2*M/(M-1)^2/(1-2*M)- 2*t*(t+1)/M^2
      )
    })

## 4d
deltas <- replicate(1e3, {
M <- sample(10:200,1)
ps <- sample(1:(M+1),M,replace=TRUE)
phi <- t(sapply(ps, function(p) c(rep(0,p-1),rep(1,M+1-p))))
theta.hat <- sum(phi)/M^2; t <- sum(diag(phi)); S <- sum(M+1 - ps)
terms <- 4*S^2/M^3*(1-(M/(M-1))^2)+ 4*S*t/M/(M-1)^2 - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))
stopifnot(terms <= 4/M^3*(1-(M/(M-1))^2)*(t*(t+1)/2 + t*M^2/2/(1-2*M))^2 - t^2*M/(M-1)^2/(1-2*M)- 2*t*(t+1)/M^2)
stopifnot(terms >=       4/M^3*(1-(M/(M-1))^2)*(M^2/2-t^2/2-M/2+t/2+M*t+t*M^2/2/(1-2*M))^2- t^2*M/(M-1)^2/(1-2*M)- 4*t/M)
terms - (4/M^3*(1-(M/(M-1))^2)*(M^2/2-t^2/2-M/2+t/2+M*t+t*M^2/2/(1-2*M))^2- t^2*M/(M-1)^2/(1-2*M)- 4*t/M)
})
hist(deltas)


## 4e lower bound on first summands -- O(1) terms
M <- 100
t <- sample(0:M,1)
phi <- t(sapply(1:M, function(m)c(rep(0+(M-m<t),m),rep(1,M-m))))
S <- sum(phi)
c(4*S^2/M^3*(1-(M/(M-1))^2)+ 4*S*t/M/(M-1)^2 - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))),
  4/M^3*(1-(M/(M-1))^2)*(M^2/2-t^2/2-M/2+t/2+M*t+t*M^2/2/(1-2*M))^2- t^2*M/(M-1)^2/(1-2*M)- 4*t/M,
  1/M^3*(1-2*M)/(M-1)^2*(M^2-t^2+2*M*t+t*M^2/(1-2*M))^2 - t^2*M/(M-1)^2/(1-2*M) - 4*t/M
  )
r <- t/M
-2*(r^4-3*r^3+5*r+1) ## 4f


## 4f upper bound on first summands -- O(1) terms
M <- 100
t <- sample(0:M,1)
phi <- t(sapply(1:M, function(m)c(rep(0,m-1),rep(0+(M-m<t),M-m+1))))
S <- sum(phi)
theta.hat <- S/M^2
4/M^3*(1-(M/(M-1))^2)*(t*(t+1)/2 + t*M^2/2/(1-2*M))^2 - t^2*M/(M-1)^2/(1-2*M)- 2*t*(t+1)/M^2
(1-2*M)/M^5*(t^2-t*M/2)^2 + t^2/2/M^2 - 2*t^2/M^2
-2*(t^4/M^4-t^3/M^3+t^2/M^2)
## 4*M*theta.hat^2*(1-(M/(M-1))^2) + 4*M/(M-1)^2*theta.hat*t - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) - (2*M+1)*t/M^2 - 1/M*(t/(M-1))^2 
## 4*M*theta.hat^2*(1-(M/(M-1))^2) + 4*M/(M-1)^2*theta.hat*t - 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi)))
## 4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2 - t^2*M/(M-1)^2/(1-2*M) #- 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi)))

-1/M*((t-2*M^2*theta.hat)/(M-1))^2 + t/M^2 + 4*M*theta.hat^2 - 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi)))
4*M*theta.hat^2*(1-(M/(M-1))^2) + 4*M/(M-1)^2*theta.hat*t - 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi))) - t^2/M/(M-1)^2 + t/M^2
4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2 - t^2*M/(M-1)^2/(1-2*M) - 2/M^2*sum(diag(phi)*(rowSums(phi)+colSums(phi))) - t^2/M/(M-1)^2 + t/M^2


concord.xy <- function(x,y) {
    m <- length(x)
    n <- length(y)
    ## stopifnot(m>1,n>1)
    if(length(x)==0 | length(y)==0)return(0)

    sum(sapply(y,function(y.i)x<y.i))
    ## sum(outer(x,y,`<`))
}

auc.obu <- function(x,y,alpha=.05) {

    I <- length(x)
    stopifnot(length(x)==length(y))
    m <- sapply(x,length)
    n <- sapply(y,length)
    I.10 <- sum(m>0)
    I.01 <- sum(n>0)
    M <- sum(m)
    N <- sum(n)

    theta.hat <- concord.xy(unlist(x),unlist(y)) / (M*N)
    V.10 <- sapply(x,function(x.i)concord.xy(x.i,unlist(y)))/N
    V.01 <- sapply(y,function(y.i)concord.xy(unlist(x),y.i))/M

    S.10 <- sum((V.10 - m*theta.hat)^2) * I.10 / ((I.10-1)*M)
    S.01 <- sum((V.01 - n*theta.hat)^2) * I.01 / ((I.01-1)*N)
    S.11 <- sum((V.10 - m*theta.hat)*(V.01 - n*theta.hat)) * I / (I-1)

    var.hat <- S.10/M + S.01/N + 2*S.11/(M*N)
    q <- qnorm(1-alpha/2)
    return(c(theta.hat=theta.hat,var.hat=var.hat,CI.lower=theta.hat-q*sqrt(var.hat),CI.upper=theta.hat+q*sqrt(var.hat)))

}


auc.jk <- function(x,y,alpha=.05) {
    concord <- function(x,y) mean(outer(x,y,`<`))
    concord.sum <- function(x,y)sum(outer(x,y,`<`))
    m <- sapply(x,length); n <- sapply(y,length)
    M <- sum(m); N <- sum(n)
    I <- length(x); stopifnot(length(x)==length(y))

    theta.hat <- concord(unlist(x),unlist(y))
    theta.del <- sapply(1:I,function(i)concord(unlist(x[-i]),unlist(y[-i])))
    pseudo <- I*theta.hat - (I-1)*theta.del
    theta.jk <- mean(pseudo)
    var.hat <- var(pseudo)/I
    q <- qnorm(1-alpha/2)


    return(c(theta.hat=theta.jk,var.hat=var.hat,CI.lower=theta.jk-q*sqrt(var.hat),CI.upper=theta.jk+q*sqrt(var.hat)))
}

## verify obu/jk variance formulas
M <- 100
x <- rnorm(M); y <- rnorm(M)
phi <- t(sapply(x,function(x.i)as.numeric(x.i<y)))
theta.hat <- sum(phi)/M^2; t <- sum(diag(phi))
S <- sum(phi)
auc.obu(x,y)['var.hat']*M*(M-1)
1/M^2*sum((colSums(phi)+rowSums(phi)-diag(phi))^2) - t/M^2 - 4*M*theta.hat^2 + 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))
auc.jk(x,y)['var.hat']*M*(M-1)
1/(M-1)^2*sum((colSums(phi)+rowSums(phi)-diag(phi))^2) - M*((t-2*M^2*theta.hat)/M/(M-1))^2

auc.jk(x,y)['var.hat']*M*(M-1) - auc.obu(x,y)['var.hat']*M*(M-1)
- t/M^2 - 4*M*theta.hat^2 + 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) + M*((t-2*M^2*theta.hat)/M/(M-1))^2 + sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2)

## 5 full objective
 t/M^2 + 4*M*theta.hat^2 - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) - M*((t-2*M^2*theta.hat)/M/(M-1))^2 
-    sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2)

4*M*theta.hat^2*(1-(M/(M-1))^2) + 4*M*theta.hat*t/(M-1)^2 - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) - t^2/M/(M-1)^2 + t/M^2
-    sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2)

4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2 - t^2*M/(M-1)^2/(1-2*M)  - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) - t^2/M/(M-1)^2 + t/M^2
-    sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2) # 5a

4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2 - t^2*M/(M-1)^2/(1-2*M)  - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) - t^2/M/(M-1)^2 + t/M^2 -  sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2) 

## O(1) terms?
4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2  - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) - (1/M^2 - 1/(M-1)^2) *sum((colSums(phi)+rowSums(phi)-diag(phi))^2)
-8/M^4*(S+t*M^2/2/(1-2*M))^2  - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) + 2/M^3*sum((colSums(phi)+rowSums(phi)-diag(phi))^2)

4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2   - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))# - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2
-    sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2)
- 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) -    sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2)
-2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))*((2*M-1)/(M*(M-1))^2 + 1/M^2) + (2*M-1)*t/(M*(M-1))^2 + (2*M-1)/(M*(M-1))^2*sum((colSums(phi)+rowSums(phi))^2) # 5b
4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2   - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) -  sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2) + (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2

## check size of cross terms
M <- 100
cross <- replicate(1e2,{
x <- rnorm(M); y <- rnorm(M)
phi <- t(sapply(x,function(x.i)as.numeric(x.i<y)))
theta.hat <- sum(phi)/M^2; t <- sum(diag(phi))
S <- sum(phi)
## (2*S^2 - sum((colSums(phi)+rowSums(phi))^2))/(2*S^2)
## c(4*S^2 - sum((colSums(phi)+rowSums(phi))^2),sum((colSums(phi)+rowSums(phi))^2) / (4*S^2 - sum((colSums(phi)+rowSums(phi))^2)))
c(sum((colSums(phi)+rowSums(phi))^2)/(4*S^2/M),4*S^2/M)
})
rowMeans(cross) # squared terms are 1/M cross terms

M <- 100
cross <- replicate(1e2,{
x <- rnorm(M); y <- rnorm(M)
phi <- t(sapply(x,function(x.i)as.numeric(x.i<y)))
theta.hat <- sum(phi)/M^2; t <- sum(diag(phi))
S <- sum(phi)
c(4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2   -2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))*((2*M-1)/(M*(M-1))^2 + 1/M^2) + (2*M-1)*t/(M*(M-1))^2 + (2*M-1)/(M*(M-1))^2*sum((colSums(phi)+rowSums(phi))^2),
4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2   -2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))*((2*M-1)/(M*(M-1))^2 + 1/M^2) + (2*M-1)*t/(M*(M-1))^2 + (2*M-1)/(M*(M-1))^2*4*S^2/M)
})
rowMeans(abs(cross)) # squared terms are 1/M cross terms

4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2   -2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))*((2*M-1)/(M*(M-1))^2 + 1/M^2) + (2*M-1)*t/(M*(M-1))^2 + (2*M-1)/(M*(M-1))^2*sum((colSums(phi)+rowSums(phi))^2) # 5b
4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2   -2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))*((2*M-1)/(M*(M-1))^2 + 1/M^2) + (2*M-1)*t/(M*(M-1))^2 + (2*M-1)/(M*(M-1))^2*4*S^2/M # apprx by S--too crude

M <- 90
terms <- replicate(1e2,{
    t <- sample(0:M,1)
    phi <- t(sapply(1:M, function(m)c(rep(0,m-1),rep(0+(M-m<t),M-m+1))))
    S <- sum(phi)
    theta.hat <- S/M^2
    c(4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2, - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))), -    sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2) )
})
rowMeans(abs(terms))
cor(t(terms))
var(t(terms))


ps.to.m <- function(ps) t(sapply(ps, function(p) c(rep(0,p-1),rep(1,M+1-p))))
M <- 6
ps.list <- expand.grid(rep(list(1:(M+1)),M))
diffs <- apply(ps.list,1,function(ps) {
    t <- sum(ps<=(1:M))
    S <- sum(M+1 - ps)
    phi <- ps.to.m(ps)
    c(t,- 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) - sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2))
    c(t,sum((colSums(phi)+rowSums(phi)-diag(phi))^2))
    ## c(t, 4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2 - t^2*M/(M-1)^2/(1-2*M)  - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) - t^2/M/(M-1)^2 + t/M^2 -    sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2))
})
by.t <- split(cbind(diffs=diffs[2,],ps.list),diffs[1,])
max.phis <- lapply(by.t,function(ps){
    ## ps.to.m(ps[which.max(ps$diffs),-1])
    lapply(split(ps[ps$diffs==max(ps$diffs),-1],1:sum(ps$diffs==max(ps$diffs))),ps.to.m)
})
min.phis <- lapply(by.t,function(ps){
    lapply(split(ps[ps$diffs==min(ps$diffs),-1],1:sum(ps$diffs==min(ps$diffs))),ps.to.m)
})



## full objective
M <- 50
t <- sample(0:M,1)
phi <- t(sapply(1:M, function(m)c(rep(0,m-1),rep(0+(M-m<t),M-m+1))))
S <- sum(phi)
theta.hat <- S/M^2
4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2   - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) -  sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2) + (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2
4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2   - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) -  sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2) ## without parts only depending on t
4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2*(1/(1-2*M)))^2 ## O(1) apprx of first term
-2/M^4*(2*S - t*M/2)^2  
## - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))
## -  sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2)
c(4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2,   - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) ,-  sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2)) ## comparing magnitues of terms -- first looks 1/M
c(- 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) ,-  sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2))
c(4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2,   -2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))*((2*M-1)/(M*(M-1))^2 + 1/M^2), + (2*M-1)*t/(M*(M-1))^2 + (2*M-1)/(M*(M-1))^2*sum((colSums(phi)+rowSums(phi))^2)) # 5b

M <- 100
terms <- replicate(1e3,{
    t <- sample(0:M,1)
    phi <- t(sapply(1:M, function(m)c(rep(0,m-1),rep(0+(M-m<t),M-m+1))))
    S <- sum(phi)
    theta.hat <- S/M^2
c(   -2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))*((2*M-1)/(M*(M-1))^2 + 1/M^2), (2*M-1)/(M*(M-1))^2*sum((colSums(phi)+rowSums(phi))^2)) # 5b
})
rowMeans(abs(terms))
cor(t(terms))
var(t(terms))

c(   -2*sum(diag(phi)*(colSums(phi)+rowSums(phi)))*((2*M-1)/(M*(M-1))^2 + 1/M^2), + (2*M-1)*t/(M*(M-1))^2 + (2*M-1)/(M*(M-1))^2*sum((colSums(phi)+rowSums(phi))^2)) # 5b


## delta parametrization

## deltas.to.m <- function(ps) t(sapply(ps, function(p) c(rep(0,p-1),rep(1,M+1-p))))
M <- 10
## deltas <- numeric()
## for (j in 1:M) deltas <- c(deltas,sample(0:(M-sum(deltas)),1))
deltas <- numeric(M)
repeat {
    deltas <<- sample(0:M,M,replace=TRUE)
    if (sum(deltas)<=M) break
    }
## deltas <- sort(deltas)
pi <- sample(1:M)
phi <- sapply(cumsum(deltas), function(c)c(rep(1,c),rep(0,M-c)))[pi,]
S <- sum(phi)
t <- sum(diag(phi))

stopifnot(sum(phi)==sum(cumsum(deltas)))
stopifnot(colSums(phi)==cumsum(deltas))
rowSums(phi)
deltas[pi]*(M+1-pi)



## objective invariant under transposing phi
M <- 20
phi <- matrix(rnorm(M^2),nrow=M)
sapply(list(phi,t(phi)),function(phi) {
    t <- sum(diag(phi)); S <- sum(phi)
    4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2   - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) -  sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2) + (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2
    })



ps.to.m <- function(ps) t(sapply(ps, function(p) c(rep(0,p-1),rep(1,M+1-p))))
M <- 3
ps.list <- expand.grid(rep(list(1:(M+1)),M))
vecs <- apply(ps.list,1,function(ps) {
    phi <- ps.to.m(ps)
    colSums(phi)+rowSums(phi)
})
sort(unique(apply(vecs,2,function(r)paste(as.character(r),collapse=''))))


M <- 10
beta0 <- rnorm(1); beta1 <- rnorm(1); beta2 <- rnorm(1); beta3 <- rnorm(1)
beta0 <- 4/M^3*(1-(M/(M-1))^2); beta3 <- M^2/2/(1-2*M); beta1 <- -(1/M^2 - 1/(M-1)^2); beta2 <- - 2/M^2
beta5 <- beta0*(beta4+1)/(beta0/4+beta1); beta6 = beta2/(beta0/4+beta1)
x <- rnorm(M); y <- rnorm(M)
phi <- t(sapply(x,function(x.i)as.numeric(x.i<y)))
theta.hat <- sum(phi)/M^2; t <- sum(diag(phi))
S <- sum(phi)
gamma <- colSums(phi)+rowSums(phi)-2*diag(phi)
4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2   - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) -  sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2) #+ (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2
## 6a
beta0*(S+beta3*t)^2 + beta1*sum((rowSums(phi)+colSums(phi)-diag(phi))^2) + beta2*sum(diag(phi)*(rowSums(phi)+colSums(phi)))
   beta0*(S+beta3*t)^2 + beta1*sum(gamma^2) + (2*beta1+beta2)*sum(diag(phi)*gamma) + (2*beta2+beta1)*t
beta0*(1/2*sum(gamma)+(beta3+1)*t)^2  + beta1*sum(gamma^2) + (2*beta1+beta2)*sum(diag(phi)*gamma) + (2*beta2+beta1)*t
beta0/4*sum(gamma)^2 + sum((beta0*(beta3+1)*t + (2*beta1+beta2)*diag(phi))*gamma) + beta0*(beta3+1)^2*t^2 + (2*beta2+beta1)*t + beta1*sum(gamma^2)
beta0/4*sum(gamma)^2 + beta1*sum((gamma+beta0/beta1/2*(beta3+1)*t+(1+beta2/2/beta1)*diag(phi))^2) + beta0*(beta3+1)^2*t^2 + (2*beta2+beta1)*t - beta1*sum( (beta0/beta1/2*(beta3+1)*t+(1+beta2/2/beta1)*diag(phi))^2)
## 6b t terms
beta0*(beta3+1)^2*t^2 + (2*beta2+beta1)*t - beta1*sum( (beta0/beta1/2*(beta3+1)*t+(1+beta2/2/beta1)*diag(phi))^2)
(beta0*(beta3+1)^2 - M*beta0^2/4/beta1*(beta3+1)^2 - beta0*(beta3+1)*(1+beta2/2/beta1))*t^2 + (2*beta2+beta1-beta1*(1+beta2/2/beta1)^2)*t
(beta0*(beta3+1)^2*(1-M*beta0/4/beta1) - beta0*(beta3+1)*(1+beta2/2/beta1))*t^2 + beta2*(1-beta2/4/beta1)*t
beta0/4/beta1*(beta3+1)*(-M*beta0 +4*beta1*beta3-M*beta0*beta3-2*beta2)*t^2 + beta2*(1-beta2/4/beta1)*t
## 6c
beta0/4*sum(gamma)^2 + beta1*sum((gamma+beta0/beta1/2*(beta3+1)*t+(1+beta2/2/beta1)*diag(phi))^2)
beta0/beta1/4*sum(gamma)^2 + sum((gamma+beta0/beta1/2*(beta3+1)*t+(1+beta2/2/beta1)*diag(phi))^2)
uu <- gamma+beta0/beta1/2*(beta3+1)*t+(1+beta2/2/beta1)*diag(phi); sum((uu-mean(uu))^2); (M-1)*var(uu)
## 6d
beta0*(S+beta3*t)^2 + beta1*sum((rowSums(phi)+colSums(phi)-diag(phi))^2) + beta2*sum(diag(phi)*(rowSums(phi)+colSums(phi)))
beta1*(beta0/beta1/4*sum(gamma)^2 + sum((gamma+beta0/beta1/2*(beta3+1)*t+(1+beta2/2/beta1)*diag(phi))^2)) + beta0*(beta3+1)^2*t^2 + (2*beta2+beta1)*t - beta1*sum( (beta0/beta1/2*(beta3+1)*t+(1+beta2/2/beta1)*diag(phi))^2)
beta1*(M-1)*var(gamma+beta0/beta1/2*(beta3+1)*t+(1+beta2/2/beta1)*diag(phi)) + beta0/4/beta1*(beta3+1)*(-M*beta0 +4*beta1*beta3-M*beta0*beta3-2*beta2)*t^2 + beta2*(1-beta2/4/beta1)*t
## full objective with variance term
c(4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2   - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) -  sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2), + (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2)
c(beta0*(S+beta3*t)^2 + beta1*sum((rowSums(phi)+colSums(phi)-diag(phi))^2) + beta2*sum(diag(phi)*(rowSums(phi)+colSums(phi))) ,+ (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2)
c(beta1*(M-1)*var(gamma+(1+beta2/2/beta1)*diag(phi)) + beta0/4/beta1*(beta3+1)*(-M*beta0 +4*beta1*beta3-M*beta0*beta3-2*beta2)*t^2 + beta2*(1-beta2/4/beta1)*t,+ (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2)
c(beta1*(M-1)*var(gamma+(1+beta2/2/beta1)*diag(phi)), + beta0/4/beta1*(beta3+1)*(-M*beta0 +4*beta1*beta3-M*beta0*beta3-2*beta2)*t^2 + beta2*(1-beta2/4/beta1)*t + (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2) ## t terms all together
beta1*(M-1)*var(gamma+(1+beta2/2/beta1)*diag(phi)) + beta0/4/beta1*(beta3+1)*(-M*beta0 +4*beta1*beta3-M*beta0*beta3-2*beta2)*t^2 + beta2*(1-beta2/4/beta1)*t + (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2
(2*M-1)/M^2/(M-1)*var(gamma + (-M^2+4*M-2)/(2*M-1)*diag(phi)) + 1/M/(2*M-1)*(t-(M^2-M-4)/2/(M-1))^2 + t/M^2/(M-1)^2/(2*M-1) - (M^2-M-4)^2/4/M/(2*M-1)/(M-1)^2 # 6h
(2*M-1)/M^2/(M-1)*var(gamma + (-M^2+4*M-2)/(2*M-1)*diag(phi)) + 1/M/(2*M-1)*(t-(M^2-M-4)/2/(M-1))^2 + t/M^2/(M-1)^2/(2*M-1) - (M^2-M-4)^2/4/M/(2*M-1)/(M-1)^2 # 6h
(2*M-1)/M^2/(M-1)*var(colSums(phi)+rowSums(phi) - M^2/(2*M-1)*diag(phi)) + 1/M/(2*M-1)*(t-(M^2-M-4)/2/(M-1))^2 + t/M^2/(M-1)^2/(2*M-1) - (M^2-M-4)^2/4/M/(2*M-1)/(M-1)^2 
(2*M-1)/M^2/(M-1)*var(colSums(phi)+rowSums(phi) - M/2*diag(phi)) + 1/M/(2*M-1)*(t-(M^2-M-4)/2/(M-1))^2 + t/M^2/(M-1)^2/(2*M-1) - (M^2-M-4)^2/4/M/(2*M-1)/(M-1)^2 # approximating diagonal term by M/2, a little crude


beta1*(M-1)*var(gamma+(1+beta2/2/beta1)*diag(phi))
var(gamma+(1+beta2/2/beta1)*diag(phi))

## verifying formula with variance term
ps.to.m <- function(ps) t(sapply(ps, function(p) c(rep(0,p-1),rep(1,M+1-p))))
M <- 5
ps.list <- expand.grid(rep(list(1:(M+1)),M))
diffs <- apply(ps.list,1,function(ps) {
    t <- sum(ps<=(1:M))
    S <- sum(M+1 - ps)
    phi <- ps.to.m(ps)
    gamma <- colSums(phi)+rowSums(phi)-2*diag(phi)
    beta0 <- 4/M^3*(1-(M/(M-1))^2); beta3 <- M^2/2/(1-2*M); beta1 <- -(1/M^2 - 1/(M-1)^2); beta2 <- - 2/M^2
    c(4/M^3*(1-(M/(M-1))^2)*(S+t*M^2/2/(1-2*M))^2   - 2/M^2*sum(diag(phi)*(colSums(phi)+rowSums(phi))) -  sum((colSums(phi)+rowSums(phi)-diag(phi))^2)*(1/M^2 - 1/(M-1)^2) + (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2,
      beta1*(M-1)*var(gamma+(1+beta2/2/beta1)*diag(phi)) + beta0/4/beta1*(beta3+1)*(-M*beta0 +4*beta1*beta3-M*beta0*beta3-2*beta2)*t^2 + beta2*(1-beta2/4/beta1)*t + (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2)
})
plot(t(diffs))
sum(diffs[1,]-diffs[2,])


diffs <- apply(ps.list,1,function(ps) {
    t <- sum(ps<=(1:M))
    S <- sum(M+1 - ps)
    phi <- ps.to.m(ps)
    gamma <- colSums(phi)+rowSums(phi)-2*diag(phi)
    beta0 <- 4/M^3*(1-(M/(M-1))^2); beta3 <- M^2/2/(1-2*M); beta1 <- -(1/M^2 - 1/(M-1)^2); beta2 <- - 2/M^2
    c(beta1*(M-1)*var(gamma+(1+beta2/2/beta1)*diag(phi)), + beta0/4/beta1*(beta3+1)*(-M*beta0 +4*beta1*beta3-M*beta0*beta3-2*beta2)*t^2 + beta2*(1-beta2/4/beta1)*t + (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2) ## t terms all together
})
rowMeans(abs(diffs))
mean(abs(diffs[1,]) - abs(diffs[2,]))
plot(abs(diffs[1,]),abs(diffs[2,]))
hist(abs(diffs[1,]-diffs[2,]))



## matrix maximizing variance of colsums+rowsums
ps.to.m <- function(ps) t(sapply(ps, function(p) c(rep(0,p-1),rep(1,M+1-p))))
M <- 7
ps.list <- expand.grid(rep(list(1:(M+1)),M))
vars <- apply(ps.list,1,function(ps) {
    t <- sum(ps<=(1:M))
    S <- sum(M+1 - ps)    
    phi <- ps.to.m(ps)
    ## c(t,var(colSums(phi)+rowSums(phi)))
    gamma <- colSums(phi)+rowSums(phi)-2*diag(phi)
    beta0 <- 4/M^3*(1-(M/(M-1))^2); beta3 <- M^2/2/(1-2*M); beta1 <- -(1/M^2 - 1/(M-1)^2); beta2 <- - 2/M^2
    ## c(t,var(gamma+(1+beta2/2/beta1)*diag(phi)))
    c(t,(2*M-1)/M^2/(M-1)*var(gamma + (-M^2+4*M-2)/(2*M-1)*diag(phi)) + 1/M/(2*M-1)*(t-(M^2-M-4)/2/(M-1))^2 + t/M^2/(M-1)^2/(2*M-1) - (M^2-M-4)^2/4/M/(2*M-1)/(M-1))
    ## c(t,var(colSums(phi)+rowSums(phi)-M^2/(2*M-1)*diag(phi)))
    })
## ps.to.m(ps.list[which.max(vars),])
lapply(split(ps.list[vars[2,]==max(vars[2,]),],1:sum(vars[2,]==max(vars[2,]))),ps.to.m)
lapply(split(ps.list[vars[2,]==min(vars[2,]),],1:sum(vars[2,]==min(vars[2,]))),ps.to.m)
lapply(split(ps.list[abs(vars[2,])==max(abs(vars[2,])),],1:sum(abs(vars[2,])==max(abs(vars[2,])))),ps.to.m)
lapply(split(ps.list[abs(vars[2,])==max(abs(vars[2,])),],1:sum(abs(vars[2,])==max(abs(vars[2,])))),function(m){phi <- ps.to.m(m);sum(diag(phi))})

by.t <- split(cbind(vars=vars[2,],ps.list),vars[1,])
max.phis <- lapply(by.t,function(ps){
    ## ps.to.m(ps[which.max(ps$vars),-1])
    lapply(split(ps[abs(ps$vars-max(ps$vars))<1e-5,-1],1:sum(abs(ps$vars-max(ps$vars))<1e-5)),ps.to.m)
})
min.phis <- lapply(by.t,function(ps){
    lapply(split(ps[abs(ps$vars-min(ps$vars))<1e-5,-1],1:sum(abs(ps$vars-min(ps$vars))<1e-5)),ps.to.m)
})

lapply(max.phis, function(l)lapply(l,function(ll)list(phi=ll,trace=sum(diag(ll)),gg=colSums(ll)+rowSums(ll)-M^2/(2*M-1)*diag(ll),hh=var(colSums(ll)+rowSums(ll)-M^2/(2*M-1)*diag(ll)))))
tr <- function(m){m[m==1] <- 2;m[m==0] <- 1;m[m==2] <- 0;m[nrow(m):1,][,ncol(m):1]}
ff <- function(m){var(colSums(m)+rowSums(m)-M^2/(2*M-1)*diag(m))}

## with(list(t=3),{
##     plot(rep(1,length(by.t[[t]]$vars)),by.t[[3]]$vars-1/M/(2*M-1)*(t-(M^2-M-4)/2/(M-1 ))^2 + t/M^2/(M-1)^2/(2*M-1) - (M^2-M-4)^2/4/M/(2*M-1)/(M-1))
##     abline(h=1/M/(2*M-1)*(t-(M^2-M-4)/2/(M-1 ))^2 + t/M^2/(M-1)^2/(2*M-1) - (M^2-M-4)^2/4/M/(2*M-1)/(M-1))
##     } )


## reducing t terms
o+ beta0/4/beta1*(beta3+1)*(-M*beta0 +4*beta1*beta3-M*beta0*beta3-2*beta2)*t^2 + beta2*(1-beta2/4/beta1)*t + (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2
beta0/4/beta1*(1+beta3)
-1/M*(M^2-4*M+2)/2/(1-2*M)
-M*beta0 +4*beta1*beta3-M*beta0*beta3-2*beta2
beta2*(1-beta2/4/beta1)
-1/M^2*(M^2+2*M-1)/(2*M-1)
beta2*(1-beta2/4/beta1) + (2*M-1)/M^2/(M-1)^2 + 1/M^2
(-M^4+2*M^3+3*M^2-4*M+1)/M^2/(M-1)^2/(2*M-1)
1/M^2/(M-1)^2/(2*M-1) - (M^2-M-4)/M/(M-1)/(2*M-1) #6e
-M/(M-1)^2/(1-2*M) - 1/M/(M-1)^2
1/M/(2*M-1)
+ beta0/4/beta1*(beta3+1)*(-M*beta0 +4*beta1*beta3-M*beta0*beta3-2*beta2)*t^2 + beta2*(1-beta2/4/beta1)*t + (2*M-1)*t/(M*(M-1))^2 - t^2*M/(M-1)^2/(1-2*M)- t^2/M/(M-1)^2 + t/M^2
t^2/M/(2*M-1) + t*(1/M^2/(M-1)^2/(2*M-1) - (M^2-M-4)/M/(M-1)/(2*M-1)) # 6f
1/M/(2*M-1)*(t-(M^2-M-4)/2/(M-1))^2 - 1/M/(2*M-1)*((M^2-M-4)/2/(M-1))^2 #+ t/M^2/(M-1)^2/(2*M-1)

(2*M-1)/M^2/(M-1)*var(

## checking trace formula for var(colSums + rowSums)
M <- 10
ones <- matrix(1,nrow=M,ncol=1)
x <- rnorm(M); y <- rnorm(M)
x <- x[order(y)]; y <- sort(y)
phi <- t(sapply(x,function(x.i)as.numeric(x.i<y)))
A <- phi + t(phi)
diag(A) <- 0
sum(diag(A%*%A))
sum((colSums(phi)+rowSums(phi))^2)
sum(diag( t(ones)%*%(A%*%A)%*%ones))
sum(A%*%A)
A

## comparing phi-dependent terms to trace terms
ps.to.m <- function(ps) t(sapply(ps, function(p) c(rep(0,p-1),rep(1,M+1-p))))
tr <- function(m){m[m==1] <- 2;m[m==0] <- 1;m[m==2] <- 0;m[nrow(m):1,][,ncol(m):1]}
M <- 5
ps.list <- expand.grid(rep(list(1:(M+1)),M))
vars <- apply(ps.list,1,function(ps) {
    t <- sum(ps<=(1:M))
    S <- sum(M+1 - ps)    
    phi <- ps.to.m(ps)
    ## c(t,var(colSums(phi)+rowSums(phi)))
    gamma <- colSums(phi)+rowSums(phi)-2*diag(phi)
    beta0 <- 4/M^3*(1-(M/(M-1))^2); beta3 <- M^2/2/(1-2*M); beta1 <- -(1/M^2 - 1/(M-1)^2); beta2 <- - 2/M^2
    ## c(t,var(gamma+(1+beta2/2/beta1)*diag(phi)))
    c(t,var(gamma + (-M^2+4*M-2)/(2*M-1)*diag(phi)))#, + 1/M/(2*M-1)*(t-(M^2-M-4)/2/(M-1))^2 + t/M^2/(M-1)^2/(2*M-1) - (M^2-M-4)^2/4/M/(2*M-1)/(M-1))
    ## c(t,var(colSums(phi)+rowSums(phi)-M^2/(2*M-1)*diag(phi)))
})
by.t <- split(cbind(vars=vars[2,],ps.list),vars[1,])
with(list(t=2),{hist(by.t[[t]]$vars);abline(v=1/M/(2*M-1)*(t-(M^2-M-4)/2/(M-1))^2 + t/M^2/(M-1)^2/(2*M-1) - (M^2-M-4)^2/4/M/(2*M-1)/(M-1),col='red')})
curve(1/M/(2*M-1)*(t-(M^2-M-4)/2/(M-1))^2 + t/M^2/(M-1)^2/(2*M-1) - (M^2-M-4)^2/4/M/(2*M-1)/(M-1),xname='t',from=0,to=M)
with(list(t=2),{plot(by.t[[t]]$vars,ylim=c(-.5,.5));points(by.t[[t]]$vars + 1/M/(2*M-1)*(t-(M^2-M-4)/2/(M-1))^2 + t/M^2/(M-1)^2/(2*M-1) - (M^2-M-4)^2/4/M/(2*M-1)/(M-1),col='red')})
tr(ps.to.m(ps.list[which.max(vars[2,]),]))
ps.to.m(ps.list[which.max(vars[2,]),])

M <- 10
j <- 3
ps <- c(rep(M-j+1,M-j-1),M-j,rep(1,j))
phi <- ps.to.m(ps)
x <- colSums(phi)+rowSums(phi) - M^2/(2*M-1)*diag(phi)
alpha <- M^2/(2*M-1)
sum(x^2)
-4*j^3 + 4*M*j^2 + j*(4*M^2-4*(M+1)*alpha+8+alpha^2) + (alpha-2)^2
(2*M-alpha)^2*j+4*(j+1)^2-4*alpha*(j+1)+alpha^2+4*j^2*(M-j-1)
ff <- function(M,j)-4*j^3 + 4*M*j^2 + j*(4*M^2-4*(M+1)*alpha+8+alpha^2) + (alpha-2)^2
## curve(-4*j^3 + 4*M*j^2 + j*(4*M^2-4*(M+1)*alpha+8+alpha^2) + (alpha-2)^2,xname='j',-10,5);abline(h=0)
## M*(M-1)*var(x)
M*sum(x^2)-sum(x)^2
-4*j^4 + j^3*(12*M-4*alpha)+j^2*(-12*M^2+8*M*alpha-alpha^2-4*alpha+8) + j*(4*M^3-4*M*(M+1)*alpha + M*(alpha^2+8) - 2*(4*M-alpha)*(2-alpha)) + (alpha-2)^2*(M-1)
curve(-4*j^4 + j^3*(12*M-4*alpha)+j^2*(-12*M^2+8*M*alpha-alpha^2-4*alpha+8) + j*(4*M^3-4*M*(M+1)*alpha + M*(alpha^2+8) - 2*(4*M-alpha)*(2-alpha)) + (alpha-2)^2*(M-1),xname='j',-2.5,10);abline(h=0)
-4*j^4 + 4*j^3*(5*M^2-3*M)/(2*M-1) + j^2/(2*M-1)^2*(-33*M^4+32*M^3+24*M^2-32*M+8) + j*(4*M^3-4*M*alpha*(M+1)+M*(alpha^2+8)-2*(4*M-alpha)*(2-alpha)) + (alpha-2)^2*(M-1)
ff <- function(j,M)-4*j^4 + 4*j^3*(5*M^2-3*M)/(2*M-1) + j^2/(2*M-1)^2*(-33*M^4+32*M^3+24*M^2-32*M+8) + j*1/(2*M-1)^2*M*(9*M^4-6*M^3-24*M^2+28*M-8) + (M^2/(2*M-1)-2)^2*(M-1)
ff(1e4,1e4+1)
M <- 1e4
curve(-4*j^4 + 4*j^3*(5*M^2-3*M)/(2*M-1) + j^2/(2*M-1)^2*(-33*M^4+32*M^3+24*M^2-32*M+8) + j*1/(2*M-1)^2*M*(9*M^4-6*M^3-24*M^2+28*M-8) + (M^2/(2*M-1)-2)^2*(M-1),xname='j',9.5e3,1.1e4);abline(h=0)
-4*j^4 + 4*j^3*(5*M^2-3*M)/(2*M-1) + j^2/(2*M-1)^2*(-33*M^4+32*M^3+24*M^2-32*M+8) + j*1/(2*M-1)^2*M*(9*M^4-6*M^3-24*M^2+28*M-8) + (M^2/(2*M-1)-2)^2*(M-1)
(-4*j^3 + 4*(3*M^2-1)/(2*M-1)*j^2 +1/(2*M-1)^2*(-9*M^4-4*M^3+28*M^2-20*M+4)*j - 1/(2*M-1)^2*(M^4-8*M^3+20*M^2-16*M+4))*(j-(M-1)) ##7a
curve(-4*j^3 + 4*(3*M^2-1)/(2*M-1)*j^2 +1/(2*M-1)^2*(-9*M^4-4*M^3+28*M^2-20*M+4)*j - 1/(2*M-1)^2*(M^4-8*M^3+20*M^2-16*M+4),ylab='',xname='j',14,18);abline(h=0)
-4*(2*M-1)^2*j^3 + 4*(3*M^2-1)*(2*M-1)*j^2 +(-9*M^4-4*M^3+28*M^2-20*M+4)*j - (M^4-8*M^3+20*M^2-16*M+4)
curve(-4*(2*M-1)^2*j^3 + 4*(3*M^2-1)*(2*M-1)*j^2 +(-9*M^4-4*M^3+28*M^2-20*M+4)*j - (M^4-8*M^3+20*M^2-16*M+4),xname='j',-1,1);abline(h=0)
-16*(2*M-1)^2*j^3+12*M*(2*M-1)*(5*M-3)*j^2+2*(-33*M^4+32*M^3+24*M^2-32*M+8)*j + M*(9*M^4-6*M^3-24*M^2+28*M-8)
curve(-16*(2*M-1)^2*j^3+12*M*(2*M-1)*(5*M-3)*j^2+2*(-33*M^4+32*M^3+24*M^2-32*M+8)*j + M*(9*M^4-6*M^3-24*M^2+28*M-8),xname='j',0,10);abline(h=0)

a <- 2;b <- 2
M <- 10
obj <- function(phi) {
    M <- nrow(phi)
    print(colSums(phi)+rowSums(phi)-M^2/(2*M-1)*diag(phi))
    cat('var: ',var(colSums(phi)+rowSums(phi)-M^2/(2*M-1)*diag(phi)),'\n')
    cat('mean: ',mean(colSums(phi)+rowSums(phi)-M^2/(2*M-1)*diag(phi)),'\n')
}
gen <- function(M=10,a=2,b=2,print=TRUE) {
    phi <- matrix(0,nrow=M,ncol=M)
    phi[(M-a+1):M,] <- 1
    phi[,(M-a+1):M] <- 1
    phi[(M-a-b+1):M,(M-a-b+1):M] <- 1
    obj(phi)
    if(print){
        print(phi)
    } else {
        return(phi)
    }
}




x <- c(rep(3/2*M,j),rep(j-M/2,M-j))
sum(x)
-M^2/2+3*M*j-j^2
sum(x^2)
M^3/4+M^2*j+2*M*j^2-j^3
M^4/4+9*M^2*j^2+j^4-3*M^3*j+M^2*j^2-6*M*j^3
sum(x)^2
M*sum(x^2)-sum(x)^2
j*(4*M^3)+j^2*(-8*M^2)+5*M*j^3 - j^4

curve(-4*j^3+15*M*j^2-16*M^2*j+4*M^3,xname='j',-20,40);abline(h=0)
abline(v=-1/8*(-7+c(-1,1)*sqrt(17))*M)
curve(j*(4*M^3)+j^2*(-8*M^2)+5*M*j^3 - j^4,xname='j',00,40)
abline(v=-1/8*(-7+c(-1,1)*sqrt(17))*M)


