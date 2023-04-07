## the distribution of the distances of a
## p-dimensional normal vector with mean mu and covariance matrix
## Sigma


.rdist <- function(n,mu,Sigma) {
    eig <- eigen(Sigma) #spectral decomposition
    P <- t(eig$vectors)
    Lambda <- diag(eig$values)
    ncps <- (diag(1/sqrt(diag(Lambda)))%*%P%*%mu)^2 #noncentrality parameters
    W.squared <- replicate(n, sapply(ncps, function(ncp) rchisq(1,df=1,ncp=ncp)))
   sqrt(diag(Lambda)%*%W.squared)
}


rdist <- function(n,mu.1,mu.2,Sigma.1,Sigma.2,weights) {
    weights.mtx <- diag(weights)
    mu <- weights.mtx%*%(mu.2-mu.1)
    Sigma <- weights.mtx%*%(Sigma.1+Sigma.2)%*%weights.mtx
    .rdist(n,mu,Sigma)
    }


require(mvtnorm)
p <- 6
## generate example normal parameters and weights
mu.1 <- runif(p)
mu.2 <- runif(p)
Sigma.1 <- matrix(rnorm(p^2),p)
Sigma.1 <- Sigma.1%*%t(Sigma.1)
Sigma.2 <- matrix(rnorm(p^2),p)
Sigma.2 <- Sigma.2%*%t(Sigma.2)
weights <- runif(p)
weights.mtx <- diag(weights)


x.1 <- rmvnorm(1e4,mu.1,Sigma.1)%*%weights.mtx
x.2 <- rmvnorm(1e4,mu.2,Sigma.2)%*%weights.mtx
squared.distances <- rowSums((x.1-x.2)^2)
hist(sqrt(squared.distances),prob=TRUE)

theoretical <- rdist(1e4,mu.1,mu.2,Sigma.1,Sigma.2,weights)
lines(density(theoretical),col=2)
