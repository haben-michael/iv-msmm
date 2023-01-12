## (8/11/22) updated to reflect simplifed rbinormal routine (eg
## redundant random effect Z removed) (8/13/22) refactoring to allow
## dgp to be binormal #1 without or #2 with censoring. interpolation
## to allow irregular grid. [outcome: did not seem to improve image to
## use random points/irregular grid. Also needed to use akima::interp
## which was going too slow on ~1e6 points. update: worked much better
## after rounding x and y values, and applying moving averages to the
## power.]

start <- Sys.time()
## require(mvtnorm)
source('../../misc.R')



## bound <- abs(qnorm(.2)) # 2
## auc.to.params <- function(theta.12,theta.11)auc.to.params.binormal.censored(theta.12,theta.11,bound) #2
## sampler <- rbinormal.censored #2
auc.to.params <- auc.to.params.binormal # binormal #1
sampler <- rbinormal #1
I <- 10
args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0) B <- 1e3 else B <- as.numeric(args[1])
## B <- 1e4
alpha <- .05
## theta.11 <- seq(1/2,1,length.out=B)
## theta.12 <- seq(0.501,.95,length.out=B) # could change to just \theta.12\in [1/2, 1]
mean.x <- 0
var.x <- var.y <- 1 #; var.z <- 0
## reject.rates <- matrix(NA,B,B)
reject <- as.list(numeric(B))
theta.12 <- runif(B,.501,.95)
theta.11 <- runif(B,theta.12,1)
## dim(reject) <- c(B,B)
for(i in 1:B) {
    ## for(j in 1:B) {
    ## reject[i,j] <- NA
    if(qnorm(theta.11[i])*qnorm(theta.12[i])<0) next
    if(theta.12[i] < 1-pnorm(sqrt(2)*qnorm(theta.11[i])) || theta.12[i] > pnorm(sqrt(2)*qnorm(theta.11[i]))) next
    ## mean.y <- delta <- sqrt(2)*qnorm(theta.12[j])
    ## cov.xy <- 1 - (mean.y / qnorm(theta.11[i]))^2 / 2
    Delta.rho <- auc.to.params(theta.12=theta.12[i],theta.11=theta.11[i])
    mean.y <- Delta.rho['Delta']
    cov.xy <- Delta.rho['rho']
    cov.xx <- cov.yy <- abs(cov.xy) ## ancillary but must be >= cov.xy for psd vcov
    if(abs(cov.xy)>mean(c(var.x,var.y))-.01) next
    ## if(abs(theta.11[i]-pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy)))>1e-5)browser()
    ## if(abs(theta.11[i]-pnorm(mean.y/sqrt(2*(1-cov.xy))))>1e-5)browser()
    ## if(abs(theta.12[j]-pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z)))>1e-5) browser()
    reject[[i]] <- replicate(1e1, {
        ## cat('.')
        m <- 1+rpois(I,5); n <- 1+rpois(I,5)
        ## data <- rbinormal(I=I,m=m,n=n,mean.x=mean.x,mean.y=mean.y,cov.xx=cov.xx,cov.xy=cov.xy,cov.yy=cov.yy,plot=FALSE)
        data <- sampler(I=I,m=m,n=n,mean.x=mean.x,mean.y=mean.y,cov.xx=cov.xx,cov.xy=cov.xy,cov.yy=cov.yy,plot=FALSE)
        x <- data$x; y <- data$y
        ## if(length(x)!=length(y)) {
        ##     print(x)
        ##     print(y)
        ##     }
        contrast <- matrix(c(1,-1),nrow=2)
        z.stat <- with(auc.cluster(x,y), sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
        abs(z.stat) > qnorm(1-alpha/2)
    })
}
elapsed <- Sys.time() - start
print(elapsed)


filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save.image(filename)












## require(mvtnorm)
## source('../../misc.R')
## bound <- abs(qnorm(.2)) # 2
## auc.to.params <- function(theta.12,theta.11)auc.to.params.binormal.censored(theta.12,theta.11,bound) #2
## ## sampler <- rbinormal.censored #2
## ## auc.to.params <- auc.to.params.binormal # binormal #1
## ## sampler <- rbinormal #1
## ## I <- 10
## ## alpha <- .05
## B <- 2e2
## theta.11 <- seq(1/2,.98,length.out=B)
## theta.12 <- seq(0.501,.95,length.out=B) # could change to just \theta.12\in [1/2, 1]
## ## mean.x <- 0
## ## var.x <- var.y <- 1 #; var.z <- 0
## ## reject.rates <- matrix(NA,B,B)
## reject <- as.list(numeric(B^2))
## dim(reject) <- c(B,B)
## for(i in 1:B)
##     for(j in 1:B) {
##         reject[i,j] <- NA
##         if(qnorm(theta.11[i])*qnorm(theta.12[j])<0) next
##         if(theta.12[j] < 1-pnorm(sqrt(2)*qnorm(theta.11[i])) || theta.12[j] > pnorm(sqrt(2)*qnorm(theta.11[i]))) next
##         ## mean.y <- delta <- sqrt(2)*qnorm(theta.12[j])
##         ## cov.xy <- 1 - (mean.y / qnorm(theta.11[i]))^2 / 2
##         Delta.rho <- auc.to.params(theta.12[j],theta.11[i])
##         mean.y <- Delta.rho['Delta']
##         cov.xy <- Delta.rho['rho']
##         ## cov.xx <- cov.yy <- abs(cov.xy) ## ancillary but must be >= cov.xy for psd vcov
##         ## if(abs(cov.xy)>mean(c(var.x,var.y))-.01) next
##         ## ## if(abs(theta.11[i]-pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy)))>1e-5)browser()
##         ## if(abs(theta.11[i]-pnorm(mean.y/sqrt(2*(1-cov.xy))))>1e-5)browser()
##         ## ## if(abs(theta.12[j]-pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z)))>1e-5) browser()
##         ## if(abs(cov.xy)>mean(c(var.x,var.y))-.01) next
##         ## ## browser()
##         print('.')
##         }
