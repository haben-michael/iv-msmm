## updated to reflect simplifed rbinormal routine (eg redundant random effect Z removed)

start <- Sys.time()
## require(mvtnorm)
source('../../misc.R')
I <- 10
alpha <- .05
B <- 2e2
theta.11 <- seq(1/2,1,length.out=B)
theta.12 <- seq(0.501,.95,length.out=B) # could change to just \theta.12\in [1/2, 1]
mean.x <- 0
var.x <- var.y <- 1 #; var.z <- 0
## reject.rates <- matrix(NA,B,B)
reject <- as.list(numeric(B^2))
dim(reject) <- c(B,B)
for(i in 1:B)
    for(j in 1:B) {
        reject[i,j] <- NA
        if(qnorm(theta.11[i])*qnorm(theta.12[j])<0) next
        if(theta.12[j] < 1-pnorm(sqrt(2)*qnorm(theta.11[i])) || theta.12[j] > pnorm(sqrt(2)*qnorm(theta.11[i]))) next
        ## c12 <- 2*var.z/(var.x+var.y)
        ## c11 <- (1+c12)*(qnorm(theta.12[j])/qnorm(theta.11[i]))^2-1
        ## cov.xy <- -(var.x+var.y)/2*c11
        ## mean.y <- qnorm(theta.11[i])*sqrt(1+c11)*sqrt(var.x+var.y)
        ## delta <- (mean.y-mean.x)/sqrt(var.x+var.y)
                ## browser()
        mean.y <- delta <- sqrt(2)*qnorm(theta.12[j])
        cov.xy <- 1 - (mean.y / qnorm(theta.11[i]))^2 / 2
        cov.xx <- cov.yy <- abs(cov.xy) ## ancillary but must be >= cov.xy for psd vcov
        if(abs(cov.xy)>mean(c(var.x,var.y))-.01) next
        ## if(abs(theta.11[i]-pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy)))>1e-5)browser()
        if(abs(theta.11[i]-pnorm(mean.y/sqrt(2*(1-cov.xy))))>1e-5)browser()
        ## if(abs(theta.12[j]-pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z)))>1e-5) browser()
        if(abs(cov.xy)>mean(c(var.x,var.y))-.01) next
        Delta.rho <- auc.to.params.binormal(theta.12[j],theta.11[i])
        mean.y <- Delta.rho['Delta']
        rho <- Delta.rho['rho']
        ## browser()
        reject[[i,j]] <- replicate(2e1, {
            ## cat('.')
            m <- 1+rpois(I,5); n <- 1+rpois(I,5)
            data <- rbinormal(I=I,m=m,n=n,mean.x=mean.x,mean.y=mean.y,cov.xx=cov.xx,cov.xy=cov.xy,cov.yy=cov.yy,plot=FALSE)
            x <- data$x; y <- data$y
            if(length(x)!=length(y)) {
                print(x)
                print(y)
                }
            contrast <- matrix(c(1,-1),nrow=2)
            z.stat <- with(auc.cluster(x,y), sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
            abs(z.stat) > qnorm(1-alpha/2)
        })
    }
elapsed <- Sys.time() - start
print(elapsed)


filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save.image(filename)













## i <- j <- 1
## while(TRUE) {
##     reject[i,j] <- NA
##         if(qnorm(theta.11[i])*qnorm(theta.12[j])<0) next
##         if(theta.12[j] < 1-pnorm(sqrt(2)*qnorm(theta.11[i])) || theta.12[j] > pnorm(sqrt(2)*qnorm(theta.11[i]))) next
##         ## c12 <- 2*var.z/(var.x+var.y)
##         ## c11 <- (1+c12)*(qnorm(theta.12[j])/qnorm(theta.11[i]))^2-1
##         ## cov.xy <- -(var.x+var.y)/2*c11
##         ## mean.y <- qnorm(theta.11[i])*sqrt(1+c11)*sqrt(var.x+var.y)
##         ## delta <- (mean.y-mean.x)/sqrt(var.x+var.y)
##                 ## browser()
##         mean.y <- delta <- sqrt(2)*qnorm(theta.12[j])
##         cov.xy <- 1 - (mean.y / qnorm(theta.11[i]))^2 / 2
##         cov.xx <- cov.yy <- abs(cov.xy) ## ancillary but must be >= cov.xy for psd vcov
##         if(abs(cov.xy)>mean(c(var.x,var.y))-.01) next
##         ## if(abs(theta.11[i]-pnorm((mean.y-mean.x)/sqrt(var.x+var.y-2*cov.xy)))>1e-5)browser()
##         if(abs(theta.11[i]-pnorm(mean.y/sqrt(2*(1-cov.xy))))>1e-5)browser()
##         ## if(abs(theta.12[j]-pnorm((mean.y-mean.x)/sqrt(var.x+var.y+2*var.z)))>1e-5) browser()
##         if(abs(cov.xy)>mean(c(var.x,var.y))-.01) next
##         reject[[i,j]] <- replicate(2e1, {
##             ## cat('.')
##             m <- 1+rpois(I,5); n <- 1+rpois(I,5)
##             data <- rbinormal(I=I,m=m,n=n,mean.x=mean.x,mean.y=mean.y,cov.xx=cov.xx,cov.xy=cov.xy,cov.yy=cov.yy,plot=FALSE)
##             x <- data$x; y <- data$y
##             if(length(x)!=length(y)) browser()
##             contrast <- matrix(c(1,-1),nrow=2)
##             z.stat <- with(auc.cluster(x,y), sqrt(I) * (t(contrast)%*%c(theta.11.hat,theta.12.hat)) / sqrt(t(contrast)%*%vcov.hat%*%contrast))
##             abs(z.stat) > qnorm(1-alpha/2)
##         })
##     }
