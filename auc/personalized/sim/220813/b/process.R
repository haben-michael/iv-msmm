filelist <- dir()
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
xyz <- lapply(filelist, function(file) {#!!!
    load(file)
    ## print(file)
    cbind(theta.11=theta.11,theta.12=theta.12,power=sapply(reject,mean))
    })
## xyz <- as.data.frame(xyz)    
xyz <- as.data.frame(do.call(rbind,xyz))

prec <- 2
xyz$theta.11 <- round(xyz$theta.11,prec)
xyz$theta.12 <- round(xyz$theta.12,prec)
## xyz <- xyz[order(xyz$theta.11,xyz$theta.12),]
xyz <- aggregate(power ~ theta.11 + theta.12, FUN=mean, data=xyz)
theta.11 <- theta.12 <- seq(.5,1,by=10^-prec)
xyz.full <- expand.grid(theta.12=theta.12,theta.11=theta.11)
xyz.full <- merge(xyz.full,xyz,all.x=TRUE)
## xyz.full <- as.data.frame(xyz.full)
power.mat <- array(xyz.full[,'power'], c(length(theta.11),length(theta.12)))
## contour(theta.11,theta.12,power.mat,nlevels=10,method='simple')


bw <- 2
power.mat.ma <- t(apply(power.mat,1, function(r)filter(r,filter=rep(1/(2*bw+1),2*bw+1),'conv')))
power.mat.ma <- apply(power.mat.ma,2, function(c)filter(c,filter=rep(1/(2*bw+1),2*bw+1),'conv'))
## contour(theta.11,theta.12,power.mat.ma,nlevels=10,method='simple')


## ## print(sum(is.na(xyz)))
## ## xyz.interp <- with(na.omit(xyz), akima::interp(x = theta.11, y = theta.12, z = power, linear=TRUE))

## ## load(filelist[1])
## ## reject.total <- as.list(numeric(B^2))
## ## dim(reject.total) <- c(B,B)
## ## for(i in 1:B) for(j in 1:B) reject.total[[i,j]] <- logical()

## ## for(file in filelist) {
## ##     load(file)
## ##     for(i in 1:B) for(j in 1:B) reject.total[[i,j]] <- c(reject.total[[i,j]],reject[[i,j]])
## ##     }

## ## print(length(reject.total[[100,100]]))

## ## power <- matrix(NA,B,B)
## ## for(i in 1:B) for(j in 1:B) power[i,j] <- mean(reject.total[[i,j]])
## with(xyz.interp, image(x,y,z))
## abline(a=0,b=1)#;abline(a=1,b=-1)

## pdf('220813.pdf')
plot(0,type='n',xlim=c(1/2,1),ylim=c(1/2,1),xlab=expression(theta[11]),ylab=expression(theta[12]),xaxs='i',yaxs='i')
contour(theta.11,theta.12,power.mat.ma,add=TRUE,nlevels=10,method='simple')
## contourplot(power)
abline(a=0,b=1)
text(.73,.75,expression(paste(H[0]:theta[11]~"="~theta[12])),srt=42)
## text(.73,.75,expression(paste(H[0]:~"|"~theta[11]-1/2~"|=|"~theta[12]-1/2~"|")),srt=43)
## dev.off()
