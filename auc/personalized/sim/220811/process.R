filelist <- dir()
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]

load(filelist[1])
reject.total <- as.list(numeric(B^2))
dim(reject.total) <- c(B,B)
for(i in 1:B) for(j in 1:B) reject.total[[i,j]] <- logical()

for(file in filelist) {
    load(file)
    for(i in 1:B) for(j in 1:B) reject.total[[i,j]] <- c(reject.total[[i,j]],reject[[i,j]])
    }

print(length(reject.total[[100,100]]))

power <- matrix(NA,B,B)
for(i in 1:B) for(j in 1:B) power[i,j] <- mean(reject.total[[i,j]])
image(theta.11,theta.12,power)
abline(a=0,b=1)#;abline(a=1,b=-1)

## png('211120.png')
plot(0,type='n',xlim=c(1/2,1),ylim=c(1/2,1),xlab=expression(theta[11]),ylab=expression(theta[12]),xaxs='i',yaxs='i')
contour(theta.11,theta.12,power,add=TRUE,nlevels=10,method='edge')
## contourplot(power)
abline(a=0,b=1)
text(.73,.75,expression(paste(H[0]:theta[11]~"="~theta[12])),srt=42)
## text(.73,.75,expression(paste(H[0]:~"|"~theta[11]-1/2~"|=|"~theta[12]-1/2~"|")),srt=43)
## dev.off()


## a = melt(power)
## b = loess(value ~ Var1+Var2, data=a, span=.08)
## smooth = predict(b, newdata=data.frame(Var1=a$Var1, Var2=a$Var2))
## smooth = array(smooth, dim=dim(power))
## contour(theta.11,theta.12,t(smooth))
## ## contourplot(smooth)
