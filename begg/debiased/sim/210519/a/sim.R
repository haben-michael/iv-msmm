args <- commandArgs(trailingOnly=TRUE)
B <- as.numeric(args[1])
## B <- 1e2
## print(B)

ns <- round(seq(1e2,2e2,len=10))
## by.n <- mclapply(ns, mc.cores=detectCores()-2,FUN=function(n) {
by.n <- lapply(ns,FUN=function(n) {
    ## cat('.')
    proj <- lapply(1:B, function(jj){
        ## s <- runif(n,0,1)
        ## s <- rgamma(n,shape=.54)
        s <- rexp(n)
        ## z <- runif(n,-1,1)
        ## z <- rbeta(n,.5,.5)-1/2
        ## z <- rt(n,df=40)
        z <- rnorm(n)
        Svec1 <- c(c(1/(s[1] - s[2]) - s[1]/sum(s^2), -1/(s[1] - s[2]) - s[2]/sum(s^2), -s[3]/sum(s^2), -s[4]/sum(s^2)), -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        Svec2 <- c(c(-s[1]/sum(s^2), -s[2]/sum(s^2), 1/(s[3] - s[4]) - s[3]/sum(s^2), -1/(s[3] - s[4]) - s[4]/sum(s^2)),  -s[4+seq_len(max(n+1-5,0))]/sum(s^2))
        SSmat <- (Svec1%*%t(Svec2) + t(Svec1%*%t(Svec2)))/2
        eig <- eigen(SSmat)
        eig$values <- eig$values[c(1,n)]
        eig$vectors <- eig$vectors[,c(1,n)]
        eig
        ## c(pos4=eig$vec[1:4,1]%*%z[1:4],neg4=eig$vec[1:4,n]%*%z[1:4],posn=eig$vec[5:n,1]%*%z[5:n],negn=eig$vec[5:n,n]%*%z[5:n])
    })
    ## boxplot(t(abs(proj)))
    ## possq <- cbind(pos4sq=proj['pos4',]^2,poscross=proj['pos4',]*proj['posn',],posnsq=proj['posn',]^2)
    ## negsq <- cbind(neg4sq=proj['neg4',]^2,negcross=proj['neg4',]*proj['negn',],negnsq=proj['negn',]^2)
    ## c(colMeans(possq),colMeans(negsq))
})
## by.n <- simplify2array(by.n)
## plot(ns,by.n['neg4sq',]-by.n['pos4sq',],type='l')
## lines(ns,by.n['negcross',]-by.n['poscross',],col=2)
## lines(ns,by.n['negnsq',]-by.n['posnsq',],col=3)
## lines(ns,by.n['neg4sq',]+by.n['negcross',]+by.n['negnsq',] - (by.n['pos4sq',]+by.n['poscross',]+by.n['posnsq',]),col=4,lty=2)
## abline(h=0,lty=2)
## mean(by.n['pos4sq',]<by.n['neg4sq',])
## mean(by.n['neg4sq',]+by.n['negcross',]+by.n['negnsq',] > by.n['pos4sq',]+by.n['poscross',]+by.n['posnsq',])
filename <- paste0('/scratch/users/habnice/save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(ns,by.n,file=filename)
