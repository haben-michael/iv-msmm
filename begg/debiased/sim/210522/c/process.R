
filelist <- dir()
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
## load(filelist[1])

res <- lapply(filelist, function(file) {
    ## file <- filelist[1]
    load(file)    
    rZ.names <- sapply(rZs, function(s)deparse(s)[2])
    do.call(rbind, lapply(1:length(rZs), function(j) cbind(distr=rZ.names[j],do.call(rbind, lapply(c('probs.full','probs.1','symdiff.1','symdiff.2','symdiff.3','symdiff.3a','symdiff.3b'), function(prob)data.frame(B=B,n=ns,prob=prob,mean=as.numeric(by.rZ[[j]][[prob]])))))))
})
skinny <- do.call(rbind,res)



means <- by(skinny, list(skinny$n,skinny$distr,skinny$prob), function(df)data.frame(distr=unique(df$distr),prob=unique(df$prob),n=unique(df$n),B=sum(df$B),mean=sum(df$B*df$mean)/sum(df$B)))
means <- do.call(rbind,means)

means <- subset(means,n!=6) ## for plotting

rZ.names <- sort(unique(means$distr))
rZ.names <- rZ.names[-4] # drop beta(.4,.4)

## png('1.png')
op <- par(mfrow=c(3,4))
for(j in 1:length(rZ.names)) {
    mmm <- subset(means,subset=means$distr==rZ.names[j],select=-distr)
    mmm$mean[mmm$prob %in% c('probs.full','probs.1')] <- mmm$mean[mmm$prob %in% c('probs.full','probs.1')]-1/2
    with(list(df=subset(mmm,subset=prob=='symdiff.1')), plot(df$n,df$mean,ylim=range(mmm$mean),main=rZ.names[j]))
    with(list(df=subset(mmm,subset=prob=='symdiff.2')), points(df$n,df$mean,col=2))
    with(list(df=subset(mmm,subset=prob=='symdiff.3')), points(df$n,df$mean,col=3))
    lll <- merge(subset(mmm,subset=prob=='symdiff.1',select=-c(B,prob)),subset(mmm,subset=prob=='symdiff.2',select=-c(B,prob)),by='n')
    with(lll,lines(n,mean.x+mean.y),main=rZ.names[j],col=4)
abline(h=0,lty=2)
stopifnot(sum(with(mmm,mean[prob=='symdiff.3'] > mean[prob=='symdiff.1'] + mean[prob=='symdiff.2']))==0)
legend('topright',pch=1,col=1:4,legend=c('1','2','3','1+2'))
}
par(op)



## png('1.png')
op <- par(mfrow=c(3,4))
for(j in 1:length(rZ.names)) {
    mmm <- subset(means,subset=means$distr==rZ.names[j],select=-distr)
    mmm$mean[mmm$prob %in% c('probs.full','probs.1')] <- mmm$mean[mmm$prob %in% c('probs.full','probs.1')]-1/2
    with(list(df=subset(mmm,subset=prob=='symdiff.3')), plot(df$n,df$mean,ylim=range(mmm$mean),main=rZ.names[j]))
    with(list(df=subset(mmm,subset=prob=='symdiff.3a')), points(df$n,df$mean,col=2))
    with(list(df=subset(mmm,subset=prob=='symdiff.3b')), points(df$n,df$mean,col=3))
    abline(h=0,lty=2)
    stopifnot(sum(with(mmm,mean[prob=='symdiff.3'] > mean[prob=='symdiff.1'] + mean[prob=='symdiff.2']))==0)
    legend('topright',pch=1,col=1:3,legend=c('3','3a','3b'))
}
par(op)




op <- par(mfrow=c(3,4))
for(j in 1:length(rZ.names)) {
    mmm <- subset(means,subset=means$distr==rZ.names[j],select=-distr)
    mmm$mean[mmm$prob %in% c('probs.full','probs.1')] <- mmm$mean[mmm$prob %in% c('probs.full','probs.1')]-1/2
    ## with(list(df=subset(mmm,subset=prob=='symdiff.3')), plot(df$n,df$mean,ylim=range(mmm$mean),main=rZ.names[j]))
    ## with(list(df=subset(mmm,subset=prob=='symdiff.3a')), points(df$n,df$mean,col=2))
    with(list(df=subset(mmm,subset=prob=='symdiff.3b')), plot(df$n,df$n*df$mean,col=3,main=rZ.names[j]))
    abline(h=0,lty=2)
    stopifnot(sum(with(mmm,mean[prob=='symdiff.3'] > mean[prob=='symdiff.1'] + mean[prob=='symdiff.2']))==0)
    legend('topleft',pch=1,col=1:3,legend=c('3','3a','3b'))
}
par(op)
