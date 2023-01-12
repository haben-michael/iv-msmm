filelist <- dir()
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
## filelist <- filelist[1]

by.Zs <- lapply(filelist, function(file) {
    ## file <- filelist[1]
    load(file)    
    ## rZ.names <- sapply(rZs, function(s)deparse(s)[2])
    ## do.call(rbind, lapply(1:length(rZs), function(j) cbind(distr=rZ.names[j],do.call(rbind, lapply(c('probs.full','probs.1','probs.2'), function(prob)data.frame(B=B,n=ns,prob=prob,mean=as.numeric(by.rZ[[j]][[prob]])))))))
    by.Z*B
})

Bs <- sapply(filelist, function(file) {
        print(file)
    ## file <- filelist[1]
    load(file)    
    ## rZ.names <- sapply(rZs, function(s)deparse(s)[2])
    ## do.call(rbind, lapply(1:length(rZs), function(j) cbind(distr=rZ.names[j],do.call(rbind, lapply(c('probs.full','probs.1','probs.2'), function(prob)data.frame(B=B,n=ns,prob=prob,mean=as.numeric(by.rZ[[j]][[prob]])))))))
    B
})

by.Z <- Reduce(`+`,by.Zs) / sum(Bs)

by.Z <- apply(by.Z,2,function(col)predict(loess(col ~ ns,span=.75)))

## load(filelist[1])
png('../../ms/comm/neg_bias_students.png')
colors <- rev(gray.colors(n.rZ))
matplot(ns,by.Z,type='l',col=colors,lty=1,xlab='no. of studies',ylab='FPR',xlim=c(1500,5000))
abline(h=.1,lty=2)
legend('bottomright',col=colors,lty=1,legend=formatC(Z.params,format='e',2),title=expression(df~-~2))
dev.off()


