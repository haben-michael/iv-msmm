filelist <- dir()
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]

res <- lapply(filelist,function(file) {
    load(file)
    data.frame(T=T,n=n,t(beta.hats[,T,])-beta[T])
})
res <- do.call(rbind,res)


biases <- by(res, list(res$T,res$n), function(df)data.frame(bias.assoc=median(df$assoc),bias.sra=median(df$sra),bias.causal=median(df$causal,na.rm=TRUE),bias.oracle=median(df$oracle)))
biases <- cbind(structure(apply(expand.grid(dimnames(biases)),2,as.numeric),dimnames=list(NULL,c('T','n'))), do.call(rbind,biases))



for(t in sort(unique(biases$T))) {
    with(subset(biases,T==t), {
        png(paste0('sim_bias_T',t,'.png'),width=480,height=480)
        plot(bias.assoc ~ n,ylim=range(subset(biases,subset=T==t,select=-c(T,n))),type='l',xlab='sample size',ylab='bias')
        lines(bias.sra ~ n,type='l',lty=2)
        lines(bias.causal ~ n,type='l',lty=3)
        lines(bias.oracle ~ n,type='l',lty=4)
        legend('topleft',lty=1:4,legend=c('associational','SRA','proposed','oracle'),bg='white')
        abline(h=0)
        dev.off()
    })
}

