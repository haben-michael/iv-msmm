require(abind)
require(xtable)
filelist <- dir()
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]

res <- lapply(filelist,function(file) {
    ## print(file)
    load(file)
     list(reps=reps,n=n,T=T,beta.true=beta[target.idx],beta.causal=beta.causal,sd.sandwich=sd.sandwich,sd.bootstrap=sd.bootstrap,cov.bootstrap=coverages[['bootstrap']],cov.sandwich=coverages[['sandwich']])
})

beta.ipw <- lapply(res,function(lst)with(lst,cbind(beta.hat=beta.causal,beta.true=beta.true,sd.sandwich=sd.sandwich,sd.bootstrap=sd.bootstrap,n=n,T=T)))
beta.ipw <- as.data.frame(do.call(rbind,beta.ipw))
grouped <- by(beta.ipw,list(beta.ipw$T,beta.ipw$n),function(df)c(T=unique(df$T),n=unique(df$n),count=nrow(df),bias=mean(df$beta.hat)-unique(df$beta.true),sd.mc=sd(df$beta.hat),sd.sandwich=mean(df$sd.sandwich),sd.bootstrap=mean(df$sd.bootstrap)))
grouped.bias <- do.call(rbind,grouped)


coverage <- lapply(res,function(lst)with(lst,cbind(cov.bootstrap=cov.bootstrap,cov.sandwich=cov.sandwich,n=n,T=T,reps=reps)))
coverage <- as.data.frame(do.call(rbind,coverage))
grouped <- by(coverage,list(coverage$T,coverage$n), function(df)c(T=unique(df$T),n=unique(df$n),cov.bootstrap=mean(df$cov.bootstrap),cov.sandwich=mean(df$cov.sandwich),reps=sum(df$reps)))
grouped.coverage <- do.call(rbind,grouped)


joined <- merge(grouped.bias,grouped.coverage)

out <- subset(joined,subset= n%in%c(200,400,600,800,1000),select=c(T,n,bias,sd.mc,sd.sandwich,sd.bootstrap,cov.sandwich,cov.bootstrap))
out <- out[order(out$T,out$n),]
out <- round(out,3)
out <- cbind(apply(out[,1:2],2,as.character),apply(out[,3:ncol(out)],2, function(x)sprintf('%.3f',x)))


out.filename <- 'sim.table.probit.tex'
colnames(out) <- c('T','n','bias','$\\sigma_{mc}$','$\\sigma_{sw}$','$\\sigma_{bs}$','coverage (sw)', 'coverage (bs)')
sink(out.filename)
print(xtable(out,display=rep("s",ncol(out)+1),align=rep('c',ncol(out)+1)), sanitize.text.function=function(x){x},include.rownames=FALSE)
sink()


