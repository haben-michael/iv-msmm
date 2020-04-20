require(abind)
require(ggplot2)
require(dplyr)
require(xtable)
filelist <- dir()
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]


res <- lapply(filelist,function(file) {
    ## print(file)
    load(file)
    list(reps=reps,n=n,J=J,beta.true=beta[target.idx],beta.causal=beta.causal,sd.sandwich=sd.sandwich,sd.bootstrap=sd.bootstrap,cov.bootstrap=coverages[['bootstrap']]['0.95'],cov.sandwich=coverages[['sandwich']]['0.95'])
})

beta.ipw <- lapply(res,function(l)with(l,cbind(beta.ipw=beta.causal,beta.true=beta.true,sd.sandwich=sd.sandwich,sd.bootstrap=sd.bootstrap,n=n,J=J)))
beta.ipw <- as.data.frame(do.call(rbind,beta.ipw))
group_by(beta.ipw,J,n) %>% summarize(count=n(),bias=median(beta.ipw)-unique(beta.true),sd.mc=sd(beta.ipw),sd.sandwich=mean(sd.sandwich,trim=.0),sd.bootstrap=mean(sd.bootstrap,trim=.0)) %>% as.data.frame()

coverage <- lapply(res,function(l)with(l,cbind(cov.bootstrap=cov.bootstrap,cov.sandwich=cov.sandwich,n=n,J=J,reps=reps)))
coverage <- as.data.frame(do.call(rbind,coverage))
group_by(coverage,J,n) %>% summarize(cov.bootstrap=mean(cov.bootstrap),cov.sandwich=mean(cov.sandwich),reps=sum(reps)) %>% as.data.frame()


tab <- inner_join(
group_by(beta.ipw,J,n) %>% summarize(count=n(),bias=median(beta.ipw)-unique(beta.true),sd.mc=sd(beta.ipw),sd.sandwich=mean(sd.sandwich,trim=.0),sd.bootstrap=mean(sd.bootstrap,trim=.0)) %>% as.data.frame(),
group_by(coverage,J,n) %>% summarize(cov.sandwich=mean(cov.sandwich),cov.bootstrap=mean(cov.bootstrap),reps=sum(reps))%>% as.data.frame()
)  %>% mutate_all(round,2)
tab


## by.J <- do.call(rbind,by.J)
## rownames(by.J) <- NULL
rownames(tab) <- NULL
tab <- subset(tab,J%in%(2:4),select=c(J,n,bias,sd.mc,sd.sandwich,sd.bootstrap,cov.bootstrap,cov.sandwich))
colnames(tab) <- c('T','n','bias','$\\sigma_{mc}$','$\\sigma_{sw}$','$\\sigma_{bs}$','coverage (sw)', 'covarage (bs)')
## sink('010620a.tex')
## print(xtable(tab), sanitize.text.function=function(x){x},floating=TRUE,floating.environment='sidewaystable')
print(xtable(tab), sanitize.text.function=function(x){x})
## sink()
