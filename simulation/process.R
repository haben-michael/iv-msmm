require(abind)
require(ggplot2)
require(dplyr)
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

coverage <- lapply(res,function(l)with(l,cbind(cov.bootstrap=cov.bootstrap,cov.sandwich=cov.sandwich,n=n,J=J)))
coverage <- as.data.frame(do.call(rbind,coverage))
group_by(coverage,J,n) %>% summarize(cov.bootstrap=mean(cov.bootstrap),cov.sandwich=mean(cov.sandwich),batches=n()) %>% as.data.frame()


