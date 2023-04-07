filelist <- dir()
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]

alpha <- .05
## file <- filelist[1] #!!!!!
coverage <- lapply(filelist, function(file) {
    ## print(file)
    load(file)
    if(isTRUE(is.na(terms)))  return(NA)
    terms <- as.data.frame(t(terms))
    obs <- terms$obs
    true <- unname(unique(terms$true))
    stopifnot(length(true)==1)
    ## stopifnot(abs(unique(true)-auc)<1e-7)
    var.hats <- as.matrix(subset(terms,select=-c(obs,true)))
    test.stats <- (obs - true) / sqrt(var.hats)
    n.covered <- colSums(abs(test.stats)<qnorm(1-alpha/2))
    length.sums <- colSums(2*sqrt(var.hats)*qnorm(1-alpha/2))
    z.stats <- obs / sqrt(var.hats)
    data.frame( p.full=p.full,p.red=p.red,auc=unique(true),pi.0=pi.0,n=n,estimator=names(n.covered),reps=reps, n.covered=n.covered,length.sums=length.sums)
    ## data.frame( p.full=p.full,p.red=p.red,adjust.size=round(adjust.size,3),n=n, estimator=names(n.covered),reps=reps, n.covered=n.covered)
})
coverage <- do.call(rbind,coverage)
## print(sum(is.na(coverage)))
coverage <- aggregate( cbind(reps,n.covered,length.sums) ~ ., sum, data=coverage)
rownames(coverage) <- NULL
coverage <- within(coverage, {
    fpr <- n.covered/reps
    mean.length <- length.sums/reps
})





## figure
## png('sim_logit.png')
out <- subset(coverage,select=c(estimator,n,pi.0,fpr))
by.pi.0 <- lapply( split(out,out$pi.0), function(df) {
    out <- reshape(subset(df,select=-pi.0),direction='wide',timevar='estimator',idvar='n')
    colnames(out) <- sub('out.','',colnames(out))
    out
})
out <- by.pi.0[[2]]
matplot(out$n,subset(out,select=-n),pch=1,type='l',lty=1:4,col=1,xlab='n',ylab='coverage')
legend('bottomright',legend=substr(colnames(out),5,nchar(colnames(out)))[-1],lty=1:4)
abline(h=1-alpha)
## dev.off()




## z.stats <- lapply(filelist, function(file) {
##     ## print(file)
##     load(file)
##     if(isTRUE(is.na(terms)))  return(NA)
##     terms <- as.data.frame(t(terms))
##     obs <- unlist(subset(terms,select=obs))
##     ## true <- unlist(subset(terms,select=true))
##     var.hats <- as.matrix(subset(terms,select=-c(obs,true)))
##     ## print(true)
##     ## test.stats <- (obs - true) / sqrt(var.hats)
##     ## n.covered <- colSums(abs(test.stats)<qnorm(1-alpha/2))
##     ## obs / sqrt(var.hats)
##     data.frame( p.full=p.full,p.red=p.red,adjust.size=round(adjust.size,3),n=n, (obs) / sqrt(var.hats) )
## })
## z.stats <- do.call(rbind,z.stats)
## ## print(sum(is.na(coverage)))
## std.error <- aggregate( cbind(delong,oracle,proposed) ~ ., sd, data=z.stats)
## std.error






## ## figure
## ## png('sim_logit.png')
## out <- subset(coverage,select=c(estimator,n,pi.0,fpr))
## by.adj.size <- lapply( split(out,out$adjust.size), function(df) {
##     out <- reshape(subset(df,select=-adjust.size),direction='wide',timevar='estimator',idvar='n')
##     colnames(out) <- sub('out.','',colnames(out))
##     out
## })
## out <- by.adj.size[[1]]
## matplot(out$n,subset(out,select=-n),pch=1,type='l',lty=1:3,col=1,xlab='n',ylab='coverage')
## legend('bottomright',legend=substr(colnames(out),5,nchar(colnames(out)))[-1],lty=1:3)
## abline(h=1-alpha)
## ## dev.off()

## table
require(xtable)
ns <- c(1e3,5e3,8e3)
fpr <- subset(coverage, n%in%ns, select=c(estimator,n,pi.0,fpr))
fpr <- xtabs(fpr ~ ., fpr)
fpr <- ftable(aperm(fpr,c(3,1,2)))
fpr <- round(fpr,3)
lengths <- subset(coverage, n%in%ns, select=c(estimator,n,pi.0,mean.length))
lengths <- xtabs(mean.length ~ ., lengths)
lengths <- ftable(aperm(lengths,c(3,1,2)))
lengths <- round(lengths,3)
out <- fpr
for(i in 1:ncol(fpr)) out[,i] <- paste0(format(fpr[,i],nsmall=3),' (',format(lengths[,i],nsmall=3),')')
names(attr(out,'row.vars'))[names(attr(out,'row.vars'))=='pi.0'] <- 'imbalance'
## save.file <- 'table_lda_auc.tex'
sink(save.file)
xtableFtable(out)
sink()
lines <- scan(save.file,what='',sep='\n')
start.idx <- grep('begin\\{tabular\\}',lines)
end.idx <- grep('end\\{tabular\\}',lines)
writeLines(lines[start.idx:end.idx],save.file)





