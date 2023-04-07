filelist <- dir()
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
alpha <- .05
## file <- filelist[1]
coverage <- lapply(filelist, function(file) {
    ## print(file)
    load(file)
    if(isTRUE(is.na(z.stats)))  return(NA)
    n.covered <- rowSums(abs(z.stats)<qnorm(1-alpha/2))
    max.dist <- apply(z.stats,1,function(knots) {
        knots <- sort(knots)
        distances <- pnorm(knots) - 1:length(knots)/length(knots)
        max(abs(distances))
    })
    data.frame( p.full=p.full,p.red=p.red,adjust.size=round(adjust.size,3),n=n,estimator=names(n.covered),reps=reps, n.covered=n.covered)
})
coverage <- do.call(rbind,coverage)
## print(sum(is.na(coverage)))

coverage <- aggregate( cbind(reps,n.covered) ~ ., sum, data=coverage)
rownames(coverage) <- NULL
coverage <- within(coverage, fpr <- n.covered/reps)

## figure
## png('sim_logit.png')
out <- subset(coverage,select=c(estimator,n,adjust.size,fpr))
by.adj.size <- lapply( split(out,out$adjust.size), function(df) {
    out <- reshape(subset(df,select=-adjust.size),direction='wide',timevar='estimator',idvar='n')
    colnames(out) <- sub('out.','',colnames(out))
    out
})
out <- by.adj.size[[1]]
matplot(out$n,subset(out,select=-n),pch=1,type='l',lty=1:3,col=1,xlab='n',ylab='coverage')
legend('bottomright',legend=substr(colnames(out),5,nchar(colnames(out)))[-1],lty=1:3)
abline(h=1-alpha)
## dev.off()

## table
require(xtable)
ns <- c(1e3,2e3,3e3)
out <- subset(coverage, n%in%ns, select=c(estimator,n,adjust.size,fpr))
out <- xtabs(fpr ~ ., out)
out <- ftable(aperm(out,c(3,1,2)))
out <- round(out,3)
sink('table_lda.tex')
xtableFtable(out)
sink()





