filelist <- dir()
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
## filelist <- filelist[1:10]
pb.alphas <- c(.05,.15)
ma.alpha <- .05
powers <- lapply(filelist, function(file) {
    ## file <- filelist[1]
    load(file)
    rows <- lapply(pb.alphas, function(pb.alpha) {
        null.idx <- stats['begg.pval',] > pb.alpha
        begg.power <- mean(1-pnorm(abs(stats['ma.stat',null.idx])) < ma.alpha/2)
        null.idx <- stats['egger.pval',] > pb.alpha
        egger.power <- mean(1-pnorm(abs(stats['ma.stat',null.idx])) < ma.alpha/2)
        unconditional.power <- mean(1-pnorm(abs(stats['ma.stat',])) < ma.alpha/2)
        ## if(begg.power>1)browser()
        data.frame(B=B,n=n,pb.alpha=pb.alpha,rZ.idx=rZ.idx,Z.param=Z.param,theta=theta,begg=begg.power,egger=egger.power,unconditional=unconditional.power)
    })
    do.call(rbind,rows)
})
powers <- do.call(rbind,powers)
powers <- powers[order(powers$Z.param),]
## params <-subset(powers,select=c('n','df'))
## powers <- subset(powers,select=c('B','unconditional','egger','begg'))
powers <- by(powers,list(powers$n,powers$pb.alpha,powers$rZ.idx,powers$Z.param,powers$theta),FUN=function(df) {
    wts <- df$B/sum(df$B)
    c(B=sum(df$B),n=unique(df$n),pb.alpha=unique(df$pb.alpha),rZ.idx=unique(df$rZ.idx),Z.param=unique(df$Z.param),theta=unique(df$theta),colSums(wts*subset(df,select=c('begg','egger','unconditional'))))
})
powers <- as.data.frame(do.call(rbind,powers))
powers <- powers[order(powers$n,powers$rZ.idx,powers$Z.param),]


## figure
powers <- subset(powers,theta==0.2 & pb.alpha==.05)

ns <- unique(powers$n)
op <- par(mfrow=c(length(ns),length(unique(powers$rZ.idx))))
for(n0 in ns) {
    for(rZ.idx0 in unique(powers$rZ.idx)) {
        with(list(powers=subset(powers,rZ.idx==rZ.idx0 & n==n0)), {
            ## browser()
            plot(powers$Z.param,powers$unconditional,ylim=range(subset(powers,select=c('begg','egger','unconditional'))),type='l')    
            lines(powers$Z.param,powers$egger,col=2)
            lines(powers$Z.param,powers$begg,col=3)
            ## legend('bottomright',lty=1,col=1:3,legend=c('unconditional','egger','begg'))
        })
    }
}
par(op)



## table

rZ.idxs <- 1:3
Z.params <- list(c(low=4,med=3,high=2.1), c(low=-.1,med=-.5,high=-.9), c(low=.9,med=.5,high=.1))
theta0 <- 0

require(xtable)
require(abind)

for(theta0 in c(0,.2)) {
    out <- subset(powers, theta==theta0 & ((rZ.idx==rZ.idxs[[1]] & Z.param %in% Z.params[[1]]) | (rZ.idx==rZ.idxs[[2]] & Z.param %in% Z.params[[2]]) | (rZ.idx==rZ.idxs[[3]] & Z.param %in% Z.params[[3]])), select=c(n,pb.alpha,rZ.idx,Z.param,begg,egger,unconditional))
    out <- reshape(out,varying=c('begg','egger','unconditional'),v.names='power',timevar='condition',times=c('begg','egger','unconditional'),direction='long',new.row.names=NULL)
    out$id <- NULL
    out$zeta <- out$Z.param
    for(j in rZ.idxs) for(k in 1:length(Z.params[[j]])) out$zeta[out$rZ.idx==j & out$Z.param==Z.params[[j]][k] ]  <- names(Z.params[[j]])[k]
    out$zeta <- factor(out$zeta,levels=c('low','med','high'))
    out <- out[order(out$n,out$pb.alpha,out$rZ.idx,out$condition,out$zeta),]
    dims <- sapply(subset(out,select=rev(c(n,pb.alpha,rZ.idx,condition))),function(x)sort(unique(x)))
    dims <- c(list(zeta=levels(out$zeta)),dims)
    out <- out$power
    dim(out) <- sapply(dims,length)
    dimnames(out) <- dims
    out <- round(ftable(out,row.vars=c(2,4,5),col.vars=c(3,1)),3)
    out <- xtableFtable(out,method='compact')
    sink(paste0('power_',theta0,'.tex'))
    ## addtorow <- list(pos=list(0, 0),
                     ## command=c("& &\\multicolumn{3}{c}{meta-analysis size} \\\\\n"
                               ## paste0('precision distribution &',paste0(ns[ns.idx],collapse='&'), '\\\\\n'))
                               ## )
    ## )
    attr(out,'col.vars')$rZ.idx <- c('t','power','beta')
    ## names(attr(out,'col.vars'))[1] <- '$f_Z$'
    names(attr(out,'col.vars')) <- c('$f_Z$','$\\zeta$')
    names(attr(out,'row.vars'))[names(attr(out,'row.vars'))=='pb.alpha'] <- '$\\alpha_0$'
    print.xtableFtable(out,floating=FALSE,latex.environment=NULL)
    sink()
}

dd






## ## with(powers, egger[rZ.idx==1 & Z.param %in% c(2.1) & theta==0.2])

## rZ.idxs <- 1:3
## Z.params <- list(c(low=4,med=3,high=2.1), c(low=-.1,med=-.5,high=-.9), c(low=.9,med=.5,high=.1))


## require(xtable)
## out <- subset(powers, (rZ.idx==rZ.idxs[[1]] & Z.param %in% Z.params[[1]]) | (rZ.idx==rZ.idxs[[2]] & Z.param %in% Z.params[[2]]) | (rZ.idx==rZ.idxs[[3]] & Z.param %in% Z.params[[3]]), select=c(n,pb.alpha,rZ.idx,Z.param,begg,egger,unconditional))
## out <- reshape(out,varying=c('begg','egger','unconditional'),v.names='power',timevar='condition',times=c('begg','egger','unconditional'),direction='long',new.row.names=NULL)
## out$id <- NULL
## out$zeta <- out$Z.param
## for(j in rZ.idxs) for(k in 1:length(Z.params[[j]])) out$zeta[out$rZ.idx==j & out$Z.param==Z.params[[j]][k] ]  <- names(Z.params[[j]])[k]


## out <- xtabs(power ~ . - Z.param, data=out)
## out <- ftable(out,row.vars=c('condition','pb.alpha','n'),col.vars=c('rZ.idx','zeta'))



## dd
## ## out <- as.matrix(out)


## ## out <- by(out, list(out$n,out$pb.alpha,out$rZ.idx), function(df)as.matrix(df))
## ## attr(out,'class') <- NULL
## ## out2 <- lapply(out, function(m)
## ##     structure(t(m[,c('begg','egger','unconditional')]),dimnames=list(c('begg','egger','unconditional'),m[,'Z.param']))
## ##     )

## ## attr(out2,'dim') <- c(attr(out2,'dim'),3,3)
## ## attr(out2,'dimnames') <- c(attr(out,'dimnames'),list(c('begg','egger','unconditional'),c('low','medium','high')))
## ## ## attr(out2,'dimnames') <- attr(out,'dimnames')    

## ## out <- split(out,with(out,list(n,pb.alpha,rZ.idx)))
## ## out[sapply(out,nrow)==0] <- NULL
## ## ## stopifnot(sum(sapply(out,nrow))==nrow(powers))

## ## dd



## require(xtable)

##                                         # set up data frame
## df <- data.frame(c(replicate(2, c("L1")), replicate(2, c("L2"))),
##                  replicate(4, "b"),
##                  replicate(4, runif(4, 1, 10)) )

##                                         # only needed if first column consists of numbers
## df[[1]] <- as.character(df[[1]])

## rle.lengths <- rle(df[[1]])$lengths
## first <- !duplicated(df[[1]])
## df[[1]][!first] <- ""

##                                         # define appearance of \multirow
## df[[1]][first] <-
##     paste0("\\midrule\\multirow{", rle.lengths, "}{*}{\\textbf{", df[[1]][first], "}}")

## strCaption <- paste0("\\textbf{Table Whatever} This table is just produced with some ",
##                      "random data and does not mean anything. Just to show you how ",
##                      "things work.")

##                                         # set up xtable output
## print(xtable(df, digits = c(0, 0, 0, 3, 1, 0, 6), # first zero "represents" row numbers which we skip later
##              align = "lllrr|rr",  # align and put a vertical line (first "l" again represents column of row numbers)
##              caption = strCaption, label = "testTable"),
##       size = "footnotesize", #Change size; useful for bigger tables "normalsize" "footnotesize"
##       include.rownames = FALSE, #Don't print rownames
##       include.colnames = FALSE, #We create them ourselves
##       caption.placement = "top", #"top", NULL
##       hline.after=NULL, #We don't need hline; we use booktabs
##       floating=TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
##       sanitize.text.function = force, # Important to treat content of first column as latex function
##       add.to.row = list(pos = list(-1,
##                                    2,
##                                    nrow(df)),
##                         command = c(paste("\\toprule \n",  # NEW row
##                                           "\\multicolumn{2}{c}{} & \\multicolumn{2}{c}{\\textbf{colLabel1}} & \\multicolumn{2}{c}{colLabel2} \\\\\n",
##                                           "\\cmidrule(l){3-4} \\cmidrule(l){5-6}\n",
##                                           " & & a1 & a2 & a3 & a4 \\\\\n", # NEW row
##                                           "\\midrule \n"
##                                           ),
##                                     paste("\\cmidrule(l){3-4} \\cmidrule(l){5-6}\n" # we may also use 'pos' and 'command' to add a midrule
##                                           ),
##                                     paste("\\bottomrule \n"  # paste is used as it is more flexible regarding adding lines
##                                           )
##                                     )
##                         )
##             )
