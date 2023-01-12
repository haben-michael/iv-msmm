## from sherlock
require(parallel)
save.path <- '/mnt/d/data'
filelist <- dir(save.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
## rZ <- function(n)rbeta(n,.2,.2)-1/2
## rZ <- function(n)rnorm(n)
rZs <- list(rnorm, function(n)runif(n,-1/2,1/2), function(n)rbeta(n,.2,.2)-1/2,function(n)rbeta(n,.3,.3)-1/2,function(n)rbeta(n,.4,.4)-1/2,function(n)rbeta(n,.5,.5)-1/2,function(n)rbeta(n,.6,.6)-1/2,function(n)rt(n,df=1),function(n)rt(n,df=1.5),function(n)rt(n,df=2),function(n)rt(n,df=2.25),function(n)rt(n,df=2.5),function(n)rt(n,df=3))

for(j in 1:length(rZs[])) {
    rZ <- rZs[[j]]
    diffs <- mclapply(filelist,mc.cores=detectCores()-2,FUN=function(file) {
        print(file)
        load(paste0(save.path,'/',file))    
        B <- length(by.n[[1]])
        zs <- lapply(ns, function(n)matrix(rZ(n*B),nrow=n))
        evec.pos <- lapply(1:length(ns), function(n) sapply(by.n[[n]],function(eig)eig$vectors[,1]))
        evec.neg <- lapply(1:length(ns), function(n) sapply(by.n[[n]],function(eig)eig$vectors[,2]))
        proj4.pos <- lapply(1:length(ns), function(n)colSums(zs[[n]][1:4,]*evec.pos[[n]][1:4,]))
        proj4.neg <- lapply(1:length(ns), function(n)colSums(zs[[n]][1:4,]*evec.neg[[n]][1:4,]))
        projrest.pos <- lapply(1:length(ns), function(n)colSums(zs[[n]][5:(ns[n]),]*evec.pos[[n]][5:(ns[n]),]))
        projrest.neg <- lapply(1:length(ns), function(n)colSums(zs[[n]][5:(ns[n]),]*evec.neg[[n]][5:(ns[n]),]))
        sq.pos <- sapply(1:length(ns), function(n)c(sq4=mean(proj4.pos[[n]]^2),cross=mean(proj4.pos[[n]]*projrest.pos[[n]]),sqrest=mean(projrest.pos[[n]]^2)))
        sq.neg <- sapply(1:length(ns), function(n)c(sq4=mean(proj4.neg[[n]]^2),cross=mean(proj4.neg[[n]]*projrest.neg[[n]]),sqrest=mean(projrest.neg[[n]]^2)))
        ## structure(sq.neg-sq.pos,colnames=ns)
        colnames(sq.neg) <- colnames(sq.pos) <- ns
        sq.neg-sq.pos
    })
    diffs <- simplify2array(diffs)

    medians <- apply(diffs,1:2,median)
    ns <- colnames(medians)
    png(paste0(j,'.png'))
    plot(ns,medians['sq4',],type='l',ylim=range(medians),main=paste0(deparse(rZ)[2]))
    abline(h=0,lty=2)
    lines(ns,medians['cross',],col=2)
    lines(ns,medians['sqrest',],col=3)
    lines(ns,colSums(medians),col=4,lty=2)
    abline(h=0,lty=2)
    dev.off()
    mean(colSums(medians)<0)
}

## ## from sherlock -- big mean rather than mean of means
## require(parallel)
## save.path <- '/mnt/d/data'
## filelist <- dir(save.path)
## filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]


## ## rZ <- function(n)rbeta(n,.2,.2)-1/2
## ## rZ <- runif(n)-1/2
## ## rZ <- function(n)rnorm(n)
## ## rZ <- function(n)rt(n,df=2.25)
## rZs <- list(rnorm, function(n)runif(n,-1/2,1/2), function(n)rbeta(n,.2,.2)-1/2,function(n)rbeta(n,.3,.3)-1/2,function(n)rbeta(n,.4,.4)-1/2,function(n)rbeta(n,.5,.5)-1/2,function(n)rbeta(n,.6,.6)-1/2,function(n)rt(n,df=1),function(n)rt(n,df=1.5),function(n)rt(n,df=2),function(n)rt(n,df=2.25),function(n)rt(n,df=2.5),function(n)rt(n,df=3))


## for(j in 1:length(rZs)) {
##     rZ <- rZs[[j]]
##     res <- mclapply(filelist,mc.cores=detectCores()-2,FUN=function(file) {
##         ## res <- lapply(filelist,FUN=function(file) {
##         print(file)
##         load(paste0(save.path,'/',file))
##         B <- length(by.n[[1]])
##         zs <- lapply(ns, function(n)matrix(rZ(n*B),nrow=n))
##         evec.pos <- lapply(1:length(ns), function(j) sapply(by.n[[j]],function(eig)eig$vectors[,1]))
##         evec.neg <- lapply(1:length(ns), function(j) sapply(by.n[[j]],function(eig)eig$vectors[,2]))
##         eval.pos <- sapply(1:length(ns), function(j) sapply(by.n[[j]], function(eig)eig$values[1]))
##         eval.neg <- sapply(1:length(ns), function(j) sapply(by.n[[j]], function(eig)eig$values[2]))
##         proj.pos <- sapply(1:length(ns), function(j)colSums(zs[[j]]*evec.pos[[j]]))
##         proj.neg <- sapply(1:length(ns), function(j)colSums(zs[[j]]*evec.neg[[j]]))
##         data.frame(B=B, n=ns, gt=colSums(proj.neg^2*eval.neg + proj.pos^2*eval.pos < 0), gt.vec=colSums(proj.neg^2 - proj.pos^2 > 0))
##     })
##     res <- do.call(rbind,res)

##     png(paste0(j,'.png'))
##     mean.gt <- by(res,res$n,function(df)sum(df$gt)/sum(df$B))
##     plot(names(mean.gt),mean.gt,main=paste0(deparse(rZ)[2]))
##     abline(h=.5,lty=2)
##     mean.gt.vec <- by(res,res$n,function(df)sum(df$gt.vec)/sum(df$B))
##     points(names(mean.gt.vec),mean.gt.vec,col=2)
##     dev.off()

## }
