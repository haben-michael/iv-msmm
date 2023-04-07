p.val.mc <- function(p,p.hat,T,reps,L=(0:(m-1))/(m-1)) {
    ## browser()
    ## if(sum(p)!=1)browser()
    if(!isTRUE(all.equal(sum(p),1)))return(NA)
    p.star <- rmultinom(n.ref.samples,n,prob=p) / n
    ## stopifnot(!(1 %in% p.star))
    ## if(1 %in% p.star)browser()
    T.star <- apply(p.star,2,T,p=p,L=L)
    T.hat <- as.numeric(T(p.hat,p,L))
    mean(T.star >= T.hat)
}


rsimplex.hp <- function(n,c=(0:3)/3,shape=1) {
    m <- length(c)
    if(m>1) {
        u.min <- (1-c[m])/sum(c[-m])
        ## if(u.min<0 || u.min>1)return(NA)
        if(u.min>1)return(NA)
        u.min <- max(0,u.min)
    } else {u.min <- 0}
    ## while(TRUE) {
    ##     u <- sort(runif(m))
    ##     u <- u / c(u%*%c)
    ##     if(u[1]<u.min && u[m]<1) {           
    ##         print(u)
    ##         break
    ##     }
    ## }
    ## browser()
    x <- replicate(n, {
        ## u <- sort(runif(m,u.min,1))
        u <- 1-rbeta(m,shape,shape)*(1-u.min)
        u <- sort(u)
        u / c(u%*%c)
    })## ,simplify=FALSE)
}


rank.to.unit <- function(x.bar,m=nrow(x.bar)) {
    M <- diag(m)
    M[row(M)==col(M)+1] <- -1
    M%*%x.bar
}


## sampling strategy in #3. not good since rsimplex.hp wasn't sampling
## all the points.
sample.p <- function(theta,c,p.resolution,shape=1) {
    m <- length(c)
    d <- (c[1:(m-1)]-c[m])/(theta-c[m])
    d.bar <- c(-diff(d),d[m-1])
    p.by.theta.bar <- rsimplex.hp(n=p.resolution,c=d.bar,shape=shape)
    if(m==2)p.by.theta.bar <- matrix(p.by.theta.bar,nrow=1)
    p.by.theta <- rank.to.unit(p.by.theta.bar)
    p.by.theta <- rbind(p.by.theta,1-colSums(p.by.theta))
}

## sampling strategy using orthog complement
p.sampler.1 <- function(n,c,shape=1) {
    d <- c - 1
    d.p <- d[d>=0]; d.m <- -d[d<0]
    if(length(d.p)*length(d.m)==0)return(NA)
    if(anyNA(c))return(NA)
    ## if(length(d.p)+length(d.m)!=length(d))return(NA)
    us <- replicate(n, {
        u.p <- rbeta(length(d.p),shape,shape)
        u.m <- rbeta(length(d.m),shape,shape)
        u.m <- u.m * c(u.p%*%d.p / (u.m%*%d.m))
        ## c <- c(-c.m,c.p)
        u <- c(u.m,u.p)
    })
    ## us <- t(us)
    us <- t(us) / rowSums(t(us))
}

p.sampler.2 <- function(n,c,tol=1e-6) {
    ## browser()
    ## print(c)
    if(anyNA(c))return(NA)
    m <- length(c)
    A <- rbind(c,-diag(m),rep(1,m))
    b <- c(1,rep(0,m),1)
    vv <- vertexenum::enumerate.vertices(A=A,b=b,warn=TRUE)
    ## idx <- vv%*%c==1 & vv%*%rep(1,m)==1
    idx <- which(abs(vv%*%c-1 + rowSums(vv)-1) < tol)
    if(length(idx)<1)return(NA)
    vv <- vv[idx,,drop=FALSE]
    ## if(is.null(dim(vv[idx,])))browser() # need weaker tol
    ## w <- matrix(rexp(n*nrow(vv[idx,])),ncol=nrow(vv[idx,]))
    w <- matrix(rexp(n*nrow(vv)),ncol=nrow(vv))
    w <- w/rowSums(w)
    ## u <- w%*%vv[idx,]
    u <- w%*%vv
}

test.coverage.asy <- function(theta.try,p.hat,c,n,alpha=.05) {
    ## browser()
    var.mle.hat <- (diag(as.numeric(p.hat)) - p.hat%*%t(p.hat)) / n
    z.stat <- (c%*%p.hat - theta.try) / sqrt(t(c)%*%var.mle.hat%*%c)
    abs(z.stat) <= qnorm(1-alpha/2)
}



test.stat.1 <- function(p.obs,p,c) {
    ## browser()
    ## cat('.')
    var.c <- t(c)%*%(diag(as.numeric(p.obs)) - p.obs%*%t(p.obs))%*%c
    abs(c%*%(p.obs-p)) / sqrt(var.c)
}


test.stat.2 <- function(p.obs,p,c,lambda=rep(.1,length(p))) {
    ## browser()
    test.stat.1.part <- test.stat.1(p.obs,p,c)
    addl <- lambda*abs(p.obs-p) / sqrt(p.obs*(1-p.obs))
    test.stat.1.part + sum(addl)
}

p.val.mc <- function(p,p.obs,n,n.ref.samples,c,T) {
    ## browser()
    ## print(p)
    if(!isTRUE(abs(sum(p)-1)<1e-8))return(NA)
    ## if(anyNA(p))return(NA)
    p.star <- rmultinom(n.ref.samples,n,prob=p) / n
    T.star <- apply(p.star,2,T,p=p,c=c)
    T.hat <- as.numeric(T(p.obs,p,c))
    mean(T.star >= T.hat)
}

test.coverage.exact <- function(theta.try,p.hat,c,n,alpha=.05,n.ref.samples=1e2,test.stat,p.sampler,theta.resolution,p.resolution,shape=1) {    
   ## browser()
    ## mc <- function(p) {
    ##     p.star <- rmultinom(n.ref.samples,n,prob=p) / n
    ##     T.star <- apply(p.star,2,T,p=p,c=c)
    ##     T.hat <- as.numeric(T(p.hat,p,c))
    ##     mean(T.star >= T.hat)
    ## }

    ## theta.grid <- seq(min(c)+.001,max(c)-.001,len=theta.resolution) # clean
    ## p.by.theta <- lapply(theta.grid, function(theta) {
    ##     sample.p(n=p.resolution,c=c/theta,shape=shape)
    ## })
    
    ## p.val.bins <- lapply(p.by.theta, function(p){
    ##     if(isTRUE(is.na(p)))return(NA)
    ##     ## apply(p,2,mc)#p.val.mc,p.hat=p.hat,T=T,reps=n.ref.samples)
    ##     apply(p,2,p.val.mc,p.obs=p.hat,n=n,n.ref.samples=n.ref.samples,c=c,T=T)
    ##     })
    ## p.val <- sapply(p.val.bins, max)
    ## idx.max <- which.max(theta.grid>=theta.try)
    ## idx.min <- rev(which(theta.grid<=theta.try))[1]
    ## p.val.try <- mean(p.val[c(idx.max,idx.min)],na.rm=TRUE)
    ## p.val.try > alpha

    ## just the p : c^t p = theta.true
    p.by.theta <- p.sampler(n=p.resolution,c=c/c(theta.try))
    ## print(p.by.theta)
    ## p.vals.simplex <- apply(p.by.theta,2,mc)
    ## print(c/theta.try)
    ## print(p.by.theta)
    p.vals.simplex <- apply(p.by.theta,1,p.val.mc,p.obs=p.hat,n=n,n.ref.samples=n.ref.samples,c=c,T=test.stat)
    max(p.vals.simplex) > alpha
    
    ## ## just at p.true
    ## p.val.mc(p.true,p.hat=p.hat,T=T,reps=n.ref.samples)
}


## collect mc params in a list
multinom.linear <- function(p.obs,c,n,theta,n.ref.samples=1e2,test.stat,p.sampler,theta.resolution=NULL,p.resolution=NULL) {
    ## browser()
    ## theta.grid <- seq(min(c)+.001,max(c)-.001,len=theta.resolution) # clean
    p.by.theta <- lapply(theta, function(theta.i) {
        p.sampler(n=p.resolution,c=c/theta.i)
    })
    ## browser()
    ## p.by.theta[[1]] <- matrix(p.true,nrow=1) #!!!!! for dbg at truth
    
    p.val.bins <- lapply(p.by.theta, function(p){
        if(isTRUE(anyNA(p)))return(NA)
        ## apply(p,2,mc)#p.val.mc,p.hat=p.hat,T=T,reps=n.ref.samples)
        apply(p,1,p.val.mc,p.obs=p.obs,n=n,n.ref.samples=n.ref.samples,c=c,T=test.stat)
    })

    max.pvals <- sapply(p.val.bins,max)

    out <- structure(list(theta=theta,p=p.by.theta,p.val.p=p.val.bins,p.val.theta=max.pvals,data.name=paste(deparse(substitute(p.obs))),p.obs=p.obs),class='multinom.linear')
    return(out)
    ## p.val <- sapply(p.val.bins, max)
    ## idx.max <- which.max(theta.grid>=theta.try)
    ## idx.min <- rev(which(theta.grid<=theta.try))[1]
    ## p.val.try <- mean(p.val[c(idx.max,idx.min)],na.rm=TRUE)
    ## p.val.try > alpha

    ## ## just the p : c^t p = theta.true
    ## p.by.theta <- p.sampler(n=p.resolution,c=c/c(theta.try))
    ## ## print(p.by.theta)
    ## ## p.vals.simplex <- apply(p.by.theta,2,mc)
    ## ## print(c/theta.try)
    ## ## print(p.by.theta)
    ## p.vals.simplex <- apply(p.by.theta,1,p.val.mc,p.obs=p.hat,n=n,n.ref.samples=n.ref.samples,c=c,T=test.stat)
    ## max(p.vals.simplex) > alpha
}


plot.multinom.linear <- function(test.output,alpha=.05){
    theta <- test.output$theta
    p <- test.output$p
    p.val.theta <- test.output$p.val.theta
    plot(theta,p.val.theta,ylim=c(0,1),xlab='theta',ylab='p-value')
    abline(h=alpha)
    }


plot3d.multinom.linear <- function(test.output,alpha=.05,palette=heat.colors(10)){
    theta <- test.output$theta
    p <- test.output$p
    p.val.p <- test.output$p.val.p
    p <- do.call(rbind,p)
    p.val.p <- do.call(base:::c,p.val.p)

    n.colors <- length(palette)
    p.val.binned <- cut(p.val.p,breaks=n.colors)
    col <- palette[p.val.binned]
    rgl::clear3d()
    rgl::plot3d(p,col=col)

    rgl::legend3d('topright',col=palette,pch=1,legend=levels(p.val.binned))
}


## instantiates class htest excpt for "statistic" entry, since
## multinom.linear doesnt currently expose the test statistics, only
## the p-values; spec:https://search.r-project.org/CRAN/refmans/EnvStats/html/htest.object.html
test.multinom.linear <- function(test.output,theta.null) {
    ## theta <- test.output$theta
    ## p <- test.output$p
    ## p.val.theta <- test.output$p.val.theta
    idx <- which.min(abs(test.output$theta-c(theta.null)))
    p.val.null <- test.output$p.val.theta[idx]
    
    out <- with(test.output,
                list(null.value=c(theta=theta.null),alternative='two-sided',method='Monte carlo exact test',estimate=c(p=as.numeric(p.obs)),data.name=data.name,statistic=NA,parameters=c(),p.value=p.val.null)
                )
    class(out) <- 'htest'
    return(out)
}
