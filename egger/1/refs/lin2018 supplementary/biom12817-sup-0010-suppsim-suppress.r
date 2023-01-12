sim.suppress.pval <- function(mu, tau, s.low, s.upp, pub.rate,
    n, n.resam, seed){
  set.seed(seed)
  pub.rate.sig <- pub.rate[1]
  pub.rate.nonsig <- pub.rate[2]
  y <- s2 <- is.sig <- NULL
  n.gen <- 0
  while(n.gen < n){
    mui.temp <- rnorm(1, mu, tau)
    si.temp <- runif(1, s.low, s.upp)
    y.temp <- rnorm(1, mui.temp, si.temp)
    s2.temp <- si.temp^2
    pval.temp <- 2*(1 - pnorm(abs(y.temp)/si.temp))
    if(pval.temp < 0.05){
      select <- (rbinom(1, 1, pub.rate.sig) == 1)
    }else{
      select <- (rbinom(1, 1, pub.rate.nonsig) == 1)
    }
    if(select){
      y <- c(y, y.temp)
      s2 <- c(s2, s2.temp)
      is.sig <- c(is.sig, pval.temp < 0.05)
      n.gen <- n.gen + 1
    }
  }

  out.metapb <- metapb.sim(y, s2, n.resam = n.resam)
  out <- c("prop.nonsig" = 1 - mean(is.sig),
    "egger.pval" = out.metapb$egger.pval,
    "rank.pval" = out.metapb$rank.pval,
    "tf.pval" = out.metapb$tf.pval,
    "reg.pval" = out.metapb$reg.pval,
    "reg.pval.resam" = out.metapb$reg.pval.resam,
    "skewness.pval" = out.metapb$skewness.pval,
    "skewness.pval.resam" = out.metapb$skewness.pval.resam,
    "combined.pval" = out.metapb$combined.pval,
    "combined.pval.resam" = out.metapb$combined.pval.resam)
  return(out)
}