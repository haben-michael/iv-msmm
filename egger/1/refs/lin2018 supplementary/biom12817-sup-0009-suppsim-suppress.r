sim.suppress.effectsize <- function(mu, tau, s.low, s.upp, m, n, n.resam, seed){
  set.seed(seed)
  mui <- rnorm(n + m, mu, tau)
  si <- runif(n + m, s.low, s.upp)
  y <- rnorm(n + m, mui, si)
  s2 <- si^2
  if(m > 0){
    rks <- rank(y, ties.method = "first")
    y <- y[rks > m]
    s2 <- s2[rks > m]
  }

  out.metapb <- metapb.sim(y, s2, n.resam = n.resam)
  out <- c("egger.pval" = out.metapb$egger.pval,
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