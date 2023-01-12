library(parallel)
library(metafor)

setwd("...") # set your working directory
source("metapb.sim.R")
source("sim.suppress.effectsize.R")
source("tnf.R")

n <- 10 # 10, 30, or 50

n.iter <- 10000
n.resam <- 1000

mu <- 1
s.low <- 1
s.upp <- 4

set.seed(12345)
iter.seeds <- sample(1000000, n.iter)

taus <- c(0, 1, 4)
ms <- c(0, floor(n/3), floor(2*n/3))

for(id.tau in 1:length(taus)){
  tau <- taus[id.tau]
  for(id.m in 1:length(ms)){
    m <- ms[id.m]
    pvals.list <- mclapply(1:n.iter, function(i){
        sim.suppress.effectsize(mu, tau, s.low, s.upp, m,
          n, n.resam, iter.seeds[i])
      }, mc.cores = 20)
    pvals <- unlist(pvals.list)
    test.names <- c("egger.pval", "rank.pval", "tf.pval",
      "reg.pval", "reg.pval.resam", "skewness.pval", "skewness.pval.resam",
      "combined.pval", "combined.pval.resam")
    rslts <- numeric(length(test.names))
    names(rslts) <- c("egger", "rank", "tf",
      "reg", "reg.resam", "skewness", "skewness.resam",
      "combined", "combined.resam")
    for(i in 1:length(test.names)){
      pvals.temp <- pvals[names(pvals) == test.names[i]]
      rslts[i] <- mean(pvals.temp < 0.1)
    }

    write.table(rslts, file = paste("suppress_effectsize_n_", n, "_tau_", tau,
      "_m_", m, "_results.txt", sep = ""), row.names = TRUE, col.names = FALSE)
  }
}
