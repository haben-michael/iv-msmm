library(parallel)
library(metafor)

setwd("...") # set your working directory
source("metapb.sim.R")
source("sim.suppress.smallstudies.R")
source("tnf.R")

n <- 10 # 10, 30, or 50

n.iter <- 10000
n.resam <- 1000

mu <- 1
s.low <- 1
s.upp <- 4
s.cutoff <- 1.5

set.seed(12345)
iter.seeds <- sample(1000000, n.iter)

taus <- c(0, 1, 4)
pub.rates <- list(c(1, 1), c(1, 0.2), c(1, 0.1), c(1, 0))

for(id.tau in 1:length(taus)){
  tau <- taus[id.tau]
  for(id.pub.rate in 1:length(pub.rates)){
    pub.rate <- pub.rates[[id.pub.rate]]
    pvals.list <- mclapply(1:n.iter, function(i){
        sim.suppress.smallstudies(mu, tau, s.low, s.upp, s.cutoff, pub.rate,
          n, n.resam, iter.seeds[i])
      }, mc.cores = 20)
    pvals <- unlist(pvals.list)
    props.nonsig <- pvals[names(pvals) == "prop.nonsig"]
    props.nonsig.smry <- c(mean(props.nonsig), sd(props.nonsig),
      quantile(props.nonsig, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
    props.nonsig.smry <- t(props.nonsig.smry)
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

    write.table(props.nonsig.smry, file =
      paste("suppress_smallstudies_n_", n, "_tau_", tau,
      "_pub.rate.sig_", pub.rate[1], "_pub.rate.nonsig_", pub.rate[2],
      "_props.nonsig.txt", sep = ""), row.names = FALSE,
      col.names = c("mean", "se", "lower", "1st", "median", "3rd", "upper"))
    write.table(rslts, file = paste("suppress_smallstudies_n_", n, "_tau_", tau,
      "_pub.rate.sig_", pub.rate[1], "_pub.rate.nonsig_", pub.rate[2],
      "_results.txt", sep = ""), row.names = TRUE, col.names = FALSE)
  }
}
