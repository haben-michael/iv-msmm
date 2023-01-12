metapb.sim <- function(y, s2, model = "RE", n.resam = 1000){
  if(length(y) != length(s2) | any(s2 < 0)) stop("error in the input data.")
  n <- length(y)
  w <- 1/s2
  mu.bar <- sum(w*y)/sum(w)
  Q <- sum(w*(y - mu.bar)^2)
  p.Q <- 1 - pchisq(Q, df = n - 1)
  I2 <- (Q - n + 1)/Q
  tau2 <- (Q - n + 1)/(sum(w) - sum(w^2)/sum(w))
  tau2 <- max(c(0, tau2))

  if(missing(model)){
    if(p.Q >= 0.05){
      model <- "FE"
    }else{
      model <- "RE"
    }
  }

  if(!is.element(model, c("FE", "RE"))){
    stop("wrong input for the argument model.")
  }

  y.egger <- y/sqrt(s2)
  x.egger <- 1/sqrt(s2)
  reg.egger <- summary(lm(y.egger ~ x.egger))$coef
  egger.stat <- as.numeric(reg.egger["(Intercept)", "Estimate"])
  egger.pval <- as.numeric(reg.egger["(Intercept)", "Pr(>|t|)"])
  res <- rma(yi = y, vi = s2, method = ifelse(model == "FE", "FE", "DL"))
  rank <- ranktest(res)
  rank.stat <- as.numeric(rank$tau)
  rank.pval <- as.numeric(rank$pval)
  tf <- tnf(res, estimator = "R0")
  tf.stat <- tf$k0
  tf.pval <- tf$p.k0

  if(model == "FE"){
    y.reg <- y/sqrt(s2)
    x.reg <- 1/sqrt(s2)
    mu.est <- sum(y/s2)/sum(1/s2)
  }else{
    y.reg <- y/sqrt(s2 + tau2)
    x.reg <- 1/sqrt(s2 + tau2)
    mu.est <- sum(y/(s2 + tau2))/sum(1/(s2 + tau2))
  }
  reg <- lm(y.reg ~ x.reg)
  reg.coef <- summary(reg)$coef
  reg.int <- reg.coef["(Intercept)", "Estimate"]
  reg.int.se <- reg.coef["(Intercept)", "Std. Error"]
  reg.int.ci.low <- reg.coef["(Intercept)", "Estimate"] -
    qt(0.975, df = n - 2)*reg.coef["(Intercept)", "Std. Error"]
  reg.int.ci.upp <- reg.coef["(Intercept)", "Estimate"] +
    qt(0.975, df = n - 2)*reg.coef["(Intercept)", "Std. Error"]
  reg.int.ci <- c("2.5%" = reg.int.ci.low, "97.5%" = reg.int.ci.upp)
  reg.pval <- reg.coef["(Intercept)", "Pr(>|t|)"]

  std.dev <- as.numeric(summary(reg)$residuals)
  cm2 <- var(std.dev)
  cm3 <- mean((std.dev - mean(std.dev))^3)
  cm4 <- mean((std.dev - mean(std.dev))^4)
  cm5 <- mean((std.dev - mean(std.dev))^5)
  cm6 <- mean((std.dev - mean(std.dev))^6)
  skewness <- cm3/(cm2^(1.5))

  var.skew <- 9 + 35/4*cm3^2/cm2^3 -6*cm4/cm2^2 + cm6/cm2^3 +
    9/4*cm3^2*cm4/cm2^5 - 3*cm3*cm5/cm2^4
  skewness.ci.low <- skewness - 1.96*sqrt(var.skew/n)
  skewness.ci.upp <- skewness + 1.96*sqrt(var.skew/n)
  skewness.ci <- c("2.5%" = skewness.ci.low, "97.5%" = skewness.ci.upp)
  var0.skew <- 6
  skewness.pval <- 2*(1 - pnorm(sqrt(n/var0.skew)*abs(skewness)))

  combined.pval <- 1 - (1 - min(c(reg.pval, skewness.pval)))^2

  reg.int.star <- skewness.star <- numeric(n.resam)
  for(i in 1:n.resam){
    if(model == "FE"){
      y.star <- rnorm(n = n, mean = mu.est, sd = sqrt(s2))
      y.reg.star <- y.star/sqrt(s2)
      x.reg.star <- 1/sqrt(s2)
    }else{
      y.star <- rnorm(n = n, mean = mu.est, sd = sqrt(s2 + tau2))
      w.star <- 1/s2
      mu.bar.star <- sum(w.star*y.star)/sum(w.star)
      Q.star <- sum(w.star*(y.star - mu.bar.star)^2)
      tau2.star <- (Q.star - n + 1)/(sum(w.star) - sum(w.star^2)/sum(w.star))
      tau2.star <- max(c(0, tau2.star))
      y.reg.star <- y.star/sqrt(s2 + tau2.star)
      x.reg.star <- 1/sqrt(s2 + tau2.star)
    }

    reg.star <- lm(y.reg.star ~ x.reg.star)
    reg.coef.star <- summary(reg.star)$coef
    reg.int.star[i] <- reg.coef.star["(Intercept)", "Estimate"]
    std.dev.star <- as.numeric(summary(reg.star)$residuals)
    cm2.star <- var(std.dev.star)
    cm3.star <- mean((std.dev.star - mean(std.dev.star))^3)
    cm4.star <- mean((std.dev.star - mean(std.dev.star))^4)
    cm5.star <- mean((std.dev.star - mean(std.dev.star))^5)
    cm6.star <- mean((std.dev.star - mean(std.dev.star))^6)
    skewness.star[i] <- cm3.star/(cm2.star^(1.5))
  }

  reg.pval.resam <- (sum(abs(reg.int.star) >= abs(reg.int)) + 1)/(n.resam + 1)
  skewness.pval.resam <-
    (sum(abs(skewness.star) >= abs(skewness)) + 1)/(n.resam + 1)
  combined.pval.resam <-
    1 - (1 - min(c(reg.pval.resam, skewness.pval.resam)))^2

  out <- NULL
  out$n <- n
  out$p.Q <- p.Q
  out$I2 <- I2
  out$tau2 <- tau2
  out$model <- model
  out$std.dev <- std.dev
  out$egger.pval <- egger.pval
  out$rank.pval <- rank.pval
  out$tf.pval <- tf.pval
  out$reg.int <- reg.int
  out$reg.int.ci <- reg.int.ci
  out$reg.pval <- reg.pval
  out$reg.pval.resam <- reg.pval.resam
  out$skewness <- skewness
  out$skewness.ci <- skewness.ci
  out$skewness.pval <- skewness.pval
  out$skewness.pval.resam <- skewness.pval.resam
  out$combined.pval <- combined.pval
  out$combined.pval.resam <- combined.pval.resam
  return(out)
}