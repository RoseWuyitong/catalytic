plot.cat_tau_glm <- function() {

}

plot.cat_freq_glm <- function(freq_fit) {
  plot(freq_fit$taus,
       freq_fit$risk_estimate_collect,
       type="l",
       col="black",
       lwd=2,
       xlab="taus",
       ylab=paste0("risk estimate (", freq_fit$risk_estimate_method,")"),
       main="Plot of Risk Estimate from Freqentist")

  text(freq_fit$optimal_tau,
       freq_fit$min_risk_estimate,
       labels = round(freq_fit$min_risk_estimate, 2),
       pos = 3,
       offset = 0.5)

  points(freq_fit$optimal_tau,
         freq_fit$min_risk_estimate,
         pch = 4,
         col = "red")

}

plot.cat_bayesian_glm <- function() {

}
