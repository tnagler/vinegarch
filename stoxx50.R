library(VineCopula)
library(foreach)
library(doParallel)
load("sp500_lret.RData")

stoxx <- lret[, 1:50]

margin_models <- lapply(stoxx, select_garch)

# transform to pseudo observaions of the copula
u <- sapply(margin_models, get_u_from_garch)

system.time(fit1 <- RVineSparseSelect(u, 1:10, NA, NA, method = "itau"))
system.time(fit2 <- RVineSparseSelect(u, 1:10, 0, 50, method = "itau"))


thres <- c(NA, NA, 0, rep(seq(0, 0.15, l = 20), 2), rep(0, 49))
trnuc <- c(NA, Inf, NA,  rep(Inf, 20), rep(NA, 20), 1:49)


foreach(
    i = seq_along(thresholds)
) %dopar% {
    fit <- RVineSparseSelect(u, 1:6, thres[i], trunc[i], method = "itau")
    data.frame(
        mBIC = RVineMBIC(u, fit),
        BIC = fit$BIC
    )
}