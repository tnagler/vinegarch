select_aparch <- function(x) {
    # model formulations to select from
    formulas <- list(
        ~ arma(1, 1) + aparch(1, 1),
        ~ arma(0, 1) + aparch(1, 1),
        ~ arma(1, 0) + aparch(1, 1),
        ~ arma(0, 0) + aparch(1, 1)
    )
    # fit all models
    models <- suppressWarnings(
        lapply(formulas, function(f) garchFit(formula = f,
                                              data = x,
                                              cond.dist = 'std',
                                              trace = FALSE))
    )
    # choose the model with smallest BIC (second entry of ics)
    models[[which.min(sapply(models, function(m) m@fit$ics[2]))]]
}

# apply the conditional marginal distribution (from aparch) to a time series
get_u_from_aparch <- function(models) {
    sapply(models, function(m)
        pstd(residuals(m, standardize = TRUE), nu = m@fit$par['shape'])
    )
}


aparch_coefs <- function(m) {
    # extract coefficients from fitte model
    out <- c(
        mu     = as.numeric(m@fit$par['mu']),
        ar1    = as.numeric(m@fit$par['ar1']),
        ma1    = as.numeric(m@fit$par['ma1']),
        omega  = as.numeric(m@fit$par['omega']),
        alpha1 = as.numeric(m@fit$par['alpha1']),
        gamma1 = as.numeric(m@fit$par['gamma1']),
        beta1  = as.numeric(m@fit$par['beta1']),
        delta  = as.numeric(m@fit$par['delta']),
        nu     = as.numeric(m@fit$par["shape"])
    )
    # set missing entries to zero
    out[is.na(out)] <- 0

    out
}

csim_aparch <- function(model, x_tm, eps_tm, sigma_tm, n_ahead = 1, u = runif(10)) {
    # extract model coefficients
    cf <- aparch_coefs(model)
    stopifnot(NCOL(u) == n_ahead)

    # initialize new series
    n_MC <- NROW(u)
    x <- mu <- sigma <- z <- eps <- matrix(NA, n_MC, n_ahead + 1)

    # start with last time point in training time series (t - 1)
    x[, 1] <- x_tm
    sigma[, 1] <- sigma_tm
    eps[, 1] <- eps_tm

    # simulate standardized residuals for full horizon
    z[, -1] <- qstd(u, nu = cf['nu'])

    ## simulate aparch series
    for (t in seq.int(n_ahead)) {
        mu[, t + 1] <- cmean_aparch(model, x[, t], eps[, t])
        sigma[, t + 1] <- csig_aparch(model, eps[, t], sigma[, t])
        eps[, t + 1] <- sigma[, t + 1] * z[, t + 1]
        x[, t + 1] <- mu[, t + 1] + eps[, t + 1]
    }

    x <- x[, -1, drop = "FALSE"]
    colnames(x) <- paste0("t_", seq.int(n_ahead))
    x
}

cmean_aparch <- function(model, x_tm, eps_tm) {
    cf <- aparch_coefs(model)
    cf['mu'] + cf['ma1'] * eps_tm + cf['ar1'] * x_tm
}

csig_aparch <- function(model, eps_tm, sigma_tm) {
    cf <- aparch_coefs(model)
    h <- cf['omega'] +
        cf['alpha1'] * (abs(eps_tm) - cf['gamma1'] * eps_tm)^cf['delta'] +
        cf['beta1'] * sigma_tm^cf['delta']
    h^(1 / cf['delta'])
}
