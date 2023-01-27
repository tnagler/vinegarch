#' Select ARMA-GARCH model (for marginal models)
#'
#' @param x time series.
#'
#' @return A model returned from [fGarch::garchFit()].
#' @export
#' @importFrom fGarch garchFit
#'
#' @examples
#' fit <- select_garch(rnorm(100))
#'
select_garch <- function(x) {
    ## model formulations to select from
    formulas <- list(
        ~ garch(1, 0),
        ~ arma(0, 1) + garch(1, 0),
        ~ arma(1, 0) + garch(1, 0),
        ~ arma(1, 1) + garch(1, 0),
        ~ garch(1, 1),
        ~ arma(0, 1) + garch(1, 1),
        ~ arma(1, 0) + garch(1, 1),
        ~ arma(1, 1) + garch(1, 1)
    )

    ## fit all models
    models <- suppressWarnings(
        lapply(formulas, fit_garch, data = x)
    )
    null_models <- which(sapply(models, is.null))
    if (length(null_models) > 0)
        models <- models[-null_models]
    if (length(models) == 0)
        stop("wasn't able to fit any model")

    ## choose the model with smallest BIC (second entry of ics)
    models[[which.min(sapply(models, function(m) m@fit$ics[2]))]]
}

fit_garch <- function(formula, data) {
    tryCatch(
        garchFit(formula = formula,
                 data = data,
                 cond.dist = 'std',
                 trace = FALSE),
        error = function(e) NULL
    )
}

#' PIT for GARCH models and its inverse.
#'
#' Applies the conditional marginal distribution (from GARCH) to the residuals.
#'
#' @param models a list of models returned from [fGarch::garchFit()].
#'
#' @export
#' @importFrom fGarch residuals pstd qstd
#'
#' @examples
#' fit <- select_garch(rnorm(100))
#' get_u_from_garch(fit)
#'
get_u_from_garch <- function(models) {
    if (inherits(models, "fGARCH"))
        models <- list(models)
    sapply(models, function(model)
        pstd(residuals(model, standardize = TRUE), nu = model@fit$par['shape'])
    )
}

#' @rdname get_u_from_garch
#' @param u a matrix of evaluation points.
#' @param time_index the index for the time point in the original data.
#' @return A matrix of conditional quantiles from the garch models (conditioned
#'     on `time == time_s[time_index]`).
#' @export
get_q_from_garch <- function(u, models, time_index) {
    stopifnot(ncol(u) == length(models))
    col_names <- colnames(u)
    q <- sapply(1:length(models), function(i) {
        x_t <- models[[i]]@fitted[time_index]
        sig_t <- models[[i]]@sigma.t[time_index]
        x_t + qstd(u[, i], nu = models[[i]]@fit$par['shape']) * sig_t
    })
    colnames(q) <- col_names
    q
}

garch_coefs <- function(m) {
    # extract coefficients from fitte model
    out <- c(
        mu     = as.numeric(m@fit$par['mu']),
        ar1    = as.numeric(m@fit$par['ar1']),
        ma1    = as.numeric(m@fit$par['ma1']),
        omega  = as.numeric(m@fit$par['omega']),
        alpha1 = as.numeric(m@fit$par['alpha1']),
        beta1  = as.numeric(m@fit$par['beta1']),
        sigma  = as.numeric(m@sigma.t[length(m@sigma.t)]),
        df     = as.numeric(m@fit$coef["shape"])
    )
    # set missing entries to zero
    out[is.na(out)] <- 0

    out
}

continue_garch <- function(u = runif(10), model) {
    # extract model coefficients
    cf <- garch_coefs(model)

    # initialize new series
    n <- length(u)
    x <- mu <- sigma <- z <- eps <- numeric(n + 1)

    # start with last time point in training time series (t - 1)
    last <- length(model@residuals)
    x[1] <- model@data[last]
    sigma[1] <- model@sigma.t[last]
    eps[1] <- model@residuals[last]

    # simulate standardized residuals for full horizon
    # (scale to unit variance)
    z[-1] <- qstd(u, nu = cf['df'])

    ## simulate GARCH series
    for (t in seq.int(n)) {
        mu[t + 1] <- cmean_garch(model, x[t], eps[t])
        sigma[t + 1] <- csig_garch(model, eps[t], sigma[t])
        # rescale with conditional variance
        eps[t + 1] <- sigma[t + 1] *  z[t + 1]
        x[t + 1] <- mu[t + 1] + eps[t + 1]
    }

    x[-1]
}

#' @export
#' @importFrom stats runif
csim_garch <- function(model, x_tm, eps_tm, sigma_tm, n_ahead = 1, u = runif(10)) {
    # extract model coefficients
    cf <- garch_coefs(model)
    stopifnot(NCOL(u) == n_ahead)

    # initialize new series
    n_MC <- NROW(u)
    x <- mu <- sigma <- z <- eps <- matrix(NA, n_MC, n_ahead + 1)

    # start with last time point in training time series (t - 1)
    x[, 1] <- x_tm
    sigma[, 1] <- sigma_tm
    eps[, 1] <- eps_tm

    # simulate standardized residuals for full horizon
    z[, -1] <- qstd(u, nu = cf['df'])

    ## simulate GARCH series
    for (t in seq.int(n_ahead)) {
        mu[, t + 1] <- cmean_garch(model, x[t], eps[t])
        sigma[, t + 1] <- csig_garch(model, eps[t], sigma[t])
        eps[, t + 1] <- sigma[, t + 1] * z[, t + 1]
        x[, t + 1] <- mu[, t + 1] + eps[, t + 1]
    }

    x <- x[, -1, drop = "FALSE"]
    colnames(x) <- paste0("t_", seq.int(n_ahead))
    x
}

#' @export
cmean_garch <- function(model, x_tm, eps_tm) {
    # extract model coefficients
    cf <- garch_coefs(model)
    cf['mu'] + cf['ma1'] * eps_tm + cf['ar1'] * x_tm
}

#' @export
csig_garch <- function(model, eps_tm, sigma_tm) {
    # extract model coefficients
    cf <- garch_coefs(model)
    h <- cf['omega'] + cf['beta1'] * sigma_tm^2  + cf['alpha1'] * eps_tm^2
    sqrt(h)
}
