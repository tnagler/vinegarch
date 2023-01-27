library(vinegarch)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

do_backtest <- function(data, thresh, trunc, n_train = 750, n_test = 250,
                        n_MC = 10^5, cores = 1) {
    data <- as.data.frame(data)
    n <- nrow(data)
    d <- ncol(data)
    VaR <- matrix(NA, nrow(data), 3)
    cmean <- eps <- csig <- matrix(NA, nrow(data), d)

    for (t in seq.int(n - 1)) {
        print(t)
        # dont start before n_train time points have passed
        if (t < n_train)
            next

        ## always refit the model after n_test time points
        if ((t - n_train) %% 125 == 0) {
            message("fit")
            # extract new training data
            train_data <- data[(t - n_train + 1):t, ]

            # select GARCH model for each margin
            garchs <- lapply(train_data, select_garch)

            # transform to pseudo observaions of the copula
            u <- sapply(garchs, get_u_from_garch)

            # fit and select vine copula model
            vine <- RVineSparseSelect(u,
                                      0:6,
                                      thresh,
                                      trunc,
                                      method = "itau",
                                      cores = cores)

            # Monte Carlo samples for the vine model
            u_MC <- RVineSim(n_MC, vine)

            # initialize with last day of training period
            eps[t, ]  <- sapply(garchs, function(m) m@residuals[n_train])
            csig[t, ] <- sapply(garchs, function(m) m@sigma.t[n_train])
        }

        # MC simulation for 1-day ahead returns
        x_ahead <- sapply(seq.int(d), function(i)
            csim_garch(
                model    = garchs[[i]],
                x_tm     = data[t, i],
                eps_tm   = eps[t, i],
                sigma_tm = csig[t, i],
                n_ahead  = 1,
                u        = u_MC[, i]
            )
        )

        # build portfolio (equally weighted)
        x_port <- rowMeans(x_ahead)

        # calcualte VaR forecast
        VaR[t + 1, ] <- quantile(x_port, 1 - c(0.9, 0.95, 0.99))

        # calculate 1-day ahead conditional mean and sd
        cmean[t + 1, ] <- sapply(seq.int(d), function(i)
            cmean_garch(garchs[[i]], data[t, i], eps[t, i])
        )
        csig[t + 1, ] <- sapply(seq.int(d), function(i)
            csig_garch(garchs[[i]], eps[t, i], csig[t, i])
        )

        # calculate residuals for time t + 1
        eps[t + 1, ] <- as.numeric(data[t + 1, ] - cmean[t + 1, ])
    }

    VaR <- as.data.frame(VaR)
    colnames(VaR) <- c("VaR90", "VaR95", "VaR99")
    cbind(index = rowMeans(data), VaR)
}

# load data
load("sp500_lret.RData")

res <- do_backtest(
    data = lret[2000:nrow(lret), 1:20],
    thresh = NA,
    trunc = NA,
    cores = 1,
    n_MC = 10^4
)

save(res, file = "back_test.RData")

load("back_test.RData")

library(tidyverse)

res %>%
    gather(key, value, -time) %>%
    ggplot(aes(x = time, y = value, color = key)) +
    geom_line()

res[complete.cases(res), ] %>%
    summarize(
        hits90 = sum(index < VaR90),
        hits95 = sum(index < VaR95),
        hits99 = sum(index < VaR99),
        n = length(VaR90),
        lr90 = binomial_lr(hits90, n, 0.10),
        lr95 = binomial_lr(hits95, n, 0.05),
        lr99 = binomial_lr(hits99, n, 0.01),
        p90 = pchisq(lr90, 1, lower.tail = FALSE),
        p95 = pchisq(lr95, 1, lower.tail = FALSE),
        p99 = pchisq(lr99, 1, lower.tail = FALSE)
    )

binomial_lr <- function(hits, n, alpha) {
    p_hat <- hits / n
    ll_alpha <- log(1 - alpha) * (n - hits) + log(alpha) * hits
    ll_hat   <- log(1 - p_hat) * (n - hits) + log(p_hat) * hits
    - 2  * (ll_alpha - ll_hat)
}

