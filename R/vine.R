select_garchvine <- function(data, thresh, trunc, meth = "itau") {
    if (!inherits(data, "data.frame"))
        data <- as.data.frame(data)
    # select GARCH model for each margin
    margin_models <- lapply(data, select_garch)

    # transform to pseudo observaions of the copula
    u <- sapply(margin_models, get_u_from_garch)

    # fit vine copula model
    vine_model <- suppressWarnings(
        RVineSparseSelect(u,
                          familyset = 0:10,
                          threshold = thresh,
                          trunclevel = trunc,
                          method = meth)
    )

    # return margins + vine copula
    list(margin_models = margin_models, vine_model = vine_model)
}

continue_garchvine <- function(model, n_ahead = 1, times) {
    # simulate from vine copula
    u <- RVineSim(m, model$vine_model)

    # use copula data to generate innovations for GARCH
    d <- ncol(u)
    out <- sapply(seq.int(d), function(i)
        continue_garch(u[, i], model$margin_models[[i]])
    )

    # return as matrix
    if (!is.matrix(out))
        out <- t(as.matrix(out))
    out
}