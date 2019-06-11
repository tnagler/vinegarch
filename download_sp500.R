library(tidyquant)

# get prices for all sp500 constituents
sp500_stocks <- tq_index("SP500")
sp500_prices <- sp500_stocks %>%
    tq_get(get = "stock.prices", from = "2011-12-31", to = "2017-01-01")

# extract log returns
sp500_returns <- sp500_prices %>%
    group_by(symbol) %>%
    tq_transmute(ohlc_fun   = Ad,
                 mutate_fun = periodReturn,
                 period     = "daily",
                 type       = "log",
                 col_rename = "return") %>%
    ungroup() %>%
    spread(symbol, return)

# remove first day (return is zero)
sp500_returns <- sp500_returns[-1, ]

# remove stocks that don't have full history
sp500_returns <- sp500_returns[, colSums(is.na(sp500_returns)) < 1]
ncol(sp500_returns)
save(sp500_returns, file = "sp500_returns.RData")
