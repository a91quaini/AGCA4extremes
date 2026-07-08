#' Fama-French daily portfolio losses
#'
#' A transformed daily loss panel from the Kenneth R. French Data Library. The
#' data contain 24 daily portfolio loss series formed from four 2-by-3
#' portfolio-sort families: size and book-to-market, size and operating
#' profitability, size and investment, and size and prior return. Losses are
#' negative daily returns, measured in percent.
#'
#' @format A data frame with 15,833 rows and 25 columns. The first column is
#' the trading date. The remaining 24 columns are daily loss series for the
#' Fama-French portfolio sorts used as a package example.
#' @source Kenneth R. French Data Library,
#' \url{https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html}
"ff_portfolio_losses"
