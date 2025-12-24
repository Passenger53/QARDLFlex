#' Unit Root Test Function (Internal)
#'
#' This internal function performs multiple unit root tests on a time series vector,
#' including ADF, PP, KPSS, ERS, and Zivot-Andrews tests. It returns a tibble with
#' test statistics and critical values for each test specification.
#'
#' @param x A numeric vector representing the time series to test for unit roots.
#'
#' @return A tibble containing test results with columns: test, selectlags, type,
#'         statistic, lag, critical values (1%, 5%, 10%), and breakpoint (for ZA test).
#'
#' @keywords internal
unitroot <- function(x){
  # Initialize results container
  res <- NULL

  # ADF Test: Test different specifications (none, drift, trend) and lag selection criteria
  types <- c("none", "drift", "trend")
  criterion <- c("BIC", "AIC")
  for (i in 1:length(types)){
    for (j in 1:length(criterion)){
      xx <- urca::ur.df(x, selectlags = criterion[j], type = types[i])
      res1 <- tibble::tibble(test = "ADF",
                             selectlags = criterion[j],
                             type = types[i],
                             statistic = xx@teststat[1],
                             lag = xx@lags,
                             `1pct` = xx@cval[1,1],
                             `5pct` = xx@cval[1,2],
                             `10pct` = xx@cval[1,3])
      res <- dplyr::bind_rows(res, res1)
    }
  }

  # PP Test: Phillips-Perron test with constant and trend models
  types <- c("constant", "trend")
  for (i in 1:length(types)){
    xx <- urca::ur.pp(x, model = types[i], type = "Z-tau")
    res1 <- tibble::tibble(test = "PP",
                           selectlags = NA,
                           type = types[i],
                           statistic = xx@teststat[1],
                           lag = xx@lag,
                           `1pct` = xx@cval[1,1],
                           `5pct` = xx@cval[1,2],
                           `10pct` = xx@cval[1,3])
    res <- dplyr::bind_rows(res, res1)
  }

  # KPSS Test: Kwiatkowski-Phillips-Schmidt-Shin test for stationarity
  xx <- urca::ur.kpss(x)
  res1 <- tibble::tibble(test = "KPSS",
                         selectlags = NA,
                         type = NA,
                         statistic = xx@teststat[1],
                         lag = xx@lag,
                         `1pct` = xx@cval[1,4],
                         `5pct` = xx@cval[1,2],
                         `10pct` = xx@cval[1,1])
  res <- dplyr::bind_rows(res, res1)

  # ERS Test: Elliott-Rothenberg-Stock test with constant and trend models
  types <- c("constant", "trend")
  for (i in 1:length(types)){
    xx <- urca::ur.ers(x, model = types[i])
    res1 <- tibble::tibble(test = "ERS",
                           selectlags = NA,
                           type = types[i],
                           statistic = xx@teststat[1],
                           lag = xx@lag,
                           `1pct` = xx@cval[1,1],
                           `5pct` = xx@cval[1,2],
                           `10pct` = xx@cval[1,3])
    res <- dplyr::bind_rows(res, res1)
  }

  # Prepare for ZA test by adding breakpoint column
  res %>% dplyr::mutate(bpoint = NA) -> res

  # ZA Test: Zivot-Andrews test with structural breaks (intercept, trend, both)
  types <- c("intercept", "trend", "both")
  for (i in 1:length(types)){
    xx <- urca::ur.za(x, model = types[i])
    res1 <- tibble::tibble(test = "ZA",
                           selectlags = NA,
                           type = types[i],
                           statistic = xx@teststat[1],
                           lag = xx@lag,
                           `1pct` = xx@cval[1],
                           `5pct` = xx@cval[2],
                           `10pct` = xx@cval[3],
                           bpoint = xx@bpoint)
    res <- dplyr::bind_rows(res, res1)
  }
  return(res)
}
