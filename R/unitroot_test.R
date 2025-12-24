#' Comprehensive Unit Root Testing for Multiple Variables
#'
#' This function performs multiple unit root tests (ADF, PP, KPSS, ERS, ZA) on all
#' numeric variables in a dataset at level, first difference, and second difference.
#' It provides a comprehensive stationarity analysis for time series variables.
#'
#' @param data A data frame containing the time series variables to test.
#'        If the first column is named "date", it will be excluded from testing.
#' @param export Logical. If TRUE, results are exported to a file. Default is TRUE.
#' @param file Character. File path for exporting results. Required if export = TRUE.
#'
#' @return A tibble with unit root test results for all variables at different
#'         differencing levels, including test statistics with significance stars.
#'         Results are organized by variable, differencing level, and test type.
#'
#' @export
#'
#' @examples
#' # Generate sample data using sim_data_generate()
#' sample_data <- sim_data_generate(
#'   b = c(0.5, 0.7, 0.3, 0.2),
#'   order = c(1, 1),
#'   size = 100,
#'   reps = 1
#' )[[1]]
#'
#' # Perform unit root tests
#' unitroot_results <- unitroot_test(sample_data, export = FALSE)
#' print(unitroot_results)
unitroot_test <- function(data, export = TRUE, file = NULL){
  # Remove date column if present (first column named "date")
  if (colnames(data)[1] == "date"){
    dd <- dplyr::select(data, -date)
  } else {
    dd <- data
  }

  # Perform unit root tests on level data for all variables
  purrr::map(dd, unitroot) -> xx
  # Add variable names to results
  for (i in 1:length(xx)){
    xx[[i]] %>% dplyr::mutate(variable = colnames(dd[i])) -> xx[[i]]
  }
  xx <- dplyr::bind_rows(xx)
  xx %>% dplyr::relocate(variable, .before = test) -> xx

  # Perform unit root tests on first differences
  dd %>%
    purrr::map_dfc(base::diff) %>%
    purrr::map(unitroot) -> xx1
  for (i in 1:length(xx1)){
    xx1[[i]] %>% dplyr::mutate(variable = colnames(dd[i])) -> xx1[[i]]
  }
  xx1 <- dplyr::bind_rows(xx1)
  xx1 %>% dplyr::relocate(variable, .before = test) -> xx1

  # Perform unit root tests on second differences
  dd %>%
    purrr::map_dfc(base::diff) %>%
    purrr::map_dfc(base::diff) %>%
    purrr::map(unitroot) -> xx2
  for (i in 1:length(xx2)){
    xx2[[i]] %>% dplyr::mutate(variable = colnames(dd[i])) -> xx2[[i]]
  }
  xx2 <- dplyr::bind_rows(xx2)
  xx2 %>% dplyr::relocate(variable, .before = test) -> xx2

  # Combine results and add differencing level information
  xx %>% dplyr::mutate(difference = "level") %>%
    dplyr::relocate(difference, .after = variable) -> xx
  xx1 %>% dplyr::mutate(difference = "1st difference") %>%
    dplyr::relocate(difference, .after = variable) -> xx1
  xx2 %>% dplyr::mutate(difference = "2nd difference") %>%
    dplyr::relocate(difference, .after = variable) -> xx2
  xx <- dplyr::bind_rows(xx, xx1, xx2)

  # Add significance stars to test statistics
  # For KPSS test: higher values indicate non-stationarity (reject stationarity)
  # For other tests: lower values indicate stationarity (reject unit root)
  xx %>%
    dplyr::mutate(statistic = base::ifelse(
      test == "KPSS",
      dplyr::case_when(
        statistic > `1pct` ~ stringr::str_c(
          base::round(statistic, 3), "***"
        ),
        statistic > `5pct` ~ stringr::str_c(
          base::round(statistic, 3), "**"
        ),
        statistic > `10pct` ~ stringr::str_c(
          base::round(statistic, 3), "*"
        ),
        TRUE ~ base::as.character(
          base::round(statistic, 3)
        )
      ),
      dplyr::case_when(
        statistic < `1pct` ~ stringr::str_c(
          base::round(statistic, 3), "***"
        ),
        statistic < `5pct` ~ stringr::str_c(
          base::round(statistic, 3), "**"
        ),
        statistic < `10pct` ~ stringr::str_c(
          base::round(statistic, 3), "*"
        ),
        TRUE ~ base::as.character(
          base::round(statistic, 3)
        )
      )
    )) -> xx

  # Export results if requested
  if (export){
    if (base::is.null(file)) {
      stop("File path must be specified when export = TRUE")
    }
    xx %>% rio::export(file)
  }

  return(xx)
}
