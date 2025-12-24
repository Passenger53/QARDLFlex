#' Generate Summary Statistics for Data Frame Variables
#'
#' This function calculates comprehensive descriptive statistics for numeric variables
#' in a data frame, including measures of central tendency, dispersion, shape, and
#' normality tests. Results can be exported to a file.
#'
#' @param data A data frame containing the variables to analyze. If the first column
#'        is named "date", it will be excluded from statistical calculations.
#' @param export Logical. If TRUE, results are exported to a file. Default is TRUE.
#' @param file Character. File path for exporting results. Required if export = TRUE.
#'
#' @return A tibble with one row per variable and columns for each statistic:
#'         variable name, mean, median, maximum, minimum, standard deviation,
#'         skewness, kurtosis, and Jarque-Bera test statistic with significance stars.
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
#' # Calculate summary statistics
#' stats <- summary_statistics(sample_data, export = FALSE)
#' print(stats)
summary_statistics <- function(data, export = TRUE, file = NULL){

  # Remove date column if present (first column named "date")
  if (colnames(data)[1] == "date"){
    dd <- dplyr::select(data, -date)
  } else {
    dd <- data
  }

  # Calculate all descriptive statistics for each variable
  # Using explicit package references with ::
  tibble::tibble(
    Mean = purrr::map_dbl(dd, base::mean, na.rm = TRUE),
    Median = purrr::map_dbl(dd, stats::median, na.rm = TRUE),
    Maximum = purrr::map_dbl(dd, base::max, na.rm = TRUE),
    Minimum = purrr::map_dbl(dd, base::min, na.rm = TRUE),
    Std.Dev. = purrr::map_dbl(dd, stats::sd, na.rm = TRUE),
    Skewness = purrr::map_dbl(dd, e1071::skewness, na.rm = TRUE),
    Kurtosis = purrr::map_dbl(dd, e1071::kurtosis, na.rm = TRUE),
    `Jarque-Bera` = purrr::map_dbl(dd, ~ tseries::jarque.bera.test(.)$statistic),
    pvalue = purrr::map_dbl(dd, ~ tseries::jarque.bera.test(.)$p.value)
  ) -> xx

  # Add variable names as first column and rearrange
  xx %>%
    dplyr::mutate(variable = colnames(dd)) %>%
    dplyr::relocate(variable, .before = Mean) -> xx

  # Add significance stars to Jarque-Bera statistic based on p-value
  # *** p < 0.01, ** p < 0.05, * p < 0.1
  xx %>%
    dplyr::mutate(`Jarque-Bera` = dplyr::case_when(
      pvalue < 0.01 ~ stringr::str_c(sprintf("%.3f", `Jarque-Bera`), "***"),
      pvalue < 0.05 ~ stringr::str_c(sprintf("%.3f", `Jarque-Bera`), "**"),
      pvalue < 0.1 ~ stringr::str_c(sprintf("%.3f", `Jarque-Bera`), "*"),
      TRUE ~ sprintf("%.3f", `Jarque-Bera`)
    )) %>%
    dplyr::select(-pvalue) -> xx

  # Format all numeric columns to 3 decimal places for consistent presentation
  xx %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ sprintf("%.3f", .x))) -> xx

  # Export results if requested
  if (export){
    if (is.null(file)) {
      stop("File path must be specified when export = TRUE")
    }
    xx %>% export(file)
  }

  return(xx)
}
