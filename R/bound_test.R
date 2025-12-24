#' Perform Bounds F Test for Cointegration
#'
#' This function performs the bounds F test for cointegration in ARDL models,
#' testing the null hypothesis of no long-run relationship between variables.
#'
#' @param data A data frame containing the time series data. If the first column
#'        is named "date", it will be automatically excluded from the analysis.
#' @param order Integer vector. ARDL lag orders for each variable, typically obtained
#'        from `order_find()`. The length must match the number of variables in the
#'        data frame (excluding date).
#' @param export Logical. If TRUE (default), exports results to a file.
#' @param file Character. File path for exporting results. Required if export = TRUE.
#'
#' @return The function displays test results in the console. If export = TRUE,
#'         results are saved to the specified file. The test evaluates whether
#'         there is a long-run cointegrating relationship among the variables.
#'
#' @export
#'
#' @examples
#' # Generate sample data
#' sample_data <- sim_data_generate(
#'   b = c(0.5, 0.7, 0.3, 0.2),
#'   order = c(1, 1),
#'   size = 100,
#'   reps = 1
#' )[[1]]
#'
#' # Find optimal lag order
#' optimal_order <- order_find(sample_data, selection = "BIC", export = FALSE)
#'
#' # Perform bounds F test
#' bound_test(sample_data, order = optimal_order, export = FALSE)
bound_test <- function(data, order, export = TRUE, file = NULL){
  # Remove date column if first column is named "date"
  if (base::colnames(data)[1] == "date"){
    dd <- dplyr::select(data, -date)
  } else {
    dd <- data
  }

  # Create formula for ARDL model: first variable ~ all other variables
  fomu <- base::colnames(dd)
  fomu <- stringr::str_c(fomu[1], "~", stringr::str_c(fomu[-1], collapse = "+"))
  fomu <- stats::as.formula(fomu)

  # Estimate ARDL model using ARDL package
  ardl <- ARDL::ardl(fomu, data = dd, order = order)

  # Perform bounds F test for cointegration (case 2: restricted intercept, no trend)
  bound <- ARDL::bounds_f_test(ardl, case = 2)

  # Export results if requested
  if (export){
    if (base::is.null(file)) {
      stop("File path must be specified when export = TRUE")
    }
    tibble::tibble(statistic = bound$statistic,
                   pvalue = bound$p.value) %>%
      rio::export(file)
  }

  # Display test results in console
  cat("\n Bounds F test with statistic",
      base::sprintf("%.3f", bound$statistic),
      "and p-value",
      base::sprintf("%.3f", bound$p.value))
}
