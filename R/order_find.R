#' Determine Optimal Lag Order for ARDL Model
#'
#' This function automatically selects or manually specifies the optimal lag order
#' for ARDL model variables. It supports both automatic lag selection based on
#' information criteria and manual specification of lag orders.
#'
#' @param data A data frame containing the time series data. The first column
#'        should be named "date" and will be excluded from analysis.
#' @param selection Character. Information criterion for automatic lag selection.
#'        Currently supports "AIC" (Akaike Information Criterion, tends to select
#'        models with better prediction accuracy) and "BIC" (Bayesian Information
#' @param auto_lag Logical. If TRUE (default), automatically selects optimal lags
#'        using specified information criterion. If FALSE, uses manually provided lags.
#' @param lags Integer vector. Manual lag orders for each variable. Required when
#'        auto_lag = FALSE. Length should match number of variables (excluding date).
#' @param export Logical. If TRUE (default), exports results to a file.
#' @param file Character. File path for exporting results. Required if export = TRUE.
#'
#' @return An integer vector containing the optimal lag orders for each variable,
#'         with variable names as the vector names.
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
#' # Automatic lag selection using BIC criterion
#' optimal_lags <- order_find(sample_data, selection = "BIC", export = FALSE)
#' print(optimal_lags)
#'
#' # Manual lag specification
#' manual_lags <- order_find(sample_data, auto_lag = FALSE, lags = c(2, 1), export = FALSE)
#' print(manual_lags)
order_find <- function(data, selection, auto_lag = TRUE, lags = NULL, export = TRUE, file = NULL){

  # Automatic lag selection using information criteria
  if (auto_lag){
    # Generate formula from data structure
    formula <- function_generate(data)
    # Remove date column for analysis
    dd <- dplyr::select(data, -date)

    # Automatically select optimal ARDL lags using specified criterion
    # max_order sets maximum lags to consider for dependent and independent variables
    xx <- ARDL::auto_ardl(formula,
                    dd,
                    max_order = c(5, rep(4, base::ncol(dd) - 1)),
                    selection = selection)
    order <- xx$best_order

    # Display results
    cat("Optimal lag orders based on", selection, "criterion:\n")
    base::print(order)

    # Export results if requested
    if (export){
      rio::export(tibble::tibble(variable = base::names(order),
                                 order = order,
                                 criterion = selection),
                  file)
    }
  } else {
    # Manual lag specification
    # Use provided lag orders and assign variable names
    order <- lags
    base::names(order) <- base::names(data)[-1]

    # Display manually specified lags
    cat("Manually specified lag orders:\n")
    base::print(order)

    # Export results if requested
    if (export){
      rio::export(tibble::tibble(variable = base::names(order),
                                 order = order,
                                 criterion = "fixed"),
                  file)
    }
  }

  # Return lag orders as named vector
  return(order)
}
