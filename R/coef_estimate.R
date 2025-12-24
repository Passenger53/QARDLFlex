#' Estimate OLS and Quantile Regression ARDL Models
#'
#' This function exports estimation results for both OLS and quantile regression
#' ARDL models at specified quantile levels. It provides comprehensive results
#' including ARDL and ECM specifications, with the option to export formatted
#' coefficient tables.
#'
#' @param data A data frame containing the time series data. If the first column
#'   is named "date", it will be automatically excluded from the analysis.
#' @param order Integer vector. ARDL lag orders for each variable, typically obtained
#'   from `order_find()`. The length must match the number of variables in the
#'   data frame (excluding date).
#' @param taus Numeric vector. Quantile levels to estimate, all values must be
#'   between 0 and 1. For example, c(0.25, 0.5, 0.75) for 25th, 50th (median),
#'   and 75th percentiles.
#' @param boots Integer. The number of bootstrap replications for quantile
#'   regression standard error estimation. Higher values provide more accurate
#'   confidence intervals but increase computation time.
#' @param export Logical. If TRUE, exports formatted coefficient tables to files.
#'   Default is FALSE.
#' @param file_ardl Character. File path for exporting ARDL coefficient table.
#'   Required if export = TRUE.
#' @param file_ecm Character. File path for exporting ECM coefficient table.
#'   Required if export = TRUE.
#'
#' @return A list containing two elements:
#'         - OLS: Full OLS estimation results from `ardl_estimate()`
#'         - QR: List of quantile regression estimation results, one for each
#'           quantile level specified in `taus`, from `qardl_estimate()`
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
#' # Estimate OLS and multiple quantile regressions
#' results <- coef_estimate(
#'   data = sample_data,
#'   order = optimal_order,
#'   taus = c(0.25, 0.5, 0.75),
#'   boots = 50,
#'   export = FALSE
#' )
#'
#' # Access OLS results
#' print(results$OLS$ARDL_coefficient)
#'
#' # Access median regression results
#' print(results$QR[[2]]$QARDL_coefficient)  # Second element corresponds to tau = 0.5
coef_estimate <- function(data, order, taus, boots, export = FALSE, file_ardl = NULL, file_ecm = NULL){
  # Estimate OLS ARDL and ECM models
  ardl_result <- ardl_estimate(data = data, order = order)

  # Initialize list to store quantile regression results
  qardl_results <- base::vector("list", base::length(taus))

  # Estimate quantile regression models for each specified quantile level
  for (i in 1:base::length(taus)){
    tau <- taus[i]
    qardl_results[[i]] <- qardl_estimate(data = data,
                                         order = order,
                                         tau = tau,
                                         boots = boots)
    cat("\n Quantile regression completed: Tau = ", base::sprintf("%.2f", taus[i]))
  }

  # Extract and rename OLS coefficient tables
  coef_ardl <- ardl_result$ARDL_coefficient
  coef_ardl %>% dplyr::rename(OLS = Value) -> coef_ardl

  coef_ecm <- ardl_result$ECM_coefficient
  coef_ecm %>% dplyr::rename(OLS = Value) -> coef_ecm

  # Merge quantile regression results with OLS results
  for (i in 1:base::length(taus)){
    # Extract and rename quantile ARDL coefficients
    coef_ardl1 <- qardl_results[[i]]$QARDL_coefficient
    base::names(coef_ardl1)[3] <- stringr::str_c("Tau", base::sprintf("%.2f", taus[i]))

    # Left join to combine OLS and quantile regression results
    coef_ardl <- dplyr::left_join(coef_ardl, coef_ardl1,
                                  by = dplyr::join_by(Variable, `Coef (Std. Error)`))

    # Extract and rename quantile ECM coefficients
    coef_ecm1 <- qardl_results[[i]]$QECM_coefficient
    base::names(coef_ecm1)[3] <- stringr::str_c("Tau", base::sprintf("%.2f", taus[i]))

    # Left join to combine OLS and quantile regression results
    coef_ecm <- dplyr::left_join(coef_ecm, coef_ecm1,
                                 by = dplyr::join_by(Variable, `Coef (Std. Error)`))
  }

  # Export formatted coefficient tables if requested
  if (export){
    if (base::is.null(file_ardl) || base::is.null(file_ecm)) {
      stop("File paths for both ARDL and ECM must be specified when export = TRUE")
    }
    rio::export(coef_ardl, file_ardl)
    rio::export(coef_ecm, file_ecm)
  }

  # Create comprehensive results list
  estimates = base::list(OLS = ardl_result,
                         QR = qardl_results)

  return(estimates)
}
