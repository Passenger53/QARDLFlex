#' Single Variable Cross-Quantile Wald Test
#'
#' This function performs Wald tests for equality of coefficients across specified
#' quantiles for each variable in the ECM model. It tests both long-run and
#' short-run coefficients separately for each variable, with special handling
#' for the dependent variable (ECM term) and variables with zero ARDL order.
#' For short-run coefficients with multiple lags, the sum of coefficients across
#' lags is tested.
#'
#' @param estimation List. Estimation results obtained from `coef_estimate()`.
#' @param order Integer vector. ARDL orders for each variable, typically obtained
#'   from `order_find()`. The length must match the number of variables in the
#'   data. The first element should be the dependent variable.
#' @param boots Integer. Number of bootstrap replications used in the estimation.
#'   This should match the value used in `coef_estimate()`.
#' @param taus Numeric vector. All quantile levels estimated, same as the
#'   `taus` parameter in `coef_estimate()`.
#' @param test_taus Numeric vector. Specific quantile levels to test for equality.
#'   Must be a subset of `taus` and contain at least two distinct quantiles.
#'   All values must be between 0 and 1.
#' @param export Logical. Whether to export the results to a file. Default is TRUE.
#' @param file Character. File path for exporting results. Required if export = TRUE.
#'
#' @return A tibble containing Wald test results for each variable, including:
#'   - Variable: Name of the variable tested
#'   - term: Type of coefficient tested ("Long", "Short", or NA for special cases)
#'   - Wald: Wald test statistic with significance stars (***p<0.01, **p<0.05, *p<0.1)
#'   - pvalue: p-value of the test
#'   - Tau: The quantile levels tested, concatenated as a string
#'
#' @export
#'
#' @examples
#' \dontrun{
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
#' # Estimate models at multiple quantiles
#' est_result <- coef_estimate(
#'   data = sample_data,
#'   order = optimal_order,
#'   taus = c(0.25, 0.5, 0.75),
#'   boots = 50,
#'   export = FALSE
#' )
#'
#' # Perform cross-quantile Wald tests for all variables
#' wald_results <- single_variable_cross_wald(
#'   estimation = est_result,
#'   order = optimal_order,
#'   boots = 50,
#'   taus = c(0.25, 0.5, 0.75),
#'   test_taus = c(0.25, 0.5, 0.75),
#'   export = FALSE
#' )
#'
#' # View results
#' print(wald_results)
#' }
single_variable_cross_wald <- function(estimation, order, boots, taus, test_taus, export = TRUE, file) {
  # Validate that test_taus has at least two distinct quantiles
  if (base::length(base::unique(test_taus)) < 2) {
    stop("The test_taus parameter must contain at least two distinct quantile levels.")
  }

  # Validate that all test_taus are in the original taus
  if (!base::all(test_taus %in% taus)) {
    stop("All quantile levels in test_taus must be present in the original taus parameter.")
  }

  # Get variable names from the order vector
  names <- base::names(order)

  # Initialize container for Wald test results
  wald_results <- NULL

  # Loop through each variable to perform cross-quantile tests
  for (i in 1:base::length(names)) {

    # Case 1: Dependent variable (first variable in the order vector)
    if (i == 1) {
      # Test ECM coefficient for the dependent variable
      variable <- "ECM"
      xx <- cross_quantile_test(
        estimation = estimation,
        order = order,
        boots = boots,
        taus = taus,
        test_taus = test_taus,
        variable = variable,
        term = NULL  # ECM term doesn't have long/short distinction
      )

      # Format results for ECM term
      xx <- tibble::tibble(
        Variable = "ECM",
        term = base::as.character(NA),  # NA for ECM term
        Wald = xx$statistic,
        pvalue = xx$p.value,
        Tau = stringr::str_c(
          base::sprintf("%.2f", test_taus),
          collapse = "="
        )
      )
      wald_results <- dplyr::bind_rows(wald_results, xx)

      # If dependent variable has ARDL order > 1, test its short-run coefficients
      if (order[i] > 1) {
        variable <- names[i]
        xx <- cross_quantile_test(
          estimation = estimation,
          order = order,
          boots = boots,
          taus = taus,
          test_taus = test_taus,
          variable = variable,
          term = "short"
        )

        xx <- tibble::tibble(
          Variable = variable,
          term = "Short",
          Wald = xx$statistic,
          pvalue = xx$p.value,
          Tau = stringr::str_c(
            base::sprintf("%.2f", test_taus),
            collapse = "="
          )
        )
        wald_results <- dplyr::bind_rows(wald_results, xx)
      }
    }
    # Case 2: Independent variables with positive ARDL order
    else if (order[i] > 0) {
      variable <- names[i]

      # Test long-run coefficients
      xx <- cross_quantile_test(
        estimation = estimation,
        order = order,
        boots = boots,
        taus = taus,
        test_taus = test_taus,
        variable = variable,
        term = "long"
      )

      xx <- tibble::tibble(
        Variable = variable,
        term = "Long",
        Wald = xx$statistic,
        pvalue = xx$p.value,
        Tau = stringr::str_c(
          base::sprintf("%.2f", test_taus),
          collapse = "="
        )
      )
      wald_results <- dplyr::bind_rows(wald_results, xx)

      # Test short-run coefficients
      xx <- cross_quantile_test(
        estimation = estimation,
        order = order,
        boots = boots,
        taus = taus,
        test_taus = test_taus,
        variable = variable,
        term = "short"
      )

      xx <- tibble::tibble(
        Variable = variable,
        term = "Short",
        Wald = xx$statistic,
        pvalue = xx$p.value,
        Tau = stringr::str_c(
          base::sprintf("%.2f", test_taus),
          collapse = "="
        )
      )
      wald_results <- dplyr::bind_rows(wald_results, xx)
    }
    # Case 3: Variables with zero ARDL order (contemporaneous effect only)
    else {
      variable <- names[i]
      xx <- cross_quantile_test(
        estimation = estimation,
        order = order,
        boots = boots,
        taus = taus,
        test_taus = test_taus,
        variable = variable,
        term = "short"  # For zero-order variables, only contemporaneous effect exists
      )

      xx <- tibble::tibble(
        Variable = variable,
        term = base::as.character(NA),  # NA for zero-order variables
        Wald = xx$statistic,
        pvalue = xx$p.value,
        Tau = stringr::str_c(
          base::sprintf("%.2f", test_taus),
          collapse = "="
        )
      )
      wald_results <- dplyr::bind_rows(wald_results, xx)
    }
  }

  # Format Wald statistics with significance stars and p-values
  wald_results %>%
    dplyr::mutate(
      Wald = base::sprintf("%.4f", Wald),
      Wald = dplyr::case_when(
        pvalue < 0.01 ~ stringr::str_c(Wald, "***"),
        pvalue < 0.05 ~ stringr::str_c(Wald, "**"),
        pvalue < 0.1 ~ stringr::str_c(Wald, "*"),
        TRUE ~ Wald
      )
    ) %>%
    dplyr::mutate(pvalue = base::sprintf("%.4f", pvalue)) -> wald_results

  # Export results if requested
  if (export) {
    if (base::is.null(file)) {
      stop("File path must be specified when export = TRUE")
    }
    rio::export(wald_results, file)
  }

  return(wald_results)
}
