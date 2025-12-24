#' Within-Quantile Test for Coefficient Equality
#'
#' This function performs Wald tests to check if coefficients of multiple
#' explanatory variables are equal within each quantile level. It tests both
#' long-run and short-run coefficients separately for each quantile. This is
#' particularly useful for NARDL-type models where variables are decomposed
#' into multiple components.
#'
#' @param estimation List. Estimation results obtained from `coef_estimate()`.
#' @param order Integer vector. ARDL orders for each variable, typically obtained
#'   from `order_find()`. The length must match the number of variables in the
#'   data. The first element should be the dependent variable.
#' @param boots Integer. Number of bootstrap replications used in the estimation.
#'   This should match the value used in `coef_estimate()`.
#' @param taus Numeric vector. All quantile levels estimated, same as the
#'   `taus` parameter in `coef_estimate()`.
#' @param variables Character vector. Names of explanatory variables to test
#'   for equality. Must contain at least two variable names and cannot include
#'   the dependent variable. All variables must have positive ARDL orders.
#' @param export Logical. Whether to export the results to a file. Default is TRUE.
#' @param file Character. File path for exporting results. Required if export = TRUE.
#'
#' @return A tibble containing Wald test results for each quantile level, including:
#'   - Variable: Concatenated variable names tested for equality
#'   - term: Type of coefficient tested ("Long" or "Short")
#'   - Wald: Wald test statistic with significance stars (***p<0.01, **p<0.05, *p<0.1)
#'   - pvalue: p-value of the test
#'   - Tau: The quantile level at which the test was performed
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' sample_data <- sim_data_generate(
#'   b = c(0.5, 0.7, 0.3, 0.4, 0.1, 0.2),
#'   order = c(1, 1, 1),
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
#' # Test if coefficients of two explanatory variables are equal at each quantile
#' test_results <- within_quantile_test(
#'   estimation = est_result,
#'   order = optimal_order,
#'   boots = 50,
#'   taus = c(0.25, 0.5, 0.75),
#'   variables = c("x1", "x2"),
#'   export = FALSE
#' )
#'
#' # View results
#' print(test_results)
#' }
within_quantile_test <- function(estimation, order, boots, taus, variables, export = TRUE, file = NULL){

  # Extract quantile regression models from estimation results
  models <- estimation$QR

  # Initialize list to store ECM models for each quantile
  ECM_models <- base::vector("list", base::length(taus))

  # Extract ECM coefficients and bootstrap samples for each quantile
  for (i in 1:base::length(taus)){
    ECM_models[[i]] <- base::list(
      boot = models[[i]]$QECM_boot,
      coef = models[[i]]$QECM_coef
    )
  }

  # Get coefficient names from the first quantile model
  coef_names <- base::names(ECM_models[[1]]$coef)

  # Validate input parameters
  if (base::length(variables) < 2){
    stop("The 'variables' parameter must contain at least two explanatory variables!")
  }

  # Ensure dependent variable is not included in test variables
  if (base::names(order)[1] %in% variables){
    stop("The 'variables' parameter can only contain explanatory variables, not the dependent variable!")
  }

  # Find coefficient indices for each variable
  find <- base::vector("list", base::length(variables))
  for (i in 1:base::length(variables)){
    find[[i]] <- base::which(stringr::str_detect(coef_names, variables[i]))
  }

  # Check that all variables have positive ARDL orders (non-zero coefficients)
  if (!base::all(purrr::map_dbl(find, base::length) > 1)){
    stop("Some variables in 'variables' have ARDL order zero and cannot be tested with Wald test!")
  }

  # Separate long-run and short-run coefficients
  find_long <- NULL
  find_short <- base::vector("list", base::length(find))

  for (i in 1:base::length(find)){
    find_long <- base::c(find_long, find[[i]][1])  # First coefficient is long-run
    find_short[[i]] <- find[[i]][-1]  # Remaining coefficients are short-run
  }

  # Initialize containers for coefficients and bootstrap samples
  coef_long <- NULL
  coef_short <- NULL
  boot_long <- base::vector("list", base::length(taus))
  boot_short <- base::vector("list", base::length(taus))

  # Extract coefficients and bootstrap samples for each quantile
  for (i in 1:base::length(taus)){
    # Long-run coefficients (first coefficient for each variable)
    coef_long <- base::rbind(coef_long, ECM_models[[i]]$coef[find_long])
    boot_long[[i]] <- ECM_models[[i]]$boot[, find_long]

    # Short-run coefficients (sum of remaining coefficients for each variable)
    ss <- NULL
    bb <- NULL
    for (j in 1:base::length(find_short)){
      ss <- base::c(ss, base::sum(ECM_models[[i]]$coef[find_short[[j]]]))
      bb <- base::cbind(bb, base::rowSums(ECM_models[[i]]$boot[, find_short[[j]]]))
    }
    coef_short <- base::rbind(coef_short, ss)
    boot_short[[i]] <- bb
  }

  # Set row and column names for better identification
  base::rownames(coef_long) <- stringr::str_c("Tau", base::sprintf("%.2f", taus))
  base::names(boot_long) <- stringr::str_c("Tau", base::sprintf("%.2f", taus))
  base::rownames(coef_short) <- stringr::str_c("Tau", base::sprintf("%.2f", taus))
  base::names(boot_short) <- stringr::str_c("Tau", base::sprintf("%.2f", taus))

  # Calculate covariance matrices for bootstrap samples
  covs_long <- base::vector("list", base::length(taus))
  covs_short <- base::vector("list", base::length(taus))

  for (i in 1:base::length(taus)){
    covs_long[[i]] <- stats::cov(boot_long[[i]])
    covs_short[[i]] <- stats::cov(boot_short[[i]])
  }

  base::names(covs_long) <- stringr::str_c("Tau", base::sprintf("%.2f", taus))
  base::names(covs_short) <- stringr::str_c("Tau", base::sprintf("%.2f", taus))

  # Create contrast matrix for Wald test (testing equality of coefficients)
  R <- base::matrix(0, nrow = base::length(variables) - 1, ncol = base::length(variables))
  R[, 1] <- 1
  for (i in 1:base::nrow(R)){
    R[i, i + 1] <- -1
  }

  # Initialize results container
  res <- NULL

  # Perform Wald tests for each quantile level
  for (i in 1:base::length(taus)){
    # Test long-run coefficients
    xx <- wald_test(b = coef_long[i, ], cov = covs_long[[i]], R = R)

    xx <- tibble::tibble(
      Variable = stringr::str_c(variables, collapse = "="),
      term = "Long",
      Wald = base::sprintf("%.4f", xx$statistic),
      pvalue = xx$p.value,
      Tau = base::sprintf("%.2f", taus[i])
    ) %>%
      dplyr::mutate(
        Wald = dplyr::case_when(
          pvalue < 0.01 ~ stringr::str_c(Wald, "***"),
          pvalue < 0.05 ~ stringr::str_c(Wald, "**"),
          pvalue < 0.1 ~ stringr::str_c(Wald, "*"),
          TRUE ~ Wald
        )
      ) %>%
      dplyr::mutate(pvalue = base::sprintf("%.4f", pvalue))

    res <- dplyr::bind_rows(res, xx)

    # Test short-run coefficients
    xx <- wald_test(b = coef_short[i, ], cov = covs_short[[i]], R = R)

    xx <- tibble::tibble(
      Variable = stringr::str_c(variables, collapse = "="),
      term = "Short",
      Wald = base::sprintf("%.4f", xx$statistic),
      pvalue = xx$p.value,
      Tau = base::sprintf("%.2f", taus[i])
    ) %>%
      dplyr::mutate(
        Wald = dplyr::case_when(
          pvalue < 0.01 ~ stringr::str_c(Wald, "***"),
          pvalue < 0.05 ~ stringr::str_c(Wald, "**"),
          pvalue < 0.1 ~ stringr::str_c(Wald, "*"),
          TRUE ~ Wald
        )
      ) %>%
      dplyr::mutate(pvalue = base::sprintf("%.4f", pvalue))

    res <- dplyr::bind_rows(res, xx)
  }

  # Export results if requested
  if (export){
    if (base::is.null(file)) {
      stop("File path must be specified when export = TRUE")
    }
    rio::export(res, file)
  }

  return(res)
}
