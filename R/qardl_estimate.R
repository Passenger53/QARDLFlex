#' Estimate QARDL and QECM Models with Bootstrapping
#'
#' This function estimates Quantile Autoregressive Distributed Lag (QARDL) and
#' Quantile Error Correction Model (QECM) with bootstrapped confidence intervals.
#' It provides comprehensive results including bootstrap samples, estimated
#' coefficients with standard errors, significance stars, and pseudo R-squared values.
#'
#' @param data A data frame containing the variables for the model. If the first
#'   column is a date variable, it will be automatically removed from analysis.
#' @param order A numeric vector specifying the ARDL orders. The length should
#'   match the number of non-date variables in the data frame. Typically obtained
#'   from `order_find()`.
#' @param tau The quantile to estimate. Default is 0.5 (median). Value must be
#'   between 0 and 1.
#' @param boots The number of bootstrap replications. Default is 200. Higher
#'   values provide more accurate confidence intervals but increase computation time.
#'
#' @return A list containing six elements:
#'         - QARDL_boot: Bootstrap samples of QARDL coefficients
#'         - QECM_boot: Bootstrap samples of QECM coefficients
#'         - QARDL_coefficient: Formatted QARDL coefficient table with standard
#'           errors and significance stars
#'         - QECM_coefficient: Formatted QECM coefficient table with standard
#'           errors and significance stars, including long-run coefficients
#'         - QARDL_coef: Original QARDL coefficient vector
#'         - QECM_coef: Original QECM coefficient vector
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
#' # Find optimal ARDL orders
#' optimal_order <- order_find(sample_data, selection = "BIC", export = FALSE)
#'
#' # Estimate QARDL and QECM models
#' result <- qardl_estimate(sample_data, optimal_order, tau = 0.5, boots = 50)
#'
#' # View QARDL coefficients
#' print(result$QARDL_coefficient)
#'
#' # View QECM coefficients
#' print(result$QECM_coefficient)
qardl_estimate <- function(data, order, tau = 0.5, boots = 200){
  # Generate lagged variables for QARDL and QECM models
  dd_qardl <- generate_ardl_variables(data, order)
  dd_qecm <- generate_ecm_variables(data, order)

  # Create formulas for QARDL and QECM models
  fm_qardl <- stringr::str_c(base::colnames(dd_qardl)[1], "~",
                             stringr::str_c(base::colnames(dd_qardl)[-1], collapse = "+")) %>%
    stats::as.formula()
  fm_qecm <- stringr::str_c(base::colnames(dd_qecm)[1], "~",
                            stringr::str_c(base::colnames(dd_qecm)[-1], collapse = "+")) %>%
    stats::as.formula()

  # Estimate quantile regression models at specified quantile
  q_ardl <- quantreg::rq(formula = fm_qardl, data = dd_qardl, tau = tau)
  q_ecm <- quantreg::rq(formula = fm_qecm, data = dd_qecm, tau = tau)

  # Generate bootstrap data for inference
  boot_data <- boot_data_generate(data = data,
                                  order = order,
                                  tau = 0.5,
                                  boots = boots)

  # Bootstrap QARDL coefficients: estimate model on each bootstrap sample
  boot_ardl <- purrr::map_dfr(boot_data,
                              ~stats::coef(quantreg::rq(formula = fm_qardl,
                                                        data = generate_ardl_variables(., order),
                                                        tau = tau)))

  # Bootstrap QECM coefficients: estimate model on each bootstrap sample
  boot_ecm <- purrr::map_dfr(boot_data,
                             ~stats::coef(quantreg::rq(formula = fm_qecm,
                                                       data = generate_ecm_variables(., order),
                                                       tau = tau)))

  # Extract variable names for coefficient tables
  names1 <- base::colnames(boot_ardl)
  names2 <- base::colnames(boot_ecm)

  # Format QARDL coefficients: calculate bootstrap standard errors, t-statistics, and p-values
  dplyr::bind_cols(Variable = names1,
                   Estimate = stats::coef(q_ardl),
                   `Std. Error` = base::sqrt(base::diag(stats::cov(boot_ardl)))
  ) %>%
    dplyr::mutate(t = Estimate / `Std. Error`) %>%
    dplyr::mutate(`Pr(>|t|)` = 2 * stats::pt(-base::abs(t), base::nrow(data) - base::length(names1))) %>%
    dplyr::mutate(Estimate = base::sprintf("%.4f", Estimate)) %>%
    dplyr::mutate(Estimate = dplyr::case_when(
      `Pr(>|t|)` < 0.01 ~ stringr::str_c(Estimate, "***"),
      `Pr(>|t|)` < 0.05 ~ stringr::str_c(Estimate, "**"),
      `Pr(>|t|)` < 0.1 ~ stringr::str_c(Estimate, "*"),
      TRUE ~ Estimate
    )) %>%
    dplyr::mutate(`Std. Error` = base::sprintf("%.4f", `Std. Error`)) %>%
    dplyr::mutate(`Std. Error` = stringr::str_c("(", `Std. Error`, ")")) %>%
    dplyr::select(Variable, Estimate, `Std. Error`) %>%
    tidyr::pivot_longer(cols = -Variable,
                        names_to = "Coef (Std. Error)",
                        values_to = "Value") -> coef_qardl

  # Fit an intercept-only model for pseudo R-squared calculation
  f1 <- stringr::str_c(dd_qardl %>% base::colnames() %>% .[1], "~1") %>% stats::as.formula()
  f2 <- stringr::str_c(dd_qecm %>% base::colnames() %>% .[1], "~1") %>% stats::as.formula()
  q_ardl1 <- quantreg::rq(formula = f1, data = dd_qardl, tau = tau)
  q_ecm1 <- quantreg::rq(formula = f2, data = dd_qecm, tau = tau)

  # Calculate pseudo R-squared: 1 - (model objective function / intercept-only objective function)
  R2_qardl <- 1 - q_ardl$rho / q_ardl1$rho
  R2_qecm <- 1 - q_ecm$rho / q_ecm1$rho

  # Add pseudo R-squared to QARDL coefficient table
  tibble::tibble(Variable = "R2",
                 `Coef (Std. Error)` = NA,
                 Value = R2_qardl %>% base::sprintf("%.4f", .)) -> R2_qardl
  coef_qardl <- dplyr::bind_rows(coef_qardl, R2_qardl)

  # Process QECM coefficients: identify ECM term and calculate long-run coefficients
  find <- stringr::str_detect(names2, "L1.") %>% base::which()
  names2[find[1]] <- "ECM"
  coef_qecm <- stats::coef(q_ecm)
  lambda <- coef_qecm[find[1]]

  # Calculate long-run coefficients: -β/λ
  coef_qecm[find[-1]] <- - coef_qecm[find[-1]] / lambda
  for (i in 2:base::length(find)){
    boot_ecm[, find[i]] <- - boot_ecm[, find[i]] / boot_ecm[, find[1]]
  }

  # Format QECM coefficients: calculate bootstrap standard errors, t-statistics, and p-values
  dplyr::bind_cols(Variable = names2,
                   Estimate = coef_qecm,
                   `Std. Error` = base::sqrt(base::diag(stats::cov(boot_ecm)))
  ) %>%
    dplyr::mutate(t = Estimate / `Std. Error`) %>%
    dplyr::mutate(`Pr(>|t|)` = 2 * stats::pt(-base::abs(t), base::nrow(data) - base::length(coef_qecm))) %>%
    dplyr::mutate(Estimate = base::sprintf("%.4f", Estimate)) %>%
    dplyr::mutate(Estimate = dplyr::case_when(
      `Pr(>|t|)` < 0.01 ~ stringr::str_c(Estimate, "***"),
      `Pr(>|t|)` < 0.05 ~ stringr::str_c(Estimate, "**"),
      `Pr(>|t|)` < 0.1 ~ stringr::str_c(Estimate, "*"),
      TRUE ~ Estimate
    )) %>%
    dplyr::mutate(`Std. Error` = base::sprintf("%.4f", `Std. Error`)) %>%
    dplyr::mutate(`Std. Error` = stringr::str_c("(", `Std. Error`, ")")) %>%
    dplyr::select(Variable, Estimate, `Std. Error`) %>%
    tidyr::pivot_longer(cols = -Variable,
                        names_to = "Coef (Std. Error)",
                        values_to = "Value") -> coef_qecm

  # Add pseudo R-squared to QECM coefficient table
  tibble::tibble(Variable = "R2",
                 `Coef (Std. Error)` = NA,
                 Value = R2_qecm %>% base::sprintf("%.4f", .)) -> R2_qecm
  coef_qecm <- dplyr::bind_rows(coef_qecm, R2_qecm)

  # Return comprehensive results
  return(
    base::list(
      QARDL_boot = boot_ardl,
      QECM_boot = boot_ecm,
      QARDL_coefficient = coef_qardl,
      QECM_coefficient = coef_qecm,
      QARDL_coef = stats::coef(q_ardl),
      QECM_coef = stats::coef(q_ecm)
    )
  )
}
