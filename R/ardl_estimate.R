#' Estimate ARDL and ECM Models
#'
#' This function estimates both ARDL (Autoregressive Distributed Lag) and ECM (Error Correction Model)
#' specifications based on the provided data and lag orders. It returns comprehensive estimation results
#' including model objects, summary statistics, and formatted coefficient tables with significance stars.
#'
#' @param data A data frame containing the time series data. If the first column
#'        is named "date", it will be automatically excluded from the analysis.
#' @param order Integer vector. ARDL lag orders for each variable, typically obtained
#'        from `order_find()`. The length must match the number of variables in the
#'        data frame (excluding date).
#'
#' @return A list containing six elements:
#'         - ARDL_estimate: The full ARDL model object from `stats::lm()`
#'         - ECM_estimate: The full ECM model object from `stats::lm()`
#'         - ARDL_summary: Summary of the ARDL model
#'         - ECM_summary: Summary of the ECM model
#'         - ARDL_coefficient: Formatted ARDL coefficient table with significance stars
#'         - ECM_coefficient: Formatted ECM coefficient table with significance stars
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
#' # Estimate ARDL and ECM models
#' results <- ardl_estimate(sample_data, order = optimal_order)
#'
#' # Access ARDL coefficients
#' print(results$ARDL_coefficient)
#'
#' # Access ECM coefficients
#' print(results$ECM_coefficient)
#'
#' # Get full ARDL model summary
#' summary(results$ARDL_estimate)
ardl_estimate <- function(data, order){
  # Generate lagged variables for ARDL model
  dd_ardl <- generate_ardl_variables(data, order)
  # Generate lagged and differenced variables for ECM model
  dd_ecm <- generate_ecm_variables(data, order)

  # Create formula for ARDL model: first variable ~ all other variables
  fm_ardl <- stringr::str_c(base::colnames(dd_ardl)[1], "~",
                            stringr::str_c(base::colnames(dd_ardl)[-1], collapse = "+")) %>%
    stats::as.formula()

  # Create formula for ECM model: first variable ~ all other variables
  fm_ecm <- stringr::str_c(base::colnames(dd_ecm)[1], "~",
                           stringr::str_c(base::colnames(dd_ecm)[-1], collapse = "+")) %>%
    stats::as.formula()

  # Estimate ARDL model using ordinary least squares
  lm_ardl <- stats::lm(formula = fm_ardl, data = dd_ardl)
  # Estimate ECM model using ordinary least squares
  lm_ecm <- stats::lm(formula = fm_ecm, data = dd_ecm)

  # Get model summaries
  sum_ardl <- base::summary(lm_ardl)
  sum_ecm <- base::summary(lm_ecm)

  # Extract coefficient tables
  coef_ardl <- sum_ardl$coefficients
  coef_ecm <- sum_ecm$coefficients

  # Get coefficient names
  names1 <- base::rownames(coef_ardl)
  names2 <- base::rownames(coef_ecm)

  # Format ARDL coefficients: add significance stars and format standard errors
  dplyr::bind_cols(tibble::tibble(Variable = names1),
                   tibble::as_tibble(coef_ardl)) %>%
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
                        values_to = "Value") -> coef_ardl

  # Add R-squared to ARDL coefficient table
  tibble::tibble(Variable = "R2",
                 `Coef (Std. Error)` = NA,
                 Value = sum_ardl$r.squared %>% base::sprintf("%.4f", .)) -> R2_ardl
  coef_ardl <- dplyr::bind_rows(coef_ardl, R2_ardl)

  # Process ECM coefficients: identify ECM term and calculate long-run coefficients
  find <- stringr::str_detect(names2, "L1.") %>% base::which()
  base::rownames(coef_ecm)[find[1]] <- "ECM"
  lambda <- coef_ecm[find[1], 1]

  # Calculate long-run coefficients for lagged independent variables
  coef_ecm[find[-1], 1] <- - coef_ecm[find[-1], 1] / lambda
  coef_ecm[find[-1], 2] <- - coef_ecm[find[-1], 2] / lambda

  # Format ECM coefficients: add significance stars and format standard errors
  dplyr::bind_cols(tibble::tibble(Variable = base::rownames(coef_ecm)),
                   tibble::as_tibble(coef_ecm)) %>%
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
                        values_to = "Value") -> coef_ecm

  # Add R-squared to ECM coefficient table
  tibble::tibble(Variable = "R2",
                 `Coef (Std. Error)` = NA,
                 Value = sum_ecm$r.squared %>% base::sprintf("%.4f", .)) -> R2_ecm
  coef_ecm <- dplyr::bind_rows(coef_ecm, R2_ecm)

  # Return comprehensive results as a named list
  return(
    base::list(
      ARDL_estimate = lm_ardl,
      ECM_estimate = lm_ecm,
      ARDL_summary = sum_ardl,
      ECM_summary = sum_ecm,
      ARDL_coefficient = coef_ardl,
      ECM_coefficient = coef_ecm
    )
  )
}
