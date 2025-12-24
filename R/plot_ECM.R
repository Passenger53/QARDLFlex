#' Plot ECM Model Coefficient Trends Across Quantiles
#'
#' This function plots the long-term and short-term coefficient trends for
#' ECM model variables across specified quantile levels. For short-term
#' coefficients, if a variable has multiple lags, the sum of coefficients
#' across all lags is plotted.
#'
#' @param estimation List. Estimation results obtained from `coef_estimate()`.
#' @param order Integer vector. ARDL orders obtained from `order_find()`.
#'   Length must match the number of variables, otherwise an error is raised.
#' @param plot_term Character. Specify whether to plot "short" (short-term)
#'   or "long" (long-term) coefficients. Only "short" and "long" are allowed.
#' @param variable Character. Name of the variable to plot. Must match one
#'   of the variable names in the data. Special cases:
#'   - If the variable is the dependent variable and plot_term = "long",
#'     the ECM coefficient is plotted by default.
#'   - If the variable is an independent variable with ARDL order of zero,
#'     its contemporaneous effect coefficient is plotted.
#' @param taus Numeric vector. Quantile levels to estimate, same as the
#'   `taus` parameter in `coef_estimate()`.
#' @param ci Logical. Whether to plot confidence intervals. Default is TRUE.
#' @param level Numeric. Confidence level. Only 0.9, 0.95, and 0.99 are allowed.
#'   Default is 0.95.
#' @param export Logical. Whether to export the plot to disk. Default is TRUE.
#' @param file Character. File path for exporting the plot. Required if export = TRUE.
#'
#' @return A list containing two elements:
#'   - plot: A ggplot object of the coefficient trend plot
#'   - plot_data: The data used for plotting
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
#' # Estimate OLS and quantile regression models
#' results <- coef_estimate(
#'   data = sample_data,
#'   order = optimal_order,
#'   taus = c(0.25, 0.5, 0.75),
#'   boots = 50,
#'   export = FALSE
#' )
#'
#' # Plot long-term coefficient for the first variable
#' plot_result <- plot_ECM(
#'   estimation = results,
#'   order = optimal_order,
#'   plot_term = "long",
#'   variable = names(optimal_order)[1],
#'   taus = c(0.25, 0.5, 0.75),
#'   ci = TRUE,
#'   level = 0.95,
#'   export = FALSE
#' )
#'
#' # Display the plot
#' print(plot_result$plot)
plot_ECM <- function(
    estimation,
    order,
    plot_term,
    variable,
    taus,
    ci = TRUE,
    level = 0.95,
    export = TRUE,
    file = NULL
){

  # Validate plot_term parameter: only "short" or "long" are allowed
  if (!(plot_term %in% c("short", "long", "none"))){
    stop("The plot_term parameter can only be 'short' (short-term) or 'long' (long-term). Please check carefully!")
  }

  # Validate that variable exists in the order vector
  if (!(variable %in% base::names(order))){
    stop("The variable parameter does not match any variable name in the data. Please check carefully!")
  }

  # Validate confidence level if confidence intervals are requested
  if (ci){
    if (!(level %in% c(0.9, 0.95, 0.99))){
      stop("The level parameter can only be set to 0.9, 0.95, or 0.99!")
    }
  }

  # Check for special case: dependent variable with lag order 1 has no short-term coefficient
  if (variable == base::names(order)[1] & plot_term == "short" & order[1] == 1){
    stop("Variable ", variable, " is the dependent variable with lag order 1, and has no short-term coefficient!")
  }

  # Extract OLS coefficients from ECM estimation
  ols <- stats::coef(estimation$OLS$ECM_estimate)

  # Initialize containers for quantile regression coefficients and bootstrap samples
  qr_coef <- NULL
  qr_boot <- NULL

  # Extract coefficients and bootstrap samples for each quantile level
  for (i in 1:base::length(estimation$QR)){
    # Combine quantile regression coefficients
    qr_coef <- base::rbind(qr_coef, estimation$QR[[i]]$QECM_coef)

    # Process bootstrap samples: add quantile and bootstrap iteration identifiers
    estimation$QR[[i]]$QECM_boot %>%
      dplyr::mutate(tau = taus[i]) %>%
      dplyr::relocate(tau, 1) %>%
      dplyr::mutate(boot = 1:dplyr::n()) %>%
      dplyr::relocate(boot, .after = tau) %>%
      tidyr::pivot_longer(cols = -c(tau, boot),
                          names_to = "var",
                          values_to = "boot_coef") -> xx
    qr_boot <- dplyr::bind_rows(qr_boot, xx)
  }

  # Format quantile regression coefficients: add tau column and reshape to long format
  qr_coef %>% tibble::as_tibble() %>% dplyr::mutate(tau = taus) %>%
    dplyr::relocate(tau, 1) -> qr_coef
  qr_coef %>%
    tidyr::pivot_longer(-tau,
                        names_to = "var",
                        values_to = "coef") -> qr_coef

  # Join bootstrap samples with point estimates
  dplyr::left_join(qr_boot, qr_coef, by = dplyr::join_by(tau, var)) -> qr

  # Classify coefficients as long-term (L1.*), short-term (D*), or none
  # Remove intercept, extract variable names, and aggregate coefficients
  qr %>%
    dplyr::filter(var != "(Intercept)") %>%
    dplyr::mutate(term = dplyr::case_when(
      base::substring(var, 1, 2) == "L1" ~ "long",
      base::substring(var, 1, 1) == "D" ~ "short",
      TRUE ~ "none"
    )) %>%
    dplyr::relocate(term, .after = tau) %>%
    dplyr::mutate(var = base::substring(var, 4, base::nchar(var))) %>%
    dplyr::summarise(boot_coef = base::sum(boot_coef),
                     coef = base::sum(coef),
                     .by = c(tau, term, boot, var)) -> qr

  # Adjust bootstrap coefficients: center around mean of point estimates
  qr %>%
    dplyr::mutate(boot_coef = boot_coef - base::mean(boot_coef) + base::mean(coef),
                  .by = c(tau, term, var)) -> qr

  # If confidence intervals are requested, calculate and plot with error bands
  if (ci){
    # Determine critical value based on confidence level
    mm <- dplyr::case_when(
      level == 0.9 ~ 1.645,
      level == 0.95 ~ 1.960,
      level == 0.99 ~ 2.576,
      TRUE ~ NA
    )

    # Calculate confidence intervals from bootstrap samples
    qr %>%
      dplyr::summarise(sd = stats::sd(boot_coef),
                       .by = c(tau, term, var, coef)) %>%
      dplyr::mutate(up = coef + mm * sd,
                    lo = coef - mm * sd) %>%
      dplyr::select(-sd) -> qr

    # Case 1: Dependent variable with long-term coefficient (ECM coefficient)
    if (variable == base::names(order)[1] & plot_term == "long"){
      base::message("Variable ", variable, " is the dependent variable. Its long-term coefficient is actually the ECM coefficient.")

      qr %>%
        dplyr::filter(var == variable & term == term) %>%
        ggplot2::ggplot(ggplot2::aes(tau, coef)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(shape = 7) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = up), alpha = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Quantiles", y = "Coefficients") -> p

    }
    # Case 2: Independent variable with zero lag order
    else if (order[base::which(base::names(order) == variable)] == 0){
      base::message("The variable has ARDL order zero. The plot_term parameter is ignored. Plotting the contemporaneous effect of ", variable, " on the dependent variable.")

      qr %>%
        dplyr::filter(var == variable) %>%
        ggplot2::ggplot(ggplot2::aes(tau, coef)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(shape = 7) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = up), alpha = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Quantiles", y = "Coefficients") -> p

    }
    # Case 3: General case
    else {
      qr %>%
        dplyr::filter(var == variable & term == plot_term) %>%
        ggplot2::ggplot(ggplot2::aes(tau, coef)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(shape = 7) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = up), alpha = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Quantiles", y = "Coefficients") -> p
    }
  }
  # If no confidence intervals, plot without error bands
  else {
    if (variable == base::names(order)[1] & plot_term == "long"){
      base::message("Variable ", variable, " is the dependent variable. Its long-term coefficient is actually the ECM coefficient.")

      qr %>%
        dplyr::filter(var == variable & term == term) %>%
        ggplot2::ggplot(ggplot2::aes(tau, coef)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(shape = 7) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Quantiles", y = "Coefficients") -> p

    } else if (order[base::which(base::names(order) == variable)] == 0){
      base::message("The variable has ARDL order zero. The plot_term parameter is ignored. Plotting the contemporaneous effect of ", variable, " on the dependent variable.")

      qr %>%
        dplyr::filter(var == variable) %>%
        ggplot2::ggplot(ggplot2::aes(tau, coef)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(shape = 7) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Quantiles", y = "Coefficients") -> p

    } else {
      qr %>%
        dplyr::filter(var == variable & term == plot_term) %>%
        ggplot2::ggplot(ggplot2::aes(tau, coef)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(shape = 7) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Quantiles", y = "Coefficients") -> p
    }
  }

  # Export the plot if requested
  if (export){
    if (base::is.null(file)) {
      stop("File path must be specified when export = TRUE")
    }
    ggplot2::ggsave(file, p)
  }

  # Print the plot
  base::print(p)

  # Return plot and data
  return(base::list(plot = p,
                    plot_data = qr))
}
