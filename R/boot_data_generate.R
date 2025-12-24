#' Generate Bootstrap Data for Quantile ARDL Model Analysis
#'
#' This function creates bootstrap datasets for quantile ARDL model analysis by
#' resampling residuals from a quantile regression model and generating new
#' time series data based on the estimated coefficients and ARDL structure.
#'
#' @param data A data frame containing the time series data. If the first column
#'        is named "date", it will be automatically excluded from the bootstrap
#'        generation process.
#' @param order Integer vector. ARDL lag orders for each variable, typically obtained
#'        from `order_find()`. The length must match the number of variables in the
#'        data frame (excluding date).
#' @param tau Numeric. The quantile level for the quantile regression, between 0 and 1.
#'        Common values are 0.25, 0.5 (median), and 0.75.
#' @param boots Integer. The number of bootstrap samples to generate.
#'
#' @return A list of length `boots` containing bootstrap datasets. Each dataset
#'         is a tibble with the same structure as the original data, including
#'         the date column (if present).
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
#' # Generate bootstrap data for median (tau = 0.5) regression
#' boot_data <- boot_data_generate(sample_data, order = optimal_order,
#'                                 tau = 0.5, boots = 10)
#'
#' # Check structure of first bootstrap sample
#' str(boot_data[[1]])
boot_data_generate <- function(data, order, tau, boots){
  # Generate ARDL variables for quantile regression
  dd_qardl <- generate_ardl_variables(data, order)

  # Create formula for quantile ARDL model
  fm_qardl <- stringr::str_c(base::colnames(dd_qardl)[1], "~",
                             stringr::str_c(base::colnames(dd_qardl)[-1], collapse = "+")) %>%
    stats::as.formula()

  # Estimate quantile ARDL model at specified quantile level
  q_ardl <- quantreg::rq(formula = fm_qardl, data = dd_qardl, tau = tau)

  # Extract coefficients and residuals from quantile regression
  b <- stats::coef(q_ardl)
  epsilon <- stats::residuals(q_ardl)

  # Select only the variables used in the model (excluding date if present)
  dd <- data %>% dplyr::select(base::names(order))

  # Initialize list to store bootstrap datasets
  boot_data <- base::vector("list", boots)
  rows <- base::nrow(data)
  max_lag <- base::max(order)

  # Set seed for reproducible bootstrap sampling
  base::set.seed(123)
  # Sample residuals with replacement to create bootstrap errors
  boot_e <- base::sample(epsilon, size = rows * boots, replace = TRUE)
  boot_e %>% base::matrix(nrow = rows) -> boot_e

  # Create location matrix to track variable indices and lag orders
  # This matrix maps each parameter to its corresponding variable and lag
  location <- NULL
  for (i in 1:base::length(order)){
    if (i == 1){
      # For dependent variable: include lags 1 through order[i]
      xx <- base::cbind(i, 1:order[i])
      base::colnames(xx) <- base::c("col", "row")
      location <- base::rbind(location, xx)
    } else {
      # For independent variables: include contemporaneous (0) and lags
      xx <- base::cbind(i, 0:order[i])
      base::colnames(xx) <- base::c("col", "row")
      location <- base::rbind(location, xx)
    }
  }

  # Generate bootstrap datasets
  for (i in 1:boots){
    # Start with original data
    boot_data[[i]] <- dd
    # Set dependent variable values to NA for periods after max lag
    boot_data[[i]][(max_lag + 1):rows, 1] <- NA
    # Convert to matrix for efficient computation
    boot_data[[i]] <- base::as.matrix(boot_data[[i]])
    # Get bootstrap errors for this replication
    ee <- boot_e[, i]

    # Generate ARDL process for each time period after max lag
    for (j in (base::max(order) + 1):rows){
      values <- base::rep(NA, base::nrow(location))

      # Collect lagged values according to location matrix
      for (k in 1:base::nrow(location)){
        rr <- location[k, 2]  # Lag order
        cc <- location[k, 1]  # Variable index
        values[k] <- boot_data[[i]][j - rr, cc]
      }

      # Add intercept and calculate new dependent variable value
      # y_t = b0 + b1*y_{t-1} + ... + b_k*x_{t-l} + bootstrap_error
      values <- base::c(1, values)  # Include intercept
      boot_data[[i]][j, 1] <- base::sum(b * values) + ee[j]
    }

    # Convert back to tibble and add date column
    boot_data[[i]] <- tibble::as_tibble(boot_data[[i]])
    boot_data[[i]] %>%
      dplyr::mutate(date = data$date) %>%
      dplyr::relocate(date, 1) -> boot_data[[i]]
  }

  return(boot_data)
}
