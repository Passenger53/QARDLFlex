#' Calculate Dynamic Multipliers for All Explanatory Variables in ARDL Model
#'
#' This function computes dynamic multipliers for all explanatory variables in an
#' ARDL model across specified quantile levels. It extends the `dynamic_multipliers()`
#' function to automatically handle multiple variables and quantiles, providing
#' a comprehensive analysis of how exogenous shocks propagate through the system
#' over time.
#'
#' @param data A data frame containing the time series data. If the first column
#'   is a date variable, it will be automatically excluded from analysis.
#' @param order Integer vector. ARDL lag orders for each variable, typically
#'   obtained from `order_find()`. The length must match the number of variables
#'   in the data frame (excluding date).
#' @param taus Numeric vector. Quantile levels to estimate. All values must be
#'   between 0 and 1. Common values include c(0.25, 0.5, 0.75) for quartiles.
#' @param steps Integer. Number of time periods after the shock to calculate
#'   multipliers for. Must be a positive integer greater than 0.
#' @param cumulative Logical. If TRUE, returns cumulative multipliers (the sum
#'   of multipliers up to each period). If FALSE (default), returns period-by-period
#'   multipliers.
#' @param impulse Logical. Specifies the type of shock:
#'   - FALSE (default): Permanent step shock (x_t = 1 for all t ≥ 0)
#'   - TRUE: Temporary impulse shock (x_0 = 1, x_t = 0 for t > 0)
#' @param export Logical. Whether to export the results to a file. Default is TRUE.
#' @param file Character. File path for exporting results. Required if export = TRUE.
#'
#' @return A tibble containing dynamic multipliers for all variables, with columns:
#'   - variable: Name of the explanatory variable
#'   - tau: Quantile level at which the multiplier was calculated
#'   - step: Time period after the shock (0 = impact period, 1 = one period after, etc.)
#'   - value: The dynamic multiplier value
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate sample data using sim_data_generate()
#' sample_data <- sim_data_generate(
#'   b = c(0.5, 0.7, 0.3, 0.2),
#'   order = c(1, 1),
#'   size = 100,
#'   reps = 1
#' )[[1]]
#'
#' # Find optimal lag order using order_find()
#' optimal_order <- order_find(sample_data, selection = "BIC", export = FALSE)
#'
#' # Calculate dynamic multipliers for all explanatory variables
#' results <- all_dynamic_multipliers(
#'   data = sample_data,
#'   order = optimal_order,
#'   taus = c(0.25, 0.5, 0.75),
#'   steps = 10,
#'   cumulative = FALSE,
#'   impulse = FALSE,
#'   export = FALSE
#' )
#'
#' # View results
#' print(results)
#' }
all_dynamic_multipliers <- function(data, order, taus, steps, cumulative = FALSE,
                                    impulse = FALSE, export = TRUE, file = NULL){

  # Generate ARDL variables from the original data
  dd <- generate_ardl_variables(data, order)

  # Create formula for regression: dependent_var ~ var1 + var2 + ...
  stringr::str_c(base::names(dd)[1], "~",
                 stringr::str_c(base::names(dd)[-1], collapse = "+")) %>%
    stats::as.formula() -> fomu

  # Estimate OLS and quantile regression models
  ols <- stats::coef(stats::lm(fomu, dd))
  qr <- stats::coef(quantreg::rq(formula = fomu, data = dd, tau = taus))
  coefs <- base::cbind(ols, qr)

  # Calculate cumulative positions for coefficient extraction
  find <- base::cumsum(order + 1)

  # Extract autoregressive coefficients (phi) for the dependent variable
  phis <- coefs[2:find[1], ]
  # Handle case when there's only one AR coefficient
  if (find[1] == 2){
    phis <- base::matrix(phis, nrow = 1)
  }

  # Initialize lists for beta coefficients and multipliers
  betas <- base::vector("list", base::length(order) - 1)
  multipliers <- base::vector("list", base::length(order) - 1)

  # Process each explanatory variable
  for (i in 1:(base::length(order) - 1)){
    # Extract coefficients for current variable
    begin <- find[i] + 1
    end <- find[i+1]
    betas[[i]] <- coefs[begin:end, ]

    # Ensure betas is a matrix even for single coefficient cases
    if (base::is.null(base::dim(betas[[i]]))){
      betas[[i]] <- base::matrix(betas[[i]], nrow = 1)
    }

    # Initialize matrix for storing multipliers for current variable
    ms <- NULL

    # Calculate dynamic multipliers for each quantile level
    for (j in 1:base::ncol(coefs)){
      phi <- phis[, j] %>% base::as.vector()
      beta <- betas[[i]][, j] %>% base::as.vector()

      # Calculate dynamic multipliers using the core function
      m0 <- dynamic_multipliers(
        phi = phi,
        beta = beta,
        steps = steps,
        cumulative = cumulative,
        impulse = impulse
      )
      ms <- base::cbind(ms, m0)
    }

    # Set column names and store results
    base::colnames(ms) <- base::colnames(coefs)
    multipliers[[i]] <- ms
  }

  # Name the multiplier list with explanatory variable names
  base::names(multipliers) <- base::names(order)[-1]

  # Reshape results from wide to long format for better analysis
  multipliers %>%
    reshape2::melt() %>%
    tibble::as_tibble() %>%
    dplyr::rename(step = Var1,
                  tau = Var2,
                  variable = L1) %>%
    dplyr::mutate(step = base::as.character(step),
                  tau = base::as.character(tau),
                  variable = base::as.character(variable)) %>%
    dplyr::mutate(step = readr::parse_number(step)) %>%
    dplyr::select(variable, tau, step, value) -> multipliers

  # Export results if requested
  if (export){
    if (base::is.null(file)) {
      stop("File path must be specified when export = TRUE")
    }
    rio::export(multipliers, file)
  }

  return(multipliers)
}
