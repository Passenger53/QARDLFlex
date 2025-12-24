#' Generate Lagged Variables for ARDL Regression
#'
#' This function creates lagged variables (including zero-order) for ARDL model estimation
#' based on specified lag orders for each variable. The resulting data frame contains
#' all necessary lagged terms for ARDL regression analysis.
#'
#' @param data A data frame containing the time series data. If the first column
#'        is named "date", it will be automatically excluded from variable generation.
#' @param order Integer vector specifying the lag order for each variable.
#'        The length must match the number of columns in the data frame (excluding date).
#'        First element corresponds to the dependent variable, subsequent elements to
#'        independent variables. All values must be non-negative integers.
#'
#' @return A tibble containing all lagged variables for ARDL regression, with column
#'         names in the format Ln.variable where n is the lag order (0 for contemporaneous).
#'         The data is ordered as: L0.X, L1.X, ..., L0.Y, L1.Y, ... for each variable.
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
#' # Generate ARDL variables
#' ardl_vars <- generate_ardl_variables(sample_data, order = optimal_order)
#' head(ardl_vars)
generate_ardl_variables <- function(data, order) {
  # Remove date column if first column is named "date"
  if (colnames(data)[1] == "date"){
    dd <- dplyr::select(data, -date)
  } else {
    dd <- data
  }

  if (is.null(names(order))){
    names(order) <- colnames(dd)
  }

  # Validate parameter compatibility
  if (base::ncol(dd) != base::length(order)) {
    stop("Length of 'order' must equal the number of columns in the data frame (excluding date column).")
  }

  # Validate lag orders are non-negative
  if (base::any(order < 0)) {
    stop("All lag orders must be non-negative integers.")
  }

  # Get variable names
  var_names <- base::names(dd)

  # Initialize result container
  result_list <- base::list()

  # Generate lagged variables for each variable
  for (i in base::seq_along(var_names)) {
    var_name <- var_names[i]
    max_lag <- order[i]

    # Create all lag orders for this variable (including lag 0)
    lag_vars <- purrr::map(0:max_lag, function(lag_order) {
      if (lag_order == 0) {
        # Zero-order (contemporaneous) term
        dd[[var_name]]
      } else {
        # Lagged term
        dplyr::lag(dd[[var_name]], lag_order)
      }
    })

    # Create appropriate column names: L0.X, L1.X, ..., Ln.X
    lag_names <- base::paste0("L", 0:max_lag, ".", var_name)
    lag_df <- base::as.data.frame(lag_vars)
    base::colnames(lag_df) <- lag_names

    result_list[[var_name]] <- lag_df
  }

  # Combine all variables maintaining the order: L0.X, L1.X..., L0.Y, L1.Y...
  final_result <- dplyr::bind_cols(result_list)

  # Return as tibble
  return(tibble::as_tibble(final_result))
}
