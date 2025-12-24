#' Generate Variables for Error Correction Model (ECM) Regression
#'
#' This function creates lagged and differenced variables for Error Correction Model
#' estimation based on ARDL lag orders. It handles cases where some variables have
#' zero lag order (contemporaneous only) differently from those with positive lag orders.
#'
#' @param data A data frame containing the time series data. If the first column
#'        is named "date", it will be automatically excluded from variable generation.
#' @param order Integer vector specifying the lag order for each variable.
#'        The length must match the number of columns in the data frame (excluding date).
#'        Values can be 0 (contemporaneous only) or positive integers.
#'
#' @return A tibble containing variables for ECM regression, including:
#'         - L1.*: First lag of all variables
#'         - D0.*: First difference of dependent variable
#'         - Dn.*: First differences of lagged variables
#'         - L0.*: Contemporaneous variables with zero lag order (if any)
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
#' # Generate ECM variables
#' ecm_vars <- generate_ecm_variables(sample_data, order = optimal_order)
#' head(ecm_vars)
generate_ecm_variables <- function(data, order){
  # Remove date column if first column is named "date"
  if (colnames(data)[1] == "date"){
    dd <- dplyr::select(data, -date)
  } else {
    dd <- data
  }

  if (is.null(names(order))){
    names(order) <- colnames(dd)
  }

  # Case 1: Some variables have zero lag order (contemporaneous only)
  if (base::any(order == 0)){
    # Identify variables with zero lag order
    find <- base::which(order == 0)
    order1 <- order[-find]
    names_zero <- base::names(order)[find]

    # Create subset of data without zero-lag variables
    data1 <- dplyr::select(data, -names_zero)
    # Reduce lag orders by 1 for remaining variables
    order1 <- order1 - 1

    # Generate ARDL variables for non-zero lag variables
    ecm_data <- generate_ardl_variables(data = data1, order = order1)

    # Calculate first differences (ECM transformation)
    ecm_data <- purrr::map_dfc(ecm_data, function(x){
      base::c(NA, base::diff(x))
    })

    # Rename columns: change L to D for differenced variables
    var_names <- base::colnames(ecm_data)
    var_names <- base::substr(var_names, 2, base::nchar(var_names))
    var_names <- stringr::str_c("D", var_names)
    base::colnames(ecm_data) <- var_names

    # Create first lag of all non-zero-lag variables
    dd1 <- dd %>% dplyr::select(-names_zero)
    dplyr::lag(dd1) -> ecm_data0
    base::colnames(ecm_data0) <- stringr::str_c("L1.", base::colnames(ecm_data0))

    # Combine lagged and differenced variables
    ecm_data <- dplyr::bind_cols(ecm_data0, ecm_data)

    # Move dependent variable difference to first column
    ecm_data <- ecm_data %>%
      dplyr::relocate(stringr::str_c("D0.", base::colnames(dd)[1]), .before = 1)

    # Add back zero-lag variables as L0.*
    xx <- dd %>% dplyr::select(names_zero)
    var_names <- base::colnames(xx)
    var_names <- stringr::str_c("L0.", var_names)
    base::colnames(xx) <- var_names

    ecm_data <- dplyr::bind_cols(ecm_data, xx)
    return(ecm_data)
  }

  # Case 2: All variables have positive lag orders
  if (base::all(order > 0)){
    # Reduce all lag orders by 1
    order1 <- order - 1

    # Generate ARDL variables
    ecm_data <- generate_ardl_variables(data = data, order = order1)

    # Calculate first differences
    ecm_data <- purrr::map_dfc(ecm_data, function(x){
      base::c(NA, base::diff(x))
    })

    # Rename columns: change L to D for differenced variables
    var_names <- base::colnames(ecm_data)
    var_names <- base::substr(var_names, 2, base::nchar(var_names))
    var_names <- stringr::str_c("D", var_names)
    base::colnames(ecm_data) <- var_names

    # Create first lag of all variables
    dplyr::lag(dd) -> ecm_data0
    base::colnames(ecm_data0) <- stringr::str_c("L1.", base::colnames(ecm_data0))

    # Combine lagged and differenced variables
    ecm_data <- dplyr::bind_cols(ecm_data0, ecm_data)

    # Move dependent variable difference to first column
    ecm_data <- ecm_data %>%
      dplyr::relocate(stringr::str_c("D0.", base::colnames(dd)[1]), 1) -> ecm_data

    return(ecm_data)
  }
}
