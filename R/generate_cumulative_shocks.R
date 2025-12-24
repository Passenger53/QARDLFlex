#' Generate Cumulative Positive and Negative Shock Variables for NARDL Models
#'
#' This function creates cumulative positive and negative shock variables from
#' original time series data, which are essential for Nonlinear Autoregressive
#' Distributed Lag (NARDL) model estimation. It calculates the cumulative sums
#' of positive and negative changes for specified variables.
#'
#' @param variable_names Character vector. Names of the variables in the data
#'   frame for which to generate cumulative shock variables. Each variable must
#'   exist in the provided data frame.
#' @param data A data frame containing the time series data. Should include
#'   all variables specified in variable_names as columns.
#'
#' @return A modified data frame with additional columns for cumulative positive
#'   and negative shocks. For each original variable, two new columns are added:
#'   - p_{variable_name}: Cumulative positive shocks (increases from previous period)
#'   - n_{variable_name}: Cumulative negative shocks (decreases from previous period)
#'   The new columns are placed immediately after their corresponding original variables.
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
#' # Generate cumulative shocks for specific variables
#' processed_data <- generate_cumulative_shocks(
#'   variable_names = c("x1", "x2"),  # Replace with actual variable names
#'   data = sample_data
#' )
#'
#' # View the modified data structure
#' str(processed_data)
#' }
generate_cumulative_shocks <- function(variable_names, data){
  # Validate that at least one variable name is provided
  if (base::length(variable_names) >= 1){
    # Process each specified variable
    for (i in 1:base::length(variable_names)){
      # Calculate first differences of the variable
      xx <- base::c(0, base::diff(data[[variable_names[i]]]))

      # Calculate cumulative positive shocks (increases from previous period)
      positive <- base::cumsum(base::ifelse(xx > 0, xx, 0))

      # Calculate cumulative negative shocks (decreases from previous period)
      negative <- base::cumsum(base::ifelse(xx < 0, xx, 0))

      # Add new columns to data frame
      data[stringr::str_c("p_", variable_names[i])] <- positive
      data[stringr::str_c("n_", variable_names[i])] <- negative

      # Reorder columns to place shock variables next to original variable
      data %>%
        dplyr::relocate(dplyr::contains(stringr::str_c("p_", variable_names[i])),
                        dplyr::contains(stringr::str_c("n_", variable_names[i])),
                        .after = variable_names[i]) -> data
    }
  }
  return(data)
}
