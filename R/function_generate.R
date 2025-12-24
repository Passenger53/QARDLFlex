#' Generate Formula from Data Frame Columns
#'
#' This function automatically creates a formula object from the column names of a data frame.
#' It uses the first non-date column as the dependent variable and the remaining columns
#' as independent variables in the formula.
#'
#' @param data A data frame containing the variables for formula generation.
#'        The function expects that the first column might be a "date" column (which is excluded),
#'        and uses the second column as the dependent variable, with subsequent columns as independent variables.
#'
#' @return A formula object where the first variable (after excluding the first column) is the
#'         dependent variable, and all remaining variables are independent variables joined by "+".
#'         The formula has the structure: `second_column ~ third_column + fourth_column + ...`
#'
#' @export
#'
#' @examples
#' # Generate sample data using sim_data_generate()
#' sample_data <- sim_data_generate(
#'   b = c(0.5, 0.7, 0.3, 0.2),
#'   order = c(1, 1),
#'   size = 100,
#'   reps = 1
#' )[[1]]
#'
#' # Create formula from the data
#' model_formula <- function_generate(sample_data)
#' print(model_formula)
#' # Output: y ~ x1 + x2  (assuming the data has columns: date, y, x1, x2)
function_generate <- function(data){
  # Extract all column names from the data frame
  xx <- colnames(data)

  # Remove the first column (typically "date" column) from consideration
  xx <- xx[-1]

  # Create formula string: first remaining variable as dependent variable,
  # all other variables as independent variables joined by "+"
  ff <- stringr::str_c(xx[1], "~", stringr::str_c(xx[-1], collapse = "+"))

  # Convert the string to a formula object
  ff <- stats::as.formula(ff)

  return(ff)
}
