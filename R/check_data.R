#' Check Data Frame for Modeling Requirements
#'
#' This function validates that a data frame meets the basic requirements for
#' quantile ARDL model analysis. It checks for: correct first column name,
#' absence of spaces in variable names, appropriate data types, and missing values.
#'
#' @param data A data frame. The data object to be checked, where the first column
#' should be a date column named "date" (all lowercase).
#'
#' @return Invisibly returns `TRUE` if all checks pass. This function primarily
#' operates through side effects, generating messages, warnings, or errors.
#' If data passes all checks, a "Data validation passed!" message is displayed.
#' If potential issues are found (non-numeric variables, missing values), warnings are issued.
#' If critical errors are found (incorrect first column name, spaces in variable names),
#' errors are thrown via `stop()`.
#'
#' @export
#'
#' @examples
#' # Create a sample data frame that meets requirements
#' good_data <- data.frame(
#'   date = seq.Date(from = as.Date("2000-01-01"), by = "month", length.out = 24),
#'   y = rnorm(24),
#'   x1 = rnorm(24),
#'   x2 = rnorm(24)
#' )
#' check_data(good_data)
check_data <- function(data) {

  # Validate input is a data frame
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame.")
  }

  # Check 1: First column name must be "date" (all lowercase)
  if (colnames(data)[1] != "date") {
    stop("First column name is '", colnames(data)[1], "', but should be 'date'! All lowercase letters required.")
  }

  # Check 2: Variable names should not contain spaces
  if (any(stringr::str_detect(colnames(data), " "))) {
    stop("Variable names contain spaces. Please rename variables without spaces!")
  }

  # Check 3: All columns except date should be numeric (numeric or integer)
  # Select all columns except date and check their classes
  x1 <- purrr::map_chr(data %>% dplyr::select(-date), class)
  # Identify non-numeric variable types
  non_numeric_types <- base::setdiff(x1, c("numeric", "integer"))
  if (length(non_numeric_types) > 0) {
    warning("Data contains non-numeric variables (",
            paste(non_numeric_types, collapse = ", "),
            "). Please check your data! For example, check if any cell contains extra spaces?")
  }

  # Check 4: Check for missing values throughout the entire dataset
  missing_count <- base::length(base::which(is.na(data)))
  if (missing_count > 0) {
    warning("Data contains ", missing_count, " missing values. Please check carefully!")
  }

  # All checks passed successfully
  message("Data validation passed!")
  return(invisible(TRUE))
}
