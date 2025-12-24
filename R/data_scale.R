#' Scale/Normalize Data for Analysis
#'
#' This function standardizes (z-score normalization) the numeric variables in a data frame.
#' It preserves the date column if present, applying scaling only to numeric variables.
#'
#' @param data A data frame containing the data to be scaled. If the first column is named "date",
#'        it will be preserved and not scaled.
#' @param scale Logical. If TRUE, scaling is applied to all numeric columns. If FALSE,
#'        the original data is returned unchanged.
#'
#' @return A data frame with the same structure as the input, but with numeric columns
#'         scaled to have mean 0 and standard deviation 1. The date column (if present)
#'         remains unchanged.
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
#' # Apply scaling
#' scaled_data <- data_scale(sample_data, scale = TRUE)
#'
#' # Check the results
#' head(scaled_data)
#' summary(scaled_data)
data_scale <- function(data, scale){
  # Only proceed if scaling is requested
  if (scale){
    # Check if first column is "date" - if yes, preserve it during scaling
    if (colnames(data)[1] == "date"){
      # Select all columns except date for scaling
      dd <- dplyr::select(data, -date)
      # Apply scaling to all columns: z-score normalization (mean=0, sd=1)
      dd %>% dplyr::mutate_all(~as.vector(base::scale(.))) -> dd
      # Recombine date column with scaled data
      data %>% dplyr::select(date) %>% dplyr::bind_cols(dd) -> data
    } else {
      # If no date column, scale all columns directly
      data %>% dplyr::mutate_all(~as.vector(base::scale(.))) -> data
    }
  }
  # Return the data (scaled if scale=TRUE, original if scale=FALSE)
  return(data)
}
