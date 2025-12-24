#' Plot Heatmap of Dynamic Multipliers Across Quantiles and Horizons
#'
#' This function creates a heatmap visualization of dynamic multipliers for a
#' specific variable across different quantile levels and time horizons. It
#' provides an intuitive way to analyze how the effect of a variable evolves
#' over time and across the conditional distribution.
#'
#' @param multipliers A data frame containing dynamic multiplier results,
#'   typically obtained from `all_dynamic_multipliers()`. Must contain columns:
#'   'variable', 'tau', 'step', and 'value'.
#' @param plot_variable Character. The name of the variable to plot. Must be one
#'   of the explanatory variables present in the multipliers data frame.
#' @param export Logical. Whether to export the plot to a file. Default is TRUE.
#' @param file Character. File path for saving the plot. Required if export = TRUE.
#'
#' @return A ggplot object representing the heatmap visualization. The plot
#'   shows time horizons on the x-axis, quantile levels on the y-axis, and
#'   multiplier values color-coded using a turbo color scale.
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
#' # Calculate dynamic multipliers for all variables
#' multiplier_results <- all_dynamic_multipliers(
#'   data = sample_data,
#'   order = optimal_order,
#'   taus = c(0.25, 0.5, 0.75),
#'   steps = 10,
#'   cumulative = FALSE,
#'   impulse = FALSE,
#'   export = FALSE
#' )
#'
#' # Create heatmap for a specific variable
#' heatmap_plot <- plot_heatmap_multipliers(
#'   multipliers = multiplier_results,
#'   plot_variable = "x1",  # Replace with actual variable name from your data
#'   export = FALSE
#' )
#'
#' # Display the heatmap
#' print(heatmap_plot)
#' }
plot_heatmap_multipliers <- function(multipliers, plot_variable, export = TRUE, file = NULL){

  # Get unique variable names from the multipliers data
  xx <- base::unique(multipliers$variable)

  # Validate that the specified plot_variable exists in the data
  if (!(plot_variable %in% xx)){
    stop("The plot_variable parameter must be one of: ",
         stringr::str_c(xx, collapse = ", "), "!")
  }

  # Create heatmap visualization using tidyplots framework
  multipliers %>%
    dplyr::filter(tau != "ols") %>%  # Exclude OLS results if present
    dplyr::mutate(tau = readr::parse_number(tau)) %>%  # Convert tau to numeric
    dplyr::filter(variable == plot_variable) -> multipliers1   # Filter for specified variable
  multipliers1 %>%
    tidyplots::tidyplot(step, tau, color = value) %>%  # Create base plot
    tidyplots::add_heatmap() %>%  # Add heatmap layer
    tidyplots::adjust_colors(new_colors = tidyplots::colors_continuous_turbo) %>%  # Set color scale
    tidyplots::remove_legend_title() %>%  # Clean up legend
    tidyplots::adjust_x_axis_title("Horizon") %>%  # Set x-axis label
    tidyplots::adjust_y_axis_title("Quantile") %>%  # Set y-axis label
    tidyplots::adjust_x_axis(breaks = base::seq(0, base::max(multipliers1$step), 2)) %>%  # X-axis breaks
    tidyplots::adjust_y_axis(breaks = base::seq(base::min(multipliers1$tau),
                                                base::max(multipliers1$tau), 0.1)) -> p  # Y-axis breaks

  # Export plot if requested
  if (export){
    if (base::is.null(file)) {
      stop("File path must be specified when export = TRUE")
    }
    tidyplots::save_plot(p, file)
  }

  # Display the plot
  base::print(p)

  # Return the plot object for further customization
  return(p)
}
