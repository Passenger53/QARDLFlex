#' Plot Line Graph of Dynamic Multipliers Across Quantiles
#'
#' This function creates a line plot visualization of dynamic multipliers for a
#' specific variable across different quantile levels and time horizons. It
#' provides an intuitive way to analyze how the effect of a variable evolves
#' over time at different points in the conditional distribution.
#'
#' @param multipliers A data frame containing dynamic multiplier results,
#'   typically obtained from `all_dynamic_multipliers()`. Must contain columns:
#'   'variable', 'tau', 'step', and 'value'.
#' @param plot_variable Character. The name of the variable to plot. Must be one
#'   of the explanatory variables present in the multipliers data frame.
#'   Currently, only single variable plotting is supported.
#' @param plot_taus Numeric vector. The quantile levels to include in the plot.
#'   Must be a subset of the quantiles estimated in the original model.
#' @param export Logical. Whether to export the plot to a file. Default is TRUE.
#' @param file Character. File path for saving the plot. Required if export = TRUE.
#'
#' @return A ggplot object representing the line plot visualization. The plot
#'   shows time horizons on the x-axis, multiplier values on the y-axis, with
#'   different quantile levels distinguished by color and linetype. When only
#'   one quantile is plotted, the legend is automatically removed for cleaner
#'   visualization.
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
#' # Create line plot for a specific variable with multiple quantiles
#' line_plot <- plot_line_multipliers(
#'   multipliers = multiplier_results,
#'   plot_variable = "x1",  # Replace with actual variable name from your data
#'   plot_taus = c(0.25, 0.5, 0.75),
#'   export = FALSE
#' )
#'
#' # Display the plot
#' print(line_plot)
#'
#' # Example with single quantile (legend automatically removed)
#' single_quantile_plot <- plot_line_multipliers(
#'   multipliers = multiplier_results,
#'   plot_variable = "x1",
#'   plot_taus = 0.5,
#'   export = FALSE
#' )
#' print(single_quantile_plot)
#' }
plot_line_multipliers <- function(multipliers, plot_variable, plot_taus, export = TRUE, file = NULL){

  # Get unique variable names from the multipliers data
  xx <- base::unique(multipliers$variable)

  # Validate that the specified plot_variable exists in the data
  if (!(plot_variable %in% xx)){
    stop("The plot_variable parameter must be one of: ",
         stringr::str_c(xx, collapse = ", "), "!")
  }

  # Extract and process available quantile levels from the multipliers data
  taus <- multipliers %>%
    dplyr::filter(tau != "ols") %>%  # Exclude OLS results
    dplyr::pull(tau) %>%  # Extract tau column
    base::unique() %>%  # Get unique values
    readr::parse_number() %>%  # Convert to numeric
    stats::na.omit()  # Remove any NA values

  # Validate that all specified plot_taus are available in the data
  if (!base::all(plot_taus %in% taus)){
    stop("All quantile levels in plot_taus must be among those estimated in the original model!")
  }

  # Create line plot visualization using tidyplots framework
  multipliers %>%
    dplyr::filter(tau != "ols") %>%  # Exclude OLS results
    dplyr::filter(variable == plot_variable) %>%  # Filter for specified variable
    dplyr::mutate(tau = readr::parse_number(tau)) %>%  # Convert tau to numeric
    dplyr::filter(tau %in% plot_taus) %>%  # Filter for specified quantiles
    dplyr::mutate(tau = base::as.character(tau),  # Convert back to character for labeling
                  tau = stringr::str_c("Tau = ", tau),  # Create descriptive labels
                  tau = base::factor(tau)) -> multipliers1 # Convert to factor for proper legend ordering
  multipliers1 %>%
    tidyplots::tidyplot(step, value, color = tau, linetype = tau) %>%  # Create base plot
    tidyplots::add_line(linewidth = 0.7) %>%  # Add lines with specified thickness
    tidyplots::add_reference_lines(y = 0,  # Add reference line at y=0
                                   linewidth = 0.5) %>%  # Set reference line thickness
    tidyplots::remove_legend_title() %>%  # Remove legend title for cleaner appearance
    tidyplots::adjust_x_axis_title("Horizon") %>%  # Set x-axis label
    tidyplots::adjust_y_axis_title("Dynamic Multiplier") %>%  # Set y-axis label
    tidyplots::adjust_x_axis(breaks = base::seq(0, base::max(multipliers1$step), 2)) -> p  # Set x-axis breaks

  # Automatically remove legend when only one quantile is plotted
  if (base::length(plot_taus) == 1){
    p %>%
      tidyplots::remove_legend() -> p  # Remove legend for cleaner single-line plots
  }

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
