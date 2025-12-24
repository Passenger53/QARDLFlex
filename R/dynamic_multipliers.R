#' Dynamic Multipliers Calculation for ARDL Models
#'
#' This function computes dynamic multipliers that trace the time path of the
#' effect of an exogenous variable on the dependent variable in an ARDL model
#' framework. It can handle both permanent (step) and temporary (impulse)
#' shocks, and provides options for cumulative multipliers.
#'
#' @param phi Numeric vector. Autoregressive coefficients in the ARDL
#'   representation, ordered from lag 1 to lag p (e.g., c(phi1, phi2, ..., phip)).
#'   The length determines the autoregressive order p.
#' @param beta Numeric vector. Distributed lag coefficients for the exogenous
#'   variable, ordered from contemporaneous effect (lag 0) to maximum lag q
#'   (e.g., c(beta0, beta1, ..., betaq)). The length determines the lag order q.
#' @param steps Integer. Number of time periods after the shock to calculate
#'   multipliers for. Must be a positive integer greater than 0.
#' @param cumulative Logical. If TRUE, returns cumulative multipliers (the sum
#'   of multipliers up to each period). If FALSE (default), returns period-by-period
#'   multipliers.
#' @param impulse Logical. Specifies the type of shock:
#'   - FALSE (default): Permanent step shock (x_t = 1 for all t ≥ 0)
#'   - TRUE: Temporary impulse shock (x_0 = 1, x_t = 0 for t > 0)
#'
#' @return A numeric vector of dynamic multipliers with names indicating the
#'   time period (h = 0, 1, ..., steps-1). The vector contains either
#'   period-specific multipliers or cumulative multipliers depending on the
#'   'cumulative' parameter.
#'
#' @export
#'
#' @examples
#' \dontrun{
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
#' # Example 1: Permanent step shock with ARDL(1,1) model
#' phi <- c(0.6)    # AR coefficient
#' beta <- c(0.8, 0.3)  # DL coefficients (beta0, beta1)
#' multipliers <- dynamic_multipliers(
#'   phi = phi,
#'   beta = beta,
#'   steps = 10,
#'   cumulative = FALSE,
#'   impulse = FALSE
#' )
#' print(multipliers)
#'
#' # Example 2: Temporary impulse shock with cumulative multipliers
#' cumulative_multipliers <- dynamic_multipliers(
#'   phi = phi,
#'   beta = beta,
#'   steps = 10,
#'   cumulative = TRUE,
#'   impulse = TRUE
#' )
#' print(cumulative_multipliers)
#'
#' # Example 3: No autoregressive component (distributed lag model only)
#' multipliers_dl <- dynamic_multipliers(
#'   phi = numeric(0),  # Empty vector for no AR component
#'   beta = c(0.5, 0.2, 0.1),
#'   steps = 8,
#'   cumulative = FALSE,
#'   impulse = FALSE
#' )
#' print(multipliers_dl)
#' }
dynamic_multipliers <- function(phi, beta, steps, cumulative = FALSE, impulse = FALSE) {

  # Input validation: steps must be positive integer
  if (steps <= 0) {
    stop("Number of steps must be a positive integer")
  }

  # Check stationarity condition for autoregressive component
  if (base::length(phi) > 0) {
    # Create characteristic polynomial: 1 - phi1*z - phi2*z^2 - ... - phi_p*z^p
    poly_coef <- base::c(1, -phi)
    # Find roots of characteristic polynomial
    roots <- base::polyroot(poly_coef)
    # Check if any root lies inside or on unit circle (non-stationarity)
    if (base::any(base::abs(roots) <= 1)) {
      base::warning("Autoregressive part may be non-stationary")
    }
  }

  # Determine model orders
  p <- base::length(phi)  # Autoregressive order
  q <- base::length(beta) - 1  # Distributed lag order (q = length(beta) - 1)
  sum_beta <- base::sum(beta)  # Sum of all distributed lag coefficients

  # Initialize multiplier vector
  multipliers <- base::numeric(steps)

  # Calculate dynamic multipliers for each period
  for (h in 1:steps) {
    t <- h - 1  # Current time period (t=0 is shock occurrence period)

    # 1. Calculate autoregressive component effect: sum(phi_j * multiplier_{t-j})
    ar_effect <- 0
    # Only calculate if AR order > 0 and we have past periods (t > 0)
    if (p > 0 && t > 0) {
      # Sum over available lags (min(p, t) ensures we don't go beyond available data)
      for (j in base::seq_len(base::min(p, t))) {
        ar_effect <- ar_effect + phi[j] * multipliers[h - j]
      }
    }

    # 2. Calculate distributed lag component effect
    if (impulse) {
      # Temporary impulse shock: x_0 = 1, x_t = 0 for t > 0
      if (t <= q && t >= 0) {
        # For impulse shock, effect is just the contemporaneous beta_t
        # beta[1] corresponds to lag 0, beta[2] to lag 1, etc.
        dl_effect <- if (t < base::length(beta)) beta[t + 1] else 0
      } else {
        dl_effect <- 0
      }
    } else {
      # Permanent step shock: x_t = 1 for all t >= 0
      if (t <= q && t >= 0) {
        # For step shock, effect is sum of coefficients up to current lag
        dl_effect <- base::sum(beta[base::seq_len(t + 1)])
      } else {
        # Beyond maximum lag, effect is sum of all coefficients
        dl_effect <- sum_beta
      }
    }

    # Total multiplier = AR effect + DL effect
    multipliers[h] <- ar_effect + dl_effect
  }

  # Calculate cumulative multipliers if requested
  if (cumulative) {
    multipliers <- base::cumsum(multipliers)
  }

  # Add period labels for clarity
  base::names(multipliers) <- base::paste("h =", 0:(steps-1))

  return(multipliers)
}
