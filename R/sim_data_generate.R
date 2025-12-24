#' Generate Simple Simulated Data for ARDL Models
#'
#' This function creates bootstrap datasets for quantile ARDL model simulation studies.
#' It generates multiple time series datasets with specified ARDL structure, coefficients,
#' and properties for Monte Carlo simulations.
#'
#' @param b Numeric vector. Coefficients for the ARDL model including intercept and lag terms.
#'        The length of this vector must match the total number of parameters implied by the
#'        `order` argument. Specifically, the length should be: sum(order + 1).
#'        For example, an ARDL(1,2) model requires ((1+1) + (2+1)) = 5 coefficients:
#'        intercept, y_lag1, x_contemporaneous, x_lag1, x_lag2.
#' @param order Integer vector. Lag orders for each variable in the system. The first element
#'        should be the lag order for the dependent variable, followed by orders for independent variables.
#' @param size Integer. Length of each time series to generate.
#' @param reps Integer. Number of bootstrap replications/datasets to generate.
#'
#' @details
#' A random seed `set.seed(123)`for generating residuals and a random seed
#'     `set.seed(234)`for generating mutually independent, multivariate normally
#'     distributed explanatory variables have been built-in.
#'
#' @return A list of length `reps` containing simulated datasets. Each dataset is a tibble
#'         with columns: date (year-month format), y (dependent variable), x1, x2, etc.
#'         (independent variables).
#'
#' @export
#'
#' @examples
#' # Generate 2 datasets with 100 observations each
#' # ARDL(1,1) model: y_t = 0.5 + 0.7*y_{t-1} + 0.3*x_{t} + 0.2*x_{t-1} + error
#' # Coefficient vector length: 1 (intercept) + 1 (y_lag1) + 2 (x_contemporaneous + x_lag1) = 4
#' sim_data <- sim_data_generate(b = c(0.5, 0.7, 0.3, 0.2),
#'                              order = c(1, 1),
#'                              size = 100,
#'                              reps = 2)
#' str(sim_data)
#'
#' # ARDL(1,2) model example: requires 5 coefficients
#' # 1 (intercept) + 1 (y_lag1) + 3 (x_contemporaneous + x_lag1 + x_lag2) = 5
#' sim_data2 <- sim_data_generate(b = c(0.5, 0.6, 0.3, 0.2, 0.1),
#'                               order = c(1, 2),
#'                               size = 50,
#'                               reps = 1)
sim_data_generate <- function(b, order, size, reps){

  # Validate coefficient vector length matches ARDL order requirements
  # Total parameters = intercept (1) + sum of (lag order + 1) for each variable
  required_length <- sum(order + 1)
  if (length(b) != required_length) {
    stop("Coefficient vector length mismatch. For order = c(",
         paste(order, collapse = ", "), "), required ", required_length,
         " coefficients, but got ", length(b), ".")
  }

  # Set seed for reproducible bootstrap errors
  set.seed(123)
  boot_e <- stats::rnorm(size * reps, 0, 1)
  boot_e <- matrix(boot_e, nrow = size)

  # Create covariance matrix for multivariate normal distribution
  # Identity matrix for independent variables
  sigma <- matrix(0, length(order), length(order))
  diag(sigma) <- 1

  # Generate initial multivariate normal data
  # Creates base dataset with specified number of variables
  set.seed(234)
  dd <- MASS::mvrnorm(size, rep(0, length(order)), sigma)
  colnames(dd) <- c("y", stringr::str_c("x", 1:(length(order) - 1)))

  # Initialize list to store bootstrap datasets
  boot_data <- vector("list", reps)
  for (i in 1:reps){
    boot_data[[i]] <- dd
  }

  rows <- size
  max_lag <- max(order)

  # Create location matrix to track variable lags
  # Matrix stores which variable (col) and which lag (row) each parameter corresponds to
  location <- NULL
  for (i in 1:length(order)){
    if (i == 1){
      # For dependent variable (y), include lags 1 through order[1]
      # Excludes contemporaneous term for dependent variable
      xx <- cbind(i, 1:order[i])
      colnames(xx) <- c("col", "row")
      location <- rbind(location, xx)
    } else {
      # For independent variables, include contemporaneous (lag 0) and lags
      xx <- cbind(i, 0:order[i])
      colnames(xx) <- c("col", "row")
      location <- rbind(location, xx)
    }
  }

  # Generate ARDL process for each bootstrap replication
  for (i in 1:reps){
    ee <- boot_e[,i]  # Error term for this replication

    # Generate ARDL process starting after max lag period
    # Initial periods are used as starting values
    for (j in (max(order) + 1):rows){
      values <- rep(NA, nrow(location))

      # Collect lagged values for ARDL equation
      for (k in 1:nrow(location)){
        rr <- location[k, 2]  # Lag order
        cc <- location[k, 1]  # Variable column
        values[k] <- boot_data[[i]][j - rr, cc]
      }

      # Add intercept and calculate ARDL value: y_t = b0 + b1*y_{t-1} + ... + error
      # Coefficient vector b: first element is intercept, rest match location matrix order
      values <- c(1, values)  # Include intercept (first element of b)
      boot_data[[i]][j, 1] <- sum(b * values) + ee[j]
    }

    # Convert to tibble and add date column for better data handling
    boot_data[[i]] <- tibble::as_tibble(boot_data[[i]])

    # Create monthly date sequence for time series data
    # Generates monthly dates starting from calculated year
    start <- 2025 - floor(size / 12)
    range <- start:2025
    stringr::str_c(rep(range, each = 12),
                   "-",
                   rep(1:12, length(range))) %>%
      lubridate::ym() -> dates
    dates <- dates[1:size]

    # Add date column and rearrange columns with date first
    boot_data[[i]] %>%
      dplyr::mutate(date = dates) %>%
      dplyr::relocate(date, 1) -> boot_data[[i]]
  }

  return(boot_data)
}
