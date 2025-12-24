#' Wald Test for Linear Hypotheses
#'
#' This function performs a Wald test for linear hypotheses on coefficient vectors.
#' It tests the null hypothesis that the linear combinations of coefficients
#' specified by the constraint matrix R are equal to the values in r.
#'
#' @param b Numeric vector. Coefficient vector of length p.
#' @param cov Numeric matrix. Variance-covariance matrix of dimensions p x p.
#' @param R Numeric matrix. Constraint matrix of dimensions q x p, where q is
#'        the number of constraints.
#' @param r Numeric vector. Constraint values of length q. Default is zero vector
#'        (testing if linear combinations equal zero).
#'
#' @return An object of class "htest" containing:
#'   - statistic: The Wald test statistic
#'   - p.value: The p-value for the test
#'   - parameter: Degrees of freedom (number of constraints)
#'   - method: Name of the test ("Wald Test")
#'
#' @export
#'
#' @examples
#' # Test whether the second and third coefficients are simultaneously zero
#' b <- c(1.5, -0.8, 0.3)  # Coefficient vector
#' cov <- matrix(c(0.1, 0.02, 0.01,
#'                0.02, 0.05, 0.03,
#'                0.01, 0.03, 0.08), nrow = 3)  # Variance-covariance matrix
#' R <- matrix(c(0, 1, 0,
#'              0, 0, 1), nrow = 2, byrow = TRUE)  # Constraint matrix
#' result <- wald_test(b, cov, R)
#' print(result)
#'
#' # Example with ARDL model coefficients
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
#' # Estimate ARDL model
#' ardl_result <- ardl_estimate(sample_data, order = optimal_order)
#'
#' # Test if two coefficients are simultaneously zero
#' b_coef <- coef(ardl_result$ARDL_estimate)
#' cov_matrix <- vcov(ardl_result$ARDL_estimate)
#' R_matrix <- matrix(c(0, 1, 0, 0,
#'                      0, 0, 1, 0), nrow = 2, byrow = TRUE)
#' wald_result <- wald_test(b_coef, cov_matrix, R_matrix)
#' print(wald_result)
#' }
wald_test <- function(b, cov, R, r = NULL) {
  # Input validation and preprocessing
  if (base::is.null(base::dim(b))) {
    b <- base::as.matrix(b)  # Convert b to column vector if it's not a matrix
  }
  p <- base::length(b)  # Number of parameters

  # Validate covariance matrix dimensions
  if (base::nrow(cov) != p || base::ncol(cov) != p) {
    stop("cov must be a ", p, "x", p, " matrix")
  }

  # Validate constraint matrix dimensions
  if (base::ncol(R) != p) {
    stop("R must have ", p, " columns to match the length of b")
  }

  q <- base::nrow(R)  # Number of constraints
  if (q == 0) {
    stop("R must have at least one row (at least one constraint is required for hypothesis testing)")
  }

  # Set default constraint values to zero if not provided
  if (base::is.null(r)) {
    r <- base::rep(0, q)
  } else if (base::length(r) != q) {
    stop("r must have length ", q, " to match the number of rows in R")
  }

  # Calculate the difference vector: R*b - r
  diff_vec <- R %*% b - base::as.matrix(r)

  # Calculate the middle matrix: M = R %*% cov %*% t(R)
  M <- R %*% cov %*% base::t(R)

  # Calculate Wald statistic: W = t(diff_vec) %*% solve(M) %*% diff_vec
  # The statistic follows a chi-square distribution with q degrees of freedom
  W <- base::t(diff_vec) %*% base::solve(M) %*% diff_vec
  W_stat <- base::as.numeric(W)  # Convert to scalar

  # Calculate p-value using chi-square distribution
  p_value <- stats::pchisq(W_stat, df = q, lower.tail = FALSE)

  # Create result object with htest class for compatibility with R's hypothesis testing framework
  result <- base::list(
    statistic = W_stat,
    p.value = p_value,
    parameter = q,  # Degrees of freedom
    method = "Wald Test"
  )

  base::class(result) <- "htest"  # Set as hypothesis test result class

  return(result)
}
