#' Cross-Quantile Test for ECM Model Coefficients
#'
#' This function performs a Wald test to compare coefficients across different
#' quantiles in the Error Correction Model (ECM) framework. It can test for
#' equality of coefficients across specified quantile levels, either individually
#' or jointly for multiple variables. The test is particularly useful for
#' examining parameter stability across the conditional distribution.
#'
#' @param estimation List. Estimation results obtained from `coef_estimate()`.
#' @param order Integer vector. ARDL orders obtained from `order_find()`.
#'   Length must match the number of variables, otherwise an error is raised.
#' @param boots Integer. Number of bootstrap replications used in estimation.
#'   This should match the value used in `coef_estimate()`.
#' @param taus Numeric vector. All quantile levels estimated, same as the
#'   `taus` parameter in `coef_estimate()`.
#' @param test_taus Numeric vector. Specific quantile levels to test for equality.
#'   Must be a subset of `taus` and contain at least two distinct quantiles.
#' @param variable Character vector. Names of variables to test. Can include
#'   "ECM", "Intercept", or regular variable names. Special cases:
#'   - If testing the dependent variable with term = "long", tests the ECM coefficient.
#'   - If testing a variable with ARDL order 0, tests the contemporaneous effect.
#'   - Cannot mix ECM/Intercept with regular variables in the same test.
#'   - If length of variable larger than 1. For example, variable = c("x1", "x2"),
#'   then this function will test H0: beta_1(tau1) = beta_1(tau2), beta_2(tau1) = beta_2(tau2)
#' @param term Character. For regular variables, specify "long" for long-term
#'   coefficients or "short" for short-term coefficients. Only applicable when
#'   testing regular variables (not ECM or Intercept).
#' @param joint Logical. For short-term coefficients, if TRUE, tests the sum
#'   of coefficients across lags. If FALSE, tests each lag separately.
#'   Default is TRUE.
#'
#' @return An object of class "htest" containing Wald test results, including:
#'   - statistic: The Wald test statistic
#'   - p.value: The p-value for the test
#'   - parameter: Degrees of freedom
#'   - method: Name of the test
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
#' # Estimate models at multiple quantiles
#' est_result <- coef_estimate(
#'   data = sample_data,
#'   order = optimal_order,
#'   taus = c(0.25, 0.5, 0.75),
#'   boots = 50,
#'   export = FALSE
#' )
#'
#' # Test if ECM coefficient is equal across quantiles 0.25, 0.5, and 0.75
#' test_result <- cross_quantile_test(
#'   estimation = est_result,
#'   order = optimal_order,
#'   boots = 50,
#'   taus = c(0.25, 0.5, 0.75),
#'   test_taus = c(0.25, 0.5, 0.75),
#'   variable = "ECM",
#'   joint = TRUE
#' )
#'
#' print(test_result)
#' }
cross_quantile_test <- function(estimation, order, boots, taus, test_taus, variable, term = NULL, joint = TRUE){
  # Validate that at least two distinct quantiles are specified for testing
  if (base::length(base::unique(test_taus)) < 2){
    stop("The test_taus parameter should contain at least two distinct quantile levels!")
  }

  # Validate variable names
  valid_vars <- base::c("ECM", "Intercept", base::names(order))
  if (base::length(base::setdiff(variable, valid_vars)) > 0){
    stop("Please check that all variable names in the 'variable' parameter are correct!")
  }

  # Validate that all test_taus are in the original taus
  if (!base::all(test_taus %in% taus)){
    stop("Some quantile levels in test_taus are not present in the original taus parameter!")
  }

  # Extract models for the specified quantiles
  models <- base::vector("list", base::length(test_taus))
  for (i in 1:base::length(test_taus)){
    find <- base::which(taus == test_taus[i])
    models[[i]] <- estimation$QR[[find]]
  }

  # Extract ECM models (coefficients and bootstrap samples)
  ECM_models <- base::vector("list", base::length(test_taus))
  for (i in 1:base::length(test_taus)){
    ECM_models[[i]] <- base::list(
      boot = models[[i]]$QECM_boot,
      coef = models[[i]]$QECM_coef
    )
  }

  # Get coefficient names and rename the second to "ECM" (error correction term)
  coef_names <- base::names(ECM_models[[1]]$coef)
  coef_names[2] <- "ECM"

  # Create matrices to store coefficients and bootstrap coefficients
  coefs <- base::matrix(NA, base::length(coef_names), base::length(test_taus))
  base::rownames(coefs) <- coef_names
  base::colnames(coefs) <- stringr::str_c("Tau", test_taus)

  # Create 3D array for bootstrap coefficients: variables x quantiles x bootstrap samples
  boot_coefs <- base::array(
    NA,
    dim = base::c(base::length(coef_names), base::length(test_taus), boots)
  )
  base::dimnames(boot_coefs)[[1]] <- coef_names
  base::dimnames(boot_coefs)[[2]] <- stringr::str_c("Tau", test_taus)

  # Fill the coefficient matrices
  for (i in 1:base::length(test_taus)){
    coefs[, i] <- ECM_models[[i]]$coef
    boot_coefs[, i, ] <- base::t(base::as.matrix(ECM_models[[i]]$boot))
  }

  # Case 1: Variable with ARDL order 0 (contemporaneous effect only)
  if (0 %in% order[variable] & base::length(variable) > 1){
    stop("When testing multiple variables, variables with ARDL order 0 are not allowed!")
  } else if (0 %in% order[variable] & base::length(variable) == 1){
    base::message("Variable ", variable, " has ARDL order 0. Long-term vs. short-term distinction is not applicable.")

    # Extract coefficients for the zero-order variable
    find <- base::which(coef_names == stringr::str_c("L0.", variable))
    b <- coefs[find, ]
    boot_b <- base::t(boot_coefs[find, , ])

    # Calculate covariance matrix
    cov <- stats::cov(boot_b)

    # Create contrast matrix for testing equality across quantiles
    R <- base::matrix(0, nrow = base::length(b) - 1, ncol = base::length(b))
    R[, 1] <- 1
    for (i in 1:base::nrow(R)){
      R[i, i + 1] <- -1
    }

    # Perform Wald test
    res <- wald_test(b = b, cov = cov, R = R)
  }
  # Case 2: Single variable tests
  else {
    # 2a: Intercept test
    if (base::length(variable) == 1 & "Intercept" %in% variable){
      b <- coefs[1, ]
      boot_b <- base::t(boot_coefs[1, , ])
      cov <- stats::cov(boot_b)

      R <- base::matrix(0, nrow = base::length(b) - 1, ncol = base::length(b))
      R[, 1] <- 1
      for (i in 1:base::nrow(R)){
        R[i, i + 1] <- -1
      }

      res <- wald_test(b = b, cov = cov, R = R)
    }
    # 2b: ECM term test
    else if (base::length(variable) == 1 & "ECM" %in% variable){
      b <- coefs[2, ]
      boot_b <- base::t(boot_coefs[2, , ])
      cov <- stats::cov(boot_b)

      R <- base::matrix(0, nrow = base::length(b) - 1, ncol = base::length(b))
      R[, 1] <- 1
      for (i in 1:base::nrow(R)){
        R[i, i + 1] <- -1
      }

      res <- wald_test(b = b, cov = cov, R = R)
    }
    # 2c: Dependent variable test (defaults to ECM coefficient)
    else if (base::length(variable) == 1 & base::colnames(estimation$OLS$ARDL_estimate$model)[1] %in% variable & term == "long"){
      base::message("The variable is the dependent variable. Testing the ECM coefficient by default.")
      b <- coefs[2, ]
      boot_b <- base::t(boot_coefs[2, , ])
      cov <- stats::cov(boot_b)

      R <- base::matrix(0, nrow = base::length(b) - 1, ncol = base::length(b))
      R[, 1] <- 1
      for (i in 1:base::nrow(R)){
        R[i, i + 1] <- -1
      }

      res <- wald_test(b = b, cov = cov, R = R)
    }
    # 2d: Invalid combinations
    else if (base::colnames(estimation$OLS$ARDL_estimate$model)[1] %in% variable & base::length(variable) > 1){
      stop("Testing equality between ECM coefficient and other variable coefficients is not allowed!")
    } else if ("ECM" %in% variable & base::length(variable) > 1){
      stop("Testing equality between ECM coefficient and other variable coefficients is not allowed!")
    }
    # 2e: Long-term coefficient test
    else if (term == "long"){
      # Find coefficients for the specified variables
      find1 <- base::which(stringr::str_detect(coef_names, stringr::str_c(variable, collapse = "|")))
      find2 <- base::which(base::substring(coef_names, 1, 2) == "L1" | coef_names == "ECM")
      find <- base::intersect(find1, find2)

      b <- base::vector("list", base::length(find))
      boot_b <- base::vector("list", base::length(find))
      R <- base::vector("list", base::length(find))

      for (i in 1:base::length(find)){
        b[[i]] <- coefs[find[i], ]
        boot_b[[i]] <- base::t(boot_coefs[find[i], , ])
        R[[i]] <- base::matrix(0, nrow = base::length(b[[i]]) - 1, ncol = base::length(b[[i]]))
        R[[i]][, 1] <- 1
        for (j in 1:base::nrow(R[[i]])){
          R[[i]][j, j + 1] <- -1
        }
      }

      b <- base::unlist(b)
      cov <- stats::cov(purrr::map_dfc(boot_b, tibble::as_tibble))
      R <- Matrix::bdiag(R)
      res <- wald_test(b = b, cov = cov, R = as.matrix(R))
    }
    # 2f: Short-term coefficient test
    else if (term == "short"){
      # Find short-term coefficients for the specified variables
      find1 <- base::which(stringr::str_detect(coef_names, stringr::str_c(variable, collapse = "|")))
      find2 <- base::which(base::substring(coef_names, 1, 1) == "D")
      find <- base::intersect(find1, find2)

      if (base::length(find) == 0){
        stop("None of the specified variables have short-term coefficients to test!")
      } else {
        # Check which variables actually have short-term coefficients
        xx <- base::unique(base::substring(coef_names[find], 4, base::nchar(coef_names[find])))
        if (base::length(variable) > base::length(xx)){
          base::warning(
            "Only ", stringr::str_c(xx, collapse = " and "), " have short-term coefficients. ",
            "Variable(s) ", stringr::str_c(base::setdiff(variable, xx), collapse = " and "),
            " in the 'variable' parameter are not meaningful!"
          )
        }
      }

      # Non-joint test: test each coefficient separately
      if (!joint){
        b <- base::vector("list", base::length(find))
        boot_b <- base::vector("list", base::length(find))
        R <- base::vector("list", base::length(find))

        for (i in 1:base::length(find)){
          b[[i]] <- coefs[find[i], ]
          boot_b[[i]] <- base::t(boot_coefs[find[i], , ])
          R[[i]] <- base::matrix(0, nrow = base::length(b[[i]]) - 1, ncol = base::length(b[[i]]))
          R[[i]][, 1] <- 1
          for (j in 1:base::nrow(R[[i]])){
            R[[i]][j, j + 1] <- -1
          }
        }

        b <- base::unlist(b)
        cov <- stats::cov(purrr::map_dfc(boot_b, tibble::as_tibble))
        R <- Matrix::bdiag(R)
        res <- wald_test(b = b, cov = cov, R = as.matrix(R))
      }

      # Joint test: test the sum of coefficients across lags
      if (joint){
        xx <- base::substring(coef_names[find], 4, base::nchar(coef_names[find]))

        # Calculate summed coefficients across lags
        if (base::length(xx) == 1){
          b <- base::as.vector(coefs[find, ])
        } else {
          b <- coefs[find, ] %>%
            tibble::as_tibble() %>%
            dplyr::mutate(name = xx) %>%
            dplyr::summarise(dplyr::across(dplyr::where(is.numeric), base::sum), .by = name) %>%
            dplyr::select(-name) %>%
            base::as.matrix() %>%
            base::t() %>%
            base::as.vector()
        }

        # Calculate covariance matrix for summed coefficients
        if (base::length(xx) == 1){
          cov <- stats::cov(base::t(boot_coefs[find, , ]))
        } else {
          # Use reshape2 to melt the 3D array, then sum across lags for each bootstrap sample
          cov <- boot_coefs[find, , ] %>%
            reshape2::melt() %>%
            dplyr::mutate(Var1 = base::as.character(Var1)) %>%
            dplyr::mutate(Var1 = base::substring(Var1, 4, base::nchar(Var1))) %>%
            dplyr::summarise(value = base::sum(value), .by = c(Var1, Var2, Var3)) %>%
            dplyr::group_split(Var1) %>%
            purrr::map(~dplyr::select(., -Var1)) %>%
            purrr::map(~tidyr::pivot_wider(., id_cols = Var3, names_from = Var2, values_from = value)) %>%
            purrr::map(~dplyr::select(., -Var3)) %>%
            dplyr::bind_cols() %>%
            stats::cov()
        }

        # Create block-diagonal contrast matrix
        R <- base::vector("list", base::length(base::unique(xx)))
        for (i in 1:base::length(R)){
          R[[i]] <- base::matrix(0, nrow = base::length(test_taus) - 1, ncol = base::length(test_taus))
          R[[i]][, 1] <- 1
          for (j in 1:base::nrow(R[[i]])){
            R[[i]][j, j + 1] <- -1
          }
        }

        R <- Matrix::bdiag(R)
        res <- wald_test(b = b, cov = cov, R = as.matrix(R))
      }
    }
  }

  return(res)
}
