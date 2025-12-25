
<img width="5333" height="2666" alt="QARDLFlex_Social_Preview" src="https://github.com/user-attachments/assets/e095f17a-ea4b-4662-aaf3-ac029a5c232a" />


## Introduction

### What is `QARDLFlex`?

A **Flex**ible implementation for the Quantile Autoregressive Distributed Lag (**QARDL**) model has been lacking. To address this gap, I developed `QARDLFlex`, an R package designed to provide a flexible and reliable toolkit for QARDL-related analysis. Although `QARDLFlex` is currently in its early development stages, it already offers a complete workflow for empirical research—from descriptive statistics and QARDL model estimation to parameter hypothesis testing.

The package is available on GitHub at <https://github.com/Passenger53/QARDLFlex>. Installation and usage details will be provided later. For users in mainland China who may experience intermittent access to GitHub, persistent retries are recommended.

### What can `QARDLFlex` Do?

Currently, `QARDLFlex` can perform the following tasks:

-   **Descriptive Statistics**: Calculate basic descriptive statistics for variables (e.g., mean, standard deviation, skewness, kurtosis), store the results in `tibble` data frames, and optionally export them to disk.

-   **Unit Root Testing**: Perform a battery of unit root tests on integrated series (including ADF, PP, KPSS, ERS, ZA), present results in `tibble`s for easy export, allowing researchers to review all outcomes and select those relevant for their study.

-   **Model Estimation**: Estimate ARDL models using both OLS and quantile regression. Results are provided in both ARDL and Error Correction Model (ECM) forms, stored as `tibble`s, and can be exported.

-   **Slope Equality Testing**: Conduct tests for the equality of slopes across different quantiles or multiple variables within a specified quantile, utilizing Bootstrap and Wald statistics.

-   **Visualization**: Calculate and visualize QARDL model coefficients and dynamic multipliers.

### Suggestions and Bug Reports

Your feedback is welcome! If you have any suggestions or encounter bugs, please email me at [shin_li\@foxmail.com](shin_li@foxmail.com).

## ARDL at a Glance

An ARDL model, estimated using ordinary least squares (OLS), is a linear model that comprises two key components: the autoregressive part (AR) and the distributed lags (DL) of the independent variables. In the AR part, the dependent variable is considered in lagged form, while in the DL part, the independent variables are included in both their levels and lagged forms. Accordingly, an ARDL model denoted as *ARDL*(3, 0, 2) would include, up to three lags of the dependent variable, one independent variable considered at the level, and another independent variable included in both its level and lagged forms up to two lags.

ARDL formula:

$$
y_t=c_0+c_1 t+\sum_{i=1}^p b_{y, i} y_{t-i}+\sum_{j=1}^k \sum_{l=0}^{q_j} b_{j, l} x_{j, t-l}+\epsilon_t
$$

An ECM utilizes the first difference of the dependent variable, regressed on the first lags of both the dependent and independent variables. The remaining regressors in the model consist of the lags of the first differences of both the dependent and independent variables.

ECM formula:

$$
\begin{aligned}
\Delta y_t= & c_0+c_1 t+\pi_y y_{t-1}+\sum_{j=1}^k \pi_j x_{j, t-1}+\sum_{i=1}^{p-1} \psi_{y, i} \Delta y_{t-i} \\
& +\sum_{j=1}^k \sum_{l=1}^{q_j-1} \psi_{j, l} \Delta x_{j, t-l}+\sum_{j=1}^k \omega_j \Delta x_{j, t}+\epsilon_t
\end{aligned}
$$

Where $\psi_{j, l}=0 \forall q_j \leq 1, \psi_{y, i}=0$ if $p=1$, and $x_{j, t-1}$ and $\Delta x_{j, t}$ cancel out becoming $x_{j, t} \forall q_j=0$.

The hypothesis test for the bounds F-test is based on the ECM and the terms $c_0$ and $c_1$ are included only in cases 2 and 4, respectively, where this essentially means that they enter the long-run relationship. The hypothesis is the following:

$$
\begin{aligned}
& \mathbf{H}_{\mathbf{0}}: \pi_y=\pi_1=\cdots=\pi_k=c_0=c_1=0 \\
& \mathbf{H}_{\mathbf{1}}: \pi_y \neq \pi_1 \neq \cdots \neq \pi_k \neq c_0 \neq c_1 \neq 0
\end{aligned}
$$

The OLS-based ARDL model described above can be easily extended to quantile regression; the details of this extension are not elaborated further here.

## Installation

### Online Installation

You can install `QARDLFlex` online using the `install_github()` function from the **devtools**​ package. If you haven't installed the **devtools**​ package yet, you need to run `install.packages("devtools")` first.

```R
# install.packages("devtools")
devtools::install_github("Passenger53/QARDLFlex")
```
### Offline Installation

If have unstable network connectivity (which may be particularly relevant for users in Mainland China), you can choose to install `QARDLFlex` locally by first downloading its archive from GitHub. However, since `QARDLFlex` depends on several other R packages, you will need to install these dependencies first. The required packages are listed below for clarity.

`ARDL`, `dplyr`, `e1071`, `ggplot2`, `lubridate`, `magrittr`, `MASS`, `Matrix`, `purrr`, `quantreg`, `readr`, `reshape2`, `rio`, `stats`, `stringr`, `tibble`, `tidyplots`, `tidyr`, `tseries`, `urca`

Typically, you can install these packages using the `install.packages()` function. For example:

```R
install.packages("urca")
```
After installing the dependencies, you need to download the `QARDLFlex` archive from <https://github.com/Passenger53/QARDLFlex>, which is usually named `"QARDLFlex-main.zip"`. Unzip this file in your working directory to obtain a folder named `"QARDLFlex-main"`. Then, run the following command for offline installation:

```R
install.packages("QARDLFlex-main",
                 repos = NULL,
                 type = "source")
```
## Load `QARDLFlex`

Run the following code to library `QARDLFlex`:

```R
library(QARDLFlex)
```
## `ARDLFlex` Usage Examples

### Data

`ARDLFlex` has several specific requirements for the input data to ensure a smoother workflow for empirical analysis. Users must ensure that the first column of the data frame is a **date** column, and its name must be "**date**". The reason for this specification is that empirical research often requires plotting time series trends. If the date column is deliberately removed at the initial stage, subsequent plotting would necessitate re-importing a dataset containing the date column, which would be cumbersome.

We can use the built-in `sim_data_generate()` function to create a simulated dataset, which allows us to examine the structure of a typical dataset required by `ARDLFlex`.

```R
sim_data <- sim_data_generate(b = c(0.5, 0.7, 0.3, 0.2),
                             order = c(1, 1),
                             size = 100,
                             reps = 2)

head(sim_data[[1]])
```

```R
# A tibble: 6 × 3
  date            y     x1
  <date>      <dbl>  <dbl>
1 2017-01-01 -0.190  0.661
2 2017-02-01 -0.347 -2.05 
3 2017-03-01  0.955 -1.50 
4 2017-04-01  1.38   1.47 
5 2017-05-01  2.33   1.46 
6 2017-06-01  4.18   0.140
```

The code above generates a sample set based on an ARDL(1,1) model (`order = c(1, 1)`), with an intercept of 0.5, a coefficient of 0.7 for the lagged dependent variable, and coefficients of 0.3 and 0.2 for the contemporaneous and first lag of the explanatory variable, respectively. The parameter `reps = 2` specifies the generation of two datasets, and `size = 100` specifies a sample size of 100 for each dataset. Note that the first column is a `date` column, while the remaining columns are numeric series without `NA`s. Column names can be set by the user, but it is not recommended to include spaces, dots, or other special characters. Names that are too simple, such as "L1", "D1", "L0", or "D0", are also not allowed, as they have special meanings in `QARDLFlex`.

Users should note that when using the `sim_data_generate()` function, the length of parameter `b` must match the total number of parameters implied by the `order` argument. Specifically, the required length is: `sum(order + 1)`. For example, an ARDL(1,2) model requires ((1+1) + (2+1)) = 5 coefficients: intercept, y_lag1, x_contemporaneous, x_lag1, and x_lag2.

Below is an example of a simulated dataset corresponding to an ARDL(1,2) model:

```R
sim_data2 <- sim_data_generate(
  b = c(0.5, 0.6, 0.3, 0.2, 0.1),
  order = c(1, 2),
  size = 50,
  reps = 1
)
head(sim_data2[[1]])
```

```R
# A tibble: 6 × 3
  date            y     x1
  <date>      <dbl>  <dbl>
1 2021-01-01 -0.776  0.661
2 2021-02-01  0.636 -2.05 
3 2021-03-01  1.65  -1.50 
4 2021-04-01  1.49   1.47 
5 2021-05-01  2.11   1.46 
6 2021-06-01  3.96   0.140
```

The **QARDLFlex** package provides a `data_scale()` function that can be used to standardize data (z-score normalization). The `scale` argument is logical: if `TRUE`, scaling is applied to all numeric columns; if `FALSE`, the original data is returned unchanged.

```R
data_scale(sim_data[[1]], scale = T) %>% 
  head()
```

```R
# A tibble: 6 × 3
  date            y     x1
  <date>      <dbl>  <dbl>
1 2017-01-01 -1.50   0.714
2 2017-02-01 -1.61  -2.24 
3 2017-03-01 -0.672 -1.64 
4 2017-04-01 -0.366  1.60 
5 2017-05-01  0.316  1.58 
6 2017-06-01  1.65   0.148
```

### Summary Statistics

Taking the dataset `sim_data2` as an example, the `summary_statistics()` function can be used to compute descriptive statistics. Since all variables were originally simulated using normal random numbers, we square one of the variables (rendering it non-normal) to demonstrate the significance of the Jarque-Bera normality test:

```R
data <- sim_data[[1]]
data %>% dplyr::mutate(x1 = x1 ^ 2) -> data1
summary_statistics(data = data1,
                   export = F,
                   file = NULL)
```

```R
# A tibble: 2 × 9
  variable Mean  Median Maximum Minimum Std.Dev. Skewness Kurtosis `Jarque-Bera`
  <chr>    <chr> <chr>  <chr>   <chr>   <chr>    <chr>    <chr>    <chr>        
1 y        1.889 1.761  6.141   -1.173  1.389    0.362    0.068    2.325        
2 x1       0.836 0.400  9.218   0.000   1.235    3.692    19.823   1948.889***  
```

Here, one asterisk, two asterisks, and three asterisks denote significance levels of 10%, 5%, and 1%, respectively.

Users can choose to directly export the results to disk, for example:

```r
summary_statistics(data = data,
                   export = T,
                   file = "Summary Statistics.xlsx")
```
### Unit Root Tests

The `unitroot_test()` function integrates multiple unit root testing methods, supports various model specifications (e.g., with drift, with trend, etc.), and provides the optimal lag order selected by both the AIC and BIC information criteria. Additionally, the function reports unit root test results for both the original series and its differenced version.

```R
ur.result <- unitroot_test(data = data,
              export = F,
              file = NULL)
head(ur.result)
```

```R
# A tibble: 6 × 11
  variable difference test  selectlags type  statistic   lag `1pct` `5pct`
  <chr>    <chr>      <chr> <chr>      <chr> <chr>     <int>  <dbl>  <dbl>
1 y        level      ADF   BIC        none  -2.084**      1  -2.6   -1.95
2 y        level      ADF   AIC        none  -2.084**      1  -2.6   -1.95
3 y        level      ADF   BIC        drift -4.562***     1  -3.51  -2.89
4 y        level      ADF   AIC        drift -4.562***     1  -3.51  -2.89
5 y        level      ADF   BIC        trend -4.774***     1  -4.04  -3.45
6 y        level      ADF   AIC        trend -4.774***     1  -4.04  -3.45
# ℹ 2 more variables: `10pct` <dbl>, bpoint <int>
```

### Find Optimal ARDL Order

The `order_find()` function automatically selects or manually specifies the optimal lag order for ARDL model variables. It supports both automatic lag selection based on information criteria and manual specification of lag orders. The `selection` argument currently supports `"AIC"` (Akaike Information Criterion, which tends to select models with better predictive accuracy) and `"BIC"` (Bayesian Information Criterion, which favors more parsimonious models).

The `auto_lag` argument is logical. If `TRUE` (default), it automatically selects the optimal lags using the specified information criterion. If `FALSE`, it uses manually provided lag orders.

The `lags` argument is an integer vector specifying the manual lag orders for each variable. It is required when `auto_lag = FALSE`. Its length should match the number of variables (excluding the date variable).

```R
# Automatic lag selection using BIC criterion
optimal_lags <- order_find(data = data,
                           selection = "BIC",
                           export = FALSE)
```

```R
Optimal lag orders based on BIC criterion:
 y x1 
 1  1 
```

```R
# Manual lag specification
manual_lags <- order_find(data = data, 
                          auto_lag = FALSE, 
                          lags = c(2, 1), 
                          export = FALSE)
```

```R
Manually specified lag orders:
 y x1 
 2  1 
```

### Bounds F-test

The `bound_test()` function performs the bounds F-test for cointegration in ARDL models, testing the null hypothesis of no long-run relationship between the variables. The `order` argument is typically obtained from `order_find()`. Its length must match the number of variables in the data frame (excluding the date variable).

```R
bound_test(data = data, 
           order = optimal_lags, 
           export = FALSE)
```

```R

 Bounds F test with statistic 13.768 and p-value 0.000
```

### OLS based ARDL Estimation

The `ardl_estimate()` function can be used to estimate an OLS-based ARDL model. It returns a list containing six elements:

-   ARDL_estimate: The full ARDL model object from stats::lm()

-   ECM_estimate: The full ECM model object from stats::lm()

-   ARDL_summary: Summary of the ARDL model

-   ECM_summary: Summary of the ECM model

-   ARDL_coefficient: Formatted ARDL coefficient table with significance stars

-   ECM_coefficient: Formatted ECM coefficient table with significance stars

```R
ardl_results <- ardl_estimate(data = data, 
                         order = optimal_lags)
print(ardl_results$ARDL_coefficient)
```

```R
# A tibble: 9 × 3
  Variable    `Coef (Std. Error)` Value    
  <chr>       <chr>               <chr>    
1 (Intercept) Estimate            0.7736***
2 (Intercept) Std. Error          (0.1608) 
3 L1.y        Estimate            0.6062***
4 L1.y        Std. Error          (0.0702) 
5 L0.x1       Estimate            0.4295***
6 L0.x1       Std. Error          (0.1005) 
7 L1.x1       Estimate            0.2422** 
8 L1.x1       Std. Error          (0.1059) 
9 R2          <NA>                0.5745   
```

```R
print(ardl_results$ECM_coefficient)
```

```R
# A tibble: 9 × 3
  Variable    `Coef (Std. Error)` Value     
  <chr>       <chr>               <chr>     
1 (Intercept) Estimate            0.7736*** 
2 (Intercept) Std. Error          (0.1608)  
3 ECM         Estimate            -0.3938***
4 ECM         Std. Error          (0.0702)  
5 L1.x1       Estimate            1.7059*** 
6 L1.x1       Std. Error          (0.3725)  
7 D0.x1       Estimate            0.4295*** 
8 D0.x1       Std. Error          (0.1005)  
9 R2          <NA>                0.3290    
```

### QARDL Estimation

The `qardl_estimate()` function can be used to estimate QARDL and QECM models with bootstrapped confidence intervals. It provides comprehensive results, including bootstrap samples, estimated coefficients with standard errors, significance stars, and pseudo R-squared values. The `data` and `order` arguments are consistent with those in the aforementioned functions. The key difference is that the `tau` argument specifies the quantile(s) to be estimated, with a default value of 0.5, and the `boots` argument specifies the number of bootstrap replications, which defaults to 200.

```R
qardl_results <- qardl_estimate(data = data, 
                                order = optimal_lags, 
                                tau = 0.5, 
                                boots = 50)
print(qardl_results$QARDL_coefficient)
```

```R
# A tibble: 9 × 3
  Variable    `Coef (Std. Error)` Value    
  <chr>       <chr>               <chr>    
1 (Intercept) Estimate            0.7325***
2 (Intercept) Std. Error          (0.1797) 
3 L1.y        Estimate            0.6228***
4 L1.y        Std. Error          (0.0824) 
5 L0.x1       Estimate            0.4927***
6 L0.x1       Std. Error          (0.1312) 
7 L1.x1       Estimate            0.2504*  
8 L1.x1       Std. Error          (0.1463) 
9 R2          <NA>                0.3511   
```

```R
print(qardl_results$QECM_coefficient)
```

```R
# A tibble: 9 × 3
  Variable    `Coef (Std. Error)` Value     
  <chr>       <chr>               <chr>     
1 (Intercept) Estimate            0.7325*** 
2 (Intercept) Std. Error          (0.1797)  
3 ECM         Estimate            -0.3772***
4 ECM         Std. Error          (0.0824)  
5 L1.x1       Estimate            1.9701*** 
6 L1.x1       Std. Error          (0.6310)  
7 D0.x1       Estimate            0.4927*** 
8 D0.x1       Std. Error          (0.1312)  
9 R2          <NA>                0.1918    
```

### Comprehensive estimation results

The `coef_estimate()` function may be more commonly used in empirical research than `ardl_estimate()` and `qardl_estimate()`, as it can simultaneously estimate both the OLS-based ARDL model and QARDL models at multiple quantiles, and it also allows exporting the results directly to disk.

```R
all_results <- coef_estimate(
  data = data,
  order = optimal_lags,
  taus = seq(0.1, 0.9, 0.1),
  boots = 50,
  export = T,
  file_ardl = "ARDL estimation results.xlsx",
  file_ecm = "ECM estimation results.xlsx"
)
```

Note that, for quantile regression, the R² value in the bottom row refers to the pseudo R².

### Plot QARDL Coefficients across Quantiles

The `plot_ECM()` function can be used to visualize the long-term and short-term coefficient trends for ECM model variables across specified quantile levels. For short-term coefficients, if a variable has multiple lags, the sum of the coefficients across all lags is plotted.

-   The `estimation` argument takes the estimation results obtained from `coef_estimate()`.

-   The `plot_term` argument specifies whether to plot `"short"` (short-term) or `"long"` (long-term) coefficients. Only the values `"short"` and `"long"` are allowed.

-   The `variable` argument specifies the name of the variable to plot. It must match one of the variable names in the data. Special cases:

    -   If the variable is the dependent variable and `plot_term = "long"`, the ECM coefficient is plotted by default.

    -   If the variable is an independent variable with an ARDL order of zero, its contemporaneous effect coefficient is plotted.

-   The `taus` argument must be the same as the `taus` parameter used in `coef_estimate()`.

-   The `ci` argument specifies whether to plot confidence intervals. The default is `TRUE`.

-   The `level` argument sets the confidence level. Only `0.9`, `0.95`, and `0.99` are allowed, with a default of `0.95`.

The `plot_ECM()` function returns a list containing two elements:

-   `plot`: a `ggplot` object of the coefficient trend plot

-   `plot_data`: the underlying data used for plotting, enabling users to create more flexible custom visualizations.

```R
plot_result <- plot_ECM(
  estimation = all_results,
  order = optimal_lags,
  plot_term = "long",
  variable = names(optimal_lags)[1],
  taus = seq(0.1, 0.9, 0.1),
  ci = TRUE,
  level = 0.95,
  export = FALSE
)
```

<img width="1344" height="960" alt="image" src="https://github.com/user-attachments/assets/6b7ce5a0-a777-4995-95cf-d46af019850e" />


### Cross-Quantile Wald Tests

The `cross_quantile_test()` function performs a Wald test to compare coefficients across different quantiles within the Error Correction Model (ECM) framework. It can test whether coefficients are equal across specified quantile levels—either individually or jointly for multiple variables. This test is particularly useful for assessing parameter stability across the conditional distribution.

-   **`estimation`**: Estimation results obtained from `coef_estimate()`.

-   **`order`**: Same as in previous functions.

-   **`boots`**: Same as in previous functions (number of bootstrap replications).

-   **`taus`**: Same as in previous functions (quantile levels used in estimation).

-   **`test_taus`**: A numeric vector specifying the quantile levels to be tested for coefficient equality. Must be a subset of `taus` and contain at least two distinct quantiles.

-   **`variable`**: A character vector listing the names of variables to test. Acceptable values include `"ECM"`, `"Intercept"`, or regular variable names. Special cases:

    -   If testing the dependent variable with `term = "long"`, the test applies to the ECM coefficient.

    -   If testing a variable with an ARDL order of zero, the test applies to its contemporaneous effect coefficient.

    -   **Note**: `"ECM"` or `"Intercept"` cannot be combined with regular variables in the same test.
 
    -   **If length of `variable` larger than 1**. For example, `variable = c("x1", "x2")`, then this function will test **H0: beta_1(tau1) = beta_1(tau2), beta_2(tau1) = beta_2(tau2)**

-   **`term`**: A character string. For regular variables, specify `"long"` to test long-term coefficients or `"short"` to test short-term coefficients. This argument is only applicable when testing regular variables (i.e., not `"ECM"` or `"Intercept"`).

-   **`joint`**: Logical. For short-term coefficients, if `TRUE` (default), the test evaluates the sum of coefficients across all lags; if `FALSE`, it tests each lagged coefficient separately.

```{.r .cell-code}
# Test if ECM coefficient is equal across quantiles 0.2, 0.5, and 0.8
test_result <- cross_quantile_test(
  estimation = all_results,
  order = optimal_lags,
  boots = 50,
  taus = seq(0.1, 0.9, 0.1),
  test_taus = c(0.2, 0.5, 0.8),
  variable = "ECM",
  joint = TRUE
)

print(test_result)
```

```

	Wald Test

data:  
= 0.72015, = 2, p-value = 0.6976
```

Compared to the `cross_quantile_test()` function, `single_variable_cross_wald()` is likely more commonly used. It performs Wald tests for equality of coefficients across specified quantiles **for each variable** in the ECM model. Specifically, it tests long-run and short-run coefficients separately for every variable, with special handling for:

-   The dependent variable (i.e., the ECM term), and

-   Variables with an ARDL order of zero (for which the contemporaneous effect is tested).

For short-run coefficients involving multiple lags, the test is conducted on the **sum of the coefficients across all lags**.

```{.r .cell-code}
# Perform cross-quantile Wald tests for all variables
wald_results <- single_variable_cross_wald(
  estimation = all_results,
  order = optimal_lags,
  boots = 50,
  taus = seq(0.1, 0.9, 0.1),
  test_taus = c(0.2, 0.5, 0.8),
  export = FALSE
)

# View results
print(wald_results)
```

```
# A tibble: 3 × 5
  Variable term  Wald   pvalue Tau           
  <chr>    <chr> <chr>  <chr>  <chr>         
1 ECM      <NA>  0.7201 0.6976 0.20=0.50=0.80
2 x1       Long  0.0423 0.9791 0.20=0.50=0.80
3 x1       Short 0.1762 0.9157 0.20=0.50=0.80
```

### Within-Quantile Wald Tests

The `within_quantile_test()` function performs Wald tests to assess whether the coefficients of multiple explanatory variables are equal **within each quantile level**. It conducts separate tests for both long-run and short-run coefficients at each quantile. This is particularly useful for NARDL-type models, where variables are decomposed into multiple components (e.g., positive and negative partial sums).

Because this test involves multiple independent variables, we simulate a new multi-variable dataset:

```{.r .cell-code}
# Generate sample data
sample_data <- sim_data_generate(
  b = c(0.5, 0.7, 0.3, 0.4, 0.1, 0.2),
  order = c(1, 1, 1),
  size = 100,
  reps = 1
)[[1]]

# Find optimal lag order
optimal_order <- order_find(sample_data, 
                            selection = "BIC", 
                            export = FALSE)
```

```
Optimal lag orders based on BIC criterion:
 y x1 x2 
 1  1  1 
```

```{.r .cell-code}
# Estimate models at multiple quantiles
est_result <- coef_estimate(
  data = sample_data,
  order = optimal_order,
  taus = c(0.25, 0.5, 0.75),
  boots = 50,
  export = FALSE
)
```

```{.r .cell-code}
# Test if coefficients of two explanatory variables are equal at each quantile
test_results <- within_quantile_test(
  estimation = est_result,
  order = optimal_order,
  boots = 50,
  taus = c(0.25, 0.5, 0.75),
  variables = c("x1", "x2"),
  export = FALSE
)

# View results
print(test_results)
```

```
# A tibble: 6 × 5
  Variable term  Wald    pvalue Tau  
  <chr>    <chr> <chr>   <chr>  <chr>
1 x1=x2    Long  0.6387  0.4242 0.25 
2 x1=x2    Short 3.4976* 0.0615 0.25 
3 x1=x2    Long  0.2009  0.6540 0.50 
4 x1=x2    Short 0.6827  0.4087 0.50 
5 x1=x2    Long  0.0735  0.7863 0.75 
6 x1=x2    Short 0.0102  0.9197 0.75 
```

### Calculate Dynamic Multipliers

The `all_dynamic_multipliers()` function computes dynamic multipliers for all explanatory variables in an ARDL model across specified quantile levels. It extends the `dynamic_multipliers()` function by automatically handling multiple variables and quantiles, offering a comprehensive analysis of how exogenous shocks propagate through the system over time. Please refer to the function’s help documentation for further details.

```{.r .cell-code}
# Generate sample data using sim_data_generate()
sample_data <- sim_data_generate(
  b = c(0.5, 0.7, 0.3, 0.2),
  order = c(1, 1),
  size = 100,
  reps = 1
)[[1]]

# Find optimal lag order using order_find()
optimal_order <- order_find(sample_data, selection = "BIC", export = FALSE)
```

```
Optimal lag orders based on BIC criterion:
 y x1 
 1  1 
```

```{.r .cell-code}
# Calculate dynamic multipliers for all explanatory variables
results <- all_dynamic_multipliers(
  data = sample_data,
  order = optimal_order,
  taus = seq(0.1, 0.9, 0.1),
  steps = 10,
  cumulative = FALSE,
  impulse = FALSE,
  export = FALSE
)

# View results
head(results)
```

```
# A tibble: 6 × 4
  variable tau    step value
  <chr>    <chr> <dbl> <dbl>
1 x1       ols       0 0.430
2 x1       ols       1 0.932
3 x1       ols       2 1.24 
4 x1       ols       3 1.42 
5 x1       ols       4 1.53 
6 x1       ols       5 1.60 
```

### Plot Dynamic Multipliers

**QARDLFlex** provides two ways to visualize dynamic multipliers. The `plot_line_multipliers()` function can be used to draw dynamic multiplier line charts, while the `plot_heatmap_multipliers()` function can be used to draw dynamic multiplier heat maps. In both functions, the `multipliers` parameter is the estimation result of the `all_dynamic_multipliers()` function. The generated plots are based on the `tidyplots` package, allowing users to further modify them as needed. In fact, given the estimation results from `all_dynamic_multipliers()` (in the form of a tibble data frame), users can easily perform their own visualizations of dynamic multipliers.

```{.r .cell-code}
line_plot <- plot_line_multipliers(
  multipliers = results,
  plot_variable = "x1", 
  plot_taus = c(0.2, 0.5, 0.8),
  export = FALSE
)
```

<img width="1344" height="960" alt="image" src="https://github.com/user-attachments/assets/73590a0c-8d0a-4f17-897b-ded1b79082cd" />


```{.r .cell-code}
heatmap_plot <- plot_heatmap_multipliers(
  multipliers = results,
  plot_variable = "x1",  # Replace with actual variable name from your data
  export = FALSE
)
```

<img width="1344" height="960" alt="image" src="https://github.com/user-attachments/assets/876cc720-eac1-4815-a987-da7a26be4b9c" />


## Further Extensions

As a small illustration, the `generate_cumulative_shocks()` function offers a method for time series decomposition. It computes the first differences of a specified series and creates two new variables:

- **`p_variable_name`**: Cumulative positive shocks (increases relative to the previous period)  
- **`n_variable_name`**: Cumulative negative shocks (decreases relative to the previous period)

These new columns are inserted immediately after their corresponding original variables in the data frame.

This decomposition allows users to construct a dataset containing the decomposed series and then estimate a QARDL model. Using the previously introduced suite of Wald test functions, users can examine whether the long-run and short-run coefficients associated with `p_variable_name` and `n_variable_name` are statistically equal—thereby generating richer empirical insights.

We will not repeat the full modeling, testing, and visualization workflow here; instead, we only demonstrate the usage of the `generate_cumulative_shocks()` function:

```{.r .cell-code}
# Generate sample data using sim_data_generate()
sample_data <- sim_data_generate(
  b = c(0.5, 0.7, 0.3, 0.2),
  order = c(1, 1),
  size = 100,
  reps = 1
)[[1]]

str(sample_data)
```

```
tibble [100 × 3] (S3: tbl_df/tbl/data.frame)
 $ date: Date[1:100], format: "2017-01-01" "2017-02-01" ...
 $ y   : num [1:100] -0.19 -0.347 0.955 1.381 2.328 ...
 $ x1  : num [1:100] 0.661 -2.053 -1.499 1.471 1.459 ...
```

```{.r .cell-code}
# Generate cumulative shocks for specific variables
processed_data <- generate_cumulative_shocks(
  variable_names = "x1",  # Replace with actual variable names
  data = sample_data
)

# View the modified data structure
str(processed_data)
```

```
tibble [100 × 5] (S3: tbl_df/tbl/data.frame)
 $ date: Date[1:100], format: "2017-01-01" "2017-02-01" ...
 $ y   : num [1:100] -0.19 -0.347 0.955 1.381 2.328 ...
 $ x1  : num [1:100] 0.661 -2.053 -1.499 1.471 1.459 ...
 $ p_x1: num [1:100] 0 0 0.554 3.524 3.524 ...
 $ n_x1: num [1:100] 0 -2.71 -2.71 -2.71 -2.73 ...
```


