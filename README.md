---
title: "The QARDLFlex Package"
author: 
  - name: "Xin Li"
    email: "shin_li@foxmail.com"
    affiliation: "School of Economics, Qingdao University"
format: 
  html:
    embed-resources: true
    keep-md: true
    toc: true
    number-sections: true
    toc-depth: 4
    number-depth: 4
    number-offset: 1
editor: visual
---






## Introduction

### What is `QARDLFlex`?

A mature and stable implementation for the Quantile Autoregressive Distributed Lag (QARDL) model has been lacking. To address this gap, I developed `QARDLFlex`, an R package designed to provide a flexible and reliable toolkit for QARDL-related analysis. Although `QARDLFlex` is currently in its early development stages, it already offers a complete workflow for empirical research—from descriptive statistics and QARDL model estimation to parameter hypothesis testing.

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



::: {.cell}

```{.r .cell-code}
# install.packages("devtools")
devtools::install_github("Passenger53/QARDLFlex")
```
:::



### Offline Installation

If have unstable network connectivity (which may be particularly relevant for users in Mainland China), you can choose to install `QARDLFlex` locally by first downloading its archive from GitHub. However, since `QARDLFlex` depends on several other R packages, you will need to install these dependencies first. The required packages are listed below for clarity.

-   `ARDL`, `dplyr`, `e1071`, `ggplot2`, `lubridate`, `magrittr`, `MASS`, `Matrix`, `purrr`, `quantreg`, `readr`, `reshape2`, `rio`, `stats`, `stringr`, `tibble`, `tidyplots`, `tidyr`, `tseries`, `urca`

Typically, you can install these packages using the `install.packages()` function. For example:



::: {.cell}

```{.r .cell-code}
install.packages("urca")
```
:::



After installing the dependencies, you need to download the `QARDLFlex` archive from <https://github.com/Passenger53/QARDLFlex>, which is usually named `"QARDLFlex-main.zip"`. Unzip this file in your working directory to obtain a folder named `"QARDLFlex-main"`. Then, run the following command for offline installation:



::: {.cell}

```{.r .cell-code}
install.packages("QARDLFlex-main",
                 repos = NULL,
                 type = "source")
```
:::



## Load `QARDLFlex`

Run the following code to library `QARDLFlex`:



::: {.cell}

```{.r .cell-code}
library(QARDLFlex)
```
:::



## `ARDLFlex` Usage Examples

### Data

`ARDLFlex` has several specific requirements for the input data to ensure a smoother workflow for empirical analysis. Users must ensure that the first column of the data frame is a **date** column, and its name must be "**date**". The reason for this specification is that empirical research often requires plotting time series trends. If the date column is deliberately removed at the initial stage, subsequent plotting would necessitate re-importing a dataset containing the date column, which would be cumbersome.

We can use the built-in `sim_data_generate()` function to create a simulated dataset, which allows us to examine the structure of a typical dataset required by `ARDLFlex`.



::: {.cell}

```{.r .cell-code}
sim_data <- sim_data_generate(b = c(0.5, 0.7, 0.3, 0.2),
                             order = c(1, 1),
                             size = 100,
                             reps = 2)

head(sim_data[[1]])
```

::: {.cell-output .cell-output-stdout}

```
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


:::
:::



The code above generates a sample set based on an ARDL(1,1) model (`order = c(1, 1)`), with an intercept of 0.5, a coefficient of 0.7 for the lagged dependent variable, and coefficients of 0.3 and 0.2 for the contemporaneous and first lag of the explanatory variable, respectively. The parameter `reps = 2` specifies the generation of two datasets, and `size = 100` specifies a sample size of 100 for each dataset. Note that the first column is a `date` column, while the remaining columns are numeric series without `NA`s. Column names can be set by the user, but it is not recommended to include spaces, dots, or other special characters. Names that are too simple, such as "L1", "D1", "L0", or "D0", are also not allowed, as they have special meanings in `QARDLFlex`.

Users should note that when using the `sim_data_generate()` function, the length of parameter `b` must match the total number of parameters implied by the `order` argument. Specifically, the required length is: `sum(order + 1)`. For example, an ARDL(1,2) model requires ((1+1) + (2+1)) = 5 coefficients: intercept, y_lag1, x_contemporaneous, x_lag1, and x_lag2.

Below is an example of a simulated dataset corresponding to an ARDL(1,2) model:



::: {.cell}

```{.r .cell-code}
sim_data2 <- sim_data_generate(
  b = c(0.5, 0.6, 0.3, 0.2, 0.1),
  order = c(1, 2),
  size = 50,
  reps = 1
)
head(sim_data2)
```

::: {.cell-output .cell-output-stdout}

```
[[1]]
# A tibble: 50 × 3
   date            y     x1
   <date>      <dbl>  <dbl>
 1 2021-01-01 -0.776  0.661
 2 2021-02-01  0.636 -2.05 
 3 2021-03-01  1.65  -1.50 
 4 2021-04-01  1.49   1.47 
 5 2021-05-01  2.11   1.46 
 6 2021-06-01  3.96   0.140
 7 2021-07-01  3.57   0.209
 8 2021-08-01  0.524 -3.04 
 9 2021-09-01 -0.605 -0.487
10 2021-10-01 -1.04  -1.09 
# ℹ 40 more rows
```


:::
:::
