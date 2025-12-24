
rm(list = ls())
gc()
options(warn = -1)
Sys.setlocale("LC_TIME", "English")

library(QARDLFlex)
library(bruceR)
library(tidyverse)

# 数据
data <- import("example_data.xlsx", as = "tibble")
data %>% select(date, VIX:gold) -> data
data %>% mutate(date = ymd(date)) -> data
check_data(data) # 检查数据是否符合要求

# 参数
taus <- seq(0.1, 0.9, 0.1) # 设定待估计的模型分位点
test_taus <- taus[c(1,5,9)] # 设定拟检验的分位点（需要是taus的子集）
criterion <- "AIC" # 设定最佳滞后阶数计算准则
auto_lag <- T # 是（TRUE）否（FALSE）基于信息准则选择最佳滞后期
lags <- c(1,1,1,1) # 若auto_lag为FALSE，则使用使用该手动设定的滞后期
boots <- 100 # 分位数回归参数估计自助法重抽样次数
scale <- F # 是（TRUE）否（FALSE）对数据进行标准化处理
steps <- 20 # 动态乘数预测步长

# 创建结果存放目录
if (dir.exists("results")){
  unlink("results", recursive = T)
}
dir.create("results")

# 基本描述性统计
basic_statistic <- summary_statistics(
  data = data,
  export = T,
  file = "results/表01-基本描述性统计.xlsx"
)
basic_statistic

# 单位根检验
ur.result <- unitroot_test(data = data,
                           export = TRUE,
                           file = "results/表02-单位根检验.xlsx")
ur.result

# 数据标准化
data %>% data_scale(scale = scale) -> data

# 寻找最佳滞后期
order <- order_find(data = data,
                    selection = criterion,
                    auto_lag = auto_lag,
                    lags = lags,
                    export = TRUE,
                    file = "results/表03-模型滞后阶数.xlsx")

# 边界协整检验
bound_test(data = data,
           order = order,
           export = TRUE,
           file = "results/表04-边界协整F检验.xlsx")

# ARDL估计（演示）
ardl_results <- ardl_estimate(data = data, order = order)
ardl_results$ARDL_coefficient
ardl_results$ECM_coefficient

# 单一分位点的QARDL估计（演示）
qardl_results <- qardl_estimate(data = data,
                                order = order,
                                tau = 0.5,
                                boots = boots)
qardl_results$QARDL_coefficient
qardl_results$QECM_coefficient
qardl_results$QARDL_coef
qardl_results$QECM_coef

# 多个分位点的QARDL估计（导出OLS和所有分位点结果）
estimation <- coef_estimate(data = data,
              order = order,
              taus = taus,
              boots = boots,
              export = TRUE,
              file_ardl = "results/表05-QARDL估计结果（ARDL形态）.xlsx",
              file_ecm = "results/表06-QARDL估计结果（ECM形态）.xlsx")

# 绘制系数图
plot_coef <- plot_ECM(
  estimation,
  order,
  plot_term = "long",
  variable = "VIX",
  taus,
  ci = T,
  level = 0.9,
  export = F)

# （演示）批量绘制所有系数图
x1 <- c("VIX", "bitcoin", "bitcoin", "GreenBond", "GreenBond", "gold")
x2 <- c("long", "long", "short", "long", "short", "long")
x3 <- c("ECM", "长期", "短期", "长期", "短期", "当期")
x3 <- str_c("results/图01-", x1, "-", x3, "系数走势.png")

plot_coefs <- list(x1, x2, x3) %>%
  pmap(function(x1, x2, x3){
    plot_ECM(
      estimation,
      order,
      plot_term = x2,
      variable = x1,
      taus,
      ci = T,
      level = 0.9,
      export = T,
      file = x3)
  })

# （演示）检验给定不同分位点（两个或多个均可）ECM模型中的系数是否相等
cross_quantile_test(
  estimation = estimation,
  order = order,
  boots = boots,
  taus = taus,
  test_taus = taus[c(1,5,9)],
  variable = c("bitcoin"),
  term = "short",
  joint = F # 不同滞后期短期系数是否加总后再检验
)

# 所有单变量长短期系数跨指定分位数相等性Wald检验
cross_single <- single_variable_cross_wald(
  estimation = estimation,
  order = order,
  boots = boots,
  taus = taus,
  test_taus = test_taus,
  export = T,
  file = "results/表07a-单变量跨指定分位数Wald检验.xlsx"
)

combn(taus, 2) %>% as_tibble() %>% as.list() %>%
  map_dfr(~ single_variable_cross_wald(
    estimation = estimation,
    order = order,
    boots = boots,
    taus = taus,
    test_taus = .,
    export = F
  )) %>%
  export("results/表07b-单变量跨所有分位数（两个）Wald检验.xlsx")

combn(taus, 3) %>% as_tibble() %>% as.list() %>%
  map_dfr(~ single_variable_cross_wald(
    estimation = estimation,
    order = order,
    boots = boots,
    taus = taus,
    test_taus = .,
    export = F
  )) %>%
  export("results/表07c-单变量跨所有分位数（三个）Wald检验.xlsx")

# 指定几个变量的分位数内系数相等性检验
within_wald <- within_quantile_test(
  estimation = estimation,
  order = order,
  boots = boots,
  taus = taus,
  variables = c("bitcoin", "GreenBond"),
  export = T,
  file = "results/表08-多变量分位数内Wald检验.xlsx"
)
within_wald

# 所有可能的双变量分位数内检验结果（演示）
order[-1][order[-1] > 0] %>% names() %>%
  combn(2) %>%
  as_tibble() %>%
  as.list() %>%
  map_dfr(
    ~ within_quantile_test(
      estimation = estimation,
      order = order,
      boots = boots,
      taus = taus,
      variables = .,
      export = F
    )
  )

# 计算动态乘数（包括OLS和分位数回归）
multipliers <- all_dynamic_multipliers(
  data = data,
  order = order,
  taus = taus,
  steps = steps,
  cumulative = FALSE,
  impulse = FALSE,
  export = TRUE,
  file = "results/表09-动态乘数计算结果.xlsx"
)
multipliers

# 绘制动态乘数折线图
plot_line_multipliers(multipliers,
                      plot_variable = c("bitcoin"),
                      plot_taus = c(0.2, 0.5, 0.8),
                      export = FALSE,
                      file = NULL)

# 绘制动态乘数热图
plot_heatmap_multipliers(multipliers,
                      plot_variable = c("bitcoin"),
                      export = FALSE,
                      file = NULL)









