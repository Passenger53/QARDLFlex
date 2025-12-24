
# 使用管道
usethis::use_pipe()
usethis::use_package("stats")
usethis::use_package("MASS")
usethis::use_package("stringr")
usethis::use_package("tibble")
usethis::use_package("lubridate")
usethis::use_package("dplyr")
usethis::use_package("purrr")
usethis::use_package("e1071")
usethis::use_package("tseries")
usethis::use_package("urca")
usethis::use_package("rio")
usethis::use_package("ARDL")
usethis::use_package("quantreg")
usethis::use_package("tidyr")
usethis::use_package("ggplot2")
usethis::use_package("reshape2")
usethis::use_package("Matrix")
usethis::use_package("readr")
usethis::use_package("tidyplots")

devtools::document()
devtools::build()

library(QARDLFlex)
help(plot_line_multipliers)



# 【使用 git 命令行上传项目到 GitHub（以 R 包为例）】 https://www.bilibili.com/video/BV1Zu4y1i7uN/?share_source=copy_web&vd_source=3b82cc67f255864d8b26ae8a0de7e0b0

# Terminal窗口：
# git init
# git add .
# git commit -m "First submit for test 2025-12-24"
# git branch -M main
# git remote add origin git@github.com:Passenger53/QARDLFlex.git
# git push -u origin main

devtools::install_github("Passenger53/QARDLFlex")




