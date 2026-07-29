pkgs <- c(
  "data.table", "dplyr", "tidyr", "purrr", "stringr", "lubridate",
  "xgboost", "mgcv", "pROC", "PRROC",
  "ggplot2", "cowplot", "ggpubr", "gridExtra", "patchwork", "scales",
  "writexl", "openxlsx", "zip"
)

new_pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(new_pkgs)) install.packages(new_pkgs)
