# 1. 清理环境并设置工作目录
rm(list=ls())

# 设置工作目录，请根据实际情况修改
setwd("E:/FILE/Masterexp/Paper/1_mother_choice/1_data/GLMM/alga/")

# 2. 加载必要的包
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(nlme)
library(glmmTMB)
library(DHARMa)
library(car)
library(ggplot2)
library(openxlsx)
library(cowplot)

# 3. 读取数据
density_M <- read.csv("density_M.csv", header = TRUE, stringsAsFactors = FALSE)
mixs_M <- read.csv("mixsrate_M.csv", header = TRUE, stringsAsFactors = FALSE)
reproductive_M <- read.csv("reproductive_M.csv", header = TRUE, stringsAsFactors = FALSE)
survival_M <- read.csv("survival_M.csv", header = TRUE, stringsAsFactors = FALSE)

# 检查数据结构
cat("检查数据结构:\n")
cat("density_M 维度:", dim(density_M), "\n")
cat("mixs_M 维度:", dim(mixs_M), "\n")
cat("reproductive_M 维度:", dim(reproductive_M), "\n")
cat("survival_M 维度:", dim(survival_M), "\n")

# 4. 定义通用函数用于数据预处理和分析（修复版）
process_dataset_glmm_fixed <- function(data, dataset_name) {
  cat(sprintf("\n=== 处理数据集: %s ===\n", dataset_name))
  
  # 添加group列（P和F的组合）
  if(!"group" %in% colnames(data)) {
    data$group <- paste0(data$P, data$F)
  }
  
  # 将数据转换为长格式
  time_cols <- names(data)[grepl("^t\\d+", names(data))]
  
  df_long <- data %>%
    pivot_longer(
      cols = all_of(time_cols),
      names_to = "time",
      values_to = "value"
    ) %>%
    mutate(
      time = as.numeric(gsub("t", "", time)),
      group = factor(group, levels = c("MC", "CC", "MM", "CM")),
      P = factor(P, levels = c("C", "M")),
      F = factor(F, levels = c("C", "M")),
      id = factor(id)
    ) %>%
    filter(!is.na(value))  # 移除缺失值
  
  # 为繁殖和生存数据创建小时单位的时间变量
  if (dataset_name %in% c("繁殖(reproductive)", "生存(survival)")) {
    df_long <- df_long %>%
      mutate(time_hours = time * 12)  # 将天转换为小时（乘以12）
  }
  
  # 确定数据类型并选择合适模型
  cat("数据特征分析:\n")
  cat(sprintf("  样本数量: %d\n", nrow(df_long)))
  cat(sprintf("  最小值: %.4f\n", min(df_long$value, na.rm = TRUE)))
  cat(sprintf("  最大值: %.4f\n", max(df_long$value, na.rm = TRUE)))
  cat(sprintf("  平均值: %.4f\n", mean(df_long$value, na.rm = TRUE)))
  cat(sprintf("  标准差: %.4f\n", sd(df_long$value, na.rm = TRUE)))
  
  # 检查是否为计数数据
  is_count_data <- all(df_long$value == floor(df_long$value), na.rm = TRUE)
  
  # 检查是否为比例数据
  is_proportion_data <- all(df_long$value >= 0 & df_long$value <= 1, na.rm = TRUE)
  
  cat(sprintf("  是否为计数数据: %s\n", ifelse(is_count_data, "是", "否")))
  cat(sprintf("  是否为比例数据: %s\n", ifelse(is_proportion_data, "是", "否")))
  
  # 初始化变量
  model <- NULL
  model_type <- "未知GLMM"
  has_random_slope <- FALSE  # 明确初始化为逻辑值
  
  # 简化的模型拟合函数
  fit_model <- function(formula_str, family_obj, use_slope = TRUE) {
    tryCatch({
      model_fit <- glmmTMB(as.formula(formula_str), 
                           data = df_long, 
                           family = family_obj,
                           control = glmmTMBControl(
                             optimizer = optim,
                             optArgs = list(method = "BFGS")
                           ))
      
      if (model_fit$fit$convergence == 0) {
        return(list(model = model_fit, success = TRUE))
      } else {
        return(list(model = model_fit, success = FALSE))
      }
    }, error = function(e) {
      return(list(model = NULL, success = FALSE))
    })
  }
  
  # 根据数据类型选择模型（简化逻辑）
  if (is_count_data) {
    cat("选择负二项GLMM模型处理计数数据\n")
    
    # 先尝试随机截距模型（更稳定）
    result <- fit_model("value ~ P * F + (1|id)", nbinom2(), use_slope = FALSE)
    
    if (result$success) {
      model <- result$model
      model_type <- "负二项GLMM"
      has_random_slope <- FALSE
    } else {
      # 尝试泊松分布
      result <- fit_model("value ~ P * F + (1|id)", poisson(), use_slope = FALSE)
      if (result$success) {
        model <- result$model
        model_type <- "泊松GLMM"
        has_random_slope <- FALSE
      }
    }
    
  } else if (is_proportion_data) {
    cat("选择Beta GLMM模型处理比例数据\n")
    
    # 调整边界值
    df_long <- df_long %>%
      mutate(value_adj = case_when(
        value == 0 ~ 0.001,
        value == 1 ~ 0.999,
        TRUE ~ value
      ))
    
    # 尝试随机截距模型
    result <- fit_model("value_adj ~ P * F + (1|id)", beta_family(), use_slope = FALSE)
    
    if (result$success) {
      model <- result$model
      model_type <- "Beta GLMM"
      has_random_slope <- FALSE
    } else {
      # 尝试高斯分布+logit链接
      result <- fit_model("value ~ P * F + (1|id)", gaussian(link = "logit"), use_slope = FALSE)
      if (result$success) {
        model <- result$model
        model_type <- "高斯GLMM(logit链接)"
        has_random_slope <- FALSE
      }
    }
    
  } else {
    cat("选择高斯GLMM模型处理连续数据\n")
    
    # 直接使用高斯分布
    result <- fit_model("value ~ P * F + (1|id)", gaussian(), use_slope = FALSE)
    
    if (result$success) {
      model <- result$model
      model_type <- "高斯GLMM"
      has_random_slope <- FALSE
    }
  }
  
  # 如果所有尝试都失败，使用简化模型
  if (is.null(model)) {
    cat("所有GLMM尝试失败，使用简化模型\n")
    model <- glmmTMB(value ~ group + (1|id), data = df_long, family = gaussian())
    model_type <- "简化GLMM"
    has_random_slope <- FALSE
  }
  
  # 添加随机效应信息到模型类型
  if (has_random_slope) {
    model_type <- paste(model_type, "(随机斜率)")
  } else {
    model_type <- paste(model_type, "(随机截距)")
  }
  
  cat(sprintf("\n最终使用模型: %s\n", model_type))
  
  # 双因素方差分析
  cat(sprintf("\n=== %s: 双因素方差分析结果 ===\n", dataset_name))
  
  anova_result <- Anova(model, type = "III")
  
  # 多重比较
  cat(sprintf("\n=== %s: 多重比较(Tukey HSD) ===\n", dataset_name))
  
  # 创建组合变量
  df_long <- df_long %>%
    mutate(P_F = interaction(P, F, sep = ""))
  
  # 拟合多重比较模型（简化：使用随机截距）
  cat("创建用于多重比较的模型...\n")
  
  # 确定响应变量
  if (model_type == "Beta GLMM (随机截距)") {
    resp_var <- "value_adj"
  } else {
    resp_var <- "value"
  }
  
  # 使用随机截距模型进行多重比较（更稳定）
  formula_mc <- as.formula(paste(resp_var, "~ P_F + (1|id)"))
  
  # 使用与原模型相同的分布
  family_used <- family(model)
  
  # 拟合多重比较模型
  model_mc <- tryCatch({
    glmmTMB(formula_mc, 
            data = df_long, 
            family = family_used)
  }, error = function(e) {
    cat("多重比较模型拟合失败，使用高斯分布替代\n")
    glmmTMB(formula_mc, data = df_long, family = gaussian())
  })
  
  # 计算边际均值
  emm_group <- emmeans(model_mc, ~ P_F, type = "response")
  
  # 重新排序P_F的levels为MC, CC, MM, CM
  current_levels <- levels(emm_group@grid$P_F)
  new_levels <- c("MC", "CC", "MM", "CM")
  
  # 简单的映射（假设顺序是M.C, C.C, M.M, C.M）
  if (length(current_levels) == 4) {
    emm_group@grid$P_F <- factor(new_levels, levels = new_levels)
    emm_group@levels$P_F <- new_levels
  }
  
  # 执行Tukey HSD检验
  tukey_result <- pairs(emm_group, adjust = "tukey")
  
  # 生成显著性字母
  cat(sprintf("\n=== %s: 生成显著性字母 ===\n", dataset_name))
  
  suppressMessages({
    cld_result <- cld(emm_group, alpha = 0.05, Letters = letters, adjust = "tukey")
  })
  
  # 提取组名和显著性字母
  sig_letters <- data.frame(
    group = cld_result$P_F,
    significance = cld_result$.group,
    stringsAsFactors = FALSE
  )
  
  # 如果没有显著差异，用"ns"表示
  if (length(unique(sig_letters$significance)) == 1) {
    sig_letters$significance <- "ns"
  } else {
    sig_letters$significance <- ifelse(sig_letters$significance == "", 
                                       "ns", 
                                       sig_letters$significance)
  }
  
  # 准备详细的方差分析结果
  cat(sprintf("\n=== %s: 准备方差分析结果 ===\n", dataset_name))
  
  anova_df_full <- as.data.frame(anova_result)
  
  # 提取需要的信息
  if ("Chisq" %in% colnames(anova_df_full)) {
    anova_df <- data.frame(
      Term = rownames(anova_df_full),
      Chisq = anova_df_full$Chisq,
      DF = anova_df_full$Df,
      p_value = anova_df_full$`Pr(>Chisq)`,
      stringsAsFactors = FALSE
    )
  } else {
    anova_df <- data.frame(
      Term = rownames(anova_df_full),
      stringsAsFactors = FALSE
    )
    
    # 尝试提取统计量和p值
    if ("F" %in% colnames(anova_df_full)) {
      anova_df$F_value <- anova_df_full$F
    }
    if ("Chisq" %in% colnames(anova_df_full)) {
      anova_df$Chisq <- anova_df_full$Chisq
    }
    if ("Pr(>F)" %in% colnames(anova_df_full)) {
      anova_df$p_value <- anova_df_full$`Pr(>F)`
    } else if ("Pr(>Chisq)" %in% colnames(anova_df_full)) {
      anova_df$p_value <- anova_df_full$`Pr(>Chisq)`
    }
  }
  
  # 移除截距项
  if ("(Intercept)" %in% anova_df$Term) {
    anova_df <- anova_df[anova_df$Term != "(Intercept)", ]
  }
  
  # 格式化数值
  for (col in colnames(anova_df)) {
    if (is.numeric(anova_df[[col]]) && col != "DF") {
      anova_df[[col]] <- format(round(anova_df[[col]], 3), nsmall = 3, scientific = FALSE)
    }
  }
  
  # 添加显著性标记
  if ("p_value" %in% colnames(anova_df)) {
    anova_df$Significance <- ifelse(
      as.numeric(anova_df$p_value) < 0.001, "***",
      ifelse(as.numeric(anova_df$p_value) < 0.01, "**",
             ifelse(as.numeric(anova_df$p_value) < 0.05, "*", ""))
    )
  }
  
  # 输出到控制台
  cat("方差分析结果:\n")
  print(anova_df)
  
  # 准备Tukey HSD结果
  cat(sprintf("\n=== %s: 准备Tukey HSD结果 ===\n", dataset_name))
  
  tukey_summary <- summary(tukey_result)
  
  # 创建数据框
  tukey_df <- data.frame(
    Comparison = tukey_summary$contrast,
    stringsAsFactors = FALSE
  )
  
  # 添加其他列
  if ("ratio" %in% colnames(tukey_summary)) {
    tukey_df$Estimate <- tukey_summary$ratio
  } else if ("estimate" %in% colnames(tukey_summary)) {
    tukey_df$Estimate <- tukey_summary$estimate
  }
  
  if ("SE" %in% colnames(tukey_summary)) {
    tukey_df$SE <- tukey_summary$SE
  }
  
  if ("df" %in% colnames(tukey_summary)) {
    tukey_df$DF <- tukey_summary$df
  } else {
    tukey_df$DF <- NA
  }
  
  if ("p.value" %in% colnames(tukey_summary)) {
    tukey_df$p_value <- tukey_summary$p.value
  }
  
  # 添加统计量列
  if ("z.ratio" %in% colnames(tukey_summary)) {
    tukey_df$Statistic <- tukey_summary$z.ratio
  } else if ("t.ratio" %in% colnames(tukey_summary)) {
    tukey_df$Statistic <- tukey_summary$t.ratio
  } else {
    tukey_df$Statistic <- NA
  }
  
  # 格式化数值列
  numeric_cols <- c("Estimate", "SE", "Statistic", "p_value")
  for (col in numeric_cols) {
    if (col %in% colnames(tukey_df) && !all(is.na(tukey_df[[col]]))) {
      tukey_df[[col]] <- as.numeric(tukey_df[[col]])
      
      if (col == "p_value") {
        tukey_df[[col]] <- format(tukey_df[[col]], digits = 3, scientific = FALSE)
      } else {
        tukey_df[[col]] <- format(round(tukey_df[[col]], 3), nsmall = 3, scientific = FALSE)
      }
    }
  }
  
  # 添加显著性标记
  tukey_df$Significance <- ifelse(
    as.numeric(tukey_df$p_value) < 0.001, "***",
    ifelse(as.numeric(tukey_df$p_value) < 0.01, "**",
           ifelse(as.numeric(tukey_df$p_value) < 0.05, "*",
                  ifelse(as.numeric(tukey_df$p_value) < 0.1, ".", "")))
  )
  
  # 输出到控制台
  cat("Tukey HSD结果:\n")
  print(tukey_df, row.names = FALSE)
  
  # 准备绘图数据
  data_plot <- df_long %>%
    group_by(group, time) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      se = sd(value, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = 'drop'
    )
  
  # 为繁殖和生存数据创建小时单位的时间变量
  if (dataset_name %in% c("繁殖(reproductive)", "生存(survival)")) {
    data_plot <- data_plot %>%
      mutate(time_hours = time * 12)  # 将天转换为小时（乘以12）
  }
  
  # 按照"MC, CC, MM, CM"的顺序设置group的因子水平
  group_order <- c("MC", "CC", "MM", "CM")
  data_plot$group <- factor(data_plot$group, levels = group_order)
  
  data_plot <- data_plot %>%
    left_join(sig_letters, by = "group") %>%
    arrange(group, time)
  
  # 创建identity列
  data_plot <- data_plot %>%
    mutate(identity = paste0(group, ",Significance:", significance))
  
  # 将identity转换为因子
  data_plot$identity <- factor(data_plot$identity, levels = unique(data_plot$identity))
  
  # 返回结果
  return(list(
    dataset_name = dataset_name,
    model_type = model_type,
    model = model,
    anova_df = anova_df,
    tukey_df = tukey_df,
    sig_letters = sig_letters,
    data_plot = data_plot,
    raw_data = data,
    df_long = df_long
  ))
}

# 5. 分别处理四个数据集（使用修复版函数）
results <- list()

# 密度数据
cat("\n\n=== 开始处理密度数据 ===\n")
results$density <- process_dataset_glmm_fixed(density_M, "密度(density)")

# 性比数据
cat("\n\n=== 开始处理性比数据 ===\n")
results$mixs <- process_dataset_glmm_fixed(mixs_M, "性比(mixs)")

# 繁殖数据
cat("\n\n=== 开始处理繁殖数据 ===\n")
results$reproductive <- process_dataset_glmm_fixed(reproductive_M, "繁殖(reproductive)")

# 生存数据
cat("\n\n=== 开始处理生存数据 ===\n")
results$survival <- process_dataset_glmm_fixed(survival_M, "生存(survival)")

# 6. 创建Excel文件，包含所有结果
cat("\n=== 创建Excel文件 ===\n")

excel_data <- list()

for (dataset_name in names(results)) {
  result <- results[[dataset_name]]
  
  # 添加ANOVA结果
  anova_sheet_name <- paste0(result$dataset_name, "_ANOVA")
  excel_data[[anova_sheet_name]] <- result$anova_df
  
  # 添加Tukey HSD结果
  tukey_sheet_name <- paste0(result$dataset_name, "_Tukey_HSD")
  excel_data[[tukey_sheet_name]] <- result$tukey_df
  
  # 添加显著性字母
  sig_sheet_name <- paste0(result$dataset_name, "_显著性字母")
  excel_data[[sig_sheet_name]] <- result$sig_letters
  
  # 添加模型信息
  model_info <- data.frame(
    数据集 = result$dataset_name,
    模型类型 = result$model_type,
    字符串AsFactors = FALSE
  )
  excel_data[[paste0(result$dataset_name, "_模型信息")]] <- model_info
  
  # 将显著性字母添加到原始数据
  data_with_sig <- result$raw_data %>%
    left_join(result$sig_letters, by = "group")
  
  # 保存带显著性字母的数据
  write.csv(data_with_sig, 
            sprintf("%s_with_significance.csv", dataset_name), 
            row.names = FALSE)
  cat(sprintf("带显著性字母的%s数据已保存为: %s_with_significance.csv\n", 
              result$dataset_name, dataset_name))
}

# 保存为Excel文件
write.xlsx(excel_data, file = "All_Statistical_Analysis_Results_GLMM_M.xlsx")
cat("所有统计分析结果Excel文件已保存为: All_Statistical_Analysis_Results_GLMM_M.xlsx\n")

# 7. 分别绘制四个数据集的图形
cat("\n=== 分别绘制四个数据集的图形 ===\n")

# 定义绘图函数 - 修改版，支持不同数据集的时间单位
create_plot_simple <- function(data_plot, dataset_name, y_label, filename) {
  
  # 根据数据集确定x轴变量和标签
  if (dataset_name %in% c("繁殖(reproductive)", "生存(survival)")) {
    # 繁殖和生存数据使用小时单位
    x_var <- "time_hours"
    x_label <- "Time(Hours)"
    # 创建小时单位的x轴刻度
    max_hours <- max(data_plot$time_hours)
    x_breaks <- seq(0, max_hours, by = 12)  # 每12小时一个刻度
    error_width <- 2  # 调整误差条宽度以适应小时单位
  } else {
    # 密度和性比数据使用天单位
    x_var <- "time"
    x_label <- "Time(Day)"
    max_day <- max(data_plot$time)
    x_breaks <- seq(0, max_day, by = 1)  # 每天一个刻度
    error_width <- 0.5
  }
  
  # 创建图形 - 根据数据集选择x轴变量
  p <- ggplot(data_plot, aes_string(x = x_var, y = "mean", 
                                    color = "identity", 
                                    shape = "identity",
                                    group = "identity")) +
    geom_line(size = 1.5) +
    geom_point(size = 3, aes_string(fill = "identity", shape = "identity"), 
               color = "black", stroke = 1) +
    geom_errorbar(aes_string(ymin = "mean - se", ymax = "mean + se"), 
                  color = "black", width = error_width) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 12),
          axis.line = element_line(),
          legend.background = element_blank(),
          legend.key.size = unit(25, "pt"),
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          legend.position = c(0.15, 0.75),
          plot.margin = unit(rep(1, 4), 'mm')) +
    # 使用适合M数据的颜色和形状映射
    scale_color_manual(values = c("#B2DF8A", "#33A02C","#E3B0A3" , "#D65D48")) +
    scale_fill_manual(values = c("#B2DF8A", "#33A02C","#E3B0A3" , "#D65D48")) +
    scale_shape_manual(values = c(21, 24, 21, 24)) +
    labs(x = x_label, y = y_label)
  
  # 设置x轴刻度
  p <- p + scale_x_continuous(breaks = x_breaks)
  
  # 根据数据集调整y轴
  if (dataset_name == "密度(density)") {
    p <- p + scale_y_continuous(breaks = seq(0, 6, 1), limits = c(0, 6))
  } else if (dataset_name == "性比(mixs)") {
    p <- p + scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1))
  } else if (dataset_name == "繁殖(reproductive)") {
    p <- p + scale_y_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1.25))
  } else if (dataset_name == "生存(survival)") {
    p <- p + scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1))
  }
  
  # 调整图例位置
  if (dataset_name == "密度(density)") {
    p <- p 
  } else if (dataset_name == "性比(mixs)") {
    p <- p + theme(legend.position = c(0.8, 0.75))
  } else if (dataset_name == "繁殖(reproductive)") {
    p <- p 
  } else if (dataset_name == "生存(survival)") {
    p <- p + theme(legend.position = c(0.8, 0.75))
  }
  
  # 保存图形
  ggsave(filename, p, width = 190, height = 100, units = "mm", dpi = 1000)
  cat(sprintf("%s图形已保存为: %s\n", dataset_name, filename))
  
  # 显示图形
  print(p)
  
  return(p)
}

# 绘制密度图形
cat("\n--- 绘制密度图形 ---\n")
p_density <- create_plot_simple(
  data_plot = results$density$data_plot,
  dataset_name = "密度(density)",
  y_label = "Population density(ind/mL)",
  filename = "density_M_GLMM.svg"
)

# 绘制性比图形
cat("\n--- 绘制性比图形 ---\n")
p_mixs <- create_plot_simple(
  data_plot = results$mixs$data_plot,
  dataset_name = "性比(mixs)",
  y_label = "Proportion of mictic females",
  filename = "mixs_M_GLMM.svg"
)

# 绘制繁殖图形
cat("\n--- 绘制繁殖图形 ---\n")
p_reproductive <- create_plot_simple(
  data_plot = results$reproductive$data_plot,
  dataset_name = "繁殖(reproductive)",
  y_label = "Reproductive rate",
  filename = "reproductive_M_GLMM.svg"
)

# 绘制生存图形
cat("\n--- 绘制生存图形 ---\n")
p_survival <- create_plot_simple(
  data_plot = results$survival$data_plot,
  dataset_name = "生存(survival)",
  y_label = "Survival rate",
  filename = "survival_M_GLMM.svg"
)

# 创建2x2的图形组合
combined_plot <- plot_grid(
  p_density, p_mixs, p_reproductive, p_survival,
  labels = c("a", "b", "c", "d"),
  ncol = 2,
  nrow = 2
)

# 保存组合图形
ggsave('population_M_2x2_GLMM.svg', combined_plot, width = 380, height = 200, units = "mm", dpi = 1000)
cat("2x2组合图形已保存为: population_M_2x2_GLMM.svg\n")

# 创建4x1的图形组合
combined_plot_vertical <- plot_grid(
  p_density, p_mixs, p_reproductive, p_survival,
  labels = c("a", "b", "c", "d"),
  ncol = 1,
  nrow = 4
)

ggsave('population_M_4x1_GLMM.svg', combined_plot_vertical, width = 190, height = 400, units = "mm", dpi = 1000)
cat("4x1组合图形已保存为: population_M_4x1_GLMM.svg\n")

# 8. 输出分析结果摘要
cat("\n=== 分析结果摘要 ===\n")

for (dataset_name in names(results)) {
  result <- results[[dataset_name]]
  cat(sprintf("1. %s:\n", result$dataset_name))
  cat(sprintf("   使用模型: %s\n", result$model_type))
  
  # 显示主要效应和交互作用
  if (nrow(result$anova_df) > 0) {
    cat("   主要统计结果:\n")
    for (i in 1:nrow(result$anova_df)) {
      term <- result$anova_df$Term[i]
      p_val <- as.numeric(result$anova_df$p_value[i])
      sig <- result$anova_df$Significance[i]
      
      if (term %in% c("P", "F", "P:F")) {
        cat(sprintf("     %s: p = %s %s\n", term, 
                    ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val)), 
                    sig))
      }
    }
  }
}

cat("\n2. 输出文件:\n")
cat("   - All_Statistical_Analysis_Results_GLMM_M.xlsx: 所有统计分析结果Excel文件\n")
cat("   - density_with_significance.csv: 带显著性字母的密度数据\n")
cat("   - mixs_with_significance.csv: 带显著性字母的性比数据\n")
cat("   - reproductive_with_significance.csv: 带显著性字母的繁殖数据\n")
cat("   - survival_with_significance.csv: 带显著性字母的生存数据\n")
cat("   - density_M_GLMM.svg: 密度曲线图\n")
cat("   - mixs_M_GLMM.svg: 性比曲线图\n")
cat("   - reproductive_M_GLMM.svg: 繁殖曲线图\n")
cat("   - survival_M_GLMM.svg: 生存曲线图\n")
cat("   - population_M_2x2_GLMM.svg: 2x2组合图形\n")
cat("   - population_M_4x1_GLMM.svg: 4x1组合图形\n")

cat("\n分析完成！\n")

# 9. 保存所有结果到RData文件
save(results, file = "M_population_analysis_results_GLMM.RData")
cat("所有分析结果已保存到: M_population_analysis_results_GLMM.RData\n")
