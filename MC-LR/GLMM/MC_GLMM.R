# 1. 清理环境并设置工作目录
rm(list=ls())

# 设置工作目录，请根据实际情况修改
setwd("E:/FILE/Masterexp/Paper/1_mother_choice/1_data/GLMM/MC/")

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
density_L <- read.csv("density_L.csv", header = TRUE, stringsAsFactors = FALSE)
mixs_L <- read.csv("mixsrate_L.csv", header = TRUE, stringsAsFactors = FALSE)
reproductive_L <- read.csv("reproductive_L.csv", header = TRUE, stringsAsFactors = FALSE)
survival_L <- read.csv("survival_L.csv", header = TRUE, stringsAsFactors = FALSE)

# 检查数据结构
cat("检查数据结构:\n")
cat("density_L 维度:", dim(density_L), "\n")
cat("mixs_L 维度:", dim(mixs_L), "\n")
cat("reproductive_L 维度:", dim(reproductive_L), "\n")
cat("survival_L 维度:", dim(survival_L), "\n")

# 4. 定义通用函数用于数据预处理和分析（全部使用GLMM）
process_dataset_glmm <- function(data, dataset_name) {
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
      group = factor(group, levels = c("LC", "CC", "CL", "LL")),
      P = factor(P, levels = c("C", "L")),
      F = factor(F, levels = c("C", "L")),
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
  
  # 根据数据类型选择GLMM模型（跳过LMM）
  if (is_count_data) {
    cat("选择负二项GLMM模型处理计数数据\n")
    
    # 尝试负二项分布
    tryCatch({
      model <- glmmTMB(value ~ P * F + (1|id),
                       data = df_long,
                       family = nbinom2())
      model_type <- "负二项GLMM"
      cat("负二项GLMM拟合成功!\n")
    }, error = function(e) {
      # 如果负二项失败，尝试泊松分布
      cat("负二项GLMM拟合失败，尝试泊松GLMM\n")
      tryCatch({
        model <- glmmTMB(value ~ P * F + (1|id),
                         data = df_long,
                         family = poisson())
        model_type <- "泊松GLMM"
        cat("泊松GLMM拟合成功!\n")
      }, error = function(e2) {
        # 如果泊松也失败，回退到高斯GLMM
        cat("泊松GLMM拟合失败，回退到高斯GLMM\n")
        model <- glmmTMB(value ~ P * F + (1|id),
                         data = df_long,
                         family = gaussian())
        model_type <- "高斯GLMM"
        cat("高斯GLMM拟合成功!\n")
      })
    })
  } else if (is_proportion_data) {
    cat("选择Beta GLMM模型处理比例数据\n")
    
    # 对于Beta分布，值必须在(0,1)之间，不能包含0或1
    # 调整边界值
    df_long <- df_long %>%
      mutate(value_adj = case_when(
        value == 0 ~ 0.001,
        value == 1 ~ 0.999,
        TRUE ~ value
      ))
    
    tryCatch({
      model <- glmmTMB(value_adj ~ P * F + (1|id),
                       data = df_long,
                       family = beta_family())
      model_type <- "Beta GLMM"
      cat("Beta GLMM拟合成功!\n")
    }, error = function(e) {
      # 如果Beta分布失败，尝试高斯分布+logit链接
      cat("Beta GLMM拟合失败，尝试高斯分布GLMM(logit链接)\n")
      tryCatch({
        model <- glmmTMB(value ~ P * F + (1|id),
                         data = df_long,
                         family = gaussian(link = "logit"))
        model_type <- "高斯GLMM(logit链接)"
        cat("高斯GLMM(logit链接)拟合成功!\n")
      }, error = function(e2) {
        # 如果高斯也失败，回退到高斯GLMM（恒等链接）
        cat("高斯GLMM(logit链接)拟合失败，回退到高斯GLMM\n")
        model <- glmmTMB(value ~ P * F + (1|id),
                         data = df_long,
                         family = gaussian())
        model_type <- "高斯GLMM"
        cat("高斯GLMM拟合成功!\n")
      })
    })
  } else {
    # 连续数据，直接尝试GLMM（跳过LMM）
    cat("选择GLMM模型处理连续数据\n")
    
    # 先尝试Gamma GLMM（对于正值数据）
    if (min(df_long$value, na.rm = TRUE) > 0) {
      cat("数据均为正值，尝试Gamma GLMM\n")
      tryCatch({
        model <- glmmTMB(value ~ P * F + (1|id),
                         data = df_long,
                         family = Gamma(link = "log"))
        model_type <- "Gamma GLMM"
        cat("Gamma GLMM拟合成功!\n")
      }, error = function(e) {
        # Gamma失败，尝试高斯GLMM
        cat("Gamma GLMM拟合失败，尝试高斯GLMM\n")
        tryCatch({
          model <- glmmTMB(value ~ P * F + (1|id),
                           data = df_long,
                           family = gaussian())
          model_type <- "高斯GLMM"
          cat("高斯GLMM拟合成功!\n")
        }, error = function(e2) {
          # 所有尝试失败，使用简化模型
          cat("所有GLMM尝试失败，使用简化模型\n")
          model <- glmmTMB(value ~ group + (1|id), 
                           data = df_long, 
                           family = gaussian())
          model_type <- "简化GLMM"
          cat("简化GLMM拟合成功!\n")
        })
      })
    } else {
      # 数据包含0或负值，使用高斯GLMM
      cat("数据包含非正值，使用高斯GLMM\n")
      tryCatch({
        model <- glmmTMB(value ~ P * F + (1|id),
                         data = df_long,
                         family = gaussian())
        model_type <- "高斯GLMM"
        cat("高斯GLMM拟合成功!\n")
      }, error = function(e) {
        cat("高斯GLMM拟合失败，使用简化模型\n")
        model <- glmmTMB(value ~ group + (1|id), 
                         data = df_long, 
                         family = gaussian())
        model_type <- "简化GLMM"
        cat("简化GLMM拟合成功!\n")
      })
    }
  }
  
  # 双因素方差分析
  cat(sprintf("\n=== %s: 双因素方差分析结果 ===\n", dataset_name))
  
  anova_result <- Anova(model, type = "III")
  
  # 多重比较
  cat(sprintf("\n=== %s: 多重比较(Tukey HSD) ===\n", dataset_name))
  
  # 创建组合变量进行多重比较
  df_long <- df_long %>%
    mutate(P_F = interaction(P, F, sep = ""))
  
  # 根据模型类型创建用于多重比较的模型
  if (model_type == "简化GLMM") {
    model_group <- model
    emm_group <- emmeans(model_group, ~ group, type = "response")
    group_var <- "group"
  } else if (model_type == "Gamma GLMM") {
    model_group <- glmmTMB(value ~ P_F + (1|id), 
                           data = df_long, 
                           family = Gamma(link = "log"))
    emm_group <- emmeans(model_group, ~ P_F, type = "response")
    group_var <- "P_F"
  } else if (model_type == "高斯GLMM") {
    model_group <- glmmTMB(value ~ P_F + (1|id), 
                           data = df_long, 
                           family = gaussian())
    emm_group <- emmeans(model_group, ~ P_F, type = "response")
    group_var <- "P_F"
  } else if (model_type == "负二项GLMM") {
    model_group <- glmmTMB(value ~ P_F + (1|id), 
                           data = df_long, 
                           family = nbinom2())
    emm_group <- emmeans(model_group, ~ P_F, type = "response")
    group_var <- "P_F"
  } else if (model_type == "泊松GLMM") {
    model_group <- glmmTMB(value ~ P_F + (1|id), 
                           data = df_long, 
                           family = poisson())
    emm_group <- emmeans(model_group, ~ P_F, type = "response")
    group_var <- "P_F"
  } else if (model_type == "Beta GLMM") {
    # 对于Beta GLMM，使用调整后的值
    if ("value_adj" %in% names(df_long)) {
      model_group <- glmmTMB(value_adj ~ P_F + (1|id), 
                             data = df_long, 
                             family = beta_family())
    } else {
      model_group <- glmmTMB(value ~ P_F + (1|id), 
                             data = df_long, 
                             family = gaussian(link = "logit"))
    }
    emm_group <- emmeans(model_group, ~ P_F, type = "response")
    group_var <- "P_F"
  } else if (model_type == "高斯GLMM(logit链接)") {
    model_group <- glmmTMB(value ~ P_F + (1|id), 
                           data = df_long, 
                           family = gaussian(link = "logit"))
    emm_group <- emmeans(model_group, ~ P_F, type = "response")
    group_var <- "P_F"
  } else {
    # 默认高斯GLMM
    model_group <- glmmTMB(value ~ P_F + (1|id), 
                           data = df_long, 
                           family = gaussian())
    emm_group <- emmeans(model_group, ~ P_F, type = "response")
    group_var <- "P_F"
  }
  
  # 重命名levels为LC, CC, CL, LL
  if (group_var == "P_F") {
    current_levels <- levels(emm_group@grid$P_F)
    # 创建正确的映射（基于C.L, C.C, L.C, L.L到LC, CC, CL, LL）
    # 假设标准顺序是C.L, C.C, L.C, L.L
    if (all(c("C.L", "C.C", "L.C", "L.L") %in% current_levels)) {
      level_mapping <- c("L.C" = "LC", "C.C" = "CC", "C.L" = "CL", "L.L" = "LL")
      emm_group@grid$P_F <- factor(level_mapping[as.character(emm_group@grid$P_F)], 
                                   levels = c("LC", "CC", "CL", "LL"))
      emm_group@levels$P_F <- c("LC", "CC", "CL", "LL")
    } else {
      # 如果顺序不同，尝试直接重命名
      levels(emm_group@grid$P_F) <- c("LC", "CC", "CL", "LL")
      emm_group@levels$P_F <- c("LC", "CC", "CL", "LL")
    }
  } else {
    # 对于group变量，确保顺序正确
    emm_group@grid$group <- factor(emm_group@grid$group, levels = c("LC", "CC", "CL", "LL"))
    emm_group@levels$group <- c("LC", "CC", "CL", "LL")
  }
  
  # 执行Tukey HSD检验
  tukey_result <- pairs(emm_group, adjust = "tukey")
  
  # 生成显著性字母
  cat(sprintf("\n=== %s: 生成显著性字母 ===\n", dataset_name))
  
  suppressMessages({
    cld_result <- cld(emm_group, alpha = 0.05, Letters = letters, adjust = "tukey")
  })
  
  # 提取组名和显著性字母
  if (group_var == "P_F") {
    sig_letters <- data.frame(
      group = cld_result$P_F,
      significance = cld_result$.group
    )
  } else {
    sig_letters <- data.frame(
      group = cld_result$group,
      significance = cld_result$.group
    )
  }
  
  # 如果没有显著差异，用"ns"表示
  # 检查是否所有组都有相同的显著性字母
  if (length(unique(sig_letters$significance)) == 1) {
    # 如果所有组都有相同的字母，说明没有显著差异
    sig_letters$significance <- "ns"
  } else {
    # 如果有不同的字母，保留原字母，但将空字符串替换为"ns"
    sig_letters$significance <- ifelse(sig_letters$significance == "", 
                                       "ns", 
                                       sig_letters$significance)
  }
  
  # 准备详细的方差分析结果
  cat(sprintf("\n=== %s: 准备方差分析结果 ===\n", dataset_name))
  
  # 对于GLMM，使用car::Anova的结果
  anova_df_full <- as.data.frame(anova_result)
  
  # 检查列名，提取适当的信息
  colnames_df <- colnames(anova_df_full)
  
  if ("Chisq" %in% colnames_df && "Pr(>Chisq)" %in% colnames_df) {
    # Type III ANOVA结果
    anova_df <- data.frame(
      Term = rownames(anova_df_full),
      Chisq = anova_df_full$Chisq,
      DF = anova_df_full$Df,
      p_value = anova_df_full$`Pr(>Chisq)`,
      stringsAsFactors = FALSE
    )
  } else if ("LR Chisq" %in% colnames_df) {
    # 似然比检验结果
    anova_df <- data.frame(
      Term = rownames(anova_df_full),
      LR_Chisq = anova_df_full$`LR Chisq`,
      DF = anova_df_full$Df,
      p_value = anova_df_full$`Pr(>Chisq)`,
      stringsAsFactors = FALSE
    )
  } else {
    # 其他格式
    anova_df <- data.frame(
      Term = rownames(anova_df_full),
      stringsAsFactors = FALSE
    )
    
    # 尝试提取统计量和p值
    if ("F" %in% colnames_df) {
      anova_df$F_value <- anova_df_full$F
    }
    if ("Chisq" %in% colnames_df) {
      anova_df$Chisq <- anova_df_full$Chisq
    }
    if ("Pr(>F)" %in% colnames_df) {
      anova_df$p_value <- anova_df_full$`Pr(>F)`
    } else if ("Pr(>Chisq)" %in% colnames_df) {
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
  
  # 将Tukey HSD结果转换为数据框
  tukey_summary <- summary(tukey_result)
  
  # 创建数据框
  tukey_df <- data.frame(
    Comparison = tukey_summary$contrast,
    stringsAsFactors = FALSE
  )
  
  # 添加其他列，但检查它们是否存在
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
  
  # 添加统计量列（z.ratio或t.ratio）
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
      
      # 格式化为三位小数
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
  
  # 按照"LC, CC, CL, LL"的顺序设置group的因子水平
  group_order <- c("LC", "CC", "CL", "LL")
  data_plot$group <- factor(data_plot$group, levels = group_order)
  
  data_plot <- data_plot %>%
    left_join(sig_letters, by = "group") %>%
    arrange(group, time)  # 按照group的顺序排列
  
  # 创建identity列，此时已经按照正确的顺序
  data_plot <- data_plot %>%
    mutate(identity = paste0(group, ",Significance:", significance))
  
  # 将identity转换为因子，保持与group相同的顺序
  data_plot$identity <- factor(data_plot$identity, levels = unique(data_plot$identity))
  
  # 返回结果
  return(list(
    dataset_name = dataset_name,
    model_type = model_type,
    anova_df = anova_df,
    tukey_df = tukey_df,
    sig_letters = sig_letters,
    data_plot = data_plot,
    raw_data = data
  ))
}

# 5. 分别处理四个数据集
results <- list()

# 密度数据
results$density <- process_dataset_glmm(density_L, "密度(density)")

# 性比数据
results$mixs <- process_dataset_glmm(mixs_L, "性比(mixs)")

# 繁殖数据
results$reproductive <- process_dataset_glmm(reproductive_L, "繁殖(reproductive)")

# 生存数据
results$survival <- process_dataset_glmm(survival_L, "生存(survival)")

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
write.xlsx(excel_data, file = "All_Statistical_Analysis_Results_L_GLMM.xlsx")
cat("所有统计分析结果Excel文件已保存为: All_Statistical_Analysis_Results_L_GLMM.xlsx\n")

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
    # 使用适合L数据的颜色和形状映射
    scale_color_manual(values = c("#B2E2DC", "#33A02C", "#E8CF92", "#D6A228")) +
    scale_fill_manual(values = c("#B2E2DC", "#33A02C", "#E8CF92", "#D6A228")) +
    scale_shape_manual(values = c(21, 24, 21, 24)) +
    labs(x = x_label, y = y_label)
  
  # 设置x轴刻度
  p <- p + scale_x_continuous(breaks = x_breaks)
  
  # 根据数据集调整y轴
  if (dataset_name == "密度(density)") {
    p <- p + scale_y_continuous(breaks = seq(0, 8, 1), limits = c(0, 8))
  } else if (dataset_name == "性比(mixs)") {
    p <- p + scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1))
  } else if (dataset_name == "繁殖(reproductive)") {
    p <- p + scale_y_continuous(breaks = seq(0, 2, 0.5), limits = c(0, 2))
  } else if (dataset_name == "生存(survival)") {
    p <- p + scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1))
  }
  
  # 调整图例位置
  if (dataset_name == "密度(density)") {
    p <- p 
  } else if (dataset_name == "性比(mixs)") {
    p <- p + theme(legend.position = c(0.15, 0.75))
  } else if (dataset_name == "繁殖(reproductive)") {
    p <- p 
  } else if (dataset_name == "生存(survival)") {
    p <- p + theme(legend.position = c(0.85, 0.75))
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
  filename = "density_L_GLMM.svg"
)

# 绘制性比图形
cat("\n--- 绘制性比图形 ---\n")
p_mixs <- create_plot_simple(
  data_plot = results$mixs$data_plot,
  dataset_name = "性比(mixs)",
  y_label = "Proportion of mictic females",
  filename = "mixs_L_GLMM.svg"
)

# 绘制繁殖图形
cat("\n--- 绘制繁殖图形 ---\n")
p_reproductive <- create_plot_simple(
  data_plot = results$reproductive$data_plot,
  dataset_name = "繁殖(reproductive)",
  y_label = "Reproductive rate",
  filename = "reproductive_L_GLMM.svg"
)

# 绘制生存图形
cat("\n--- 绘制生存图形 ---\n")
p_survival <- create_plot_simple(
  data_plot = results$survival$data_plot,
  dataset_name = "生存(survival)",
  y_label = "Survival rate",
  filename = "survival_L_GLMM.svg"
)

# 创建2x2的图形组合
combined_plot <- plot_grid(
  p_density, p_mixs, p_reproductive, p_survival,
  labels = c("a", "b", "c", "d"),
  ncol = 2,
  nrow = 2
)

# 保存组合图形
ggsave('population_L_2x2_GLMM.svg', combined_plot, width = 380, height = 200, units = "mm", dpi = 1000)
cat("2x2组合图形已保存为: population_L_2x2_GLMM.svg\n")

# 创建4x1的图形组合
combined_plot_vertical <- plot_grid(
  p_density, p_mixs, p_reproductive, p_survival,
  labels = c("a", "b", "c", "d"),
  ncol = 1,
  nrow = 4
)

ggsave('population_L_4x1_GLMM.svg', combined_plot_vertical, width = 190, height = 400, units = "mm", dpi = 1000)
cat("4x1组合图形已保存为: population_L_4x1_GLMM.svg\n")

# 8. 输出分析结果摘要
cat("\n=== 分析结果摘要 ===\n")

for (dataset_name in names(results)) {
  result <- results[[dataset_name]]
  cat(sprintf("1. %s:\n", result$dataset_name))
  cat(sprintf("   使用模型: %s\n", result$model_type))
}

cat("\n2. 输出文件:\n")
cat("   - All_Statistical_Analysis_Results_L_GLMM.xlsx: 所有统计分析结果Excel文件\n")
cat("   - density_with_significance.csv: 带显著性字母的密度数据\n")
cat("   - mixs_with_significance.csv: 带显著性字母的性比数据\n")
cat("   - reproductive_with_significance.csv: 带显著性字母的繁殖数据\n")
cat("   - survival_with_significance.csv: 带显著性字母的生存数据\n")
cat("   - density_L_GLMM.svg: 密度曲线图\n")
cat("   - mixs_L_GLMM.svg: 性比曲线图\n")
cat("   - reproductive_L_GLMM.svg: 繁殖曲线图\n")
cat("   - survival_L_GLMM.svg: 生存曲线图\n")
cat("   - population_L_2x2_GLMM.svg: 2x2组合图形\n")
cat("   - population_L_4x1_GLMM.svg: 4x1组合图形\n")

cat("\n分析完成！\n")

# 9. 保存所有结果到RData文件
save(results, file = "L_population_analysis_results_GLMM.RData")

cat("所有分析结果已保存到: L_population_analysis_results_GLMM.RData\n")
