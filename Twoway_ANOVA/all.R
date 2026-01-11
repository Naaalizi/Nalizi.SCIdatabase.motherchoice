# 清理环境
rm(list = ls())

# 加载必要的包
library(carData)
library(car)
library(rcompanion)
library(agricolae)
library(DescTools)
library(dplyr)

# 设置工作目录
setwd("D:/SSData/variance_analysis/twoway_ANOVA/")

# 定义提取Scheirer-Ray-Hare结果的函数
extract_srh_df <- function(formula, data, param_name, treatment) {
  srh <- scheirerRayHare(formula, data = data)
  
  df_out <- data.frame(
    Parameter = param_name,
    Treatment = treatment,
    Effect = rownames(srh),
    Df = srh$Df,
    Sum_Sq = round(srh$`Sum Sq`, 3),
    Statistic = round(srh$H, 3),
    p_value = round(srh$`p.value`, 3),
    Test_Type = "Scheirer-Ray-Hare",
    stringsAsFactors = FALSE
  )
  
  rownames(df_out) <- NULL
  return(df_out)
}

# 定义提取ANOVA结果的函数
extract_anova_df <- function(aov_model, param_name, treatment) {
  anova_summary <- summary(aov_model)
  anova_df <- as.data.frame(anova_summary[[1]])
  
  df_out <- data.frame(
    Parameter = param_name,
    Treatment = treatment,
    Effect = rownames(anova_df),
    Df = anova_df$Df,
    Sum_Sq = round(anova_df$`Sum Sq`, 3),
    Mean_Sq = round(anova_df$`Mean Sq`, 3),
    Statistic = round(anova_df$`F value`, 3),
    p_value = round(anova_df$`Pr(>F)`, 3),
    Test_Type = "ANOVA",
    stringsAsFactors = FALSE
  )
  
  rownames(df_out) <- NULL
  return(df_out)
}

# 定义主要分析函数，自动选择检验方法
analyze_dataset <- function(file_path, response_var, param_name, treatment) {
  # 读取数据
  data <- read.csv(file_path)
  
  # 确保因子变量是因子类型
  data$P <- as.factor(data$P)
  data$F <- as.factor(data$F)
  
  # 构建公式
  formula_str <- paste(response_var, "~ P * F")
  formula_obj <- as.formula(formula_str)
  
  # 执行正态性检验（检验残差）
  model_lm <- lm(formula_obj, data = data)
  shapiro_result <- shapiro.test(residuals(model_lm))
  
  # 执行方差齐性检验
  levene_result <- leveneTest(formula_obj, data = data)
  
  cat("\n==========================================\n")
  cat("分析数据集:", file_path, "\n")
  cat("参数:", param_name, "\n")
  cat("处理:", treatment, "\n")
  cat("Shapiro-Wilk正态性检验p值:", round(shapiro_result$p.value, 4), "\n")
  cat("Levene方差齐性检验p值:", round(levene_result$`Pr(>F)`[1], 4), "\n")
  
  # 判断条件：p值>0.05表示满足假设
  normality_ok <- shapiro_result$p.value > 0.05
  homogeneity_ok <- levene_result$`Pr(>F)`[1] > 0.05
  
  cat("满足正态性假设:", ifelse(normality_ok, "是", "否"), "\n")
  cat("满足方差齐性假设:", ifelse(homogeneity_ok, "是", "否"), "\n")
  
  # 自动选择检验方法
  if (normality_ok && homogeneity_ok) {
    cat("使用双因素ANOVA\n")
    aov_model <- aov(formula_obj, data = data)
    result_df <- extract_anova_df(aov_model, param_name, treatment)
    
    # 进行事后检验
    tukey_result <- TukeyHSD(aov_model)
    posthoc <- list(Tukey = tukey_result)
    
  } else {
    cat("使用Scheirer-Ray-Hare检验\n")
    result_df <- extract_srh_df(formula_obj, data, param_name, treatment)
    
    # 进行事后检验（使用Wilcoxon秩和检验）
    data$wilcox.group <- paste0(data$P, data$F)
    posthoc <- pairwise.wilcox.test(x = data[[response_var]], 
                                    g = data$wilcox.group,
                                    p.adjust.method = "BH")
  }
  
  cat("检验完成\n")
  
  # 返回结果
  return(list(
    result_df = result_df,
    posthoc = posthoc,
    shapiro_p = shapiro_result$p.value,
    levene_p = levene_result$`Pr(>F)`[1],
    test_used = ifelse(normality_ok && homogeneity_ok, "ANOVA", "Scheirer-Ray-Hare"),
    normality_ok = normality_ok,
    homogeneity_ok = homogeneity_ok
  ))
}

# 初始化结果存储列表
all_results <- list()
posthoc_results <- list()
assumption_results <- list()

#####################################################
# Alga部分
#####################################################

# 1. growthr_M
growthr_M_result <- analyze_dataset("growthrate_M.csv", "rate", "Growth Rate", "Alga")
all_results[[1]] <- growthr_M_result$result_df
posthoc_results[["growthr_M"]] <- growthr_M_result$posthoc
assumption_results[["growthr_M"]] <- data.frame(
  Dataset = "growthr_M",
  Parameter = "Growth Rate",
  Treatment = "Alga",
  Shapiro_p = growthr_M_result$shapiro_p,
  Levene_p = growthr_M_result$levene_p,
  Normality_ok = growthr_M_result$normality_ok,
  Homogeneity_ok = growthr_M_result$homogeneity_ok,
  Test_used = growthr_M_result$test_used,
  stringsAsFactors = FALSE
)

# 2. rm_M
rm_M_result <- analyze_dataset("rm_M.csv", "rate", "r", "Alga")
all_results[[2]] <- rm_M_result$result_df
posthoc_results[["rm_M"]] <- rm_M_result$posthoc
assumption_results[["rm_M"]] <- data.frame(
  Dataset = "rm_M",
  Parameter = "r",
  Treatment = "Alga",
  Shapiro_p = rm_M_result$shapiro_p,
  Levene_p = rm_M_result$levene_p,
  Normality_ok = rm_M_result$normality_ok,
  Homogeneity_ok = rm_M_result$homogeneity_ok,
  Test_used = rm_M_result$test_used,
  stringsAsFactors = FALSE
)

# 3. G0_M
G0_M_result <- analyze_dataset("G0_M.csv", "rate", "G0", "Alga")
all_results[[3]] <- G0_M_result$result_df
posthoc_results[["G0_M"]] <- G0_M_result$posthoc
assumption_results[["G0_M"]] <- data.frame(
  Dataset = "G0_M",
  Parameter = "G0",
  Treatment = "Alga",
  Shapiro_p = G0_M_result$shapiro_p,
  Levene_p = G0_M_result$levene_p,
  Normality_ok = G0_M_result$normality_ok,
  Homogeneity_ok = G0_M_result$homogeneity_ok,
  Test_used = G0_M_result$test_used,
  stringsAsFactors = FALSE
)

# 4. R0_M
R0_M_result <- analyze_dataset("R0_M.csv", "rate", "R0", "Alga")
all_results[[4]] <- R0_M_result$result_df
posthoc_results[["R0_M"]] <- R0_M_result$posthoc
assumption_results[["R0_M"]] <- data.frame(
  Dataset = "R0_M",
  Parameter = "R0",
  Treatment = "Alga",
  Shapiro_p = R0_M_result$shapiro_p,
  Levene_p = R0_M_result$levene_p,
  Normality_ok = R0_M_result$normality_ok,
  Homogeneity_ok = R0_M_result$homogeneity_ok,
  Test_used = R0_M_result$test_used,
  stringsAsFactors = FALSE
)

# 5. T_M
T_M_result <- analyze_dataset("T_M.csv", "rate", "T", "Alga")
all_results[[5]] <- T_M_result$result_df
posthoc_results[["T_M"]] <- T_M_result$posthoc
assumption_results[["T_M"]] <- data.frame(
  Dataset = "T_M",
  Parameter = "T",
  Treatment = "Alga",
  Shapiro_p = T_M_result$shapiro_p,
  Levene_p = T_M_result$levene_p,
  Normality_ok = T_M_result$normality_ok,
  Homogeneity_ok = T_M_result$homogeneity_ok,
  Test_used = T_M_result$test_used,
  stringsAsFactors = FALSE
)

# 6. ex_M
ex_M_result <- analyze_dataset("ex_M.csv", "rate", "ex", "Alga")
all_results[[6]] <- ex_M_result$result_df
posthoc_results[["ex_M"]] <- ex_M_result$posthoc
assumption_results[["ex_M"]] <- data.frame(
  Dataset = "ex_M",
  Parameter = "ex",
  Treatment = "Alga",
  Shapiro_p = ex_M_result$shapiro_p,
  Levene_p = ex_M_result$levene_p,
  Normality_ok = ex_M_result$normality_ok,
  Homogeneity_ok = ex_M_result$homogeneity_ok,
  Test_used = ex_M_result$test_used,
  stringsAsFactors = FALSE
)

# 7. offspring
off_M_result <- analyze_dataset("offspring_M.csv", "offspring", "Offspring", "Alga")
all_results[[7]] <- off_M_result$result_df
posthoc_results[["off_M"]] <- off_M_result$posthoc
assumption_results[["off_M"]] <- data.frame(
  Dataset = "off_M",
  Parameter = "Offspring",
  Treatment = "Alga",
  Shapiro_p = off_M_result$shapiro_p,
  Levene_p = off_M_result$levene_p,
  Normality_ok = off_M_result$normality_ok,
  Homogeneity_ok = off_M_result$homogeneity_ok,
  Test_used = off_M_result$test_used,
  stringsAsFactors = FALSE
)

# 8. lifespan
lifespan_M_result <- analyze_dataset("lifespan_M.csv", "time", "Lifespan", "Alga")
all_results[[8]] <- lifespan_M_result$result_df
posthoc_results[["lifespan_M"]] <- lifespan_M_result$posthoc
assumption_results[["lifespan_M"]] <- data.frame(
  Dataset = "lifespan_M",
  Parameter = "Lifespan",
  Treatment = "Alga",
  Shapiro_p = lifespan_M_result$shapiro_p,
  Levene_p = lifespan_M_result$levene_p,
  Normality_ok = lifespan_M_result$normality_ok,
  Homogeneity_ok = lifespan_M_result$homogeneity_ok,
  Test_used = lifespan_M_result$test_used,
  stringsAsFactors = FALSE
)

# 9. body size
body_M_result <- analyze_dataset("body_M.csv", "V", "Body Size", "Alga")
all_results[[9]] <- body_M_result$result_df
posthoc_results[["body_M"]] <- body_M_result$posthoc
assumption_results[["body_M"]] <- data.frame(
  Dataset = "body_M",
  Parameter = "Body Size",
  Treatment = "Alga",
  Shapiro_p = body_M_result$shapiro_p,
  Levene_p = body_M_result$levene_p,
  Normality_ok = body_M_result$normality_ok,
  Homogeneity_ok = body_M_result$homogeneity_ok,
  Test_used = body_M_result$test_used,
  stringsAsFactors = FALSE
)

# 10. spine
sp_M_result <- analyze_dataset("spine_M.csv", "spine", "Spine", "Alga")
all_results[[10]] <- sp_M_result$result_df
posthoc_results[["sp_M"]] <- sp_M_result$posthoc
assumption_results[["sp_M"]] <- data.frame(
  Dataset = "sp_M",
  Parameter = "Spine",
  Treatment = "Alga",
  Shapiro_p = sp_M_result$shapiro_p,
  Levene_p = sp_M_result$levene_p,
  Normality_ok = sp_M_result$normality_ok,
  Homogeneity_ok = sp_M_result$homogeneity_ok,
  Test_used = sp_M_result$test_used,
  stringsAsFactors = FALSE
)

############################################################################
# MC-LR部分
############################################################################

# 1. growth_L
growth_L_result <- analyze_dataset("growthrate_MC.csv", "rate", "Growth Rate", "MC-LR")
all_results[[11]] <- growth_L_result$result_df
posthoc_results[["growth_L"]] <- growth_L_result$posthoc
assumption_results[["growth_L"]] <- data.frame(
  Dataset = "growth_L",
  Parameter = "Growth Rate",
  Treatment = "MC-LR",
  Shapiro_p = growth_L_result$shapiro_p,
  Levene_p = growth_L_result$levene_p,
  Normality_ok = growth_L_result$normality_ok,
  Homogeneity_ok = growth_L_result$homogeneity_ok,
  Test_used = growth_L_result$test_used,
  stringsAsFactors = FALSE
)

# 2. ir_L
ir_L_result <- analyze_dataset("r_MC.csv", "rate", "r", "MC-LR")
all_results[[12]] <- ir_L_result$result_df
posthoc_results[["ir_L"]] <- ir_L_result$posthoc
assumption_results[["ir_L"]] <- data.frame(
  Dataset = "ir_L",
  Parameter = "r",
  Treatment = "MC-LR",
  Shapiro_p = ir_L_result$shapiro_p,
  Levene_p = ir_L_result$levene_p,
  Normality_ok = ir_L_result$normality_ok,
  Homogeneity_ok = ir_L_result$homogeneity_ok,
  Test_used = ir_L_result$test_used,
  stringsAsFactors = FALSE
)

# 3. R_L
R_L_result <- analyze_dataset("R0_MC.csv", "rate", "R0", "MC-LR")
all_results[[13]] <- R_L_result$result_df
posthoc_results[["R_L"]] <- R_L_result$posthoc
assumption_results[["R_L"]] <- data.frame(
  Dataset = "R_L",
  Parameter = "R0",
  Treatment = "MC-LR",
  Shapiro_p = R_L_result$shapiro_p,
  Levene_p = R_L_result$levene_p,
  Normality_ok = R_L_result$normality_ok,
  Homogeneity_ok = R_L_result$homogeneity_ok,
  Test_used = R_L_result$test_used,
  stringsAsFactors = FALSE
)

# 4. T_L
T_L_result <- analyze_dataset("T_MC.csv", "rate", "T", "MC-LR")
all_results[[14]] <- T_L_result$result_df
posthoc_results[["T_L"]] <- T_L_result$posthoc
assumption_results[["T_L"]] <- data.frame(
  Dataset = "T_L",
  Parameter = "T",
  Treatment = "MC-LR",
  Shapiro_p = T_L_result$shapiro_p,
  Levene_p = T_L_result$levene_p,
  Normality_ok = T_L_result$normality_ok,
  Homogeneity_ok = T_L_result$homogeneity_ok,
  Test_used = T_L_result$test_used,
  stringsAsFactors = FALSE
)

# 5. G_L
G_L_result <- analyze_dataset("G_MC.csv", "rate", "G", "MC-LR")
all_results[[15]] <- G_L_result$result_df
posthoc_results[["G_L"]] <- G_L_result$posthoc
assumption_results[["G_L"]] <- data.frame(
  Dataset = "G_L",
  Parameter = "G",
  Treatment = "MC-LR",
  Shapiro_p = G_L_result$shapiro_p,
  Levene_p = G_L_result$levene_p,
  Normality_ok = G_L_result$normality_ok,
  Homogeneity_ok = G_L_result$homogeneity_ok,
  Test_used = G_L_result$test_used,
  stringsAsFactors = FALSE
)

# 6. ex_L
ex_L_result <- analyze_dataset("ex_MC.csv", "rate", "ex", "MC-LR")
all_results[[16]] <- ex_L_result$result_df
posthoc_results[["ex_L"]] <- ex_L_result$posthoc
assumption_results[["ex_L"]] <- data.frame(
  Dataset = "ex_L",
  Parameter = "ex",
  Treatment = "MC-LR",
  Shapiro_p = ex_L_result$shapiro_p,
  Levene_p = ex_L_result$levene_p,
  Normality_ok = ex_L_result$normality_ok,
  Homogeneity_ok = ex_L_result$homogeneity_ok,
  Test_used = ex_L_result$test_used,
  stringsAsFactors = FALSE
)

# 7. spine_L
spine_L_result <- analyze_dataset("spine_MC.csv", "spine", "Spine", "MC-LR")
all_results[[17]] <- spine_L_result$result_df
posthoc_results[["spine_L"]] <- spine_L_result$posthoc
assumption_results[["spine_L"]] <- data.frame(
  Dataset = "spine_L",
  Parameter = "Spine",
  Treatment = "MC-LR",
  Shapiro_p = spine_L_result$shapiro_p,
  Levene_p = spine_L_result$levene_p,
  Normality_ok = spine_L_result$normality_ok,
  Homogeneity_ok = spine_L_result$homogeneity_ok,
  Test_used = spine_L_result$test_used,
  stringsAsFactors = FALSE
)

# 8. body size_L
body_L_result <- analyze_dataset("body_MC.csv", "V", "Body Size", "MC-LR")
all_results[[18]] <- body_L_result$result_df
posthoc_results[["body_L"]] <- body_L_result$posthoc
assumption_results[["body_L"]] <- data.frame(
  Dataset = "body_L",
  Parameter = "Body Size",
  Treatment = "MC-LR",
  Shapiro_p = body_L_result$shapiro_p,
  Levene_p = body_L_result$levene_p,
  Normality_ok = body_L_result$normality_ok,
  Homogeneity_ok = body_L_result$homogeneity_ok,
  Test_used = body_L_result$test_used,
  stringsAsFactors = FALSE
)

# 9. lifespan_L
life_L_result <- analyze_dataset("lifespan_MC.csv", "time", "Lifespan", "MC-LR")
all_results[[19]] <- life_L_result$result_df
posthoc_results[["life_L"]] <- life_L_result$posthoc
assumption_results[["life_L"]] <- data.frame(
  Dataset = "life_L",
  Parameter = "Lifespan",
  Treatment = "MC-LR",
  Shapiro_p = life_L_result$shapiro_p,
  Levene_p = life_L_result$levene_p,
  Normality_ok = life_L_result$normality_ok,
  Homogeneity_ok = life_L_result$homogeneity_ok,
  Test_used = life_L_result$test_used,
  stringsAsFactors = FALSE
)

# 10. offspring_L
off_L_result <- analyze_dataset("offspring_MC.csv", "offspring", "Offspring", "MC-LR")
all_results[[20]] <- off_L_result$result_df
posthoc_results[["off_L"]] <- off_L_result$posthoc
assumption_results[["off_L"]] <- data.frame(
  Dataset = "off_L",
  Parameter = "Offspring",
  Treatment = "MC-LR",
  Shapiro_p = off_L_result$shapiro_p,
  Levene_p = off_L_result$levene_p,
  Normality_ok = off_L_result$normality_ok,
  Homogeneity_ok = off_L_result$homogeneity_ok,
  Test_used = off_L_result$test_used,
  stringsAsFactors = FALSE
)

#####################################################
# 合并所有结果到一个数据框
#####################################################

# 处理合并函数
process_results <- function(all_results) {
  # 找出所有结果中存在的列名
  all_cols <- unique(unlist(lapply(all_results, colnames)))
  
  # 创建空的合并数据框
  final_df <- data.frame()
  
  # 处理每个结果
  for (i in seq_along(all_results)) {
    df <- all_results[[i]]
    
    # 如果是Scheirer-Ray-Hare结果
    if ("Statistic" %in% colnames(df) && df$Test_Type[1] == "Scheirer-Ray-Hare") {
      # 添加缺失的Mean_Sq列
      df$Mean_Sq <- NA
      
      # 重新排列列顺序
      srh_df <- data.frame(
        Parameter = df$Parameter,
        Treatment = df$Treatment,
        Effect = df$Effect,
        Df = df$Df,
        Sum_Sq = df$Sum_Sq,
        Mean_Sq = df$Mean_Sq,
        Statistic = df$Statistic,
        p_value = df$p_value,
        Test_Type = df$Test_Type,
        stringsAsFactors = FALSE
      )
      final_df <- rbind(final_df, srh_df)
    }
    # 如果是ANOVA结果
    else if ("Statistic" %in% colnames(df) && df$Test_Type[1] == "ANOVA") {
      # 已经包含所有需要的列
      final_df <- rbind(final_df, df)
    }
  }
  
  # 移除临时行名
  rownames(final_df) <- NULL
  
  return(final_df)
}

# 合并假设检验结果
assumptions_df <- do.call(rbind, assumption_results)
rownames(assumptions_df) <- NULL

# 处理所有结果
final_results <- process_results(all_results)

# 查看假设检验结果
cat("\n=== 假设检验结果汇总 ===\n")
print(assumptions_df)

# 查看统计检验结果
cat("\n=== 处理后结果的前几行 ===\n")
print(head(final_results, 10))

# 查看数据结构
cat("\n=== 结果数据结构 ===\n")
str(final_results)

# 查看汇总统计
cat("\n=== 结果汇总 ===\n")
print(table(final_results$Test_Type))
print(table(final_results$Treatment))

# 保存结果到CSV文件
write.csv(final_results, "two_way_analysis_results.csv", row.names = FALSE)
cat("\n结果已保存到: two_way_analysis_results.csv\n")

# 保存假设检验结果
write.csv(assumptions_df, "assumption_test_results.csv", row.names = FALSE)
cat("假设检验结果已保存到: assumption_test_results.csv\n")

# 按处理类型分开保存
alga_results <- subset(final_results, Treatment == "Alga")
mc_lr_results <- subset(final_results, Treatment == "MC-LR")

write.csv(alga_results, "Alga_results.csv", row.names = FALSE)
write.csv(mc_lr_results, "MC_LR_results.csv", row.names = FALSE)

cat("\nAlga结果已保存到: Alga_results.csv")
cat("\nMC-LR结果已保存到: MC_LR_results.csv\n")

# 按检验类型分开保存
anova_results <- subset(final_results, Test_Type == "ANOVA")
srh_results <- subset(final_results, Test_Type == "Scheirer-Ray-Hare")

write.csv(anova_results, "ANOVA_results_all.csv", row.names = FALSE)
write.csv(srh_results, "Scheirer_Ray_Hare_results_all.csv", row.names = FALSE)

cat("\nANOVA结果已保存到: ANOVA_results_all.csv")
cat("\nScheirer-Ray-Hare结果已保存到: Scheirer_Ray_Hare_results_all.csv\n")

# 生成汇总表格
summary_table <- final_results %>%
  group_by(Parameter, Treatment, Test_Type) %>%
  summarize(
    Effects = n(),
    Significant_Effects = sum(p_value < 0.05, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\n=== 分析汇总表格 ===\n")
print(summary_table)

write.csv(summary_table, "analysis_summary.csv", row.names = FALSE)
cat("\n分析汇总表格已保存到: analysis_summary.csv\n")

# 生成格式化的报告表格
generate_report_table <- function(results_df) {
  report_df <- data.frame(
    Parameter = results_df$Parameter,
    Treatment = results_df$Treatment,
    Effect = results_df$Effect,
    Df = results_df$Df,
    Statistic = ifelse(results_df$Test_Type == "ANOVA", 
                       paste0("F = ", round(results_df$Statistic, 3)),
                       paste0("H = ", round(results_df$Statistic, 3))),
    p_value = ifelse(results_df$p_value < 0.001, "< 0.001",
                     ifelse(results_df$p_value < 0.01, round(results_df$p_value, 4),
                            round(results_df$p_value, 3))),
    Significance = ifelse(results_df$p_value < 0.001, "***",
                          ifelse(results_df$p_value < 0.01, "**",
                                 ifelse(results_df$p_value < 0.05, "*", "ns"))),
    Test_Type = results_df$Test_Type,
    stringsAsFactors = FALSE
  )
  
  return(report_df)
}

report_table <- generate_report_table(final_results)
cat("\n=== 格式化报告表格（前20行） ===\n")
print(head(report_table, 20))

write.csv(report_table, "formatted_report_table.csv", row.names = FALSE)
cat("\n格式化报告表格已保存到: formatted_report_table.csv\n")

# 生成交互作用效应的单独表格
interaction_effects <- final_results %>%
  filter(Effect %in% c("P:F", "P*F"))

cat("\n=== 交互作用效应汇总 ===\n")
print(interaction_effects)

write.csv(interaction_effects, "interaction_effects.csv", row.names = FALSE)
cat("\n交互作用效应结果已保存到: interaction_effects.csv\n")

# 生成主效应表格
main_effects <- final_results %>%
  filter(Effect %in% c("P", "F"))

cat("\n=== 主效应汇总 ===\n")
print(main_effects)

write.csv(main_effects, "main_effects.csv", row.names = FALSE)
cat("\n主效应结果已保存到: main_effects.csv\n")

# 保存事后检验结果（改进版本）
save_posthoc_summary <- function(posthoc_results) {
  posthoc_summary <- data.frame()
  
  for (name in names(posthoc_results)) {
    ph_result <- posthoc_results[[name]]
    
    if (inherits(ph_result, "pairwise.htest")) {
      # 对于Wilcoxon检验
      if (!is.null(ph_result$p.value)) {
        # 将p值矩阵转换为长格式
        p_mat <- ph_result$p.value
        p_df <- data.frame(
          Dataset = character(0),
          Group1 = character(0),
          Group2 = character(0),
          p_value = numeric(0),
          Test_Type = character(0),
          stringsAsFactors = FALSE
        )
        
        # 提取矩阵中的所有组合
        for (i in 1:nrow(p_mat)) {
          for (j in 1:ncol(p_mat)) {
            p_val <- p_mat[i, j]
            if (!is.na(p_val)) {
              row_name <- rownames(p_mat)[i]
              col_name <- colnames(p_mat)[j]
              p_df <- rbind(p_df, data.frame(
                Dataset = name,
                Group1 = row_name,
                Group2 = col_name,
                p_value = round(p_val, 4),
                Test_Type = "Wilcoxon",
                stringsAsFactors = FALSE
              ))
            }
          }
        }
        
        posthoc_summary <- rbind(posthoc_summary, p_df)
      }
    } else if (inherits(ph_result, "list") && "Tukey" %in% names(ph_result)) {
      # 对于Tukey检验
      # 尝试提取交互作用项P:F
      if ("P:F" %in% names(ph_result$Tukey)) {
        tukey_interaction <- as.data.frame(ph_result$Tukey$`P:F`)
        tukey_interaction$Comparison <- rownames(tukey_interaction)
        tukey_interaction$Dataset <- name
        tukey_interaction$Test_Type <- "TukeyHSD"
        
        # 重命名列
        colnames(tukey_interaction) <- c(
          "diff", "lwr", "upr", "p_adj",
          "Comparison", "Dataset", "Test_Type"
        )
        
        posthoc_summary <- rbind(posthoc_summary, tukey_interaction)
      }
      # 如果没有P:F交互作用，尝试提取其他交互作用项
      else if (any(grepl(":", names(ph_result$Tukey)))) {
        interaction_name <- names(ph_result$Tukey)[grepl(":", names(ph_result$Tukey))][1]
        tukey_interaction <- as.data.frame(ph_result$Tukey[[interaction_name]])
        tukey_interaction$Comparison <- rownames(tukey_interaction)
        tukey_interaction$Dataset <- name
        tukey_interaction$Test_Type <- "TukeyHSD"
        
        colnames(tukey_interaction) <- c(
          "diff", "lwr", "upr", "p_adj",
          "Comparison", "Dataset", "Test_Type"
        )
        
        posthoc_summary <- rbind(posthoc_summary, tukey_interaction)
      }
      # 如果没有交互作用项，尝试提取主效应
      else if ("P" %in% names(ph_result$Tukey)) {
        tukey_P <- as.data.frame(ph_result$Tukey$P)
        tukey_P$Comparison <- rownames(tukey_P)
        tukey_P$Dataset <- name
        tukey_P$Test_Type <- "TukeyHSD"
        
        colnames(tukey_P) <- c(
          "diff", "lwr", "upr", "p_adj",
          "Comparison", "Dataset", "Test_Type"
        )
        
        posthoc_summary <- rbind(posthoc_summary, tukey_P)
      }
    }
  }
  
  return(posthoc_summary)
}

posthoc_summary <- save_posthoc_summary(posthoc_results)
if (nrow(posthoc_summary) > 0) {
  write.csv(posthoc_summary, "posthoc_results_summary.csv", row.names = FALSE)
  cat("\n事后检验结果摘要已保存到: posthoc_results_summary.csv\n")
} else {
  cat("\n没有有效的事后检验结果\n")
}

cat("\n=== 分析完成 ===")
cat("\n所有结果已保存到CSV文件中。")
cat("\n请查看以下文件：")
cat("\n1. two_way_analysis_results.csv - 完整结果")
cat("\n2. assumption_test_results.csv - 假设检验结果")
cat("\n3. Alga_results.csv - Alga处理结果")
cat("\n4. MC_LR_results.csv - MC-LR处理结果")
cat("\n5. analysis_summary.csv - 分析汇总")
cat("\n6. formatted_report_table.csv - 格式化报告表格")
cat("\n7. interaction_effects.csv - 交互作用效应")
cat("\n8. main_effects.csv - 主效应结果")
cat("\n9. posthoc_results_summary.csv - 事后检验结果摘要\n")
