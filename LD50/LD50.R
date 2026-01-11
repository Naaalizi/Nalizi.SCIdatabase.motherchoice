rm(list = ls())

library(ggplot2)
library(dplyr)

# 设置工作目录
setwd("E:/FILE/Masterexp/Paper/1_mother_choice/1_data/picture/LD50/")

# 读取数据
data_raw <- read.csv("LD50_alga.csv")

# 查看数据结构
print("数据概览:")
print(head(data_raw))
print(str(data_raw))
print(summary(data_raw))

# 检查列名，确保与数据匹配
# 根据您提供的数据，列名应该是: alga, individual, deathrate
# 如果不是，请修改为正确的列名

# 重新命名列以便清晰
colnames(data_raw) <- c("alga", "individual", "deathrate")

# 如果需要，可以计算死亡率的平均值和标准差用于后续分析
# 但为了使用所有数据点进行拟合，我们将使用原始重复数据

# 使用所有数据点进行线性回归
# x: 死亡率(deathrate), y: 藻浓度(alga)
my_fit <- lm(alga ~ deathrate, data = data_raw)

# 获取模型摘要
fit_summary <- summary(my_fit)
print(fit_summary)

# 计算当死亡率为50%时的藻浓度预测值
y_val <- predict(my_fit, newdata = data.frame(deathrate = 50))
y_val <- round(y_val, 3)

# 计算R²
R_square <- fit_summary$r.squared
R_square <- round(R_square, 3)

# 提取回归系数
intercept <- coef(my_fit)[1]
slope <- coef(my_fit)[2]

# 四舍五入系数
intercept_rounded <- round(intercept, 3)
slope_rounded <- round(slope, 3)

# 创建公式表达式
# 注意：这里使用bquote创建表达式，以便在图中正确显示公式
formula_text <- bquote(y == .(intercept_rounded) + .(slope_rounded) ~ x)

# 创建图形
LD50_alga <- ggplot(data_raw, aes(x = deathrate, y = alga)) +
  # 绘制所有数据点
  geom_point( colour = "blue", size = 2) +
  
  # 添加线性回归线
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black",
    size = 1.2,
    formula = y ~ x
  ) +
  
  # 添加垂直线和水平线标记LD50
  annotate( "segment",  x = 50, y = 0, xend = 50, yend = y_val,lty = 5, size = 1, colour = "#426671"
  ) +
  annotate( "segment", x = 0, y = y_val,xend = 50, yend = y_val, lty = 5, size = 1, colour = "#426671"
  ) +
  
  # 标记LD50点
  annotate( "point", x = 50, y = y_val,color = "red", size = 3 ) +
  # 添加公式文本
  annotate( "text",   x = 25, y = 17,  parse = TRUE,
    size = 5,  label = paste0("y == ", intercept_rounded, " + ", slope_rounded, " * x")
  ) +
  
  # 添加R²文本
  annotate( "text",  x = 25, y = 15,
    parse = TRUE,
    size = 5,
    label = paste0("R^2 == ", R_square) ) +
  
  # 添加LD50坐标文本
  annotate( "text", x = 60, y = y_val-1, size = 5, label = paste0("(50, ", y_val, ")") ) +
  
  # 设置主题
  theme_bw() +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.line = element_line(),
    plot.margin = unit(rep(1, 4), 'mm') ) +
  
  # 设置坐标轴
  scale_x_continuous( breaks = seq(0, 100, 10),limits = c(0, 100)  ) +
  scale_y_continuous(breaks = seq(0, 20, 5),limits = c(0,20)) +
  
  # 设置坐标轴标签
  ylab( expression( atop( "Concentration of non-toxic ",
                      italic("M.aeruginosa") ~ "(" ~ 10^4 ~ "cells/mL)") ) ) +
  xlab("Death Rate(%)")

# 显示图形
print(LD50_alga)

# 保存图形
ggsave("alga_LD50_all_points.svg", width =180, height= 90 , units ="mm",dpi=1000)




###################################################MC
# 读取数据
data_raw <- read.csv("LD50_MC.csv")

# 查看数据结构
print("数据概览:")
print(head(data_raw))
print(str(data_raw))
print(summary(data_raw))

# 检查列名，确保与数据匹配


# 重新命名列以便清晰
colnames(data_raw) <- c("MC",  "deathrate")

# 如果需要，可以计算死亡率的平均值和标准差用于后续分析
# 但为了使用所有数据点进行拟合，我们将使用原始重复数据

# 使用所有数据点进行线性回归
# x: 死亡率(deathrate), y: MC
my_fit <- lm(MC ~ deathrate, data = data_raw)

# 获取模型摘要
fit_summary <- summary(my_fit)
print(fit_summary)

# 计算当死亡率为50%时的藻浓度预测值
y_val <- predict(my_fit, newdata = data.frame(deathrate = 50))
y_val <- round(y_val, 3)

# 计算R²
R_square <- fit_summary$r.squared
R_square <- round(R_square, 3)

# 提取回归系数
intercept <- coef(my_fit)[1]
slope <- coef(my_fit)[2]

# 四舍五入系数
intercept_rounded <- round(intercept, 3)
slope_rounded <- round(slope, 3)

# 创建公式表达式
# 注意：这里使用bquote创建表达式，以便在图中正确显示公式
formula_text <- bquote(y == .(intercept_rounded) + .(slope_rounded) ~ x)

# 创建图形
LD50_MC <- ggplot(data_raw, aes(x = deathrate, y = MC)) +
  # 绘制所有数据点
  geom_point( colour = "blue", size = 2) +
  
  # 添加线性回归线
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black",
    size = 1.2,
    formula = y ~ x
  ) +
  
  # 添加垂直线和水平线标记LD50
  annotate( "segment",  x = 50, y = 0, xend = 50, yend = y_val,lty = 5, size = 1, colour = "#426671"
  ) +
  annotate( "segment", x = 0, y = y_val,xend = 50, yend = y_val, lty = 5, size = 1, colour = "#426671"
  ) +
  
  # 标记LD50点
  annotate( "point", x = 50, y = y_val,color = "red", size = 3 ) +
  # 添加公式文本
  annotate( "text",   x = 25, y = 500,  parse = TRUE,
            size = 5,  label = paste0("y == ", intercept_rounded, " + ", slope_rounded, " * x")
  ) +
  
  # 添加R²文本
  annotate( "text",  x = 25, y = 450,
            parse = TRUE,
            size = 5,
            label = paste0("R^2 == ", R_square) ) +
  
  # 添加LD50坐标文本
  annotate( "text", x = 63, y = y_val-40, size = 5, label = paste0("(50, ", y_val, ")") ) +
  
  # 设置主题
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.line = element_line(),
        plot.margin = unit(rep(1, 4), 'mm') ) +
  
  # 设置坐标轴
  scale_x_continuous( breaks = seq(0, 100, 10),limits = c(0, 100)  ) +
  scale_y_continuous(breaks = seq(0, 600, 100),limits = c(0,600)) +
  
  # 设置坐标轴标签
  ylab( " \nConcentration of MC-LR (ng/mL)") +
  xlab("Death Rate(%)")

# 显示图形
print(LD50_MC)

# 保存图形
ggsave("MC_LD50_all_points.svg", width =180, height= 90 , units ="mm",dpi=1000)

