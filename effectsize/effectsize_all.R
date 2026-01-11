rm(list=ls())
library(effsize)

setwd("D:/SSData/effectsize/growthrate/")

growthrate <- read.csv("effectsize_of_growthrate.csv")

d1 <- cohen.d(growthrate$MM,growthrate$CM,data=growthrate,conf.level = 0.95)
d2 <- cohen.d(growthrate$LL,growthrate$CL,data=growthrate,conf.level = 0.95)

# 打印效应量和置信区间
print(d1)
print(d2)

d1$estimate
d1$conf.int["lower"]
d1$conf.int["upper"]

# 定义两个变量的均值和置信区间界限
mean1 <- d1$estimate
ci_lower1 <- d1$conf.int["lower"]
ci_upper1 <- d1$conf.int["upper"]

mean2 <- d2$estimate
ci_lower2 <- d2$conf.int["lower"]
ci_upper2 <- d2$conf.int["upper"]

# 计算每个均值的标准误差
se1 <- (ci_upper1 - ci_lower1) / (2 * 1.96)
se2 <- (ci_upper2 - ci_lower2) / (2 * 1.96)

# 计算Z得分
z_score <- (mean1 - mean2) / sqrt(se1^2 + se2^2)

# 计算p值
p_value <- 2 * (1 - pnorm(abs(z_score)))

# 输出结果
print(paste("Z Score:", z_score))
print(paste("P Value:", p_value))

library(ggplot2)
library(ggsignif)

df <- data.frame('group'=c('Alga','MC-LR'),'mean'=c(mean1,mean2),'se'=c(se1,se2))

# 将group转换为因子，并设置因子水平
df$group <- factor(df$group, 
                   levels = c('Alga','MC-LR'),
                   labels = c(
                     expression(atop("non-toxic", italic("M. aeruginosa"))),
                     "MC-LR"
                   ))

pr <- ggplot(df, aes(x = group, y = mean, shape = group)) +
  geom_point(size = 2) +
  geom_signif(  
    comparisons = list(c(1, 2)),  # 使用数字索引
    y_position = 6.5, 
    map_signif_level = FALSE, 
    annotations = c("p<0.001"),
    textsize = 2.5,
    size = 0.2,
    tip_length = c(0.9, 0.1)) +
  geom_hline(aes(yintercept = 0), size = 0.5, colour = "dark grey", linetype = "dashed") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), size = 0.3, width = 0.1) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),  # 调整行间距
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(rep(1,4), 'mm'),
        legend.position = "none") +
  scale_y_continuous(breaks = seq(-4, 8, 2), expand = c(0, 0), limits = c(-4, 8)) + 
  scale_shape_manual(values = c(16, 1)) +
  labs(x = "Group", y = c("Maternal effect size\nof growth rate")) +
  scale_x_discrete(labels = function(x) parse(text = as.character(x)))

pr

setwd("D:/SSData/effectsize/instrincrate/")

instrincrate <- read.csv("effectsize_of_instrincrate.csv")

d1 <- cohen.d(instrincrate$MM,instrincrate$CM,data=instrincrate,conf.level = 0.95)
d2 <- cohen.d(instrincrate$LL,instrincrate$CL,data=instrincrate,conf.level = 0.95)

# 打印效应量和置信区间
print(d1)
print(d2)

d1$estimate
d1$conf.int["lower"]
d1$conf.int["upper"]

# 定义两个变量的均值和置信区间界限
mean1 <- d1$estimate
ci_lower1 <- d1$conf.int["lower"]
ci_upper1 <- d1$conf.int["upper"]

mean2 <- d2$estimate
ci_lower2 <- d2$conf.int["lower"]
ci_upper2 <- d2$conf.int["upper"]

# 计算每个均值的标准误差
se1 <- (ci_upper1 - ci_lower1) / (2 * 1.96)
se2 <- (ci_upper2 - ci_lower2) / (2 * 1.96)

# 计算Z得分
z_score <- (mean1 - mean2) / sqrt(se1^2 + se2^2)

# 计算p值
p_value <- 2 * (1 - pnorm(abs(z_score)))

# 输出结果
print(paste("Z Score:", z_score))
print(paste("P Value:", p_value))

library(ggplot2)
library(ggsignif)

df <- data.frame('group'=c('Alga','MC-LR'),'mean'=c(mean1,mean2),'se'=c(se1,se2))

# 将group转换为因子，并设置因子水平
df$group <- factor(df$group, 
                   levels = c('Alga','MC-LR'),
                   labels = c(
                     expression(atop("non-toxic", italic("M. aeruginosa"))),
                     "MC-LR"
                   ))

prm <- ggplot(df, aes(x = group, y = mean, shape = group)) +
  geom_point(size = 2) +
  geom_signif(  
    comparisons = list(c(1, 2)),  # 使用数字索引
    y_position = 6, 
    map_signif_level = FALSE, 
    annotations = c("p<0.001"),
    textsize = 2.5,
    size = 0.2,
    tip_length = c(0.9, 0.12)) +
  geom_hline(aes(yintercept = 0), size = 0.5, colour = "dark grey", linetype = "dashed") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), size = 0.3, width = 0.1) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),  # 调整行间距
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(rep(1,4), 'mm'),
        legend.position = "none") +
  scale_y_continuous(breaks = seq(-2, 8, 2), expand = c(0, 0), limits = c(-3, 7)) + 
  scale_shape_manual(values = c(16, 1)) +
  labs(x = "Group", y = c("Maternal effect size\nof instrinc rate")) +
  scale_x_discrete(labels = function(x) parse(text = as.character(x)))

prm

setwd("D:/SSData/effectsize/R0/")

R0 <- read.csv("effectsize_of_R0.csv")

d1 <- cohen.d(R0$MM,R0$CM,data=R0,conf.level = 0.95)
d2 <- cohen.d(R0$LL,R0$CL,data=R0,conf.level = 0.95)

# 打印效应量和置信区间
print(d1)
print(d2)

d1$estimate
d1$conf.int["lower"]
d1$conf.int["upper"]

# 定义两个变量的均值和置信区间界限
mean1 <- d1$estimate
ci_lower1 <- d1$conf.int["lower"]
ci_upper1 <- d1$conf.int["upper"]

mean2 <- d2$estimate
ci_lower2 <- d2$conf.int["lower"]
ci_upper2 <- d2$conf.int["upper"]

# 计算每个均值的标准误差
se1 <- (ci_upper1 - ci_lower1) / (2 * 1.96)
se2 <- (ci_upper2 - ci_lower2) / (2 * 1.96)

# 计算Z得分
z_score <- (mean1 - mean2) / sqrt(se1^2 + se2^2)

# 计算p值
p_value <- 2 * (1 - pnorm(abs(z_score)))

# 输出结果
print(paste("Z Score:", z_score))
print(paste("P Value:", p_value))

library(ggplot2)
library(ggsignif)

df <- data.frame('group'=c('Alga','MC-LR'),'mean'=c(mean1,mean2),'se'=c(se1,se2))

# 将group转换为因子，并设置因子水平
df$group <- factor(df$group, 
                   levels = c('Alga','MC-LR'),
                   labels = c(
                     expression(atop("non-toxic", italic("M. aeruginosa"))),
                     "MC-LR"
                   ))

pR0 <- ggplot(df, aes(x = group, y = mean, shape = group)) +
  geom_point(size = 2) +
  geom_signif(  
    comparisons = list(c(1, 2)),  # 使用数字索引
    y_position = 5.5, 
    map_signif_level = TRUE, 
    annotations = c("p<0.001"),
    textsize = 2.5,
    size = 0.2,
    tip_length = c(0.95, 0.15)) +
  geom_hline(aes(yintercept = 0), size = 0.5, colour = "dark grey", linetype = "dashed") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), size = 0.3, width = 0.1) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),  # 调整行间距
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(rep(1,4), 'mm'),
        legend.position = "none") +
  scale_y_continuous(breaks = seq(-2, 6, 2), expand = c(0, 0), limits = c(-3, 7)) + 
  scale_shape_manual(values = c(16, 1)) +
  labs(x = "Group", y = c("Maternal effect size\nof net reproduction rate")) +
  scale_x_discrete(labels = function(x) parse(text = as.character(x)))

pR0

library(ggpubr)

ggarrange(prm, pR0,
          labels = c("A", "B"), label.x = 0, label.y = 1,
          font.label = list(size = 10),
          align = "v",
          widths = c(1, 1),
          heights = c(1, 1),
          ncol = 2, nrow = 1)

ggsave("effectsize.svg", width = 160, height = 70, units = "mm", dpi = 1000)
