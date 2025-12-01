# 构造数据框
df <- data.frame(
  Position = c(150, 800, 16000),
  AltCount = c(5, 20, 2),
  TotalCount = c(500, 1000, 400)
)

# 设置背景错误率（empirical error rate）
background_error_rate <- 0.001

# 用 binom.test 来计算每个位点是否高于背景错误率
binom.test(
  x = 5,
  n = 500,
  p = background_error_rate,
)
df$p_value <- mapply(function(alt, total) {
  binom.test(x = alt, n = total, p = background_error_rate, alternative = "greater")$p.value
}, df$AltCount, df$TotalCount)

# 添加突变频率
df$MutationFreq <- df$AltCount / df$TotalCount
df
# 查看结果
print(df)




# ? beta-binomial --------------------------------------------------------------------


# 测序深度（每个位点的总reads数）
depth <- c(1000, 1200, 950, 1100, 1050, 1300, 900, 1000, 1150, 1020)

# 变异计数（低频替代等位基因数）
alt_counts <- c(1, 2, 0, 3, 1, 2, 1, 0, 2, 1)

# 变异频率
allelic_ratios <- alt_counts / depth
print(allelic_ratios)
df <- data.frame(success = alt_counts, size = depth)

df
library(VGAM)

fit <- vglm(cbind(success, size - success) ~ 1, betabinomial, data = df)
summary(fit)
coef_values <- Coef(fit)
alpha <- exp(coef_values["logitmu"]) * exp(coef_values["logitphi"])
beta <- exp(coef_values["logitphi"])

coef_values <- Coef(fit)
alpha <- exp(coef_values["logitmu"]) * exp(coef_values["logitphi"])
beta <- exp(coef_values["logitphi"])

# 计算平均成功概率（估计的错误率）
error_rate <- alpha / (alpha + beta)
print(paste("Estimated background error rate:", error_rate))
