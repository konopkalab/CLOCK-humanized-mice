# STATISTICS
- [01_Statistics.R](01_Statistics.R)

```{R}
library(multcomp)
library(emmeans)
library(lme4)
library(lmerTest)

#set your working directory
setwd("your working directory")
list <- list.files(path = "your working directory")

tab <- read.table(file = list[i], sep = "\t", header = T, fill = TRUE)
#this step is only required if the range of your value is between 0 and 1
tab <- log(tab)

#General linear mixed model (GLMM) and general linear model (GLM), the model depends on the data for testing
model <- lmer(depentent ~ fixed + sex + (1|random/image), data=tab)
model <- lmer(depentent ~ fixed + sex +time + (1|random/neuron/segment/spine), data=tab)
model <- lm(depentent ~ fixed, data=tab)

main <- summary(model)
posthoc <- emmeans(model, pairwise ~ fixed, adjust = "tukey")
anova <- anova(model)

setwd("file saving directory")
capture.output(main, posthoc, anova, file = "statistical_result.txt")

#repeated measures ANOVA 
setwd("your working directory")
tab <- read.table(file = list[i], sep = "\t", header = T, fill = TRUE)
res.aov <- anova_test(data = tab, dv = observation, wid = id, within = c("time", "sex"))
res.aov <- anova_test(data = tab, observation ~ time * fixed * sex, wid = id)
aov <- get_anova_table(res.aov)
#posthoc test
pwc <- tab %>%
  pairwise_t_test(observation ~ fixed, paired = F, p.adjust.method = "bonferroni")

setwd("file saving directory")
capture.output(aov, pwc, file = "statistical_result.txt")

```
