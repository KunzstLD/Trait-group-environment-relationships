library(glmmTMB)
library(data.table)
library(dplyr)

# Data
dat <- read.csv("/home/kunzst/Dokumente/Literature/Ecological_Literature/TerBraak_2019/Tutorial/Data/whittakerrevisitdata.csv")
str(dat)
summary(dat)

# Y is response variable
# for logit, must be matrix of successes and failures
dat$y <- with(dat, cbind(value,100-value))
dat$site <- factor(dat$site)

# Model
formula.MLM3.linear <- y ~ env * trait + (1 + trait| site) + (1 + env|species)
MLM3 <- glmmTMB(formula.MLM3.linear, family = betabinomial, data= dat)
summary(MLM3)
