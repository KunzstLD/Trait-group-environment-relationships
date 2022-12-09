# Approaches to study TPG - environment relationships ----
# Q: How do species/taxa with similar trait expressions/trait profiles 
# respond to environmental gradients?
# Kleyer et al.: 
library(ade4)
library(MASS)
library(vegan)
library(ecodist)
library(maptools)
library(rpart)
library(splines)
library(gam)
library(pgirmess)
library(utils)
library(combinat)
library(mvpart)
library(cluster)
library(fpc)
library(clusterSim)
library(lmtest)
library(Hmisc)


install.packages(
  c(
    "ade4",
    "vegan",
    "ecodist",
    "maptools",
    "gam",
    "pgirmess",
    "combinat",
    "mvpart",
    "fpc",
    "clusterSim",
    "lmtest"
  )
)

# TODO Check Piliere 

## "CWM"-like approach ----
# E.g., CWM-RDA -> seems to be a bit different than what is proposed in Beauchard et al. 2017
# group species according to trait profile (gr. 1, gr. 2, ...)
# then record abundance of traits 
# then use sum or frequency and relate to environmental variable 
# This is an analysis on community level and thus does not reflect Q

## Cluster regression ----
# 1) - Cluster traits table -> species with similar trait expressions 
#    - Optimal number of clusters 
#    - Resulting clusterings and their clusters are bootstrapped to check cluster stability

# 2) Group responses to env. variables are modeled
# As response, kleyer seemed to use some kind of frequency on how often each group is present
# at a certain plot (for each group, nr. of species was counted in each site/plot weighted by abundance) 
# These counts act as dependent variable regressed from environmental variables
# Logistic regression (?) and model selection (Lasso?) 






## Double CCA ----
# Unimodal env. gradient
# How to incorporate functional groups a-priori?
# Looks like first are species and env. correlated, later traits come into play
# (same for RLQ?)

## RLQ ----
# Unimodal env. gradient
# How to incorporate functional groups a-priori?
