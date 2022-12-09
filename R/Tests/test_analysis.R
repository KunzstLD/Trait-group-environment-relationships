############################################################################
# Goal: Similarity in the response of traits and TPGs to (pesticide-stress)
#
# Approaches: 
#   - Clustering -> TPG -> Use relative abundance or summed abundance per group
#   - HMSC: Use as trait the toxicity?
# Q: How to make the comparison to single traits?
# What do traits do in HMSC?
# For species niches (Beta,ß) HMSC assumes multivariate normal variation
# The expected value (µ) is species-specific by considering species traits
# (see HMSC materials)
# -> How much traits explain species niches + how much do traits explain of variation in species abundances?

#   - Piliere Approach: Fit many BRT: Abund. per trait or TPG vs. environm. variables
#   - Correlation of responses of traits and TPGs
# TODO: Add trait on toxicity; TU sensitive?
############################################################################

# HMSC Example
# 
# https://cran.r-project.org/web/packages/Hmsc/vignettes/vignette_3_multivariate_high.pdf
library(Hmsc)
library(corrplot)
library(parallel)
library(MASS)
set.seed(1)

# Simulate phylogeny and trait data
ns = 50
phy = ape::rcoal(n=ns, tip.label = sprintf('species_%.3d',1:ns), br = "coalescent")
plot(phy, show.tip.label = FALSE, no.margin = TRUE)

C = ape::vcv(phy, model = "Brownian", corr = TRUE)
spnames = colnames(C)
traits = matrix(NA, ncol = 2, nrow = ns)
for (i in 1:2) {
  traits[, i] = matrix(mvrnorm(
    n = 1,
    mu = rep(0, ns),
    Sigma = C
  ))
}
rownames(traits) = spnames
colnames(traits) = c("habitat.use", "thermal.optimum")
traits = as.data.frame(traits)

par(fig = c(0,0.6,0,0.8), mar=c(6,0,2,0))
plot(phy, show.tip.label = FALSE)
par(fig = c(0.6,0.9,0.025,0.775), mar=c(6,0,2,0), new=T)
plot.new()
fields::image.plot(t(traits),axes=FALSE,legend.width = 3,legend.shrink=1,
           col = colorRampPalette(c("blue","white","red"))(200))
text(x=1.1, y=0.72, srt = 90, "H", cex=0.9, pos = 4)
text(x=1.4, y=0.72, srt = 90, "T", cex=0.9, pos = 4)

# Simulate env. and species data
n = 200
habitat = factor(sample(x = c("forest","open"), size = n, replace=TRUE))
climate = rnorm(n)

nc = 4
mu = matrix(0,nrow=nc,ncol=ns)
mu[1, ] = -traits$thermal.optimum^2/4-traits$habitat.use
mu[2, ] = 2*traits$habitat.use
mu[3, ] = traits$thermal.optimum/2
mu[4, ] = -1/4
beta = mu + 0.25*matrix(rnorm(n = ns*nc), ncol=ns)
X = cbind(rep(1,ns), as.numeric(habitat=="forest"), climate, climate*climate)
L = X%*%beta

Y = L + mvrnorm(n=n, mu=rep(0,ns), Sigma=diag(ns))
colnames(Y) = spnames

# HMSC model 
XData = data.frame(climate = climate, habitat = habitat)
XFormula = ~habitat + poly(climate,degree = 2,raw = TRUE)
TrFormula = ~habitat.use + thermal.optimum
studyDesign = data.frame(sample = sprintf('sample_%.3d',1:n), stringsAsFactors=TRUE)
rL = HmscRandomLevel(units = studyDesign$sample)
rL$nfMax = 15
m = Hmsc(Y = Y, XData = XData, XFormula = XFormula,
         TrData = traits, 
         TrFormula = TrFormula,
         phyloTree = phy,
         studyDesign = studyDesign, ranLevels = list(sample = rL))



