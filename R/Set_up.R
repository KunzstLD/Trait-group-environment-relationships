# ______________________________________________________________________
# Set up script
# loading libraries, paths & functions
# ______________________________________________________________________

## Libraries

# Data processing
library(data.table)
library(dplyr)
library(readxl)
library(zeallot)
library(purrr)

# Retrieval of taxonomic information and habitat information
# library(taxize)
# library(worms)

# Retrieval of ecotoxicological test  information
library(standartox)
# library(chemlook)

# Plotting & tables
library(ggplot2)
library(reactable)
library(htmltools)
library(patchwork)
library(ggsci)
library(GGally)
library(broom)

# Distance matrices and clustering
library(ade4)
library(cluster)
library(NbClust)
library(dendextend)
library(mclust)

# Boosted regression trees
library(mlr3)
library(mlr3learners)
library(mlr3viz)
library(mlr3verse)
library(xgboost)
library(paradox)
library(DALEX)
library(DALEXtra)

# LM and GAMs
library(DHARMa)
library(mgcv)
library(gratia)

# RDA
library(vegan)
library(corrplot)

## Paths
path_base <- "C:/Users/sk193/Documents/PhD/Trait-group-environment-relationships"
path_in <- file.path(path_base, "Data")
path_cache <- file.path(path_base, "Cache")
path_src <- file.path(path_base, "R")
path_out <- file.path(path_base, "Output")
path_paper <- file.path(path_base, "Paper")

## Functions
source(file.path(path_src, "Functions.R"))

## North American harmonised trait dataset
# TODO: necessary to load this?
noa_traits <- readRDS(file.path(path_in, "NoA", "Traits_US_LauraT_pp_harmonized.rds"))
noa_traits[, taxon := coalesce(species, genus, family, order)]
noa_traits <- noa_traits[, .SD,
    .SDcols = patterns("species|genus|family|order|taxon|feed|resp|size|locom|volt")
]
