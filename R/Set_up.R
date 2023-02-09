# ______________________________________________________________________
# Set up script
# loading libraries, paths & functions
# ______________________________________________________________________

## Libraries

# Data processing
library(data.table)
library(dplyr)
library(readxl)
library(purrr)
library(zeallot)
library(purrr)

# Retrieval of taxonomic information and habitat information
library(taxize)
library(worms)

# Retrieval of ecotoxicological test  information
library(standartox)
library(chemlook)

# Plotting
library(ggplot2)
library(reactable)
library(htmltools)
library(patchwork)
library(ggsci)

# Distance matrices and clustering
library(ade4)
library(cluster)
library(NbClust)
library(dendextend)

# Boosted regression trees
library(mlr3)
library(mlr3learners)
library(mlr3viz)
library(mlr3verse)
library(xgboost)
library(paradox)
library(DALEX)
library(DALEXtra)

## Paths
path_in <- "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Data"
path_cache <- "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Cache"
path_scr <- "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/R"
path_repo <- "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships"
path_out <- "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Output"
path_paper <- "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Paper"

## Functions
source("/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/R/Functions.R")

## North American harmonised trait dataset
# TODO: necessary to load this?
noa_traits <- readRDS(file.path(path_in, "NoA", "Traits_US_LauraT_pp_harmonized.rds"))
noa_traits[, taxon := coalesce(species, genus, family, order)]
noa_traits <- noa_traits[, .SD,
    .SDcols = patterns("species|genus|family|order|taxon|feed|resp|size|locom|volt")
]
