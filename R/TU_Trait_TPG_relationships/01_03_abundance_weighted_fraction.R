# ________________________________________________
# Abundance weighted fraction of TPGs ----
# Same families as in the CWM approach should be covered
## _______________________________________________

# Abundance weighted fraction
# TODO:
# - Abundance data, connect occurring taxa
# with groups (start with family dataset)
# - weight by abundance
trait_family <- readRDS(file.path(path_cache, "trait_family_tpg.rds"))
abund <- readRDS(file.path(path_cache, "total_abund_CEOPT.rds"))
abund <- abund[Region == "California", .(
    site,
    species,
    genus,
    family,
    order,
    taxon,
    taxonomic_level,
    abundance
)]

# Merge TPG
abund[trait_family, group := i.group, on = "family"]

# 5 families which are not classified to a certain TPG
# (California)
unique(abund[is.na(group), .(family, order)])

trait_matrix

trait_matrix[family %in% unique(abund[is.na(group), family]), ]
