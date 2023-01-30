# ________________________________________________
# Functional groups ----
# Preprocessing script
# Idea: Generate general clusters based on the North American trait database:
# on both genus and on family-level
# compare both classifications to individual traits
# TODO: Switch S_org -> assign it in the beginning
## _______________________________________________

# North American trait databases (CONUS + Vieiera) Trait DB
noa_trait_conus <-
  readRDS(file.path(path_in, "NoA", "CONUS_trait_db_preproc.rds"))

# Bob Zuelig
noa_zuelig <-
  readRDS(file.path(path_in, "NoA", "trait_db_enhanced_preproc.rds"))
setnames(noa_zuelig,
         names(noa_zuelig),
         tolower(names(noa_zuelig)))

# S_org
tab_spear <- readRDS(file.path(path_cache, "tab_spear.rds"))

# Add values for certain Diptera families
# Der Wert -0.35 gilt dennoch mit Ausnahme von:
# Chironomidae -0.39
# Culicidae  -0.29
# Simuliidae -0.46
tab_spear[
  Ordnung == "Diptera" & Rang == "group",
  Sensitivität := -0.35
]
tab_spear <- rbind(
  tab_spear,
  data.table(
    Ordnung = "Diptera",
    Familie = "Culicidae",
    Gattung = NA_character_,
    Taxon = "Culiciadae", Sensitivität = -0.29
  ),
  fill = TRUE
)

# Data preproc ----
# Subset to relevant orders
# Aggregate conus to genus lvl
# Add Bob Zuelig DB
noa_trait_conus <- noa_trait_conus[order %in% c("Coleoptera",
                                                "Diptera",
                                                "Ephemeroptera",
                                                "Odonata",
                                                "Plecoptera",
                                                "Trichoptera"), ]
noa_zuelig <-
  noa_zuelig[order %in% c("Coleoptera",
                          "Diptera",
                          "Ephemeroptera",
                          "Odonata",
                          "Plecoptera",
                          "Trichoptera"),]

## Aggregation to genus level ----
noa_aggr_conus <- direct_agg(
  trait_data = noa_trait_conus,
  non_trait_cols = c(
    "species",
    "genus",
    "family",
    "order",
    "taxon",
    "taxonomic_level"
  ),
  method = median,
  on = "genus"
)
noa_aggr_conus <- na.omit(noa_aggr_conus)

# Subset Zuelig DB to taxa that are not present
# in the aggregated North American datasets
col_names <-
  names(noa_aggr_conus)[names(noa_aggr_conus) %in% names(noa_zuelig)]
noa_zuelig_subset <-
  noa_zuelig[!is.na(genus) & !genus %in% noa_aggr_conus$genus,
    .SD,
    .SDcols = col_names
  ]

# Final genus trait dataset
trait_genera <- rbind(noa_aggr_conus[, .SD, .SDcols = col_names],
                      noa_zuelig_subset[, .SD, .SDcols = col_names])

# Add S_org
trait_genera[tab_spear, sensitivity_organic := i.Sensitivität,
  on = c(genus = "Taxon")
]

# Add S_org from coarser taxonomic levels to finer taxonomic levels
# S_org from family level to genus level
genera_subset <- trait_genera[is.na(sensitivity_organic), ]
family_assign_sensitivity <-
  tab_spear[Taxon %in% unique(trait_genera[is.na(sensitivity_organic), family]), ]
genera_subset[family_assign_sensitivity, sensitivity_organic := i.Sensitivität,
  on = c("family" = "Familie")
]
trait_genera[genera_subset, sensitivity_organic := i.sensitivity_organic,
  on = "genus"
]

# S_org from order level to genus level
genera_subset_order <- trait_genera[is.na(sensitivity_organic), ]
order_assign_sensitivity <-
  tab_spear[Taxon %in% unique(trait_genera[is.na(sensitivity_organic), order]), ]
genera_subset_order[order_assign_sensitivity, sensitivity_organic := i.Sensitivität,
  on = c("order" = "Ordnung")
]
trait_genera[genera_subset_order, sensitivity_organic := i.sensitivity_organic,
  on = "genus"
]
# saveRDS(trait_genera, file.path(path_cache, "trait_genera.rds"))

# Overview range of taxonomic information ----
# Overview: 712 genera of 95 families
trait_genera[, .N, by = "family"]
trait_genera[, .N, by = "order"]

# By comparison in North American trait DB:
# 799 genera of 112 families
# In total 133 families
noa_trait_conus[is.na(species) & !is.na(genus),]
noa_trait_conus[is.na(species) & !is.na(genus), .N, by = "family"]
noa_trait_conus[, uniqueN(family)]

# Bob Zuelig DB:
# 534 genera of 65 families (all families in North American trait DB)
# 21 Genera not in North American trait DB
# overall: 799+21 = 820 genera of 112 families
# i.e., we cover 87 % of all genera listed in both trait DBs
# from 85 % (88 of 112) of the families that have information up to genus level.
# Overall, 71 % (88 of 133) off all families 
# are represented with the genus dataset
noa_zuelig[, uniqueN(genus)]
noa_zuelig[, uniqueN(family)]
noa_zuelig[!is.na(genus) & !genus %in% noa_trait_conus$genus,]
noa_zuelig[!family %in% noa_trait_conus$family,]

# Aggregation to family level ----
# Add genera from Zuelig that are not contained 
# in the North American trait database
noa_trait_conus_zuelig <- rbind(noa_trait_conus,
  noa_zuelig[!genus %in% noa_trait_conus$genus, .SD, .SDcols = col_names],
  fill = TRUE
)
trait_family <- direct_agg(
  trait_data = noa_trait_conus_zuelig,
  non_trait_cols = c(
    "species",
    "genus",
    "family",
    "order",
    "taxon",
    "taxonomic_level"
  ),
  method = median,
  on = "family"
)

# Assign S_org on family level
trait_family[tab_spear, sensitivity_organic := i.Sensitivität,
  on = c("family" = "Familie")
]

# S_org from order level to family level
family_subset <- trait_family[is.na(sensitivity_organic), ]
order_assign_sensitivity <-
  tab_spear[Taxon %in% unique(trait_family[is.na(sensitivity_organic), order]), ]
family_subset[order_assign_sensitivity, sensitivity_organic := i.Sensitivität,
  on = c("order" = "Ordnung")
]
trait_family[family_subset, sensitivity_organic := i.sensitivity_organic,
  on = "family"
]

# 101 families of 133 possible (76 %)
trait_family <- na.omit(trait_family)
# saveRDS(trait_family, file.path(path_cache, "trait_family.rds"))
