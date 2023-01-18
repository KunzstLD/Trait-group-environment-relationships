# __________________________________________________________________________________________________
# Functional groups ----
# - Use cluster analysis to delineate groups with similar trait profiles
# - Then use the abundance per TPG

# Idea: Generate general clusters based on the North American trait database:
# on both genus and on family-level
# compare both classifications to individual traits
## __________________________________________________________________________________________________

# North American trait databases (CONUS + Vieiera) Trait DB
noa_trait_conus <- readRDS(file.path(path_in, "NoA", "CONUS_trait_db_preproc.rds"))

# Bob Zuelig
noa_zuelig <- readRDS(file.path(path_in, "NoA", "trait_db_enhanced_preproc.rds"))
setnames(noa_zuelig,
         names(noa_zuelig),
         tolower(names(noa_zuelig)))
# S_org
# Not much information on species level anyways -> add on genus lvl
tab_spear <-
  fread(file.path(path_in, "sensitivity_values_SPEAR.csv"))
tab_spear[, Taxon := sub(" Gen\\. sp\\.| sp\\.| ssp\\.", "", Taxon)]
tab_spear[Familie == "Stratiomyiidae", `:=`(Gattung = NA_character_,
                                            Familie = "Stratiomyidae",
                                            Taxon = "Stratiomyidae")]

# Data preproc ----
# Subset to relevant orders
# Aggregate conus to genus lvl
# Add Bob Zuelig DB
noa_trait_conus <- noa_trait_conus[order %in% c("Coleoptera",
                                                "Diptera",
                                                "Ephemeroptera",
                                                "Odonata",
                                                "Plecoptera",
                                                "Trichoptera"),]
noa_zuelig <-
  noa_zuelig[order %in% c("Coleoptera",
                                         "Diptera",
                                         "Ephemeroptera",
                                         "Odonata",
                                         "Plecoptera",
                                         "Trichoptera"), ]

# Subset to taxa that are not present in CONUS
noa_zuelig_subset <-
  noa_zuelig[!genus %in% noa_aggr_conus$genus, .SD, .SDcols = col_names]

# Aggregation
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
# Identify taxa where information for one grouping feature is missing 
# need to add trait names, do not consider now
# noa_aggr_conus[, row_sum := apply(.SD, 1, function(y) sum(y, na.rm = TRUE)), .SDcols = trait_names[-20]]
# noa_aggr_conus[row_sum == 4, ]
noa_aggr_conus <- na.omit(noa_aggr_conus)
col_names <- names(noa_aggr_conus)[names(noa_aggr_conus) %in% names(noa_zuelig_subset)]

trait_genera <- rbind(noa_aggr_conus[, .SD, .SDcols = col_names],
                      noa_zuelig_subset[, .SD, .SDcols = col_names])

# Add S_org
trait_genera[tab_spear, sensitivity_organic := Sensitivität,
             on = c(genus = "Taxon")]

# S_org from family level to genus level
genera_subset <- trait_genera[is.na(sensitivity_organic), ] 
family_assign_sensitivity <-
  tab_spear[Taxon %in% unique(trait_genera[is.na(sensitivity_organic), family]),]
genera_subset[family_assign_sensitivity, sensitivity_organic := i.Sensitivität,
              on = c("family" = "Familie")]
trait_genera[genera_subset, sensitivity_organic := i.sensitivity_organic, on = "genus"]

# S_org from order level to genus level
family_subset <- trait_genera[is.na(sensitivity_organic), ]
order_assign_sensitivity <-
  tab_spear[Taxon %in% unique(trait_genera[is.na(sensitivity_organic), order]),]
family_subset[order_assign_sensitivity, sensitivity_organic := i.Sensitivität,
              on = c("order" = "Ordnung")]
trait_genera[family_subset, sensitivity_organic := i.sensitivity_organic, on = "family"]

# Few NA's left for sensitivity organic (Diptera), need to be removed because we don't
# have any information for these taxa 
trait_genera <- trait_genera[!is.na(sensitivity_organic), ]

# Overview range of taxonomic information ----
# Overview: 782 genera of 88 families
# Mostly Dipterans (279, then Trichoptera 140, ...Odonata 59)
trait_genera[, .N, by = "family"]
trait_genera[, .N, by = "order"]

# By comparison in CONUS:
# 799 genera of 112 families
noa_trait_conus[is.na(species) & !is.na(genus), ]
noa_trait_conus[is.na(species) & !is.na(genus), .N, by = "family"]

# And Bob Zuelig DB:
# 183 genera, 110 not in CONUS
# 65 families -> all in CONUS
# overall: 799+110 = 909 genera of 112 families
# i.e., we cover 86 % of all genera listed in both trait DBs 
# and 78.5 % of families
noa_zuelig[, .N, by = "family"]
noa_zuelig[!genus %in% noa_trait_conus$genus, ]
noa_zuelig[!family %in% noa_trait_conus$family, ]


