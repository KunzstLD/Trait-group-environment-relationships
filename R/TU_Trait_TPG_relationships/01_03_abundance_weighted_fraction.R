# ________________________________________________
# Abundance weighted fraction of TPGs ----
# Same families as in the CWM approach should be covered
# (California)
## _______________________________________________

# Abundance weighted fraction
trait_family <- readRDS(file.path(path_cache, "trait_family_tpg.rds"))
trait_genera <- readRDS(file.path(path_cache, "trait_genera_tpg.rds"))
abund <- readRDS(file.path(path_cache, "total_abund_CEOPT.rds"))
unique(abund[, .(taxon, taxonomic_level)]) %>% 
.[, .N, by = "taxonomic_level"]
trait_genera$order %>% unique


# Merge TPGs
# family level
abund[trait_family, group_family := i.group, on = "family"]
# abund[, .N, by = c("Region", "group")] %>%
#   .[order(group), ]

# genus level
abund[trait_genera, group_genus := i.group, on = "genus"]

# Families which are not classified to a certain TPG
# throw out for now
# unique(abund[is.na(group), .(family, order), by = "Region"])
# also not used in CWM approach?

# Some taxa have a relatively high abundance (up to > 8000)
# use sqrt tranformation
# abund[, .(min(abundance), max(abundance)), by = "Region"]
abund[, abundance := sqrt(abundance)]

# Abundance weighted fraction
abund_family <- abund[!is.na(group_family),]
abund_genus <- abund[!is.na(group_genus),]
setnames(abund_family, "group_family", "group")
setnames(abund_genus, "group_genus", "group")
abund_comb <- list(
  "family_lvl" = abund_family,
  "genus_lvl" = abund_genus
)

# Taxonomic composition
# family level
abund_family[, .(
  Region,
  group,
  taxon,
  species,
  genus,
  family,
  order,
  abundance
)] %>%
  .[abundance > 0, ] %>%
  .[order(-abundance), ] %>%
  .[, sum_abund_reg_gr_ord := sum(abundance),
    by = c("Region", "group", "order")
  ] %>%
  .[, sum_abund_reg_gr := sum(abundance), by = c("Region", "group")] %>%
  .[, prop_ord := sum_abund_reg_gr_ord / sum_abund_reg_gr] %>%
  .[, .(Region, group, order, prop_ord)] %>%
  unique(.) %>%
  .[order(Region, group, -prop_ord), ] %>% 
  saveRDS(., file.path(path_cache, "tpg_taxonomic_composition_family.rds"))

abund_genus[, .(
  Region,
  group,
  taxon,
  species,
  genus,
  family,
  order,
  abundance
)] %>%
  .[abundance > 0, ] %>%
  .[order(-abundance), ] %>%
  .[, sum_abund_reg_gr_ord := sum(abundance),
    by = c("Region", "group", "order")
  ] %>%
  .[, sum_abund_reg_gr := sum(abundance), by = c("Region", "group")] %>%
  .[, prop_ord := sum_abund_reg_gr_ord / sum_abund_reg_gr] %>%
  .[, .(Region, group, order, prop_ord)] %>%
  unique(.) %>%
  .[order(Region, group, -prop_ord), ] %>% 
  saveRDS(., file.path(path_cache, "tpg_taxonomic_composition_genus.rds"))

# Calc weighted fraction
abund_comb <- lapply(abund_comb, function(x) {
  
  x[, tot_abund_site := sum(abundance), by = "site"]
  x[, `:=`(
    weighted_fraction = sum(abundance),
    weighted_fraction_rel = sum(abundance) / tot_abund_site
  ),
  by = c(
    "site",
    "group"
  )
  ]

  # TPGs as columns for BRT analyses
  x[, group := paste0("T", group)]

  unique(x[, .(
    Region,
    site,
    group,
    weighted_fraction_rel
  )])
})
# Use of relative fraction, slightly better
# performance than absolute (with family level trait data)
# abund_weighted_frac <- unique(abund[, .(
#   Region,
#   site,
#   group,
#   weighted_fraction
# )])

# All 15 groups observed in the sampling campaign but not in every region!
lapply(abund_comb, function(x) x[, uniqueN(group), by = "Region"])
# TPG11 and TPG14 do not occur in all regions
abund_comb$family_lvl[, .N, by = c("group")] %>% .[order(N), ]
abund_comb$family_lvl[group == "T11", Region] %>% unique()
abund_comb$family_lvl[group == "T14", Region] %>% unique()

# TPG7_genus for genus lvl trait data
abund_comb$genus_lvl[, .N, by = c("group")] %>% .[order(N), ]

# Trait groups
# Relative fraction
trait_groups_rel <- lapply(abund_comb, function(x) {
  split(x,
    f = x$Region
  ) %>%
    lapply(., function(y) {
      dcast(y[, .(site, group, weighted_fraction_rel)],
        site ~ group,
        value.var = "weighted_fraction_rel"
      )
    })
})
trait_groups_rel$genus_lvl$Midwest[abund[Region == "Midwest", ], STAID := i.STAID,
  on = "site"
]
trait_groups_rel$genus_lvl$Midwest[,STAID := sub("\"([0-9]{1,})\"", "\\1", STAID)]
trait_groups_rel$genus_lvl$Midwest[,STAID := paste0("T", STAID)]

trait_groups_rel$family_lvl$Midwest[abund[Region == "Midwest", ], STAID := i.STAID,
  on = "site"
]
trait_groups_rel$family_lvl$Midwest[,STAID := sub("\"([0-9]{1,})\"", "\\1", STAID)]
trait_groups_rel$family_lvl$Midwest[,STAID := paste0("T", STAID)]

# Number of TPGs represented at each region
lapply(trait_groups_rel$genus_lvl, function(x) dim(x[, .SD, .SDcols = patterns("T.*")]))
lapply(trait_groups_rel$family_lvl, function(x) dim(x[, .SD, .SDcols = patterns("T.*")]))
saveRDS(trait_groups_rel, file.path(path_cache, "trait_groups_rel.rds"))

# Taxa that fall in certain group (absolute)
# trait_groups <- split(abund_weighted_frac,
#   f = abund_weighted_frac$Region
# ) %>%
#   lapply(., function(y) {
#     dcast(y[, .(site, group, weighted_fraction)],
#       site ~ group,
#       value.var = "weighted_fraction"
#     )
#   })
# trait_groups$Midwest[abund[Region == "Midwest", ], STAID := i.STAID,
#   on = "site"
# ]
# trait_groups$Midwest[,STAID := sub("\"([0-9]{1,})\"", "\\1", STAID)]
# trait_groups$Midwest[,STAID := paste0("T", STAID)]
# saveRDS(trait_groups, file.path(path_cache, "trait_groups.rds"))
