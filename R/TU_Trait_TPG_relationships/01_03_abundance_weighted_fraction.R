# ________________________________________________
# Abundance weighted fraction of TPGs ----
# Same families as in the CWM approach should be covered
# (California)
## _______________________________________________

# Abundance weighted fraction
trait_family <- readRDS(file.path(path_cache, "trait_family_tpg.rds"))
abund <- readRDS(file.path(path_cache, "total_abund_CEOPT.rds"))

# Merge TPGs
abund[trait_family, group := i.group, on = "family"]
# abund[, .N, by = c("Region", "group")] %>%
#   .[order(group), ]

# Families which are not classified to a certain TPG
# throw out for now
# unique(abund[is.na(group), .(family, order), by = "Region"])
# Also not covered in CWM approach?
abund <- abund[!is.na(group),]

# Calc weighted fraction
abund[, tot_abund_site := sum(abundance), by = "site"]
abund[, `:=`(
    weighted_fraction = sum(abundance),
    weighted_fraction_rel = sum(abundance) / tot_abund_site
),
by = c(
    "site",
    "group"
)
]

# TPGs as columns for BRT analyses
abund[, group := paste0("T", group)]
abund_weighted_frac <- unique(abund[, .(
  Region,
  site,
  group,
  weighted_fraction
)])
abund_weighted_frac_rel <- unique(abund[, .(
  Region,
  site,
  group,
  weighted_fraction_rel
)])

# Taxonomic composition
abund[, .(
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
  saveRDS(., file.path(path_cache, "tpg_taxonomic_composition.rds"))

# Trait groups
# Taxa that fall in certain group (absolute fraction)
trait_groups <- split(abund_weighted_frac,
  f = abund_weighted_frac$Region
) %>%
  lapply(., function(y) {
    dcast(y[, .(site, group, weighted_fraction)],
      site ~ group,
      value.var = "weighted_fraction"
    )
  })
trait_groups$Midwest[abund[Region == "Midwest", ], STAID := i.STAID,
  on = "site"
]
trait_groups$Midwest[,STAID := sub("\"([0-9]{1,})\"", "\\1", STAID)]
trait_groups$Midwest[,STAID := paste0("T", STAID)]
# saveRDS(trait_groups, file.path(path_cache, "trait_groups.rds"))

# Number of TPGs represented at each region
lapply(trait_groups, function(y) dim(y[, .SD, .SDcols = patterns("T.*")]))

# Relative fraction
trait_groups_rel <- split(abund_weighted_frac_rel,
  f = abund_weighted_frac_rel$Region
) %>%
  lapply(., function(y) {
    dcast(y[, .(site, group, weighted_fraction_rel)],
      site ~ group,
      value.var = "weighted_fraction_rel"
    )
  })
trait_groups_rel$Midwest[abund[Region == "Midwest", ], STAID := i.STAID,
  on = "site"
]
trait_groups_rel$Midwest[,STAID := sub("\"([0-9]{1,})\"", "\\1", STAID)]
trait_groups_rel$Midwest[,STAID := paste0("T", STAID)]
# saveRDS(trait_groups_rel, file.path(path_cache, "trait_groups_rel.rds"))