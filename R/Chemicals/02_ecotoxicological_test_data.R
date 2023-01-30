# ____________________________________________
# Match ecotoxicological test information ----
# ____________________________________________

## Load data ----
# Pesticide cmax concentration for pesticides and degradates 
# Values below detection limit subtituted with 0
rsqa_cmax <-
  read_xlsx(
    file.path(
      path_in,
      "Chemicals",
      "s2437_rsqa_compound_info_Aug2020_forStefanKunz3.xlsx"
    ),
    sheet = 3,
    skip = 6
  )

# Meta information including CAS-RN
# Is included in the header of the table, nedds some preprocessing
meta_rsqa_cmax <- read_xlsx(
  file.path(
    path_in,
    "Chemicals",
    "s2437_rsqa_compound_info_Aug2020_forStefanKunz3.xlsx"
  ),
  sheet = 3,
  range = "D1:GC7"
)
setDT(meta_rsqa_cmax)

# - Preprocess meta_rsqa_cmax
# str(meta_rsqa_cmax)
meta_rsqa_cmax <- melt(meta_rsqa_cmax, id.vars = "...1") %>% 
  dcast(., variable ~ ..., value.var = "value")
meta_rsqa_cmax[, variable := NULL]
setnames(meta_rsqa_cmax, 
         c("Region/Chemname", 
           "Parent/Degradate",
           "Pesticide use group"), 
         c("Chemname",
           "Parent_Degradate", 
           "Pesticide_use_group"))

# 179 Pesticides + degradates
# length(unique(meta_rsqa_cmax$CASRN))

# Query standartox ----
# Most sensitive organism:
# It's not important to use LC50 data for the observed taxa.
# Rather all available LC50 values for invertebrates are queried and then 
# the lowest test value that can be found for a given substance is taken. 
# This then represents the most sensitive organism (i.e. most sensitive is for each substance).

# Stx query:
# Query for all substances including degradates
# In the end, we may exclude herbicide degradates
# ecotox <- stx_query(
#   cas = meta_rsqa_cmax$CASRN,
#   concentration_type = "active ingredient",
#   endpoint = "XX50",
#   ecotox_grp = "invertebrate",
#   habitat = "freshwater",
#   duration = c(48, 96),
#   taxa = NULL,
#   effect = "Mortality"
# )
# saveRDS(ecotox, file.path(data_cache, "ecotox.rds"))

# Initially, data for 86 pesticides could be obtained
ecotox <- readRDS(file.path(path_cache, "ecotox.rds"))
# ecotox_aggr <- ecotox$aggregated
# ecotox_aggr[, tax_all := NULL]
ecotox_subset <- ecotox$filtered
# Through Standartox test data for 84 substances could be obtained 
# length(unique(ecotox$filtered$cas))

## Concentration unit ---- 
# ppb; ug/l
# ppb is equivalent to ug/L
# Hmisc::describe(ecotox_subset[,.(concentration_unit)]) 
# g/m2?; mg/kg -> probably soil tox. tests
ecotox_subset[concentration_unit == "ul/l",
              `:=`(concentration = concentration *
                     1000,
                   concentration_unit = "ug/l")] # (10^-6 * 10^3)/10^3

# Value for Methomyl too high?
ecotox_subset[cas == "16752-77-5", concentration := 0.1]

# After filtering 83 substances remain 
ecotox_subset <- ecotox_subset[concentration_unit %in% c("ppb", "ug/l"), ]
# ecotox_subset$cas %>% unique()


## Most sensitive Taxon ---- 
# We identify now the taxon with the lowest LC50 for a given pesticide
# Then the multiple tests for one taxon are aggregated with the geometric mean.
# If we have multiple most sensitive taxa for a given pesticide, 
# the taxon with the lowest geometric mean LC50 is taken.
# unique(ecotox_subset[, .(tax_taxon, min(concentration)), by = "cname"]) 
lc50 <- ecotox_subset[,.(gmean = gmean(concentration), 
                         cas),
              by = c("cname", "tax_taxon")] %>%
  .[order(gmean), .SD, by = "cname"] %>% 
  .[, .SD[1, ], by = "cname"] %>% 
  .[order(cname), ]
setnames(
  lc50,
  c("cname",
    "tax_taxon",
    "gmean"),
  c("pesticide",
    "taxon_most_sensitive",
    "conc_ug_l")
)

# Merge back meta information
lc50[meta_rsqa_cmax,
     `:=`(Parent_Degradate = i.Parent_Degradate,
          Pesticide_use_group = i.Pesticide_use_group),
     on = c(cas = "CASRN")]
lc50[, source := "standartox"]

## Missing test information ----
# PPDB or Registration documents mentioned by Lisa Nowell
# meta_rsqa_cmax[!CASRN %in% lc50$cas & 
#                  !(Pesticide_use_group == "H" & Parent_Degradate == "Degradate"), ]
# ecotox_missing <- meta_rsqa_cmax[!CASRN %in% lc50$cas,]
# saveRDS(ecotox_missing, file.path(data_cache, "ecotox_missing.rds"))

# Query PPDB ----
# According to the PPDB: 
# "Where data for for several species are available, data for the most sensitive is given"
ppdb_queried <- fread(file.path(path_in, "Chemicals", "PPDB_queried.csv"))
setnames(ppdb_queried, "Source/Quality", "taxon_most_sensitive")
ppdb_queried[meta_rsqa_cmax, `:=`(Parent_Degradate = i.Parent_Degradate,
                               Pesticide_use_group = i.Pesticide_use_group), 
             on = c(cas = "CASRN")]
ppdb_queried[, source := "PPDB"]
# 103 substances covered, mostly Pesticides, 9 Degradates
lc50 <- rbind(lc50, ppdb_queried[, .(pesticide, 
                                     cas, 
                                     taxon_most_sensitive,
                                     conc_ug_l, 
                                     Pesticide_use_group,
                                     Parent_Degradate,
                                     source)])
# Hmisc::describe(lc50)
# saveRDS(lc50, file.path(data_cache, "lc50.rds"))

# Remaining pesticides with no information ----
# Focus on insecticides and fungicides because these are relevant 
# for aquatic insects
# Manually queried through 
# https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/aquatic-life-benchmarks-and-ecological-risk
# meta_rsqa_cmax[!CASRN %in% lc50$cas &
#                  !(Parent_Degradate ==  "Degradate" &
#                      Pesticide_use_group == "H"), ] %>%
#   fwrite(., file.path(path_repo, "ecotox_info_missing.csv"))
# TODO: Search for I and F Degradates -> check Registration documents
ecotox_info_missing <- fread(file.path(path_repo, "ecotox_info_missing.csv"))
ecotox_info_missing[Availability_aq_life_bm == "Yes", source := "USEPA_Reg_doc"]
#ecotox_info_missing[pesticide == "Fipronil amide", ]

# Final list of ecotox test values for now
lc50 <-
  rbind(lc50, ecotox_info_missing[Availability_aq_life_bm == "Yes", .(
    pesticide,
    cas,
    conc_ug_l,
    taxon_most_sensitive,
    Parent_Degradate,
    Pesticide_use_group,
    source
  )])

# Assign ecotox information for fentin (triphenyl tin) from fentin hydorixde (only salt
# of fentin registered for use in the US)
# Lowest value for Daphnia magna, 48h acute test
lc50 <- rbind(
  lc50,
  data.table(
    pesticide = "Fentin",
    taxon_most_sensitive = "Daphnia magna",
    conc_ug_l = 0.0165,
    cas = "668-34-8",
    Parent_Degradate = "Pesticide",
    Pesticide_use_group = "F",
    source = "PPDB"
  )
)
# saveRDS(lc50, file.path(path_cache, "lc50.rds"))

# Add infos for a few degradates from seven parent compounds that reach the maxTU
# Information obtained from the scientific literature

# Complete list of all pesticides + degradates detected
meta_rsqa_cmax[lc50, `:=`(
  lc50_ug_l = i.conc_ug_l,
  taxon_most_sensitive = i.taxon_most_sensitive,
  ecotox_data_av = "Yes",
  source = i.source
),
on = c(CASRN = "cas")]
meta_rsqa_cmax[is.na(ecotox_data_av), ecotox_data_av := "No"]
setcolorder(meta_rsqa_cmax, c("CASRN", "Chemname"))

# Complete list of pesticides + degradates with missing ecotox test info
# meta_rsqa_cmax[ecotox_data_av == "No", ] %>% 
#   .[order(-Parent_Degradate, Pesticide_use_group), 
#     .(CASRN, Chemname, Parent_Degradate, Pesticide_use_group, lc50_ug_l)] %>% 
#   fwrite(., file.path(path_repo, "compl_list_ecotox_info_missing.csv"))

# Save for overview table
saveRDS(meta_rsqa_cmax, file.path(path_cache, "meta_rsqa_cmax.rds"))