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
# All available LC50 values for invertebrates are queried and then
# the lowest test value that can be found for a given substance is taken.

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

# Few taxa are actually marine and need to be filtered out
# marine_taxa <- wormsbynames(unique(ecotox_subset[, tax_taxon]))
# saveRDS(marine_taxa, file.path(path_cache, "marine_taxa.rds"))
marine_taxa <- readRDS(file.path(path_cache, "marine_taxa.rds"))
setDT(marine_taxa)

# 9 Taxa potentially only marine or brackish
saline_taxa <- marine_taxa[, .(
  scientificname,
  status,
  isMarine,
  isBrackish,
  isFreshwater,
  isTerrestrial,
  isExtinct
)] %>%
  .[!is.na(scientificname), ] %>%
  .[
    (isMarine == 1 | isBrackish == 1) & (isFreshwater == 0 | is.na(isFreshwater)),
    scientificname
  ]

# RM Hydropsyche from the list (belongs to freshwater)
saline_taxa <- saline_taxa[saline_taxa != "Hydropsyche"]
ecotox_subset <- ecotox_subset[!tax_taxon %in% saline_taxa, ]

## Manual corrections and unit conversions to standartox ----

# Few taxa that are not freshwater remain
# Apis mellifera
# Eisenia fetida
ecotox_subset <- ecotox_subset[!tax_taxon %in% c("Apis mellifera", "Eisenia fetida"), ]

# Concentration unit
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

# Some pesticides queried from EPA ecotox DB have a unit error
# Dimethoate (60-51-5) is in mg/L given in the original source
# This value also appears in the Standartox list, but also 1.29 ug/L
# TODO: ask Andi
ecotox_subset[
  cas == "60-51-5" & tax_taxon == "Chironomus dilutus",
  `:=`(
    concentration = 1290,
    concentration_unit_orig = "mg/L"
  )
]

# Methamidophos
ecotox_subset[
  cas == "10265-92-6" & tax_taxon == "Daphnia magna",
  `:=`(
    concentration = 34000,
    concentration_unit_orig = "mg/L"
  )
]

# Propiconazole has a strange study as source from the 1970s 
# Removed all tests from this study
ecotox_subset <- ecotox_subset[!(cas == "60207-90-1" & tax_taxon %in% c(
  "Baetis rhodani",
  "Gammarus lacustris",
  "Physa fontinalis",
  "Hydropsyche siltalai",
  "Heptagenia sulphurea"
)), ]

# After filtering 60 substances remain
ecotox_subset <- ecotox_subset[concentration_unit %in% c("ppb", "ug/l"), ]
# unique(ecotox_subset$cas)

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
# PPDB and EPA Registration documents 
# meta_rsqa_cmax[!CASRN %in% lc50$cas &
#   !(Pesticide_use_group == "H" & Parent_Degradate == "Degradate"), ]
# ecotox_missing <- meta_rsqa_cmax[!CASRN %in% lc50$cas,]
# saveRDS(ecotox_missing, file.path(path_cache, "ecotox_missing.rds"))

# Query PPDB ----
# According to the PPDB: 
# "Where data for for several species are available, data for the most sensitive is given"
ppdb_queried <- fread(file.path(path_in, "Chemicals", "PPDB_queried.csv"))
setnames(ppdb_queried, "Source/Quality", "taxon_most_sensitive")
ppdb_queried[meta_rsqa_cmax, `:=`(Parent_Degradate = i.Parent_Degradate,
                               Pesticide_use_group = i.Pesticide_use_group), 
             on = c(cas = "CASRN")]
ppdb_queried[, source := "PPDB"]
# 82 substances covered
lc50 <- rbind(lc50, ppdb_queried[, .(
  pesticide,
  cas,
  taxon_most_sensitive,
  conc_ug_l,
  Pesticide_use_group,
  Parent_Degradate,
  source
)])

# Remaining pesticides with no information ----
# Focus on insecticides and fungicides because these are relevant for aquatic insects
# Manually queried: 
# https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/aquatic-life-benchmarks-and-ecological-risk
# new_missing <- meta_rsqa_cmax[!CASRN %in% lc50$cas &
#                  !(Parent_Degradate ==  "Degradate" &
#                      Pesticide_use_group == "H"), ] %>%
#   fwrite(., file.path(path_repo, "ecotox_info_missing.csv"))
ecotox_info_missing <- fread(file.path(path_repo, "ecotox_info_missing.csv"))
ecotox_info_missing[Availability_aq_life_bm == "Yes", source := "USEPA_Reg_doc"]

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

# Assign ecotox information for fentin (triphenyl tin) 
# from fentin hydorixde (only salt of fentin registered for use in the US)
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

# Add data from literature search (probably won't play a role) for important pc
ecotox_important_pc <- fread(file.path(path_repo, "pot_important_degradates_ecotox_info_missing.csv"))
ecotox_important_pc <- ecotox_important_pc[!is.na(LC50_ug_l), ]
ecotox_important_pc[Degradate == "Chlorpyrifos oxon", LC50_ug_l := gmean(LC50_ug_l)]
ecotox_important_pc <- ecotox_important_pc[!duplicated(Degradate), ]
setnames(
  ecotox_important_pc,
  c(
    "LC50_ug_l",
    "most_sensitive_taxon",
    "Source",
    "Notes",
    "Cas_degradates",
    "Degradate"
  ),
  c(
    "lc50_ug_l",
    "taxon_most_sensitive",
    "source",
    "notes",
    "cas",
    "pesticide"
  )
)

# Add lc50 data obtained from literature
lc50 <- rbind(lc50, ecotox_important_pc[, .(
  cas,
  lc50_ug_l,
  taxon_most_sensitive,
  source,
  notes,
  cas,
  pesticide
)],
fill = TRUE
)
lc50[is.na(Parent_Degradate), Parent_Degradate := "Degradate"]
# saveRDS(lc50, file.path(path_cache, "lc50.rds"))

# Complete list of all pesticides + degradates detected
meta_rsqa_cmax[lc50, `:=`(
  lc50_ug_l = i.conc_ug_l,
  taxon_most_sensitive = i.taxon_most_sensitive,
  ecotox_data_av = "Yes",
  source = i.source,
  notes = i.notes
),
on = c("CASRN" = "cas")
]
meta_rsqa_cmax[is.na(ecotox_data_av), ecotox_data_av := "No"]
setcolorder(meta_rsqa_cmax, c("CASRN", "Chemname"))

# Complete list of pesticides + degradates with missing ecotox test info
# meta_rsqa_cmax[ecotox_data_av == "No", ] %>% 
#   .[order(-Parent_Degradate, Pesticide_use_group), 
#     .(CASRN, Chemname, Parent_Degradate, Pesticide_use_group, lc50_ug_l)] %>% 
#   fwrite(., file.path(path_repo, "compl_list_ecotox_info_missing.csv"))

# Save for overview on lc50 values
# meta_rsqa_cmax[!is.na(lc50_ug_l), ] %>%
#   .[order(lc50_ug_l), .(
#     CASRN,
#     Chemname,
#     Parent_Degradate,
#     Pesticide_use_group,
#     lc50_ug_l = round(lc50_ug_l, digits = 4),
#     taxon_most_sensitive,
#     source
#   )] %>%
#   fwrite(., file.path(path_out, "pesticides_lc50.csv"))

# Save for further use in R
saveRDS(meta_rsqa_cmax, file.path(path_cache, "meta_rsqa_cmax.rds"))
