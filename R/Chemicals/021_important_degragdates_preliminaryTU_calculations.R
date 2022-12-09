# _________________________________________________________________________________________________
# Identify important degradates ----
# For which substances should we try to obtain ecotox information from the literature 
# We use the maxTU of the most sensitive freshwater organism 
# _________________________________________________________________________________________________

## Load concentration and ecotox data ----
# Load and preprcoess Cmax data
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
setDT(rsqa_cmax)
# rsqa_cmax[duplicated(TSITE_NO_WQ), ]

# Meta data site concentrations
# (also contains lc50 values)
meta_rsqa_cmax <- readRDS(file.path(path_cache, "meta_rsqa_cmax.rds"))

# Few columns containing pesticide conc. are character which should be numeric
# Filter(is.character, rsqa_cmax)
pest_char <- c("cis-Cyhalothric acid", "MCPA", "Metribuzin DK")
rsqa_cmax[, (pest_char) := lapply(.SD, as.numeric),.SDcols = pest_char]
rsqa_cmax <- melt(
  rsqa_cmax,
  id.vars = c(
    "SHORT_NAME",
    "TSITE_NO_WQ",
    "TSITE_NO_ECO",
    "Region/Chemname"
  ),
  variable.name = "pesticide",
  value.name = "cmax"
)
rsqa_cmax[meta_rsqa_cmax, `:=`(cas = i.CASRN,
                               Parent_Degradate = i.Parent_Degradate),
          on = c(pesticide = "Chemname")]

# Concentrations are in ng/L, convert to uq/L
rsqa_cmax[, cmax := cmax/1000]

# Few NAs in cmax - why?
# rsqa_cmax[is.na(cmax), ]

# Some sites have no pesticides detected
# sites_zero_cmax <- rsqa_cmax[, sum(cmax, na.rm = TRUE), by = "TSITE_NO_WQ"] %>% 
#   .[V1 == 0, TSITE_NO_WQ]
# rsqa_cmax <- rsqa_cmax[!TSITE_NO_WQ %in% sites_zero_cmax, ]

# Add total nr of sites
rsqa_cmax[, n_sites := uniqueN(TSITE_NO_WQ)]
# saveRDS(rsqa_cmax, file.path(path_cache, "rsqa_cmax_preproc.rds"))

# Add and join ecotox information
lc50 <- readRDS(file.path(path_cache, "lc50.rds"))
rsqa_cmax[lc50, lc50_most_sensitive := conc_ug_l, on = "cas"]

## TU Calculation ----
# max_TU of the most sensitive test organism
# c_i/lc50_i,j

# Handle degradates that lack ecotox information:
# Create list of parent compounds to the degradates that lack ecotox information
# Check which parent compounds occur widely and have a high toxicity (logTU > -3)
# -> then degradates likely also play a role
# Load data file with parent compound information
parent_compound_info <-
  read_xlsx(
    file.path(
      path_in,
      "Chemicals",
      "s2437_rsqa_compound_info_Aug2020_forStefanKunz2.xlsx"
    ),
    sheet = 2,
    skip = 1
  )
setDT(parent_compound_info)
pot_important_degradates <- unique(rsqa_cmax[cmax > 0 & is.na(lc50_most_sensitive), .(cas, pesticide)])
setnames(pot_important_degradates,
         c("cas",
           "pesticide"),
         c("cas_degradates",
           "degradate"))
pot_important_degradates[parent_compound_info,
                         parent_compound := `i.parent group full`,
                         on = c(cas_degradates = "CASRN2")]

# Calculate TUs where possible
rsqa_cmax[cmax > 0 &
            !is.na(lc50_most_sensitive), log_tu := log(cmax / lc50_most_sensitive,
                                                   base = 10), 
          by = "TSITE_NO_WQ"] 

# Check occurrence per site where logTU > -3 
# Overall 37 compounds, however, 16 only occur on two or less sites
important_parent_compounds <-
  rsqa_cmax[log_tu > -3 &
              Parent_Degradate %in% c("Pesticide",
                                      "Pesticide/Degradate"),
            .(cas, uniqueN(TSITE_NO_WQ) /
                n_sites),
            by = "pesticide"] %>%
  .[order(-V2),] %>%
  unique(.)
setnames(important_parent_compounds, "V2", "occr_pc_TU_larger_minus3")
important_parent_compounds[rsqa_cmax, Parent_Degradate := i.Parent_Degradate, 
                           on = "cas"]

# Check degradates that occur at many sites
degradates_occr <-
  rsqa_cmax[Parent_Degradate %in% c("Degradate", "Pesticide/Degradate") &
              cmax > 0 &
              is.na(log_tu), .(cas, uniqueN(TSITE_NO_WQ) / n_sites), by = "pesticide"] %>%
  .[order(-V2), ] %>%
  unique(.)
setnames(degradates_occr, "V2", "occr_degradates")

# Add to degradates list with missing ecotox data 
# - logTU parent compound
# - occur parent compound at x sites
# - occur of degradate at x sites
pot_important_degradates[important_parent_compounds,
                         `:=`(occr_pc_TU_larger_minus3 = i.occr_pc_TU_larger_minus3,
                              cas_parent_compound = i.cas,
                              Parent_Degradate = i.Parent_Degradate),
                         on = c(parent_compound = "pesticide")]
pot_important_degradates <-
  pot_important_degradates[!is.na(occr_pc_TU_larger_minus3),] %>%
  .[order(-occr_pc_TU_larger_minus3),]

# Merge back occr of degradates for those with parent compounds that had high TU
pot_important_degradates[degradates_occr, 
                         occr_degradates := i.occr_degradates, 
                         on = c(cas_degradates = "cas")]
pot_important_degradates[, `:=`(
  occr_degradates = round(occr_degradates * 100, digits = 1),
  occr_pc_TU_larger_minus3 = round(occr_pc_TU_larger_minus3 * 100, digits = 1)
)]
pot_important_degradates[meta_rsqa_cmax,
                         Pesticide_use_group := i.Pesticide_use_group,
                         on = c(cas_parent_compound = "CASRN")]
pot_important_degradates[, Parent_Degradate := NULL]

# Check if one of these degradates occurs at a site where the parent does
# not occur
# -> all seem to occur on sites where the parent doesn't occur
rsqa_cmax[cmax > 0 & cas %in% pot_important_degradates$cas_parent_compound, ]
sites_imp_pc <- rsqa_cmax[log_tu > -3 &
                            cas %in% pot_important_degradates$cas_parent_compound, TSITE_NO_WQ]
rsqa_cmax[cas %in% pot_important_degradates$cas_degradates, ] %>% 
  .[!(TSITE_NO_WQ %in% sites_imp_pc), pesticide] %>% 
  unique(.)

# 
setnames(
  pot_important_degradates,
  names(pot_important_degradates),
  c(
    "Cas_degradates",
    "Degradate",
    "Parent_compound",
    "Occurrences parent compound with TU > -3 [% sites]",
    "Cas_parent_compound",
    "Occurrences degradates [% sites]",
    "Pesticide_use_group"
  )
)
setcolorder(
  pot_important_degradates,
  c(
    "Cas_degradates",
    "Degradate",
    "Parent_compound",
    "Cas_parent_compound"
  )
)
pot_important_degradates[, N_sites := uniqueN(rsqa_cmax$TSITE_NO_WQ)]
# List with potentially important degradtes
# We checked for Fipronil, Atrazine, Diazinon, 
# Chlorpyrifos, Diuron, Tebupirimfos, Bifenthrin degradates 
# (these parent compounds reach the maximum TU at some sites)
# fwrite(
#   pot_important_degradates,
#   file.path(
#     path_repo,
#     "pot_important_degradates_ecotox_info_missing.csv"
#   )
# )
# pot_important_degradates[`Occurrences degradates [% sites]` > 5, ]

# Calculation maximum TU values
rsqa_cmax[, max_log_tu := max(log_tu, na.rm = TRUE), by = "TSITE_NO_WQ"]

# Search for those parent compounds that have the maximum TU at certain sites
# Try to obtain ecotox information for the degradates of these parent compounds 
tox_parent_compounds <- rsqa_cmax[max_log_tu == log_tu,] %>%
  .[pesticide %in% pot_important_degradates$Parent_compound, .N / n_sites *
      100, by = "pesticide"] %>%
  unique

pot_important_degradates[Parent_compound %in% tox_parent_compounds$pesticide, ]




