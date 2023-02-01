# ________________________
# MaxTU calculations ----
# Also provide overview over toxicity data
# ________________________

## Add ecotox info ----
# Final list of ecotox_test data  
meta_rsqa_cmax <- readRDS(file.path(path_cache, "meta_rsqa_cmax.rds"))

# Add data from literature search (probably won't play a role) for important pc
ecotox_important_pc <- fread(file.path(path_repo, "pot_important_degradates_ecotox_info_missing.csv"))
ecotox_important_pc <- ecotox_important_pc[!is.na(LC50_ug_l), ]
ecotox_important_pc[Degradate == "Chlorpyrifos oxon", LC50_ug_l := gmean(LC50_ug_l)]
ecotox_important_pc <- ecotox_important_pc[!duplicated(Degradate), ]

# Add lc50 data obtained from literature 
meta_rsqa_cmax[ecotox_important_pc,
               `:=`(
                 lc50_ug_l = i.LC50_ug_l,
                 taxon_most_sensitive = i.most_sensitive_taxon,
                 source = i.Source,
                 ecotox_data_av = "Yes",
                 notes = i.Notes
               ),
               on = c(CASRN = "Cas_degradates")]
# saveRDS(meta_rsqa_cmax, file.path(path_cache, "meta_rsqa_cmax.rds"))

# Cmax data preprocessed
# Add lc50 of most sensitive organism
rsqa_cmax <- readRDS(file.path(path_cache, "rsqa_cmax_preproc.rds"))
rsqa_cmax[meta_rsqa_cmax[ecotox_data_av == "Yes", ],
  lc50_ug_l := i.lc50_ug_l,
  on = c(cas = "CASRN")
]
# These sites (PN and Midwest) have no chemical information
# rsqa_cmax[TSITE_NO_WQ %in% c("T12073525",
#                              "T03611200",
#                              "T03318800",
#                              "T393247089260701"),]

# Occr. pesticides overall
rsqa_cmax[cmax > 0, n_occr_pesticide := .N, by = "pesticide"]
rsqa_cmax[, n_occr_pesticide := n_occr_pesticide / n_sites]

## Calculate max TU ----
rsqa_cmax[cmax > 0 &
            !is.na(lc50_ug_l), log_tu := log(cmax / lc50_ug_l,
                                             base = 10),
          by = "TSITE_NO_WQ"]
rsqa_cmax[, max_log_tu := max(log_tu, na.rm = TRUE), by = "TSITE_NO_WQ"]
setnames(rsqa_cmax, "Region/Chemname", "Region")
# PCA_chem.R

# Unexpectedly low and high TU values
check_low_tu <- rsqa_cmax[max_log_tu <= -5, .(
  Region,
  TSITE_NO_WQ,
  cas,
  pesticide, 
  cmax,
  lc50_ug_l,
  log_tu,
  max_log_tu
)] %>%
  .[max_log_tu == log_tu, ] %>% 
  .[order(Region), ] 


# check_high_tu <- rsqa_cmax[max_log_tu > 1, .(
#   Region,
#   TSITE_NO_WQ,
#   cas,
#   pesticide,
#   cmax,
#   lc50_ug_l,
#   log_tu,
#   max_log_tu
# )] %>%
#   .[max_log_tu == log_tu, ] %>% 
#   .[order(Region), ]
# saveRDS(check_high_tu, file.path(path_cache, "high_tu.rds"))

# ## Check ecotox info ----

# # Full list
# meta_rsqa_cmax[CASRN %in% check_low_tu[lc50_ug_l >= 10000, cas], ]
# meta_rsqa_cmax[CASRN %in% check_high_tu[, cas], ]

# # merge ecotox info 
# check_low_tu[meta_rsqa_cmax, `:=`(
#   taxon_most_sensitive = i.taxon_most_sensitive,
#   source = i.source
# ),
# on = c("cas" = "CASRN")
# ]
# saveRDS(check_low_tu, file.path(path_cache, "low_tu.rds"))
# check_high_tu[meta_rsqa_cmax, `:=`(
#   taxon_most_sensitive = i.taxon_most_sensitive,
#   source = i.source
# ),
# on = c("cas" = "CASRN")
# ]
# saveRDS(check_high_tu, file.path(path_cache, "high_tu.rds"))
# fwrite(
#   check_low_tu[order(Region, cas), ],
#   file.path(path_out, "pesticides_tu_below_minus5.csv")
# )

# Check substances were lc is below 1
# meta_rsqa_cmax[lc50_ug_l <= 1, .(
#   CASRN,
#   Chemname,
#   Parent_Degradate,
#   Pesticide_use_group,
#   lc50_ug_l,
#   taxon_most_sensitive,
#   source
# )] %>% 
# saveRDS(., file.path(path_cache, "pesticides_lc_50_below_1.rds"))
# fwrite(., file.path(path_out, "pesticides_lc50_below_1_µg_L.csv"))


## standartox
# ecotox <- readRDS(file.path(path_cache, "ecotox.rds"))
# # ecotox_aggr <- ecotox$aggregated
# # ecotox_aggr[, tax_all := NULL]
# ecotox_subset <- ecotox$filtered
# ecotox_subset[concentration_unit == "ul/l",
#               `:=`(concentration = concentration *
#                      1000,
#                    concentration_unit = "ug/l")] # (10^-6 * 10^3)/10^3

# ecotox_subset <- ecotox_subset[concentration_unit %in% c("ppb", "ug/l"), ]
# ecotox_subset[
#   cas %in% check_low_tu[lc50_ug_l >= 10000, cas],
#   .(
#     cname, 
#     cas,
#     exposure,
#     endpoint,
#     concentration,
#     concentration_unit,
#     concentration_orig,
#     concentration_unit_orig,
#     tax_taxon
#   )
# ]

# ## Ceriodaphnia dubia has a super low lc50 value
# ## with a strange unit
# ecotox_subset[
#   cas %in% check_high_tu[, cas],
#   .(
#     cname, 
#     cas,
#     exposure,
#     endpoint,
#     concentration,
#     concentration_unit,
#     concentration_orig,
#     concentration_unit_orig,
#     tax_taxon
#   )
# ]

# # PPDB
# ppdb_queried <- fread(file.path(path_in, "Chemicals", "PPDB_queried.csv"))
# setnames(ppdb_queried, "Source/Quality", "taxon_most_sensitive")
# ppdb_queried[meta_rsqa_cmax, `:=`(Parent_Degradate = i.Parent_Degradate,
#                                Pesticide_use_group = i.Pesticide_use_group), 
#              on = c(cas = "CASRN")]
# ppdb_queried[, source := "PPDB"]
# ppdb_queried[cas %in% check_low_tu[lc50_ug_l >= 30000, cas], ]

# For 2 sites where only 2-Hydroxyatrazine has been detected we cannot calculate a log TU
# rsqa_cmax[cmax > 0 & !TSITE_NO_WQ %in% rsqa_cmax[max_log_tu == log_tu, TSITE_NO_WQ], ]
max_tu <- rsqa_cmax[max_log_tu == log_tu, ]
max_tu[, n_sites_tu := .N]
max_tu[, n_max_tu := .N, by = "pesticide"] 
max_tu[, n_max_tu := n_max_tu/n_sites_tu]
max_tu <- max_tu[, .(
  SHORT_NAME,
  TSITE_NO_WQ,
  Region,
  pesticide,
  cmax,
  cas,
  Parent_Degradate,
  n_occr_pesticide,
  max_log_tu,
  n_max_tu
)]

# Sites with no pesticides detected or no where no ecotox info could be assigned
# Assign -5, likewise assign -5 for all pesticides with max_log_tu < -5 
max_tu <- rbind(max_tu,
  unique(rsqa_cmax[is.infinite(max_log_tu), .(
    TSITE_NO_WQ,
    Region,
    max_log_tu
  )]),
  fill = TRUE
)
max_tu[is.infinite(max_log_tu), max_log_tu := -5]
max_tu[max_log_tu < -5, max_log_tu := -5]
# saveRDS(max_tu, file.path(path_cache, "max_tu.rds"))

# Table for Co-Authors ecotox data
# meta_rsqa_cmax[rsqa_cmax[!is.na(n_occr_pesticide), ], 
#                n_occr_pesticide := i.n_occr_pesticide, 
#                on = c(CASRN = "cas")]
# meta_rsqa_cmax[source == "PPDB", source := "ppdb"]
# meta_rsqa_cmax[max_tu, n_max_tu := i.n_max_tu, on = c(CASRN = "cas")]
# meta_rsqa_cmax[, .(
#   CASRN,
#   Chemname,
#   Parent_Degradate,
#   Pesticide_use_group,
#   occr_pesticides = round(n_occr_pesticide * 100, digits = 2),
#   lc50_ug_l = round(lc50_ug_l, digits = 2),
#   max_tu_at_sites = round(n_max_tu * 100, digits = 2),
#   taxon_most_sensitive,
#   ecotox_data_av,
#   source,
#   notes
# )] %>% 
#   .[order(-ecotox_data_av, -source), ] %>% 
#   fwrite(., file.path(path_repo, "overview_pesticides_ecotox_data.csv"))

# Plot
# Overview over log max TU
ggplot(max_tu, aes(x = Region, y = max_log_tu)) +
  ggdist::stat_halfeye() +
  geom_boxplot(
    width = .2,
    outlier.shape = NA
  )+
  labs(y = "Max Log TU") +
  # stat_summary(fun = "mean", color = "forestgreen") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    axis.text.y = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    legend.position = "none",
  )
ggsave(
  filename = file.path(path_paper,
                       "Graphs",
                       "max_tu_overview.png"),
  width = 50,
  height = 30,
  units = "cm"
)

# Graphical overview of sites and max TU
ggplot(max_tu, aes(x = max_log_tu)) +
  geom_density() +
  facet_wrap(. ~ Region) +
  theme_bw()
# Hmisc::describe(max_tu)
