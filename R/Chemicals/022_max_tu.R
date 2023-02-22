# ________________________
# MaxTU calculations ----
# Also provide overview over toxicity data
# ________________________

## Add ecotox info ----
# Final list of ecotox_test data  
meta_rsqa_cmax <- readRDS(file.path(path_cache, "meta_rsqa_cmax.rds"))

# Cmax data preprocessed
# Add lc50 of most sensitive organism
rsqa_cmax <- readRDS(file.path(path_cache, "rsqa_cmax_preproc.rds"))
rsqa_cmax[meta_rsqa_cmax[ecotox_data_av == "Yes", ],
  lc50_ug_l := i.lc50_ug_l,
  on = c(cas = "CASRN")
]

# Occr. pesticides overall
rsqa_cmax[cmax > 0, n_occr_pesticide := .N, by = "pesticide"]
rsqa_cmax[, n_occr_pesticide := n_occr_pesticide / n_sites]

## Calculate max TU ----
rsqa_cmax[cmax > 0 &
  !is.na(lc50_ug_l), log_tu := log(cmax / lc50_ug_l,
  base = 10
),
by = "TSITE_NO_WQ"
]
rsqa_cmax[, max_log_tu := max(log_tu, na.rm = TRUE), by = "TSITE_NO_WQ"]
setnames(rsqa_cmax, "Region/Chemname", "Region")

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

# Sites with no pesticides detected or no where no ecotox info could be assigned:
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
