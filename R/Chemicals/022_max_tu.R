# ________________________
# MaxTU calculations ----
# Also provide overview over toxicity data
# ________________________

## Add ecotox info ----
# Final list of ecotox_test data  
meta_rsqa_cmax <- readRDS(file.path(path_cache, "meta_rsqa_cmax.rds"))
# meta_rsqa_cmax[ecotox_data_av == "Yes", ]

# Cmax data preprocessed
# Add lc50 of most sensitive organism
rsqa_cmax <- readRDS(file.path(path_cache, "rsqa_cmax_preproc.rds"))
rsqa_cmax[meta_rsqa_cmax[ecotox_data_av == "Yes", ],
  lc50_ug_l := i.lc50_ug_l,
  on = c(cas = "CASRN")
]
# Change region name
setnames(rsqa_cmax, "Region/Chemname", "Region")

# Load degradates and parent compound info
# For the degradates that occur, subset to the 
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

# Consider all Degradates by assigning lc50 values from the parent
# Then correlate with the maxTU calculated without these degradates
meta_rsqa_cmax[parent_compound_info, parent_compound := `i.parent group full`, on = c(CASRN = "CASRN2")]

# Select degradates without lc50 values and merge via column parent_compound 
degradates_missing <- merge.data.table(meta_rsqa_cmax[Parent_Degradate == "Degradate" & is.na(lc50_ug_l), ], 
                                       meta_rsqa_cmax[Parent_Degradate == "Pesticide" & !is.na(lc50_ug_l), 
                                                      .(parent_compound, lc50_ug_l)], 
                                       by = "parent_compound",
                                       suffixes = c("", "_parent"))

meta_rsqa_cmax[degradates_missing, lc50_ug_l_parent := i.lc50_ug_l_parent, on = "CASRN"]
meta_rsqa_cmax[, lc50_ug_l_parent := coalesce(lc50_ug_l, lc50_ug_l_parent)]
rsqa_cmax[meta_rsqa_cmax, lc50_ug_l_parent := i.lc50_ug_l_parent, on = c(cas="CASRN")]

meta_rsqa_cmax[Parent_Degradate == "Degradate" & is.na(lc50_ug_l), .(Parent_Degradate, 
                                                                     parent_compound,
                                                                     Chemname,
                                                                     CASRN,
                                                                     lc50_ug_l, 
                                                                     lc50_ug_l_parent)]
# Test:
# meta_rsqa_cmax[CASRN == "196611-19-5", ]
# meta_rsqa_cmax[parent_compound=="Metolachlor",]  
# Critical degradates?
# rsqa_cmax[pesticide %in% c("Diazinon oxon", "Tebupirimfos oxon") & cmax > 0, TSITE_NO_WQ] 
meta_rsqa_cmax[!is.na(lc50_ug_l), ]
rsqa_cmax[!is.na(lc50_ug_l), uniqueN(pesticide)]

# Occr. pesticides overall
rsqa_cmax[cmax > 0, n_occr_pesticide := .N, by = "pesticide"]
rsqa_cmax[, n_occr_pesticide := n_occr_pesticide / n_sites]

## Calculate max & sum TU ----
rsqa_cmax[cmax > 0 &
            !is.na(lc50_ug_l), `:=`(
              log_tu = log(cmax / lc50_ug_l, base = 10),
              log_tu_parent = log(cmax / lc50_ug_l_parent, base = 10)
            ), by = "TSITE_NO_WQ"]

# Take the maximum
rsqa_cmax[, `:=`(
  max_log_tu = max(log_tu, na.rm = TRUE),
  max_log_tu_parent = max(log_tu_parent, na.rm = TRUE)
), by = "TSITE_NO_WQ"]

# Sum TU
rsqa_cmax[cmax > 0 &
            !is.na(lc50_ug_l), tu := cmax / lc50_ug_l, by = "TSITE_NO_WQ"]
rsqa_cmax[, sum_log_tu := log(sum(tu, na.rm = TRUE), base = 10), by = "TSITE_NO_WQ"]

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
  sum_log_tu,
  max_log_tu_parent,
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
# Correlation max and sum tu
max_tu[, cor(x=max_log_tu, y=sum_log_tu, method="pearson")]

# Correlation max log tu parent and max log tu parent
max_tu[, cor(x=max_log_tu, y=max_log_tu_parent, method="pearson")]

# Assign -5 for sites with very low estimated toxicity
max_tu[is.infinite(max_log_tu), max_log_tu := -5]
max_tu[max_log_tu < -5, max_log_tu := -5]
max_tu[is.infinite(sum_log_tu), sum_log_tu := -5]
max_tu[sum_log_tu < -5, sum_log_tu := -5]

# Rename regions
max_tu[, Region := fcase(
  Region == "CSQA", "California",
  Region == "MSQA", "Midwest",
  Region == "NESQA", "Northeast",
  Region == "PNSQA", "PN",
  Region == "SESQA", "Southeast"
)]
# saveRDS(max_tu, file.path(path_cache, "max_tu.rds"))

## SI Table ecotox data ----
meta_rsqa_cmax[rsqa_cmax[!is.na(n_occr_pesticide), ], 
               n_occr_pesticide := i.n_occr_pesticide, 
               on = c(CASRN = "cas")]
meta_rsqa_cmax[source == "PPDB", source := "ppdb"]
meta_rsqa_cmax[max_tu, n_max_tu := i.n_max_tu, on = c(CASRN = "cas")]
meta_rsqa_cmax[, Pesticide_use_group := fcase(
  Pesticide_use_group == "H", "Herbicide",
  Pesticide_use_group == "I", "Insecticide",
  Pesticide_use_group == "F", "Fungicide"
)] %>%
  .[, Pesticide_use_group := factor(Pesticide_use_group,
    levels = c("Insecticide", "Fungicide", "Herbicide")
  )] %>%
  .[source %in% c("Literature", "Syngenta/ECHA"), source := "Scientific literature"] %>%
  .[
    order(-ecotox_data_av, -source, Pesticide_use_group),
    .(
      CASRN,
      Chemname,
      Parent_Degradate,
      Pesticide_use_group,
      lc50_ug_l = round(lc50_ug_l, digits = 2),
      taxon_most_sensitive,
      source
    )
  ] %>%
  setnames(
    .,
    c(
      "Chemname",
      "Parent_Degradate",
      "Pesticide_use_group",
      "lc50_ug_l",
      "taxon_most_sensitive",
      "source"
    ),
    c(
      "Pesticide name",
      "Parent or Degradate",
      "Pesticide use group",
      "LC50 [µg/L]",
      "Most sensitive taxon",
      "Source"
    )
  ) %>%
  fwrite(., file.path(path_paper, "Tables", "overview_pesticides_ecotox_data.csv"), sep = ";")

# Dominating pesticides
max_tu <- readRDS(file.path(path_cache, "max_tu.rds"))
# max_tu[TSITE_NO_WQ %in% c(
#   "T06804000",
#   "T451244123050200",
#   "T452538122213700",
#   "T05458800",
#   "T410133082465301"
# ), ]

# Overview over log max TU
summary_max_tu <- max_tu[, lapply(summary_funs, function(f) f(max_log_tu)),
  by = "Region"
]
max_tu[, quantile(max_log_tu, seq(0, 1, 0.25)), by = "Region"]
# max_tu[max_log_tu==-5, .N, by='Region']

# Pesticides that contributed most often to maxTU 
# SI Table
max_tu[, .(cas, .N), by = c("pesticide", "Region")] %>%
  unique() %>%
  .[!is.na(pesticide), ] %>%
  dcast(., ... ~ Region, value.var = "N") %>%
  .[, c("California", "Midwest", "Northeast", "PN", "Southeast") := lapply(.SD, function(x) fifelse(is.na(x), 0, x)),
    .SDcols = c("California", "Midwest", "Northeast", "PN", "Southeast")
  ] %>%
  .[meta_rsqa_cmax, Pesticide_use_group := i.Pesticide_use_group,
    on = c("cas" = "CASRN")
  ] %>%
  .[, Pesticide_use_group := factor(Pesticide_use_group, levels = c("I", "F", "H"))] %>%
  .[, sum_occr := California + Midwest + Northeast + PN + Southeast] %>%
  .[order(Pesticide_use_group, -sum_occr), ] %>%
  .[, .(pesticide, Pesticide_use_group, California, Midwest, Northeast, PN, Southeast)] %>%
  .[, Pesticide_use_group := fcase(
    Pesticide_use_group == "I", "Insecticide",
    Pesticide_use_group == "F", "Fungicide",
    Pesticide_use_group == "H", "Herbicide"
  ), ] %>%
  setnames(
    .,
    c("pesticide", "Pesticide_use_group", "PN"),
    c("Pesticide", "Pesticide use group", "Northwest")
  ) %>%
  fwrite(., file.path(path_paper, "Tables", "N_most_toxic_pesticides.csv"),
    sep = ";"
  )

# Dominant pesticides: range of TU values
# Imidaclopird
max_tu[
  pesticide %in% c(
    "Imidacloprid",
    "Fipronil",
    "Fipronil sulfone",
    "Chlorpyrifos"
  ), range(max_log_tu),
  by = "pesticide"
]

# Pesticide with the hightest maxTU value
max_tu[max_log_tu >= 1.5, ]

# Atrazine: Occurrences as most toxic pesticide
max_tu[pesticide == "Atrazine", .N, by = "Region"]
Hmisc::describe(max_tu[pesticide == "Atrazine", .(max_log_tu)])

## Plotting ----
summary_max_tu[, x_coord := c(1,2,3,4,5)]
ggplot(max_tu, aes(x = Region, y = max_log_tu)) +
  ggdist::stat_halfeye(
    adjust = .6,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .2,
    outlier.shape = NA
  ) +
  geom_point(
    size = 2,
    alpha = .4,
    position = position_jitter(
      seed = 1, width = .05
    )
  ) +
  geom_text(
    data = summary_max_tu,
    aes(
      x = Region,
      y = median,
      label = round(median, digits = 2)
    ),
    size = 6,
    nudge_x = -.25
  ) +
  scale_x_discrete(labels = c(
    "California (n = 85)",
    "Midwest (n = 100)",
    "Northeast (n = 94)",
    "Northwest (n = 87)",
    "Southeast (n = 76)"
  )) +
  labs(x = "", y = "maxTU") +
  coord_flip() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 18
    ),
    axis.text.y = element_text(
      family = "Roboto Mono",
      size = 18
    ),
    legend.position = "none",
  )
ggsave(
  filename = file.path(path_paper,
                       "Figures",
                       "max_tu_overview_flipped.png"),
  width = 27,
  height = 20,
  units = "cm"
)
