# ________________________________________________
# Relationship between fraction of EPT taxa and toxicity 
# ________________________________________________

# Abundances
abund <- readRDS(file.path(path_cache, "total_abund_CEOPT.rds"))
abund <- abund[, .(
    Region,
    site, 
    STAID,
    taxon,
    species,
    genus,
    family,
    order,
    taxonomic_level,
    abundance)]

# Calc. abundance per site, but merge chemical information per STAID
# Midwest has more sites than STAID, e.g. T391114085205801
abund[Region == "Midwest", STAID := sub("\"([0-9]{1,})\"", "\\1", STAID)]
abund[Region == "Midwest", STAID := paste0("T", STAID)]

# Few sites have two IDs (for eco and wq)
# because chemical and eco information were taken from slightly
# different positions (mainly because of accessability issues)
abund[STAID == "T03611200", STAID := "T03611100"]
abund[site == "T12073525", site := "T12073425"]
# save this version of the abundance data with improved STAID labels 
# and the site IDs changed for later use
saveRDS(abund, file.path(path_cache, "total_abund_CEOPT_corrected.rds"))

# Toxicity
max_tu <- readRDS(file.path(path_cache, "max_tu.rds"))
setnames(max_tu, "TSITE_NO_WQ", "site")

# Calculate Fraction EPT of the whole assemblage per site
abund[, abund_order_site := sum(abundance), by = c("site", "order")]
abund[, total_abund := sum(abundance), by = "site"]
abund[order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera"),
    abund_EPT := sum(abundance),
    by = "site"
]
abund_subs <- unique(abund[
    !is.na(abund_EPT),
    .(Region, site, STAID, abund_EPT, total_abund)
])
abund_subs[, frac_EPT := abund_EPT/total_abund]
# Hmisc::describe(abund_subs$frac_EPT)

# Merge toxicity
abund_subs[max_tu[Region != "Midwest", ],
    max_log_tu := i.max_log_tu,
    on = "site"
]
abund_subs[max_tu[Region == "Midwest"],
    max_log_tu := i.max_log_tu,
    on = c("STAID" = "site")
]

# Few sites don't have chemical information
abund_subs <- abund_subs[!is.na(max_log_tu), ]
# Few STAIDs are duplicates (same chemical info, but different ecology)
# abund_subs[!is.na(STAID) & duplicated(STAID),]

# LM ----
# Significant in California and Northwest
# Negative relationship in California, Midwest, 
# Northwest, Southeast
# Positive relationship in Northeast (tough not significant) -> why?
# TODO: 
# Check assumptions of linear models!
EPT_regr <- abund_subs[, sregr_dt(
    x = .SD,
    form = "max_log_tu ~ frac_EPT"
), by = "Region"]
setnames(EPT_regr, "Pr(>|t|)", "p_value")
EPT_regr[p_value <= 0.05, ]

# Plotting ----
EPT_regr[, `:=`(
    coord_x = rep(0.75, 10),
    coord_y = rep(0.8, 10)
)]
EPT_regr[Region == "PN", Region := "Northwest"]
abund_subs[Region == "PN", Region := "Northwest"]

ggplot(abund_subs, aes(x = frac_EPT, y = max_log_tu)) +
    facet_wrap(as.factor(Region) ~.) +
    geom_point() +
    geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "steelblue"
    ) +
    geom_text(data = EPT_regr[id == "frac_EPT", ], aes(
        x = coord_x,
        y = coord_y,
        label = paste0(
            "slope = ", round(Estimate, digits = 2), "\n",
            "p = ", round(p_value, digits = 4)
        )
    )) +
    labs(
        x = "Fraction EPT taxa on the whole assemblage",
        y = "Max logTU"
    ) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text.x = element_text(
            family = "Roboto Mono",
            size = 14
        ),
        axis.text.y = element_text(
            family = "Roboto Mono",
            size = 14
        ),
        strip.text = element_text(
            family = "Roboto Mono",
            size = 14
        )
    )
ggsave(file.path(path_paper, "Graphs", paste0("EPT_toxicity.png")),
    width = 50,
    height = 30,
    units = "cm"
)
