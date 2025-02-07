# _____________________________________________________________
# Env factors preprocessing
# Relevant factors: temp, nutrients, habitat quality
# TODO: Update paper results
# _____________________________________________________________

# Info from Ian:
#  full site list: more sites than I'm currently are using 
# (few have only water chemistry data - no Ecology).
#
# Merge info:
# TSTAID and use "Has.Invert = y".
#
# Columns L - AZ: summary variables for Habitat:
# Probably the best variables to describe Habitat quality are:
# "SubstrateD84.M" (size of Substrate at 84th point)
# "WetWidthIQR.STD"  (Standard deviation of channel wetted width),
# DMax42d.M (Max Depth over
# the previous 42 days based on stage records - but can be a lot of missing values), 
# other "Riffle.FRC" percent riffle habitat as a fraction.

# Columns BA - BT: basic water quality parameters:
# Temperature, ph, conductivity, etc.
# summarized as max., median, and min
# Ian suggests to the median values over the last 4 weeks.
# Two sets of variables that record temperature:
# - WQ parameters over the last 4 weeks,
# - Columns AV-AZ, based on continuous recorders
# (42 days prior to the ecological sampling)
# lost all the data from the California study, 
# best to use the ones mentioned above in the WQ variables.

# Columns BU - end: Summaries for Nutrients
# Again max, median, and min over the last 4 weeks,
# suggest using the median based variables.

# TODO:
# - Which variables to use for habitat and nutrients? -> Discuss with Ralf
# - Check in which regions most responsive traits/TPGs were found

# Env. parameter datasets
env_param <- fread(file.path(
    path_in,
    "NoA",
    "Env_parameters",
    "National Environ n Eco Metrics combined_v6_Germany.csv"
))
setnames(env_param, "TSTAID", "site")

# Subset to relevant columns
env_param <- env_param[, .(
    site,
    `Has Inverts`,
    SubstrateD84.M,
    WetWidthIQR.STD,
    DMax_42day.M,
    Riffle.FRC,
    Temp.median,
    NH3_4wk.median,
    orthoP_4wk.median,
    NO3NO2_4wk.median
)]

# Riffle frc for midwest needs to be divided by 100. There's one NA, which we omit
env_param[Riffle.FRC > 1, Riffle.FRC := Riffle.FRC / 100]

# CWM traits and abundance weighted fraction TPGs (family level)
data_cwm <- readRDS(file.path(path_cache, "data_cwm_final.rds"))
trait_groups_rel <- readRDS(file.path(path_cache, "trait_groups_rel_final.rds"))
trait_groups_rel_fam <- trait_groups_rel$family_lvl
trait_groups_rel_genus <- trait_groups_rel$genus_lvl

# Sites that contian biological information
sites_bio <- lapply(data_cwm, function(x) {
    if ("STAID" %in% names(x)) unique(x$STAID) else unique(x$site)
}) %>% Reduce(c, .)

# Subset to sites with biological info
env_param <- env_param[`Has Inverts` == "y" & site %in% sites_bio, ]

# Merge with biological data & toxicity data
# Responding CWMtraits: gatherer, predator, sens. organic
data_cwm_env <- lapply(data_cwm, function(x) {
    on <- if ("STAID" %in% names(x)) c("STAID" = "site") else "site"
    x[env_param, `:=`(
        SubstrateD84.M = i.SubstrateD84.M,
        WetWidthIQR.STD = i.WetWidthIQR.STD,
        DMax_42day.M = i.DMax_42day.M,
        Riffle.FRC = i.Riffle.FRC,
        Temp.median = i.Temp.median,
        NH3_4wk.median = i.NH3_4wk.median,
        orthoP_4wk.median = i.orthoP_4wk.median,
        NO3NO2_4wk.median = i.NO3NO2_4wk.median
    ), on = on]
})
saveRDS(data_cwm_env, file.path(path_cache, "data_cwm_env.rds"))

data_tpg_env_fam <- lapply(trait_groups_rel_fam, function(x) {
    on <- if ("STAID" %in% names(x)) c("STAID" = "site") else "site"
    x[env_param, `:=`(
        SubstrateD84.M = i.SubstrateD84.M,
        WetWidthIQR.STD = i.WetWidthIQR.STD,
        DMax_42day.M = i.DMax_42day.M,
        Riffle.FRC = i.Riffle.FRC,
        Temp.median = i.Temp.median,
        NH3_4wk.median = i.NH3_4wk.median,
        orthoP_4wk.median = i.orthoP_4wk.median,
        NO3NO2_4wk.median = i.NO3NO2_4wk.median
    ), on = on]
})

# subset to most responsive TPGs
# Family-level: TPG1_fam, TPG2_fam, TPG5_fam,
# TPG8_fam, TPG9_fam, TPG10_fam, and TPG12_fam
data_tpg_env_fam <- lapply(data_tpg_env_fam, function(x) {
    melt(x[, .SD,
        .SDcols = c(
            "site",
            "SubstrateD84.M",
            "WetWidthIQR.STD",
            "DMax_42day.M",
            "Riffle.FRC",
            "Temp.median",
            "NH3_4wk.median",
            "orthoP_4wk.median",
            "NO3NO2_4wk.median",
            "max_log_tu",
            "T1",
            "T2",
            "T5",
            "T8",
            "T9",
            "T10",
            "T12",
            "T15"
        )
    ], measure.vars = c(
        "T1",
        "T2",
        "T5",
        "T8",
        "T9",
        "T10",
        "T12",
        "T15"
    ),
    variable.name = "tpg")
})

# Add family label
lapply(data_tpg_env_fam, function(x) x[, tpg := paste0(sub("T", "TPG", tpg), "_fam")])
saveRDS(data_tpg_env_fam, file.path(path_cache, "data_tpg_env_fam.rds"))

# Genus lvl data
data_tpg_env_genus <- lapply(trait_groups_rel_genus, function(x) {
    on <- if ("STAID" %in% names(x)) c("STAID" = "site") else "site"
    x[env_param, `:=`(
        SubstrateD84.M = i.SubstrateD84.M,
        WetWidthIQR.STD = i.WetWidthIQR.STD,
        DMax_42day.M = i.DMax_42day.M,
        Riffle.FRC = i.Riffle.FRC,
        Temp.median = i.Temp.median,
        NH3_4wk.median = i.NH3_4wk.median,
        orthoP_4wk.median = i.orthoP_4wk.median,
        NO3NO2_4wk.median = i.NO3NO2_4wk.median
    ), on = on]
})

# subset to most responsive TPGs
# Genus-level: TPG4_genus, TPG10_genus, and TPG12_genus
data_tpg_env_genus <- lapply(data_tpg_env_genus, function(x) {
    melt(x[, .SD,
        .SDcols = c(
            "site",
            "SubstrateD84.M",
            "WetWidthIQR.STD",
            "DMax_42day.M",
            "Riffle.FRC",
            "Temp.median",
            "NH3_4wk.median",
            "orthoP_4wk.median",
            "NO3NO2_4wk.median",
            "max_log_tu",
            "T4",
            "T10",
            "T12"
        )
    ], measure.vars = c(
        "T4",
        "T10",
        "T12"
    ),
    variable.name = "tpg")
})

# Add genus label
lapply(data_tpg_env_genus, function(x) x[, tpg := paste0(sub("T", "TPG", tpg), "_genus")])
saveRDS(data_tpg_env_genus, file.path(path_cache, "data_tpg_env_genus.rds"))

## Variable selection ----
# Habitat vars
# SubstrateD84.M: 84th percentile of particles on bed surface in the reach [m]
# WetWidthIQR.STD: Standardized inner quarterile range of wetted widths: 
# (75th percentile of transect wetted widths - 25th percentile)/50th percentile
# DMax_42day.M: Maximum depth during the 42 days prior to ecological sampling [m]
# Riffle.FRC: Fraction of transects in the reach classified as riffle [0 to 1]
# 84 NAs for DMax_42day.M
env_param[, .(
    SubstrateD84.M,
    WetWidthIQR.STD,
    DMax_42day.M,
    Riffle.FRC
)]  %>% .[, lapply(.SD, function(x) sum(is.na(x)))]
Hmisc::describe(env_param[, .(
    SubstrateD84.M,
    WetWidthIQR.STD,
    DMax_42day.M,
    Riffle.FRC
)])

# Correlation with consistent CWM traits
# small size, gills, gatherer, filterer, predator, bi/multivoltinism, semivoltinism, swimming, and  S_org

pl_habitat <- list()
titles_cwm_traits <- c(
    "Size large",
    "Size small",
    "Gills",
    "Feed Gatherer",
    "Feed Filterer",
    "Feed Predator",
    "Bi/multivoltinism",
    "Semivoltinism",
    "Swimming",
    "Sensitivity Organic"
)
for (region in names(data_cwm_env)) {
    x <- data_cwm_env[[region]]
    pl_habitat[[region]] <- Map(
        function(i, title) {
            x[
                trait == i, ggpairs(.SD, title = title),
                .SDcols = c(
                    "SubstrateD84.M",
                    "WetWidthIQR.STD",
                    "Riffle.FRC",
                    "cwm_val"
                )
            ]
        },
        c(
            "size_large",
            "size_small",
            "resp_gil",
            "feed_gatherer",
            "feed_filter",
            "feed_predator",
            "volt_bi_multi",
            "volt_semi",
            "locom_swim",
            "sensitivity_organic"
        ),
        titles_cwm_traits
    )
}
pl_habitat$Southeast

# "Strongest" correlations TPGs (graphically):
pl_habitat_tpgs <- list()
titles_tpgs <- c(
        "TPG1_fam",
        "TPG2_fam",
        "TPG5_fam",
        "TPG8_fam",
        "TPG9_fam",
        "TPG10_fam",
        "TPG12_fam"
)
for (region in names(data_tpg_env_fam)) {
    x <- data_tpg_env_fam[[region]]
    pl_habitat_tpgs[[region]] <- Map(
        function(i, title) {
            x[
                tpg == i, ggpairs(.SD, title = title),
                .SDcols = c(
                    "SubstrateD84.M",
                    "WetWidthIQR.STD",
                    "Riffle.FRC",
                    "value"
                )
            ]
        },
        c(
            "TPG1_fam",
            "TPG2_fam",
            "TPG5_fam",
            "TPG8_fam",
            "TPG9_fam",
            "TPG10_fam",
            "TPG12_fam"
        ),
        titles_tpgs
    )
}
pl_habitat_tpgs$California

# According to correlation
# Riffle FRC seems to be most strongly correlated to TPGs
cor_habitat_tpgs <- list()
for (region in names(data_tpg_env_fam)) {
    cor_habitat_tpgs[[region]] <- data_tpg_env_fam[[region]][, lapply(.SD, function(x) {
        cor(x, value, use = "complete.obs")
    }),
    .SDcols = c(
        "SubstrateD84.M",
        "WetWidthIQR.STD",
        "Riffle.FRC"
    ), by = "tpg"
    ]
}
rbindlist(cor_habitat_tpgs, idcol = "region") %>%
    melt(., id.vars = c("region", "tpg")) %>%
    .[, max_val := max(value), by = c("region", "tpg")] %>% 
    .[value == max_val, ] %>% 
    .[order(region), ]

# Nutrients [mg/L]
# CWM Traits
pl_nutrients <- list()
for (region in names(data_cwm_env)) {
    x <- data_cwm_env[[region]]
    pl_nutrients[[region]] <- Map(
        function(i, title) {
            x[
                trait == i, ggpairs(.SD, title = title),
                .SDcols = c(
                    "NH3_4wk.median",
                    "orthoP_4wk.median",
                    "NO3NO2_4wk.median",
                    "cwm_val"
                )
            ]
        },
        c(
            "size_large",
            "size_small",
            "resp_gil",
            "feed_gatherer",
            "feed_filter",
            "feed_predator",
            "volt_bi_multi",
            "volt_semi",
            "locom_swim",
            "sensitivity_organic"
        ),
        titles_cwm_traits
    )
}
pl_nutrients$Southeast

# For TPGs
# most strongly/oftencorrelating: ortho P and NO3NO2
cor_nutrients_tpgs <- list()
for (region in names(data_tpg_env_fam)) {
    cor_nutrients_tpgs[[region]] <- data_tpg_env_fam[[region]][, lapply(.SD, function(x) {
        cor(x, value, use = "complete.obs")
    }),
    .SDcols = c(
        "NH3_4wk.median",
        "orthoP_4wk.median",
        "NO3NO2_4wk.median"
    ), by = "tpg"
    ]
}
rbindlist(cor_nutrients_tpgs, idcol = "region") %>%
    melt(., id.vars = c("region", "tpg")) %>%
    .[, max_val := max(value), by = c("region", "tpg")] %>% 
    .[value == max_val, ] %>% 
    .[order(region), ] %>% 
    .[, .N, by = "variable"]

# Correlations & Range between NH3 and P
# Take Phosphate, seems to correlate quite often with TPGs 
lapply(data_tpg_env_fam, function(x) x[, cor(NH3_4wk.median, orthoP_4wk.median)])
lapply(data_tpg_env_fam, function(x) {
    x[, .(
        r_NH3 = range(NH3_4wk.median),
        r_P = range(orthoP_4wk.median)
    )]
})

# EPT ----
ept  <- readRDS(file.path(path_cache, "ept.rds"))
ept[Region != "Midwest", STAID := site]

# Site assigned to STAID which was previously just used for Midwest
# to enable merge with env_param
ept[env_param, `:=`(
    Riffle.FRC = i.Riffle.FRC,
    Temp.median = i.Temp.median,
    orthoP_4wk.median = i.orthoP_4wk.median
),
on = c("STAID" = "site")
]
saveRDS(ept, file.path(path_cache, "ept_env.rds"))

# SPEAR ----
result_spear <- readRDS(file.path(path_cache, "spear_preprocessed.rds"))
result_spear$SPEAR_Midwest[data_cwm$Midwest, STAID := i.STAID, on = "site"]

# Merge with env data
result_spear  <- lapply(result_spear, function(x) {
    on <- if ("STAID" %in% names(x)) c("STAID" = "site") else "site"
    x[env_param, `:=`(
        Riffle.FRC = i.Riffle.FRC,
        Temp.median = i.Temp.median,
        orthoP_4wk.median = i.orthoP_4wk.median
    ),
    on = on
    ]
})

# Merge Toxicity
max_tu <- readRDS(file.path(path_cache, "max_tu.rds"))
setnames(max_tu, "TSITE_NO_WQ", "site")
lapply(result_spear, function(x) {
    on <- if ("STAID" %in% names(x)) c("STAID" = "site") else "site"
    x[max_tu, max_log_tu := i.max_log_tu, on = on]
})
saveRDS(result_spear, file.path(path_cache, "spear_env.rds"))