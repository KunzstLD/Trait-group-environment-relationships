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

# CWM traits and abundance weighted fraction TPGs (family level)
data_cwm <- readRDS(file.path(path_cache, "data_cwm_final.rds"))
trait_groups_rel <- readRDS(file.path(path_cache, "trait_groups_rel_final.rds"))
trait_groups_rel <- trait_groups_rel$family_lvl

# Sites that contian biological information
sites_bio <- lapply(data_cwm, function(x) {
    if ("STAID" %in% names(x)) unique(x$STAID) else unique(x$site)
}) %>% Reduce(c, .)

# Subset to sites with biological info
env_param <- env_param[`Has Inverts` == "y" & site %in% sites_bio, ]

# Subset to relevant columns
env_param <- env_param[, .(
    site,
    SubstrateD84.M,
    WetWidthIQR.STD,
    DMax_42day.M,
    Riffle.FRC,
    Temp.median,
    NH3_4wk.median,
    orthoP_4wk.median,
    NO3NO2_4wk.median
)]

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
data_tpg_env <- lapply(trait_groups_rel, function(x) {
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
data_tpg_env <- lapply(data_tpg_env, function(x) {
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
            "T12",
            "T5",
            "T2",
            "T10",
            "T1",
            "T8"
        )
    ], measure.vars = c(
        "T12",
        "T5",
        "T2",
        "T10",
        "T1",
        "T8"
    ),
    variable.name = "tpg")
})

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
pl_habitat <- list()
for (region in names(data_cwm_env)) {
    x <- data_cwm_env[[region]]
    titles <- c("Feed Predator", "Feed Gatherer", "Sensitivity Organic")
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
        c("feed_predator", "feed_gatherer", "sensitivity_organic"),
        titles
    )
}
pl_habitat$Southeast
# "Strongest" correlations CWM:
# California: Riffle fraction highest corr.
# Midwest: Wetwidth (pred & gatherer), Riffle fraction (sens. organic)
# Northeast: SubstrateD84 (~ 0.1)
# Northwest: Wetwidth (predator, ~0.1), SubstrateD84 (gatherer, > 0.3), 
# RiffleFrc (sens. organic, > 0.3)
# Southeast: SubstrateD84 (predator & sens. organic, ~0.1), 
# RiffleFRC (gatherer, ~ 0.2)

# "Strongest" correlations TPGs (graphically):
#   T12, T5 in 4 regions
#   T2, T10, T1, T8 in 3 regions
pl_habitat_tpgs <- list()
for (region in names(data_cwm_env)) {
    x <- data_tpg_env[[region]]
    titles <- c("T12", "T5", "T2", "T10", "T1", "T8")
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
        c("T12", "T5", "T2", "T10", "T1", "T8"),
        titles
    )
}
pl_habitat_tpgs$California

# According to correlation
# Riffle FRC seems to be most strongly correlated to TPGs
cor_habitat_tpgs <- list()
for (region in names(data_tpg_env)) {
    cor_habitat_tpgs[[region]] <- data_tpg_env[[region]][, lapply(.SD, function(x) {
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
# TODO: NO3 + NO2?
pl_nutrients <- list()
for (region in names(data_cwm_env)) {
    x <- data_cwm_env[[region]]
    titles <- c("Feed Predator", "Feed Gatherer", "Sensitivity Organic")
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
        c("feed_predator", "feed_gatherer", "sensitivity_organic"),
        titles
    )
}
pl_nutrients$Southeast
# "Strongest" correlations:
# California: NO3NO2 (predator, sens.organic, 0.2 - 0.3),
# orthoP (gatherer, ~0.25)
# Midwest: NO3NO2 (predator, sens.organic, 0.1), 
# orthoP (gatherer, 0.1)
# Northeast: NO3NO2 (~0.2)
# Northwest: NH3 (predator, 0.1 - 0.33)
# Southeast: NO3NO2 (predator, ~0.2), 
# orthoP (gatherer, sens.organic, 0.1 - 0.2) 

# For TPGs
# most strongly correlating: ortho P and NO3NO2
cor_nutrients_tpgs <- list()
for (region in names(data_tpg_env)) {
    cor_nutrients_tpgs[[region]] <- data_tpg_env[[region]][, lapply(.SD, function(x) {
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


# Testing Interactions ----
interactions <- list(
    "max_log_tu * Riffle.FRC",
    "max_log_tu * Temp.median",
    "max_log_tu * NO3NO2_4wk.median",
    "max_log_tu * Riffle.FRC * Temp.median",
    "max_log_tu * Riffle.FRC * NO3NO2_4wk.median",
    "max_log_tu * Temp.median * NO3NO2_4wk.median"
)
trait <- list(
    "feed_gatherer",
    "feed_predator",
    "sensitivity_organic"
)
collect_form <- expand.grid(trait, "~", interactions)
titles_traits <- collect_form$Var1
collect_form <- as.data.table(collect_form) %>%
    .[, paste(Var1, Var2, Var3)] %>%
    as.list(.)

# Calc. interactions CWM
cwm_interactions <- lapply(data_cwm_env, function(dt) {
    dcast(dt,
        max_log_tu + Riffle.FRC + Temp.median + NO3NO2_4wk.median + site ~ trait,
        value.var = "cwm_val"
    ) %>%
        fun_interactions(
            x = .,
            formulas = collect_form,
            titles = titles_traits,
            naming_category = "cwm_trait"
        )
}) %>%
    rbindlist(., idcol = "region") %>%
    setnames(., "Pr(>|t|)", "p_value")

# Search for significant interactions with maxTU
cwm_interactions[p_value <= 0.05 & id != "(Intercept)", ] %>%
    .[id %like% ":", ] %>%
    .[id %like% "max\\_log\\_tu", ]

# Overview cwm traits interactions : 
# California: gatherer ~ toxicity * riffle
# Midwest: predator ~ toxicity * temp
# Northeast:
#   sens. organic ~ toxicity * riffle
#   sens. organic ~ toxicity * nutrients
#   sens. organic ~ toxicity * temp.
# Northwest:
#   predator ~ toxicity * riffle
#   gatherer ~ toxicity * nutrients
#   gatherer ~ toxicity * nutrients * temp
# Southeast: predator ~ toxicity * riffle

# Calc. interactions TPGs:
tpg <- list(
    "T12",
    "T5",
    "T2",
    "T10",
    "T1",
    "T8"
)
collect_form_tpg <- expand.grid(tpg, "~", interactions)
titles_tpg <- collect_form_tpg$Var1
collect_form_tpg <- as.data.table(collect_form_tpg) %>%
    .[, paste(Var1, Var2, Var3)] %>%
    as.list(.)

tpg_interactions <- lapply(data_tpg_env, function(dt) {
    dcast(dt,
        max_log_tu + Riffle.FRC + Temp.median + NO3NO2_4wk.median + site ~ tpg,
        value.var = "value"
    ) %>%
        fun_interactions(
            x = .,
            formulas = collect_form_tpg,
            titles = titles_tpg,
            naming_category = "tpg"
        )
}) %>%
    rbindlist(., idcol = "region") %>% 
    setnames("Pr(>|t|)", "p_value")

# Search for significant interactions with maxTU
tpg_interactions[p_value <= 0.05 & id != "(Intercept)", ] %>%
    .[id %like% ":", ] %>%
    .[id %like% "max\\_log\\_tu", ]
