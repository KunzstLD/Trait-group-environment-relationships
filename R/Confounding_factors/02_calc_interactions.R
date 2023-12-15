# _____________________________________________________________
# Calculating Interactions ----
# TODO: 
# -> need to decide how much of the results should be shown
# -> CWM interactions only in California for size large -> is size large relevant there?
# _____________________________________________________________

# Load data
data_cwm_env <- readRDS(file.path(path_cache, "data_cwm_env.rds"))
data_tpg_env_fam <- readRDS(file.path(path_cache, "data_tpg_env_fam.rds"))
data_tpg_env_genus <- readRDS(file.path(path_cache, "data_tpg_env_genus.rds"))

# Interactions fitted for cwm traits, tpgs, ept and spear 
interactions <- list(
    "max_log_tu * Riffle.FRC",
    "max_log_tu * Temp.median",
    "max_log_tu * orthoP_4wk.median",
    "max_log_tu * Riffle.FRC * Temp.median",
    "max_log_tu * Riffle.FRC * orthoP_4wk.median",
    "max_log_tu * Temp.median * orthoP_4wk.median"
)

# CWM ---- 
trait <- list(
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
)
collect_form_cwm <- expand.grid(trait, "~", interactions)
titles_traits <- collect_form_cwm$Var1
collect_form_cwm <- as.data.table(collect_form_cwm) %>%
    .[, paste(Var1, Var2, Var3)] %>%
    as.list(.)

# LF & interactions
data_cwm_env_lf <- lapply(data_cwm_env, function(dt) {
    dcast(dt,
        max_log_tu + Riffle.FRC + Temp.median + orthoP_4wk.median + site ~ trait,
        value.var = "cwm_val"
    )
})
cwm_interactions <- list()
for (i in names(data_cwm_env_lf)) {
    cwm_interactions[[i]] <- fun_interactions(
        x = data_cwm_env_lf[[i]],
        form = collect_form_cwm
    )
}
cwm_interactions <- rbindlist(cwm_interactions, idcol = "Region")
# saveRDS(cwm_interactions, file.path(path_cache, "cwm_interactions.rds"))

# TPG family level ----
tpg_fam <- list(
    "TPG1_fam",
    "TPG2_fam",
    "TPG5_fam",
    "TPG8_fam",
    "TPG9_fam",
    "TPG10_fam",
    "TPG12_fam"
)
collect_form_tpg <- expand.grid(tpg_fam, "~", interactions)
collect_form_tpg <- as.data.table(collect_form_tpg) %>%
    .[, paste(Var1, Var2, Var3)] %>%
    as.list(.)

# LF & interactions
data_tpg_env_fam_lf <- lapply(data_tpg_env_fam, function(dt) {
    dcast(dt,
        max_log_tu + Riffle.FRC + Temp.median + orthoP_4wk.median + site ~ tpg,
        value.var = "value"
    )
})
tpg_fam_interactions <- list()
for (i in names(data_tpg_env_fam_lf)) {
    tpg_fam_interactions[[i]] <- fun_interactions(
        x = data_tpg_env_fam_lf[[i]],
        form = collect_form_tpg
    )
}
tpg_fam_interactions <- rbindlist(tpg_fam_interactions, idcol = "Region")
# saveRDS(tpg_fam_interactions, file.path(path_cache, "tpg_fam_interactions.rds"))

# Search for significant interactions with maxTU
# Note for paper: T10 and T1 seem to correlate frequently with other variables
tpg_fam_interactions[fo %like% "TPG1_.+", ] %>%
    .[p.value <= 0.05 & term != "(Intercept)", ] %>%
    .[term %like% ":", ] %>% 
    .[term %like% "max_log_tu", ]

# TPG genus level ---- 
tpg_genus <- list(
    "TPG4_genus",
    "TPG10_genus",
    "TPG12_genus"
)
collect_form_tpg_genus <- expand.grid(tpg_genus, "~", interactions)
collect_form_tpg_genus <- as.data.table(collect_form_tpg_genus) %>%
    .[, paste(Var1, Var2, Var3)] %>%
    as.list(.)

# LF & interactions
data_tpg_env_genus_lf <- lapply(data_tpg_env_genus, function(dt) {
    dcast(dt,
        max_log_tu + Riffle.FRC + Temp.median + orthoP_4wk.median + site ~ tpg,
        value.var = "value"
    )
})
tpg_genus_interactions <- list()
for (i in names(data_tpg_env_genus_lf)) {
    tpg_genus_interactions[[i]] <- fun_interactions(
        x = data_tpg_env_genus_lf[[i]],
        form = collect_form_tpg_genus
    )
}
tpg_genus_interactions <- rbindlist(tpg_genus_interactions, idcol = "Region")
saveRDS(
    tpg_genus_interactions,
    file.path(path_cache, "tpg_genus_interactions.rds")
)

# Search for significant interactions with maxTU
tpg_genus_interactions[p.value <= 0.05 & term != "(Intercept)", ] %>%
    .[term %like% ":", ] %>%
    .[term %like% "max_log_tu", ]


# EPT ----
ept_env <- readRDS(file.path(path_cache, "ept_env.rds"))
collect_form_ept <- expand.grid("frac_EPT", "~", interactions)
collect_form_ept <- as.data.table(collect_form_ept) %>%
    .[, paste(Var1, Var2, Var3)] %>%
    as.list(.)

# Data combined, easy to calc interactions
# (could have used this approach for the other data as well?)
ept_interactions <- ept_env[, fun_interactions(.SD, form = collect_form_ept), by = "Region"]
saveRDS(
    ept_interactions,
    file.path(path_cache, "ept_interactions.rds")
)

# Search for significant interactions with maxTU
ept_interactions[p.value <= 0.05 & term != "(Intercept)", ] %>%
    .[term %like% ":", ] %>%
    .[term %like% "max_log_tu", ]

# SPEAR ----
spear_env <- readRDS(file.path(path_cache, "spear_env.rds"))
collect_form_spear <- expand.grid("SPEAR_Pestizide", "~", interactions)
collect_form_spear <- as.data.table(collect_form_spear) %>%
    .[, paste(Var1, Var2, Var3)] %>%
    as.list(.)

spear_interactions <- list()
for (i in names(spear_env)) {
    spear_interactions[[i]] <- fun_interactions(
        x = spear_env[[i]],
        form = collect_form_spear
    )
}
spear_interactions <- rbindlist(spear_interactions, idcol = "Region")
spear_interactions[, Region := sub("SPEAR_", "", Region)]
saveRDS(
    spear_interactions,
    file.path(path_cache, "spear_interactions.rds")
)

# Search for significant interactions with maxTU
# Only few two-way interactions with temp., riffle fraction, and phospate
spear_interactions[p.value <= 0.05 & term != "(Intercept)", ] %>%
    .[term %like% ":", ] %>%
    .[term %like% "max_log_tu", ]