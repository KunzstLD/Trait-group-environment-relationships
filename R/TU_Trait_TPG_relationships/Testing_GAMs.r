# GAMs
# Intro: https://noamross.github.io/gams-in-r-course/chapter2
library(mgcv)
# TODO: test if GAMs yield better results than lm's
# Alternative: GLM with the right error distribution? Zero inflated?

# Most_important CWM traits 
data_cwm_final <- readRDS(file.path(path_cache, "data_cwm_final.rds"))
most_imp_cwm <- readRDS(file.path(path_cache, "most_important_traits_cwm.rds"))
# lapply(data_cwm_final, function(x) x[trait == "sensitivity_organic", range(cwm_val)])

# CWM ----
cwm_combined <- rbindlist(data_cwm_final, idcol = "region", fill = TRUE) %>%
    .[trait %in% "size_large", ] %>%
    .[, region := fifelse(region == "PN", "Northwest", region)]

gam_cwm <- list()
gam_res_trait <- list()
for (i in unique(cwm_combined$region)) {
    for (j in unique(cwm_combined$trait)) {
        gam_res_trait[[j]] <- gam(
            max_log_tu ~ s(cwm_val),
            data = cwm_combined[region == i & trait == j, ]
        )
    }
    gam_cwm[[i]]  <- gam_res_trait 
}
for (i in names(gam_cwm)) {
    png(file = file.path(path_out, "Graphs", paste0("gam_cwm_", i, ".png")))
    par(mfrow = c(1, 1))
    lapply(
        gam_cwm[[i]],
        function(x) plot(x, residuals = TRUE, pch = 1, 
        ylab = "max_log_TU", 
        main = "trait: large size")
    )
    dev.off()
}
for (i in names(gam_cwm)) {
    png(file = file.path(path_out, "Graphs", paste0("gam_cwm_check_", i, ".png")))
    par(mfrow = c(2, 2))
    lapply(
        gam_cwm[[i]],
        function(x) gam.check(x)
    )
    dev.off()
}

# Fit appropriate GAMs
# - Ausreißer?
# - Varianzhomogenität und Normalverteilung relevant wenn Gaussian Smoother verwendet
# - Schätzungen der Koeffizienten davon nicht beeinträchtigt, sondern SEs und entsprechend CIs.

# Normal GAM 
# test_gam <- gam(max_log_tu ~ s(cwm_val),
#     family = gaussian,
#     data = cwm_combined[region == "Northeast", ],
#     method = "REML"
# )
# summary(test_gam)
# gam.check(test_gam)
# plot(test_gam, residuals = TRUE, pch = 1, ylab = "s(max_log_TU)")
# coef(test_gam)

# Censored GAM
# Create a control object with desired settings
# test_dat <- cwm_combined[region == "Northeast", .(max_log_tu, cwm_val)]
# setDF(test_dat)
# tu_old <- test_dat$max_log_tu
# tu_new <- cbind(tu_old, tu_old)
# tu_new[tu_old<= -5, 2] <- -Inf
# test_dat$tu_new <- tu_new

# test_gam_c <- gam(tu_new ~ s(cwm_val),
#     family = cnorm,
#     data = test_dat,
#     method = "REML"
# )
# summary(test_gam_c)
# gam.check(test_gam_c)
# plot(test_gam_c, residuals = TRUE, pch = 1, ylab = "s(max_log_TU)")
# coef(test_gam_c)

# test_dat$resid <- resid(test_gam_c)
# test_dat$predictions <- predict(test_gam_c)

# ggplot(test_dat, aes(x = predictions, y = tu_new[, 1])) +
# geom_point()
# test_dat[test_dat$max_log_tu == -5, c("max_log_tu", "predictions", "resid")]



# Most important TPGs ----
trait_groups_rel_final <- readRDS(file.path(path_cache, "trait_groups_rel_final.rds"))
tpgs_names <- grep("T[0-9]{1,}" ,names(trait_groups_rel_final$family_lvl$Midwest), value = TRUE)
trait_groups_rel_final$family_lvl <- lapply(
    trait_groups_rel_final$family_lvl,
    function(x) {
        setnames(
            x,
            tpgs_names,
            paste0(tpgs_names, "_fam"),
            skip_absent = TRUE
        )
    }
)
tpg_comb_family <- rbindlist(trait_groups_rel_final$family_lvl, idcol = "region", fill = TRUE) %>%
    .[, .(region, T12_fam, T5_fam, max_log_tu)] %>%
    melt(.,
        id.vars = c("max_log_tu", "region"),
        variable.name = "tpg",
        value.name = "tpg_val"
    ) %>%
    .[, region := fifelse(region == "PN", "Northwest", region)]

gam_tpg <- list()
gam_res_tpg <- list()
for (i in unique(tpg_comb_family$region)) {
    for (j in unique(tpg_comb_family$tpg)) {
        gam_res_tpg[[j]] <- gam(
            max_log_tu ~ s(tpg_val),
            data = tpg_comb_family[region == i & tpg == j, ]
        )
    }
    gam_tpg[[i]]  <- gam_res_tpg 
}
for (i in names(gam_tpg)) {
    png(file = file.path(path_out, "Graphs", paste0("gam_tpg_", i, ".png")))
    par(mfrow = c(2, 1))
    lapply(
        gam_tpg[[i]],
        function(x) plot(x, residuals = TRUE, pch = 1, ylab = "max_log_TU")
    )
    dev.off()
}

# EPT ----
# Except for California, GAMs seem to fit a linear 
# relationship (default)
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
# saveRDS(abund, file.path(path_cache, "total_abund_CEOPT_corrected.rds"))

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

lm_ept <- list()
gam_ept <- list()
for (i in unique(abund_subs$Region)) {
    lm_ept[[i]] <- lm(max_log_tu ~ frac_EPT,
        data = abund_subs[Region == i, ]
    )

    gam_ept[[i]] <- gam(
        max_log_tu ~ s(frac_EPT),
        data = abund_subs[Region == i, ])
}
lapply(lm_ept, function(x) summary(x)[["r.squared"]]*100)
lapply(gam_ept, function(x) summary(x)) # EDFs often 1, i.e. linear fit!

png(file = file.path(path_out, "Graphs", paste0("gam_EPT.png")))
par(mfrow = c(2, 3))
lapply(gam_ept, function(x) plot(x, residuals = TRUE, pch = 1, ylab = "max_log_TU"))
dev.off()

# SPEAR ----
result_spear <- load_data(
    path = path_out,
    pattern = ".*exported.*.csv",
    format = "csv",
    name_rm_pattern = "_exported"
)
lapply(result_spear, function(x) {
    setnames(
        x, 
        c("Name_2", "Name_3"),
        c("Region", "site")
    )
})
lapply(result_spear, function(x) x[, Name_1 := NULL])

max_tu <- readRDS(file.path(path_cache, "max_tu.rds"))
setnames(max_tu, "TSITE_NO_WQ", "site")

result_spear$SPEAR_Midwest[abund, STAID := i.STAID, on = "site"]
lapply(result_spear, function(x) {
    on <- if ("STAID" %in% names(x)) c("STAID" = "site") else "site"
    x[max_tu, max_log_tu := i.max_log_tu, on = on]
})
result_spear <- lapply(result_spear, function(x) x[!is.na(max_log_tu), ])

gam_spear <- list()
for (i in names(result_spear)) {
    gam_spear[[i]] <- gam(
        max_log_tu ~ s(SPEAR_Pestizide),
        data = result_spear[[i]])
}
png(file = file.path(path_out, "Graphs", paste0("gam_SPEAR.png")))
par(mfrow = c(2, 3))
lapply(gam_spear, function(x) plot(x, residuals = TRUE, pch = 1, ylab = "max_log_TU"))
dev.off()

gam.check(gam_spear$SPEAR_California)

# zero-inflated model
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#zero-inflation-k-inflation-or-deficits
# fit seems to be okay, can we still improve it?
library("glmmTMB")

fit_zinormal <- glmmTMB(
    max_log_tu+5 ~ SPEAR_Pestizide,
    data = result_spear[["SPEAR_California"]],
    ziformula = ~1,
    family = gaussian
)
summary(fit_zinormal)

library("DHARMa")
sim_zinormal <- simulateResiduals(fittedModel = fit_zinormal)
plot(sim_zinormal)
testZeroInflation(sim_zinormal)

# Comparison to lm
fit_normal <- lm(max_log_tu + 5 ~ SPEAR_Pestizide,
    data = result_spear[["SPEAR_California"]]
)
sim_normal <- simulateResiduals(fittedModel = fit_normal)
plot(sim_normal)
testZeroInflation(sim_normal)

