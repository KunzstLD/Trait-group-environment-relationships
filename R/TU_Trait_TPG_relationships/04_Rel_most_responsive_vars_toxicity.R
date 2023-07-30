# ________________________________________________
# Post analysis
# Relationship of the most impportant traits and TPGs with max logTU
# Investigated with LM 

# TODO: 
# - Add other most consistent TPG_fam
# - Table with Coefficients for SI/What to report for GAMs (Brown Paper)?
# - Graph for EPT & SPEAR as well 
# ________________________________________________
data_cwm_final <- readRDS(file.path(path_cache, "data_cwm_final.rds"))
most_imp_cwm <- readRDS(file.path(path_cache, "most_important_traits_cwm.rds"))
most_imp_cwm[, .N, by = "trait_TPG"] %>% 
.[order(-N)]

# CWM
# size large: California, Midwest, Northwest
# traits that were most important in two regions might be added later
cwm_combined <- rbindlist(data_cwm_final, idcol = "region", fill = TRUE) %>%
    .[trait == "size_large", ] %>%
    .[, region := fifelse(region == "PN", "Northwest", region)]

# bivar. relationship
biv_plot <- list()
for(i in unique(cwm_combined$region)){
    biv_plot[[i]] <- ggpairs(cwm_combined[region == i, .(max_log_tu, cwm_val)],
    title = i)
}
biv_plot

# Most consistent responding traits
# Decide on a model
lm_cwm <- list()
sim_res_lm <- list()
gam_cwm <- list()
gam_assump <- list()
for (i in unique(cwm_combined$region)) {
    lm_cwm[[i]] <- lm(max_log_tu ~ cwm_val,
        data = cwm_combined[region == i, ]
    )

    sim_res_lm[[i]] <- simulateResiduals(fittedModel = lm_cwm[[i]])

    gam_cwm[[i]] <- gam(
        max_log_tu ~ s(cwm_val),
        data = cwm_combined[region == i, ],
        method = "REML"
    )

    gam_assump[[i]] <- gam.check(gam_cwm[[i]])
}

# Diagnostics
# GAM
# Suggests linear model for California, Midwest, and maybe for Southeast
# GAMs for the rest
for(i in names(gam_cwm)){
    plot(gam_cwm[[i]], main = i, residuals = TRUE, pch = 1, ylab = "max_log_TU")
}

# LM
# Midwest Model ok
# California shows some deviations
for(i in names(sim_res_lm)){
    plot(sim_res_lm[[i]], main = i)
}

# TPGs ----
## Family level ----
# "TPG1_fam", 
# "TPG2_fam", TPG5_fam,
#  TPG8_fam, TPG9_fam, TPG10_fam, and TPG12_fam
trait_groups_rel_final <- readRDS(file.path(path_cache, "trait_groups_rel_final.rds"))
tpgs_names <- grep("T[0-9]{1,}" ,names(trait_groups_rel_final$family_lvl$Midwest), value = TRUE)
trait_groups_rel_final$family_lvl <- lapply(
    trait_groups_rel_final$family_lvl,
    function(x) {
        setnames(
            x,
            tpgs_names,
            paste0("TPG", sub("T", "", tpgs_names), "_fam"),
            skip_absent = TRUE
        )
    }
)
tpg_comb_family <- rbindlist(trait_groups_rel_final$family_lvl, idcol = "region", fill = TRUE) %>%
    .[, .(
        region,
        TPG1_fam,
        TPG2_fam,
        TPG5_fam,
        TPG8_fam,
        TPG9_fam,
        TPG10_fam,
        TPG12_fam,
        max_log_tu
    )] %>%
    melt(.,
        id.vars = c("max_log_tu", "region"),
        variable.name = "tpg",
        value.name = "tpg_val"
    ) %>%
    .[, region := fifelse(region == "PN", "Northwest", region)]

biv_plot_tpg_fam <- list()
for (i in unique(tpg_comb_family$region)) {
    biv_plot_tpg_fam[[i]] <- ggpairs(
        tpg_comb_family[
            region == i & tpg == "TPG1_fam",
            .(max_log_tu, tpg_val)
        ],
        title = i
    )
}
biv_plot_tpg_fam

# Choose a model
lm_tpg_fam <- list()
sim_res_lm_tpg_fam <- list()
gam_tpg_fam <- list()
gam_assump_tpg_fam <- list()
for (i in unique(tpg_comb_family$region)) {
    for (j in unique(tpg_comb_family$tpg)) {
        lm_tpg_fam[[paste0(i, "_", j)]] <- lm(max_log_tu ~ tpg_val,
            data = tpg_comb_family[region == i & tpg == j, ]
        )

        sim_res_lm_tpg_fam[[paste0(i, "_", j)]] <- simulateResiduals(fittedModel = lm_tpg_fam[[paste0(i, "_", j)]])

        gam_tpg_fam[[paste0(i, "_", j)]] <- gam(
            max_log_tu ~ s(tpg_val),
            data = tpg_comb_family[region == i & tpg == j, ],
            method = "REML"
        )

        gam_assump_tpg_fam[[paste0(i, "_", j)]] <- gam.check(gam_tpg_fam[[paste0(i, "_", j)]])
    }
}

# Plots
# - California all LM except for TPG10_fam 
# - Midwest All LM except for TPG5_fam & TPG12_fam
# - Northeast all GAM except TPG1_fam & TPG9_fam
# - Northwest all GAM except TPG10_fam & TPG12_fam
# - Southeast all LM except TPG1_fam
# linear model California and Midwest, maybe Southeast as well
par(mfrow = c(5,7))
for(i in names(gam_tpg_fam)){
    plot(gam_tpg_fam[[i]], main = i, residuals = TRUE, pch = 1, ylab = "max_log_TU")
}

# LM Diagnostics
for(i in names(sim_res_lm_tpg_fam)){
    plot(sim_res_lm_tpg_fam[[i]], main = i)
}


## Genus level ----
# most consistent (3 regions) TPG4_genus, TPG10_genus, and TPG12_genus
trait_groups_rel_final$genus_lvl <- lapply(
    trait_groups_rel_final$genus_lvl,
    function(x) {
        setnames(
            x,
            tpgs_names,
            paste0(tpgs_names, "_genus"),
            skip_absent = TRUE
        )
    }
)
tpg_comb_genus <- rbindlist(trait_groups_rel_final$genus_lvl, idcol = "region", fill = TRUE) %>%
    .[, .(region, T4_genus, T10_genus, T12_genus, max_log_tu)] %>%
    melt(.,
        id.vars = c("max_log_tu", "region"),
        variable.name = "tpg",
        value.name = "tpg_val"
    ) %>%
    .[region == "PN", region := "Northwest"]

# Choose a model
lm_tpg_gen <- list()
sim_res_lm_tpg_gen <- list()
gam_tpg_gen <- list()
gam_assump_tpg_gen <- list()
for (i in unique(tpg_comb_genus$region)) {
    for (j in unique(tpg_comb_genus$tpg)) {
        lm_tpg_gen[[paste0(i, "_", j)]] <- lm(max_log_tu ~ tpg_val,
            data = tpg_comb_genus[region == i & tpg == j, ]
        )

        sim_res_lm_tpg_gen[[paste0(i, "_", j)]] <-
            simulateResiduals(fittedModel = lm_tpg_gen[[paste0(i, "_", j)]])

        tryCatch(
            {
                gam_tpg_gen[[paste0(i, "_", j)]] <- gam(
                    max_log_tu ~ s(tpg_val),
                    data = tpg_comb_genus[region == i &
                        tpg == j, ],
                    method = "REML"
                )
            },
            error = function(e) {
                cat("ERROR :", conditionMessage(e), "\n")
            }
        )

        tryCatch(
            {
                gam_assump_tpg_gen[[paste0(i, "_", j)]] <- gam.check(gam_tpg_gen[[paste0(i, "_", j)]])
            },
            error = function(e) {
                cat("ERROR :", conditionMessage(e), "\n")
            }
        )
    }
}

# Diagnostics
# GAM
# Suggests linear model California and Midwest, maybe Southeast as well
par(mfrow = c(5,3))
for(i in names(gam_tpg_gen)){
    plot(gam_tpg_gen[[i]], main = i, residuals = TRUE, pch = 1, ylab = "max_log_TU")
}

# LM
# Midwest Model ok
# California shows some deviations
for(i in names(sim_res_lm_tpg_gen)){
    plot(sim_res_lm_tpg_gen[[i]], main = i)
}


# Summary plots ----
# cwm_combined <- cwm_combined[, .(region, trait, cwm_val, max_log_tu)]
# saveRDS(cwm_combined, file.path(path_cache, "cwm_combined.rds"))
cwm_combined <- readRDS(file.path(path_cache, "cwm_combined.rds"))

# tpg_comb_family <- tpg_comb_family[, .(region, tpg, tpg_val, max_log_tu)]
# tpg_comb <- rbindlist(
#     list(
#         "tpg_fam" = tpg_comb_family,
#         "tpg_gen" = tpg_comb_genus
#     ),
#     use.names = TRUE,
#     idcol = "approach"
# )
# tpg_comb[, id := paste0(region, "_", tpg)]
# saveRDS(tpg_comb, file.path(path_cache, "tpg_combined.rds"))
tpg_comb <- readRDS(file.path(path_cache, "tpg_combined.rds"))
tpg_comb[tpg %like% "_genus", tpg := sub("T", "TPG", tpg)]

## CWM plot ----
cwm_plot <- ggplot(cwm_combined, aes(x = cwm_val, y = max_log_tu)) +
    facet_wrap(.~as.factor(region)) +
    geom_point(alpha = 0.5, shape = 1, size = 2) +
    geom_smooth(
        data = ~ .x[region %in% c("California", "Midwest", "Southeast"), ],
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "steelblue"
    ) +
    geom_smooth(
        data = ~ .x[region %in% c("Northeast", "Northwest"), ],
        method = "gam",
        formula = y ~ s(x),
        se = TRUE,
        color = "steelblue"
    ) +
    labs(x = "CWM values size large", 
         y = "Max logTU") +
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
ggsave(file.path(path_paper, "Graphs", "gam_lm_cwm_toxicity.png"),
    width = 35,
    height = 30,
    units = "cm"
)

## TPG plot ----
tpg_plot <- ggplot(tpg_comb, aes(x = tpg_val, y = max_log_tu)) +
    facet_grid(as.factor(tpg) ~ as.factor(region)) +
    geom_point(alpha = 0.5, shape = 1, size = 2) +
    geom_smooth(
        data = ~ .x[approach == "tpg_fam" &
            (
                (region == "California" & id != "California_TPG10_fam") |
                (region == "Midwest" & !id %in% c("Midwest_TPG5_fam", "Midwest_TPG12_fam")) |
                (region == "Southeast" & id != "Southeast_TPG1_fam") |
                (region == "Northeast" & id %in% c("Northeast_TPG1_fam", "Northeast_TPG9_fam")) |
                (region == "Northwest" & id %in% c("Northwest_TPG10_fam", "Northwest_TPG12_fam"))
            ), ],
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "gold"
    ) +
    geom_smooth(
        data = ~ .x[approach == "tpg_fam" &
            (
                (region == "California" & id == "California_TPG10_fam") |
                (region == "Midwest" & id %in% c("Midwest_TPG5_fam", "Midwest_TPG12_fam")) | 
                (region == "Southeast" & id == "Southeast_TPG1_fam") | 
                (region == "Northeast" & !id %in% c("Northeast_TPG1_fam", "Northeast_TPG9_fam"))|
                (region == "Northwest" & !id %in% c("Northwest_TPG10_fam", "Northwest_TPG12_fam"))
            ), ],
        method = "gam",
        formula = y ~ s(x),
        se = TRUE,
        color = "gold"
    ) +
    geom_smooth(
        data = ~ .x[(region %in% c("California", "Midwest", "Southeast") | id %in% c(
            "Northeast_T4_genus",
            "Northwest_T4_genus",
            "Northwest_T12_genus"
        )) & approach == "tpg_gen", ],
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "forestgreen"
    ) +
    geom_smooth(
        data = ~ .x[id %in% c(
            "Northeast_T10_genus",
            "Northeast_T12_genus",
            "Northwest_T10_genus"
        ) & approach == "tpg_gen", ],
        method = "gam",
        formula = y ~ s(x),
        se = TRUE,
        color = "forestgreen"
    ) +
   scale_x_continuous(limits = c(0, 0.35)) +
    labs(x = "Abundance weighted-fraction TPG", 
         y = "Max logTU") +
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
ggsave(file.path(path_paper, "Graphs", "gam_lm_tpg_toxicity.png"),
    width = 35,
    height = 40,
    units = "cm"
)

# How many values are above 0.3 for TPG values?
tpg_comb[, .(.N, tpg_val, tpg), by = "region"] %>% 
.[tpg_val > 0.3, (.N/N)*100, by = "tpg"] %>% 
unique()


# Summary table ----
## CWM table ----
cwm_table <- rbind(
    lapply(
        lm_cwm[c("California", "Midwest", "Southeast")],
        function(x) {
            broom::tidy(x)
        }
    ) %>%
        rbindlist(., idcol = "Region") %>% 
        setnames(., "statistic", "t"),
    rbind(
        lapply(
            gam_cwm[c("Northeast", "Northwest")],
            function(x) {
                broom::tidy(x, parametric = TRUE)
            }
        ) %>%
            rbindlist(., idcol = "Region")%>% 
            setnames(., "statistic", "t"),
        lapply(
            gam_cwm[c("Northeast", "Northwest")],
            function(x) {
                broom::tidy(x, parametric = FALSE)
            }
        ) %>%
            rbindlist(., idcol = "Region") %>%
            setnames(., 
            c("statistic", "p.value"), 
            c("F", "p.value_gam")),
        fill = TRUE
    ),
    fill = TRUE
)
cwm_table[term == "cwm_val", term := "cwm_val_size_large"]
cwm_table[term == "s(cwm_val)", term := "s(cwm_val_size_large)"]
# saveRDS(cwm_table, file.path(path_cache, "cwm_tabl_publ.rds"))

## TPG table ----
tpg_table_fam <- rbind(
    lapply(
        lm_tpg_fam[
            !names(lm_tpg_fam) %in% c(
                "California_TPG10_fam",
                "Midwest_TPG5_fam",
                "Midwest_TPG12_fam",
                "Southeast_TPG1_fam",
                "Northeast_TPG2_fam",
                "Northeast_TPG5_fam",
                "Northeast_TPG8_fam",
                "Northeast_TPG10_fam",
                "Northeast_TPG12_fam",
                "Northwest_TPG1_fam",
                "Northwest_TPG2_fam",
                "Northwest_TPG5_fam",
                "Northwest_TPG8_fam",
                "Northwest_TPG9_fam"
            )
        ],
        function(x) {
            broom::tidy(x)
        }
    ) %>%
        rbindlist(., idcol = "Region_TPG") %>% 
        setnames(., "statistic", "t"),
    rbind(
        lapply(
            gam_tpg_fam[
                names(lm_tpg_fam) %in% c(
                    "California_TPG10_fam",
                    "Midwest_TPG5_fam",
                    "Midwest_TPG12_fam",
                    "Southeast_TPG1_fam",
                    "Northeast_TPG2_fam",
                    "Northeast_TPG5_fam",
                    "Northeast_TPG8_fam",
                    "Northeast_TPG10_fam",
                    "Northeast_TPG12_fam",
                    "Northwest_TPG1_fam",
                    "Northwest_TPG2_fam",
                    "Northwest_TPG5_fam",
                    "Northwest_TPG8_fam",
                    "Northwest_TPG9_fam"
                )
            ], function(x) {
                broom::tidy(x, parametric = TRUE)
            }
        ) %>%
            rbindlist(., idcol = "Region_TPG") %>% 
            setnames(., "statistic", "t"),
        lapply(
            gam_tpg_fam[
                names(lm_tpg_fam) %in% c(
                    "California_TPG10_fam",
                    "Midwest_TPG5_fam",
                    "Midwest_TPG12_fam",
                    "Southeast_TPG1_fam",
                    "Northeast_TPG2_fam",
                    "Northeast_TPG5_fam",
                    "Northeast_TPG8_fam",
                    "Northeast_TPG10_fam",
                    "Northeast_TPG12_fam",
                    "Northwest_TPG1_fam",
                    "Northwest_TPG2_fam",
                    "Northwest_TPG5_fam",
                    "Northwest_TPG8_fam",
                    "Northwest_TPG9_fam"
                )
            ], function(x) {
                broom::tidy(x, parametric = FALSE)
            }
        ) %>%
            rbindlist(., idcol = "Region_TPG") %>% 
            setnames(., 
            c("statistic", "p.value"), 
            c("F", "p.value_gam")),
        fill = TRUE
    ) %>% 
    .[order(Region_TPG), ],
    fill = TRUE
) 
tpg_table_fam[(p.value > 0.05 | p.value_gam > 0.05) & term != "(Intercept)"]
tpg_table_fam[, `:=`(
    estimate = round(estimate, digits = 2),
    std.error = round(std.error, digits = 2),
    t = round(t, digits = 2),
    F = round(F, digits = 2),
    p.value = fifelse(
        round(p.value, digits = 3) == 0,
        "<0.01",
        paste(round(p.value, digits = 3))
    ),
    p.value_gam = fifelse(
        round(p.value_gam, digits = 3) == 0,
        "<0.01",
        paste(round(p.value_gam, digits = 3))
    )
)]
# saveRDS(tpg_table_fam, file.path(path_cache, "tpg_tabl_publ.rds"))

tpg_table_genus <- rbind(
    lapply(
        lm_tpg_gen[
            names(lm_tpg_gen) %in% c(
                "Northeast_T4_genus",
                "Northwest_T4_genus",
                "Northwest_T12_genus"
            ) |
                names(lm_tpg_gen) %like% "California.*|Midwest.*|Southeast.*"
        ],
        function(x) {
            broom::tidy(x)
        }
    ) %>%
        rbindlist(., "Region_TPG",
            use.names = TRUE,
            fill = TRUE
        ) %>%
        setnames(., "statistic", "t"),
    rbind(
        lapply(
            gam_tpg_gen[
                names(gam_tpg_gen) %in% c(
                    "Northeast_T10_genus",
                    "Northeast_T12_genus",
                    "Northwest_T10_genus"
                )
            ],
            function(x) {
                broom::tidy(x, parametric = TRUE)
            }
        ) %>%
            rbindlist(., "Region_TPG",
                use.names = TRUE,
                fill = TRUE
            ) %>%
            setnames(
                .,
                c("statistic", "p.value"),
                c("F", "p.value_gam")
            ),
        lapply(
            gam_tpg_gen[
                names(gam_tpg_gen) %in% c(
                    "Northeast_T10_genus",
                    "Northeast_T12_genus",
                    "Northwest_T10_genus"
                )
            ],
            function(x) {
                broom::tidy(x, parametric = FALSE)
            }
        ) %>%
            rbindlist(., "Region_TPG",
                use.names = TRUE,
                fill = TRUE
            ) %>%
            setnames(
                .,
                c("statistic", "p.value"),
                c("F", "p.value_gam")
            ),
        fill = TRUE
    ) %>%  
    .[order(Region_TPG), ],
    fill = TRUE
)
# tpg_table_genus[(p.value > 0.05 | p.value_gam > 0.05) & term != "(Intercept)"]
tpg_table_genus[, `:=`(
    estimate = round(estimate, digits = 2),
    std.error = round(std.error, digits = 2),
    t = round(t, digits = 2),
    F = round(F, digits = 2),
    p.value = fifelse(
        round(p.value, digits = 3) == 0,
        "<0.01",
        paste(round(p.value, digits = 3))
    ),
    p.value_gam = fifelse(
        round(p.value_gam, digits = 3) == 0,
        "<0.01",
        paste(round(p.value_gam, digits = 3))
    )
)]
# saveRDS(tpg_table_genus, file.path(path_cache, "tpg_tabl_genus_publ.rds"))

# Tables for publication
cwm_table <- readRDS(file.path(path_cache, "cwm_tabl_publ.rds"))
cwm_table[, .(
    Region,
    term,
    estimate = round(estimate, digits = 2),
    std.error = round(std.error, digits = 2),
    t = round(t, digits = 2),
    p.value = round(p.value, digits = 3),
    edf = round(edf, digits = 2),
    F = round(F, digits = 2),
    p.value_gam = round(p.value_gam, digits = 3)
)] %>%
    .[, `:=`(
        p.value = fifelse(p.value == 0, "< 0.001", paste(p.value)),
        p.value_gam = fifelse(p.value_gam == 0, "< 0.001", paste(p.value_gam))
    )] %>%
    .[order(Region), ] %>%  
fwrite(., file.path(path_paper, "Tables", "CWM_log_tu_results.csv"))

tpg_table_fam <- readRDS(file.path(path_cache, "tpg_tabl_publ.rds"))
tpg_table_fam[, Region := sub("^(.+)(_)(.+)(_)(.+)", "\\1", Region_TPG)]
tpg_table_fam[, TPG := sub("^(.+)(_)(.+)(_)(.+)", "\\3\\4\\5", Region_TPG)]

tpg_table_genus <- readRDS(file.path(path_cache, "tpg_tabl_genus_publ.rds"))
tpg_table_genus[, Region := sub("^(.+)(_)(.+)(_)(.+)", "\\1", Region_TPG)]
tpg_table_genus[, TPG := sub("^(.+)(_)(.+)(_)(.+)", "\\3\\4\\5", Region_TPG)]

tpg_table_final <- rbind(
    tpg_table_fam[, .(
        Region,
        TPG,
        term,
        estimate,
        std.error,
        t,
        p.value,
        edf,
        F,
        p.value_gam
    )],
    tpg_table_genus[, .(
        Region,
        TPG,
        term,
        estimate,
        std.error,
        t,
        p.value,
        edf,
        F,
        p.value_gam
    )]
) 
tpg_table_final[, edf := round(edf, digits = 2)] %>%
    fwrite(., file.path(path_paper, "Tables", "TPG_log_tu_results.csv"))


# Mclust approach
# TODO 
# - fix variance
# - do this for every region
# - do wilcox.test subseqnetly
# Load necessary package
# library(mclust)

# cwm_combined[region == "Northwest", {
#     mclust_fit <- Mclust(max_log_tu, G = 2)
#     .(
#         bimod_mean = mclust_fit$parameters$mean,
#         bimod_var = mclust_fit$parameters$variance
#     )
# }]
# # MWU/wilcox test
# lapply(cwm_combined[region == "Northeast", ])
# wilcox.test(
#     x = x[max_log_tu == -5, cwm_val],
#     y = x[max_log_tu != 5, cwm_val]
# )
# cwm_combined[, cens := fifelse(max_log_tu == -5, 1, 0)]
# cwm_combined[region == "Southeast",
#     wilcox.test(
#         x = .SD[cens == 1, cwm_val],
#         y = .SD[cens == 0, cwm_val]
#     )]