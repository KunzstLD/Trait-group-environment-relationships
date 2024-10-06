# ________________________________________________
# Post analysis
# Relationship of the most impportant traits and TPGs with max logTU
# Investigated with LM 

# TODO: 
# - Add other most consistent TPG_fam
# - Table with Coefficients for SI/What to report for GAMs (Brown Paper)?
# - Graph for EPT & SPEAR as well 
# - Calculate for all TPGs?
# ________________________________________________
data_cwm_final <- readRDS(file.path(path_cache, "data_cwm_final.rds"))
most_imp_cwm <- readRDS(file.path(path_cache, "most_important_traits_cwm.rds"))
most_imp_cwm[, .N, by = "trait_TPG"] %>% 
.[order(-N)]

# CWM ----
# size large: California, Midwest, Northwest
# Other traits: most important traits per region
cwm_combined <- rbindlist(data_cwm_final, idcol = "region", fill = TRUE) %>%
    .[trait %in% c("size_large", "feed_predator", "sensitivity_organic", "size_small", "feed_parasite"), ] %>%
    .[, region := fifelse(region == "PN", "Northwest", region)]

# bivar. relationship
# TODO: Plots need be saved correctly
biv_plot <- list()
for(i in unique(cwm_combined$region)){
  for(j in unique(cwm_combined$trait)){
    biv_plot[[paste(i, "_", j)]] <- ggpairs(cwm_combined[region == i & trait==j, .(max_log_tu, cwm_val)],
                                            title = paste(i, "_", j)) 
  }
}
biv_plot


# Most consistent responding traits
# Decide on a model
lm_cwm <- list()
sim_res_lm <- list()
gam_cwm <- list()
gam_assump <- list()
for (i in unique(cwm_combined$region)) {
  for(j in unique(cwm_combined$trait)){
    lm_cwm[[paste0(i, "_", j)]] <- lm(max_log_tu ~ cwm_val,
                                      data = cwm_combined[region == i & trait==j, ]
    )
    
    sim_res_lm[[paste0(i, "_", j)]] <- simulateResiduals(fittedModel = lm_cwm[[paste0(i, "_", j)]])
    
    gam_cwm[[paste0(i, "_", j)]] <- gam(
      max_log_tu ~ s(cwm_val),
      data = cwm_combined[region == i & trait==j, ],
      method = "REML"
    )
    
    gam_assump[[paste0(i, "_", j)]] <- gam.check(gam_cwm[[paste0(i, "_", j)]]) 
  }
}

# Diagnostics GAMs, visual inspection
par(mfrow = c(5,5))
for(i in names(gam_cwm)){
    plot(gam_cwm[[i]], main = i, residuals = TRUE, pch = 1, ylab = "max_log_TU")
}

# EDF for final decision between LM or GAM 
edf_cwm <- list()
for(i in names(gam_cwm)){
  edf_ <- broom::tidy(gam_cwm[[i]], parametric=FALSE)
  edf_cwm[[i]] <- data.table(edf = edf_$edf, gam_or_lm = ifelse(edf_$edf > 1.01, "gam", "lm"))
}
edf_cwm <- rbindlist(edf_cwm, idcol = "region_cwm")

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
        TPG15_fam, # Added afterwards, not consistently most important, but fifth most important in Midwest
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

# EDF for decision between LM or GAM 
edf_tpg_fam <- list()
for(i in names(gam_tpg_fam)){
  edf_ <- broom::tidy(gam_tpg_fam[[i]], parametric=FALSE)
  edf_tpg_fam[[i]] <- data.table(edf = edf_$edf, gam_or_lm = ifelse(edf_$edf > 1.01, "gam", "lm"))
}
edf_tpg_fam <- rbindlist(edf_tpg_fam, idcol = "region_tpg_fam")

# Plots
# - California all LM except for TPG10_fam 
# - Midwest All LM except for TPG5_fam & TPG12_fam
# - Northeast all GAM except TPG1_fam & TPG9_fam
# - Northwest all GAM except TPG10_fam & TPG12_fam
# - Southeast all LM except TPG1_fam
# linear model California and Midwest, maybe Southeast as well
par(mfrow = c(6,7))
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

# Added T9_genus and T15_genus, as they were the most important TPGs in some regions
tpg_comb_genus <- rbindlist(trait_groups_rel_final$genus_lvl, idcol = "region", fill = TRUE) %>%
    .[, .(region, T4_genus, T10_genus, T9_genus, T12_genus, T15_genus, max_log_tu)] %>%
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
names(lm_tpg_gen)
names(gam_tpg_gen)

# EDF
edf_tpg_genus <- list()
for(i in names(gam_tpg_gen)){
  edf_ <- broom::tidy(gam_tpg_gen[[i]], parametric=FALSE)
  edf_tpg_genus[[i]] <- data.table(edf = edf_$edf, gam_or_lm = ifelse(edf_$edf > 1.01, "gam", "lm"))
}
edf_tpg_genus <- rbindlist(edf_tpg_genus, idcol = "region_tpg_gen")

# Diagnostics
# GAM
par(mfrow = c(5,3))
for(i in names(gam_tpg_gen)){
    plot(gam_tpg_gen[[i]], main = i, residuals = TRUE, pch = 1, ylab = "max_log_TU")
}

# LM
for(i in names(sim_res_lm_tpg_gen)){
    plot(sim_res_lm_tpg_gen[[i]], main = i)
}


# Summary plots ----

# Load the data
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
# Subset to the most consistent CWM trait: size large
cwm_plot <- ggplot(cwm_combined[trait=="size_large", ], aes(x = cwm_val, y = max_log_tu)) +
    facet_wrap(.~as.factor(region)) +
    geom_point(alpha = 0.5, shape = 1, size = 2) +
    geom_smooth(
        data = ~ .x[region %in% c("California", "Midwest"), ],
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "steelblue"
    ) +
    geom_smooth(
        data = ~ .x[region %in% c("Northeast", "Northwest", "Southeast"), ],
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
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(
            family = "Roboto Mono",
            size = 16
        ),
        axis.text.y = element_text(
            family = "Roboto Mono",
            size = 16
        ),
        strip.text = element_text(
            family = "Roboto Mono",
            size = 16
        )
    )
ggsave(file.path(path_paper, "Graphs", "gam_lm_cwm_toxicity.png"),
    width = 35,
    height = 20,
    units = "cm"
)

## TPG plot ----
# Removed T9 and T15 again, as they were not among the most consistent responding TPGs
lm_tpg_fam_ls <- edf_tpg_fam[gam_or_lm=="lm",][["region_tpg_fam"]]
gam_tpg_fam_ls <- edf_tpg_fam[gam_or_lm=="gam",][["region_tpg_fam"]]
lm_tpg_gen_ls <- edf_tpg_genus[gam_or_lm=="lm",][["region_tpg_gen"]]
gam_tpg_gen_ls <- edf_tpg_genus[gam_or_lm=="gam",][["region_tpg_gen"]]


tpg_plot <- ggplot(tpg_comb[!tpg %in% c("TPG9_genus", "TPG15_genus", "TPG15_fam")], 
                   aes(x = tpg_val, y = max_log_tu)) +
    facet_grid(as.factor(tpg) ~ as.factor(region)) +
    geom_point(alpha = 0.5, shape = 1, size = 2) +
    geom_smooth(
        data = ~ .x[id %in% lm_tpg_fam_ls],
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "gold"
    ) +
    geom_smooth(
        data = ~ .x[id %in% gam_tpg_fam_ls, ],
        method = "gam",
        formula = y ~ s(x),
        se = TRUE,
        color = "gold"
    ) +
    geom_smooth(
        data = ~ .x[id %in% c(lm_tpg_gen_ls, "California_T10_genus", "Midwest_T15_genus"), ],
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "forestgreen"
    ) +
    geom_smooth(
        data = ~ .x[id %in% gam_tpg_gen_ls, ],
        method = "gam",
        formula = y ~ s(x),
        se = TRUE,
        color = "forestgreen"
    ) +
    scale_x_continuous(limits = c(0, 0.5), 
                       breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
    labs(x = "Abundance weighted-fraction TPG", 
         y = "Max logTU") +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(
            family = "Roboto Mono",
            size = 16
        ),
        axis.text.y = element_text(
            family = "Roboto Mono",
            size = 16
        ),
        strip.text = element_text(
            family = "Roboto Mono",
            size = 16
        ),
        panel.spacing.x = unit(6, "mm")
    )
ggsave(file.path(path_paper, "Graphs", "gam_lm_tpg_toxicity.png"),
    width = 38,
    height = 45,
    units = "cm"
)

# How many values are above 0.3 for TPG values?
tpg_comb[, .(.N, tpg_val, tpg), by = "region"] %>% 
.[tpg_val > 0.5, (.N/N)*100, by = "tpg"] %>% 
unique()


# Summary table ----
## CWM table ----
lm_cwm_ls <- edf_cwm[gam_or_lm=="lm",][["region_cwm"]]
gam_cwm_ls <- edf_cwm[gam_or_lm=="gam",][["region_cwm"]]

cwm_table <- rbind(
    lapply(
        lm_cwm[lm_cwm_ls],
        function(x) {
            broom::tidy(x)
        }
    ) %>%
        rbindlist(., idcol = "Region") %>% 
        setnames(., "statistic", "t"),
    rbind(
        lapply(
            gam_cwm[gam_cwm_ls],
            function(x) {
                broom::tidy(x, parametric = TRUE)
            }
        ) %>%
            rbindlist(., idcol = "Region")%>% 
            setnames(., "statistic", "t"),
        lapply(
            gam_cwm[gam_cwm_ls],
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

# Explained variance & deviance
cwm_r2 <- lapply(lm_cwm[lm_cwm_ls], function(x) {
  broom::glance(x)
}) |>
  rbindlist(, idcol = 'Region') |>
  _[, c('Region', 'r.squared')] |> 
  _[, r.squared := round(r.squared, digits=3)]

cwm_dev <- lapply(gam_cwm[gam_cwm_ls], function(x) {
  data.table('dev_explained' = summary(x)$dev.expl)
}) |>
  rbindlist(, idcol = 'Region') |> 
  _[, dev_explained := round(dev_explained, digits=3)]

cwm_table[cwm_r2, r_squared := r.squared, on = 'Region']
cwm_table[cwm_dev, dev_explained := dev_explained, on = 'Region']
cwm_table[,r_squared_dev_expl := coalesce(r_squared, dev_explained)]
cwm_table[,p_p_gam := coalesce(p.value, p.value_gam)]
cwm_table[, cwm_trait := sub("(.+)_(.+_.+)", "\\2", Region)]
cwm_table[, Region := sub("(.+)_(.+_.+)", "\\1", Region)]
# saveRDS(cwm_table, file.path(path_cache, "cwm_tabl_publ.rds"))

## TPG table ----

### Family level ----
lm_tpg_fam_ls <- edf_tpg_fam[gam_or_lm=="lm",][["region_tpg_fam"]]
gam_tpg_fam_ls <- edf_tpg_fam[gam_or_lm=="gam",][["region_tpg_fam"]]

tpg_table_fam <- rbind(
    lapply(
        lm_tpg_fam[
            names(lm_tpg_fam) %in% lm_tpg_fam_ls],
        function(x) {
            broom::tidy(x)
        }
    ) %>%
        rbindlist(., idcol = "Region_TPG") %>% 
        setnames(., "statistic", "t"),
    rbind(
        lapply(
            gam_tpg_fam[
                names(gam_tpg_fam) %in% gam_tpg_fam_ls
            ], function(x) {
                broom::tidy(x, parametric = TRUE)
            }
        ) %>%
            rbindlist(., idcol = "Region_TPG") %>% 
            setnames(., "statistic", "t"),
        lapply(
            gam_tpg_fam[
                names(gam_tpg_fam) %in% gam_tpg_fam_ls], function(x) {
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
# tpg_table_fam[(p.value > 0.05 | p.value_gam > 0.05) & term != "(Intercept)"]
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

# Explained variance/deviance
tpg_fam_r2 <- lapply(lm_tpg_fam[names(lm_tpg_fam) %in% lm_tpg_fam_ls], function(x) {
  broom::glance(x)
}) |>
  rbindlist(, idcol = 'Region_TPG') |>
  _[, c('Region_TPG', 'r.squared')] |>
  _[, r.squared := round(r.squared, digits = 3)]

tpg_fam_dev <- lapply(gam_tpg_fam[names(gam_tpg_fam) %in% gam_tpg_fam_ls], function(x) {
  data.table('dev_explained' = summary(x)$dev.expl)
}) |>
  rbindlist(, idcol = 'Region_TPG') |>
  _[, dev_explained := round(dev_explained, digits = 3)]

tpg_table_fam[tpg_fam_r2, r_squared := r.squared, on = 'Region_TPG']
tpg_table_fam[tpg_fam_dev, dev_explained := dev_explained, on = 'Region_TPG']
tpg_table_fam[,r_squared_dev_expl := coalesce(r_squared, dev_explained)]
tpg_table_fam[,p_p_gam := coalesce(p.value, p.value_gam)]
tpg_table_fam[, TPG := sub("(.+)_(.+_.+)", "\\2", Region_TPG)]
tpg_table_fam[, Region := sub("(.+)_(.+_.+)", "\\1", Region_TPG)]
# saveRDS(tpg_table_fam, file.path(path_cache, "tpg_tabl_publ.rds"))


### Genus level ----
lm_tpg_gen_ls <- edf_tpg_genus[gam_or_lm=="lm",][["region_tpg_gen"]]
gam_tpg_gen_ls <- edf_tpg_genus[gam_or_lm=="gam",][["region_tpg_gen"]]

tpg_table_genus <- rbind(
    lapply(
        lm_tpg_gen[
            names(lm_tpg_gen) %in% lm_tpg_gen_ls
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
                names(gam_tpg_gen) %in% gam_tpg_gen_ls
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
                names(gam_tpg_gen) %in% gam_tpg_gen_ls
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

# Explained Variance/Deviance
tpg_gen_r2 <- lapply(lm_tpg_gen[names(lm_tpg_gen) %in% lm_tpg_gen_ls], function(x) {
  broom::glance(x)
}) |>
  rbindlist(, idcol = 'Region_TPG') |>
  _[, c('Region_TPG', 'r.squared')] |>
  _[, r.squared := round(r.squared, digits = 3)]

tpg_gen_dev <- lapply(gam_tpg_gen[names(gam_tpg_gen) %in% gam_tpg_gen_ls], function(x) {
  data.table('dev_explained' = summary(x)$dev.expl)
}) |>
  rbindlist(, idcol = 'Region_TPG') |>
  _[, dev_explained := round(dev_explained, digits = 3)]

# Postprocessing
tpg_table_genus[tpg_gen_r2, r_squared := r.squared, on = 'Region_TPG']
tpg_table_genus[tpg_gen_dev, dev_explained := dev_explained, on = 'Region_TPG']
tpg_table_genus[,r_squared_dev_expl := coalesce(r_squared, dev_explained)]
tpg_table_genus[,p_p_gam := coalesce(p.value, p.value_gam)]
tpg_table_genus[, TPG := sub("(.+)_(.+_.+)", "\\2", Region_TPG)]
tpg_table_genus[, Region := sub("(.+)_(.+_.+)", "\\1", Region_TPG)]
# saveRDS(tpg_table_genus, file.path(path_cache, "tpg_tabl_genus_publ.rds"))


# Tables for publication ----
## CWM ----
cwm_table <- readRDS(file.path(path_cache, "cwm_tabl_publ.rds"))
cwm_table[, .(
    Region,
    cwm_trait,
    term,
    estimate = round(estimate, digits = 2),
    std.error = round(std.error, digits = 2),
    t = round(t, digits = 2),
    p.value = round(p.value, digits = 3),
    edf = round(edf, digits = 2),
    F = round(F, digits = 2),
    p.value_gam = round(p.value_gam, digits = 3),
    R2 = round(r_squared, digits = 3), 
    Deviance = round(dev_explained, digits = 3)
)] %>%
    .[, `:=`(
        p.value = fifelse(p.value == 0, "< 0.001", paste(p.value)),
        p.value_gam = fifelse(p.value_gam == 0, "< 0.001", paste(p.value_gam))
    )] %>%
    .[order(Region), ] %>%  
fwrite(., file.path(path_paper, "Tables", "CWM_log_tu_results.csv"))

## TPGs ----
tpg_table_fam <- readRDS(file.path(path_cache, "tpg_tabl_publ.rds"))
tpg_table_genus <- readRDS(file.path(path_cache, "tpg_tabl_genus_publ.rds"))

tpg_table_final <- rbind(
    tpg_table_fam[, c(
        "Region",
        "TPG",
        "term",
        "estimate",
        "std.error",
        "t",
        "p_p_gam",
        "edf",
        "F",
        "r_squared_dev_expl"
    )],
    tpg_table_genus[, c(
        "Region",
        "TPG",
        "term",
        "estimate",
        "std.error",
        "t",
        "p_p_gam",
        "edf",
        "F",
        "r_squared_dev_expl"
    )]
) 
tpg_table_final[, Taxonomic_resolution := ifelse(grepl("_fam", TPG), "family_level", "genus_level")]
tpg_table_final = tpg_table_final[order(Taxonomic_resolution, Region), ]
tpg_table_final[, edf := round(edf, digits = 2)] %>%
  .[, .SD, .SDcols = !names(tpg_table_final) %in% "Taxonomic_resolution"] %>%
    fwrite(., file.path(path_paper, "Tables", "TPG_log_tu_results.csv"))