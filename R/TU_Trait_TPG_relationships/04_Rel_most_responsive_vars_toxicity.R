# ________________________________________________
# Post analysis
# Relationship of the most impportant traits and TPGs with max logTU
# Investigated with LM 
# ________________________________________________
data_cwm_final <- readRDS(file.path(path_cache, "data_cwm_final.rds"))
most_imp_cwm <- readRDS(file.path(path_cache, "most_important_traits_cwm.rds"))
# lapply(data_cwm_final, function(x) x[trait == "sensitivity_organic", range(cwm_val)])

# CWM
cwm_combined <- rbindlist(data_cwm_final, idcol = "region", fill = TRUE) %>%
    .[trait %in% c("feed_predator", "feed_gatherer", "sensitivity_organic"), ] %>%
    .[, region := fifelse(region == "PN", "Northwest", region)]

# Three most consistent responding traits
lm_cwm <- cwm_combined[,
    {
        model <- lm(max_log_tu ~ cwm_val)
        p_value <- summary(model)$coefficients["cwm_val", "Pr(>|t|)"]
        .(
            coef_trait = coef(model)["cwm_val"],
            p_value
        )
    },
    by = c("region", "trait")
]
lm_cwm[, `:=`(
    coord_x = rep(c(
        0.35,
        0.3,
        -0.5
    ), 5),
    coord_y = rep(-8, 15)
)]
# Plot
ggplot(cwm_combined, aes(x = cwm_val, y = max_log_tu)) +
    facet_grid(as.factor(region) ~ as.factor(trait), scales = "free_x") +
    geom_point() +
    geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "steelblue"
    ) +
    geom_text(data = lm_cwm, aes(
        x = coord_x,
        y = coord_y,
        label = paste0(
            "slope = ", round(coef_trait, digits = 2), "; ",
            "p = ", round(p_value, digits = 4)
        )
    )) +
    labs(
        x = "Community-weighted mean traits",
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
ggsave(file.path(path_paper, "Graphs", paste0("rel_traits_logTu.png")),
    width = 50,
    height = 30,
    units = "cm"
)

# TPGs
# Two most consistent TPGs
# Family level
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

# LM
lm_tpg <- tpg_comb_family[,
    {
        model <- lm(max_log_tu ~ tpg_val)
        p_value <- summary(model)$coefficients["tpg_val", "Pr(>|t|)"]
        .(
            coef_tpg = coef(model)["tpg_val"],
            p_value
        )
    },
    by = c("region", "tpg")
]
lm_tpg[, `:=`(
    coord_x = rep(0.11, 10),
    coord_y = rep(-7, 10)
)]

# Plot
ggplot(tpg_comb_family, aes(x = tpg_val, y = max_log_tu)) +
    facet_grid(as.factor(region) ~ as.factor(tpg),
        labeller = as_labeller(c(
            "California" = "California",
            "Midwest" = "Midwest",
            "Northeast" = "Northeast",
            "Northwest" = "Northwest",
            "Southeast" = "Southeast",
            "T12_fam" = "TPG12_fam",
            "T5_fam" = "TPG5_fam"
        ))
    ) +
    geom_point() +
    geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "steelblue"
    ) +
    geom_text(data = lm_tpg, aes(
        x = coord_x,
        y = coord_y,
        label = paste0(
            "slope = ", round(coef_tpg, digits = 2), "; ",
            "p = ", round(p_value, digits = 4)
        )
    )) +
    labs(x = "Abundance weighted fraction TPG", y = "Max logTU") +
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
ggsave(file.path(path_paper, "Graphs", paste0("rel_tpgs_logTu_family.png")),
        width = 50,
        height = 30,
        units = "cm"
    )

# Genus level
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
rbindlist(trait_groups_rel_final$genus_lvl, idcol = "region", fill = TRUE) %>%
    .[, .(region, T1_genus, max_log_tu)] %>%  
    melt(.,
        id.vars = c("max_log_tu", "region"),
        variable.name = "tpg",
        value.name = "tpg_val"
    ) %>%
    ggplot(., aes(x = tpg_val, y = max_log_tu)) +
    facet_grid(as.factor(region) ~ as.factor(tpg)) +
    geom_point() +
    geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "steelblue"
    ) +
    labs(x = "Abundance weighted fraction TPG", y = "max logTU") +
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
ggsave(file.path(path_paper, "Graphs", paste0("rel_tpgs_logTu_genus.png")),
        width = 50,
        height = 30,
        units = "cm"
    )
