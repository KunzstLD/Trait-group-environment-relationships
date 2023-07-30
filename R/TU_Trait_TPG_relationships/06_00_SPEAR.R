# ________________________________________________
# SPEAR preproc
# ________________________________________________
abund <- readRDS(file.path(path_cache, "total_abund_CEOPT_corrected.rds"))
abund <- abund[, .(Region, site, STAID, taxon, species, genus, family, order, abundance)]

# Preparation for Indicate tool
# Split into a file for each region and order
abund_ls <- split(abund, f = abund$Region)
names(abund_ls)[[4]] <- "Northwest" 
abund_ls <- lapply(abund_ls, function(x) x[order(site), ])

# SPEAR cannot match all taxa
# Therefore, use lower taxonomic resolution
# Taxa that could actually merged with the trait database 
# used for the SPEAR indicator
taxa_used_spear <- load_data(
    path = path_out,
    pattern = ".+_taxa_used_.+.csv",
    format = "csv",
    name_rm_pattern = "_taxa_used"
)
lapply(
    taxa_used_spear,
    function(x) {
        data.table(
            "taxa_overall" = x[, uniqueN(`Taxa in Monitoring-Daten`)],
            "taxa_used_SPEAR" = x[`Taxa in Merkmals-Datenbank` != "", .N]
        ) %>%
            .[, .(taxa_overall,
                taxa_used_SPEAR,
                frac_used = taxa_used_SPEAR / taxa_overall
            )]
    }
)

# Taxonomy added for taxa that are missing
lapply(taxa_used_spear, function(x) {
    x[abund,
        `:=`(
            Art = i.species,
            Genus = i.genus,
            Familie = i.family,
            Ordnung = i.order
        ),
        on = c("Taxa in Monitoring-Daten" = "taxon")
    ]
})

# Try next lower taxonomic resolution for taxa that 
# could not be matched by the indicate Tool
lapply(taxa_used_spear, function(x) {
    x[`Taxa in Merkmals-Datenbank` == "",
        Taxa_new := fcase(
            !is.na(Art), Genus,
            is.na(Art) & !is.na(Genus), Familie,
            is.na(Art) & is.na(Genus) & !is.na(Familie), Ordnung
        )
    ] %>%
        .[, Taxa_new := coalesce(Taxa_new, `Taxa in Monitoring-Daten`)]
})
Map(
    function(x, y) {
        x[y, Taxa_new := i.Taxa_new,
            on = c("taxon" = "Taxa in Monitoring-Daten")
        ]
    },
    abund_ls,
    taxa_used_spear
)

# There are a few taxa on genus or family level
# that still cannot be matched with the Indicate Trait database
# Use the next lower taxonomic resolution (i.e. family or order level)
# See script 06_01_SPEAR_taxa_join_via_order.R
source(file.path(path_scr, "TU_Trait_TPG_relationships", "06_01_SPEAR_taxa_join_via_order.R"))

# Save
lapply(
    names(abund_ls),
    function(region) {
        fwrite(
            abund_ls[[region]][, .(Region, site, Taxa_new, abundance)],
            file.path(path_out, paste0("SPEAR_abund_", region, ".csv"))
        )
    }
)

# SPEAR results ----
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

# EQ Distribution (5 classes based on WRRL)
lapply(result_spear, function(x)
    x[, EQ_Pestizide := factor(EQ_Pestizide, levels = c(
        "I: High",
        "II: Good",
        "III: Moderate",
        "IV: Poor",
        "V: Bad"
    ))])
lapply(result_spear, function(x) {
    x[, total_n := .N] %>%
        .[, .(EQ_frac = round((.N / total_n) * 100, digits = 2)), by = "EQ_Pestizide"] %>%
        .[order(EQ_Pestizide), ] %>%
        unique()
})

## Correlation of SPEAR values with maxTU ----
# Merge max tu values
max_tu <- readRDS(file.path(path_cache, "max_tu.rds"))
setnames(max_tu, "TSITE_NO_WQ", "site")

result_spear$SPEAR_Midwest[abund, STAID := i.STAID, on = "site"]
lapply(result_spear, function(x) {
    on <- if ("STAID" %in% names(x)) c("STAID" = "site") else "site"
    x[max_tu, max_log_tu := i.max_log_tu, on = on]
})
result_spear <- lapply(result_spear, function(x) x[!is.na(max_log_tu), ])
result_spear <- rbindlist(result_spear, idcol = "Region", fill = TRUE)
result_spear[, Region := sub("SPEAR_", "", Region)]
result_spear[Region == "PN", Region := "Northwest"]

# GAMs & LMs ----
lm_spear <- list()
sim_res_spear <- list()
gam_spear <- list()
gam_assump_spear <- list()
for (i in unique(result_spear$Region)) {
    lm_spear[[i]] <- lm(max_log_tu ~ SPEAR_Pestizide,
        data = result_spear[Region == i, ]
    )

    sim_res_ept[[i]] <- simulateResiduals(fittedModel = lm_spear[[i]])

    gam_spear[[i]] <- gam(
        max_log_tu ~ s(SPEAR_Pestizide),
        data = result_spear[Region == i, ],
        method = "REML"
    )

    gam_assump_spear[[i]] <- gam.check(gam_spear[[i]])
}

## Diagnostics ----
# GAM
# All linear except Northeast
par(mfrow = c(2,3))
for(i in names(gam_spear)){
    plot(gam_spear[[i]], main = i, residuals = TRUE, pch = 1, ylab = "max_log_TU")
}
lapply(gam_spear, summary)

# LM
for(i in names(sim_res_ept)){
    plot(sim_res_ept[[i]], main = i)
}


# Plotting ----
ggplot(result_spear, aes(x = SPEAR_Pestizide, y = max_log_tu)) +
    facet_wrap(as.factor(Region) ~.) +
    geom_point(alpha = 0.5, shape = 1, size = 2) +
    geom_smooth(
        data = ~ .x[Region != "Northeast", ],
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "steelblue"
    ) +
    geom_smooth(
        data = ~ .x[Region == "Northeast", ],
        method = "gam",
        formula = y ~ s(x),
        se = TRUE,
        color = "steelblue"
    ) +
    labs(
        x = "SPEAR pesticides",
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
ggsave(file.path(path_paper, "Graphs", "SPEAR_toxicity.png"),
    width = 35,
    height = 40,
    units = "cm"
)
