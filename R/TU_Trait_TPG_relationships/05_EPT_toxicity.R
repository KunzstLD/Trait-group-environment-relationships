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
# saveRDS(abund, file.path(path_cache, "total_abund_CEOPT_corrected.rds"))

# Toxicity
max_tu <- readRDS(file.path(path_cache, "max_tu.rds"))
setnames(max_tu, "TSITE_NO_WQ", "site")

# Calculate Fraction EPT of the whole assemblage per site
# apply sqrt transformation
abund[, abundance := sqrt(abundance)]
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
abund_subs[Region == "PN", Region := "Northwest"]

# Few sites don't have chemical information
abund_subs <- abund_subs[!is.na(max_log_tu), ]
# Few STAIDs are duplicates (same chemical info, but different ecology)
# abund_subs[!is.na(STAID) & duplicated(STAID),]

# Final EPT data
# saveRDS(abund_subs, file.path(path_cache, "ept.rds"))

# GAMs & LMs ----
# - Significant coeff. in California and Northwest
# - Negative relationship in California, Midwest, 
# Northwest, Southeast
# - Positive relationship in Northeast (tough not significant) -> why?
lm_ept <- list()
sim_res_ept <- list()
gam_ept <- list()
gam_assump_ept <- list()
for (i in unique(abund_subs$Region)) {
    lm_ept[[i]] <- lm(max_log_tu ~ frac_EPT,
        data = abund_subs[Region == i, ]
    )

    sim_res_ept[[i]] <- simulateResiduals(fittedModel = lm_ept[[i]])

    gam_ept[[i]] <- gam(
        max_log_tu ~ s(frac_EPT),
        data = abund_subs[Region == i, ],
        method = "REML"
    )

    gam_assump_ept[[i]] <- gam.check(gam_ept[[i]])
}

## Diagnostics ----

# EDF
edf_ept <- list()
for(i in names(gam_ept)){
  edf_ <- broom::tidy(gam_ept[[i]], parametric=FALSE)
  edf_ept[[i]] <- data.table(edf = edf_$edf, gam_or_lm = ifelse(edf_$edf > 1.01, "gam", "lm"))
}
edf_ept <- rbindlist(edf_ept, idcol = "Region")

# GAM
# All linear except Midwest & California
par(mfrow = c(2,3))
for(i in names(gam_ept)){
    plot(gam_ept[[i]], main = i, residuals = TRUE, pch = 1, ylab = "max_log_TU")
}
lapply(gam_ept, summary)

# LM
for(i in names(sim_res_ept)){
    plot(sim_res_ept[[i]], main = i)
}


# Plotting ----
ggplot(abund_subs, aes(x = frac_EPT, y = max_log_tu)) +
    facet_wrap(as.factor(Region) ~.) +
    geom_point(alpha = 0.5, shape = 1, size = 2) +
    geom_smooth(
        data = ~ .x[!Region %in% c("California", "Midwest"), ],
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "steelblue"
    ) +
    labs(
        x = "Fraction EPT taxa per site",
        y = "Max logTU"
    ) + 
    geom_smooth(
        data = ~ .x[Region %in% c("California", "Midwest"), ],
        method = "gam",
        formula = y ~ s(x),
        se = TRUE,
        color = "steelblue"
    ) +
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
ggsave(file.path(path_paper, "Graphs", "EPT_toxicity.png"),
    width = 35,
    height = 20,
    units = "cm"
)

# Summary Table EPT ----
ept_table <- rbind(
    lapply(
        lm_ept[c("Northeast", "Northwest", "Southeast")],
        function(x) {
            broom::tidy(x)
        }
    ) %>%
        rbindlist(., idcol = "Region") %>%
        setnames(., "statistic", "t"),
    rbind(
        lapply(
            gam_ept[c("California","Midwest")],
            function(x) {
                broom::tidy(x, parametric = TRUE)
            }
        ) %>%
            rbindlist(., idcol = "Region") %>%
            setnames(., "statistic", "t"),
        lapply(
            gam_ept[c("California","Midwest")],
            function(x) {
                broom::tidy(x, parametric = FALSE)
            }
        ) %>%
            rbindlist(., idcol = "Region") %>%
            setnames(
                .,
                c("statistic", "p.value"),
                c("F", "p.value_gam")
            ),
        fill = TRUE
    ),
    fill = TRUE
)
ept_gof <- lapply(lm_ept, broom::glance) %>%
    rbindlist(., idcol = "Region")
ept_table[, `:=`(
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
    ), 
    edf = round(edf, digits = 2),
    ref.df = round(ref.df, digits = 2)
)]

# Add R2 or deviance
ept_table[ept_gof[!Region %in% c("California", "Midwest"), ],
    r.squared := i.r.squared,
    on = "Region"
]
ept_dev <- lapply(gam_ept[c("California", "Midwest")], function(x) {
  data.table("dev_explained" = summary(x)$dev.expl)
}) |>
  rbindlist(, idcol = "Region") |> 
  _[, dev_explained := round(dev_explained, digits=3)]

ept_table[ept_dev, dev_explained := i.dev_explained, on="Region"]
ept_table[,r_squared_dev_expl := coalesce(r.squared, dev_explained)]
ept_table[, r_squared_dev_expl := round(r_squared_dev_expl, digits = 3)]
ept_table[,p_p_gam := coalesce(p.value, p.value_gam)]
setcolorder(
    ept_table,
    c("Region", "term", "estimate", "std.error", "t", "p_p_gam", "r_squared_dev_expl")
)
fwrite(ept_table, file.path(path_paper, "Tables", "EPT_log_tu_results.csv"))
