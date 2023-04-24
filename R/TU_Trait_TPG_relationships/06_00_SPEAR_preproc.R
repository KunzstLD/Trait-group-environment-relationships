# ________________________________________________
# SPEAR preproc
# ________________________________________________
abund <- readRDS(file.path(path_cache, "total_abund_CEOPT_corrected.rds"))
abund <- abund[, .(Region, site, STAID, taxon, species, genus, family, order, abundance)]

# Preparation for Indicate tool
# Split into a file for each region and order
abund_ls <- split(abund, f = abund$Region)
names(abund_ls)[[4]] <- "Northwest" 
abund_ls <- lapply(abund_ls, function(x) x[order(Site), ])

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
# - Not all taxa can be merged with the trait data
    # TODO: Check how many could not be matched!
# - Strange transformation of abundances log(4x+1)
# - Currently obtained results are thus preliminary
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

# Simple linear regression
lm_spear <- list()
for (i in names(result_spear)) {
    lm_spear[[i]] <- lm(max_log_tu ~ SPEAR_Pestizide, data = result_spear[[i]])
}

regrres_spear <- lapply(lm_spear, function(x) lm_summary_to_dt(lm_obj = x)) %>%
    rbindlist(., id = "Region")
setnames(regrres_spear, "Pr(>|t|)", "p_value")
regrres_spear[, Region := sub("SPEAR_", "", Region)]

# Regression results
regrres_spear[p_value <= 0.05, ]

# TODO: check assumptions, same for EPT