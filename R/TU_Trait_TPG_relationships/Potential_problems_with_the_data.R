# _______________________________________
# Potential problemse with the datasets
# _______________________________________
data_cwm <- readRDS(file.path(path_cache, "data_cwm_final.rds"))

# Collinearity ----

# LF
lapply(data_cwm, function(x) {
    dcast(x, site + max_log_tu ~ trait,
        value.var = "cwm_val"
    )
}) %>% 
lapply(., function(x) {
    x[, cor(.SD), .SDcols = !c("site", "max_log_tu")] %>%
        as.data.table(., keep.rownames = TRUE) %>%
        melt(., id.vars = "rn") %>%
        .[value != 1 & value > 0.5, ]
})


cor(test_midwest[, .SD, .SDcols = !c("site", "max_log_tu")]) %>%
    as.data.table(., keep.rownames = TRUE) %>%
    melt(., id.vars = "rn") %>%
    .[value != 1 & value > 0.5, ]