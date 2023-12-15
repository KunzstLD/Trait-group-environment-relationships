# Create Table for publication (and extended for SI) ----
cwm_interactions <- readRDS(file.path(path_cache, "cwm_interactions.rds"))
tpg_fam_interactions <- readRDS(file.path(path_cache, "tpg_fam_interactions.rds"))
tpg_genus_interactions <- readRDS(file.path(path_cache, "tpg_genus_interactions.rds"))
spear_interactions <- readRDS(file.path(path_cache, "spear_interactions.rds"))
ept_interactions <- readRDS(file.path(path_cache, "ept_interactions.rds"))

# Analysis for paper
rbindlist(
    list(
    "cwm" = cwm_interactions,
    "tpg_fam" = tpg_fam_interactions,
    "tpg_genus" = tpg_genus_interactions,
    "spear" = spear_interactions
    ),
    idcol = "response_cat"
) %>% 
.[p.value <= 0.05 & term != "(Intercept)", ] %>%
.[term %like% ":", ] %>% 
.[term %like% "max_log_tu", ] %>%
#.[fo %like% "size_large|TPG1_fam|SPEAR", ] %>%  
.[, response := sub("(.+)( ~)(.+)", "\\1", fo)]  %>%
.[response_cat == "tpg_genus", ] %>% 
.[order(response_cat, term, Region), .(Region, response_cat, term, estimate, adj.r.squared)] #%>% 
.[, .N, by = "term"]


.[, .N, by = c("term", "response_var")]

lookup_interactions <- data.table(
    interactions = c(
        "max_log_tu * Riffle.FRC",
        "max_log_tu * Temp.median",
        "max_log_tu * orthoP_4wk.median",
        "max_log_tu * Riffle.FRC * Temp.median",
        "max_log_tu * Riffle.FRC * orthoP_4wk.median",
        "max_log_tu * Temp.median * orthoP_4wk.median"
    ),
    id = c(1:6)
)

# Filter all for significant interactions & create table ----
## EPT ----
# No significant interactions for EPT 
ept_interactions[p.value <= 0.05 & term != "(Intercept)", ] %>%
        .[term %like% ":", ]

## SPEAR ----
# spear_table <- fread(file.path(path_paper, "Tables", "SPEAR_log_tu_results.csv"))
rbind(
    # spear_interactions[p.value <= 0.05 & term != "(Intercept)", ] %>%
    #     .[term %like% ":", ] %>%
    #     .[term %like% "max_log_tu", ] %>%
    #     .[, .(Region,
    #      `Significant term` = term,
    #         Formula = fo,
    #         Estimate = round(estimate, digits = 3),
    #         `P value` = fifelse(
    #             round(p.value, digits = 3) == 0,
    #             "<0.01",
    #             paste(round(p.value, digits = 3))
    #         ),
    #         `Adj. R squared` = round(adj.r.squared, digits = 2)
    #     )] %>%
    #     .[, Response := sub("(.+)( ~)(.+)", "\\1", Formula)] %>%
    #     .[, Formula := sub("(.+~) ", "", Formula)] %>%
    #     .[lookup_interactions, `Formula ID` := i.id, on = c(Formula = "interactions")] %>%
    #     .[, .(
    #         Region,
    #         Response,
    #         `Significant term`,
    #         Estimate,
    #         `P value`,
    #         `Adj. R squared`,
    #         `Formula ID`
    #     )],

    ## CWM ----
    # Subset for size large and sensitivity_organic
    cwm_interactions[p.value <= 0.05 & term != "(Intercept)", ] %>%
        .[term %like% ":", ] %>%
        .[term %like% "max_log_tu", ] %>%
        .[!fo %like% "size_large.+", ] %>%
        .[, .(Region,
            `Significant term` = term,
            Formula = fo,
            Estimate = round(estimate, digits = 3),
            `P value` = fifelse(
                round(p.value, digits = 3) == 0,
                "<0.01",
                paste(round(p.value, digits = 3))
            ),
            `Adj. R squared` = round(adj.r.squared, digits = 2)
        )] %>%
        .[, Response := sub("(.+)( ~)(.+)", "\\1", Formula)] %>%
        .[, Formula := sub("(.+~) ", "", Formula)] %>%
        .[lookup_interactions, `Formula ID` := i.id, on = c(Formula = "interactions")] %>%
        .[, .(
            Region,
            Response,
            `Significant term`,
            Estimate,
            `P value`,
            `Adj. R squared`,
            `Formula ID`
        )],

    ## TPGs ----
    # Subset to most important TPGs
    list(
        "tpg_fam" = tpg_fam_interactions[!fo %like% "TPG1_fam.+", ],
        "tpg_genus" = tpg_genus_interactions[!fo %like% "TPG4_genus.+", ]
    ) %>%
        lapply(., function(x) {
            x[p.value <= 0.05 & term != "(Intercept)", ] %>%
                .[term %like% ":", ] %>%
                .[term %like% "max_log_tu", ] %>%
                .[, .(Region,
                    `Significant term` = term,
                    Formula = fo,
                    Estimate = round(estimate, digits = 3),
                    `P value` = fifelse(
                        round(p.value, digits = 3) == 0,
                        "<0.01",
                        paste(round(p.value, digits = 3))
                    ),
                    `Adj. R squared` = round(adj.r.squared, digits = 2)
                )] %>%
                .[, Response := sub("(.+)( ~)(.+)", "\\1", Formula)] %>%
                .[, Formula := sub("(.+~) ", "", Formula)] %>%
                .[lookup_interactions, `Formula ID` := i.id, on = c(Formula = "interactions")] %>%
                .[, .(
                    Response,
                    `Significant term`,
                    `Formula ID`,
                    Estimate,
                    `P value`,
                    `Adj. R squared`,
                    Region
                )]
        }) %>%
        rbindlist()
) %>% 
.[Region == "PN", Region := "Northwest"] %>% 
.[, `Significant term` := sub("ortho", "", `Significant term`)] %>% 
#.[, Response := factor(Response, levels = c("size_large", "TPG1_fam", "TPG4_genus", "SPEAR_Pestizide"))] %>% 
.[order(`Significant term`, Region), .(
                    Response,
                    `Significant term`,
                    `Formula ID`,
                    Estimate,
                    `P value`,
                    `Adj. R squared`, 
                    Region)] %>% 
                    print()
fwrite(., file = file.path(path_paper, "Tables", "SI_significant_interactions.csv"))
