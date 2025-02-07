
# Create Table for publication (and extended for SI) ----
cwm_interactions <- readRDS(file.path(path_cache, "cwm_interactions.rds"))
tpg_fam_interactions <- readRDS(file.path(path_cache, "tpg_fam_interactions.rds"))
tpg_genus_interactions <- readRDS(file.path(path_cache, "tpg_genus_interactions.rds"))
spear_interactions <- readRDS(file.path(path_cache, "spear_interactions.rds"))
ept_interactions <- readRDS(file.path(path_cache, "ept_interactions.rds"))

# Analysis for paper
# rbindlist(
#     list(
#     "cwm" = cwm_interactions,
#     "tpg_fam" = tpg_fam_interactions,
#     "tpg_genus" = tpg_genus_interactions,
#     "spear" = spear_interactions
#     ),
#     idcol = "response_cat"
# ) %>% 
# .[p.value <= 0.05 & term != "(Intercept)", ] %>%
# .[term %like% ":", ] %>% 
# .[term %like% "max_log_tu", ] %>%
# #.[fo %like% "size_large|TPG1_fam|SPEAR", ] %>%  
# .[, response := sub("(.+)( ~)(.+)", "\\1", fo)]  %>%
# .[response_cat == "tpg_genus", ] %>% 
# .[order(response_cat, term, Region), .(Region, response_cat, term, estimate, adj.r.squared)] #%>% 
# .[, .N, by = "term"]

#_______________________________________________________________________________
# Main effects ----
lookup_formula <-
  data.table(
    Formula_id = c(1, 2, 3),
    Formula = c(
      "max_log_tu * Riffles * Temperature",
      "max_log_tu * Riffles * Phosphate",
      "max_log_tu * Temperature * Phosphate"
    )
  )

main_effects_tbl <- rbind(
  # EPT
  ept_interactions[(p.value <= 0.05) & (!term %like% "Intercept|.*:.*|max_log_tu"), ] %>%
  .[, Response := sub("(.+)( ~)(.+)", "\\1", fo)] %>%
  .[, Type := "EPT"] %>%
  .[, .(Region,
        `Significant term` = term,
        Response,
        Estimate = round(estimate, digits = 3),
        `P value` = fifelse(
          round(p.value, digits = 3) == 0,
          "<0.01",
          paste(round(p.value, digits = 3))
        ),
        `Adj. R squared` = round(adj.r.squared, digits = 2), 
        Formula = fo,
        Type)],
  # SPEAR
  spear_interactions[(p.value <= 0.05) & (!term %like% "Intercept|.*:.*|max_log_tu"), ]%>%
    .[, Response := sub("(.+)( ~)(.+)", "\\1", fo)] %>%
    .[, Type := "SPEAR"] %>%
    .[, .(Region,
          `Significant term` = term,
          Response,
          Estimate = round(estimate, digits = 3),
          `P value` = fifelse(
            round(p.value, digits = 3) == 0,
            "<0.01",
            paste(round(p.value, digits = 3))
          ),
          `Adj. R squared` = round(adj.r.squared, digits = 2), 
          Formula = fo,
          Type)],
  
  # CWM
  cwm_interactions[(p.value <= 0.05) & (!term %like% "Intercept|.*:.*|max_log_tu"), ] %>%
    .[, Response := sub("(.+)( ~)(.+)", "\\1", fo)] %>%
    .[, Type := "Traits"] %>% 
    .[, .(Region,
          `Significant term` = term,
          Response,
          Estimate = round(estimate, digits = 3),
          `P value` = fifelse(
            round(p.value, digits = 3) == 0,
            "<0.01",
            paste(round(p.value, digits = 3))
          ),
          `Adj. R squared` = round(adj.r.squared, digits = 2), 
          Formula = fo,
          Type)],
  
  # TPG Fam
  tpg_fam_interactions[(p.value <= 0.05) & (!term %like% "Intercept|.*:.*|max_log_tu"), ]%>%
    .[, Response := sub("(.+)( ~)(.+)", "\\1", fo)] %>%
    .[, Type := "TPG"] %>%
    .[, .(Region,
          `Significant term` = term,
          Response,
          Estimate = round(estimate, digits = 3),
          `P value` = fifelse(
            round(p.value, digits = 3) == 0,
            "<0.01",
            paste(round(p.value, digits = 3))
          ),
          `Adj. R squared` = round(adj.r.squared, digits = 2), 
          Formula = fo,
          Type)],

  # TPG genus
  tpg_genus_interactions[(p.value <= 0.05) & (!term %like% "Intercept|.*:.*|max_log_tu"), ]%>%
    .[, Response := sub("(.+)( ~)(.+)", "\\1", fo)] %>%
    .[, Type := "TPG"] %>%
    .[, .(Region,
          `Significant term` = term,
          Response,
          Estimate = round(estimate, digits = 3),
          `P value` = fifelse(
            round(p.value, digits = 3) == 0,
            "<0.01",
            paste(round(p.value, digits = 3))
          ),
          `Adj. R squared` = round(adj.r.squared, digits = 2), 
          Formula = fo,
          Type)]
)

# CWM Traits
main_effects_tbl[Type == "Traits", .(`Significant term`, Response, Region)] |> 
  _[order(Response), ] |> 
  unique()

# All 10 traits were associated at least with one other variable
main_effects_tbl[Type == "Traits", Response] |> unique()

# How many have associations across multiple regions? 
main_effects_tbl[Type == "Traits", .(`Significant term`, Response, Region)] |> 
  _[order(Response), ] |> 
  unique() |> 
  _[, uniqueN(Region), by = "Response"] |> 
  _[order(V1), ]


# TPG fam
main_effects_tbl[Response %like% ".*_fam", .(`Significant term`, Response, Region)] |> 
  _[order(Response), ] |> 
  unique()

# 6 out of 7 TPGs were associated at least with one other variable
main_effects_tbl[Response %like% ".*_fam", Response] |> unique()

# How many have associations across multiple regions? 
main_effects_tbl[Response %like% ".*_fam", .(`Significant term`, Response, Region)] |> 
  unique() |> 
  _[, uniqueN(Region), by = "Response"] |> 
  _[order(V1), ]
main_effects_tbl[Response == "TPG8_fam",]

# TPG genus
main_effects_tbl[Response %like% ".*_genus", .(`Significant term`, Response, Region, Estimate)] |> 
  _[order(Response), ] |> 
  unique()

# SPEAR & EPT 
main_effects_tbl[Type == "EPT", .(`Significant term`, Response, Region)] 

main_effects_tbl[Type == "SPEAR", .(`Significant term`, Response, Region)]

# Postprocessing 
main_effects_tbl[Region == "PN", Region := "Northwest"] 
main_effects_tbl[, `Significant term` := sub("Temp.median", "Temperature",`Significant term`)]
main_effects_tbl[, `Significant term` := sub("orthoP_4wk.median", "Phosphate",`Significant term`)]
main_effects_tbl[, `Significant term` := sub("Riffle.FRC", "Riffles",`Significant term`)]
main_effects_tbl[, Formula := sub("Temp.median", "Temperature", Formula)]
main_effects_tbl[, Formula := sub("orthoP_4wk.median", "Phosphate",Formula)]
main_effects_tbl[, Formula := sub("Riffle.FRC", "Riffles",Formula)]

main_effects_tbl$Response <- factor(
  main_effects_tbl$Response,
  levels = c(
    "size_large",
    "size_small",
    "feed_gatherer",
    "feed_filter",
    "feed_predator",
    "volt_semi",
    "volt_bi_multi",
    "locom_swim",
    "sensitivity_organic",
    "resp_gil",
    "TPG1_fam",
    "TPG2_fam",
    "TPG8_fam",
    "TPG9_fam",
    "TPG10_fam",
    "TPG12_fam",
    "TPG4_genus",
    "TPG10_genus",
    "frac_EPT",
    "SPEAR_Pestizide"
  )
)
main_effects_tbl[, Formula := sub("(.+)(~ )(.+)", "\\3", Formula)]
main_effects_tbl[lookup_formula, Formula_id := i.Formula_id, on = "Formula"]

# Table for SI 
main_effects_tbl[, .(`Significant term`,
                     Response,
                     Region,
                     Formula_id,
                     Estimate,
                     `P value`,
                     `Adj. R squared`)] |>
  _[order(Response, `Significant term`, Region), ] |> 
  fwrite(file.path(path_paper, "Tables", "SI_significant_main_effects.csv"))

#_______________________________________________________________________________
# Filter all for significant interactions & create table ----

## EPT ----
# No significant interactions for EPT 
ept_interactions[p.value <= 0.05 & term != "(Intercept)", ] %>%
  .[term %like% ":", ]

interaction_table <- rbind(
  
    ## Spear ----
    spear_interactions[p.value <= 0.05 & term != "(Intercept)", ] %>%
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
        .[, Type := "SPEAR"] %>%
        .[, .(
            Response,
            `Significant term`,
            Estimate,
            `P value`,
            `Adj. R squared`,
             Formula,
             Region,
             Type
        )],

    ## CWM ----
    # Subset for size large and sensitivity_organic
    cwm_interactions[p.value <= 0.05 & term != "(Intercept)", ] %>%
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
        .[, Type := "Traits"] %>%
        .[, .(
            Response,
            `Significant term`,
            Estimate,
            `P value`,
            `Adj. R squared`,
            Formula,
            Region,
            Type
        )],

    ## TPGs ----
    # Subset to most important TPGs
    list(
        "tpg_fam" = tpg_fam_interactions, # [fo %like% "TPG1_fam.+", ]
        "tpg_genus" = tpg_genus_interactions # [fo %like% "TPG4_genus.+", ]
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
                .[, Type := "TPGs"] %>%
                .[, .(
                    Response,
                    `Significant term`,
                    Formula,
                    Estimate,
                    `P value`,
                    `Adj. R squared`,
                    Region,
                    Type
                )]
        }) %>%
        rbindlist()
) 

# Postprocessing
interaction_table[Region == "PN", Region := "Northwest"] 
interaction_table[, `Significant term` := sub("Temp.median", "Temperature",`Significant term`)]
interaction_table[, `Significant term` := sub("orthoP_4wk.median", "Phosphate",`Significant term`)]
interaction_table[, `Significant term` := sub("Riffle.FRC", "Riffles",`Significant term`)]
interaction_table[, Formula := sub("Temp.median", "Temperature", Formula)]
interaction_table[, Formula := sub("orthoP_4wk.median", "Phosphate",Formula)]
interaction_table[, Formula := sub("Riffle.FRC", "Riffles",Formula)]
interaction_table[lookup_formula, Formula_id := i.Formula_id, on = "Formula"]

# Table for SI
interaction_table$Response <- factor(
  interaction_table$Response,
  levels = c(
    "size_large",
    "size_small",
    "feed_gatherer",
    "feed_filter",
    "feed_predator",
    "volt_semi",
    "volt_bi_multi",
    "sensitivity_organic",
    "resp_gil",
    "TPG1_fam",
    "TPG2_fam",
    "TPG8_fam",
    "TPG9_fam",
    "TPG10_fam",
    "TPG12_fam",
    "TPG4_genus",
    "TPG10_genus",
    "TPG12_genus",
    "SPEAR_Pestizide"
  )
)
interaction_table[, .(`Significant term`,
                     Response,
                     Region,
                     Formula_id,
                     Estimate,
                     `P value`,
                     `Adj. R squared`)] |>
  _[order(Response, `Significant term`, Region), ] |> 
  fwrite(file = file.path(path_paper, "Tables", "SI_significant_interactions.csv"))

# ______________________________________________________________________________

# Interpretation ----
# How many interactions occur per trait, tpg, spear, ept? 6
# max_log_tu:Riffle.FRC:Temp.median
# max_log_tu:Riffle.FRC:orthoP_4wk.median
# max_log_tu:Temp.median:orthoP_4wk.median
# max_log_tu:Riffle.FRC
# max_log_tu:Temp.median
# max_log_tu:orthoP_4wk.median
# (excluding similar lower order interactions)
cwm_interactions[term %like% "max_log_tu:", ] |> 
  _[, .(Region, term, fo)] |> 
  _[fo %like% "size_large", ] |> 
  _[, uniqueN(term), by = "Region"]

## CWM traits ----
# 9 out of 10 have interactions (not locom swimming)
# Most interactions in California (9) and Northeast (8), overall 23 (possible? 6*10*5)
interaction_table[Type == "Traits", .N, by = "Region"] |> 
  _[, sum(N)]

# Traits that have multiple interactions?
interaction_table[Type == "Traits",]
interaction_table[Type == "Traits", .N, by = "Response"] |> 
  _[order(N), ]

interaction_table[Type == "Traits" & Response == "feed_filter",]
interaction_table[Type == "Traits" & Response == "size_small",]
interaction_table[Type == "Traits" & Response == "volt_bi_multi",]

# Interactions that occurred often
interaction_table[Type == "Traits", .N, by = "Significant term"]
interaction_table[Type == "Traits" & `Significant term` == "max_log_tu:Temperature",]

# Interactions that occurred across regions
interaction_table[Type == "Traits", .N, by = c("Response", "Significant term")] |> 
  _[order(N), ]

interaction_table[Type == "Traits" & Response == "feed_filter", ]
interaction_table[Type == "Traits" & Response == "volt_bi_multi", ]

## TPG groups ----
# 6 out of 7
# Most interactions in Northeast, than California and Midwest
interaction_table[Response %like% "_fam", Response] |> unique()
interaction_table[Response %like% "_fam", .N, by = "Region"] |> 
  _[, sum(N)]

# TPGs fam that have multiple interactions?
interaction_table[Response %like% "_fam", .N, by = "Response"]
interaction_table[Response %like% "_fam" & Response == "TPG9_fam",]
interaction_table[Response %like% "_fam" & Response == "TPG8_fam",]

# Interactions that occurred often
interaction_table[Response %like% "_fam", .N, by = "Significant term"]
interaction_table[Response %like% "_fam" & `Significant term` == "max_log_tu:Temperature",]
interaction_table[Response %like% "_fam" & `Significant term` == "max_log_tu:Riffles",]

# Interactions that occurred across regions
interaction_table[Response %like% "_fam", .N, by = c("Response", "Significant term")] |> 
  _[order(N), ]

interaction_table[Response %like% "TPG9_fam",]
interaction_table[Response %like% "TPG2_fam",]

# Genus level
interaction_table[Response %like% "_genus", Response] |> unique()
interaction_table[Response %like% "_genus", ]
interaction_table[Response %like% "_genus", .N, by = "Response"]

# Interactions that occurred across regions -> none
interaction_table[Response %like% "_genus", .N, by = c("Response", "Significant term")] |> 
  _[order(N), ]

## SPEAR ----
interaction_table[Type == "SPEAR",]
interaction_table[Type == "SPEAR", .N, by = "Region"] |> 
  _[, sum(N)]

