# __________________________________________________________________________________________________
# Post analysis for Discussion 
# __________________________________________________________________________________________________


# Do TPGs that are mainly dominated by one order have a higher predictive performance than those that contain multiple orders? ----

## Family level ---- 
abund_family <- readRDS(file.path(path_cache, "tpg_taxonomic_composition_family.rds"))
abund_family[, Region := ifelse(Region=="PN", "Northwest", Region)]
abund_family[, Region_TPG := paste0(Region, "_TPG", group, "_fam")]

abund_family[Region=='California' & group==5,]

# TPGs with results for ind. LM/GAM models
# Only calculated for the most responsive TPGs
tpg_table_fam <- readRDS(file.path(path_cache, "tpg_tabl_publ.rds"))
abund_family[tpg_table_fam, r_squared_dev_expl := i.r_squared_dev_expl, on = "Region_TPG"]

# Most responsive TPGs and proportion of orders 
abund_family[!is.na(r_squared_dev_expl), ]

# California:
# Over 80 of one order (impurity score): 
# - TPG5 (0.22)
# - TPG8 (0.12)
# - TPG10 (0.18)
abund_family[Region_TPG %in% c("California_TPG5_fam", "California_TPG10_fam", "California_TPG15_fam", "California_TPG8_fam", "California_TPG12_fam")]

# Midwest
# Over 80 of one order (impurity score): 
# - TPG1 (0.1)
# - TPG8 (0.21)
# - TPG10 (0.12)
# - TPG15 (0.09)
abund_family[Region_TPG %in% c("Midwest_TPG1_fam", "Midwest_TPG2_fam", "Midwest_TPG8_fam", "Midwest_TPG10_fam", "Midwest_TPG15_fam")]

# Northeast
# Over 80 of one order (impurity score):
# - TPG5 (0.45)
# - TPG6 (0.09)
# - TPG8 (0.21)
# - TPG9 (0.05)
abund_family[Region_TPG %in% c("Northeast_TPG1_fam", "Northeast_TPG5_fam", "Northeast_TPG6_fam", "Northeast_TPG8_fam", "Northeast_TPG9_fam")]

# Northwest
# Over 80 of one order (impurity score):
# - TPG1 (0.12)
# - TPG5 (0.3)
# - TPG9 (0.04)
abund_family[Region_TPG %in% c("Northwest_TPG1_fam", "Northwest_TPG2_fam", "Northwest_TPG5_fam", "Northwest_TPG9_fam", "Northwest_TPG12_fam")]

# Southeast
# Over 80 of one order (impurity score):
# - TPG2 (0.22)
# - TPG9 (0.04)
# - TPG12 (0.09)
abund_family[Region_TPG %in% c("Southeast_TPG1_fam", "Southeast_TPG2_fam", "Southeast_TPG9_fam", "Southeast_TPG10_fam", "Southeast_TPG12_fam")]


## Genus level ----
abund_genus <- readRDS(file.path(path_cache, "tpg_taxonomic_composition_genus.rds"))
abund_genus[, Region := ifelse(Region=="PN", "Northwest", Region)]
abund_genus[, Region_TPG := paste0(Region, "_TPG", group, "_genus")]

# California:
# Over 80 of one order (impurity score): 
# - TPG2 (0.07)
# - TPG4 (0.4)
# - TPG5 (0.10)
abund_genus[Region_TPG %in% c("California_TPG1_genus", "California_TPG2_genus", "California_TPG4_genus", "California_TPG5_genus", "California_TPG6_genus")]

# Midwest
# Over 80 of one order (impurity score): 
# - TPG1 (0.08)
# - TPG3 (0.25)
# - TPG4 (0.11)
abund_genus[Region_TPG %in% c("Midwest_TPG1_genus", "Midwest_TPG3_genus", "Midwest_TPG4_genus", "Midwest_TPG9_genus", "Midwest_TPG11_genus")]

# Northeast
# Over 80 of one order (impurity score):
# - TPG13 
# - TPG4 
# - TPG10 
abund_genus[Region_TPG %in% c("Northeast_TPG4_genus", "Northeast_TPG8_genus", "Northeast_TPG10_genus", "Northeast_TPG12_genus", "Midwest_TPG13_genus")]

# Northwest
# Over 80 of one order (impurity score):
# - TPG3 
# - TPG7 
# - TPG10 
# - TPG15
abund_genus[Region_TPG %in% c("Northwest_TPG3_genus", "Northwest_TPG7_genus", "Northwest_TPG10_genus", "Northwest_TPG12_genus", "Northwest_TPG15_genus")]

# Southeast
# Over 80 of one order (impurity score):
# - TPG5
# - TPG6 
abund_genus[Region_TPG %in% c("Southeast_TPG5_genus", "Southeast_TPG6_genus", "Southeast_TPG8_genus", "Southeast_TPG10_genus", "Southeast_TPG12_genus")]

# Taxonomic variation within TPGs ----

# Mean proportion of one order in the TPGs (genus level)
abund_genus[, mean(prop_ord), by = "Region_TPG"]
abund_genus[, mean(prop_ord), by = "Region_TPG"] |> 
  _[, mean(V1)]

# Per Region
abund_genus[, mean(prop_ord), by = "Region"] 
# abund_genus[, sd(prop_ord), by = "Region_TPG"] 


# Mean proportion of one order in the TPGs (family level)
abund_family[, mean(prop_ord), by = "Region_TPG"]
abund_family[, mean(prop_ord), by = "Region_TPG"] |> 
  _[, mean(V1)]

abund_family[, mean(prop_ord), by = "Region"] 
# abund_family[, sd(prop_ord), by = "Region_TPG"] 


# T-test
dcast(abund_genus, group + order ~ Region, value.var = "prop_ord", fill=TRUE) 
dcast(abund_family, group + order ~ Region, value.var = "prop_ord", fill=TRUE) 


# Sensitivity organic per order ----
trait_matrix <- readRDS(file.path(path_cache, "trait_matrix_final.rds"))
trait_names <-
  grep("feed.*|resp.*|locom.*|size.*|volt.*|sensitivity.*",
       names(trait_matrix),
       value = TRUE)
trait_matrix <- trait_matrix[, .SD, .SDcols = c(
  "Region",
  "taxon",
  "species",
  "genus",
  "family",
  "order",
  trait_names
)]

trait_matrix_ls <- split(trait_matrix, f = trait_matrix$Region)
trait_matrix_lf <- melt(
  trait_matrix,
  id.vars = c(
    "Region",
    "taxon",
    "species",
    "genus",
    "family",
    "order"
  ),
  variable.name = "trait"
)
trait_matrix_lf[trait=="sensitivity_organic", mean(value), by = "order"]


# Best correlations with pesticide toxicity
tpg_table_fam[, TPG := sub("(.+)_(.+_.+)", "\\2", Region_TPG)]
tpg_table_fam[term != "(Intercept)", .(Region_TPG, term, estimate, r_squared, dev_explained, r_squared_dev_expl, TPG)] |>
  _[order(-r_squared_dev_expl), ] |> 
  _[r_squared_dev_expl >= 0.1, ] |> 
  _[, .N, by = "TPG"]

tpg_table_fam[term != "(Intercept)", .(Region_TPG, term, estimate, r_squared, dev_explained, r_squared_dev_expl, TPG)] |> 
  _[r_squared_dev_expl >= 0.1, ] |> 
  _[TPG == "TPG5_fam", ]


# Number of S_org values that had to be assgined at lower taxonomic resolution ----

abund <- readRDS(file.path(path_cache, "total_abund_CEOPT.rds"))

tab_spear <- readRDS(file.path(path_cache, "tab_spear.rds"))

# Add values for certain Diptera families
# Der Wert -0.35 gilt dennoch mit Ausnahme von:
# Chironomidae -0.39
# Culicidae  -0.29
# Simuliidae -0.46
tab_spear[
  Ordnung == "Diptera" & Rang == "group",
  Sensitivität := -0.35
]
tab_spear <- rbind(
  tab_spear,
  data.table(
    Ordnung = "Diptera",
    Familie = "Culicidae",
    Gattung = NA_character_,
    Taxon = "Culiciadae", Sensitivität = -0.29
  ),
  fill = TRUE
)

# Species info but only genus info for S_org
# 37 species (8 % of all Taxa)
abund[(genus %in% tab_spear[Rang=="genus", ]$Taxon) & !is.na(species), species] |> unique() |> length()

# Genus info but only family info for S_org
# 275 genera (56 %)
abund[(family %in% tab_spear[Rang=="family", ]$Taxon) & is.na(species) & !is.na(genus), genus] |> unique() |> length()

# Family info but only order for S_org
abund[(order %in% tab_spear[Rang=="group", ]$Taxon) & is.na(species) & is.na(genus) & !is.na(family), family] |> unique() |> length()




