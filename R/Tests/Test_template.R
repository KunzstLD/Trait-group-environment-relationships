northeast_taxonomy <- readRDS("/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Cache/northeast_taxonomy.rds")
taxonomy <- copy(northeast_taxonomy)

# AQ insects
aq_insects <- c(
  "Coleoptera",
  "Ephemeroptera",
  "Diptera",
  "Hemiptera",
  "Heteroptera",
  "Lepidoptera",
  "Megaloptera",
  "Neuroptera",
  "Odonata",
  "Plecoptera",
  "Trichoptera"
)

# Overview taxon richness
taxonomy[, total_N := .N]
taxonomy[order %in% aq_insects, aq_ins_total_N := .N]
taxonomy[, taxon_richn := .N, by = order]

# Add for how many trait information is available
# Which are the most abundant?
noa_traits <-
  readRDS(
    "/home/kunzst/Dokumente/Projects/Trait_DB/Convergence-trait-profiles/Data/Traits_US_LauraT_pp_harmonized.rds"
  )
noa_traits[, taxon := coalesce(species, genus, family, order)]
noa_traits <- noa_traits[, .SD,
                         .SDcols = patterns("species|genus|family|order|taxon|feed|resp|size|locom|volt")]
setcolorder(
  x = noa_traits,
  neworder = c(
    grep("feed.*", names(noa_traits), value = TRUE),
    grep("resp.*", names(noa_traits), value = TRUE),
    grep("locom.*", names(noa_traits), value = TRUE),
    grep("size.*", names(noa_traits), value = TRUE),
    grep("volt.*", names(noa_traits), value = TRUE)
  )
)
# Add taxonomic level in trait DB
noa_traits[, taxonomic_level := fcase(
  !is.na(species),
  "species",
  is.na(species) & !is.na(genus),
  "genus",
  is.na(species) & is.na(genus) & !is.na(family),
  "family",
  is.na(species) &
    is.na(genus) & is.na(family) & !is.na(order),
  "order"
)]

# Trait matrix ----
trait_names <-
  grep("feed|resp|size|locom|volt", names(noa_traits), value = TRUE)

# Establish complete trait matrix
trait_matrix <- noa_traits[taxon %in% taxonomy$taxon, ]
trait_matrix <- rbind(
  trait_matrix,
  taxonomy[
    !taxon %in% noa_traits$taxon,
    .(taxon, species, genus, family, order, taxonomic_level)
  ],
  fill = TRUE
)

# Add identifier for aq. insects and non-insects
trait_matrix[, type := fcase(
  order %in% aq_insects, "insect",
  !order %in% aq_insects, "non_insect"
)]

# Go until family level, not below
trait_matrix <-
  trait_matrix[taxonomic_level %in% c("species", "genus", "tribe", "subfamily", "family"),]

# Load trait data from Bob Zuelig 
# locom sessil and feed parasite not defined
# Some weird names in OTU which are not used now but will be manually checked later 
# (e.g., taxa with "/" and "genus")
noa_trait_enhanced <-
  readRDS(
    "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Data/NoA/trait_db_enhanced_preproc.rds"
  )
noa_trait_enhanced_lf <- melt(
  noa_trait_enhanced,
  id.vars = c("OTU_new",
              "Phylum",
              "Class",
              "Order",
              "Family",
              "Genus"),
  variable.name = "trait"
)

# Trait matrix into long format for merge
trait_matrix_lf <-
  melt(
    trait_matrix,
    id.vars = c(
      "species",
      "genus",
      "family",
      "order",
      "taxon",
      "taxonomic_level",
      "type"
    ),
    variable.name = "trait"
  )

# Merge Bob Zuelig traits when no information is there
trait_matrix_lf[noa_trait_enhanced_lf,
                value := fifelse(is.na(value), i.value, value),
                on = c(taxon = "OTU_new",
                       trait = "trait")]
trait_matrix <-
  dcast(trait_matrix_lf, ... ~ trait, value.var = "value")

# Update column for complete or incomplete trait info
trait_matrix[, ind_trait_info := apply(.SD, 1, function(y) {
  sum(!is.na(y))
}),
.SDcols = trait_names
]
trait_matrix[, trait_info := fcase(
  ind_trait_info == length(trait_names),
  "complete",
  ind_trait_info != length(trait_names) & ind_trait_info > 0,
  "incomplete",
  ind_trait_info == 0,
  "no_information"
)]
# trait_matrix[trait_info == "complete", .N]

## Abundance and trait data `r params$region`
# Load abundance data
ca_abund <- fread(
  "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Data/NoA/CA_NoAmbig_Species Matrix n82.csv"
)

# Conversion to long format
# Some integer columns which will be coerced to double
ca_abund_lf <- melt(
  ca_abund,
  id.vars = c(
    "TSTAID",
    "CollectionDate",
    "SMCOD",
    "BU_ID"
  ),
  variable.name = "taxon",
  value.name = "abundance"
)
ca_abund_lf[taxon %like% "Concha", taxon := "Conchapelopia"]
abund <- copy(ca_abund_lf)

# Give the same name to all site columns 
setnames(
  abund,
  old = c("TSTAID", "SCODE"),
  new = c("site", "site"),
  skip_absent = TRUE
)

# Merge taxonomy
abund[taxonomy,
      `:=`(
        species = i.species,
        genus = i.genus,
        family = i.family,
        order = i.order,
        taxonomic_level = i.taxonomic_level
      ),
      on = "taxon"
]

# Calculate at how many sites a given taxon occurs
# Create ind. column
abund[, nr_sites := uniqueN(site)]
abund[, occr_at_sites := sum(ifelse(abundance > 0, 1, 0)), by = taxon]
abund[, perc_occr_at_sites := round((occr_at_sites/nr_sites) * 100, digits = 2)]

# Calculate total abundance
abund[, total_abund := round(sum(abundance)), by = taxon]
abund_unq <- unique(abund[, .(taxon, order, total_abund)])
abund_unq[, percentile := round(percent_rank(total_abund), digits = 2), by = order]
abund[abund_unq, percentile := i.percentile, on = "taxon"]
abund[, percentile_category_order := fcase(
  percentile >= 0.75, "75th_and_above",
  percentile < 0.75 & percentile >= 0.5, "between_75th_50th",
  percentile < 0.5 & percentile >= 0.25, "between_50th_25th",
  percentile < 0.25 & percentile >= 0.05, "between_25th_05th",
  percentile < 0.05, "below_5th"
)]
abund <- abund[taxonomic_level %in% c("species", "genus", "tribe", "subfamily", "family"), ]
abund[, log_total_abund := round(log(total_abund + 1), digits = 2)]

# Add identifier of complete or incomplete trait information
abund[trait_matrix,
      `:=`(
        trait_info = i.trait_info,
        type = i.type
      ),
      on = "taxon"
]

# Add to trait matrix percentile rank
trait_matrix[abund,
             percentile_category_order := i.percentile_category_order,
             on = "taxon"
]

## Data gaps `r params$region`
# For which traits per order is information missing?
compl_overview_insect <- completeness_trait_data(
  x = trait_matrix[type == "insect", ],
  non_trait_cols = c(
    "species",
    "genus",
    "family",
    "order",
    "taxon",
    "taxonomic_level",
    "ind_trait_info",
    "trait_info",
    "type"
  )
)
compl_overview_insect <- as.data.table(compl_overview_insect)
compl_overview_insect[, V2 := sub("\\^", "", V2)]
compl_overview_insect[, Type := "insect"]
setnames(compl_overview_insect,
         old = c("V1", "V2"),
         new = c("Available information [%]", "Grouping feature")
)

compl_overview_non_insect <- completeness_trait_data(
  x = trait_matrix[type == "non_insect", ],
  non_trait_cols = c(
    "species",
    "genus",
    "family",
    "order",
    "taxon",
    "taxonomic_level",
    "ind_trait_info",
    "trait_info",
    "type"
  )
)
compl_overview_non_insect <- as.data.table(compl_overview_non_insect)
compl_overview_non_insect[, V2 := sub("\\^", "", V2)]
compl_overview_non_insect[, Type := "non_insect"]
setnames(compl_overview_non_insect,
         old = c("V1", "V2"),
         new = c("Available information [%]", "Grouping feature")
)

### Trait aggregation & trait assignment `r params$region`

# Subfamilies, tribes, and genera for which traits on family level can be assigned
# missing_tribe_subfam_genus <- unique(trait_matrix[taxonomic_level %in% c("tribe", "subfamily", "genus") &
#   trait_info == "no_information", .(taxon, family)])
missing_tribe_subfam_genus <-
  unique(trait_matrix[taxonomic_level %in% c("tribe", "subfamily", "genus") &
                        trait_info %in% c("no_information", "incomplete"), .(taxon, family)])

# Get available trait information on family level for subfamilies, tribes and genera
trait_tribe_subfam_genus <- trait_matrix[family %in% unique(missing_tribe_subfam_genus$family) &
                                           trait_info == "complete" & taxonomic_level == "family", ] %>%
  melt(., id.vars = c(
    "species",
    "genus",
    "family",
    "order",
    "taxon",
    "taxonomic_level",
    "ind_trait_info",
    "trait_info",
    "type",
    "percentile_category_order"
  ))

# Merge back tribe, subfamily, or genera to family column
# Intermediate dataset that will be merged back to the complete trait matrix
# During merging the merge() function is not to ensure that all taxa within
# a family are merged back, e.g. multiple genera of Baetidae ("Acentrella", "Procloeon", ...)
trait_tribe_subfam_genus <- merge(
  x = trait_tribe_subfam_genus[, .SD, .SDcols = !"taxon"],
  y = missing_tribe_subfam_genus,
  by = "family",
  allow.cartesian = TRUE
)

# Taxa on tribe, subfamily or genus level for which no family level trait information could be assigned
families_remain_tribe_subfam <-
  unique(missing_tribe_subfam_genus[!family %in% unique(trait_tribe_subfam_genus$family), family])

# Families for which trait information could be obtained via trait aggregation
family_incomp_no_info <- unique(trait_matrix[taxonomic_level == "family" &
                                               trait_info %in% c("incomplete", "no_information"), family])

# Aggregate NOA trait data to family level
# and use these for genera, subfamilies, tribes, or families with missing trait information
# (e.g., vector families_remain_tribe_subfam)
# Add an ID to trait matrix for as unique identifier
trait_matrix[, unique_id := paste0("id_", 1:nrow(trait_matrix))]
noa_traits[, unique_id := paste0("id_", 1:nrow(noa_traits))]
family_lvl_aggr <- direct_agg(
  trait_data = noa_traits[, .SD,
                          .SDcols = c(trait_names,
                                      "species",
                                      "genus",
                                      "family",
                                      "order",
                                      "unique_id")],
  non_trait_cols = c("species",
                     "genus",
                     "family",
                     "order",
                     "unique_id"),
  method = median
)
family_lvl_aggr_final <- na.omit(family_lvl_aggr[, .SD,
                                                 .SDcols = c("family", trait_names, "unique_id")])
family_lvl_aggr_final[family_lvl_aggr,
                      order := i.order,
                      on = "family"
]

# Add type and trait_info
family_lvl_aggr_final[trait_matrix,
                      type := i.type,
                      on = "family"]
family_lvl_aggr_final[, trait_info := "complete"]

# Long format
family_lvl_aggr_final <- melt(
  family_lvl_aggr_final,
  id.vars = c("family",
              "order",
              "type",
              "trait_info",
              "unique_id")
)

# First, get subfamilies, tribes, and genera without or incomplete trait information
# These are genera/tribes/subfamilies to which trait information on family level could not be assigned 
# -> instead, aggregated trait information on family level is assigned (hence, in the conditional statement
# a few families are not considered)
trait_tribe_subfam_genus_aggr <-
  trait_matrix[taxon %in% missing_tribe_subfam_genus[!family %in% unique(trait_tribe_subfam_genus$family), taxon] &
                 family %in% unique(family_lvl_aggr_final$family), ] %>%
  melt(
    .,
    id.vars = c(
      "species",
      "genus",
      "family",
      "order",
      "taxon",
      "taxonomic_level",
      "ind_trait_info",
      "trait_info",
      "type",
      "unique_id",
      "percentile_category_order"
    )
  )
trait_tribe_subfam_genus_aggr[family_lvl_aggr_final,
                              value := i.value,
                              on = c("family", "variable")]

# Bring together tribes, subfamilies, and genera with traits assigned from family lvl
# + tribes, subfamilies, and genera with aggregated traits 
# trait_tribe_subfam_genus
trait_tribe_subfam_genus <- rbind(trait_tribe_subfam_genus,
                                  trait_tribe_subfam_genus_aggr, 
                                  fill = TRUE)

# Not all missing trait information can be covered
trait_tribe_subfam_genus <- trait_tribe_subfam_genus[!is.na(value), ]

# With the assigned trait information from family level, change trait_info column to complete
# trait_tribe_subfam_genus[is.na(value),]
trait_tribe_subfam_genus[!trait_info %in% "complete", trait_info := "complete"]

# Now convert the trait matrix to lf for effective merge
trait_matrix_cp <- copy(trait_matrix)
trait_matrix_cp <- melt(trait_matrix_cp, id.vars = c(
  "species",
  "genus", 
  "family",
  "order",
  "taxon",
  "taxonomic_level",
  "trait_info",
  "type",
  "unique_id",
  "ind_trait_info",
  "percentile_category_order"
))

# Merge back for subfamilies, tribes, and genera
# Merge conditional on NA values, don't want to "overwrite" existing trait information
# with aggregated trait information
# Merge on taxon column!
trait_matrix_cp[trait_tribe_subfam_genus,
                value := fifelse(is.na(value), i.value, value),
                on = c("taxon", "variable")]
# trait_matrix_cp[taxon == "Microtendipes", ]
# Check: Microtendipes, resp is NA, feeding is filter 1, locom is crawler 1, size NA, volt is 
# bi_multi 1
# trait_matrix_cp[taxon == "Microtendipes", ]
# tribe_subfam_genus_assign[taxon == "Microtendipes", ]
# 49 still missing
# trait_matrix_cp[, uniqueN(taxon)]
# trait_matrix_cp[is.na(value), uniqueN(taxon)]

# N Matches expectations
# trait_matrix_cp[, uniqueN(taxon)]-trait_matrix_cp[is.na(value), uniqueN(taxon)]

# Merge aggregated family level traits back
# First, create intermediate dataset
trait_matrix_cp_fam <-
  trait_matrix_cp[taxonomic_level == "family" &
                    family %in% family_incomp_no_info, ]
trait_matrix_cp_fam[family_lvl_aggr_final,
                    value := fifelse(is.na(value), i.value, value),
                    on = c("family", "variable")]
trait_matrix_cp[trait_matrix_cp_fam, 
                value := fifelse(is.na(value), i.value, value),
                on = c("unique_id", "variable")]

#### Traits on genus level that could be used for species
# Get species with incomplete or no trait information
# Distinguish between incomplete and no trait information
species_incompl_no_info <- trait_matrix[
  taxonomic_level == "species" &
    trait_info %in% c("no_information", "incomplete"),
  .(taxon, genus)
]
# trait_matrix[taxonomic_level == "genus" &
#                genus %in% unique(species_incompl_no_info$genus) &
#                trait_info == "complete", .(taxon, family, order, trait_info)]

# Those genera with complete trait profiles that can be used are transformed to lf
genus_assign_species <-
  trait_matrix[genus %in% species_incompl_no_info$genus &
                 trait_info == "complete" & 
                 is.na(species),] %>%
  melt(
    .,
    id.vars = c(
      "species",
      "genus",
      "family",
      "order",
      "taxon",
      "taxonomic_level",
      "ind_trait_info",
      "trait_info",
      "type",
      "unique_id",
      "percentile_category_order"
    ),
    variable.name = "trait"
  )

# With this approach, trait profiles from genera are taken
# For those species with incomplete trait profiles, we could also take just trait information when missing
# Caution, merge just works because there is always one species per genus
genus_assign_species[species_incompl_no_info,
                     taxon := i.taxon,
                     on = "genus"]

# Also consider Bob Zuellig's DB (has data on genus level)
setnames(noa_trait_enhanced_lf,
         names(noa_trait_enhanced_lf),
         tolower(names(noa_trait_enhanced_lf)))
setnames(noa_trait_enhanced_lf,
         "otu_new",
         "taxon")

# Reduce can be really handy
# Reduce()
# add <- function(x) Reduce(`+`, x)
# add(list(1, 2, 3, 4, 5, 6))

genus_assign_species <- rbind(genus_assign_species,
                              noa_trait_enhanced_lf[!genus %in% genus_assign_species$genus &
                                                      !is.na(genus),],
                              fill = TRUE)

# Now get the trait matrix with the species with incomplete or no trait information
tm_species_incomplete <- trait_matrix[taxonomic_level == "species" &
                                        trait_info %in% c("no_information", "incomplete"),] %>% melt(
                                          .,
                                          id.vars = c(
                                            "species",
                                            "genus",
                                            "family",
                                            "order",
                                            "taxon",
                                            "taxonomic_level",
                                            "ind_trait_info",
                                            "trait_info",
                                            "type",
                                            "unique_id",
                                            "percentile_category_order"
                                          ),
                                          variable.name = "trait"
                                        )
tm_species_incomplete <- tm_species_incomplete[is.na(value), ]
tm_species_incomplete[genus_assign_species,
                      value := i.value,
                      on = c("genus", "trait")
]

# sum(!is.na(tm_species_incomplete$value))
trait_matrix_cp[tm_species_incomplete,
                value := fifelse(is.na(value), i.value, value),
                on = c("species", variable = "trait")]
# trait_matrix_cp[, uniqueN(taxon)]
# trait_matrix_cp[is.na(value), uniqueN(taxon)]

# Check completeness and try to add trait information for remaining species 
# with family level trait information
species_incompl_no_info2 <- trait_matrix_cp %>%
  dcast(., ... ~ variable, value.var = "value") %>%
  .[, ind_trait_info := apply(.SD, 1, function(y) {
    sum(!is.na(y))
  }),
  .SDcols = trait_names] %>%
  .[, trait_info := fcase(
    ind_trait_info == length(trait_names),
    "complete",
    ind_trait_info != length(trait_names) & ind_trait_info > 0,
    "incomplete",
    ind_trait_info == 0,
    "no_information"
  )] %>% 
  .[trait_info %in% c("incomplete", "no_information") & 
      taxonomic_level == "species", .(family, taxon)]

# Use both assigned or aggregated family level traits
family_assign_species <-
  trait_matrix_cp[family %in% unique(species_incompl_no_info2$family) &
                    taxonomic_level == "family" &
                    trait_info == "complete",]

# Merge back taxon
# problem of duplicates
family_assign_species <-
  merge(
    x = family_assign_species[, .SD, .SDcols = !"taxon"],
    y = species_incompl_no_info2,
    by = "family",
    allow.cartesian = TRUE
  )

# Assign family level traits and aggregated family level traits to some species
tm_species_incomplete2 <-
  trait_matrix_cp[taxon %in% species_incompl_no_info2$taxon &
                    trait_info %in% c("no_information", "incomplete"),]
tm_species_incomplete2 <- tm_species_incomplete2[is.na(value), ]
tm_species_incomplete2[family_assign_species,
                       value := i.value,
                       on = c("taxon", "variable")
]

# Now merge to family level aggregated traits 
tm_species_incomplete2[family_lvl_aggr_final,
                       value := fifelse(is.na(value), i.value, value),
                       on = c("family", "variable")]

# And merge back the family level traits to species
trait_matrix_cp[tm_species_incomplete2,
                value := fifelse(is.na(value), i.value, value),
                on = c("species", "variable")]
# trait_matrix_cp[, uniqueN(taxon)]
# trait_matrix_cp[is.na(value), uniqueN(taxon)]

# Calculate completeness of information at the end (this time in lf)
trait_matrix_cp[, trait_info := fcase(
  sum(!is.na(value)) == length(trait_names),
  "complete",
  sum(!is.na(value)) != length(trait_names) &
    sum(!is.na(value)) > 0,
  "incomplete",
  sum(!is.na(value)) == 0,
  "no_information"
), by = taxon]
trait_matrix_cp <- dcast(trait_matrix_cp, ... ~ variable, value.var = "value")

# Normalize 
# Theoretically necessary - when assigned and aggregated trait information is mixed
normalize_by_rowSum(
  trait_matrix_cp,
  non_trait_cols = c(
    "species",
    "genus",
    "family",
    "order",
    "taxon",
    "taxonomic_level",
    "trait_info",
    "type",
    "unique_id",
    "ind_trait_info",
    "percentile_category_order"
  )
)