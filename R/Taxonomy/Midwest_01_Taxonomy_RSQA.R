# __________________________________
# Data processing MSQA (midwest)
# __________________________________
midwest_qa <- read_xlsx(file.path(path_in, "NoA", "MSQA_Opt4_Comb_MV.xlsx"))
setDT(midwest_qa)

# 102 sites, 317 taxa
# dim(midwest_qa)
# names(midwest_qa)


# ____________________________________________________________________
## Establishing taxonomy ----
# ____________________________________________________________________
taxa_names <-
    names(midwest_qa)[!names(midwest_qa) %in% c(
        "SCODE",
        "SUID",
        "STAID",
        "REACH",
        "DATE",
        "SampleID",
        "SMCOD"
    )]
midwest_taxa <- as.data.table(taxa_names)
setnames(midwest_taxa, "taxa_names", "taxon")

# Save those taxa (potential duplicates)
saveRDS(midwest_taxa[grepl("\\/.+", taxon), ],
    file = file.path(path_cache, "potential_duplicates_midwest.rds")
)

# Cleaning taxa names (ask I.Waite)
# - "(...)"
# - 3 names 
# "/" -> for now remove the second name
midwest_taxa[, taxon := sub(" \\(.+\\)", "", taxon)]
midwest_taxa[, taxon := sub("\\/.+", "", taxon)]
midwest_taxa[, taxon := sub("(.+)(\\s)(.+)(\\s)(.+)", "\\1\\2\\3", taxon)]

midwest_taxa[noa_traits,
    `:=`(
        species = i.species,
        genus = i.genus,
        family = i.family,
        order = i.order
    ),
    on = "taxon"
]

# Use taxize to add the remaining 71 
ID_gbif <-
    get_gbifid(midwest_taxa[
        is.na(species) & is.na(genus) & is.na(family) & is.na(order),
        taxon
    ])
# saveRDS(ID_gbif, file = file.path(data_cache, "midwest_ID_gbif.rds"))

ID_merge <-
    data.table(
        "taxon" = midwest_taxa[is.na(species) &
            is.na(genus) &
            is.na(family) &
            is.na(order), taxon],
        "query" = ID_gbif
    )
tax_classf <- classification(
  id = ID_gbif,
  db = "gbif"
)
tax_classf <- rbind(tax_classf)
setDT(tax_classf)

# Merge back original names of the queried taxa
tax_classf[ID_merge,
    merge_taxa := i.taxon,
    on = "query"
]
tax_classf <- dcast(tax_classf, query + merge_taxa ~ rank, value.var = "name")

# Merge to CA taxonomy
midwest_taxa[tax_classf,
    `:=`(
        species = i.species,
        genus = i.genus,
        family = i.family,
        order = i.order,
        class = i.class,
        phylum = i.phylum
    ),
    on = c(taxon = "merge_taxa")
]

# Add taxonomy for remaining taxa by hand
midwest_taxa[taxon == "Hirudinea", `:=`(
    subclass = "Hirudinea",
    class = "Clitellata",
    phylum = "Hirudinea"
)]
midwest_taxa[taxon == "Tanytarsini", `:=`(
    tribe = "Tanytarsini",
    subfamily = "Chironominae",
    family = "Chironomidae",
    order = "Diptera"
)]
midwest_taxa[taxon == "Orthocladiinae", `:=`(
    subfamily = "Orthocladiinae",
    family = "Chironomidae",
    order = "Diptera"
)]
midwest_taxa[taxon == "Megadrile", `:=`(
    order = "Megadrile",
    class = "Clitellata",
    phylum = "Annelida"
)]
midwest_taxa[taxon == "Acari", `:=`(
    subclass = "Acari",
    class = "Arachnida", 
    phylum = "Arthropoda"
)]
midwest_taxa[taxon == "Anisoptera", `:=`(
    suborder = "Anisoptera",
    order = "Odonata"
)]
midwest_taxa[taxon == "Pentaneurini", `:=`(
    tribe = "Pentaneurini",
    subfamily = "Chironominae",
    family = "Chironomidae",
    order = "Diptera"
)]
midwest_taxa[taxon == "Thienemannimyia group", `:=`(
    genus = "Thienemannimyia",
    tribe = "Pentaneurini",
    subfamily = "Chironominae",
    family = "Chironomidae",
    order = "Diptera"
)]
midwest_taxa[taxon == "Ceratopogoninae", `:=`(
    subfamily = "Ceratopogoninae",
    family = "Ceratopogonidae",
    order = "Diptera"
)]
midwest_taxa[taxon == "Oligochaeta", `:=`(
    subclass = "Oligochaeta",
    class = "Clitellata",
    phylum = "Annelida"
)]
midwest_taxa[family == "Naididae", order := "Haplotaxida"]

# Add taxonomic level
midwest_taxa[, taxonomic_level := fcase(
    !is.na(species),
    "species",
    is.na(species) & !is.na(genus),
    "genus",
    is.na(species) & is.na(genus) & !is.na(tribe),
    "tribe",
    is.na(species) & is.na(genus) & is.na(tribe) &
        !is.na(subfamily),
    "subfamily",
    is.na(species) & is.na(genus) & is.na(tribe) &
        is.na(subfamily) &
        !is.na(family),
    "family",
    is.na(species) & is.na(genus) & is.na(tribe) &
        is.na(subfamily) & is.na(family) & !is.na(suborder),
    "suborder",
    is.na(species) & is.na(genus) & is.na(tribe) &
        is.na(subfamily) & is.na(family) & is.na(suborder) &
        !is.na(order),
    "order",
    is.na(species) & is.na(genus) & is.na(tribe) &
        is.na(subfamily) & is.na(family) & is.na(suborder) &
        is.na(order) & !is.na(subclass),
    "subclass",
    is.na(species) & is.na(genus) & is.na(tribe) &
        is.na(subfamily) & is.na(family) & is.na(suborder) &
        is.na(order) & is.na(subclass) & !is.na(class),
    "class",
    is.na(species) & is.na(genus) & is.na(tribe) &
        is.na(subfamily) & is.na(family) & is.na(suborder) &
        is.na(order) & is.na(subclass) & is.na(class) & !is.na(phylum),
    "phylum"
)]

# Change column order
setcolorder(
  midwest_taxa,
  c(
    "taxon",
    "species",
    "genus",
    "tribe",
    "subfamily",
    "family",
    "suborder",
    "order",
    "subclass",
    "class",
    "phylum",
    "taxonomic_level"
  )
)


# Remove duplicates (for now)
midwest_taxa <- midwest_taxa[!duplicated(taxon), ]

# Save taxonomy
saveRDS(midwest_taxa, file = file.path(path_cache, "midwest_taxonomy.rds"))