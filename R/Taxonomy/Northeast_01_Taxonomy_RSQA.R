# __________________________________
# Data processing MSQA (midwest)
# __________________________________
northeast_qa <- fread(file.path(path_in, "NoA", "NE_NoAmbig_Species_Matrix n92_Primer.csv"))
setDT(northeast_qa)

# 92 sites, 298 taxa
# dim(northeast_qa)
# names(northeast_qa)
# str(northeast_qa)

# ____________________________________________________________________
## Establishing taxonomy ----
# ____________________________________________________________________
taxa_names <-
    names(northeast_qa)[names(northeast_qa) != "TSTAID"]
northeast_taxa <- as.data.table(taxa_names)
setnames(northeast_taxa, "taxa_names", "taxon")

northeast_taxa[taxon == "Conchapelopia (Helopelopia)", taxon := "Conchapelopia"] 
northeast_taxa[taxon == "Cricotopus (Isocladius)", taxon := "Cricotopus"] 

# Add Taxonomy from NOA traits
northeast_taxa[noa_traits,
    `:=`(
        species = i.species,
        genus = i.genus,
        family = i.family,
        order = i.order
    ),
    on = "taxon"
]

# Use Taxonomy from Midwest and CA
northeast_taxa[is.na(species) & is.na(genus) & is.na(family) & is.na(order), ]
midwest_taxa <- readRDS(file.path(path_cache, "midwest_taxonomy.rds"))
ca_taxa <- readRDS(file.path(path_cache, "ca_taxonomy.rds"))

northeast_taxa[midwest_taxa, `:=`(
    species = i.species,
    genus = i.genus,
    tribe = i.tribe,
    subfamily = i.subfamily,
    family = i.family,
    suborder = i.suborder,
    order = i.order,
    subclass = i.subclass,
    class = i.class,
    phylum = i.phylum
),
on = "taxon"
]
northeast_taxa[ca_taxa, `:=`(
    species = i.species,
    genus = i.genus,
    subfamily = i.subfamily,
    family = i.family,
    order = i.order,
    class = i.class,
    phylum = i.phylum
),
on = "taxon"
]


# Use taxize to add the remaining taxa
ID_gbif <-
    get_gbifid(northeast_taxa[
        is.na(species) & is.na(genus) & is.na(family) & is.na(order) & 
        is.na(tribe) & is.na(subfamily) & is.na(suborder) & is.na(subclass) &
        is.na(class) & is.na(phylum),
        taxon
    ])
# saveRDS(ID_gbif, file = file.path(data_cache, "northeast_ID_gbif.rds"))

ID_merge <-
    data.table(
        "taxon" = northeast_taxa[
            is.na(species) & is.na(genus) & is.na(family) & is.na(order) &
                is.na(tribe) & is.na(subfamily) & is.na(suborder) & is.na(subclass) &
                is.na(class) & is.na(phylum),
            taxon
        ],
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
northeast_taxa[tax_classf,
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

northeast_taxa[
            is.na(species) & is.na(genus) & is.na(family) & is.na(order) &
                is.na(tribe) & is.na(subfamily) & is.na(suborder) & is.na(subclass) &
                is.na(class) & is.na(phylum),
            taxon
        ]

# Two taxa remain
northeast_taxa[taxon == "Oligochaeta", `:=`(
    subclass = "Oligochaeta",
    class = "Clitellata",
    phylum = "Annelida"
)]
northeast_taxa[taxon == "Diamesinae", `:=`(
    subfamily = "Diamesinae",
    family = "Chironomidae",
    order = "Diptera"
)]

# Add taxonomic level
northeast_taxa[, taxonomic_level := fcase(
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
    northeast_taxa,
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

# Save
saveRDS(northeast_taxa, file = file.path(path_cache, "northeast_taxonomy.rds"))