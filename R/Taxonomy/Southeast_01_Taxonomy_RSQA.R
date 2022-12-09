# __________________________________
# Data processing PN
# __________________________________
southeast_qa <- fread(file.path(data_in, "NoA", "SESQA_Inverts_transp_TSTAID_n75_clean.csv"))
setDT(southeast_qa)

# 75 sites, 282 taxa
# dim(southeast_qa)
# names(southeast_qa)
# str(southeast_qa)

# ____________________________________________________________________
## Establishing taxonomy ----
# ____________________________________________________________________
taxa_names <-
    names(southeast_qa)[!names(southeast_qa) %in% c(
        "TSTAID",
        "SHORT_NAME",
        "SiteType2USE",
        "UrbanCenter_tier"
    )]
southeast_taxa <- as.data.table(taxa_names)
setnames(southeast_taxa, "taxa_names", "taxon")

# Add Taxonomy from NOA traits
southeast_taxa[noa_traits,
    `:=`(
        species = i.species,
        genus = i.genus,
        family = i.family,
        order = i.order
    ),
    on = "taxon"
]

# Use Taxonomy from other datasets
# southeast_taxa[is.na(species) & is.na(genus) & is.na(family) & is.na(order), ]
midwest_taxa <- readRDS(file.path(data_cache, "midwest_taxonomy.rds"))
ca_taxa <- readRDS(file.path(data_cache, "ca_taxonomy.rds"))
northeast_taxa <- readRDS(file.path(data_cache, "northeast_taxonomy.rds"))
pn_taxa <- readRDS(file.path(data_cache, "pn_taxonomy.rds"))

southeast_taxa[midwest_taxa, `:=`(
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
southeast_taxa[ca_taxa, `:=`(
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
southeast_taxa[northeast_taxa, `:=`(
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
southeast_taxa[pn_taxa, `:=`(
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

# Use taxize to add the remaining taxa
ID_gbif <-
    get_gbifid(southeast_taxa[
        is.na(species) & is.na(genus) & is.na(family) & is.na(order) & 
        is.na(tribe) & is.na(subfamily) & is.na(suborder) & is.na(subclass) &
        is.na(class) & is.na(phylum),
        taxon
    ])
# saveRDS(ID_gbif, file = file.path(data_cache, "southeast_ID_gbif.rds"))

ID_merge <-
    data.table(
        "taxon" = southeast_taxa[
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
southeast_taxa[tax_classf,
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

# Add taxonomic level
southeast_taxa[, taxonomic_level := fcase(
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
    southeast_taxa,
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
saveRDS(southeast_taxa, file = file.path(data_cache, "southeast_taxonomy.rds"))