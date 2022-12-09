# __________________________________
# Data processing PN
# __________________________________
pn_qa <- fread(file.path(path_in, "NoA", "PNSQA Species Matrix_Primer.csv"))
setDT(pn_qa)

# 87 sites, 263 taxa
# dim(pn_qa)
# names(pn_qa)
# str(pn_qa)

# ____________________________________________________________________
## Establishing taxonomy ----
# ____________________________________________________________________
taxa_names <-
    names(pn_qa)[names(pn_qa) != "TSTAID"]
pn_taxa <- as.data.table(taxa_names)
setnames(pn_taxa, "taxa_names", "taxon")

# Taxa with "/", "group" and "complex" (potential subspecies?) cleaned for now
pn_taxa[, taxon := sub("\\/.+", "", taxon)]
pn_taxa[, taxon := sub(" group", "", taxon)]
pn_taxa[, taxon := sub(" complex", "", taxon)]

# Add Taxonomy from NOA traits
pn_taxa[noa_traits,
    `:=`(
        species = i.species,
        genus = i.genus,
        family = i.family,
        order = i.order
    ),
    on = "taxon"
]

# Use Taxonomy from Midwest and CA
# pn_taxa[is.na(species) & is.na(genus) & is.na(family) & is.na(order), ]
midwest_taxa <- readRDS(file.path(data_cache, "midwest_taxonomy.rds"))
ca_taxa <- readRDS(file.path(data_cache, "ca_taxonomy.rds"))
northeast_taxa <- readRDS(file.path(data_cache, "northeast_taxonomy.rds"))

pn_taxa[midwest_taxa, `:=`(
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
pn_taxa[ca_taxa, `:=`(
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
pn_taxa[northeast_taxonomy, `:=`(
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
    get_gbifid(pn_taxa[
        is.na(species) & is.na(genus) & is.na(family) & is.na(order) & 
        is.na(tribe) & is.na(subfamily) & is.na(suborder) & is.na(subclass) &
        is.na(class) & is.na(phylum),
        taxon
    ])
# saveRDS(ID_gbif, file = file.path(data_cache, "pn_ID_gbif.rds"))

ID_merge <-
    data.table(
        "taxon" = pn_taxa[
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
pn_taxa[tax_classf,
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

# One taxa remains
pn_taxa[
    taxon == "Tanypodinae",
    `:=`(
        subfamily = "Tanypodinae",
        family = "Chironomidae",
        order = "Diptera"
    )
]

# Add taxonomic level
pn_taxa[, taxonomic_level := fcase(
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
    pn_taxa,
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

# Few duplicates that will be removed
pn_taxa <- pn_taxa[!duplicated(taxon), ]

# Save
saveRDS(pn_taxa, file = file.path(data_cache, "pn_taxonomy.rds"))
