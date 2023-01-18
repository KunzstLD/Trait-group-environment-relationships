# ______________________________________________________________________
# Collection of functions used during data processing & analysis 
# ______________________________________________________________________

# Check for completeness of trait dataset ----
create_pattern_ind <- function(x, non_trait_cols) {
  if (missing(non_trait_cols)) {
    trait_names_pattern <- sub("\\_.*|\\..*", "", names(x)) %>%
      unique() %>%
      paste0("^", .)
  } else{
    pat <- paste0(non_trait_cols, collapse = "|")
    # get trait names & create pattern for subset
    trait_names_pattern <-
      grep(pat, names(x), value = TRUE, invert = TRUE) %>%
      sub("\\_.*|\\..*", "", .) %>%
      unique() %>%
      paste0("^", .)
  }
  trait_names_pattern
}

completeness_trait_data <- function(x, non_trait_cols) {
  trait_names_pattern <- create_pattern_ind(
    x = x,
    non_trait_cols = non_trait_cols
  )
  
  # test how complete trait sets are
  output <- matrix(ncol = 2, nrow = length(trait_names_pattern))
  for (i in seq_along(trait_names_pattern)) {
    # vector containing either 0 (no NAs) or a number (> 0) meaning that all
    # entries for this trait contained NA
    vec <-
      x[, apply(.SD, 1, function(y) {
        base::sum(is.na(y))
      }),
      .SDcols = names(x) %like% trait_names_pattern[[i]]
      ]
    
    # How complete is the dataset for each individual trait?
    output[i, ] <-
      c(
        (length(vec[vec == 0]) / nrow(x)) %>% `*`(100) %>% round(),
        trait_names_pattern[[i]]
      )
  }
  output <- as.data.table(output)
  return(output)
}

# Normalization of trait scores ----
# All trait states of one trait are divided by their row sum
# Hence, trait affinities are represented as "%" or ratios
# this works for traits named: "groupingfeature_trait"
normalize_by_rowSum <- function(x, 
                                non_trait_cols, 
                                na.rm = TRUE) {
  # get trait names & create pattern for subset
  trait_names_pattern <- create_pattern_ind(x = x,
                                            non_trait_cols = non_trait_cols)
  
  # loop for normalization (trait categories for each trait sum up to 1)
  for (cols in trait_names_pattern) {
    # get row sum for a specific trait
    x[, rowSum := apply(.SD, 1, sum, na.rm = na.rm),
      .SDcols = names(x) %like% cols]
    
    # get column names for assignment
    col_name <- names(x)[names(x) %like% cols]
    
    # divide values for each trait state by
    # the sum of trait state values
    x[, (col_name) := lapply(.SD, function(y) {
      y / rowSum
    }),
    .SDcols = names(x) %like% cols]
  }
  # del rowSum column
  x[, rowSum := NULL]
  return(x)
}


# Trait aggregation ----
# Direct aggregation to family-level
# Works for data.tables
direct_agg <- function(trait_data,
                       non_trait_cols,
                       method,
                       na.rm = TRUE,
                       on = "family") {
  if (!on %in% c("genus", "family")) {
    warning(paste("Aggregation to", on, "is not implemented yet"))
  }
  # get names of trait columns
  pat <- paste0(non_trait_cols, collapse = "|")
  trait_col <-
    grep(pat, names(trait_data), value = TRUE, invert = TRUE)
  
  # aggregate to family-level
  # & merge information on order back
  # subset so that no NA values occur in data
  # (otherwise all NA entries are viewed as a group &
  # aggregated as well)
  if (on == "genus") {
    agg_data <- trait_data[!is.na(genus),
                           lapply(.SD, method, na.rm = na.rm),
                           .SDcols = trait_col,
                           by = on]
    agg_data[trait_data, `:=`(family = i.family,
                              order = i.order), on = on]
    return(agg_data)
  }
  if (on == "family") {
    agg_data <- trait_data[!is.na(family),
                           lapply(.SD, method, na.rm = na.rm),
                           .SDcols = trait_col,
                           by = on]
    agg_data[trait_data, order := i.order, on = on]
    return(agg_data)
  }
}

# Statistical (helper) functions ----

# test <- rnorm(n = 10000, mean = 2, sd = 0.1)
# Problem of overflow, will result in Inf for large lists
# gmean <- function(x){
#   (Reduce("*", x))^(1/length(x))
# }
# TODO: Understand the relation to the logarithm
gmean <- function(x) {
  exp(mean(log(x)))
}

## Cluster analysis ----
mycluster_hc <- function(x, k) {
  list(cluster = cutree(hclust(as.dist(x),
                               method = "ward.D2"),
                        k = k))
}

## Community weighted mean trait ----
# relatively slow, could improve the loop?
# TODO: transform abundances?
calc_cwm <- function(abund,
                     trait,
                     trait_names) {
  taxa_names <- names(abund)[names(abund) != c("site")]
  cwm <- list()
  for (i in trait_names) {
    trait_sub <- trait[, .SD, .SDcols = c("taxon", i)]
    trait_sub <- trait_sub[match(taxa_names, taxon), ]
    
    cwm[[i]] <-
      abund[, apply(.SD, 1, function(x)
        weighted.mean(trait_sub[[i]],
                      w = x)),
        .SDcols = !c("site"),
        by = "site"]
  }
  cwm
}

# General helper functions ----

## Google search ----
# For google search from Rstudio
query_google <- function(x) {
  invisible(lapply(x,
                   function(y) {
                     utils::browseURL(url = paste0("https://google.com/search?q=", y))
                   }))
}

## Load batches of similarly named datasets ----
# So far, only RDS files
load_data <- function(path, pattern, name_rm_pattern = NULL) {
  files <- list.files(path = path, pattern = pattern)
  data <- lapply(files, function(y) readRDS(file = file.path(path, y)))
  name <- sub("\\.rds", "", files)
  if (!is.null(name_rm_pattern)) {
    name <- sub(name_rm_pattern, "", name)
  }
  data <- setNames(data, name)
  data
}

## Removing strange notations ----
# E.g., in taxonomic nomenclature
rm_strange_notations <- function(x){
  sub("\\?| group", "", x)
}

# Plotting ----
## Dendrograms ----
fun_dendrog_pl <- function(hc,
                           optimal_nog,
                           labels,
                           hang_height = 0.001) {
  hc %>%
    as.dendrogram() %>%
    color_branches(k = optimal_nog) %>%
    hang.dendrogram(hang_height = hang_height) %>%
    set("labels_cex", 0.55) %>%
    dendextend::ladderize() %>%
    set("labels", labels)
}
