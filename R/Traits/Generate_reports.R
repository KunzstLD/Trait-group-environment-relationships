# https://bookdown.org/yihui/rmarkdown-cookbook/parameterized-reports.html
region_names <- c(
    "California",
    "Midwest",
    "Northeast",
    "PN",
    "Southeast"
)

# For loop for generating R markdown reports
for (region in region_names) {
  rmarkdown::render(
    file.path(path_scr, "Traits", "Template.Rmd"),
    params = list(region = region),
    output_file = paste0(region, ".html"),
    output_dir = "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Taxonomy_trait_information_overview"
  )
}
