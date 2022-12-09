# Creates plots for visualising data gaps in trait information
# Script is meant to be run within the main RMD file
trait_matrix_lf <- melt(
  trait_matrix,
  id.vars = c(
    "species",
    "genus",
    "family",
    "order",
    "taxon",
    "taxonomic_level",
    "ind_trait_info",
    "trait_info",
    "type"
  ),
  variable.name = "trait",
  value.name = "affinity"
)
trait_matrix_lf[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)]
trait_matrix_lf[, value := fcase(is.na(affinity),
                              "missing",
                              !is.na(affinity),
                              "available")]
trait_matrix_lf[, value := as.factor(value)]

# Plot
missing_traits_pl1 <- trait_matrix_lf[trait_info %in% c("incomplete", "no_information"), ] %>%
  .[order %in% aq_insects[1:5], ] %>%
  ggplot(
    .,
    aes(
      x = taxon,
      y = trait,
      fill = value
    )
  ) +
  geom_tile() +
  scale_fill_manual(values = c(
    "missing" = "grey",
    "available" = "steelblue"
  )) +
  facet_grid(grouping_feature ~ order,
    scales = "free",
    space = "free"
  ) +
  labs(
    x = "Family",
    y = "Trait",
    fill = "Trait value"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(
    #   family = "Roboto Mono",
    #   size = 16,
    #   angle = 90,
    #   hjust = 1,
    #   vjust = 0.2
    # ),
    axis.text.y = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    legend.title = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    legend.text = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    strip.text = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    plot.title = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    panel.grid = element_blank()
  )

missing_traits_pl2 <- trait_matrix_lf[trait_info %in% c("incomplete", "no_information"), ] %>%
  .[order %in% aq_insects[6:11], ] %>%
  ggplot(
    .,
    aes(
      x = taxon,
      y = trait,
      fill = value
    )
  ) +
  geom_tile() +
  scale_fill_manual(values = c(
    "missing" = "grey",
    "available" = "steelblue"
  )) +
  facet_grid(grouping_feature ~ order,
    scales = "free",
    space = "free"
  ) +
  labs(
    x = "Family",
    y = "Trait",
    fill = "Trait value"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(
    #   family = "Roboto Mono",
    #   size = 16,
    #   angle = 90,
    #   hjust = 1,
    #   vjust = 0.2
    # ),
    axis.text.y = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    legend.title = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    legend.text = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    strip.text = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    plot.title = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    panel.grid = element_blank()
  )

missing_traits_pl_non_insects <- trait_matrix_lf[trait_info %in% c("incomplete", "no_information") &
  type == "non_insect", ] %>%
  ggplot(
    .,
    aes(
      x = taxon,
      y = trait,
      fill = value
    )
  ) +
  geom_tile() +
  scale_fill_manual(values = c(
    "missing" = "grey",
    "available" = "steelblue"
  )) +
  facet_grid(grouping_feature ~ order,
    scales = "free",
    space = "free"
  ) +
  labs(
    x = "Family",
    y = "Trait",
    fill = "Trait value"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(
    #   family = "Roboto Mono",
    #   size = 16,
    #   angle = 90,
    #   hjust = 1,
    #   vjust = 0.2
    # ),
    axis.text.y = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    legend.title = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    legend.text = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    strip.text = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    plot.title = element_text(
      family = "Roboto Mono",
      size = 13
    ),
    panel.grid = element_blank()
  )