# Shared plotting conventions ---------------------------------------------

# Species colours follow the aphids themselves: pea aphids are pink, bean
# aphids are black. Because bean is near-black, avoid dark-grey ink for other
# plot elements (fitted curves, reference lines); use the species colour with
# a different linetype instead.
species_colours <- c(pea = "#d55181", bean = "#262626")

species_labels <- c(pea = "Pea aphid", bean = "Bean aphid")

theme_model <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey35"),
      legend.position = "bottom",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}
