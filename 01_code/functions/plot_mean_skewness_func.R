
PlotMeanSkewEffect <- function(
    df,
    title,
    xlab = "Change in mean (%)",
    ylab = "Change in skewness (%)",
    labels = TRUE) {
  
  df <- df |>
    mutate(
      species_plot = str_replace(species, "_", " "),
      species_plot = paste0(substr(word(species_plot, 1), 1, 1), ". ",word(species_plot, 2)),
      effect_size = sqrt(mean^2 + skewness^2),
      #label = if_else(rank(-effect_size) <= 5,species_plot,NA_character_ ), 
      ci_excludes_zero = (mean_lower > 0 | mean_upper < 0) |
                        (skewness_lower > 0 | skewness_upper < 0),
      label = if_else(ci_excludes_zero,
                      species_plot,
                      NA_character_))
  
  p <- ggplot(df, aes(x = mean, y = skewness)) +
    
    # Zero reference lines
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5, 
      colour = "grey20") +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey20") +
    
    # Horizontal 90% CI for mean
    geom_segment(aes(x = mean_lower, xend = mean_upper, y = skewness,yend = skewness),
                linewidth = 0.3,
                colour = "grey50",
                alpha = 0.5) +
    
    # Vertical 90% CI for skewness
    geom_segment(aes(x = mean, xend = mean, y = skewness_lower, yend = skewness_upper),
                linewidth = 0.3,
                colour = "grey50",
                alpha = 0.5) +
    
    # Re-draw CIs in bold when CIs dont cross 0
    geom_segment(
      data = df |> filter(ci_excludes_zero),
      aes(
        x = mean_lower,
        xend = mean_upper,
        y = skewness,
        yend = skewness, 
        colour = fished),
      linewidth = 0.35,
      alpha = 1 ) +

    geom_segment(
      data = df |> filter(ci_excludes_zero),
      aes(
        x = mean,
        xend = mean,
        y = skewness_lower,
        yend = skewness_upper, colour = fished),
      linewidth = 0.35,
      alpha = 1) +
    
    # posterior median
    geom_point(
      aes(size = max_length, alpha = ci_excludes_zero, colour = fished), shape = 16) +
    
    
    scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1), guide = "none")+
    
    scale_colour_manual(values = c("#660154", "#E69F00"))+
    
    scale_size_continuous(range = c(0.9, 2.8),  breaks = c(10, 40, 80, 120))+
    
    labs(title = title,
          x = xlab,
          y = ylab, 
         colour = "Fished", 
         size = "Maximum \nlength (cm)") +
    
    theme_classic()+
    theme(axis.title = element_text(size = 11.5, lineheight = 1.2), 
          axis.text = element_text(size = 9), 
          legend.title = element_text(size = 13, lineheight = 1.05), 
          legend.text = element_text(size = 12), 
          legend.key.height = unit(0.8, "cm"),
          axis.line = element_line(linewidth = 0.6), 
          axis.ticks = element_line(linewidth = 0.6))
  
  if (labels) {
    
    p <- p +
      # Only label selected species
      geom_text_repel(
        aes(label = label),
        size = 1.4,
        max.overlaps = Inf,
        box.padding = 0.5,
        point.padding = 0.3,
        min.segment.length = 0,
        segment.size = 0.2,
        seed = 123
      ) 
  }
  
  p
}