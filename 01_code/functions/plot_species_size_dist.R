
PlotSpeciesDistribution <- function(
    prediction,
    show_uncertainty = FALSE,
    cp = "#660154") {
    
    curve_df <- prediction$curve
    
    species_label <- gsub("_", " ", prediction$species)
   
    # -----------------------------------------
    # Panel 1: time
    # -----------------------------------------
    
    time_df <- curve_df |>
      filter(category == "Time") |>
      mutate(
        year = factor(
          year,
          levels = c(1990, 2000, 2010, 2020)),
        MPA_label = factor(
          MPA_label,
          levels = c("Non-MPA", "MPA")))
    
    p_time <- ggplot(
      time_df,
      aes(
        x = length,
        y = density,
        colour = year,  
        linetype = MPA_label,
        group = interaction(year, MPA_label)))
    
    if (show_uncertainty) {
      p_time <- p_time +
        geom_ribbon( data = time_df,
          aes( x= length,
            ymin = density_lower,
            ymax = density_upper,
            fill = year, 
            group = interaction(year, MPA_label)),
          alpha = 0.08,
          colour = NA,
          inherit.aes = FALSE )
    }
    
    p_time <- p_time +
      geom_line(linewidth = 0.7) +
      scale_x_continuous(
        limits = c(0, prediction$x_max),
        expand = expansion(mult = c(0, 0.01)) ) +
      scale_colour_manual(values = c("#191970","#56B4E9","#E69F00","#CC79A7"))+
      scale_fill_manual(values = c("#191970","#56B4E9","#E69F00","#CC79A7"))+
      scale_linetype_manual(
        values = c("Non-MPA" = "dashed", "MPA" = "solid") )+
      labs(
        x = "Length (cm)",
        y = "Probability density",
        colour = "Year",
        fill = "Year",
        linetype = NULL ) +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "right", 
        axis.title = element_text(size = 14, lineheight = 1.2), 
        axis.text = element_text(size = 13), 
        legend.title = element_text(size = 14, lineheight = 1.05), 
        legend.text = element_text(size = 13))
    
    
    # -----------------------------------------
    # Panel 2: temperature
    # -----------------------------------------
    
    temp_df <- curve_df |>
      filter(category == "Temperature") 
    
    temp_levels <- temp_df |>
      distinct(temp_label, scenario_order) |>
      arrange(scenario_order) |>
      pull(temp_label)
    
    temp_df <- temp_df |>
      mutate(
        temp_label = factor(
          temp_label,
          levels = temp_levels ))
    
    p_temp <- ggplot(
      temp_df,
      aes(
        x = length,
        y = density,
        colour = temp_label,
        group = temp_label ) )
    
    if (show_uncertainty) {
      p_temp <- p_temp +
        geom_ribbon(data = temp_df,
          aes( x = length, 
            ymin = density_lower,
            ymax = density_upper,
            fill = temp_label , 
             group = temp_label),
          alpha = 0.10,
          colour = NA, 
          inherit.aes = FALSE)
    }
    
    p_temp <- p_temp +
      geom_line(linewidth = 0.9) +
      scale_x_continuous(
        limits = c(0, prediction$x_max),
        expand = expansion(mult = c(0, 0.01))) +
      scale_colour_manual(values = c("#633372FF", "#1F6E9CFF","#F4C40FFF","#D8443CFF"))+
      scale_fill_manual(values = c("#633372FF", "#1F6E9CFF","#F4C40FFF","#D8443CFF"))+
      labs(
        x = "Length (cm)",
        y = "Probability density",
        colour = "Temperature", 
        fill = "Temperature") +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "right", 
        axis.title = element_text(size = 14, lineheight = 1.2), 
        axis.text = element_text(size = 13), 
        legend.title = element_text(size = 14, lineheight = 1.05), 
        legend.text = element_text(size = 13))

    # -----------------------------------------
    # combine the two panels
    # -----------------------------------------
    
    combined_plot <-
      (p_time / p_temp) +
      patchwork::plot_annotation(tag_levels = "A", title = species_label)&theme(plot.title = element_text(face = "italic", size = 16))
    
    combined_plot
  }