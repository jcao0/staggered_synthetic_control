# Function to plot RMSE for a given case
plot_rmse_case <- function(plot_data, case_name, save_file = NULL) {
    # Define colors, line types, and point shapes for each method
    
    # ---- Subset the data for this case ----
    df <- subset(plot_data, case == case_name)
    T <- as.numeric(sub(".*T([0-9]+).*", "\\1", case_name))
    N <- as.numeric(sub("N([0-9]+)_.*", "\\1", case_name))
    r <- as.numeric(sub(".*_r([0-9]+)", "\\1", case_name))
    # n_sims <- unique(df$n_sims)
    
    # ---- Set method order ----
    present_methods <- unique(df$method)
    method_order <- c()
    if ("ssc" %in% present_methods) method_order <- c(method_order, "ssc")
    method_order <- c(method_order, setdiff(present_methods, "ssc"))
    df$method <- factor(df$method, levels = method_order)
    
    # ---- Colors, line types, shapes ----
    method_colors    <- c("ssc" = rgb(1,0.4,0.3), 
                            "gsc" = rgb(0.55,0.45,0.75), 
                            "asy" = rgb(0,0.6,0.6))
    method_linetypes <- c("ssc" = "solid", "gsc" = "dashed", "asy" = "dotdash")
    method_shapes    <- c("ssc" = 16, "gsc" = 17, "asy" = 15)
    
    # Capitalized legend labels
    method_labels <- toupper(levels(df$method))

    p <- ggplot(
        df,
        aes(x = event_time, y = rmse, color = method, linetype = method, shape = method)
    ) +
        geom_hline(yintercept = 0, linetype = "longdash", linewidth = 0.6, color = "gray40") +
        geom_line(na.rm = TRUE, linewidth = 1) +
        geom_point(na.rm = TRUE, size = 3) +
        scale_color_manual(values = method_colors, labels = method_labels) +
        scale_linetype_manual(values = method_linetypes, labels = method_labels) +
        scale_shape_manual(values = method_shapes, labels = method_labels) +
        scale_x_continuous( 
        breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)
        ) +
        labs(
        # title = paste("T =", T, ", N =", N, ", r =", r),
        x = "Event Time",
        y = "Root Mean Squared Error (RMSE)",
        color = NULL,
        linetype = NULL,
        shape = NULL
        ) +
        theme_bw() +
        theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        legend.position = c(0.98, 0.02),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = NA, color = NA),
        legend.key.width = unit(1.5, "cm"), # wider legend keys
        legend.key.height = unit(0.5, "cm")
        )
    
    # Save plot if file path is provided
    if (!is.null(save_file)) {
        ggsave(filename = save_file, plot = p, width = 10, height = 10, dpi = 100)
    }
    
    return(p)
}

# Function to combine multiple plots into a single figure
combine_plots <- function(plot_list,
                          panel_titles = NULL,
                          shared_legend = FALSE,
                          save_file = NULL) {
  
  n <- length(plot_list)
  parameter_panel_titles <- panel_titles
  

  # ---- Default panel captions ----
  if (is.null(parameter_panel_titles)) {
    labs <- names(plot_list)
    panel_titles <- sapply(labs, function(case_name) {
                T <- as.numeric(sub(".*T([0-9]+).*", "\\1", case_name))
                N <- as.numeric(sub("N([0-9]+)_.*", "\\1", case_name))
                r <- as.numeric(sub(".*_r([0-9]+)", "\\1", case_name))    
   })
   panel_titles <- sapply(seq_len(n), function(i) {
    paste0( "(" , letters[i], ") ",
            "T = ", as.numeric(sub(".*T([0-9]+).*", "\\1", names(plot_list)[i])), 
            ", N = ", as.numeric(sub("N([0-9]+)_.*", "\\1", names(plot_list)[i])),
            ", r = ", as.numeric(sub(".*_r([0-9]+)", "\\1", names(plot_list)[i])))
   })
  }

  # --- Align y-axes ----
  ymax <- max(unlist(lapply(plot_list, function(p) {
    ggplot_build(p)$data[[2]]$y
  })),na.rm = TRUE)
  plot_list <- lapply(plot_list, function(p) {
    p + scale_y_continuous(limits = c(0, ymax + 0.1*ymax))
  })

  # ---- Attach captions ----
  plots_labeled <- lapply(seq_len(n), function(i) {
    plot_list[[i]] +
      labs(caption = panel_titles[i])
  })
  
  # ---- Combine horizontally for each r----
  if (is.null(parameter_panel_titles)) {
    labs <- names(plot_list)
    r_list <- sapply(labs, function(case_name) {
                r <- as.numeric(sub(".*_r([0-9]+)", "\\1", case_name))
    })
    n_r <- length(unique(r_list))
    combined <- wrap_plots(plots_labeled, nrow = n_r)
  } else {
    combined <- wrap_plots(plots_labeled, nrow = 1)
  }
  
  # ---- Shared legend (optional) ----
  if (shared_legend) {
    combined <- combined + plot_layout(guides = "collect") + theme(legend.position = "top")
  }
  
  # ---- Global styling ----
  combined <- combined &
    theme(
      plot.caption = element_text(
        hjust = 0.5,
        vjust = -1,
        size  = 12
      ),
      plot.margin = margin(5.5, 5.5, 15, 5.5),
    )
  
  # ---- Save plot if file path is provided ----
    if (!is.null(save_file)) {
        # ggsave(filename = save_file, plot = combined, width = 10, height = 10, dpi = 100)
        ggsave(filename = save_file, plot = combined, width = 12, height = 8, dpi = 100)
    }

  return(combined)
}
