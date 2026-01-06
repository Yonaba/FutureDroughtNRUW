root_folder <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir=root_folder)
Sys.setenv(TZ = "UTC")

library(ggplot2)
library(scales)
library(ggpubr)

obs_df <- read.csv("data/csv/ref_climate_data.csv", stringsAsFactors = FALSE)
raw_df <- read.csv("data/csv/raw_climate_model_data.csv", stringsAsFactors = FALSE)
bc_df  <- read.csv("data/csv/bc_climate_model_data.csv", stringsAsFactors = FALSE)

obs_df$date <- as.Date(obs_df$date, format = "%m/%d/%Y")
raw_df$date <- as.Date(raw_df$date)
bc_df$date  <- as.Date(bc_df$date)

# ---- Filter Historical Period (1981–2010) ----
obs_df <- obs_df[obs_df$date >= as.Date("1981-01-01") & obs_df$date <= as.Date("2014-12-31"), ]
raw_df <- raw_df[raw_df$scenario == "historical" & raw_df$date >= as.Date("1981-01-01") & raw_df$date <= as.Date("2014-12-31"), ]
bc_df  <- bc_df[bc_df$scenario == "historical" & bc_df$date >= as.Date("1981-01-01") & bc_df$date <= as.Date("2014-12-31"), ]

month_levels <- month.abb  # Jan, Feb, ..., Dec
obs_df$month <- factor(month.abb[as.integer(format(obs_df$date, "%m"))], levels = month_levels, ordered = TRUE)
raw_df$month <- factor(month.abb[as.integer(format(raw_df$date, "%m"))], levels = month_levels, ordered = TRUE)
bc_df$month  <- factor(month.abb[as.integer(format(bc_df$date, "%m"))], levels = month_levels, ordered = TRUE)

# Observed: melt into long format
obs_long <- reshape(obs_df[, c("station", "date", "month", "pr", "tasmin", "tasmax", "pet")],
                    direction = "long",
                    varying = list(c("pr", "tasmin", "tasmax", "pet")),
                    v.names = "value",
                    timevar = "variable",
                    times = c("pr", "tasmin", "tasmax", "pet"))

obs_monthly <- aggregate(value ~ station + variable + month, data = obs_long, FUN = mean)
raw_monthly <- aggregate(value ~ station + modelname + variable + month, data = raw_df, FUN = mean)
bc_monthly <- aggregate(value ~ station + modelname + variable + month, data = bc_df, FUN = mean)

model_names <- sort(unique(raw_df$modelname))
palette_colors <- scales::hue_pal()(length(model_names))
names(palette_colors) <- model_names

# ---- Plot per Variable ----
variables <- c("pr", "tasmin", "tasmax", "pet")
label.vars <- c("Rainfall", "Tmin", "Tmax", "PET")
ylabel.vars <- c("pr (mm)", "tasmin (°C)", "tasmax (°C)", "pet (mm)")
names(label.vars) <- names(ylabel.vars) <- variables

for (st in unique(obs_monthly$station)) {
  #st = "OUAHIGOUYA"
  plot_list <- list()
  plot_index <- 1

  obs_var <- obs_monthly[obs_monthly$station == st, ]
  raw_var <- raw_monthly[raw_monthly$station == st, ]
  bc_var  <- bc_monthly[bc_monthly$station == st, ]
  
  message("Station: ", st)
  for (var in variables) {
   
    #var = "pr"
    message("....... ", var)
    
    # Filter data for variable
    obs_plot <- obs_var[obs_var$variable == var, ]
    raw_plot <- raw_var[raw_var$variable == var, ]
    bc_plot  <- bc_var[bc_var$variable == var, ]
    
    obs_plot$source <- "OBSERVED"
    raw_plot$source <- raw_plot$modelname
    bc_plot$source  <- bc_plot$modelname
    
    sorted_models <- sort(model_names)
    legend_order <- c(sorted_models, "OBSERVED")
    
    raw_combined <- rbind(
      raw_plot[, c("month", "value", "source")],
      obs_plot[, c("month", "value", "source")]
    )
    bc_combined <- rbind(
      bc_plot[, c("month", "value", "source")],
      obs_plot[, c("month", "value", "source")]
    )
    
    raw_combined$source <- factor(raw_combined$source, levels = legend_order)
    bc_combined$source  <- factor(bc_combined$source, levels = legend_order)
    
    sorted_models <- sort(model_names)
    full_sources <- c(sorted_models, "OBSERVED")
    full_colors <- c(palette_colors[sorted_models], OBSERVED = "black")
    full_linetypes <- c(rep("solid", length(sorted_models)), OBSERVED = "dotted")
    names(full_colors) <- names(full_linetypes) <- full_sources
    
    p1 <- ggplot(data = raw_combined) +
      geom_line(aes(x = month, y = value, group = source, color = source, linetype = source), linewidth = 0.5) +
      scale_color_manual(values = full_colors) +
      scale_linetype_manual(values = full_linetypes) +
      labs(title = paste0(label.vars[var], " (", st, ") - raw"),
           x = NULL, y = paste0(ylabel.vars[var]), color = "Model", linetype = "Model") +
      guides(color = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3)) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title = element_text(size = 12, color = "black"),
        axis.title.y = element_text(margin = margin(r = 8)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")
      )
    
    
    p2 <- ggplot(data = bc_combined) +
      geom_line(aes(x = month, y = value, group = source, color = source, linetype = source), linewidth = 0.5) +
      scale_color_manual(values = full_colors) +
      scale_linetype_manual(values = full_linetypes) +
      labs(title = paste0(label.vars[var], " (", st, ") - bias corrected"),
           x = NULL, y = paste0(ylabel.vars[var]), color = "Model", linetype = "Model") +
      guides(color = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3)) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title = element_text(size = 12, color = "black"),
        axis.title.y = element_text(margin = margin(r = 8)),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")
      )
    
    p2

    #combined_plot <- ggarrange(plotlist = list(p1, p2), ncol = 1, align = "v",common.legend = T, legend = "bottom", align = "hv")
    plot_list[[plot_index]] <- p1
    plot_list[[plot_index + 4]] <- p2
    plot_index <- plot_index + 1
  }
  
  final_plot <- ggarrange(plotlist = plot_list, ncol = 4, nrow = 2, 
                          common.legend = T, legend = "bottom", align = "hv")
  
  ggsave(
    filename = paste0("output/figures/2_bias_correction_eval_",st,".png"),
    plot = final_plot,
    width = 20, height = 8, dpi = 350, bg = "white", scale = 1
  ) 

}



# ---- Density Distribution Plots: Raw vs. Bias Corrected vs. Observed ----
library(ggplot2)
library(ggpubr)

# ---- Density Grid Plot per Station ----

for (st in unique(obs_df$station)) {
  #st = "OUAGADOUGOU"
  message("Density layout grid - Station: ", st)
  
  # Prepare OBSERVED: wide to long
  obs_station <- obs_df[obs_df$station == st, ]
  obs_long <- reshape(
    obs_station[, c("date", "pr", "tasmin", "tasmax", "pet")],
    direction = "long",
    varying = list(c("pr", "tasmin", "tasmax", "pet")),
    v.names = "value",
    timevar = "variable",
    times = c("pr", "tasmin", "tasmax", "pet")
  )
  obs_long$source <- "OBSERVED"
  obs_long$modelname <- "OBSERVED"
  
  # Prepare RAW and BC data
  raw_station <- raw_df[raw_df$station == st & raw_df$scenario == "historical", ]
  bc_station  <- bc_df[bc_df$station == st & bc_df$scenario == "historical", ]
  raw_station$source <- "RAW"
  bc_station$source  <- "BIAS_CORRECTED"
  
  # Combine all
  combined_data <- rbind(
    obs_long[, c("variable", "value", "source", "modelname")],
    raw_station[, c("variable", "value", "source", "modelname")],
    bc_station[, c("variable", "value", "source", "modelname")]
  )
  
  # Clean NAs
  combined_data <- combined_data[!is.na(combined_data$value), ]
  
  # Set variable labels
  combined_data$variable <- factor(
    combined_data$variable,
    levels = variables,
    labels = label.vars[variables]
  )
  
  # Custom colors
  color_map <- c("OBSERVED" = "black", "RAW" = "red", "BIAS_CORRECTED" = "blue")
  
  # Build plot rows per model
  model_panels <- list()
  
  for (model in model_names) {
    #model = model_names[1]
    message("....... Model: ", model)
    
    model_data <- combined_data[
      combined_data$modelname %in% c("OBSERVED", model), ]
    
    row_plots <- list()
    
    for (v in levels(combined_data$variable)) {
      #v = "Rainfall"
      plot_df <- model_data[model_data$variable == v, ]
      
      p <- ggplot(plot_df, aes(x = value, color = source)) +
        stat_ecdf(linewidth = 0.5) +
        scale_color_manual(values = color_map) +
        labs(
          title = paste0(model, " - ", v),
          x = "", y = "",
          color = "Source"
        ) +
        theme_bw(base_size = 12) +
        theme(
          plot.title = element_text(size = 12, face = "bold", color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.title.x = element_text(margin = margin(t = 8)),
          axis.title.y = element_text(margin = margin(r = 8)),
          legend.position = "none"
        )
      #p
      row_plots[[v]] <- p
    }
    
    # Align 4 variable plots side by side for this model
    model_panels[[model]] <- ggarrange(
      plotlist = row_plots,
      ncol = 4, align = "hv"
    )
  }
  
  # Stack models vertically
  final_grid <- ggarrange(
    plotlist = model_panels,
    ncol = 1,
    common.legend = T,
    legend = "bottom"
  )
  
  # ggsave(
  #   filename = paste0("output/figures/2_ecdf_bias_correction_eval_", st, ".png"),
  #   plot = final_grid, limitsize = F,
  #   width = 20, height = 4 * length(model_names), dpi = 350, bg = "white", scale = 1
  # )
  
  ggsave(
    filename = paste0("output/figures/2_ecdf_bias_correction_eval_", st, ".pdf"),
    plot = final_grid,
    device = cairo_pdf,limitsize = F,
    width = 20, height = 4 * length(model_names),
    dpi = 600, bg = "white"
  )
}

print("All done!")
