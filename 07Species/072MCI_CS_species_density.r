
R
library(data.table)

dt_CS=fread("04Cattle_Sheep_matched_MS_95.bed")
dt_CS=dt_CS[,c(4:12,16:24)]
target_names=c("Cattle_chr","Cattle_pos","C_Ci","C_k","C_u","C_SampleN","C_Prop","C_Prela","C_state","Sheep_chr","Sheep_pos","S_Ci","S_k","S_u","S_SampleN","S_Prop","S_Prela","S_state")
setnames(dt_CS, target_names)

###density of prop
library(ggplot2)
dt_density <- dt_CS[!is.na(as.numeric(C_Prop)) & !is.na(as.numeric(S_Prop))]
dt_density[, `:=`(C_Prop = as.numeric(C_Prop), S_Prop = as.numeric(S_Prop))]
median_cattle <- median(dt_density$C_Prop)
median_sheep  <- median(dt_density$S_Prop)

plot_density_data <- melt(
  dt_density, 
  measure.vars = c("C_Prop", "S_Prop"),
  variable.name = "Species", 
  value.name = "Prop_Value"
)
plot_density_data[, Species := fifelse(Species == "C_Prop", "Cattle", "Sheep")]

nc_blue   <- "#377EB8"  
nc_orange <- "#EA580C" 

p_fig8a_clean <- ggplot(plot_density_data, aes(x = Prop_Value, color = Species)) +
  
  geom_density(linewidth = 1.2, adjust = 3.5, bw = "nrd0", key_glyph = "path") +
  
  scale_color_manual(values = c("Cattle" = nc_blue, "Sheep" = nc_orange)) +

  scale_x_continuous(limits = c(0.95, 1.0), breaks = seq(0.95, 1.0, by = 0.01), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) + 

  annotate("text", x = 0.9515, y = Inf, 
           label = paste0("Cattle Median: ", round(median_cattle, 3), "\nSheep Median: ", round(median_sheep, 3)),
           family = "sans", fontface = "bold.italic", size = 4.2, color = "#0F172A", hjust = 0, vjust = 1.4) +
  
  labs(
    x = "Sample prevalence (Prop)",
    y = "Density"
  ) +
  
  theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.title.x = element_text(face = "bold", size = 11, color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 11, color = "black", margin = margin(r = 10)),
    axis.text.x  = element_text(color = "black", size = 9.5), 
    axis.text.y  = element_text(color = "black", size = 9.5),
    axis.line    = element_line(linewidth = 0.6, color = "black"),
    axis.ticks   = element_line(linewidth = 0.6, color = "black"),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold", color = "black"),
    
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(), 
    legend.box.margin = margin(t = 5, b = 0),
    plot.margin = margin(t = 10, r = 15, b = 5, l = 10)
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 1.5,       
        #fill = "transparent",  
        alpha = 1              
      )
    )#,
    #fill = "none" 
  )

ggsave("Figure8A_Prop_Density_Clean_nofill.pdf", 
       plot = p_fig8a_clean, 
       width = 6.2, 
       height = 4.6, 
       dpi = 600,            
       useDingbats = FALSE)


