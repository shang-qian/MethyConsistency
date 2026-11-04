
workdir=/Shang/03Project/03Analysis
mkdir -p $workdir/04MCI/Cattle/01Global
cd $workdir/04MCI/Cattle/01Global

###All cattle
module load r/4.4.0
R

library(data.table)
library(matrixStats)
library(ggplot2)
library(viridis)
library(scales)
library(mgcv) 
library(ggridges)

final_result<- fread("/Shang/03Project/03Analysis/04MCI/Cattle/01Global/01_All_Cattle_Consistency_Final_List_ML_MS#.txt")
final_result <- final_result[!is.na(All_Ci)]
setDT(final_result)

###01Fig0 Validation SD vs MCI
ml_cols <- grep("^ML_", names(final_result), value = TRUE)
n_samples <- length(ml_cols)
final_result[, All_SD := sqrt(rowVars(as.matrix(.SD), na.rm = TRUE)), .SDcols = ml_cols]
setDT(final_result)


p5_nc_submit <- ggplot(
  final_result[!is.na(All_Ci) & !is.na(All_SD)], 
  aes(x = All_Ci, y = All_SD)
) +
  stat_binhex(bins = 100) +
  scale_fill_viridis_c(
    option = "plasma", 
    trans = "log10", 
    labels = label_comma(),
    name = "Loci Count" 
  ) +
  
  geom_smooth(
    method = "gam", 
    color = "white", 
    linewidth = 1.0, 
    se = FALSE
  ) +
  
  scale_x_continuous(
    limits = c(0.2, 1.0), 
    breaks = seq(0.2, 1.0, by = 0.1),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = c(0.0, 0.55), 
    breaks = seq(0.0, 0.5, by = 0.1),
    expand = c(0.01, 0.01)
  ) +
  
  labs(
    x = "Methylation Consistency Index (MCI)",
    y = "Standard Deviation across Samples (SD)"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(1.8, "mm"),
    
    axis.text.x = element_text(color = "black", size = 11, margin = margin(t = 5)),
    axis.text.y = element_text(color = "black", size = 11, margin = margin(r = 5)), 
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 12)), 
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 12)),
    
    panel.grid.major = element_line(color = "grey95", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(color = "black", size = 9),
    legend.position = "right",
    
    plot.margin = margin(t = 10, r = 10, b = 10, l = 16, unit = "mm")
  ) +
  
  guides(
    fill = guide_colourbar(
      frame.colour = "black", 
      frame.linewidth = 0.6,
      ticks.colour = "black",
      ticks.linewidth = 0.5,
      barheight = unit(5, "cm"),
      barwidth = unit(0.4, "cm")
    )
  )

ggsave(
  filename = "Figure0_NC_MCI_vs_SD_Validation.pdf", 
  plot = p5_nc_submit, 
  width = 6.8, 
  height = 5.2, 
  useDingbats = FALSE
)




###02 MP ~ MCI
final_result[, Prop_Group := cut(All_Prop, 
                                 breaks = c(-0.01, 0.20, 0.50, 0.80, 0.95, 1.01),
                                 labels = c("0-20%", "20-50%", "50-80%", "80-95%", "95-100%"))]
counts_df <- final_result[, .(N = .N), by = Prop_Group]

label_mapping <- structure(
  paste0(counts_df$Prop_Group, "\n(MS#: ", comma(counts_df$N),")"),
  names = as.character(counts_df$Prop_Group)
)

final_result[, Prop_Group_Labeled := factor(label_mapping[as.character(Prop_Group)], 
                                            levels = label_mapping[order(names(label_mapping))])]

p_nc_ridges_two_lines <- ggplot(final_result, aes(x = All_Ci, y = Prop_Group_Labeled, fill = after_stat(x))) +
  geom_density_ridges_gradient(
    scale = 1.1, 
    rel_min_height = 0.005, 
    linewidth = 0.35, 
    color = "white"
  ) +
  scale_fill_viridis_c(
    option = "magma", 
    name = "MCI Value",
    breaks = seq(0.2, 1.0, by = 0.2)
  ) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "grey40", alpha = 0.6, linewidth = 0.4) +
  
  labs(
    title = "Methylation Consistency Across Sample Presence Strata",
    x = "Methylation consistency index (MCI)",
    y = "MS prevalence"
  ) +
  
  scale_x_continuous(
    limits = c(0.2, 1.0), 
    breaks = seq(0.2, 1.0, by = 0.1), 
    expand = expansion(mult = c(0.05, 0.05)) 
  ) +
  
  theme_ridges(font_size = 12, grid = TRUE) + 
  theme(
    axis.text.x = element_text(color = "black", size = 10, hjust = 0.5),

    axis.text.y = element_text(color = "black", size = 9.5, vjust = 0, hjust = 1, lineheight = 0.85),     
    
    axis.title.x = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(t = 10)), 
    axis.title.y = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(r = 10)),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    panel.background = element_rect(fill = "grey98"),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),

    plot.margin = margin(t = 5, r = 5, b = 5, l = 22, unit = "mm")
  )

ggsave("02Fig1_NC_MCI_Ridges_TwoLines_cattle.pdf", p_nc_ridges_two_lines, width = 8.0, height = 5.5, useDingbats = FALSE)



#  0-20%  (Hyper-Specific MS) Ci> =0.9
supp_data_0_20 <- Final_out[All_Prop >= 0 & All_Prop <= 0.20 & All_Ci >=0.9]
setDT(supp_data_0_20)
float_cols <- c("All_Ci", "All_k", "All_u", "All_Prop")
supp_data_0_20[, (float_cols) := lapply(.SD, round, digits = 3), .SDcols = float_cols]
fwrite(
  supp_data_0_20, 
  file = "02Supplementary_Data_1_Hyper_Specific_MS_0_20_Ci90.csv", 
  row.names = FALSE, 
  sep = ","
)
#95-100% (Core MS) Ci>=0.9
supp_data_95_100 <- Final_out[All_Prop >= 0.95 & All_Prop <= 1.00&All_Ci>=0.9]
setDT(supp_data_95_100)
float_cols <- c("All_Ci", "All_k", "All_u", "All_Prop")
supp_data_95_100[, (float_cols) := lapply(.SD, round, digits = 3), .SDcols = float_cols]
fwrite(
  supp_data_95_100, 
  file = "02Supplementary_Data_2_Core_Housekeeping_MS_95_100_Ci90.csv", 
  row.names = FALSE, 
  sep = ","
)


####High Consistency Ci and Prop, 阈值：Prop >= 0.95, Ci >= 0.9, k的分布和u的分布
setDT(supp_data_95_100)
p_fig2_nc_clean_ticks <- ggplot(supp_data_95_100, aes(x = All_u)) +
  geom_histogram(
    bins = 50, 
    fill = "#E64B35CC", 
    color = "white", 
    linewidth = 0.12,
    boundary = 0,
    closed = "left"
  ) +

  scale_y_log10(
    limits = c(1, 1000000),
    breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
    labels = comma,
    expand = expansion(mult = c(0, 0.04))
  ) + 
  
  scale_x_continuous(
    limits = c(0.0, 1.0), 
    breaks = seq(0.0, 1.0, by = 0.1),
    expand = c(0.005, 0.005)
  ) +
  
  labs(
    #title = "Distribution of Mean Methylation (u) within Highly Consistent Sites",
    #subtitle = "Analysis of Core Housekeeping Layer (Sample Prevalence >= 95%)",
    x = "u (Mean methylation level)",
    y = expression("Number of MS (Ci " >= " 0.9 and Prop " %in% " [0.95, 1])")
  ) +
  
  annotation_logticks(
    sides = "l", 
    scaled = TRUE, 
    colour = "black", 
    linewidth = 0.4,
    short = unit(1.0, "mm"),  
    mid = unit(1.6, "mm"),    
    long = unit(2.2, "mm")  
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    
    axis.ticks.x = element_line(color = "black", linewidth = 0.6),
    axis.ticks.y = element_blank(), 
    axis.ticks.length.x = unit(1.8, "mm"),
    
    axis.text.x = element_text(color = "black", size = 11, margin = margin(t = 5)),
    axis.text.y = element_text(color = "black", size = 11, margin = margin(r = 6)), 
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 12)), 
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 12)),
    
    panel.grid.major.y = element_line(color = "grey94", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(), 
    panel.grid.minor.x = element_blank(),
    
    plot.margin = margin(t = 8, r = 8, b = 8, l = 16, unit = "mm")
  )

ggsave(
  "03Figure2_NC_u_Log_Histogram_CleanTicks.pdf", 
  p_fig2_nc_clean_ticks, 
  width = 6.8, 
  height = 4.6, 
  useDingbats = FALSE
)

#######output High and low ML
supp_data_95_100[, ML_stat := cut(All_u, breaks = c(-Inf, 0.5, Inf),  labels = c("Low", "High"))]
fwrite(supp_data_95_100, "03_All_Cattle_Consistency_High_u_Ci90_95_100.txt", sep = "\t", quote = FALSE)




