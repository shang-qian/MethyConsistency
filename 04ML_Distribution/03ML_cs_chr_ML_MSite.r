####Chromosome-based methylation level
workdir=/Shang/03Project/03Analysis
mkdir -p $workdir/03ML_dist/01All/Cattle/chr
cd $workdir/03ML_dist/01All/Cattle/chr

module load r/4.4.0
R

#data.table input files
library(data.table)

dt_covF=fread("/Shang/03Project/03Analysis/02Samples/Cattle/02Final/03Final.ML.joint237.txt")
dt_BTF_ML=fread("../01Cattle_Breed_tissue_average_ML_All.txt")

dt_filtered <- copy(dt_BTF_ML)
dt_filtered <- dt_filtered[!Breed %in% c("Cattle")]
dt_filtered <- dt_filtered[!Tissue %in% c("Milk somatic cells")]
dt_filtered <- dt_filtered[!Tissue %in% c("Macrophages peritoneal")]

valid_breeds <- dt_filtered[, .N, by = Breed][N >= 3, Breed]
dt_filtered <- dt_filtered[Breed %in% valid_breeds]

valid_tissues <- dt_filtered[, .N, by = Tissue][N >= 3, Tissue]
dt_filtered <- dt_filtered[Tissue %in% valid_tissues]

dt_filtered[, Combo_Count := .N, by = .(Breed, Tissue)]

dt_filtered <- dt_filtered[Combo_Count >= 3]

#dt_filtered[Tissue == "Cerebral Cortex", Tissue := "Cortex"]
#dt_filtered[Tissue == "Macrophages peritoneal", Tissue := "pMac"]
#dt_filtered[Breed == "Chinese red steppes", Breed := "Red steppes"]



library(ggplot2)
chr_map <- data.table(
  chr = c(
    paste0("NC_0373", 28:56, ".1")#, # Chr1-9,X
    #"NC_082638.1" #,                  # ChrY
  ),
  Chr_Name = c(
    paste0("Chr", 1:29) #, "ChrX", "ChrY"  #, "MT"
  )
)

###04 ML along chromsome: Muscle tissue- 5breeds
window_size <- 500000  # 500Kb 
window_step <- 100000  # 100Kb 
ml_cols <- grep("^ML_", names(dt_covF), value = TRUE)

message("Step 3: Generating 500Kb Sliding Windows (100Kb shifts)...")

chr_bounds <- dt_covF[, .(max_pos = max(CpG_pos, na.rm = TRUE)), by = chr]
dt_windows <- chr_bounds[, {
  starts <- seq(1, max_pos, by = window_step)
  ends <- starts + window_size - 1
  .(win_start = starts, win_end = ends)
}, by = chr]

dt_windows[, win_id := .I]
setkey(dt_windows, chr, win_start, win_end)

dt_cpg_loc <- dt_covF[, .(chr, start = CpG_pos, end = CpG_pos)]
dt_cpg_loc[, (ml_cols) := dt_covF[, ml_cols, with = FALSE]]
setkey(dt_cpg_loc, chr, start, end)

message("--> Mapping CpG sites into overlapping windows...")
#dt_joined <- foverlaps(dt_cpg_loc, dt_windows, type = "within", nomatch = NULL)
dt_window <- dt_cpg_loc[dt_windows, 
  {
    lapply(.SD, function(col) mean(col, na.rm = TRUE))
  }, 
  on = .(chr, start <= win_end, end >= win_start), 
  by = .EACHI, 
  .SDcols = ml_cols]

dt_window[, chr := dt_windows$chr]
dt_window[, win_start := dt_windows$win_start]
dt_window[, win_end := dt_windows$win_end]

dt_window[, c("start", "end") := NULL]

setcolorder(dt_window, c("chr", "win_start", "win_end", ml_cols))

dt_win_long <- melt(dt_window, id.vars = c("chr", "win_start", "win_end"), 
                   variable.name = "ColName", value.name = "Win_ML")
dt_win_long <- dt_win_long[!is.na(Win_ML)]

dt_win_final <- merge(dt_win_long, dt_filtered[, .(ColName, Breed, Tissue)], by = "ColName")
dt_win_final <- merge(dt_win_final, chr_map, by = "chr", all.x = TRUE)
dt_win_final <- dt_win_final[!is.na(Chr_Name)]
dt_win_final[, Chr_Name := factor(Chr_Name, levels = c(paste0("Chr", 1:29), "ChrX", "ChrY"))]

dt_win_final[, Pos_Mb := (win_start + win_end) / 2 / 1e6]

dt_win_trend <- dt_win_final[, .(Group_ML = mean(Win_ML)), by = .(Chr_Name, Pos_Mb, Breed, Tissue)]
dt_win_trend <- dt_win_trend[Tissue=="Muscle"]

dt_win_trend[, Breed := factor(Breed, levels = c( "Holstein","Brangus","Wagyu","Brangus-Angus","Red steppes"))]
academic_palette <- c(
  "Holstein"   = "#377EB8",  
  "Wagyu"      = "#984EA3",  
  "Brangus"    = "#D95F02",  
  "Brangus-Angus"  = "#D9A441",    
  "Red steppes" = "#A65628"  
)


p_cen_vline <- ggplot() +
  geom_line(data = dt_win_trend, aes(x = Pos_Mb, y = Group_ML, color = Breed), 
            linewidth = 0.6, alpha = 0.85) +
  
  facet_wrap(~Chr_Name, ncol = 4, scales = "free_x") +
  #scale_color_brewer(palette = "Paired") +
  scale_color_manual(values = academic_palette) + 
  #scale_fill_manual(values = academic_palette1) +
  
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0.01, 0.01)) +
  
  theme_bw(base_size = 20) +
  labs(
    title = "Cattle DNA Methylation Landscapes Across Breeds",
    subtitle = "Genome-wide sliding-window analysis (500Kb windows shifted by 100Kb steps)",
    x = "Genomic position (Mb)", 
    y = "Mean methylation level",
    color = "Breed"
  ) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 13, face = "italic", hjust = 0.5, color = "gray30", margin = margin(b = 15)),

    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),

    strip.background = element_rect(fill = "gray95", color = "gray80"),
    strip.text = element_text(size = 14, face = "bold"),

    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.4),

    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 20)
  ) +

  guides(color = guide_legend(nrow = 1, byrow = TRUE, 
                              override.aes = list(linewidth = 3, alpha = 1)))

ggsave("04_Cattle_500Kb_Sliding_ML_Muscle_Breeds.pdf", p_cen_vline, 
       width = 22, height = 28, dpi = 300, limitsize = FALSE)



###05. 500Kb Texel Tissues
message("--> Aggregating group trends by Tissue...")
dt_win_trend <- dt_win_final[, .(Group_ML = mean(Win_ML)), by = .(Chr_Name, Pos_Mb, Breed, Tissue)]
dt_win_trend <- dt_win_trend[Breed=="Holstein"]

dt_win_trend[, Tissue := factor(Tissue, levels = c("Blood", "Sperm","Mammary","Cortex","Breast","Lung","Muscle","Rumen","Liver","Jejunum","ESCs"))]
academic_palette1 <- c(
 "Muscle"   = "#E41A1C",  
  "Blood"    = "#377EB8",  
  "Sperm"    = "#4DAF4A",  
  "Liver"    = "#984EA3",  
  
  "Mammary"  = "#D95F02",  
  "Lung"     = "#D9A441",  
  "Rumen"    = "#1B9E77",  
  "Cortex"   = "#7570B3",  
  "Breast"   = "#E7298A",  
  
  "Jejunum"  = "#66A61E",  
  "ESCs"     = "#A6761D"  
)

p_cen_vline <- ggplot() +
  geom_line(data = dt_win_trend, aes(x = Pos_Mb, y = Group_ML, color = Tissue), 
            linewidth = 0.6, alpha = 0.85) +
  
  facet_wrap(~Chr_Name, ncol = 4, scales = "free_x") +
  scale_color_manual(values = academic_palette1) + 
  #scale_fill_manual(values = academic_palette1) +
  
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0.01, 0.01)) +
  
  theme_bw(base_size = 20) +
  labs(
    title = "Cattle DNA Methylation Landscapes Across Tissues",
    subtitle = "Genome-wide sliding-window analysis (500Kb windows shifted by 100Kb steps)",
    x = "Genomic position (Mb)", 
    y = "Mean methylation level",
    color = "Tissue"
  ) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 13, face = "italic", hjust = 0.5, color = "gray30", margin = margin(b = 15)),

    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),

    strip.background = element_rect(fill = "gray95", color = "gray80"),
    strip.text = element_text(size = 14, face = "bold"),

    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.4),

    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 18),
    legend.spacing.x = unit(0.2, "cm"),     
    legend.key.width = unit(1.2, "cm")    
  ) +

  guides(color = guide_legend(nrow = 2, byrow = TRUE, 
                              override.aes = list(linewidth = 3, alpha = 1)))

ggsave("05_Cattle_500Kb_Sliding_ML_Holstein_Tissues.pdf", p_cen_vline, 
       width = 22, height = 28, dpi = 300, limitsize = FALSE)


################### 06. Density: Muscle ################################################
window_size <- 500000  # 500Kb 
window_step <- 100000  # 100Kb 

message("--> Computing rolling mean across all samples...")

dt_window <- dt_cpg_loc[dt_windows, 
  {
    lapply(.SD, function(col) sum(!is.na(col)))
  }, 
  on = .(chr, start <= win_end, end >= win_start), 
  by = .EACHI, 
  .SDcols = ml_cols]

dt_window[, chr := dt_windows$chr]
dt_window[, win_start := dt_windows$win_start]
dt_window[, win_end := dt_windows$win_end]

dt_window[, c("start", "end") := NULL]

setcolorder(dt_window, c("chr", "win_start", "win_end", ml_cols))

dt_win_long <- melt(dt_window, id.vars = c("chr", "win_start", "win_end"), 
                    variable.name = "ColName", value.name = "Win_ML")
dt_win_long <- dt_win_long[!is.na(Win_ML)]

dt_win_final <- merge(dt_win_long, dt_filtered[, .(ColName, Breed, Tissue)], by = "ColName")
dt_win_final <- merge(dt_win_final, chr_map, by = "chr", all.x = TRUE)
dt_win_final <- dt_win_final[!is.na(Chr_Name)]
dt_win_final[, Chr_Name := factor(Chr_Name, levels = c(paste0("Chr", 1:29), "ChrX", "ChrY"))]

dt_win_final[, Pos_Mb := (win_start + win_end) / 2 / 1e6]

dt_win_trend <- dt_win_final[, .(Group_ML = mean(Win_ML)), by = .(Chr_Name, Pos_Mb, Breed, Tissue)]

dt_win_trend <- dt_win_trend[Tissue=="Muscle"]

dt_win_trend[, Breed := factor(Breed, levels = c("Holstein","Brangus","Wagyu","Brangus-Angus","Red steppes"))]

message("Step 4: Generating Landscape with Vertical Lines...")

p_cen_vline <- ggplot() +
  geom_line(data = dt_win_trend, aes(x = Pos_Mb, y = Group_ML, color = Breed), 
            linewidth = 0.6, alpha = 0.85) +
  
  facet_wrap(~Chr_Name, ncol = 4, scales = "free_x") +
  scale_color_manual(values = academic_palette) + 
  
  scale_y_continuous(breaks = seq(0, 20000, 5000), limits = c(0, 20000), expand = c(0.01, 0.01)) +
  
  theme_bw(base_size = 20) +
  labs(
    title = "Cattle DNA Methylation sites Landscapes Across Breeds",
    subtitle = "Genome-wide sliding-window analysis (500Kb windows shifted by 100Kb steps)",
    x = "Genomic position (Mb)", 
    y = "Mean methylation sites number",
    color = "Breed"
  ) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 13, face = "italic", hjust = 0.5, color = "gray30", margin = margin(b = 15)),

    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),

    strip.background = element_rect(fill = "gray95", color = "gray80"),
    strip.text = element_text(size = 14, face = "bold"),

    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.4),

    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 20)
  ) +

  guides(color = guide_legend(nrow = 1, byrow = TRUE, 
                              override.aes = list(linewidth = 3, alpha = 1)))

ggsave("06_Cattle_500Kb_Sliding_MSite_Muscle_Breeds.pdf", p_cen_vline, 
       width = 22, height = 28, dpi = 300, limitsize = FALSE)

###05. 500Kb Holstein Tissues Msite
message("--> Aggregating group trends by Tissue...")
dt_win_trend <- dt_win_final[, .(Group_ML = mean(Win_ML)), by = .(Chr_Name, Pos_Mb, Breed, Tissue)]
dt_win_trend <- dt_win_trend[Breed=="Holstein"]

dt_win_trend[, Tissue := factor(Tissue, levels = c("Blood", "Sperm","Mammary","Cortex","Breast","Lung","Muscle","Rumen","Liver","Jejunum","ESCs","pMac"))]

p_cen_vline <- ggplot() +

  geom_line(data = dt_win_trend, aes(x = Pos_Mb, y = Group_ML, color = Tissue), 
            linewidth = 0.6, alpha = 0.85) +
  
  facet_wrap(~Chr_Name, ncol = 4, scales = "free_x") +
  scale_color_manual(values = academic_palette1) + 
  
  scale_y_continuous(breaks = seq(0, 20000, 5000), limits = c(0, 20000), expand = c(0.01, 0.01)) +
  
  theme_bw(base_size = 20) +
  labs(
    title = "Cattle DNA Methylation site Landscapes Across Tissues",
    subtitle = "Genome-wide sliding-window analysis (500Kb windows shifted by 100Kb steps)",
    x = "Genomic position (Mb)", 
    y = "Mean methylation site number",
    color = "Tissue"
  ) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 13, face = "italic", hjust = 0.5, color = "gray30", margin = margin(b = 15)),

    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),

    strip.background = element_rect(fill = "gray95", color = "gray80"),
    strip.text = element_text(size = 14, face = "bold"),

    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.4),

    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 18),
    legend.spacing.x = unit(0.2, "cm"),      
    legend.key.width = unit(1.2, "cm")      
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, 
                              override.aes = list(linewidth = 3, alpha = 1)))

ggsave("07_Cattle_500Kb_Sliding_MSite_Holstein_Tissues.pdf", p_cen_vline, 
       width = 22, height = 28, dpi = 300, limitsize = FALSE)

