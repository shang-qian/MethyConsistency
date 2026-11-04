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
    paste0("NC_0373", 28:56, ".1")#, # Chr1-29,#X
    #"NC_082638.1" #,                  # ChrY
  ),
  Chr_Name = c(
    paste0("Chr", 1:29) #, "ChrX", "ChrY"  #, "MT"
  )
)

ml_cols <- grep("^ML_", names(dt_covF), value = TRUE)
dt_chr_summary <- dt_covF[, lapply(.SD, mean, na.rm = TRUE), by = chr, .SDcols = ml_cols]
dt_chr_long <- melt(dt_chr_summary, id.vars = "chr", 
                    variable.name = "ColName", value.name = "Chr_Mean_ML")

dt_plot_final <- merge(dt_chr_long, dt_filtered[, .(ColName, Breed, Tissue)], by = "ColName")
dt_plot_final <- merge(dt_plot_final, chr_map, by = "chr", all.x = TRUE)

dt_plot_final=dt_plot_final[!is.na(Chr_Name)]

chr_levels <- c(paste0("Chr", 1:29)) #, "ChrX", "ChrY")
dt_plot_final[, Chr_Name := factor(Chr_Name, levels = intersect(chr_levels, unique(Chr_Name)))]


target_breed_order <- c("Holstein","Hereford", "Angus","Wagyu", "Brangus","Brangus-Angus","Red steppes") 
target_tissue_order <- c("Sperm","ESCs", "Cortex", "Blood","Breast","Lung","Rumen", "Jejunum","Liver","Mammary","Spleen","Caruncle","Muscle") #"pMac"
dt_plot_final[, Breed := factor(Breed, levels = intersect(target_breed_order, unique(Breed)))]
dt_plot_final[, Tissue := factor(Tissue, levels = intersect(target_tissue_order, unique(Tissue)))]

dt_plot_final <- dt_plot_final[!is.na(Chr_Mean_ML)]
academic_palette <- c(
  "Hereford" = "#E41A1C", 
  "Holstein" = "#377EB8",  
  "Angus"    = "#4DAF4A",  
  "Wagyu"    = "#984EA3", 
  "Brangus"  = "#D95F02", 
  "Brangus-Angus" = "#D9A441",  
  "Red steppes" = "#A65628"  
)

####01. Global ML across Chromosomes
dt_global_comp <- dt_plot_final[!is.na(Chr_Mean_ML)]
p_global_chr <- ggplot(dt_global_comp, aes(x = Chr_Name, y = Chr_Mean_ML, fill = Chr_Name)) +
  geom_violin(trim = TRUE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, color = "black", outlier.shape = NA, alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  labs(
    title = "Methylation Level Comparison Across All Chromosomes",
    x = "Chromosome",
    y = "Mean methylation level"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),    
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  )
ggsave("01_Global_Chr_Comparison.pdf", p_global_chr, width = 15, height = 8)


academic_palette_14 <- c(
  "Muscle"   = "#E41A1C",  
  "Blood"    = "#377EB8", 
  "Sperm"    = "#4DAF4A", 
  "Liver"    = "#984EA3",  
  "Spleen"   = "#A65628",  
  "Mammary"  = "#D95F02", 
  "Lung"     = "#D9A441",  
  "Rumen"    = "#1B9E77",  
  "Cortex"   = "#7570B3",  
  "Breast"   = "#E7298A",  
  "Jejunum"  = "#66A61E",  
  "Caruncle" = "#E6AB02",  
  "ESCs"     = "#A6761D" 
)

###02.Tissue comparisons
dt_fingerprint <- dt_global_comp[, .(Group_Mean = mean(Chr_Mean_ML)), by = .(Chr_Name, Tissue)]
p_tissue <- ggplot(dt_fingerprint, aes(x = Chr_Name, y = Group_Mean, group = Tissue, color = Tissue)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 3) +
  
  scale_color_manual(values = academic_palette_14) + 
  
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  labs(
    title = "Chromosome Methylation Levels by Tissue",
    x = "Chromosome",
    y = "Mean methylation level"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),    
    legend.position = "right",
    legend.text = element_text(size = 20, face = "bold"),
    legend.title = element_text(size =20 ,face = "bold")
  )
ggsave("02_Tissue_Chr_Comparison_cattle.pdf", p_tissue, width = 15, height = 8)


###03.Breed comparisons
dt_fingerprint <- dt_global_comp[, .(Group_Mean = mean(Chr_Mean_ML)), by = .(Chr_Name, Breed)]
p_tissue <- ggplot(dt_fingerprint, aes(x = Chr_Name, y = Group_Mean, group = Breed, color = Breed)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 3) +
  #scale_color_brewer(palette = "Paired") +
  scale_color_manual(values = academic_palette) + 
  
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  labs(
    title = "Chromosome Methylation Levels by Breed",
    x = "Chromosome",
    y = "Mean methylation level"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),    
    legend.position = "right",
    legend.text = element_text(size = 20, face = "bold"),
    legend.title = element_text(size =20 ,face = "bold")
  )
ggsave("02_Breed_Chr_Comparison_cattle.pdf", p_tissue, width = 15, height = 8)


