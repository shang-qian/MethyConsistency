####Chromosome-based methylation level
workdir=/Shang/03Project/03Analysis
mkdir -p $workdir/03ML_dist/01All/Cattle
cd $workdir/03ML_dist/01All/Cattle

module load r/4.4.0
R

#data.table input files
library(data.table)
#237 samples ML
dt_covF=fread("/Shang/03Project/03Analysis/02Samples/Cattle/02Final/03Final.ML.joint237.txt")
###01.Cattle_sample_information
SSBTfile="/Shang/03Project/01Script/03Sample/10Breed_Tissue/01Species_sample_breed_tissue.txt"
SSBT <- fread(SSBTfile,header=T)
dt_BT=SSBT[Species=="Cattle"][, .(BioSample,Breed,Tissue)]
setnames(dt_BT,c("SampleID","Breed","Tissue"))
dt_BT[, ColName := paste0("ML_", SampleID)]
valid_cols <- intersect(dt_BT$ColName, names(dt_covF))
dt_BTF <- dt_BT[ColName %in% valid_cols]
ml_cols <- grep("^ML_", names(dt_covF), value = TRUE)
#na.rm = TRUE remove na
sample_means <- colMeans(dt_covF[, ..ml_cols], na.rm = TRUE)
dt_means <- data.table(
  ColName = names(sample_means),
  Mean_ML = as.numeric(sample_means)
)
dt_BTF_ML <- merge(dt_BTF, dt_means, by = "ColName", all.x = TRUE)
dt_BTF_ML <- dt_BTF_ML[order(Mean_ML)]
fwrite(dt_BTF_ML, paste0("01Cattle_Breed_tissue_average_ML_All.txt"),sep = "\t", na = "NA",quote=F)


####01.Mean+-SD Figure for average ML
library(ggplot2)
dt_filtered <- copy(dt_BTF_ML)
dt_filtered <- dt_filtered[!Breed %in% c("Cattle")]
dt_filtered <- dt_filtered[!Tissue %in% c("Milk somatic cells")]
# Breed count >= 3 
valid_breeds <- dt_filtered[, .N, by = Breed][N >= 3, Breed]
dt_filtered <- dt_filtered[Breed %in% valid_breeds]
# Tissue count >=3
valid_tissues <- dt_filtered[, .N, by = Tissue][N >= 3, Tissue]
dt_filtered <- dt_filtered[Tissue %in% valid_tissues]
# combination >= 3 
dt_filtered[, Combo_Count := .N, by = .(Breed, Tissue)]
#dt_filtered[order(Combo_Count)]
dt_filtered <- dt_filtered[Combo_Count >= 3]

#dt_filtered[Tissue == "Cerebral Cortex", Tissue := "Cortex"]
#dt_filtered[Tissue == "Macrophages peritoneal", Tissue := "pMac"]
#dt_filtered[Breed == "Chinese red steppes", Breed := "Red steppes"]

cat("Raw total sample#:", nrow(dt_BTF_ML), "\n")
cat("Retained sample#:", nrow(dt_filtered), "\n")
cat("Retained Breed#:", length(unique(dt_filtered$Breed)), "\n")
cat("Retained Tissue#:", length(unique(dt_filtered$Tissue)), "\n")

library(ggplot2)
target_breed_order <- c("Holstein","Hereford", "Angus","Wagyu", "Brangus","Brangus-Angus","Red steppes") 
dt_filtered[, Breed := factor(Breed, levels = target_breed_order)]
dt_filtered[, Tissue := reorder(Tissue, Mean_ML, FUN = median)]
target_tissue_order <- c("Sperm","ESCs", "Cortex", "Blood","Breast","Lung","Rumen", "Jejunum","pMac","Liver","Mammary","Spleen","Caruncle","Muscle")
dt_filtered[, Tissue := factor(Tissue, levels = target_tissue_order)]

academic_palette <- c(
  "Hereford" = "#E41A1C",  
  "Holstein" = "#377EB8",  
  "Angus"    = "#4DAF4A",  
  "Wagyu"    = "#984EA3",  
  "Brangus"  = "#D95F02",  
  "Brangus-Angus" = "#D9A441",  
  "Red steppes"   = "#A65628"  
)



###plot ML distribution across breeds and tissues
p_nc_ultimate <- ggplot(dt_filtered, aes(x = Breed, y = Mean_ML, color = Breed, fill = Breed)) +
  stat_summary(fun = mean, 
               geom = "crossbar", 
               width = 0.4, 
               color = "#000000", 
               linewidth = 0.6,
               show.legend = FALSE) +
  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult = 1), 
               geom = "errorbar", 
               width = 0.15, 
               color = "#000000", 
               linewidth = 0.5,
               show.legend = FALSE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), 
             size = 2.0, 
             alpha = 0.9, 
             show.legend = TRUE) + 
  
  facet_wrap(~ Tissue, nrow = 1, scales = "free_x") +

  scale_color_manual(values = academic_palette) + 
  scale_fill_manual(values = academic_palette) +
  
  scale_y_continuous(breaks = seq(0.3, 1.0, 0.1), 
                     expand = expansion(mult = c(0.04, 0.06))) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 26, color = "black"),
    axis.title = element_text(face = "bold", size = 22, color = "black"),
    
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 18, color = "black", face = "bold"),
    axis.text.y = element_text(size = 18, color = "black", face = "bold"),
    
    axis.line = element_line(color = "black", linewidth = 1.0),
    axis.ticks = element_line(color = "black", linewidth = 1.0),
    axis.ticks.length = unit(0.2, "cm"),
    
    strip.background = element_rect(fill = "gray93", color = "gray80", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 16, color = "black"),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray93", linetype = "dashed", linewidth = 0.5),
    
    panel.spacing = unit(0.05, "lines"),
    legend.position = "bottom",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16, face = "bold"),
    legend.background = element_blank(),
    legend.box = "horizontal",
    
    plot.margin = margin(t = 25, r = 25, b = 15, l = 25, unit = "pt")
  ) +
  
  labs(
    title = "Global Mean Methylation Level across Breeds and Tissues",
    x = " ",
    y = "Mean methylation level per sample",
    color = "Breed"
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, 
                              override.aes = list(size = 6, alpha = 1)),
         fill = "none")
ggsave("02Cattle_avgML_TissueFacet_Breed_All_sample.pdf", 
       p_nc_ultimate, width = 17, height = 9, dpi = 300, limitsize = FALSE)
fwrite(dt_filtered, paste0("02Cattle_Breed_tissue_average_ML_All_199_FigureSourceData.txt"),sep = "\t", na = "NA",quote=F)


###########Density distribution for individual samples
library(ggridges)
set.seed(62)
dt_sub <- dt_covF[sample(.N, 1000000)]
dt_long_all <- melt(dt_sub, measure.vars = ml_cols, variable.name = "ColName", value.name = "ML")
dt_long_all <- merge(dt_long_all[!is.na(ML)], dt_filtered, by = "ColName")

p_ridges_samples <- ggplot(dt_long_all, aes(x = ML, y = Breed, fill = Breed)) +
  geom_vline(xintercept = seq(0, 1, 0.25), 
             color = "gray88",          
             linetype = "solid",     
             linewidth = 0.35) +
  geom_density_ridges(aes(group = ColName), 
                      scale = 0.95, 
                      alpha = 0.35, 
                      color = "black",
                      linewidth = 0.35,
                      rel_min_height = 0.005,
                      from = 0,
                      to =1,
                      quantile_lines = TRUE, quantiles = 2, 
                      vline_width = 0.2, vline_color = "gray30", vline_linetype = "dashed") +
  
  facet_wrap(~Tissue, scales = "free_y", ncol =5) +
  
  #scale_fill_brewer(palette = "Paired") +
  scale_color_manual(values = academic_palette) + 
  scale_fill_manual(values = academic_palette) +
  
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22, color = "black", margin = margin(b = 15)),
    axis.title.x = element_text(face = "bold", size = 18, color = "black", margin = margin(t = 12)),
    axis.title.y = element_text(face = "bold", size = 18, color = "black", margin = margin(r = 12)),
    
    axis.text.x = element_text(size = 13, color = "black", face = "plain"),
    axis.text.y = element_text(size = 13, color = "black", face = "bold"),
    
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.15, "cm"),
    
    strip.background = element_rect(fill = "gray93", color = "gray80", linewidth = 0.6),
    strip.text = element_text(face = "bold", size = 15, color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray93", linetype = "solid", linewidth = 0.5),
    panel.spacing.x = unit(1.0, "lines"),
    panel.spacing.y = unit(1.2, "lines"),
    
    legend.position = "bottom",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.background = element_blank(),
    legend.box = "horizontal",
    
    plot.margin = margin(t = 20, r = 20, b = 15, l = 20, unit = "pt")
  ) +
  
  labs(
    title = "Genome-wide DNA Methylation Density Profiles across Samples",
    x = "Methylation level", 
    y = "Breeds",
    fill = "Breed"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, 
                             override.aes = list(alpha = 1)), color = "none")

ggsave("03Cattle_ML_Density_Individual_Samples_NC.pdf", 
       p_ridges_samples, width = 18, height = 12, dpi = 300, limitsize = FALSE)


