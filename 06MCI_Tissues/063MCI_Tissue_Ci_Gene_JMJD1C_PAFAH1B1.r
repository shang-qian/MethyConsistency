workdir=/Shang/03Project/03Analysis
mkdir -p $workdir/05MCI_Tissue/Cattle/01Specific_Ci_u/03Tissue_example
cd $workdir/05MCI_Tissue/Cattle/01Specific_Ci_u/03Tissue_example

R

library(rtracklayer)
library(GenomicFeatures)
library(ggplot2)
library(data.table)
library(viridis)
library(GenomicRanges)

dt_C_Tissue=fread("/Shang/03Project/03Analysis/03ML_dist/01All/Cattle/01Cattle_Breed_tissue_average_ML_All.txt",header=T)
dt_filtered <- copy(dt_C_Tissue)
dt_filtered <- dt_filtered[!Breed %in% c("Cattle")]
dt_filtered <- dt_filtered[!Tissue %in% c("Milk somatic cells")]
valid_breeds <- dt_filtered[, .N, by = Breed][N >= 3, Breed]
dt_filtered <- dt_filtered[Breed %in% valid_breeds]

valid_tissues <- dt_filtered[, .N, by = Tissue][N >= 3, Tissue]
dt_filtered <- dt_filtered[Tissue %in% valid_tissues]
dt_filtered[, Combo_Count := .N, by = .(Breed, Tissue)]
dt_filtered <- dt_filtered[Combo_Count >= 3]
dt_C_Tissue <- copy(dt_filtered)
dt_S_Tissue <- dt_C_Tissue


dt_Sperm_cattle <- fread("/Shang/03Project/03Analysis/05MCI_Tissue/Cattle/01Specific_Ci_u/02Tissue_function/01_cattle_Sperm_1specific_MS_with_Gene_annotation.txt")

######01. Single tissue-specific gene: JMJD1C
dt_S_VGLL4=dt_Sperm_cattle[grepl("JMJD1C", Gene_Name)]
target_row <- dt_S_VGLL4[TProp == 1][which.max(TCi)]

target_chrom       <- target_row$chr
target_pos         <- target_row$CpG_pos
target_tci         <- round(target_row$TCi, 2)
target_prop        <- round(target_row$TProp, 2)
target_spec_tissue <- target_row$Tissue 

ml_cols <- grep("^ML_", names(target_row), value = TRUE)
target_row[, (ml_cols) := lapply(.SD, as.character), .SDcols = ml_cols]

plot_data <- melt(target_row, 
                  id.vars = c("chr", "CpG_pos"), 
                  measure.vars = ml_cols,
                  variable.name = "ColName", 
                  value.name = "ML_value")

plot_data <- plot_data[ML_value != "." & ML_value != "" & !is.na(ML_value)]
plot_data[, ML_value := as.numeric(ML_value)]

plot_data <- merge(plot_data, dt_S_Tissue[, .(ColName, Sample_Tissue = Tissue)], by = "ColName", all.x = TRUE)
plot_data <- plot_data[!is.na(Sample_Tissue)]
plot_data <- plot_data[plot_data$Sample_Tissue!="Macrophages peritoneal"]
plot_data <- plot_data[plot_data$Sample_Tissue!="ESCs"]
print(unique(plot_data$Sample_Tissue))
tissue_colors <- c(
  "Sperm"       = "#E41A1C",  
  "Spleen"      = "#377EB8",  
  "Liver"       = "#4DAF4A",  
  "Mammary"     = "#984EA3",  
  "Caruncle"    = "#16A085",  
  "Muscle"      = "#D9A441",  
  "Cerebral Cortex" = "#A65628",  
  "Rumen"       = "#D95F02",   
  "Blood"       = "#B03A2E", 
  "Jejunum"     = "#7D6608" 
)

ann_caption <- paste0(
  "Specificity Profile | Target Tissue: ", target_spec_tissue, 
  "  MCI: ", target_tci, 
  "  MP: ", target_prop)
meta_summary <- plot_data[, .(N = .N), by = .(Sample_Tissue)]
meta_summary[, Tissue_Ci   := as.numeric(target_row[1, paste0("T_", Sample_Tissue, "_Ci"), with = FALSE])]
meta_summary[, Tissue_Prop := as.numeric(target_row[1, paste0("T_", Sample_Tissue, "_Prop"), with = FALSE])]
meta_summary[, X_Label := paste0(
  Sample_Tissue, "\n",
  "(n=", N, ", ",
  "MCI:", round(Tissue_Ci, 2), "\n",
  "MP:", round(Tissue_Prop, 2), ")"
)]

plot_data <- merge(plot_data, meta_summary[, .(Sample_Tissue, X_Label)], by = "Sample_Tissue")
p_site_scatter <- ggplot(plot_data, aes(x = X_Label, y = ML_value, color = Sample_Tissue)) +
  geom_jitter(width = 0.22, height = 0, size = 2.8, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.35, color = "grey20", linewidth = 0.4, alpha = 0.6) +
  scale_color_manual(values = tissue_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    title = paste0("Methylation Profile for JMJD1C Representative Site"),
    subtitle = paste0("Genomic Coordinates: ", target_chrom, " at pos ", target_pos),
    x = "Tissues",
    y = "Methylation level across individuals",
    caption = ann_caption
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(family = "sans", face = "bold", size = 12.5, hjust = 0.5),
    plot.subtitle = element_text(family = "sans", face = "plain", size = 9.5, hjust = 0.5, color = "grey40", margin = margin(b = 12)),
    axis.title = element_text(family = "sans", face = "bold", size = 10, color = "black"),
    axis.text.x = element_text(family = "sans", face = "bold", size = 9, color = "black", angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(family = "sans", color = "black", size = 9),
    legend.position = "none",
    plot.margin = margin(t = 12, r = 12, b = 12, l = 12)
  )
ggsave("07FigureC_JMJD1C_Site_Validation_Scatter_cattle.pdf", 
       plot = p_site_scatter, 
       width = 6.5, 
       height = 5, 
       dpi = 300, 
       useDingbats = FALSE)



####02. Duo tissue-specific sheep Gene: 
dt_S_Tissue=fread("/Shang/03Project/03Analysis/03ML_dist/01All/Sheep/01Sheep_Breed_tissue_average_ML_All.txt",header=T)
dt_filtered <- copy(dt_S_Tissue)
dt_filtered <- dt_filtered[!Breed %in% c("Unknow", "unknow", "NA", NA)]
dt_filtered <- dt_filtered[!Tissue %in% c("Unknow", "unknow", "NA", "Brain", NA)]

valid_breeds <- dt_filtered[, .N, by = Breed][N >= 3, Breed]
dt_filtered <- dt_filtered[Breed %in% valid_breeds]

valid_tissues <- dt_filtered[, .N, by = Tissue][N >= 3, Tissue]
dt_filtered <- dt_filtered[Tissue %in% valid_tissues]
dt_filtered[, Combo_Count := .N, by = .(Breed, Tissue)]
dt_filtered <- dt_filtered[Combo_Count >= 3]
dt_S_Tissue <- copy(dt_filtered)

dt_Thalamus_Sperm_sheep <- fread("/Shang/03Project/03Analysis/05MCI_Tissue/Sheep/01Specific_Ci_u/02Tissue_function/03_Sheep_Sperm_vs_Thalamus_2specific_MS_with_Gene_annotation.txt")
dt_S_VAPA=dt_Thalamus_Sperm_sheep[grepl("PAFAH1B1", Gene_Name)]
target_row <- dt_S_VAPA[TProp == 1][which.min(TCi)]

target_chrom       <- target_row$chr
target_pos         <- target_row$CpG_pos
target_tci         <- round(target_row$TCi, 4)
target_prop        <- round(target_row$TProp, 4)
target_spec_tissue <- target_row$Tissue  # 该位点理论上特异的组织（例如 Sperm）

ml_cols <- grep("^ML_", names(target_row), value = TRUE)
target_row[, (ml_cols) := lapply(.SD, as.character), .SDcols = ml_cols]

plot_data <- melt(target_row, 
                  id.vars = c("chr", "CpG_pos"), 
                  measure.vars = ml_cols,
                  variable.name = "ColName", 
                  value.name = "ML_value")

plot_data <- plot_data[ML_value != "." & ML_value != "" & !is.na(ML_value)]
plot_data[, ML_value := as.numeric(ML_value)]
plot_data <- merge(plot_data, dt_S_Tissue[, .(ColName, Sample_Tissue = Tissue)], by = "ColName", all.x = TRUE)

plot_data <- plot_data[!is.na(Sample_Tissue)]
print(unique(plot_data$Sample_Tissue))

sheep_tissue_colors <- c(
  "Sperm"         = "#E41A1C",  
  "Spleen"         = "#377EB8",  
  "Liver"           = "#4DAF4A",  
  "Mammary"          = "#984EA3",  
  "Caruncle"    = "#16A085",  
  "Muscle"      = "#D9A441",  
  "Ovary"              = "#A65628",  
  "Thalamus"  = "#F781BF",  
  "Hypothalamus"          = "#D95F02",  
  "Skin"          = "#2C3E50"  
)


ann_caption <- paste0(
  "Specificity Profile | Target Tissue: ", target_spec_tissue, 
  "  MCI: ", target_tci, 
  "  MP: ", target_prop)

meta_summary <- plot_data[, .(N = .N), by = .(Sample_Tissue)]
meta_summary[, Tissue_Ci   := as.numeric(target_row[1, paste0("T_", Sample_Tissue, "_Ci"), with = FALSE])]
meta_summary[, Tissue_Prop := as.numeric(target_row[1, paste0("T_", Sample_Tissue, "_Prop"), with = FALSE])]
meta_summary[, X_Label := paste0(
  Sample_Tissue, "\n",
  "(n=", N, ", ",
  "MCI:", round(Tissue_Ci, 2), "\n",
  "MP:", round(Tissue_Prop, 2), ")"
)]

plot_data <- merge(plot_data, meta_summary[, .(Sample_Tissue, X_Label)], by = "Sample_Tissue")

p_site_scatter <- ggplot(plot_data, aes(x = X_Label, y = ML_value, color = Sample_Tissue)) +
  geom_jitter(width = 0.22, height = 0, size = 2.8, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.35, color = "grey20", linewidth = 0.4, alpha = 0.6) +
  scale_color_manual(values = sheep_tissue_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    title = paste0("Methylation Profile for PAFAH1B1 Representative Site"),
    subtitle = paste0("Genomic Coordinates: ", target_chrom, " at pos ", target_pos),
    x = "Tissues",
    y = "Methylation level across individuals",
    caption = ann_caption
  ) +

  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(family = "sans", face = "bold", size = 12.5, hjust = 0.5),
    plot.subtitle = element_text(family = "sans", face = "plain", size = 9.5, hjust = 0.5, color = "grey40", margin = margin(b = 12)),
    axis.title = element_text(family = "sans", face = "bold", size = 10, color = "black"),
  
    axis.text.x = element_text(family = "sans", face = "bold", size = 9, color = "black", angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(family = "sans", color = "black", size = 9),

    legend.position = "none",
    plot.margin = margin(t = 12, r = 12, b = 12, l = 12)
  )

ggsave("08FigureD_PAFAH1B1_Site_Validation_Scatter_sheep2.pdf", 
       plot = p_site_scatter, 
       width = 6.5, 
       height = 5, 
       dpi = 300, 
       useDingbats = FALSE)

