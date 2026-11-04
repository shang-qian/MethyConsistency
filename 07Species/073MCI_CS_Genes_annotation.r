
workdir=/Shang/03Project/03Analysis
mkdir -p $workdir/06Species/01Species
cd $workdir/06Species/01Species

R

library(rtracklayer)
library(GenomicFeatures)
library(ggplot2)
library(data.table)
library(viridis)
library(GenomicRanges)
library(patchwork) 
library(scales)

#####Function gene annotation
TE_data <- import("/Shang/Reference/Cattle/ARS-UCD_v2/GCF_002263795.3_ARS-UCD2.0_repeats.gff3")
te_dt <- as.data.table(TE_data)
chr_map <- data.table(
  chr = c(paste0("NC_0373", 28:57, ".1"), "NC_082638.1" ),
  Chr_Name = c(paste0(1:29), "X", "Y" )
)
te_dt <- merge(te_dt, chr_map, by.x = "seqnames", by.y = "Chr_Name", all.x = TRUE)

target_classes <- c("LINE", "SINE", "LTR", "DNA")
##"+"
te_final_dt_pos <- te_dt[ strand == "+" & class %in% target_classes & !is.na(chr) ] # 确保只保留 1-29, X, Y 染色体上的位点
##"-"
te_final_dt_neg <- te_dt[ strand == "-" & class %in% target_classes & !is.na(chr) ] # 确保只保留 1-29, X, Y 染色体上的位点

te_region_gr_pos <- GRanges(
  seqnames = te_final_dt_pos$chr,
  ranges = IRanges(start = te_final_dt_pos$start, end = te_final_dt_pos$end),
  strand = te_final_dt_pos$strand,
  TE_Class = te_final_dt_pos$class,
  TE_Family = te_final_dt_pos$family
)
te_region_gr_neg <- GRanges(
  seqnames = te_final_dt_neg$chr,
  ranges = IRanges(start = te_final_dt_neg$start, end = te_final_dt_neg$end),
  strand = te_final_dt_neg$strand,
  TE_Class = te_final_dt_neg$class,
  TE_Family = te_final_dt_neg$family
)

print(te_final_dt_pos[, .N, by = class])
print(te_final_dt_neg[, .N, by = class])


gff_path <- "/Shang/Reference/Cattle/ARS-UCD_v2/GCF_002263795.3_ARS-UCD2.0_genomic.gff.gz"
txdb <- makeTxDbFromGFF(gff_path, format="gff")


promoters_gr <- promoters(txdb, upstream=2000, downstream=2000)
exons_gr <- exons(txdb)
introns_gr <- unlist(intronsByTranscript(txdb))
genes_gr <- genes(txdb)
intergenic_gr <- gaps(genes_gr) 
###Pos + strand
intergenic_gr_pos <- intergenic_gr[strand(intergenic_gr) == "+"]
promoters_gr_pos <- promoters_gr[strand(promoters_gr)=="+"]
exons_gr_pos <- exons_gr[strand(exons_gr)=="+"]
introns_gr_pos <- introns_gr[strand(introns_gr)=="+"]
###Neg - strand
intergenic_gr_neg <- intergenic_gr[strand(intergenic_gr) == "-"]
promoters_gr_neg <- promoters_gr[strand(promoters_gr)=="-"]
exons_gr_neg <- exons_gr[strand(exons_gr)=="-"]
introns_gr_neg <- introns_gr[strand(introns_gr)=="-"]

regions_list_pos <- list(
  "Promoter" = promoters_gr_pos,   "Exon" = exons_gr_pos,   "Intron" = introns_gr_pos,   "Intergenic" = intergenic_gr_pos,
  "TE-DNA" = te_region_gr_pos[te_region_gr_pos$TE_Class == "DNA"],
  "TE-SINE" = te_region_gr_pos[te_region_gr_pos$TE_Class == "SINE"],
  "TE-LINE" = te_region_gr_pos[te_region_gr_pos$TE_Class == "LINE"],
  "TE-LTR" = te_region_gr_pos[te_region_gr_pos$TE_Class == "LTR"])
regions_list_neg <- list(
  "Promoter" = promoters_gr_neg,   "Exon" = exons_gr_neg,   "Intron" = introns_gr_neg,   "Intergenic" = intergenic_gr_neg,
  "TE-DNA" = te_region_gr_neg[te_region_gr_neg$TE_Class == "DNA"],
  "TE-SINE" = te_region_gr_neg[te_region_gr_neg$TE_Class == "SINE"],
  "TE-LINE" = te_region_gr_neg[te_region_gr_neg$TE_Class == "LINE"],
  "TE-LTR" = te_region_gr_neg[te_region_gr_neg$TE_Class == "LTR"])

###input HighLowK-HighCi
dt_HH <- fread("/Shang/03Project/03Analysis/06Species/01Species/Sites/05_CattleSheep_MS_Site_Counts.txt") 
sites_gr <- GRanges(seqnames = dt_HH$Cattle_chr, 
                    ranges = IRanges(start = dt_HH$Cattle_pos, end = dt_HH$Cattle_pos),type=dt_HH$C_stat)

low_sites  <- sites_gr[sites_gr$type == "Low"]
high_sites <- sites_gr[sites_gr$type == "High"]

####calculate 
nc_max_pos <- c(
  "NC_037328.1" = 158534110, # Chr 1
  "NC_037329.1" = 136231102, # Chr 2
  "NC_037330.1" = 121005158, # Chr 3
  "NC_037331.1" = 120000601, # Chr 4
  "NC_037332.1" = 120089316, # Chr 5
  "NC_037333.1" = 117806340, # Chr 6
  "NC_037334.1" = 110682743, # Chr 7
  "NC_037335.1" = 113319770, # Chr 8
  "NC_037336.1" = 105454467, # Chr 9
  "NC_037337.1" = 103308737, # Chr 10
  "NC_037338.1" = 106982474, # Chr 11
  "NC_037339.1" = 87216183,  # Chr 12
  "NC_037340.1" = 83472345,  # Chr 13
  "NC_037341.1" = 82403003,  # Chr 14
  "NC_037342.1" = 85007780,  # Chr 15
  "NC_037343.1" = 81013979,  # Chr 16
  "NC_037344.1" = 73167244,  # Chr 17
  "NC_037345.1" = 65820629,  # Chr 18
  "NC_037346.1" = 63449741,  # Chr 19
  "NC_037347.1" = 71974595,  # Chr 20
  "NC_037348.1" = 69862954,  # Chr 21
  "NC_037349.1" = 60773035,  # Chr 22
  "NC_037350.1" = 52498615,  # Chr 23
  "NC_037351.1" = 62317253,  # Chr 24
  "NC_037352.1" = 42350435,  # Chr 25
  "NC_037353.1" = 51992305,  # Chr 26
  "NC_037354.1" = 45612108,  # Chr 27
  "NC_037355.1" = 45940150,  # Chr 28
  "NC_037356.1" = 51098607,  # Chr 29
  "NC_037357.1" = 139009144, # Chr X
  "NC_082638.1" = 59476289,  # Chr Y
  "NC_001941.1" = 16340      # Chr MT 
)
genome_full_pos <- GRanges(seqnames = names(nc_max_pos),
                            ranges = IRanges(start = 1, end = as.numeric(nc_max_pos)), strand = "+" )
genome_full_neg <- GRanges(seqnames = names(nc_max_pos),
                            ranges = IRanges(start = 1, end = as.numeric(nc_max_pos)), strand = "-" )

# Exon -> Intron -> Intergenic ---
# A. Exon 
final_exon_gr_pos <- reduce(regions_list_pos$Exon)
final_exon_gr_neg <- reduce(regions_list_neg$Exon)
# B. Intron 
all_genes_pos <- reduce(genes(txdb)[strand(genes(txdb)) == "+"])
final_intron_gr_pos <- setdiff(all_genes_pos, final_exon_gr_pos)
all_genes_neg <- reduce(genes(txdb)[strand(genes(txdb)) == "-"])
final_intron_gr_neg <- setdiff(all_genes_neg, final_exon_gr_neg)
# C. Intergenic
final_inter_gr_pos <- setdiff(genome_full_pos, all_genes_pos)
Feature_list_pos <- list(Exon = final_exon_gr_pos, Intron = final_intron_gr_pos, Intergenic = final_inter_gr_pos)
final_inter_gr_neg <- setdiff(genome_full_neg, all_genes_neg)
Feature_list_neg <- list(Exon = final_exon_gr_neg, Intron = final_intron_gr_neg, Intergenic = final_inter_gr_neg)

element_names <- c("Promoter", "TE-DNA", "TE-SINE", "TE-LINE", "TE-LTR")
element_list_pos <- lapply(element_names, function(n) {
  reduce(regions_list_pos[[n]][strand(regions_list_pos[[n]]) == "+"])
})
element_list_neg <- lapply(element_names, function(n) {
  reduce(regions_list_neg[[n]][strand(regions_list_neg[[n]]) == "-"])
})
names(element_list_pos) <- element_names
names(element_list_neg) <- element_names

calc_refined_stats <- function(gr_list, sites_bg, sites_hi, sites_lo) {
  rbindlist(lapply(names(gr_list), function(n) {
    gr <- gr_list[[n]]
    len_bp <- sum(as.numeric(width(gr)))
    c_bg <- sum(countOverlaps(sites_bg, gr) > 0)
    c_hi <- sum(countOverlaps(sites_hi, gr) > 0)
    c_lo <- sum(countOverlaps(sites_lo, gr) > 0)
    data.table(Feature = n, Total_MS = c_bg, High_ML = c_hi, Low_ML = c_lo,
               Length_Mb = len_bp, High_Density = c_hi/(len_bp/1e6), Low_Density = c_lo/(len_bp/1e6))
  }))
}

feature_df_pos <- calc_refined_stats(Feature_list_pos, sites_gr, high_sites, low_sites)
feature_df_neg <- calc_refined_stats(Feature_list_neg, sites_gr, high_sites, low_sites)
element_df_pos <- calc_refined_stats(element_list_pos, sites_gr, high_sites, low_sites)
element_df_neg <- calc_refined_stats(element_list_neg, sites_gr, high_sites, low_sites)
final_df_pos <- rbind(feature_df_pos, element_df_pos)
final_df_neg <- rbind(feature_df_neg, element_df_neg)
fwrite(final_df_pos, "01_All_cattle_HighLowML_HighCi_statistics_genomic_features_pos.txt", sep = "\t", quote = FALSE)
fwrite(final_df_neg, "01_All_cattle_HighLowML_HighCi_statistics_genomic_features_neg.txt", sep = "\t", quote = FALSE)


plot_density <- melt(final_df_pos, 
                     id.vars = "Feature", 
                     measure.vars = c("High_Density", "Low_Density"),
                     variable.name = "Group", 
                     value.name = "Density")
feature_order <- c("Promoter", "Exon", "Intron", "TE-DNA", "TE-SINE", "TE-LINE", "TE-LTR", "Intergenic")
plot_density[, Feature := factor(Feature, levels = feature_order)]

academic_colors <- c(
  "Promoter"   = "#D63031",   
  "Exon"       = "#FAB1A0",   
  "Intron"     = "#FDCB6E",   
  "TE-DNA"     = "#00A8FF",   
  "TE-SINE"    = "#9C88FF",   
  "TE-LINE"    = "#006266",   
  "TE-LTR"     = "#E1B12C",   
  "Intergenic" = "#7F8C8D"    
)

p_density_nc_compact <- ggplot(plot_density, aes(x = Feature, y = Density, fill = Feature)) +
  geom_bar(stat = "identity", width = 0.6, color = NA, alpha = 0.9) +

  facet_wrap(~Group, scales = "free_y", ncol = 2, 
             labeller = as_labeller(c("High_Density" = "High-ML Density", 
                                      "Low_Density" = "Low-ML Density"))) +

  geom_text(aes(label = sprintf("%.2f", Density)), 
            vjust = -0.5, 
            size = 3.8, 
            fontface = "bold", color = "grey20") +
  
  scale_fill_manual(values = academic_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  
  theme_classic(base_size = 14) +
  labs(
    title = "Genomic Enrichment Density of Methylation",
    subtitle = "Normalized by Sites per Mb in Positive strand",
    x = NULL,
    y = "Density (Sites / Mb)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 11, color = "grey30", margin = margin(b = 15)),

    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13, hjust = 0, color = "black"),

    axis.text.x = element_text(angle = 30, hjust = 1, color = "black", face = "bold", size = 11),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    axis.line = element_line(color = "black", linewidth = 0.5),

    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave("02_Figure2_MS_Density_Facetted_Pos_cattle.pdf", 
       plot = p_density_nc_compact, 
       width = 12, 
       height = 5, 
       device = cairo_pdf)

###Negative strand
plot_density <- melt(final_df_neg, 
                     id.vars = "Feature", 
                     measure.vars = c("High_Density", "Low_Density"),
                     variable.name = "Group", 
                     value.name = "Density")

feature_order <- c("Promoter", "Exon", "Intron", "TE-DNA", "TE-SINE", "TE-LINE", "TE-LTR", "Intergenic")
plot_density[, Feature := factor(Feature, levels = feature_order)]

p_density_nc_compact <- ggplot(plot_density, aes(x = Feature, y = Density, fill = Feature)) +
  geom_bar(stat = "identity", width = 0.6, color = NA, alpha = 0.9) +

  facet_wrap(~Group, scales = "free_y", ncol = 2, 
             labeller = as_labeller(c("High_Density" = "High-ML Density", 
                                      "Low_Density" = "Low-ML Density"))) +

  geom_text(aes(label = sprintf("%.2f", Density)), 
            vjust = -0.5, 
            size = 3.8, 
            fontface = "bold", color = "grey20") +

  scale_fill_manual(values = academic_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 

  theme_classic(base_size = 14) +
  labs(
    title = "Genomic Enrichment Density of Methylation",
    subtitle = "Normalized by Sites per Mb in Positive strand",
    x = NULL,
    y = "Density (Sites / Mb)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 11, color = "grey30", margin = margin(b = 15)),

    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13, hjust = 0, color = "black"),

    axis.text.x = element_text(angle = 30, hjust = 1, color = "black", face = "bold", size = 11),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    axis.line = element_line(color = "black", linewidth = 0.5),

    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave("02_Figure2_MS_Density_Facetted_Neg_cattle.pdf", 
       plot = p_density_nc_compact, 
       width = 12, 
       height = 5, 
       device = cairo_pdf)



###03. dt_HH gene annotation
dt_HH[, Genomic_Type_Pos := as.character(".")]
dt_HH[, Genomic_Type_Neg := as.character(".")]

# Exon
hits_exon_pos <- findOverlaps(sites_gr, final_exon_gr_pos)
hits_exon_neg <- findOverlaps(sites_gr, final_exon_gr_neg)
dt_HH[unique(queryHits(hits_exon_pos)), Genomic_Type_Pos := "Exon"]
dt_HH[unique(queryHits(hits_exon_neg)), Genomic_Type_Neg := "Exon"]
# Intron
hits_intron_pos <- findOverlaps(sites_gr, final_intron_gr_pos)
hits_intron_neg <- findOverlaps(sites_gr, final_intron_gr_neg)
dt_HH[unique(queryHits(hits_intron_pos)), Genomic_Type_Pos := "Intron"]
dt_HH[unique(queryHits(hits_intron_neg)), Genomic_Type_Neg := "Intron"]
# Intergenic
hits_inter_pos <- findOverlaps(sites_gr, final_inter_gr_pos)
hits_inter_neg <- findOverlaps(sites_gr, final_inter_gr_neg)
dt_HH[unique(queryHits(hits_inter_pos)), Genomic_Type_Pos := "Intergenic"]
dt_HH[unique(queryHits(hits_inter_neg)), Genomic_Type_Neg := "Intergenic"]

missing_count <- dt_HH[Genomic_Type_Neg == ".", .N]
total_count   <- nrow(dt_HH)

stats_summary_pos <- dt_HH[, .(
  Count = .N, 
  Mean_Beta = mean(C_Ci, na.rm = TRUE),
  Mean_u =mean(C_u,na.rm=TRUE)
), by = Genomic_Type_Pos]
stats_summary_neg <- dt_HH[, .(
  Count = .N, 
  Mean_Beta = mean(C_Ci, na.rm = TRUE),
  Mean_u =mean(C_u,na.rm=TRUE)
), by = Genomic_Type_Neg]

print(stats_summary_pos)
print(stats_summary_neg)

genes_gr <- genes(txdb)
mcols(genes_gr)$gene_name <- names(genes_gr)
genes_gr_pos <- genes_gr[strand(genes_gr)=="+"]
genes_gr_neg <- genes_gr[strand(genes_gr)=="-"]

hits_gene_pos <- findOverlaps(sites_gr, genes_gr_pos)
hits_gene_neg <- findOverlaps(sites_gr, genes_gr_neg)

hits_dt_pos <- data.table(
  site_idx = queryHits(hits_gene_pos),
  gene_name = mcols(genes_gr_pos)$gene_name[subjectHits(hits_gene_pos)])
hits_dt_neg <- data.table(
  site_idx = queryHits(hits_gene_neg),
  gene_name = mcols(genes_gr_neg)$gene_name[subjectHits(hits_gene_neg)])

hits_collapsed_pos <- hits_dt_pos[, .(
  Collapsed_Names = paste(unique(gene_name), collapse = ";")
), by = site_idx]
hits_collapsed_neg <- hits_dt_neg[, .(
  Collapsed_Names = paste(unique(gene_name), collapse = ";")
), by = site_idx]

dt_HH[, Gene_Name_Pos := as.character(NA)]
dt_HH[, Gene_Name_Neg := as.character(NA)]

dt_HH[hits_collapsed_pos$site_idx, Gene_Name_Pos := hits_collapsed_pos$Collapsed_Names]
dt_HH[hits_collapsed_neg$site_idx, Gene_Name_Neg := hits_collapsed_neg$Collapsed_Names]


message("Sites annotated with Gene Name in Pos strand: ", sum(!is.na(dt_HH$Gene_Name_Pos)))
message("Sites annotated with Gene Name in Neg strand: ", sum(!is.na(dt_HH$Gene_Name_Neg)))

fwrite(dt_HH, "03_Cattle_MS_Sites_with_GenomicType_GeneAnnotation_Pos_Neg_strand.txt", sep = "\t", quote = FALSE, na = "NA")


