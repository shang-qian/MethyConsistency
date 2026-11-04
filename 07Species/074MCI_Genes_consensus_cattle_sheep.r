workdir=/Shang/03Project/03Analysis
mkdir -p $workdir/06Species/01Species
cd $workdir/06Species/01Species

R

library(rtracklayer)
library(GenomicFeatures)
library(ggplot2)
library(data.table)
library(ggplot2)
library(viridis)
library(GenomicRanges)


###input cattle and sheep shared annotation lists
dt_cattle_gene <- fread("03_Cattle_MS_Sites_with_GenomicType_GeneAnnotation_Pos_Neg_strand.txt") 
dt_sheep_gene <- fread("03_Sheep_MS_Sites_with_GenomicType_GeneAnnotation_Pos_Neg_strand.txt") 

anno_cols <- c("Genomic_Type_Pos", "Genomic_Type_Neg", "Gene_Name_Pos", "Gene_Name_Neg")
key_cols <- c("Cattle_chr", "Cattle_pos", "Sheep_chr", "Sheep_pos")

dt_c <- copy(dt_cattle_gene)
setnames(dt_c, anno_cols, paste0("C_", anno_cols))

dt_s <- dt_sheep_gene[, c(key_cols, anno_cols), with = FALSE]
setnames(dt_s, anno_cols, paste0("S_", anno_cols))
dt_merged <- merge(dt_c, dt_s, by = key_cols, all = TRUE, sort = FALSE)
dim(dt_merged)   
fwrite(dt_merged, "05-1Cattle_Sheep_MS_Sites_GeneAnnotation_merged.txt",
       sep = "\t", quote = FALSE, na = "NA")


gene_cols <- c("C_Gene_Name_Pos", "C_Gene_Name_Neg", "S_Gene_Name_Pos", "S_Gene_Name_Neg")
has_gene <- dt_merged[, Reduce(`|`, lapply(.SD, function(x) !is.na(x) & x != "")),
                      .SDcols = gene_cols]
dt_gene <- dt_merged[has_gene]
fwrite(dt_gene, "05-2Cattle_Sheep_MS_Sites_GeneAnnotation_withGene.txt",
       sep = "\t", quote = FALSE, na = "NA")



split_genes <- function(...) {
  x <- c(...)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(character(0))
  unique(trimws(unlist(strsplit(x, ";"))))
}

dt_gene[, Shared_Genes := mapply(function(cp, cn, sp, sn) {
  shared <- intersect(split_genes(cp, cn), split_genes(sp, sn))
  if (length(shared) == 0) NA_character_ else paste(sort(shared), collapse = ";")
}, C_Gene_Name_Pos, C_Gene_Name_Neg, S_Gene_Name_Pos, S_Gene_Name_Neg)]


fwrite(dt_gene, "05-3Cattle_Sheep_MS_Sites_GeneAnnotation_sharedGene.txt",
       sep = "\t", quote = FALSE, na = "NA")

dt_shared <- dt_gene[!is.na(Shared_Genes)]
fwrite(dt_shared, "05-4Cattle_Sheep_MS_Sites_GeneAnnotation_sharedGene20601.txt",
       sep = "\t", quote = FALSE, na = "NA")




######KEGG genes
# ------------------------------------------
dt_gene_High <- dt_shared
gene_counts_dt <- dt_gene_High[, .(Shared_Genes)]
gene_counts_flat <- gene_counts_dt[, .(Shared_Genes = unlist(strsplit(as.character(Shared_Genes), ";")))]
final_gene_stats <- gene_counts_flat[, .(.N), by = Shared_Genes][order(-N)]

setnames(final_gene_stats, c("GeneID", "Site_Count"))
print(head(final_gene_stats, 20))
fwrite(final_gene_stats, "06CattleSheep_Gene_HighML_Site_Counts.txt", sep = "\t", quote = FALSE)

###02.MS count/gene length in Sheep
gff_path <- "/Shang/Reference/Sheep/ARS-UI_Ramb_v3/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.gff.gz"
txdb <- makeTxDbFromGFF(gff_path, format="gff")

genes_gr <- genes(txdb)
gene_lens <- width(genes_gr) 

gene_len_dt <- data.table(
  GeneID = mcols(genes_gr)$gene_id, 
  Length_bp = gene_lens               
)

gene_len_dt <- gene_len_dt[, .(GeneLength = max(Length_bp, na.rm = TRUE)), by = GeneID]

# Normalization
setkey(final_gene_stats, GeneID)
setkey(gene_len_dt, GeneID)
final_gene_stats <- gene_len_dt[final_gene_stats, nomatch = NULL]

## sites count / (gene length / 1000)
final_gene_stats[, Normalized_Density := Site_Count / (GeneLength / 1000)]

final_gene_stats <- final_gene_stats[order(-Normalized_Density)]
final_gene_stats[, Normalized_Density := round(Normalized_Density, 3)]

print(head(final_gene_stats, 20))

fwrite(final_gene_stats, "07Sheep_Gene_Normalized_Site_Densities.txt", 
       sep = "\t", quote = FALSE, row.names = FALSE)


GeneList_filter=final_gene_stats[Normalized_Density>=0.2&Site_Count>=5]
fwrite(GeneList_filter, "08Sheep_Gene_Normalized_Site_Densities0.2_Count5.txt", 
       sep = "\t", quote = FALSE, row.names = FALSE)


