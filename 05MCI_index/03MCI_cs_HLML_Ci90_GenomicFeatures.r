
workdir=/Shang/03Project/03Analysis
mkdir -p $workdir/04MCI/Cattle/01Global/06Function
cd $workdir/04MCI/Cattle/01Global/06Function

R

library(rtracklayer)
library(GenomicFeatures)
library(ggplot2)
library(data.table)
library(viridis)
library(GenomicRanges)
library(patchwork) 

#####02 Function gene annotation
TE_data <- import("/Shang/Reference/Cattle/ARS-UCD_v2/GCF_002263795.3_ARS-UCD2.0_repeats.gff3")
te_dt <- as.data.table(TE_data)

chr_map <- data.table(
  chr = c(paste0("NC_0373", 28:57, ".1"), "NC_082638.1" ),
  Chr_Name = c(paste0(1:29), "X", "Y" )
)
te_dt <- merge(te_dt, chr_map, by.x = "seqnames", by.y = "Chr_Name", all.x = TRUE)
target_classes <- c("LINE", "SINE", "LTR", "DNA")

te_final_dt <- te_dt[
  strand == "+" & 
  class %in% target_classes & 
  !is.na(chr) 
]

te_region_gr <- GRanges(
  seqnames = te_final_dt$chr,
  ranges = IRanges(start = te_final_dt$start, end = te_final_dt$end),
  strand = te_final_dt$strand,
  TE_Class = te_final_dt$class,
  TE_Family = te_final_dt$family
)

print(te_final_dt[, .N, by = class])


gff_path <- "/cluster/pixstor/mckaysd-lab/Reference/Cattle/ARS-UCD_v2/GCF_002263795.3_ARS-UCD2.0_genomic.gff.gz"
txdb <- makeTxDbFromGFF(gff_path, format="gff")

promoters_gr <- promoters(txdb, upstream=2000, downstream=2000)
exons_gr <- exons(txdb)
introns_gr <- unlist(intronsByTranscript(txdb))

genes_gr <- genes(txdb)
intergenic_gr <- gaps(genes_gr) 
intergenic_gr <- intergenic_gr[strand(intergenic_gr) == "+"]

regions_list <- list(
  "Promoter" = promoters_gr,
  "Exon" = exons_gr,
  "Intron" = introns_gr,
  "TE-DNA" = te_region_gr[te_region_gr$TE_Class == "DNA"],
  "TE-SINE" = te_region_gr[te_region_gr$TE_Class == "SINE"],
  "TE-LINE" = te_region_gr[te_region_gr$TE_Class == "LINE"],
  "TE-LTR" = te_region_gr[te_region_gr$TE_Class == "LTR"],
  "Intergenic" = intergenic_gr
)

###input HighLowK-HighCi
dt_HH <- fread("/Shang/03Project/03Analysis/04MCI/Cattle/01Global/03_All_Cattle_Consistency_High_u_Ci90_95_100.txt") 
sites_gr <- GRanges(seqnames = dt_HH$chr, 
                    ranges = IRanges(start = dt_HH$CpG_pos, end = dt_HH$CpG_pos),type=dt_HH$ML_stat)

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
genome_full_plus <- GRanges(seqnames = names(nc_max_pos),
                            ranges = IRanges(start = 1, end = as.numeric(nc_max_pos)),
                            strand = "+")

#Exon -> Intron -> Intergenic ---
# A. Exon
final_exon_gr <- reduce(regions_list$Exon[strand(regions_list$Exon) == "+"])

# B. Intron 
all_genes_plus <- reduce(genes(txdb)[strand(genes(txdb)) == "+"])
final_intron_gr <- setdiff(all_genes_plus, final_exon_gr)

# C. Intergenic 
final_inter_gr <- setdiff(genome_full_plus, all_genes_plus)
vessel_list <- list(Exon = final_exon_gr, Intron = final_intron_gr, Intergenic = final_inter_gr)

element_names <- c("Promoter", "TE-DNA", "TE-SINE", "TE-LINE", "TE-LTR")
element_list <- lapply(element_names, function(n) {
  reduce(regions_list[[n]][strand(regions_list[[n]]) == "+"])
})
names(element_list) <- element_names

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

vessel_df <- calc_refined_stats(vessel_list, sites_gr, high_sites, low_sites)
element_df <- calc_refined_stats(element_list, sites_gr, high_sites, low_sites)
final_df <- rbind(vessel_df, element_df)
fwrite(final_df, "01_All_cattle_HighLowML_HighCi_statistics_genomic_features.txt", sep = "\t", quote = FALSE)

