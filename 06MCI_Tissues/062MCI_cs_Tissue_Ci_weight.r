
workdir=/Shang/03Project/20251103MethRate/03Analysis
mkdir -p $workdir/05MCI_Tissue/Cattle/01Specific_Ci_u/01Tissue
cd $workdir/05MCI_Tissue/Cattle/01Specific_Ci_u/01Tissue

###Tissue
module load r/4.4.0
R

library(data.table)
library(matrixStats)
library(ggplot2)

dt_g1 <- fread("01Cattle_Tissue_specific_Group1.txt")
dt_g2 <- fread("01Cattle_Tissue_specific_Group2.txt")
dt_g3 <- fread("01Cattle_Tissue_specific_Group3.txt")

ml_cols <- grep("^ML_", names(dt_g1), value = TRUE)
b_cols  <- grep("^B_", names(dt_g1), value = TRUE)
t_cols  <- grep("^T_", names(dt_g1), value = TRUE)

core_cols <- c("chr", "CpG_pos")

# ---- Group 1  ----
sub_g1 <- dt_g1[, c(core_cols, ml_cols, b_cols, t_cols, 
                    "Group", "Specific_Tissue", "Specific_u", "Specific_k", 
                    "Specific_MSampleC", "Specific_Ci", "Specific_Prop"), with = FALSE]
setnames(sub_g1, 
         old = c("Specific_Tissue", "Specific_u", "Specific_k", "Specific_MSampleC", "Specific_Ci", "Specific_Prop"),
         new = c("Tissue", "TMean_u", "Tk", "TMsample", "TCi", "TProp"))
sub_g1[, `:=`(Spec_Hoyer = NA_real_, Spec_Tau = NA_real_, Spec_entropy = NA_real_, Spec_W3Norm = NA_real_)]

# ---- Group 2  ----
sub_g2 <- dt_g2[, c(core_cols, ml_cols, b_cols, t_cols, 
                    "Group", "Tissue", "TMean_u", "Tk", "TMsample", "TCi", "TProp"), with = FALSE]
sub_g2[, `:=`(Spec_Hoyer = NA_real_, Spec_Tau = NA_real_, Spec_entropy = NA_real_, Spec_W3Norm = NA_real_)]
# ---- Group 3  ----
sub_g3 <- dt_g3[, c(core_cols, ml_cols, b_cols, t_cols, 
                    "Group", "Tissue", "TMean_u", "Tk", "TMSample", "TCi", "TProp",
                    "Spec_Hoyer", "Spec_Tau", "Spec_entropy", "Spec_W3Norm"), with = FALSE]
setnames(sub_g3, old = "TMSample", new = "TMsample")

dt_G123_Combined <- rbind(sub_g1, sub_g2, sub_g3, use.names = TRUE)

final_col_order <- c(
  "chr", "CpG_pos", 
  ml_cols, 
  b_cols, 
  t_cols, 
  "Spec_Hoyer", "Spec_Tau", "Spec_entropy", "Spec_W3Norm", 
  "TMean_u", "Tk", "TMsample", "TCi", "TProp", "Tissue", "Group"
)

setcolorder(dt_G123_Combined, final_col_order)

print(table(dt_G123_Combined$Group))

setorder(dt_G123_Combined, chr, CpG_pos, Tissue)
dt_G123_Combined[, Specific_MSorder := rleid(chr, CpG_pos)]
dt_G123_Combined[, Specific_TissueOrder := .N, by = Specific_MSorder]
print(table(dt_G123_Combined$Specific_TissueOrder))

fwrite(dt_G123_Combined, "02Cattle_Tissue_specific_combined_raw.txt", sep = "\t", quote = FALSE, na = ".")
dt = copy(dt_G123_Combined)
setDT(dt)

dt[, Tissue_Role := "Background"]
dt[Specific_TissueOrder %in% 1:2, Tissue_Role := "True_Specific"]

dt[Specific_TissueOrder >= 3, Tissue_Role := {
  base_median = median(TMean_u)
  ifelse(abs(TMean_u - base_median) >= 0.4, "True_Specific", "Background")
}, by = Specific_MSorder]

dt_clean <- dt[Tissue_Role == "True_Specific"]
setorder(dt_clean, chr, CpG_pos, Tissue)
dt_clean[, Specific2_MSorder := rleid(chr, CpG_pos)]
dt_clean[, Specific2_Torder := .N, by = Specific2_MSorder]
print(table(dt_clean$Specific2_Torder))
fwrite(dt_clean, "03Cattle_Tissue_specific_combined_clean2tissues.txt", sep = "\t", quote = FALSE, na = ".")

