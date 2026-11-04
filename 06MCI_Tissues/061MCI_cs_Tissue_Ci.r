
workdir=/Shang/03Project/03Analysis
mkdir -p $workdir/05MCI_Tissue/Cattle/01Specific_Ci_u/01Tissue
cd $workdir/05MCI_Tissue/Cattle/01Specific_Ci_u/01Tissue

###Tissue
module load r/4.4.0
R

library(data.table)
library(matrixStats)
library(ggplot2)
library(fdrtool)


dt_Ci <- fread("/Shang/03Project/03Analysis/04MCI/Cattle/01Global/01_All_Cattle_Consistency_Final_List_ML_MS#.txt")

####Tissue-specific
ci_cols <- grep("^T_.*_Ci$", names(dt_Ci), value = TRUE)
dt_Ci[, Ci_count := rowSums(!is.na(.SD)), .SDcols = ci_cols]
dt_Ci[Ci_count == 1, Group := "G1"]
dt_Ci[Ci_count == 2, Group := "G2"]
dt_Ci[Ci_count >= 3, Group := "G3"]

dt_group1 <- dt_Ci[Ci_count == 1]
dt_group2 <- dt_Ci[Ci_count == 2]
dt_group3 <- dt_Ci[Ci_count >= 3]
rm(dt_Ci)
gc() 
###01.Group1
ci_cols      <- grep("^T_.*_Ci$", names(dt_group1), value = TRUE)
prop_cols    <- gsub("_Ci$", "_Prop", ci_cols)
msample_cols <- gsub("_Ci$", "_MSampleC", ci_cols)
u_cols       <- gsub("_Ci$", "_u", ci_cols)
k_cols       <- gsub("_Ci$", "_k", ci_cols)

active_col_idx <- max.col(!is.na(dt_group1[, ci_cols, with = FALSE]), ties.method = "first")
tissues_names  <- gsub("^T_|_Ci$", "", ci_cols)
dt_group1[, Specific_Tissue := tissues_names[active_col_idx]]
row_seq <- seq_len(nrow(dt_group1))
mat_idx <- matrix(c(row_seq, active_col_idx), ncol = 2)
dt_group1[, Specific_Ci       := as.matrix(.SD)[mat_idx], .SDcols = ci_cols]
dt_group1[, Specific_Prop     := as.matrix(.SD)[mat_idx], .SDcols = prop_cols]
dt_group1[, Specific_MSampleC := as.matrix(.SD)[mat_idx], .SDcols = msample_cols]
dt_group1[, Specific_u        := as.matrix(.SD)[mat_idx], .SDcols = u_cols]
dt_group1[, Specific_k        := as.matrix(.SD)[mat_idx], .SDcols = k_cols]
dt_G1_TSpec=dt_group1[Specific_Ci>=0.9&Specific_Prop>=0.95]
fwrite(dt_G1_TSpec, "01Cattle_Tissue_specific_Group1.txt", sep = "\t", quote = FALSE, na = ".")


####2. Group2
ci_cols      <- grep("^T_.*_Ci$", names(dt_group2), value = TRUE)
prop_cols    <- gsub("_Ci$", "_Prop", ci_cols)
msample_cols <- gsub("_Ci$", "_MSampleC", ci_cols)
u_cols       <- gsub("_Ci$", "_u", ci_cols)
k_cols       <- gsub("_Ci$", "_k", ci_cols)
tissues_names  <- gsub("^T_|_Ci$", "", ci_cols)

non_na_mat <- !is.na(dt_group2[, ci_cols, with = FALSE])
row_idx <- rep(seq_len(nrow(dt_group2)), each = 2)

col_idx_all <- which(t(non_na_mat), arr.ind = TRUE)[, "row"] 

active_idx_1 <- col_idx_all[seq(1, length(col_idx_all), by = 2)]
active_idx_2 <- col_idx_all[seq(2, length(col_idx_all), by = 2)]

row_seq <- seq_len(nrow(dt_group2))
mat_idx_1 <- matrix(c(row_seq, active_idx_1), ncol = 2)
mat_idx_2 <- matrix(c(row_seq, active_idx_2), ncol = 2)

dt_group2[, Tissue_1 := tissues_names[active_idx_1]]
dt_group2[, Tissue_2 := tissues_names[active_idx_2]]

dt_group2[, Ci_1       := as.matrix(.SD)[mat_idx_1], .SDcols = ci_cols]
dt_group2[, Ci_2       := as.matrix(.SD)[mat_idx_2], .SDcols = ci_cols]
dt_group2[, Prop_1     := as.matrix(.SD)[mat_idx_1], .SDcols = prop_cols]
dt_group2[, Prop_2     := as.matrix(.SD)[mat_idx_2], .SDcols = prop_cols]
dt_group2[, MSampleC_1 := as.matrix(.SD)[mat_idx_1], .SDcols = msample_cols]
dt_group2[, MSampleC_2 := as.matrix(.SD)[mat_idx_2], .SDcols = msample_cols]
dt_group2[, u_1        := as.matrix(.SD)[mat_idx_1], .SDcols = u_cols]
dt_group2[, u_2        := as.matrix(.SD)[mat_idx_2], .SDcols = u_cols]
dt_group2[, k_1        := as.matrix(.SD)[mat_idx_1], .SDcols = k_cols]
dt_group2[, k_2        := as.matrix(.SD)[mat_idx_2], .SDcols = k_cols]

dt_G2_candidate1 <- dt_group2[Ci_1 >= 0.9 & Prop_1 >= 0.95 & Ci_2 >= 0.9 & Prop_2 >= 0.95 &abs(u_1-u_2)>=0.5 &abs(k_1-k_2)>=0.5]
dt_G2_candidate2 <- dt_group2[Ci_1 >= 0.9 & Prop_1 >= 0.95 & abs(u_1-u_2)>=0.5 & abs(k_1-k_2)>=0.5]
dt_G2_candidate3 <- dt_group2[Ci_2 >= 0.9 & Prop_2 >= 0.95 & abs(u_1-u_2)>=0.5 & abs(k_1-k_2)>=0.5]
dt_G2_TSpec <- rbind(dt_G2_candidate1, dt_G2_candidate2, dt_G2_candidate3)
tissues <- c("Tissue_1","Tissue_2")
ci_cols    <- c("Ci_1","Ci_2")
u_cols     <- c("u_1","u_2")
prop_cols  <- c("Prop_1","Prop_2")
k_cols     <- c("k_1","k_2")
MSample_cols <- c("MSampleC_1","MSampleC_2")

dt_G2_Long_All <- melt(
  dt_G2_TSpec,
  id.vars = c("chr", "CpG_pos", "Group"), 
  measure.vars = list(tissues, u_cols, ci_cols, prop_cols,k_cols,MSample_cols),
  value.name = c("Tissue", "TMean_u", "TCi","TProp","Tk","TMsample")
#  variable.name = "Tissue"
)
dt_G2_Long_All[, variable := NULL]
dt_G2_Long_Strict <- unique(dt_G2_Long_All[ TCi >= 0.9 &  TProp >= 0.95 ])

setorder(dt_G2_Long_Strict, chr, CpG_pos, Tissue)
dt_G2_Long_Strict[, Locus_Global_ID := rleid(chr, CpG_pos)]
dt_G2_Long_Strict[, Passed_Tissue_Count := .N, by = Locus_Global_ID]
ML_cols <- grep("^ML_.*$", names(dt_group2), value = TRUE)
B_cols  <- grep("^B_.*$", names(dt_group2), value = TRUE)
T_cols  <- grep("^T_.*$", names(dt_group2), value = TRUE)
extra_cols <- c("chr", "CpG_pos", c(ML_cols, B_cols, T_cols))
dt_G2_Long_Annotated <- dt_group2[, ..extra_cols][dt_G2_Long_Strict, on = .(chr, CpG_pos)]
setorder(dt_G2_Long_Annotated, Locus_Global_ID, Tissue)
dt_G2_Final <- dt_G2_Long_Annotated
setorder(dt_G2_Final, chr, CpG_pos, Tissue)
fwrite(dt_G2_Final, "01Cattle_Tissue_specific_Group2.txt", sep = "\t", quote = FALSE, na = ".")

###03 Group3
ci_cols <- grep("^T_.*_Ci$", colnames(dt_group3), value = TRUE)
max_ci_vector <- rowMaxs(as.matrix(dt_group3[, ..ci_cols]), na.rm = TRUE)
dt_G3_filtered <- dt_group3[max_ci_vector >= 0.9]

u_cols <- grep("^T_.*_u$", colnames(dt_G3_filtered), value = TRUE)
X <- as.matrix(dt_G3_filtered[, ..u_cols])
N_valid <- rowSums(!is.na(X))
l1 <- rowSums(X, na.rm = TRUE)
l2 <- sqrt(rowSums(X^2, na.rm = TRUE) + 1e-20)

# Spec_Hoyer
Spec_Hoyer <- (sqrt(N_valid) - (l1 / l2)) / (sqrt(N_valid) - 1)
Spec_Hoyer[is.nan(Spec_Hoyer) | is.infinite(Spec_Hoyer)] <- 0

# Spec_Tau
X_max <- rowMaxs(X, na.rm = TRUE)
X_ratio <- X / X_max
X_ratio[is.na(X_ratio)] <- 1 
Tau <- rowSums(1 - X_ratio, na.rm = TRUE) / (N_valid - 1)
Tau[is.nan(Tau) | is.infinite(Tau)] <- 0

# Spec_entropy
row_sums_X <- rowSums(X, na.rm = TRUE) + 1e-20
p <- X / row_sums_X 
Entropy <- -rowSums(p * log(p + 1e-20), na.rm = TRUE) / log(N_valid)
Spec_entropy <- 1 - Entropy
Spec_entropy[is.nan(Spec_entropy) | is.infinite(Spec_entropy)] <- 0

dt_G3_filtered[, `:=`(
  Spec_Hoyer   = Spec_Hoyer,
  Spec_Tau     = Tau,
  Spec_entropy = Spec_entropy
)]

vars <- sapply(dt_G3_filtered[, .(Spec_Hoyer, Spec_Tau, Spec_entropy)], var, na.rm = TRUE)
weights <- vars / sum(vars)

dt_G3_filtered[, Spec_W3 := weights[1] * Spec_Hoyer + weights[2] * Spec_Tau + weights[3] * Spec_entropy]
normalize01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
dt_G3_filtered[, Spec_W3Norm := normalize01(Spec_W3)]


# 3.2. abs u,k >=0.5 ############
u_cols <-grep("^T_.*_u$", colnames(dt_G3_filtered), value = TRUE)
X_all <- as.matrix(dt_G3_filtered[, ..u_cols])
dt_G3_filtered[, `:=`(
  Global_Max_u = rowMaxs(X_all, na.rm = TRUE),
  Global_Min_u = rowMins(X_all, na.rm = TRUE)
)]
dt_G3_filtered[, Delta_u := Global_Max_u - Global_Min_u]
dt_G3_sig <- dt_G3_filtered[Delta_u >= 0.5&Spec_W3Norm >= 0.7]


################################################
#4.List tissue row
################################################
tissues <- gsub("^T_|_u$", "", u_cols)
u_cols      <- paste0("T_", tissues, "_u")
ci_cols     <- paste0("T_", tissues, "_Ci")
prop_cols   <- paste0("T_", tissues, "_Prop")
k_cols   <- paste0("T_", tissues, "_k")
Msample_cols   <- paste0("T_", tissues, "_MSampleC")
dt_G3_Long_All <- melt(
  dt_G3_sig,
  id.vars = c("chr","CpG_pos","Group","Spec_Hoyer","Spec_Tau","Spec_entropy","Spec_W3Norm","Global_Max_u","Global_Min_u","Delta_u"), 
  measure.vars = list(u_cols, ci_cols, prop_cols, k_cols, Msample_cols),
  value.name = c("TMean_u", "TCi", "TProp","Tk","TMSample"),
  variable.name = "Tissue"
)
dt_G3_Long_All[, Tissue := tissues[as.integer(Tissue)]]
dt_G3_Long_Strict <- dt_G3_Long_All[ TCi >= 0.9 &  TProp >= 0.95 ]

setorder(dt_G3_Long_Strict, chr, CpG_pos, Tissue)
dt_G3_Long_Strict[, Locus_Global_ID := rleid(chr, CpG_pos)]
dt_G3_Long_Strict[, Passed_Tissue_Count := .N, by = Locus_Global_ID]

ML_cols <- grep("^ML_.*$", names(dt_G3_sig), value = TRUE)
B_cols  <- grep("^B_.*$", names(dt_G3_sig), value = TRUE)
T_cols  <- grep("^T_.*$", names(dt_G3_sig), value = TRUE)
extra_cols <- c("chr", "CpG_pos", c(ML_cols, B_cols, T_cols))

dt_G3_Long_Annotated <- dt_G3_sig[, ..extra_cols][dt_G3_Long_Strict, on = .(chr, CpG_pos)]
setorder(dt_G3_Long_Annotated, Locus_Global_ID, Tissue)
dt_G3_Long_Annotated[, `:=`(
  Locus_Max_u = max(TMean_u, na.rm = TRUE),
  Locus_Min_u = min(TMean_u, na.rm = TRUE)
), by = Locus_Global_ID]

dt_G3_Long_Annotated[, Locus_Delta_u := Locus_Max_u - Locus_Min_u]
dt_G3_Final <- dt_G3_Long_Annotated[(Passed_Tissue_Count>1&Locus_Delta_u>=0.5)|(Passed_Tissue_Count==1&Spec_W3Norm>=0.8)]

setorder(dt_G3_Final, chr, CpG_pos, Tissue)
dt_G3_Final[, Locus_Global_ID := rleid(chr, CpG_pos)]
dt_G3_Final[, Passed_Tissue_Count := .N, by = Locus_Global_ID]

fwrite(dt_G3_Final, "01Cattle_Tissue_specific_Group3.txt", sep = "\t", quote = FALSE, na = ".")

