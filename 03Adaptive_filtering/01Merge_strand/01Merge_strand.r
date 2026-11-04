####Biscuit strand and merge R script.

args <- commandArgs(trailingOnly = TRUE)
#args=c("/cluster/pixstor/mckaysd-lab/WORKING/Shang/03Project/20251103MethRate/03Analysis/01Coverage/Cattle_Angus_PRJNA626458/SRR24955917_100G/03Bwameth/04Meth_sites/02SRR24955917_100G.CG.meth.bwameth.chr.txt",
#"/cluster/pixstor/mckaysd-lab/WORKING/Shang/03Project/20251103MethRate/03Analysis/01Coverage/Cattle_Angus_PRJNA626458/SRR24955917_100G/03Bwameth/05Merge_sites",
#"SRR24955917_100G_bwameth")

library(data.table)

####
inputfile=args[1]
dt <- fread(inputfile, header = FALSE)
setnames(dt, c("chr", "pos", "TRead", "MRead", "ML", "strand"))

# CpG + and -
dt[, CpG_pos := ifelse(strand == "-", pos - 1, pos)]

# 合并 + 按规则输出坐标
merged <- dt[, {
  if (.N == 2 && all(c("+","-") %in% strand)) {
    list(
      pos = pos[strand == "+"][1],
      TRead = sum(TRead),
      MRead = sum(MRead),
      n=.N
    )
  } else {
    # - Strand
    list(
      pos = pos[1],
      TRead = sum(TRead),
      MRead = sum(MRead),
      n=.N
    )
  }
}, by = .(chr, CpG_pos)]

# 重新计算 ML
merged[, ML := MRead / TRead]

# 排序
setorder(merged, chr, pos)

# 输出
fwrite(merged[, .(chr, pos, TRead, MRead, ML, n)], 
       paste0(args[2],"/",args[3],".merged.txt"), sep = "\t")



