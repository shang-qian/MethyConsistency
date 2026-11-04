#!/bin/bash

###Project information
#### 01 Bismark
workdir=/Shang/03Project/03Analysis/02Samples
species=Sheep
for file in $(ls $workdir/${species}/*_*/*/01*/04*/02*.CG.meth.Bismark.chr.txt)
do 
ls $file 
sample02=$(basename $file .CG.meth.Bismark.chr.txt)
sample=$(echo ${sample02#02})
echo $sample
cd $workdir/${species}/*_*/$sample
mkdir -p 06Merge_sites
outdir=$workdir/${species}/*_*/$sample/06Merge_sites
cd $outdir 
scriptdir=/Shang/03Project/20251103MethRate/01Script
sbatch $scriptdir/03Sample/06Merge_strand/02Merge_strand.sub $file $outdir ${sample}.bismark
done 

###02Bwameth
workdir=/Shang/03Project/03Analysis/02Samples
species=Sheep
for file in $(ls $workdir/${species}/*_*/*/03*/04*/02*.CG.meth.bwameth.chr.txt) 
do 
ls $file 
sample02=$(basename $file .CG.meth.bwameth.chr.txt)
sample=$(echo ${sample02#02})
echo $sample
cd $workdir/${species}/*_*/$sample
mkdir -p 06Merge_sites
outdir=$workdir/${species}/*_*/$sample/06Merge_sites
cd $outdir 
scriptdir=/Shang/03Project/01Script
sbatch $scriptdir/03Sample/06Merge_strand/02Merge_strand.sub $file $outdir ${sample}.bwameth
done 
####03 Biscuit
workdir=/Shang/03Project/03Analysis/02Samples
species=Sheep
for file in $(ls $workdir/${species}/*_*/*/04*/04*/02*.CG.meth.biscuit.chr.txt)
do 
ls $file 
sample02=$(basename $file .CG.meth.biscuit.chr.txt)
sample=$(echo ${sample02#02})
echo $sample
cd $workdir/${species}/*_*/$sample
mkdir -p 06Merge_sites
outdir=$workdir/${species}/*_*/$sample/06Merge_sites
cd $outdir 
sheepCpG=/Shang/Reference/Sheep/ARS-UI_Ramb_v3/biscuit_index/sheep_cpg_all.txt
scriptdir=/Shang/03Project/20251103MethRate/01Script
sbatch $scriptdir/03Sample/06Merge_strand/02Merge_Biscuit.sub $file $sheepCpG $outdir ${sample}.biscuit
done 

###04Encode
workdir=/Shang/03Project/03Analysis/02Samples
species=Sheep
for file in $(ls $workdir/${species}/*_*/*/05*/04*/02*.CG.meth.chr.ENCODE.txt) 
do 
ls $file 
sample02=$(basename $file .CG.meth.chr.ENCODE.txt)
sample=$(echo ${sample02#02})
echo $sample
cd $workdir/${species}/*_*/$sample
mkdir -p 06Merge_sites
outdir=$workdir/${species}/*_*/$sample/06Merge_sites
cd $outdir 
scriptdir=/Shang/03Project/01Script
sbatch $scriptdir/03Sample/06Merge_strand/02Merge_strand.sub $file $outdir ${sample}.encode
done 


