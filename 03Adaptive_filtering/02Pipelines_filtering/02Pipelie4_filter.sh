#!/bin/bash

workdir=/Shang/03Project/03Analysis/02Samples
###01. Cattle
species=Cattle
for file in $(ls $workdir/${species}/*_*/*/06Merge_sites/*.bismark.merged.txt)
do 
ls $file 
sample=$(basename $file .bismark.merged.txt)
echo $sample
cd $workdir/${species}/*_*/$sample
mkdir -p 08Pipeline_filter
outdir=$workdir/${species}/*_*/$sample/08Pipeline_filter
cd $outdir 

inputdir=$workdir/${species}/*_*/$sample/07Strand_filter
ls $inputdir/${sample}_06bismark.strand.filtering.txt
bismark=$inputdir/${sample}_06bismark.strand.filtering.txt
bwameth=$inputdir/${sample}_06bwameth.strand.filtering.txt
biscuit=$inputdir/${sample}_06biscuit.strand.filtering.txt
encode=$inputdir/${sample}_06encode.strand.filtering.txt

refcpg=/Shang/Reference/Cattle/ARS-UCD_v2/biscuit_index/cattle_cpg_all.txt
###02 filtering for integrated pipelines
scriptdir=/Shang/03Project/01Script
sbatch $scriptdir/03Sample/07Distribution_filter/00TPipeline_methods4.sub $sample $bismark $bwameth $biscuit $encode $outdir $refcpg
done 



###02. Sheep 
workdir=/Shang/03Project/03Analysis/02Samples
###02. Sheep
species=Sheep
for file in $(ls $workdir/${species}/*_*/*/06Merge_sites/*.bismark.merged.txt)
do 
ls $file 
sample=$(basename $file .bismark.merged.txt)
echo $sample
cd $workdir/${species}/*_*/$sample
mkdir -p 08Pipeline_filter
outdir=$workdir/${species}/*_*/$sample/08Pipeline_filter
cd $outdir 

inputdir=$workdir/${species}/*_*/$sample/07Strand_filter
ls $inputdir/${sample}_06bismark.strand.filtering.txt
bismark=$inputdir/${sample}_06bismark.strand.filtering.txt
bwameth=$inputdir/${sample}_06bwameth.strand.filtering.txt
biscuit=$inputdir/${sample}_06biscuit.strand.filtering.txt
encode=$inputdir/${sample}_06encode.strand.filtering.txt

refcpg=/Shang/Reference/Sheep/ARS-UI_Ramb_v3/biscuit_index/sheep_cpg_all.txt
###02 filtering for integrated pipelines
scriptdir=/Shang/03Project/20251103MethRate/01Script
sbatch $scriptdir/03Sample/07Distribution_filter/00TPipeline_methods4.sub $sample $bismark $bwameth $biscuit $encode $outdir $refcpg
done 


