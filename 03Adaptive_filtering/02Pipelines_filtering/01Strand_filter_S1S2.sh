#!/bin/bash

workdir=/Shang/03Project/03Analysis/02Samples
###01. Cattle strand filtering
species=Cattle
for file in $(ls $workdir/${species}/*_*/*/06Merge_sites/*.bismark.merged.txt)
do 
ls $file 
sample=$(basename $file .bismark.merged.txt)
echo $sample
cd $workdir/${species}/*_*/$sample
mkdir -p 07Strand_filter
outdir=$workdir/${species}/*_*/$sample/07Strand_filter
cd $outdir 

inputdir=$workdir/${species}/*_*/$sample/06Merge_sites
ls $inputdir/$sample.bismark.merged.txt
bismark=$inputdir/$sample.bismark.merged.txt
bwameth=$inputdir/$sample.bwameth.merged.txt
biscuit=$inputdir/$sample.biscuit.merged.txt
encode=$inputdir/$sample.encode.merged.txt
###01 filtering for single strand of each pipeline
scriptdir=/Shang/03Project/01Script
sbatch $scriptdir/03Sample/07Distribution_filter/00Strand_S1S2.sub $sample $bismark $bwameth $biscuit $encode $outdir
done 


###02. Sheep strand filtering
workdir=/Shang/03Project/20251103MethRate/03Analysis/02Samples
###02. Sheep
species=Sheep
for file in $(ls $workdir/${species}/*_*/*/06Merge_sites/*.bismark.merged.txt)
do 
ls $file 
sample=$(basename $file .bismark.merged.txt)
echo $sample
cd $workdir/${species}/*_*/$sample
mkdir -p 07Strand_filter
outdir=$workdir/${species}/*_*/$sample/07Strand_filter
cd $outdir 

inputdir=$workdir/${species}/*_*/$sample/06Merge_sites
ls $inputdir/$sample.bismark.merged.txt
bismark=$inputdir/$sample.bismark.merged.txt
bwameth=$inputdir/$sample.bwameth.merged.txt
biscuit=$inputdir/$sample.biscuit.merged.txt
encode=$inputdir/$sample.encode.merged.txt
###01 filtering for single strand of each pipeline
scriptdir=/Shang/03Project/01Script
sbatch $scriptdir/03Sample/07Distribution_filter/00Strand_S1S2.sub $sample $bismark $bwameth $biscuit $encode $outdir
done 


