#!/bin/bash

for sample in SAMN39476118 SAMN39476119
do 
workdir=/Shang/03Project/03Analysis/02Samples
species=Cattle
breed=Angus
proID=PRJNA1054159

inputdir=/Shang/02Samples/$species/${breed}_${proID}/$sample
scriptdir=/Shang/03Project/01Script/03Sample/01Align_Sbatch
ls $inputdir


##Alignment:01Bismark
mkdir -p $workdir/${species}/${breed}_${proID}/$sample/01Bismark/01Bowtie2
cd $workdir/${species}/${breed}_${proID}/$sample/01Bismark/01Bowtie2
sbatch $scriptdir/01Bismark.sub $sample $inputdir $workdir/${species}/${breed}_${proID}/$sample/01Bismark/01Bowtie2

##02bwameth
mkdir -p $workdir/${species}/${breed}_${proID}/$sample/03Bwameth/01Bwamem2
cd $workdir/${species}/${breed}_${proID}/$sample/03Bwameth/01Bwamem2
sbatch $scriptdir/03Bwameth.sub $sample $inputdir $workdir/${species}/${breed}_${proID}/$sample/03Bwameth/01Bwamem2

###03biscuit
mkdir -p $workdir/${species}/${breed}_${proID}/$sample/04Biscuit/01Align
cd $workdir/${species}/${breed}_${proID}/$sample/04Biscuit/01Align
sbatch $scriptdir/04Biscuit.sub $sample $inputdir $workdir/${species}/${breed}_${proID}/$sample/04Biscuit/01Align

###04 WGBS-pipeline
mkdir -p $workdir/${species}/${breed}_${proID}/$sample/05ENCODE_pip/01Align
cd $workdir/${species}/${breed}_${proID}/$sample/05ENCODE_pip/01Align
sbatch $scriptdir/05WGBSpipe.sub $sample $inputdir $workdir/${species}/${breed}_${proID}/$sample/05ENCODE_pip/01Align $sample

done


