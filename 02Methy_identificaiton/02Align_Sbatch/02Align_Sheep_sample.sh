#!/bin/bash

####
for i in SAMEA104494974  SAMEA104495017 
do
workdir=/Shang/03Project/03Analysis/02Samples
species=Sheep
breed=Ram
proID=PRJEB39178
sample=$i

inputdir=/Shang/02Samples/$species/${breed}_${proID}/$sample
scriptdir=/Shang/03Project/01Script/03Sample/02Alignment/sbatch
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

###04WGBS-pipeline
mkdir -p $workdir/${species}/${breed}_${proID}/$sample/05ENCODE_pip/01Align
cd $workdir/${species}/${breed}_${proID}/$sample/05ENCODE_pip/01Align
sbatch $scriptdir/05WGBSpipe.sub $sample $inputdir $workdir/${species}/${breed}_${proID}/$sample/05ENCODE_pip/01Align $sample
done


