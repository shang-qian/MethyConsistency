#!/bin/bash

####
for sample in SAMN39476118 SAMN39476119
do 
workdir=/Shang/03Project/03Analysis/02Samples
species=Cattle
breed=Angus
proID=PRJNA1054159
######01Bismark
mkdir -p $workdir/${species}/${breed}_${proID}/$sample/01Bismark/03Call
cd $workdir/${species}/${breed}_${proID}/$sample/01Bismark/03Call
resultdir=$workdir/${species}/${breed}_${proID}/$sample/01Bismark/03Call
inputdir=$workdir/${species}/${breed}_${proID}/$sample/01Bismark/01Bowtie2
dedupbam=$(ls $inputdir/*bismark_bt2_pe.deduplicated.bam)
REF=/Shang/Reference/Cattle/ARS-UCD_v2/Bismark_bowtie2
scriptdir=/Shang/03Project/01Script/03Sample/03Caller/sbatch_sample
sbatch $scriptdir/01Bismark_call.sub $dedupbam $REF $resultdir
######02Bwameth
mkdir -p $workdir/${species}/${breed}_${proID}/$sample/03Bwameth/03Call
cd $workdir/${species}/${breed}_${proID}/$sample/03Bwameth/03Call
resultdir=$workdir/${species}/${breed}_${proID}/$sample/03Bwameth/03Call
inputdir=$workdir/${species}/${breed}_${proID}/$sample/03Bwameth/01Bwamem2
reference=/Shang/Reference/Cattle/ARS-UCD_v2/bwameth_index/cattle_v2_genome.fa 
scriptdir=/Shang/03Project/01Script/03Sample/03Caller/sbatch_sample
sbatch $scriptdir/02Bwameth_call.sub $sample $inputdir $reference $resultdir
done



