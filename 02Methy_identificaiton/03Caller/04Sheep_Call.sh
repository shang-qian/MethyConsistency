#!/bin/bash

####
######01bwameth
for i in SAMEA104494974 #SAMEA104495017  SAMEA104495022  SAMEA104495023  SAMEA104495048  SAMEA104495052  SAMEA104495058  SAMEA104495065
do
workdir=/Shang/03Project/03Analysis/02Samples
species=Sheep
breed=Ram
proID=PRJEB39178
sample=$i
##Bismark
mkdir -p $workdir/${species}/${breed}_${proID}/$sample/01Bismark/03Call
cd $workdir/${species}/${breed}_${proID}/$sample/01Bismark/03Call
resultdir=$workdir/${species}/${breed}_${proID}/$sample/01Bismark/03Call
inputdir=$workdir/${species}/${breed}_${proID}/$sample/01Bismark/01Bowtie2
dedupbam=$(ls $inputdir/*bismark_bt2_pe.deduplicated.bam)
REF=/Shang/Reference/Sheep/ARS-UI_Ramb_v3/Bismark_bowtie2  
scriptdir=/Shang/03Project/01Script/03Sample/03Caller/sbatch_sample
sbatch $scriptdir/01Bismark_call.sub $dedupbam $REF $resultdir
##Bwameth
mkdir -p $workdir/${species}/${breed}_${proID}/$sample/03Bwameth/03Call
cd $workdir/${species}/${breed}_${proID}/$sample/03Bwameth/03Call
resultdir=$workdir/${species}/${breed}_${proID}/$sample/03Bwameth/03Call
inputdir=$workdir/${species}/${breed}_${proID}/$sample/03Bwameth/01Bwamem2
reference=/Shang/Reference/Sheep/ARS-UI_Ramb_v3/bwameth_index/sheep_v3_genome.fa
scriptdir=/Shang/03Project/20251103MethRate/01Script/03Sample/03Caller/sbatch_sample
sbatch $scriptdir/02Bwameth_call.sub $sample $inputdir $reference $resultdir
done

