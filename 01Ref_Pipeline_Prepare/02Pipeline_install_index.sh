#####Reference genome index for 4 pipelines
cd /Shang/01Envs/WGBS_pip

###1. bismark
#Cattle:
bismark_genome_preparation --bowtie2 --parallel 8 /Shang/03Project/01Ref/Cattle/ARS-UCD2/Bismark_bowtie2  
#Sheep:
bismark_genome_preparation --bowtie2 --parallel 8 /Shang/03Project/01Ref/Sheep/ARS-UI_Ramb_v3/Bismark_bowtie2  

    
###2. Bwameth
#Cattle: 
cd /Shang/03Project/01Ref/Cattle/ARS-UCD2/bwameth_index
bwameth.py index cattle_v2_genome.fa
#sheep:
cd /Shang/03Project/01Ref/Sheep/ARS-UI_Ramb_v3/bwameth_index
bwameth.py index sheep_v3_genome.fa


###3.Biscuit
#Cattle:
cd /Shang/03Project/01Ref/Cattle/ARS-UCD2/biscuit_index
biscuit index cattle_v2_genome.fa
#sheep:
cd /cluster/pixstor/mckaysd-lab/Reference/Sheep/ARS-UI_Ramb_v3/biscuit_index
biscuit index sheep_v3_genome.fa


###4.WGBS-pipeline
#download ENCODE_wgbs-pipeline
git clone https://github.com/ENCODE-DCC/wgbs-pipeline
#Cattle:
cd /cluster/pixstor/mckaysd-lab/Reference/Cattle/ARS-UCD_v2/Encode_index
module load samtools
samtools faidx cattle_v2_genome.fa
#Sheep:
cd /cluster/pixstor/mckaysd-lab/Reference/Sheep/ARS-UI_Ramb_v3/Encode_index
module load samtools
samtools faidx sheep_v3_genome.fa

