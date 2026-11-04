
workdir=/Shang/03Project/03Analysis
mkdir -p $workdir/06Species/00Chain
cd $workdir/06Species/00Chain


module load minimap2

minimap2 \
-x asm20 \
-c \
--cs \
-t 250 \
--secondary=no \
$cattle_ref \
$sheep_ref \
> 02sheep_to_cattle.paf

conda activate sheep_final
transanno minimap2chain \
    03sheep_to_cattle.paf \
    --output 04sheep_to_cattle.chain


workdir=/Shang/03Project/03Analysis
mkdir -p $workdir/06Species/01Species/Sites
cd $workdir/06Species/01Species/Sites



####cattle and sheep
sheep_H_Ci=/Shang/03Project/03Analysis/04MCI/Sheep/01Global/03_All_Sheep_Consistency_High_u_Ci90_95_100.txt 
cattle_H_Ci=/Shang/03Project/03Analysis/04MCI/Cattle/01Global/03_All_Cattle_Consistency_High_u_Ci90_95_100.txt

for threshold in 0.99 #0.95 # 0.96 0.97 0.98 0.99 1
do 
echo $threshold 
awk -v thre=$threshold 'NR>1&&$7>=thre {print $1"\t"$2"\t"$2+1"\t"$0}' $sheep_H_Ci |\
awk 'BEGIN{OFS="\t"}
{
    $2=sprintf("%.0f",$2)
    $3=sprintf("%.0f",$3)
    print
}' > 01Sheep_C90_Prop_${threshold}.bed

awk -v thre=$threshold 'NR>1&&$7>=thre {print $1"\t"$2"\t"$2+1"\t"$0}' $cattle_H_Ci |\
awk 'BEGIN{OFS="\t"}
{
    $2=sprintf("%.0f",$2)
    $3=sprintf("%.0f",$3)
    print
}' > 02Cattle_C90_Prop_${threshold}.bed
done 

transanno liftbed \
    --chain ../../00Chain/04sheep_to_cattle.chain \
    --output 03Sheep2cattle_C9_P99.bed \
    --failed 03unmapped_99.bed \
    01Sheep_C90_Prop_0.99.bed

module load bedtools2
bedtools intersect -loj -a 02Cattle_C90_Prop_0.99.bed  -b 03Sheep2cattle_C9_P99.bed |awk '$13!="." {print $0}' > 04Cattle_Sheep_matched_MS_99.bed

