for chr in {12..22}
do
echo $chr

for slice in {1..119} 
do
rm /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr$chr/$slice/chr${chr}Pixel.*.glm.linear
done
done


workDir=/vast/scratch/users/jackson.v/retThickness/GWAS

for chr in {9..21}
do
echo $chr

for slice in {1..119} 
do

 for fi in $workDir/results/chr${chr}/${slice}/chr${chr}Pixel*.glm.linear 
  do
  gzip $fi
 done

  rsync -av $workDir/results/chr${chr}/${slice}/*.gz /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${chr}/${slice}/
  rsync -av $workDir/clumpedResults/chr${chr}/${slice}/*.clumped /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr${chr}/${slice}/

done
done




cat <<- EOF > $workDir/scripts/compress.sh
#!/bin/bash

#SBATCH -J compress
#SBATCH -o $workDir/logs/compress_%A_%a.log
#SBATCH -t 12:0:0
#SBATCH --mem=500MB
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119%15


slice=\$SLURM_ARRAY_TASK_ID


cd $workDir

for chr in {5..1}
do


 for fi in $workDir/results/chr\${chr}/\${slice}/chr\${chr}Pixel*.glm.linear 
  do
  gzip \$fi
 done
  
  mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr\${chr}/
  mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr\${chr}/\${slice}/
  mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr\${chr}/
  mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr\${chr}/\${slice}/

  rsync -av $workDir/results/chr\${chr}/\${slice}/*.gz /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr\${chr}/\${slice}/
  rsync -av $workDir/clumpedResults/chr\${chr}/\${slice}/*.clumped /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr\${chr}/\${slice}/

done

EOF




for chr in 12 13 14 15 16 17 18 19 20 21
for chr in 1 2 3 4 5 6 7 8 9 10 11
for chr in 6 7 8 9 10 11
for chr in 1 2 3 4 5 
do
 sbatch $workDir/scripts/plinkAssoc_chr${chr}.sh
 sleep 10
done


