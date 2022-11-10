

# wget "https://www.cog-genomics.org/static/bin/plink2_src_220603.zip"

workDir=/vast/scratch/users/jackson.v/retThickness/GWAS
dataDir=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedWhiteBritUnrelatedData

cd $workDir

# mkdir -p $workDir/plink
# cd $workDir/plink
# wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20220814.zip"
# unzip plink2_linux_x86_64_20220814.zip

mkdir -p $workDir/pheno
mkdir -p $workDir/filteredGeno
mkdir -p $workDir/results
mkdir -p $workDir/clumpedResults
mkdir -p $workDir/logs
mkdir -p $workDir/scripts


rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/phenotypes* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt .

## convert data
# for chr in 1 2 3 5 7 9 10 19
for chr in 4 6 8 11 12 13 14 15 16 17 18 20 21 22
do

./plink/plink2 \
  --bgen $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr${chr}.bgen ref-first \
  --sample $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr${chr}.sample \
  --make-bed \
  --threads 2 \
  --mac 500 \
  --new-id-max-allele-len 30 truncate \
  --out $workDir/filteredGeno/ukbb_minMAC500_minInfo0.8_chr${chr} \
  --set-all-var-ids @:#_\$r_\$a

done


for chr in {1..22}
do
  cat <<- EOF > $workDir/scripts/plinkAssoc_chr${chr}.sh
#!/bin/sh

#SBATCH -J plink-${chr}
#SBATCH -o $workDir/logs/plinkChr${chr}_%A_%a.log
#SBATCH -t 48:0:0
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119%20


slice=\$SLURM_ARRAY_TASK_ID


cd $workDir


mkdir -p $workDir/results/chr${chr}/\${slice}

./plink/plink2 \
  --bfile $workDir/filteredGeno/ukbb_minMAC500_minInfo0.8_chr${chr} \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs.txt \
  --covar  $workDir/pheno/covariates_doubleIDs.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=+a1count,+a1freq \
  --threads 2 \
  --out $workDir/results/chr${chr}/\${slice}/chr${chr}Pixel


## load plink v1.9
module load plink

nPix=\$(wc -l pixels.txt | cut -d " " -f 1)

for pix in \$(seq 1 \$nPix)
do

  read pixel y x < <(sed -n \${pix}p pixels.txt)

  if [ \$y -eq \$slice ]
  then
    echo 'Clumping Pixel '\$pixel''

    plink \
    --bfile $workDir/filteredGeno/ukbb_minMAC500_minInfo0.8_chr${chr}\
    --clump $workDir/results/chr${chr}/\${slice}/chr${chr}Pixel.\${pixel}.glm.linear \
    --clump-snp-field ID \
    --clump-p1 5e-8 \
    --clump-r2 0.1 \
    --clump-p2 5e-5 \
    --clump-kb 500 \
    --threads 2 \
    --out $workDir/clumpedResults/chr${chr}/chr${chr}Pixel.\${pixel}
  done

EOF

done


for chr in 1 2 3 5 7 9 10 19
# for chr in 4 6 8 11 12 13 14 15 16 17 18 20 21 22
do
 sbatch $workDir/scripts/plinkAssoc_chr${chr}.sh
 sleep 1
done


mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsOct2022
mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseClumpedResultsOct2022

rsync -av $workDir/results/* /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsOct2022
rsync -av $workDir/clumpedResults/* /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseClumpedResultsOct2022

