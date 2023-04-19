

# wget "https://www.cog-genomics.org/static/bin/plink2_src_220603.zip"

workDir=/vast/scratch/users/jackson.v/retThickness/GWAS
dataDir=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedEURData

cd $workDir

mkdir -p $workDir/plink
cd $workDir/plink
wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20221024.zip"
unzip plink2_linux_x86_64_20221024.zip

mkdir -p $workDir/pheno
mkdir -p $workDir/results
mkdir -p $workDir/clumpedResults
mkdir -p $workDir/logs
mkdir -p $workDir/scripts

mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/EURpixelWiseResultsDec2022/
mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/EURpixelWiseResultsDec2022/results/
mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/EURpixelWiseResultsDec2022/clumpedResults/
  
  
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/phenotypes* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt .



for chr in {1..22}
do
  cat <<- EOF > $workDir/scripts/plinkAssoc_chr${chr}.sh
#!/bin/bash

#SBATCH -J plink-${chr}
#SBATCH -o $workDir/logs/plinkChr${chr}_%A_%a.log
#SBATCH -t 48:0:0
#SBATCH --mem=6G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119%20


slice=\$SLURM_ARRAY_TASK_ID


cd $workDir


mkdir -p $workDir/results/chr${chr}/\${slice}
mkdir -p $workDir/clumpedResults/chr${chr}/\${slice}

./plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs_EUR.txt \
  --covar  $workDir/pheno/covariates_doubleIDs_EUR.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=-chrom,-ref,-alt,-test,-nobs,-err \
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
    --bfile $dataDir/plinkBin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
    --clump $workDir/results/chr${chr}/\${slice}/chr${chr}Pixel.\${pixel}.glm.linear \
    --clump-snp-field ID \
    --clump-p1 5e-8 \
    --clump-r2 0.1 \
    --clump-p2 5e-5 \
    --clump-kb 500 \
    --threads 2 \
    --out $workDir/clumpedResults/chr${chr}/\${slice}/chr${chr}Pixel.\${pixel}
fi
done

 for fi in $workDir/results/chr${chr}/\${slice}/chr${chr}Pixel*.glm.linear 
  do
  gzip \$fi
 done
  
  mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${chr}/
  mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${chr}/\${slice}/
  mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr${chr}/
  mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr${chr}/\${slice}/

  rsync -av $workDir/results/chr${chr}/\${slice}/*.gz /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${chr}/\${slice}/
  rsync -av $workDir/clumpedResults/chr${chr}/\${slice}/*.clumped /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr${chr}/\${slice}/

EOF

 sbatch $workDir/scripts/plinkAssoc_chr${chr}.sh
 sleep 10

done




## re-run clumping with different prameters
mkdir -p $workDir/tmp
for slice in {1..119}
do
awk -v y="$slice" ' $2==y ' /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt > $workDir/tmp/${slice}_pixels.txt
done

for chr in {1..22}
do
  cat <<- EOF > $workDir/scripts/plinkClump_chr${chr}.sh
#!/bin/bash

#SBATCH -J clump-${chr}
#SBATCH -o $workDir/logs/clumpChr${chr}_%A_%a.log
#SBATCH -t 4:0:0
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119%15


slice=\$SLURM_ARRAY_TASK_ID


cd $workDir


mkdir -p $workDir/results/chr${chr}/\${slice}
mkdir -p $workDir/clumpedResults/chr${chr}/\${slice}

## load plink v1.9
module load plink


nPix=\$(wc -l  $workDir/tmp/\${slice}_pixels.txt | cut -d " " -f 1)


for pix in \$(seq 1 \$nPix)
do

read pixel y x < <(sed -n \${pix}p $workDir/tmp/\${slice}_pixels.txt)

    echo 'Clumping Pixel '\$pixel''

    zcat  /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${chr}/\${slice}/chr${chr}Pixel.\${pixel}.glm.linear.gz > $workDir/tmp/chr${chr}_\${pixel}.txt

    plink \
    --bfile $dataDir/plinkBin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
    --clump $workDir/tmp/chr${chr}_\${pixel}.txt \
    --clump-snp-field ID \
    --clump-p1 5e-8 \
    --clump-r2 0.001 \
    --clump-p2 5e-5 \
    --clump-kb 5000 \
    --threads 2 \
    --out $workDir/clumpedResults/chr${chr}/\${slice}/chr${chr}Pixel.\${pixel}_thresh0.001

rm $workDir/tmp/chr${chr}_\${pixel}.txt 

done

  rsync -av $workDir/clumpedResults/chr${chr}/\${slice}/*_thresh0.001.clumped /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr${chr}/\${slice}/

EOF

 sbatch $workDir/scripts/plinkClump_chr${chr}.sh
 sleep 10

done





mkdir -p $workDir/tmp
for slice in {1..119}
do
awk -v y="$slice" ' $2==y ' /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt > $workDir/tmp/${slice}_pixels.txt
done

for chr in {1..22}
do
  cat <<- EOF > $workDir/scripts/plinkClump_chr${chr}.sh
#!/bin/bash

#SBATCH -J clump-${chr}
#SBATCH -o $workDir/logs/clumpChr${chr}_overlap_%A_%a.log
#SBATCH -t 4:0:0
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119%15


slice=\$SLURM_ARRAY_TASK_ID


cd $workDir


mkdir -p $workDir/results/chr${chr}/\${slice}
mkdir -p $workDir/clumpedResults/chr${chr}/\${slice}

## load plink v1.9
module load plink


nPix=\$(wc -l  $workDir/tmp/\${slice}_pixels.txt | cut -d " " -f 1)


for pix in \$(seq 1 \$nPix)
do

read pixel y x < <(sed -n \${pix}p $workDir/tmp/\${slice}_pixels.txt)

    echo 'Clumping Pixel '\$pixel''

    zcat  /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${chr}/\${slice}/chr${chr}Pixel.\${pixel}.glm.linear.gz > $workDir/tmp/chr${chr}_\${pixel}.txt

    plink \
    --bfile $dataDir/plinkBin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
    --clump $workDir/tmp/chr${chr}_\${pixel}.txt \
    --clump-snp-field ID \
    --clump-p1 5e-8 \
    --clump-r2 0.001 \
    --clump-p2 5e-5 \
    --clump-kb 5000 \
    --clump-allow-overlap \
    --threads 2 \
    --out $workDir/clumpedResults/chr${chr}/\${slice}/chr${chr}Pixel.\${pixel}_thresh0.001_withOverlap

rm $workDir/tmp/chr${chr}_\${pixel}.txt 

done

  rsync -av $workDir/clumpedResults/chr${chr}/\${slice}/*_thresh0.001_withOverlap.clumped /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr${chr}/\${slice}/

EOF

 sbatch $workDir/scripts/plinkClump_chr${chr}.sh
 sleep 10

done






chr=X

  cat <<- EOF > $workDir/scripts/plinkAssoc_chr${chr}.sh
#!/bin/bash

#SBATCH -J plink-${chr}
#SBATCH -o $workDir/logs/plinkChr${chr}_%A_%a.log
#SBATCH -t 48:0:0
#SBATCH --mem=6G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119%20


slice=\$SLURM_ARRAY_TASK_ID


cd $workDir


mkdir -p $workDir/results/chr${chr}/\${slice}
mkdir -p $workDir/clumpedResults/chr${chr}/\${slice}

./plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs_EUR.txt \
  --covar  $workDir/pheno/covariates_doubleIDs_EUR.txt \
  --not-covar sex \
  --vif 500  \
  --xchr-model 2 \
  --covar-variance-standardize \
  --glm hide-covar cols=-chrom,-ref,-alt,-test,-nobs,-err \
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
    --bfile $dataDir/plinkBin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
    --clump $workDir/results/chr${chr}/\${slice}/chr${chr}Pixel.\${pixel}.glm.linear \
    --clump-snp-field ID \
    --clump-p1 5e-8 \
    --clump-r2 0.001 \
    --clump-p2 5e-5 \
    --clump-kb 5000 \
    --clump-allow-overlap \
    --threads 2 \
    --out $workDir/clumpedResults/chr${chr}/\${slice}/chr${chr}Pixel.\${pixel}_thresh0.001_withOverlap
fi
done

 for fi in $workDir/results/chr${chr}/\${slice}/chr${chr}Pixel*.glm.linear 
  do
  gzip \$fi
 done
  
  mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${chr}/
  mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${chr}/\${slice}/
  mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr${chr}/
  mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr${chr}/\${slice}/

  rsync -av $workDir/results/chr${chr}/\${slice}/*.gz /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${chr}/\${slice}/
  rsync -av $workDir/clumpedResults/chr${chr}/\${slice}/*.clumped /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr${chr}/\${slice}/

EOF

 sbatch $workDir/scripts/plinkAssoc_chr${chr}.sh




























cd $workDir


# mkdir -p $workDir/results/chr${chr}{slice}
# mkdir -p $workDir/clumpedResults/chr${chr}/\${slice}

# ## load plink v1.9
# module load plink


# nPix=$(wc -l  $workDir/tmp/${slice}_pixels.txt | cut -d " " -f 1)


# for pix in $(seq 1 $nPix)
# do

# read pixel y x < <(sed -n ${pix}p $workDir/tmp/${slice}_pixels.txt)

#     echo 'Clumping Pixel '$pixel''

#     zcat  /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr${chr}/${slice}/chr${chr}Pixel.${pixel}.glm.linear.gz > $workDir/tmp/chr${chr}_${pixel}.txt

#     plink \
#     --bfile $dataDir/plinkBin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
#     --clump $workDir/tmp/chr${chr}_${pixel}.txt \
#     --clump-snp-field ID \
#     --clump-p1 5e-8 \
#     --clump-r2 0.001 \
#     --clump-p2 5e-5 \
#     --clump-kb 5000 \
#     --threads 2 \
#     --out $workDir/clumpedResults/chr${chr}/${slice}/chr${chr}Pixel.${pixel}_thresh0.001

# rm $workDir/tmp/chr${chr}_${pixel}.txt 

# done

#   rsync -av $workDir/clumpedResults/chr${chr}/${slice}/*_thresh0.001.clumped /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/clumpedResults/chr${chr}/${slice}/
