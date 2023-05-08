

workDir=/vast/scratch/users/jackson.v/retThickness/AFRGWAS
dataDir=/vast/scratch/users/jackson.v/retThickness/filteringGeneticData/cleanedAFRData

mkdir -p $workDir
cd $workDir

mkdir -p $workDir/plink
cd $workDir/plink
wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20221024.zip"
unzip plink2_linux_x86_64_20221024.zip

mkdir -p $workDir/pheno
mkdir -p $workDir/results
mkdir -p $workDir/logs
mkdir -p $workDir/scripts

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/phenotypes* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/chr*sentinelsIDonly_clumpThresh0.001_withOverlap.txt $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt .


  cat <<- EOF > $workDir/scripts/plinkAssoc.sh
#!/bin/bash

#SBATCH -J plink
#SBATCH -o $workDir/logs/plink_%A_%a.log
#SBATCH -t 0:10:0
#SBATCH --mem=2G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119%20


slice=\$SLURM_ARRAY_TASK_ID


cd $workDir

for anc in AFR CSA
do

dataDir=/vast/scratch/users/jackson.v/retThickness/filteringGeneticData/cleaned\${anc}Data

for chr in {1..22}
do

mkdir -p $workDir/results/chr\${chr}/
mkdir -p $workDir/results/chr\${chr}/\${slice}

./plink/plink2 \
  --pfile \$dataDir/plink2Bin/\${anc}_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs_\${anc}.txt \
  --covar  $workDir/pheno/covariates_doubleIDs_\${anc}.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=-ref,-alt,-test,-nobs,-err,+a1freq \
  --extract $workDir/pheno/chr\${chr}sentinelsIDonly_clumpThresh0.001_withOverlap.txt \
  --threads 2 \
  --out $workDir/results/chr\${chr}/\${slice}/chr\${chr}Pixel\${anc}


done

chr=X

mkdir -p $workDir/results/chr\${chr}/
mkdir -p $workDir/results/chr\${chr}/\${slice}

./plink/plink2 \
  --pfile \$dataDir/plink2Bin/\${anc}_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs_\${anc}.txt \
  --covar  $workDir/pheno/covariates_doubleIDs_\${anc}.txt \
  --not-covar sex \
  --vif 500  \
  --xchr-model 2 \
  --covar-variance-standardize \
  --glm hide-covar cols=-ref,-alt,-test,-nobs,-err,+a1freq \
  --extract $workDir/pheno/chr\${chr}sentinelsIDonly_clumpThresh0.001_withOverlap.txt \
  --threads 2 \
  --out $workDir/results/chr\${chr}/\${slice}/chr\${chr}Pixel\${anc}

done

EOF

sbatch $workDir/scripts/plinkAssoc.sh


## format results : formattingNonEurResults.R
## plotting result comparisons : 











## testing
# mkdir -p $workDir/results/chr${chr}/

# ./plink/plink2 \
#   --pfile $dataDir/plink2Bin/AFR_minMaf0.005_minInfo0.8_chr${chr} \
#   --pheno $workDir/pheno/phenotypesSlice${slice}_doubleIDs_AFR.txt \
#   --covar  $workDir/pheno/covariates_doubleIDs_AFR.txt \
#   --vif 500  \
#   --covar-variance-standardize \
#   --glm hide-covar cols=-chrom,-ref,-alt,-test,-nobs,-err \
#   --extract $workDir/pheno/chr${chr}sentinelsIDonly_clumpThresh0.001_withOverlap.txt \
#   --threads 2 \
#   --out $workDir/results/chr${chr}/${slice}/chr${chr}PixelAFR


# ./plink/plink2 \
#   --pfile $dataDir/plink2Bin/AFR_minMaf0.005_minInfo0.8_chr${chr} \
#   --pheno $workDir/pheno/phenotypesSlice${slice}_doubleIDs_AFR.txt \
#   --covar  $workDir/pheno/covariates_doubleIDs_AFR.txt \
#   --vif 500  \
#   --covar-variance-standardize \
#   --glm hide-covar cols=-chrom,-ref,-alt,-test,-nobs,-err \
#   --extract $workDir/pheno/chr${chr}sentinelsIDonly_clumpThresh0.001_withOverlap.txt \
#   --threads 2 \
#   --out $workDir/results/chr${chr}/${slice}/chr${chr}PixelCSA

