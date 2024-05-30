

# wget "https://www.cog-genomics.org/static/bin/plink2_src_220603.zip"

workDir=/vast/scratch/users/jackson.v/retThickness/GWAS
dataDir=$workDir/geneticData/cleanedEURData

cd $workDir

mkdir -p $workDir/plink
cd $workDir/plink
# version now unavailable...
# wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20221024.zip"
# unzip plink2_linux_x86_64_20221024.zip

wget "https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240105.zip"
unzip plink2_linux_x86_64_20240105.zip

mkdir -p $workDir/sensitivityAnalysesPheno
mkdir -p $workDir/sensitivityAnalysesresults
mkdir -p $workDir/logs
mkdir -p $workDir/scripts

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/phenotypes* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output/phenotypes* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output/pixels.txt .



cat <<- EOF > $workDir/scripts/sensitivityAnalyses.sh
#!/bin/bash

#SBATCH -J sensitivity
#SBATCH -o $workDir/logs/sensitivityAnalyses_%A_%a.log
#SBATCH -t 4:0:0
#SBATCH --mem=6G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-232


locus=\$SLURM_ARRAY_TASK_ID

read locusID chr SNP pixel slice  < <(sed -n \${locus}p loci.txt)


cd $workDir


if [ \$chr -eq 23 ]; then
    
    chr="X"

## smoking as a covariate

./plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs_smoking_EUR.txt \
  --pheno-name \$pixel \
  --covar  $workDir/pheno/covariates_doubleIDs_smoking_EUR.txt \
  --not-covar sex \
  --vif 500  \
  --xchr-model 2 \
  --covar-variance-standardize \
  --glm hide-covar cols=-chrom,-test,-nobs,-err \
  --threads 2 \
  --out $workDir/results/chr\${chr}Pixel\${pixel}_smoking


## no surgery

./plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs_noSurgery_EUR.txt \
  --pheno-name \$pixel \
  --covar  $workDir/pheno/covariates_doubleIDs_noSurgery_EUR.txt \
  --not-covar sex \
  --vif 500  \
  --xchr-model 2 \
  --covar-variance-standardize \
  --glm hide-covar cols=-chrom,-test,-nobs,-err \
  --threads 2 \
  --out $workDir/results/chr\${chr}Pixel\${pixel}_noSurgery


else


## smoking as a covariate

./plink/plink2 \
--pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
--pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs_smoking_EUR.txt \
--pheno-name \$pixel \
--covar  $workDir/pheno/covariates_doubleIDs_smoking_EUR.txt \
--vif 500  \
--covar-variance-standardize \
--glm hide-covar cols=-chrom,-test,-nobs,-err \
--threads 2 \
--out $workDir/results/chr\${chr}Pixel\${pixel}_smoking


## no surgery

./plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs_noSurgery_EUR.txt \
  --pheno-name \$pixel \
  --covar  $workDir/pheno/covariates_doubleIDs_noSurgery_EUR.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=-chrom,-test,-nobs,-err \
  --threads 2 \
  --out $workDir/results/chr\${chr}Pixel\${pixel}_noSurgery

fi

 rsync -av $workDir/results/chr\${chr}*\$pixel.glm.linear /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output

EOF

 sbatch $workDir/scripts/sensitivityAnalyses.sh





cat <<- EOF > $workDir/scripts/sensitivityAnalysesTrans.sh
#!/bin/bash

#SBATCH -J sensitivity
#SBATCH -o $workDir/logs/sensitivityAnalysesTrans_%A_%a.log
#SBATCH -t 2:0:0
#SBATCH --mem=6G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-232


locus=\$SLURM_ARRAY_TASK_ID

read locusID chr SNP pixel slice  < <(sed -n \${locus}p loci.txt)


cd $workDir


if [ \$chr -eq 23 ]; then
    
    chr="X"


./plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs_EUR.txt \
  --pheno-name \$pixel \
  --pheno-quantile-normalize \
  --covar  $workDir/pheno/covariates_doubleIDs_EUR.txt \
  --not-covar sex \
  --vif 500  \
  --xchr-model 2 \
  --covar-variance-standardize \
  --glm hide-covar cols=-chrom,-test,-nobs,-err \
  --threads 2 \
  --out $workDir/results/chr\${chr}Pixel\${pixel}_phenoTrans



else



./plink/plink2 \
--pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
--pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs_EUR.txt \
--pheno-name \$pixel \
--pheno-quantile-normalize \
--covar  $workDir/pheno/covariates_doubleIDs_EUR.txt \
--vif 500  \
--covar-variance-standardize \
--glm hide-covar cols=-chrom,-test,-nobs,-err \
--threads 2 \
--out $workDir/results/chr\${chr}Pixel\${pixel}_phenoTrans



fi

 rsync -av $workDir/results/chr\${chr}*\$pixel_phenoTrans.glm.linear /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output

EOF

 sbatch $workDir/scripts/sensitivityAnalysesTrans.sh










# ./plink/plink2 \
#   --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
#   --pheno $workDir/pheno/phenotypesSlice${slice}_doubleIDs_smoking_EUR.txt \
#   --pheno-name $pixel \
#   --covar  $workDir/pheno/covariates_doubleIDs_smoking_EUR.txt \
#   --vif 500  \
#   --covar-variance-standardize \
#   --glm hide-covar cols=-chrom,-test,-nobs,-err \
#   --threads 2 \
#   --out $workDir/results/chr${chr}Pixel${pixel}_smoking

# ./plink/plink2 \
#   --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
#   --pheno $workDir/pheno/phenotypesSlice${slice}_doubleIDs_noSurgery_EUR.txt \
#   --pheno-name $pixel \
#   --covar  $workDir/pheno/covariates_doubleIDs_noSurgery_EUR.txt \
#   --vif 500  \
#   --covar-variance-standardize \
#   --glm hide-covar cols=-chrom,-test,-nobs,-err \
#   --threads 2 \
#   --out $workDir/results/chr${chr}Pixel${pixel}_noSurgery








