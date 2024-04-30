workDir=/vast/scratch/users/jackson.v/retThickness/GWAS
dataDir=$workDir/geneticData/cleanedEURData

cd $workDir

mkdir -p $workDir/plink
cd $workDir/plink
wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20221024.zip"
unzip plink2_linux_x86_64_20221024.zip

mkdir -p $workDir/pheno
mkdir -p $workDir/results
mkdir -p $workDir/logs
mkdir -p $workDir/scripts

  
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/covariates*EUR.txt $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/phenotypes*EUR.txt $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt .



cat <<- EOF > $workDir/scripts/conditionalAnalyses.sh
#!/bin/bash

#SBATCH -J conditional
#SBATCH -o $workDir/logs/conditionalAnalyses_%A_%a.log
#SBATCH -t 2:0:0
#SBATCH --mem=6G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-224


locus=\$SLURM_ARRAY_TASK_ID

read locusID chr SNP pixel slice  < <(sed -n \${locus}p loci.txt)


cd $workDir


if [ \$chr -eq 23 ]; then
    
    chr="X"


./plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs_EUR.txt \
  --pheno-name \$pixel \
  --covar  $workDir/pheno/covariates_doubleIDs_EUR.txt \
  --not-covar sex \
  --condition \$SNP \
  --vif 500  \
  --xchr-model 2 \
  --covar-variance-standardize \
  --glm hide-covar cols=-chrom,-test,-nobs,-err \
  --threads 2 \
  --out $workDir/results/chr\${chr}Pixel\${pixel}_conditional_\$SNP



else


./plink/plink2 \
--pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
--pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs_EUR.txt \
--pheno-name \$pixel \
--covar  $workDir/pheno/covariates_doubleIDs_EUR.txt \
--condition \$SNP \
--vif 500  \
--covar-variance-standardize \
--glm hide-covar cols=-chrom,-test,-nobs,-err \
--threads 2 \
--out $workDir/results/chr\${chr}Pixel\${pixel}_conditional_\$SNP


fi

 rsync -av $workDir/results/chr\${chr}Pixel\${pixel}_conditional_\${SNP}*.glm.linear /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output

EOF

 sbatch $workDir/scripts/conditionalAnalyses.sh


#  rsync -av $workDir/results/chr*Pixel*_conditional_*.glm.linear /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output
