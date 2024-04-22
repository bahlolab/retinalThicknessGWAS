

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

mkdir -p $workDir/sensitivityAnalysesResultsFPC
mkdir -p $workDir/logs
mkdir -p $workDir/scripts

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output/FPCphenotypes* $workDir/pheno



cat <<- EOF > $workDir/scripts/sensitivityAnalysesFPCs.sh
#!/bin/bash

#SBATCH -J sensitivity
#SBATCH -o $workDir/logs/sensitivityAnalysesFPC_%A_%a.log
#SBATCH -t 3:0:0
#SBATCH --mem=6G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-23

cd $workDir



if [ \$SLURM_ARRAY_TASK_ID -eq 23 ]; then
    
    chr="X"
    mkdir -p $workDir/sensitivityAnalysesResultsFPC/chr\${chr}/


    ## smoking as a covariate

./plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/FPCphenotypes_doubleIDs_smoking_EUR.txt \
  --pheno-name fpc1-fpc6 \
  --covar  $workDir/pheno/covariates_doubleIDs_smoking_EUR.txt \
  --not-covar sex \
  --vif 500  \
  --xchr-model 2 \
  --covar-variance-standardize \
  --glm hide-covar cols=-chrom,-test,-nobs,-err \
  --threads 2 \
  --out $workDir/results/chr\${chr}_smoking


## no surgery

./plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/FPCphenotypes_doubleIDs_noSurgery_EUR.txt \
  --pheno-name fpc1-fpc6 \
  --covar  $workDir/pheno/covariates_doubleIDs_noSurgery_EUR.txt \
  --not-covar sex \
  --vif 500  \
  --xchr-model 2 \
  --covar-variance-standardize \
  --glm hide-covar cols=-chrom,-test,-nobs,-err \
  --threads 2 \
  --out $workDir/results/chr\${chr}_noSurgery


    else
  
    chr=\$SLURM_ARRAY_TASK_ID
    mkdir -p $workDir/sensitivityAnalysesResultsFPC/chr\${chr}/


## smoking as a covariate

./plink/plink2 \
--pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
--pheno $workDir/pheno/FPCphenotypes_doubleIDs_smoking_EUR.txt \
  --pheno-name fpc1-fpc6 \
--covar  $workDir/pheno/covariates_doubleIDs_smoking_EUR.txt \
--vif 500  \
--covar-variance-standardize \
--glm hide-covar cols=-chrom,-test,-nobs,-err \
--threads 2 \
--out $workDir/results/chr\${chr}_smoking


## no surgery

./plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/FPCphenotypes_doubleIDs_noSurgery_EUR.txt \
  --pheno-name fpc1-fpc6 \
  --covar  $workDir/pheno/covariates_doubleIDs_noSurgery_EUR.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=-chrom,-test,-nobs,-err \
  --threads 2 \
  --out $workDir/results/chr\${chr}_noSurgery

fi

 rsync -av $workDir/results/chr\${chr}*fpc*.glm.linear /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/sensitivityAnalysesGWAS/output

EOF

 sbatch $workDir/scripts/sensitivityAnalysesFPCs.sh










