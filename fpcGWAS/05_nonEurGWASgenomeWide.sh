
workDir=/vast/scratch/users/jackson.v/retThickness/fpcGWASnoExclusions/nonEURGWAS

mkdir -p $workDir
cd $workDir

mkdir -p $workDir/plink
cd $workDir/plink
# version no longer available...
# wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20221024.zip"
# unzip plink2_linux_x86_64_20221024.zip
wget "https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240105.zip"
unzip plink2_linux_x86_64_20240105.zip

mkdir -p $workDir/pheno
mkdir -p $workDir/results
mkdir -p $workDir/logs
mkdir -p $workDir/scripts

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/*.txt $workDir/pheno


  cat <<- EOF > $workDir/scripts/plinkAssoc.sh
#!/bin/bash
#SBATCH -J plink
#SBATCH -o $workDir/logs/plink_%A_%a.log
#SBATCH -t 0:6:0
#SBATCH --mem=2G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-23

cd $workDir
dataDir=/vast/scratch/users/jackson.v/retThickness/GWAS/geneticData/

if [[ \$SLURM_ARRAY_TASK_ID -eq 23 ]]
then
  chr=X

  mkdir -p $workDir/results/chr\${chr}/

  for anc in AFR CSA
  do

  ./plink/plink2 \
    --pfile \$dataDir/cleaned\${anc}Data/plink2Bin/\${anc}_minMaf0.005_minInfo0.8_chr\${chr} \
    --pheno $workDir/pheno/FPCphenotypes_doubleIDs_\${anc}.txt \
    --covar  $workDir/pheno/covariates_doubleIDs_\${anc}.txt \
    --not-covar sex \
    --pheno-name fpc1-fpc6 \
    --vif 500  \
    --xchr-model 2 \
    --covar-variance-standardize \
    --glm hide-covar cols=+a1count,+a1freq \
    --threads 2 \
    --out $workDir/results/chr\${chr}/chr\${chr}_\${anc}

  done

rsync -av $workDir/results/chr\${chr}/\${chr}_* /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr\${chr}/

else
  chr=\$SLURM_ARRAY_TASK_ID

  for anc in AFR CSA
do


mkdir -p $workDir/results/chr\${chr}/

./plink/plink2 \
  --pfile \$dataDir/cleaned\${anc}Data/plink2Bin/\${anc}_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/FPCphenotypes_doubleIDs_\${anc}.txt \
  --covar  $workDir/pheno/covariates_doubleIDs_\${anc}.txt \
   --pheno-name fpc1-fpc6 \
   --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=+a1count,+a1freq \
  --threads 2 \
  --out $workDir/results/chr\${chr}/chr\${chr}_\${anc}


done

rsync -av $workDir/results/chr\${chr}/chr\${chr}_* /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr\${chr}/

fi


EOF

sbatch $workDir/scripts/plinkAssoc.sh

for i in {1..22} X
do
rsync -av /vast/scratch/users/jackson.v/retThickness/fpcGWASnoExclusions/nonEURGWAS/results/chr${i}/chr${i}* /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr${i}/
done

## format results : 06_formattingNonEurResults.R
## plotting result comparisons : GWASfollowUp/ancestryComparisons.R








