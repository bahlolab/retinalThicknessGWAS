

wget "https://www.cog-genomics.org/static/bin/plink2_src_220603.zip"

workDir=/vast/scratch/users/jackson.v/retThickness/GWAS
dataDir=/vast/scratch/users/jackson.v/retThickness/filteringGeneticData/cleanedWhiteBritUnrelatedData


# mkdir $workDir/plink
# cd $workDir/plink
# wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20220603.zip"
# unzip plink2_linux_x86_64_20220603.zip

mkdir $workDir/pheno
mkdir $workDir/results
mkdir $workDir/logs
mkdir $workDir/scripts

cd $workDir

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/phenotypes* $workDir/pheno





./plink/plink2 \
  --bgen $dataDir/imputed/ukbb_minMaf0.001_minInfo0.5_chr22.bgen ref-first \
  --sample $dataDir/imputed/ukbb_minMaf0.001_minInfo0.5_chr22.sample \
  --pheno $workDir/pheno/phenotypesSlice64_doubleIDs.txt \
  --covar  $workDir/pheno/covariates_doubleIDs.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm \
  --out slice64Test


chr=5
  cat <<- EOF > $workDir/scripts/plinkAssoc_chr${chr}.sh
#!/bin/sh

#SBATCH -J plink-${chr}
#SBATCH -o $workDir/logs/plinkChr${chr}_%A_%a.log
#SBATCH -t 4:0:0
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119


slice=\$SLURM_ARRAY_TASK_ID

mkdir -p $workDir/results/\$slice

cd $workDir

./plink/plink2 \
  --bgen $dataDir/imputed/ukbb_minMaf0.001_minInfo0.5_chr${chr}.bgen ref-first \
  --sample $dataDir/imputed/ukbb_minMaf0.001_minInfo0.5_chr${chr}.sample \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs.txt \
  --covar  $workDir/pheno/covariates_doubleIDs.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm \
  --hide-covar \
  --out $workDir/results/\$slice/pixel

EOF

sbatch  $workDir/scripts/plinkAssoc_chr${chr}.sh

sbatch -a 119 $workDir/scripts/plinkAssoc_chr${chr}.sh

sleep 5
done
