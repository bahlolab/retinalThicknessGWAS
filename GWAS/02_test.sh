

wget "https://www.cog-genomics.org/static/bin/plink2_src_220603.zip"

workDir=/vast/scratch/users/jackson.v/retThickness/GWAS
dataDir=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedWhiteBritUnrelatedData
dataDir=/vast/scratch/users/jackson.v/retThickness/filteringGeneticData/cleanedWhiteBritUnrelatedData/

cd $workDir

# mkdir $workDir/plink
# cd $workDir/plink
# wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20220603.zip"
# unzip plink2_linux_x86_64_20220603.zip

mkdir $workDir/pheno
mkdir $workDir/results
mkdir $workDir/logs
mkdir $workDir/scripts


rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/phenotypes* $workDir/pheno





./plink/plink2 \
  --bgen $dataDir/imputed/ukbb_minMaf0.001_minInfo0.5_chr22.bgen ref-first \
  --sample $dataDir/imputed/ukbb_minMaf0.001_minInfo0.5_chr22.sample \
  --pheno $workDir/pheno/phenotypesSlice64_doubleIDs.txt \
  --covar  $workDir/pheno/covariates_doubleIDs.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=+a1count,+a1freq \
  --threads 2 \
  --out slice64Test/chr5Pixel


  cat <<- EOF > $workDir/scripts/plinkAssoc_allChr.sh
#!/bin/sh

#SBATCH -J plinkAssoc
#SBATCH -o $workDir/logs/plink_%A_%a.log
#SBATCH -t 48:0:0
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119


slice=\$SLURM_ARRAY_TASK_ID


cd $workDir

for chr in {22..1}
do

mkdir -p $workDir/results/chr\${chr}/\${slice}

./plink/plink2 \
  --bgen $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr\${chr}.bgen ref-first \
  --sample $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr\${chr}.sample \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs.txt \
  --covar  $workDir/pheno/covariates_doubleIDs.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar  cols=+a1count,+a1freq \
  --threads 2 \
  --out $workDir/results/chr\${chr}/\${slice}/chr\${chr}Pixel

done

EOF

# sbatch  $workDir/scripts/plinkAssoc_allChr.sh

sbatch -a 64 $workDir/scripts/plinkAssoc_allChr.sh




chr=5
  cat <<- EOF > $workDir/scripts/plinkAssoc_chr${chr}.sh
#!/bin/sh

#SBATCH -J plink-${chr}
#SBATCH -o $workDir/logs/plinkChr${chr}_%A_%a.log
#SBATCH -t 6:0:0
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119


slice=\$SLURM_ARRAY_TASK_ID


cd $workDir


mkdir -p $workDir/results/chr${chr}/\${slice}

./plink/plink2 \
  --bgen $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr${chr}.bgen ref-first \
  --sample $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr${chr}.sample \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs.txt \
  --covar  $workDir/pheno/covariates_doubleIDs.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=+a1count,+a1freq \
  --threads 2 \
  --out $workDir/results/chr${chr}/\${slice}/chr${chr}Pixel


EOF

sbatch -a 64 $workDir/scripts/plinkAssoc_chr${chr}.sh

sbatch  $workDir/scripts/plinkAssoc_chr${chr}.sh
sbatch -a 1-63 $workDir/scripts/plinkAssoc_chr${chr}.sh
sbatch -a 65-119 $workDir/scripts/plinkAssoc_chr${chr}.sh



# for chr in {6..17}
chr=2
# do
  cat <<- EOF > $workDir/scripts/plinkAssoc_chr${chr}.sh
#!/bin/sh

#SBATCH -J plink-${chr}
#SBATCH -o $workDir/logs/plinkChr${chr}_%A_%a.log
#SBATCH -t 6:0:0
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119


slice=\$SLURM_ARRAY_TASK_ID


cd $workDir


mkdir -p $workDir/results/chr${chr}/\${slice}

./plink/plink2 \
  --bgen $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr${chr}.bgen ref-first \
  --sample $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr${chr}.sample \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs.txt \
  --covar  $workDir/pheno/covariates_doubleIDs.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=+a1count,+a1freq \
  --threads 2 \
  --out $workDir/results/chr${chr}/\${slice}/chr${chr}Pixel


EOF

# sbatch -a 64 $workDir/scripts/plinkAssoc_chr${chr}.sh
# sbatch $workDir/scripts/plinkAssoc_chr${chr}.sh

sleep 1
# done
