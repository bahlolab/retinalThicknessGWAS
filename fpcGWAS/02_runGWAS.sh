

# wget "https://www.cog-genomics.org/static/bin/plink2_src_220603.zip"

workDir=/vast/scratch/users/jackson.v/retThickness/fpcGWAS
dataDir=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedWhiteBritUnrelatedData

cd $workDir

mkdir -p $workDir/plink
cd $workDir/plink
wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20220603.zip"
unzip plink2_linux_x86_64_20220603.zip

mkdir -p $workDir/pheno
mkdir -p $workDir/results
mkdir -p $workDir/logs
mkdir -p $workDir/scripts
mkdir -p $workDir/clumpedResults
mkdir -p $workDir/temp

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWAS/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWAS/output/FPCphenotypes* $workDir/pheno

cd $workDir

# ## test chr 22
# ./plink/plink2 \
#   --bgen $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr22.bgen ref-first \
#   --sample $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr22.sample \
#   --pheno $workDir/pheno/FPCphenotypes_doubleIDs.txt \
#   --covar  $workDir/pheno/covariates_doubleIDs.txt \
#   --vif 500  \
#   --covar-variance-standardize \
#   --glm hide-covar cols=+a1count,+a1freq \
#   --threads 2 \
#   --out results/chr22

# /vast/scratch/users/jackson.v/retThickness/GWAS/plink/plink2 \
#   --bgen $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr22.bgen ref-first \
#   --sample $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr22.sample \
#   --pheno $workDir/pheno/FPCphenotypes_doubleIDs.txt \
#   --make-bed \
#   --rm-dup force-first \
#   --threads 2 \
#   --out $workDir/temp/ukbb_minMaf0.001_minInfo0.8_chr22

# mkdir $workDir/clumpedResults/chr22

# module load plink

# for i in {1..10}
#   do

#   plink \
#       --bfile $workDir/temp/ukbb_minMaf0.001_minInfo0.8_chr22 \
#       --clump $workDir/results/chr22.fpc${i}.glm.linear \
#       --clump-snp-field ID \
#       --clump-p1 5e-8 \
#       --clump-r2 0.2 \
#       --clump-p2 5e-5 \
#       --clump-kb 500 \
#       --threads 2 \
#       --out $workDir/clumpedResults/chr22/chr22.${i}

#   done


cat <<- EOF > $workDir/scripts/plinkAssoc_allChr.sh
#!/bin/sh

#SBATCH -J plinkAssoc
#SBATCH -o $workDir/logs/plink_%A_%a.log
#SBATCH -t 16:0:0
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-22

chr=\$SLURM_ARRAY_TASK_ID


cd $workDir

mkdir -p $workDir/results/chr\${chr}/
mkdir -p $workDir/clumpedResults/chr\${chr}/

./plink/plink2 \
  --bgen $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr\${chr}.bgen ref-first \
  --sample $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr\${chr}.sample \
  --mac 500 \
  --pheno $workDir/pheno/FPCphenotypes_doubleIDs.txt  \
  --covar  $workDir/pheno/covariates_doubleIDs.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=+a1count,+a1freq \
  --threads 2 \
  --out $workDir/results/chr\${chr}/chr\${chr}


  ./plink/plink2 \
    --bgen $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr\${chr}.bgen ref-first \
    --sample $dataDir/imputed/ukbb_minMaf0.001_minInfo0.8_chr\${chr}.sample \
    --mac 500 \
    --pheno $workDir/pheno/FPCphenotypes_doubleIDs.txt  \
    --make-bed \
    --rm-dup force-first \
    --threads 2 \
    --out $workDir/temp/ukbb_minMAC500_minInfo0.8_chr\${chr}

  mkdir $workDir/clumpedResults/chr\${chr}

  module load plink


for i in {1..100}
do

plink \
    --bfile $workDir/temp/ukbb_minMAC500_minInfo0.8_chr\${chr} \
    --clump  $workDir/results/chr\${chr}/chr\${chr}.fpc\${i}.glm.linear \
    --mac 500 \
    --clump-snp-field ID \
    --clump-p1 5e-8 \
    --clump-r2 0.2 \
    --clump-p2 5e-5 \
    --clump-kb 500 \
    --threads 2 \
    --out $workDir/clumpedResults/chr\${chr}/chr\${chr}.\${i}

done

EOF

sbatch $workDir/scripts/plinkAssoc_allChr.sh





mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/fpcGWAS/results
mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/fpcGWAS/clumpedResults
rsync -av $workDir/results/* /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/fpcGWAS/results
rsync -av $workDir/clumpedResults/* /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/fpcGWAS/ClumpedResults
