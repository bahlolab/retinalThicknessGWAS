

# wget "https://www.cog-genomics.org/static/bin/plink2_src_220603.zip"

workDir=/vast/scratch/users/jackson.v/retThickness/fpcGWASnoExclusions
dataDir=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedEURData

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
mkdir -p $workDir/clumpedResults
mkdir -p $workDir/temp

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/FPCphenotypes* $workDir/pheno

cd $workDir

mkdir -p $workDir/results/chr22
mkdir -p $workDir/clumpedResults/chr22


 ## test chr 22
 ./plink/plink2 \
   --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr22 \
   --pheno $workDir/pheno/FPCphenotypes_doubleIDs_EUR.txt \
   --pheno-name fpc1-fpc10 \
   --covar  $workDir/pheno/covariates_doubleIDs_EUR.txt \
   --vif 500  \
   --covar-variance-standardize \
   --glm hide-covar cols=+a1count,+a1freq \
   --threads 2 \
   --out results/chr22/chr22EUR

 module load plink

 for i in {1..10}
   do

   plink \
       --bfile $dataDir/plinkBin/EUR_minMaf0.005_minInfo0.8_chr22 \
       --clump $workDir/results/chr22/chr22EUR.fpc${i}.glm.linear \
       --clump-snp-field ID \
    --clump-p1 5e-8 \
    --clump-r2 0.001 \
    --clump-p2 5e-5 \
    --clump-kb 5000 \
    --clump-allow-overlap \
    --threads 2 \
       --out $workDir/clumpedResults/chr22/chr22EUR.${i}

   done


cat <<- EOF > $workDir/scripts/plinkAssoc_allChr.sh
#!/bin/sh

#SBATCH -J plinkAssoc
#SBATCH -o $workDir/logs/plink_%A_%a.log
#SBATCH -t 6:0:0
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-22

chr=\$SLURM_ARRAY_TASK_ID


cd $workDir

mkdir -p $workDir/results/chr\${chr}/
mkdir -p $workDir/clumpedResults/chr\${chr}/

./plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/FPCphenotypes_doubleIDs_EUR.txt \
  --pheno-name fpc1-fpc10 \
  --covar  $workDir/pheno/covariates_doubleIDs_EUR.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=+a1count,+a1freq \
  --threads 2 \
  --out $workDir/results/chr\${chr}/chr\${chr}EUR
  

module load plink


for i in {1..10}
do

plink \
    --bfile $dataDir/plinkBin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
    --clump  $workDir/results/chr\${chr}/chr\${chr}EUR.fpc\${i}.glm.linear \
    --clump-snp-field ID \
    --clump-p1 5e-8 \
    --clump-r2 0.001 \
    --clump-p2 5e-5 \
    --clump-kb 5000 \
    --threads 2 \
    --out $workDir/clumpedResults/chr\${chr}/chr\${chr}EUR.\${i}

done

rsync -av $workDir/results/chr\${chr}/* /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr\${chr}/
rsync -av $workDir/clumpedResults/chr\${chr}/* /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/clumpedResults/chr\${chr}/

EOF

sbatch $workDir/scripts/plinkAssoc_allChr.sh




cat <<- EOF > $workDir/scripts/plinkClump_allChr.sh
#!/bin/sh

#SBATCH -J plinkAssoc
#SBATCH -o $workDir/logs/clump_%A_%a.log
#SBATCH -t 6:0:0
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-22

chr=\$SLURM_ARRAY_TASK_ID


cd $workDir

mkdir -p $workDir/clumpedResults/chr\${chr}/

module load plink

for i in {1..10}
do

plink \
    --bfile $dataDir/plinkBin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
    --clump  /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr\${chr}/chr\${chr}EUR.fpc\${i}.glm.linear \
    --clump-snp-field ID \
    --clump-p1 5e-8 \
    --clump-r2 0.001 \
    --clump-p2 5e-5 \
    --clump-kb 5000 \
    --clump-allow-overlap \
    --threads 2 \
    --out $workDir/clumpedResults/chr\${chr}/chr\${chr}EUR.\${i}_thresh0.001_withOverlap

done

rsync -av $workDir/clumpedResults/chr\${chr}/* /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/clumpedResults/chr\${chr}/

EOF

sbatch $workDir/scripts/plinkClump_allChr.sh






#  module load plink

# for chr in {1..22}
# do
#  for i in {1..10}
#    do

#    plink \
#        --bfile $dataDir/plinkBin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
#        --clump /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr${chr}/chr${chr}EUR.fpc${i}.glm.linear \
#        --clump-snp-field ID \
#     --clump-p1 5e-8 \
#     --clump-r2 0.001 \
#     --clump-p2 5e-5 \
#     --clump-kb 5000 \
#     --clump-allow-overlap \
#     --threads 2 \
#        --out $workDir/clumpedResults/chr${chr}/chr${chr}EUR.${i}_thresh0.001_withOverlap

#    done
# rsync -av $workDir/clumpedResults/chr${chr}/* /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/clumpedResults/chr${chr}/

# done



