

workDir=/vast/scratch/users/jackson.v/retThickness/GWAS

mkdir -p $workDir
mkdir -p $workDir/geneticData
mkdir -p $workDir/scripts
mkdir -p $workDir/fpcResults
mkdir -p $workDir/ldscInFilesSensitivity
mkdir -p $workDir/ldscOutFilesSensitivity
mkdir -p $workDir/logs

mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldscSensitivity/
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldscSensitivity/output/
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldscSensitivity/output/plots/

cd $workDir
# module load anaconda3 
# # download ldsc
# git clone https://github.com/bulik/ldsc.git
# cd ldsc
# conda env create --file environment.yml
# source activate ldsc

## download ld scores in $workDir/ldscInFiles
# 
# cd  $workDir/ldscInFiles

# wget https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.tar.gz $workDir/ldscInFiles
# wget https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz $workDir/ldscInFiles
# tar -xzf UKBB.ALL.ldscore.tar.gz
# # tar -xvf baselineLF_v2.2.UKB.tar.gz

# ## use PanUKBB LDscore
# ## copy M_5_50 files, so can run with rsid versions...
# cp /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.l2.M /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.M
# cp /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.l2.M_5_50 /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.M_5_50
# cp /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.l2.M /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.rsid.l2.M
# cp /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.l2.M_5_50 /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.rsid.l2.M_5_50

## copy genetic data
# rsync -av /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedEURData.tar.gz $workDir/geneticData
# cd $workDir/geneticData
# tar -xzf cleanedEURData.tar.gz


## copy scripts and other files
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/sensitivityAnalyses/reformatLDSCfPCs.R  $workDir/scripts/
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/sensitivityAnalyses/plotLDSCgenCorrs.R  $workDir/scripts/
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/finalResultsEUR/chr*SNPinfo.txt $workDir


# copy fPC results
cd $workDir/fpcResults

for i in {1..6} 
do

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr*/chr*EUR.fpc${i}.glm.linear .
done


## reformat fpc results
## create ldsc input files
## run ldsc

mkdir -p  $workDir/ldscInFiles/fPCs
mkdir -p  $workDir/ldscOutFiles/fPCs


module load R/4.1.3

cd $workDir
## run R script to reformat data
$workDir/scripts/reformatLDSCfPCs.R 


## run LDSC

cat <<- EOF > $workDir/scripts/runLDSCfPCsGenCorr.sh
#!/bin/bash

#SBATCH -J ldscFPCsGenCoorr
#SBATCH -o $workDir/logs/ldscFPCsGenCorr_%A_%a.log
#SBATCH -t 6:0:0
#SBATCH --mem=6G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH --array=1-6

fpc=\$SLURM_ARRAY_TASK_ID

module load anaconda3 

cd $workDir

## run ldsc

cd $workDir/ldsc
conda env create --file environment.yml
source activate ldsc

cd $workDir/

  echo 'Running ldsc for FPC '\$fpc''
  ## munge sumstats
  $workDir/ldsc/munge_sumstats.py \
    --sumstats $workDir/ldscInFiles/fPCs/fPC\${fpc}_ldscFormat.txt \
    --out $workDir/ldscInFiles/fPCs/fPC\${fpc} 

   ## munge sumstats
  $workDir/ldsc/munge_sumstats.py \
    --sumstats $workDir/ldscInFiles/fPCs/fPC\${fpc}_smoking_ldscFormat.txt \
    --out $workDir/ldscInFiles/fPCs/fPC\${fpc}_smoking 

  ## munge sumstats
  $workDir/ldsc/munge_sumstats.py \
    --sumstats $workDir/ldscInFiles/fPCs/fPC\${fpc}_noSurgery_ldscFormat.txt \
    --out $workDir/ldscInFiles/fPCs/fPC\${fpc}_noSurgery 

## run correlations
echo 'Running correlations for FPC '\$fpc''

./ldsc/ldsc.py \
  --rg $workDir/ldscInFiles/fPCs/fPC\${fpc}.sumstats.gz,$workDir/ldscInFiles/fPCs/fPC\${fpc}_smoking.sumstats.gz \
  --ref-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid \
  --w-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid \
  --out $workDir/ldscOutFiles/fPCs/fPC\${fpc}.smoking

 ./ldsc/ldsc.py \
  --rg $workDir/ldscInFiles/fPCs/fPC\${fpc}.sumstats.gz,$workDir/ldscInFiles/fPCs/fPC\${fpc}_noSurgery.sumstats.gz \
  --ref-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid \
  --w-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid \
  --out $workDir/ldscOutFiles/fPCs/fPC\${fpc}.noSurgery


EOF

sbatch $workDir/scripts/runLDSCfPCsGenCorr.sh
    
## extract and plot results

module load R/4.1.3
$workDir/scripts/plotLDSCgenCorrs.R

rsync -av /vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs/fPC*smoking.log /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldscSensitivity/output/
rsync -av /vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs/smokingGenCorr.txt /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldscSensitivity/output/
rsync -av /vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPC*noSurgery.log /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldscSensitivity/output/
rsync -av /vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/fPCs/noSurgeryGenCorr.txt /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldscSensitivity/output/





