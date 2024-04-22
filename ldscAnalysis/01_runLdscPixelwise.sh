

workDir=/vast/scratch/users/jackson.v/retThickness/GWAS

mkdir -p $workDir
mkdir -p $workDir/geneticData
mkdir -p $workDir/scripts
mkdir -p $workDir/results
mkdir -p $workDir/ldscInFiles
mkdir -p $workDir/ldscOutFiles
mkdir -p $workDir/logs

mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldsc/
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldsc/output/
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldsc/output/plots/

cd $workDir
module load anaconda3 
# download ldsc
git clone https://github.com/bulik/ldsc.git
cd ldsc
conda env create --file environment.yml
source activate ldsc

## download ld scores in $workDir/ldscInFiles
# 
cd  $workDir/ldscInFiles

# wget https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.tar.gz $workDir/ldscInFiles
wget https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz $workDir/ldscInFiles
tar -xzf UKBB.ALL.ldscore.tar.gz
# tar -xvf baselineLF_v2.2.UKB.tar.gz

## use PanUKBB LDscore
## copy M_5_50 files, so can run with rsid versions...
cp /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.l2.M /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.M
cp /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.l2.M_5_50 /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.M_5_50
cp /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.l2.M /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.rsid.l2.M
cp /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.l2.M_5_50 /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.rsid.l2.M_5_50

## copy genetic data
# rsync -av /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedEURData.tar.gz $workDir/geneticData
# cd $workDir/geneticData
# tar -xzf cleanedEURData.tar.gz


## copy scripts and other files
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt  $workDir
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/ldscAnalysis/reformatLDSC.R  $workDir/scripts/
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/finalResultsEUR/chr*SNPinfo.txt $workDir


## copy pixelwise results
# cd $workDir/
# cat <<- EOF > $workDir/scripts/copyData.sh
# #!/bin/sh

# #SBATCH -J copy
# #SBATCH -o $workDir/logs/copyData_%A_%a.log
# #SBATCH -t 12:0:0
# #SBATCH --mem=4G
# #SBATCH --mail-type=FAIL,END
# #SBATCH --mail-user=jackson.v@wehi.edu.au
# #SBATCH -a 1-23

# if [ \$SLURM_ARRAY_TASK_ID -eq 23 ]; then
#     chr="X"
# else
#     chr=\$SLURM_ARRAY_TASK_ID
# fi
# cd $workDir/results
# rsync -av /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr\${chr}.tar.gz .
# tar -xzf chr\${chr}.tar.gz

# EOF

# sbatch $workDir/scripts/copyData.sh




## reformat pixelwise results
## create ldsc input files
## run ldsc

mkdir -p  $workDir/ldscInFiles/
mkdir -p  $workDir/ldscOutFiles/

  cat <<- EOF > $workDir/scripts/ldsc_reformat.sh
#!/bin/bash

#SBATCH -J ldscReformat
#SBATCH -o $workDir/logs/ldscReformat_%A_%a.log
#SBATCH -t 6:0:0
#SBATCH --mem=6G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119%60


slice=\$SLURM_ARRAY_TASK_ID

module load R/4.1.3

cd $workDir
mkdir -p $workDir/ldscInFiles/\${slice}
mkdir -p $workDir/ldscOutFiles/\${slice}

## run R script to reformat data
$workDir/scripts/reformatLDSC.R \
--slice \$slice


EOF

sbatch $workDir/scripts/ldsc_reformat.sh
    




  cat <<- EOF > $workDir/scripts/runLDSC.sh
#!/bin/bash

#SBATCH -J ldsc
#SBATCH -o $workDir/logs/ldsc_%A_%a.log
#SBATCH -t 48:0:0
#SBATCH --mem=6G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119%60


slice=\$SLURM_ARRAY_TASK_ID

module load anaconda3 

cd $workDir

## run ldsc
nPix=\$(wc -l pixels.txt | cut -d " " -f 1)

cd $workDir/ldsc
conda env create --file environment.yml
source activate ldsc

cd $workDir/
for pix in \$(seq 1 \$nPix)
do

read pixel y x < <(sed -n \${pix}p pixels.txt)

if [ \$y -eq \$slice ]
then
    echo 'Running ldsc for Pixel '\$pixel''

    ## munge sumstats
    $workDir/ldsc/munge_sumstats.py \
      --sumstats $workDir/ldscInFiles/\${slice}/\${pixel}_ldscFormat.txt \
      --out $workDir/ldscInFiles/\${slice}/\${pixel} 

    ./ldsc/ldsc.py \
    --h2 $workDir/ldscInFiles/\${slice}/\${pixel}.sumstats.gz \
    --ref-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid \
    --w-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid \
    --out $workDir/ldscOutFiles/\${slice}/\${pixel}.univariate_h2

    ./ldsc/ldsc.py \
    --h2 $workDir/ldscInFiles/\${slice}/\${pixel}.sumstats.gz \
    --ref-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.rsid \
    --w-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid \
    --out $workDir/ldscOutFiles/\${slice}/\${pixel}.stratified_h2

  fi
done

EOF

sbatch $workDir/scripts/runLDSC.sh
    


rsync -av /vast/scratch/users/jackson.v/retThickness/GWAS/ldscOutFiles/* /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/ldsc/output/





# slice=1
# pixel=1_39
# cd $workDir/ldsc
# conda env create --file environment.yml
# source activate ldsc


# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2


# ./munge_sumstats.py \
#       --sumstats $workDir/ldscInFiles/${slice}/${pixel}_ldscFormat.txt \
#       --out $workDir/ldscInFiles/${slice}/${pixel} 

# ./ldsc.py \
# --h2 $workDir/ldscInFiles/${slice}/${pixel}.sumstats.gz \
# --ref-ld-chr $workDir/ldscInFiles/baselineLF_v2.2.UKB/baselineLF2.2.UKB.@ \
# --w-ld-chr $workDir/ldscInFiles/baselineLF_v2.2.UKB/weights.UKB.@ \
# --out $workDir/ldscOutFiles/${pixel}_h2



# ./ldsc/ldsc.py \
# --h2 $workDir/ldscInFiles/${slice}/${pixel}.sumstats.gz \
# --ref-ld-chr $workDir/ldscInFiles/baselineLF_v2.2.UKB/weights.UKB.@ \
# --w-ld-chr $workDir/ldscInFiles/baselineLF_v2.2.UKB/weights.UKB.@ \
# --out $workDir/ldscOutFiles/${pixel}_h2


# cp /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.l2.M /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.M
# cp /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.l2.M_5_50 /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.M_5_50
# cp /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.l2.M /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.rsid.l2.M
# cp /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.l2.M_5_50 /vast/scratch/users/jackson.v/retThickness/GWAS/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.rsid.l2.M_5_50

# ./ldsc/ldsc.py \
# --h2 $workDir/ldscInFiles/${slice}/${pixel}.sumstats.gz \
# --ref-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid \
# --w-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid \
# --out $workDir/ldscOutFiles/${pixel}_h2

# ./ldsc/ldsc.py \
# --h2 $workDir/ldscInFiles/${slice}/${pixel}.sumstats.gz \
# --ref-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.25LDMS.rsid \
# --w-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid \
# --out $workDir/ldscOutFiles/${pixel}.25LDMS_h2

# ./ldsc/ldsc.py \
# --h2 $workDir/ldscInFiles/${slice}/${pixel}.sumstats.gz \
# --ref-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.8LDMS.rsid \
# --w-ld $workDir/ldscInFiles/UKBB.ALL.ldscore/UKBB.EUR.rsid \
# --out $workDir/ldscOutFiles/${pixel}.8LDMS_h2

