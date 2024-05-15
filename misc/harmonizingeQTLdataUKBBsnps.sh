## harmonisation of Ratnapriya et al. data with UKBB data used in RT GWAS

## set up directories
cd /vast/scratch/users/jackson.v/
mkdir -p harmonisingData
mkdir -p harmonisingData/inVcf
mkdir -p harmonisingData/outVcf
mkdir -p harmonisingData/tmp
mkdir -p harmonisingData/scripts
mkdir -p harmonisingData/logs

cd /vast/scratch/users/jackson.v/harmonisingData

for chr in {1..22}
do

## copy original vcfs
rsync -av /stornext/Bioinf/data/lab_bahlo/public_datasets/RetinaData/Ratnapriya_Data/Genotypes_2021/reimputed_Ratnapriya_et_al_Genotypedata/chr${chr}.dose.vcf.gz ./inVcf

done

# UKK data already on scratch
# in: /vast/scratch/users/jackson.v/retThickness/GWAS/geneticData/cleanedEURData/plinkBin/

# copy R script 
rsync -av /stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/misc/harmonizingeQTLdataUKBBsnps.R ./scripts

## run R script to harmonise
cat <<- EOF > ./scripts/harmoniseVCFs.sh
#!/bin/bash

#SBATCH -J conditional
#SBATCH -o /vast/scratch/users/jackson.v/harmonisingData/logs/harmoniseVcfs_%A_%a.log
#SBATCH -t 4:0:0
#SBATCH --mem=60G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-22


chr=\$SLURM_ARRAY_TASK_ID

module load R/4.3.1

cd /vast/scratch/users/jackson.v/harmonisingData/scripts/

Rscript harmonizingeQTLdataUKBBsnps.R --chr \$chr

EOF

sbatch ./scripts/harmoniseVCFs.sh

## copy to stornext
rsync -av ./outVcf/* /stornext/Bioinf/data/lab_bahlo/public_datasets/RetinaData/Ratnapriya_Data/Genotypes_2021/harmonized_ukb/
