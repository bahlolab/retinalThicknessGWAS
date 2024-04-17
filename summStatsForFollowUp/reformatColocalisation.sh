mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/colocInFiles/
mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/colocInFiles/summStats/

workDir=/vast/scratch/users/jackson.v/retThickness/GWAS

mkdir -p $workDir
mkdir -p $workDir/scripts
mkdir -p $workDir/results
mkdir -p $workDir/colocInFiles
mkdir -p $workDir/logs


cd $workDir


## copy scripts and other files
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt  $workDir
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/summStatsForFollowUp/*.R  $workDir/scripts/
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


## run makeRefFie.R
module load R/4.1.3

$workDir/scripts/makeRefFile.R

  cat <<- EOF > $workDir/scripts/colocFormat_gwSig.sh
#!/bin/bash

#SBATCH -J colocGWsig
#SBATCH -o $workDir/logs/colocGWsig_%A_%a.log
#SBATCH -t 2:0:0
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au


end=\$(wc -l loci.txt | cut -d " " -f 1)

module load R/4.1.3
cd $workDir

for locus in \$(seq 1 \$end)
do

    echo 'Reformatting locus '\$locus''
    read locusID chr SNP pixel slice  < <(sed -n \${locus}p loci.txt)

    ## run R script to reformat data
    $workDir/scripts/reformatColocalisationGWsig.R \
    --chr \$chr \
    --slice \$slice \
    --pixel \$pixel


done

EOF

sbatch $workDir/scripts/colocFormat_gwSig.sh

rsync -av $workDir/colocInFiles/*_colocFormat.txt /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/colocInFiles/summStats/
rsync -av $workDir/loci.txt /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/colocInFiles/
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/pixelWise_eQTLs.csv /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/colocInFiles/



