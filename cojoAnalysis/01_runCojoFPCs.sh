

mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output

workDir=/vast/scratch/users/jackson.v/retThickness/GWAS

mkdir -p $workDir
mkdir -p $workDir/geneticData
mkdir -p $workDir/scripts
mkdir -p $workDir/fpcResults
mkdir -p $workDir/cojoInFiles
mkdir -p $workDir/cojoOutFiles
mkdir -p $workDir/logs


cd $workDir

## download gcta
# wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip
# unzip gcta-1.94.1-linux-kernel-3-x86_64.zip

## copy genetic data
rsync -av /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleanedEURData.tar.gz $workDir/geneticData
cd $workDir/geneticData
tar -xzf cleanedEURData.tar.gz

## recode X to chr 23
module load plink

## make chr update file
echo -e "X\t23" > $workDir/chrXupdate.txt
plink \
    --bfile $workDir/geneticData/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chrX \
    --update-chr $workDir/chrXupdate.txt \
    --make-bed \
    --out $workDir/geneticData/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chr23


## copy scripts and other files
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/cojoAnalysis/reformatCojoFPCs.R  $workDir/scripts/


## copy fpc results
cd $workDir/fpcResults

for i in {1..6} 
do

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr*/chr*EUR.fpc${i}.glm.linear .
done

## run R script to reformat data
module load R/4.1.3
cd $workDir
$workDir/scripts/reformatCojoFPCs.R 



## create cojo input files
## run cojo


  cat <<- EOF > $workDir/scripts/cojofPCs.sh
#!/bin/bash

#SBATCH -J cojo-FPCs
#SBATCH -o $workDir/logs/cojoFPCs_%A_%a.log
#SBATCH -t 6:0:0
#SBATCH --mem=24G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-6


fpc=\$SLURM_ARRAY_TASK_ID

mkdir -p $workDir/cojoOutFiles/fpc\$fpc

for chr in \$(seq 1 22)
do

    echo 'Running cojo for Chr '\$chr''


    $workDir/gcta-1.94.1-linux-kernel-3-x86_64/gcta64 \
        --bfile $workDir/geneticData/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
        --chr \$chr \
        --cojo-file $workDir/cojoInFiles/chr\${chr}.fpc\${fpc}_cojoFormat.txt \
        --cojo-slct \
        --cojo-p 8.33e-9 \
        --out $workDir/cojoOutFiles/fpc\${fpc}/chr\${chr}.fpc\${fpc}_cojoOut

done

    echo 'Running cojo for Chr X'


    $workDir/gcta-1.94.1-linux-kernel-3-x86_64/gcta64 \
        --bfile $workDir/geneticData/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chr23 \
        --chr 23 \
        --cojo-file $workDir/cojoInFiles/chrX.fpc\${fpc}_cojoFormat.txt \
        --cojo-slct \
        --cojo-p 8.33e-9 \
        --out $workDir/cojoOutFiles/fpc\${fpc}/chrX.fpc\${fpc}_cojoOut



EOF

sbatch $workDir/scripts/cojofPCs.sh
    



 rsync -av $workDir/cojoOutFiles/fpc*/chr*.fpc*_cojoOut.jma.cojo /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output

