

mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output

workDir=/vast/scratch/users/jackson.v/retThickness/GWAS

mkdir -p $workDir
mkdir -p $workDir/geneticData
mkdir -p $workDir/scripts
mkdir -p $workDir/results
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
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt  $workDir
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/cojoAnalysis/reformatCojo.R  $workDir/scripts/
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/cojoAnalysis/reformatCojoGWsig.R  $workDir/scripts/
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/finalResultsEUR/chr*SNPinfo.txt $workDir


## copy pixelwise results
cd $workDir/
cat <<- EOF > $workDir/scripts/copyData.sh
#!/bin/sh

#SBATCH -J copy
#SBATCH -o $workDir/logs/copyData_%A_%a.log
#SBATCH -t 12:0:0
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-23

if [ \$SLURM_ARRAY_TASK_ID -eq 23 ]; then
    chr="X"
else
    chr=\$SLURM_ARRAY_TASK_ID
fi
cd $workDir/results
rsync -av /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/pixelWiseResultsDec2022/results/chr\${chr}.tar.gz .
tar -xzf chr\${chr}.tar.gz

EOF

sbatch $workDir/scripts/copyData.sh



## full pixelwise cojo takes too long...
## reformat pixelwise results
## create cojo input files
## run cojo

# for chr in {1..22} X
# do

# mkdir -p  $workDir/cojoInFiles/chr${chr}
# mkdir -p  $workDir/cojoOutFiles/chr${chr}

#   cat <<- EOF > $workDir/scripts/cojo_chr${chr}.sh
# #!/bin/bash

# #SBATCH -J cojo-${chr}
# #SBATCH -o $workDir/logs/cojoChr${chr}_%A_%a.log
# #SBATCH -t 48:0:0
# #SBATCH --mem=6G
# #SBATCH --mail-type=FAIL,END
# #SBATCH --mail-user=jackson.v@wehi.edu.au
# #SBATCH -a 1-119%20


# slice=\$SLURM_ARRAY_TASK_ID

# module load R/4.1.3
# cd $workDir
# mkdir -p $workDir/cojoInFiles/chr${chr}/\${slice}
# mkdir -p $workDir/cojoOutFiles/chr${chr}/\${slice}

# ## run R script to reformat data
# $workDir/scripts/reformatCojo.R \
# --chr $chr \
# --slice \$slice

# ## run cojo
# nPix=\$(wc -l pixels.txt | cut -d " " -f 1)

# for pix in \$(seq 1 \$nPix)
# do

# read pixel y x < <(sed -n \${pix}p pixels.txt)

# if [ \$y -eq \$slice ]
# then
#     echo 'Running cojo for Pixel '\$pixel''


#     $workDir/gcta-1.94.1-linux-kernel-3-x86_64/gcta64 \
#         --bfile $workDir/geneticData/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
#         --chr $chr \
#         --cojo-file $workDir/cojoInFiles/chr${chr}/\${slice}/chr${chr}Pixel.\${pixel}_cojoFormat.txt \
#         --cojo-slct \
#         --cojo-p 1.72e-12 \
#         --out $workDir/cojoOutFiles/chr${chr}/\${slice}/chr${chr}Pixel.\${pixel}_cojoOut

#   fi
# done

# EOF

# sbatch $workDir/scripts/cojo_chr${chr}.sh
    
# done

# ## run with more memory and time...
# for chr in {1..8} 
# do

# sbatch --mem=24G --time=48:0:0 $workDir/scripts/cojo_chr${chr}.sh

# done


## run cojo for genome-wide sig only
## 224 loci, split into array jobs of 5

mkdir -p cojoGWsigInFiles
mkdir -p cojoGWsigOutFiles

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/cojoAnalysis/reformatCojoGWsig.R  $workDir/scripts/


  cat <<- EOF > $workDir/scripts/cojo_gwSig.sh
#!/bin/bash

#SBATCH -J cojoGWsig
#SBATCH -o $workDir/logs/cojoGWsig_%A_%a.log
#SBATCH -t 6:0:0
#SBATCH --mem=30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-45


start=\$((\$SLURM_ARRAY_TASK_ID * 5 - 4))

if [ \$SLURM_ARRAY_TASK_ID -eq 45 ]; then
    end=\$(wc -l loci.txt | cut -d " " -f 1)
else
end=\$((\$SLURM_ARRAY_TASK_ID * 5))
fi

echo 'Running cojo for loci '\$start' to '\$end''

module load R/4.1.3
cd $workDir

for locus in \$(seq \$start \$end)
do

    echo 'Reformatting locus '\$locus''
    read locusID chr SNP pixel slice  < <(sed -n \${locus}p loci.txt)


    ## run R script to reformat data
    $workDir/scripts/reformatCojoGWsig.R \
    --chr \$chr \
    --slice \$slice \
    --pixel \$pixel

    echo 'Running cojo for locus '\$locusID' on chr '\$chr' at pixel '\$pixel''

  $workDir/gcta-1.94.1-linux-kernel-3-x86_64/gcta64 \
    --bfile $workDir/geneticData/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
    --chr \$chr \
    --cojo-file $workDir/cojoGWsigInFiles/chr\${chr}Pixel.\${pixel}_cojoFormat.txt \
    --cojo-slct \
    --cojo-p 1.72e-12 \
    --out $workDir/cojoGWsigOutFiles/chr\${chr}Pixel.\${pixel}_cojoOut

done

EOF

sbatch $workDir/scripts/cojo_gwSig.sh


 rsync -av $workDir/cojoGWsigOutFiles/chr*Pixel.*_cojoOut.jma.cojo /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output



# for locus in $(seq 221 224); do

#     read locusID chr SNP pixel slice  < <(sed -n ${locus}p loci.txt)

#     ## run R script to reformat data
#     $workDir/scripts/reformatCojoGWsig.R \
#     --chr $chr \
#     --slice $slice \
#     --pixel $pixel


#   $workDir/gcta-1.94.1-linux-kernel-3-x86_64/gcta64 \
#     --bfile $workDir/geneticData/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chr${chr} \
#     --chr $chr \
#     --cojo-file $workDir/cojoGWsigInFiles/chr${chr}Pixel.${pixel}_cojoFormat.txt \
#     --cojo-slct \
#     --cojo-p 1.72e-12 \
#     --out $workDir/cojoGWsigOutFiles/chr${chr}Pixel.${pixel}_cojoOut

# done