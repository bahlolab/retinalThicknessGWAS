
workDir=/vast/scratch/users/jackson.v/retThickness/GWAS

mkdir -p $workDir
mkdir -p $workDir/geneticData
mkdir -p $workDir/scripts
mkdir -p $workDir/fpcResults
mkdir -p $workDir/results

mkdir -p $workDir/popcornInFiles
mkdir -p $workDir/popcornOutFiles
mkdir -p $workDir/logs


## copy genetic data
rsync -av /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/*.tar.gz $workDir/geneticData
cd $workDir/geneticData

# for file in *.tar.gz
for file in cleanedAFRData.tar.gz  cleanedCSAData.tar.gz
do
    tar -xzf $file
done

# copy GWAS results
cd $workDir/fpcResults

for i in {1..6}
do
 rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/results/chr*/chr*.fpc${i}.glm.linear .
done

## load Popcorn
module load python

cd $workDir
git clone https://github.com/brielin/Popcorn.git
cd Popcorn
pip install .

python -m popcorn compute -v 1 --gen_effect \
    --bfile1 $workDir/geneticData/cleanedAFRData/plinkBin/AFR_minMaf0.005_minInfo0.8_chr22 \
    --bfile2 $workDir/geneticData/cleanedCSAData/plinkBin/CSA_minMaf0.005_minInfo0.8_chr22  \
    $workDir/popcornOutFiles/chr22TestScores_AFR_CSA.txt



   cat <<- EOF > $workDir/scripts/popcornScores.sh
#!/bin/bash
#SBATCH -J popcornScores
#SBATCH -o $workDir/logs/popcornScores_%A_%a.log
#SBATCH -t 24:0:0
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-23

cd $workDir

if [[ \$SLURM_ARRAY_TASK_ID -eq 23 ]]
then
  chr=X
else
    chr=\$SLURM_ARRAY_TASK_ID
fi


module load python

cd $workDir/Popcorn

python -m popcorn compute -v 1 --gen_effect \
    --bfile1 $workDir/geneticData/cleanedAFRData/plinkBin/AFR_minMaf0.005_minInfo0.8_chr\${chr} \
    --bfile2 $workDir/geneticData/cleanedCSAData/plinkBin/CSA_minMaf0.005_minInfo0.8_chr\${chr}  \
    $workDir/popcornOutFiles/chr\${chr}_AFR_CSA_scores.txt

python -m popcorn compute -v 1 --gen_effect \
    --bfile1 $workDir/geneticData/cleanedAFRData/plinkBin/AFR_minMaf0.005_minInfo0.8_chr\${chr} \
    --bfile2 $workDir/geneticData/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chr\${chr}  \
    $workDir/popcornOutFiles/chr\${chr}_AFR_EUR_scores.txt

python -m popcorn compute -v 1 --gen_effect \
    --bfile1 $workDir/geneticData/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
    --bfile2 $workDir/geneticData/cleanedCSAData/plinkBin/CSA_minMaf0.005_minInfo0.8_chr\${chr}  \
    $workDir/popcornOutFiles/chr\${chr}_EUR_CSA_scores.txt

EOF

sbatch $workDir/scripts/popcornScores.sh 

