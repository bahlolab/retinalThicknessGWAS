
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

## copy R script
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWASfollowUp/reformatCrossAncestryCorrelationsFPCs.R $workDir/scripts/
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

## update python packages
pip install --upgrade numpy==1.14.2
pip install --upgrade scipy==1.0.1
pip install --upgrade pandas==0.22.0
pip install --upgrade pysnptools==0.3.9
pip install --upgrade bottleneck==1.0.0
pip install --upgrade statsmodels==0.8.0


## fix bug
sed -i -e 's/np.bool/np.bool_/' popcorn/jackknife.py

# ## generate popcorn scores
#    cat <<- EOF > $workDir/scripts/popcornScores.sh
# #!/bin/bash
# #SBATCH -J popcornScores
# #SBATCH -o $workDir/logs/popcornScores_%A_%a.log
# #SBATCH -t 24:0:0
# #SBATCH --mem=48G
# #SBATCH --mail-type=FAIL,END
# #SBATCH --mail-user=jackson.v@wehi.edu.au
# #SBATCH -a 1-23

# cd $workDir

# if [[ \$SLURM_ARRAY_TASK_ID -eq 23 ]]
# then
#   chr=X
# else
#     chr=\$SLURM_ARRAY_TASK_ID
# fi


# module load python

# cd $workDir/Popcorn

# python -m popcorn compute -v 1 --gen_effect \
#     --maf 0.01  --SNPs_to_store 20000 \
#     --bfile1 $workDir/geneticData/cleanedAFRData/plinkBin/AFR_minMaf0.005_minInfo0.8_chr\${chr} \
#     --bfile2 $workDir/geneticData/cleanedCSAData/plinkBin/CSA_minMaf0.005_minInfo0.8_chr\${chr}  \
#     $workDir/popcornOutFiles/chr\${chr}_maf0.01_AFR_CSA_scores.txt

# python -m popcorn compute -v 1 --gen_effect \
#     --maf 0.01  --SNPs_to_store 20000 \
#     --bfile1 $workDir/geneticData/cleanedAFRData/plinkBin/AFR_minMaf0.005_minInfo0.8_chr\${chr} \
#     --bfile2 $workDir/geneticData/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chr\${chr}  \
#     $workDir/popcornOutFiles/chr\${chr}_maf0.01_AFR_EUR_scores.txt

# python -m popcorn compute -v 1 --gen_effect \
#     --maf 0.01  --SNPs_to_store 20000 \
#     --bfile1 $workDir/geneticData/cleanedEURData/plinkBin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
#     --bfile2 $workDir/geneticData/cleanedCSAData/plinkBin/CSA_minMaf0.005_minInfo0.8_chr\${chr}  \
#     $workDir/popcornOutFiles/chr\${chr}_maf0.01_EUR_CSA_scores.txt

# EOF

# sbatch $workDir/scripts/popcornScores.sh 

## cocatenate popcorn scores
# cat $workDir/popcornOutFiles/chr*_AFR_CSA_scores.txt > $workDir/popcornOutFiles/AFR_CSA_scores.txt
# cat $workDir/popcornOutFiles/chr*_AFR_EUR_scores.txt > $workDir/popcornOutFiles/AFR_EUR_scores.txt
# cat $workDir/popcornOutFiles/chr*_EUR_CSA_scores.txt > $workDir/popcornOutFiles/EUR_CSA_scores.txt

for comp in AFR_CSA AFR_EUR EUR_CSA
do

anc1=$(echo $comp | cut -d'_' -f1)
anc2=$(echo $comp | cut -d'_' -f2)

    cat <<- EOF > $workDir/scripts/popcornScores_${comp}.sh
#!/bin/bash
#SBATCH -J popcornScores
#SBATCH -o $workDir/logs/popcornScores_${comp}_%A_%a.log
#SBATCH -t 24:0:0
#SBATCH --mem=48G
#SBATCH --mail-type=FAIL,END
#SBATCH -a 1-23

cd $workDir

if [[ \$SLURM_ARRAY_TASK_ID -eq 23 ]]
then
  chr=X
else
    chr=\$SLURM_ARRAY_TASK_ID
fi


cd $workDir

module load python

cd $workDir/Popcorn

python -m popcorn compute -v 1 --gen_effect \
    --maf 0.01  --SNPs_to_store 20000 \
    --bfile1 $workDir/geneticData/cleaned${anc1}Data/plinkBin/${anc1}_minMaf0.005_minInfo0.8_chr\${chr} \
    --bfile2 $workDir/geneticData/cleaned${anc2}Data/plinkBin/${anc2}_minMaf0.005_minInfo0.8_chr\${chr}  \
    $workDir/popcornOutFiles/chr\${chr}_maf0.01_${comp}_scores.txt


EOF

sbatch $workDir/scripts/popcornScores_${comp}.sh

done


## cocatenate popcorn scores
cat $workDir/popcornOutFiles/chr*_maf0.01_AFR_CSA_scores.txt > $workDir/popcornOutFiles/AFR_CSA_maf0.01_scores.txt
cat $workDir/popcornOutFiles/chr*_maf0.01_AFR_EUR_scores.txt > $workDir/popcornOutFiles/AFR_EUR_maf0.01_scores.txt
cat $workDir/popcornOutFiles/chr*_maf0.01_EUR_CSA_scores.txt > $workDir/popcornOutFiles/EUR_CSA_maf0.01_scores.txt

## use MAF > 0.05 and exclude MHC
## first exclude MHC from chr 6
for anc in AFR CSA EUR
do

    ./plink/plink2 --bfile $workDir/geneticData/cleaned${anc}Data/plinkBin/${anc}_minMaf0.005_minInfo0.8_chr6 \
    --chr 6 \
    --from-bp 28477797 \
    --to-bp 33448354 \
    --make-just-bim \
    --threads 2 \
    --out $workDir/geneticData/cleaned${anc}Data/plinkBin/${anc}_MHCvars

    ./plink/plink2 --bfile $workDir/geneticData/cleaned${anc}Data/plinkBin/${anc}_minMaf0.005_minInfo0.8_chr6 \
    --exclude $workDir/geneticData/cleaned${anc}Data/plinkBin/${anc}_MHCvars.bim \
    --make-bed \
    --threads 2 \
    --out $workDir/geneticData/cleaned${anc}Data/plinkBin/${anc}_minMaf0.005_minInfo0.8_chr6_noMHC

done


for comp in AFR_CSA AFR_EUR EUR_CSA
do

anc1=$(echo $comp | cut -d'_' -f1)
anc2=$(echo $comp | cut -d'_' -f2)

    cat <<- EOF > $workDir/scripts/popcornScores_maf0.05_noMHC_${comp}.sh
#!/bin/bash
#SBATCH -J popcornScores_maf0.05_noMHC
#SBATCH -o $workDir/logs/popcornScores_maf0.05_noMHC_${comp}_%A_%a.log
#SBATCH -t 24:0:0
#SBATCH --mem=48G
#SBATCH --mail-type=FAIL,END
#SBATCH -a 1-23

module load python

cd $workDir/Popcorn

if [[ \$SLURM_ARRAY_TASK_ID -eq 6 ]]
then

python -m popcorn compute -v 1 --gen_effect \
    --maf 0.05  --SNPs_to_store 20000 \
    --bfile1 $workDir/geneticData/cleaned${anc1}Data/plinkBin/${anc1}_minMaf0.005_minInfo0.8_chr6_noMHC \
    --bfile2 $workDir/geneticData/cleaned${anc2}Data/plinkBin/${anc2}_minMaf0.005_minInfo0.8_chr6_noMHC  \
    $workDir/popcornOutFiles/chr6_maf0.05_noMHC_${comp}_scores.txt



else

if [[ \$SLURM_ARRAY_TASK_ID -eq 23 ]]
then
  chr=X
else
    chr=\$SLURM_ARRAY_TASK_ID
fi

python -m popcorn compute -v 1 --gen_effect \
    --maf 0.05  --SNPs_to_store 20000 \
    --bfile1 $workDir/geneticData/cleaned${anc1}Data/plinkBin/${anc1}_minMaf0.005_minInfo0.8_chr\${chr} \
    --bfile2 $workDir/geneticData/cleaned${anc2}Data/plinkBin/${anc2}_minMaf0.005_minInfo0.8_chr\${chr}  \
    $workDir/popcornOutFiles/chr\${chr}_maf0.05_noMHC_${comp}_scores.txt



fi


EOF

sbatch $workDir/scripts/popcornScores_maf0.05_noMHC_${comp}.sh

done

## concatenate scores
cat $workDir/popcornOutFiles/chr*_maf0.05_noMHC_AFR_CSA_scores.txt > $workDir/popcornOutFiles/AFR_CSA_maf0.05_noMHC_scores.txt
cat $workDir/popcornOutFiles/chr*_maf0.05_noMHC_AFR_EUR_scores.txt > $workDir/popcornOutFiles/AFR_EUR_maf0.05_noMHC_scores.txt
cat $workDir/popcornOutFiles/chr*_maf0.05_noMHC_EUR_CSA_scores.txt > $workDir/popcornOutFiles/EUR_CSA_maf0.05_noMHC_scores.txt


## format GWAS results

## run correlation analyses
for comp in AFR_CSA AFR_EUR EUR_CSA
do

anc1=$(echo $comp | cut -d'_' -f1)
anc2=$(echo $comp | cut -d'_' -f2)

    cat <<- EOF > $workDir/scripts/correlations_maf0.01_${comp}.sh
#!/bin/bash
#SBATCH -J correlations
#SBATCH -o $workDir/logs/maf0.01Correlations_${comp}_%A_%a.log
#SBATCH -t 24:0:0
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL,END
#SBATCH -a 1-6

FPC=\$SLURM_ARRAY_TASK_ID

cd $workDir

module load python

cd $workDir/Popcorn

python -m popcorn fit -v 1 --gen_effect \
    --maf 0.01  \
    --sfile1 $workDir/popcornInFiles/allChr.fpc\${FPC}_${anc1}popcornFormat.txt \
    --sfile2 $workDir/popcornInFiles/allChr.fpc\${FPC}_${anc2}popcornFormat.txt  \
    --cfile $workDir/popcornOutFiles/${comp}_maf0.01_scores.txt \
    $workDir/popcornOutFiles/${comp}_FPC\${FPC}_correlations.txt

EOF

sbatch $workDir/scripts/correlations_maf0.01_${comp}.sh

done


## run correlation analyses
for comp in AFR_CSA AFR_EUR EUR_CSA
do

anc1=$(echo $comp | cut -d'_' -f1)
anc2=$(echo $comp | cut -d'_' -f2)

    cat <<- EOF > $workDir/scripts/correlations_${comp}.sh
#!/bin/bash
#SBATCH -J correlations
#SBATCH -o $workDir/logs/correlations_${comp}_%A_%a.log
#SBATCH -t 24:0:0
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-6

FPC=\$SLURM_ARRAY_TASK_ID

cd $workDir

module load python

cd $workDir/Popcorn

python -m popcorn fit -v 1 --gen_effect \
    --maf 0.05  \
    --sfile1 $workDir/popcornInFiles/allChr.fpc\${FPC}_${anc1}popcornFormat.txt \
    --sfile2 $workDir/popcornInFiles/allChr.fpc\${FPC}_${anc2}popcornFormat.txt  \
    --cfile $workDir/popcornOutFiles/${comp}_maf0.05_noMHC_scores.txt \
    $workDir/popcornOutFiles/${comp}_FPC\${FPC}_correlations.txt

EOF

sbatch $workDir/scripts/correlations_${comp}.sh

done





## final attempt: genetic impact instead of genetic effect...
for comp in AFR_CSA AFR_EUR EUR_CSA
do

anc1=$(echo $comp | cut -d'_' -f1)
anc2=$(echo $comp | cut -d'_' -f2)

    cat <<- EOF > $workDir/scripts/popcornScores_maf0.05_noMHC_geneticImpact_${comp}.sh
#!/bin/bash
#SBATCH -J popcornScores_maf0.05_noMHC_geneticImpact
#SBATCH -o $workDir/logs/popcornScores_maf0.05_noMHC_geneticImpact_${comp}_%A_%a.log
#SBATCH -t 24:0:0
#SBATCH --mem=48G
#SBATCH --mail-type=FAIL,END
#SBATCH -a 1-23

module load python

cd $workDir/Popcorn

if [[ \$SLURM_ARRAY_TASK_ID -eq 6 ]]
then

python -m popcorn compute -v 1  \
    --maf 0.05  --SNPs_to_store 10000 \
    --bfile1 $workDir/geneticData/cleaned${anc1}Data/plinkBin/${anc1}_minMaf0.005_minInfo0.8_chr6_noMHC \
    --bfile2 $workDir/geneticData/cleaned${anc2}Data/plinkBin/${anc2}_minMaf0.005_minInfo0.8_chr6_noMHC  \
    $workDir/popcornOutFiles/chr6_maf0.05_noMHC_geneticImpact_${comp}_scores.txt



else

if [[ \$SLURM_ARRAY_TASK_ID -eq 23 ]]
then
  chr=X
else
    chr=\$SLURM_ARRAY_TASK_ID
fi

python -m popcorn compute -v 1  \
    --maf 0.05  --SNPs_to_store 10000 \
    --bfile1 $workDir/geneticData/cleaned${anc1}Data/plinkBin/${anc1}_minMaf0.005_minInfo0.8_chr\${chr} \
    --bfile2 $workDir/geneticData/cleaned${anc2}Data/plinkBin/${anc2}_minMaf0.005_minInfo0.8_chr\${chr}  \
    $workDir/popcornOutFiles/chr\${chr}_maf0.05_noMHC_geneticImpact_${comp}_scores.txt



fi


EOF

sbatch $workDir/scripts/popcornScores_maf0.05_noMHC_geneticImpact_${comp}.sh

done

## concatenate scores
cat $workDir/popcornOutFiles/chr*_maf0.05_noMHC_geneticImpact_AFR_CSA_scores.txt > $workDir/popcornOutFiles/AFR_CSA_maf0.05_noMHC_geneticImpact_scores.txt
cat $workDir/popcornOutFiles/chr*_maf0.05_noMHC_geneticImpact_AFR_EUR_scores.txt > $workDir/popcornOutFiles/AFR_EUR_maf0.05_noMHC_geneticImpact_scores.txt
cat $workDir/popcornOutFiles/chr*_maf0.05_noMHC_geneticImpact_EUR_CSA_scores.txt > $workDir/popcornOutFiles/EUR_CSA_maf0.05_noMHC_geneticImpact_scores.txt


## run correlation analyses
for comp in AFR_CSA AFR_EUR EUR_CSA
do

anc1=$(echo $comp | cut -d'_' -f1)
anc2=$(echo $comp | cut -d'_' -f2)

    cat <<- EOF > $workDir/scripts/correlations_geneticImpact_${comp}.sh
#!/bin/bash
#SBATCH -J correlations_geneticImpact
#SBATCH -o $workDir/logs/correlations_geneticImpact_${comp}_%A_%a.log
#SBATCH -t 24:0:0
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-6

FPC=\$SLURM_ARRAY_TASK_ID

cd $workDir

module load python

cd $workDir/Popcorn

python -m popcorn fit -v 1  \
    --maf 0.05  \
    --sfile1 $workDir/popcornInFiles/allChr.fpc\${FPC}_${anc1}popcornFormat.txt \
    --sfile2 $workDir/popcornInFiles/allChr.fpc\${FPC}_${anc2}popcornFormat.txt  \
    --cfile $workDir/popcornOutFiles/${comp}_maf0.05_noMHC_geneticImpact_scores.txt \
    $workDir/popcornOutFiles/${comp}_FPC\${FPC}_geneticImpact_correlations.txt

EOF

sbatch $workDir/scripts/correlations_geneticImpact_${comp}.sh

done

