
workDir=/vast/scratch/users/jackson.v/retThickness/fpcGWASnoExclusions/nonEURGWAS

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

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/*.txt $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/GWAS/sentinels/chr*sentinelsIDonly_clumpThresh0.001_withOverlap.txt $workDir/pheno


  cat <<- EOF > $workDir/scripts/plinkAssoc.sh
#!/bin/bash
#SBATCH -J plink
#SBATCH -o $workDir/logs/plink_%A.log
#SBATCH -t 0:2:0
#SBATCH --mem=2G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au


cd $workDir

for anc in AFR CSA
do

dataDir=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/cleanedGeneticFiles/cleaned\${anc}Data

for chr in {1..22}
do

mkdir -p $workDir/results/chr\${chr}/

./plink/plink2 \
  --pfile \$dataDir/plink2Bin/\${anc}_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/FPCphenotypes_doubleIDs_\${anc}.txt \
  --covar  $workDir/pheno/covariates_doubleIDs_\${anc}.txt \
   --pheno-name fpc1-fpc6 \
   --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=+a1count,+a1freq \
  --extract $workDir/pheno/chr\${chr}sentinelsIDonly_clumpThresh0.001_withOverlap.txt \
  --threads 2 \
  --out $workDir/results/chr\${chr}/chr\${chr}_\${anc}


done

chr=X

mkdir -p $workDir/results/chr\${chr}/

./plink/plink2 \
  --pfile \$dataDir/plink2Bin/\${anc}_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/FPCphenotypes_doubleIDs_\${anc}.txt \
  --covar  $workDir/pheno/covariates_doubleIDs_\${anc}.txt \
  --not-covar sex \
  --pheno-name fpc1-fpc6 \
   --vif 500  \
  --xchr-model 2 \
  --covar-variance-standardize \
  --glm hide-covar cols=+a1count,+a1freq \
  --extract $workDir/pheno/chr\${chr}sentinelsIDonly_clumpThresh0.001_withOverlap.txt \
  --threads 2 \
  --out $workDir/results/chr\${chr}/chr\${chr}_\${anc}

done

EOF

sbatch $workDir/scripts/plinkAssoc.sh


## format results : 06_formattingNonEurResults.R
## plotting result comparisons : GWASfollowUp/ancestryComparisons.R








