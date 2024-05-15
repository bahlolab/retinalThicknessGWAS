workDir=/vast/scratch/users/jackson.v/retThickness/GWAS
dataDir=$workDir/geneticData/cleanedEURData

cd $workDir

mkdir -p $workDir/plink
cd $workDir/plink
wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20221024.zip"
unzip plink2_linux_x86_64_20221024.zip

mkdir -p $workDir/pheno
mkdir -p $workDir/results
mkdir -p $workDir/logs
mkdir -p $workDir/scripts

  
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWASnoExclusions/output/FPCphenotypes* $workDir/pheno

# Print first 3 columns of the file, but omitting the first line
awk ' NR > 1 { print $1,$2,$3}' /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWASfollowUp/output/annotations/annotatedSentinelsFPCs.txt > $workDir/fPCloci.txt

# set n as the number of lines in $workDir/fPCloci.txt
n=$(wc -l $workDir/fPCloci.txt | awk '{print $1}')

cat <<- EOF > $workDir/scripts/conditionalAnalysesFPCs.sh
#!/bin/bash

#SBATCH -J conditional
#SBATCH -o $workDir/logs/conditionalAnalysesFPCs_%A_%a.log
#SBATCH -t 2:0:0
#SBATCH --mem=6G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-$n


locus=\$SLURM_ARRAY_TASK_ID

read chr pos SNP   < <(sed -n \${locus}p fPCloci.txt)


cd $workDir

if [ \$chr == "X" ]; then
    
./plink/plink2 \
  --pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/FPCphenotypes_doubleIDs_EUR.txt \
  --pheno-name fpc1-fpc10 \
  --covar  $workDir/pheno/covariates_doubleIDs_EUR.txt \
  --not-covar sex \
  --condition \$SNP \
  --vif 500  \
  --xchr-model 2 \
  --covar-variance-standardize \
  --glm hide-covar cols=-chrom,-test,-nobs,-err,+a1count,+a1freq  \
  --threads 2 \
  --out $workDir/results/chr\${chr}_conditional_\$SNP



else


./plink/plink2 \
--pfile $dataDir/plink2Bin/EUR_minMaf0.005_minInfo0.8_chr\${chr} \
  --pheno $workDir/pheno/FPCphenotypes_doubleIDs_EUR.txt \
  --pheno-name fpc1-fpc10 \
  --covar  $workDir/pheno/covariates_doubleIDs_EUR.txt \
--condition \$SNP \
--vif 500  \
--covar-variance-standardize \
--glm hide-covar cols=-chrom,-test,-nobs,-err,+a1count,+a1freq  \
--threads 2 \
--out $workDir/results/chr\${chr}_conditional_\$SNP


fi

 rsync -av $workDir/results/chr\${chr}_conditional_\${SNP}*.glm.linear /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output

EOF

 sbatch $workDir/scripts/conditionalAnalysesFPCs.sh


#  rsync -av $workDir/results/chr*Pixel*_conditional_*.glm.linear /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/cojoAnalysis/output
