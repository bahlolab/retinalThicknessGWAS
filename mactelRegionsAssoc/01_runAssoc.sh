
workDir=/vast/scratch/users/jackson.v/retThickness/macTelLociAssocs
dataDir=/vast/scratch/users/jackson.v/retThickness/GWAS/filteredGeno

mkdir -p $workDir
cd $workDir

# mkdir -p $workDir/plink
# cd $workDir/plink
# wget "https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20220814.zip"
# unzip plink2_linux_x86_64_20220814.zip

mkdir -p $workDir/pheno
mkdir -p $workDir/results
mkdir -p $workDir/logs
mkdir -p $workDir/scripts


rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/covariates* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/phenotypes* $workDir/pheno
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/pixels.txt .
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/rawData/macTelLoci.txt .
rsync -av /vast/scratch/users/jackson.v/retThickness/plink .


read Rsid chr BP from to < <(sed -n 1p macTelLoci.txt)
slice=1

./plink/plink2 \
  --bfile $dataDir/ukbb_minMAC500_minInfo0.8_chr${chr} \
  --pheno $workDir/pheno/phenotypesSlice${slice}_doubleIDs.txt \
  --covar  $workDir/pheno/covariates_doubleIDs.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=+a1count,+a1freq \
  --chr $chr \
  --from-bp $from \
  --to-bp $to \
  --threads 2 \
  --out $workDir/results/${Rsid}/${slice}/${Rsid}Pixel




nSNPs=$(wc -l macTelLoci.txt)

for snp in $(seq 1 $nSNPs)
do

  read Rsid chr BP from to < <(sed -n ${snp}p macTelLoci.txt)


cat <<- EOF > $workDir/scripts/plinkAssoc_${Rsid}.sh
#!/bin/sh

#SBATCH -J plink-${Rsid}
#SBATCH -o $workDir/logs/plinkChr${Rsid}_%A_%a.log
#SBATCH -t 3:0:0
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-119%20


slice=\$SLURM_ARRAY_TASK_ID


cd $workDir


mkdir -p $workDir/results/${Rsid}/\${slice}

./plink/plink2 \
  --bfile $dataDir/ukbb_minMAC500_minInfo0.8_chr${chr} \
  --pheno $workDir/pheno/phenotypesSlice\${slice}_doubleIDs.txt \
  --covar  $workDir/pheno/covariates_doubleIDs.txt \
  --vif 500  \
  --covar-variance-standardize \
  --glm hide-covar cols=+a1count,+a1freq \
  --chr $chr \
  --from-bp $from \
  --to-bp $to \
  --threads 2 \
  --out $workDir/results/${Rsid}/\${slice}/${Rsid}Pixel


EOF

 sbatch $workDir/scripts/plinkAssoc_${Rsid}.sh
 sleep 1
done

rsync -av $workDir/results/* /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/output



nSNPs=$(wc -l macTelLoci.txt 

for snp in $(seq 1 $nSNPs)
do

  read Rsid chr BP from to < <(sed -n ${snp}p macTelLoci.txt)

    rm /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/logs/${Rsid}_logs.txt
    for fi in $workDir/logs/plinkChr${Rsid}_*.log
    do

        echo "$fi $(tail -n 1 $fi)" >> /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/logs/${Rsid}_logs.txt
    done
done


# Rsid=rs2160387
# jobID=8499830

# Rsid=rs1047891
# jobID=8499819

# Rsid=rs10995566
# jobID=8490772

# Rsid=rs9820465
# jobID=8490769

# Rsid=rs677622 
# jobID=8498526

# Rsid=rs17279437
# jobID=8490775

#     for fi in $workDir/logs/plinkChr${Rsid}_${jobID}_*.log
#     do

#         echo "$fi $(tail -n 1 $fi)" >> /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/macTelLociAssocs/logs/${Rsid}_${jobID}_logs.txt
#     done
