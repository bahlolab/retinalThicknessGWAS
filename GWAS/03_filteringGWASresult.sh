
module load R/4.1.2

cd /vast/scratch/users/jackson.v/retThickness/GWAS/
mkdir /vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWAS/filteringGwasResults1.R ./scripts

for chr in {12..22}
do

echo "filtering chromosome $chr"

mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr$chr/


## Run SNP cleaning scripts
./scripts/filteringGwasResults1.R  \
  --chr $chr \
  --directory /vast/scratch/users/jackson.v/retThickness/GWAS/clumpedResults 

done

