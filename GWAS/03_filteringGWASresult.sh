
module load R/4.1.2

cd /vast/scratch/users/jackson.v/retThickness/GWAS/
mkdir /vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWAS/filteringGwasResults1.R ./scripts
rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWAS/filteringGwasResults2.sh ./scripts

for chr in {12..22}
do

echo "filtering chromosome $chr"

mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr$chr/

## Run SNP cleaning scripts
./scripts/filteringGwasResults1.R  \
  --chr $chr \
  --directory /vast/scratch/users/jackson.v/retThickness/GWAS/clumpedResults 

done

## extract SNPs

sbatch ./scripts/filteringGwasResults2.sh

## run collating results
collatingResults_v1.r

## extracts results for sentinels for all pixels
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/sentinelResults/

rsync  -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWAS/extractSentinelResults.sh ./scripts


## pixelwise plots
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/
for chr in {1..22}
do
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/chr${chr}
done
##run 
plottingSentinelsPixelwise.R

## region plots
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/regionPlots/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/scripts
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/logs/

for chr in {1..22}
do
mkdir -p /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/chr${chr}
mkdir -p  /vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/chr${chr}
done

## run
runningLocusZoom.sh