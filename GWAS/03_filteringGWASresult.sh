
module load R/4.1.2

cd /vast/scratch/users/jackson.v/retThickness/GWAS/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/

rsync -av /stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWAS/filteringGwasResults1.R ./scripts
rsync -av /stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWAS/filteringGwasResults2.sh ./scripts
rsync -av /stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWAS/collatingResults.* ./scripts

for chr in {1..22}
do

echo "filtering chromosome $chr"

mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr$chr/

## Run SNP cleaning scripts
./scripts/filteringGwasResults1.R  \
  --chr $chr \
  --directory /vast/scratch/users/jackson.v/retThickness/GWAS/clumpedResults 

done

chr=X
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/GWsigResults/chr$chr/

## Run SNP cleaning scripts
./scripts/filteringGwasResults1.R  \
  --chr $chr \
  --directory /vast/scratch/users/jackson.v/retThickness/GWAS/clumpedResults 


## extract SNPs
sbatch ./scripts/filteringGwasResults2.sh

## run collating results
## this defines signals - clumps of SNPs, and lists of pixels, assigned to each signal.
sbatch ./scripts/collatingResults.sh

## extracts results for sentinels for all pixels
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/sentinelResults/
rsync  -av /stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWAS/extractSentinelResults.sh ./scripts
sbatch ./scripts/extractSentinelResults.sh 

## generate pixelwise and region plots
mkdir -p /stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/
mkdir -p /stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/regionPlots/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/scripts
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/logs/

rsync  -av /stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWAS/runningPlots*.sh ./scripts
rsync  -av /stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWAS/creatingLocusZoomInput.R ./scripts

for chr in {1..22}
do
mkdir -p /stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/chr${chr}
mkdir -p  /vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/chr${chr}

## run
./scripts/runningPlots.sh ${chr}

done

chr=X

mkdir -p /stornext/Bioinf/data/lab_bahlo/projects/misc/retinalThickness/GWAS/output/sentinels/plots/chr${chr}
mkdir -p  /vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/chr${chr}

## run
./scripts/runningPlotsX.sh 
