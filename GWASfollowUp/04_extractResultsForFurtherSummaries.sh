

## extract pixelwise results for all FfPC sentinels, and plot.
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsFPCsentinels/logs/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsFPCsentinels/scripts/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsFPCsentinels/results/

mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsKnownLoci/logs/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsKnownLoci/scripts/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsKnownLoci/results/

cd  /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsFPCsentinels/

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWASfollowUp/extractFPCSentinelsPixelwiseResults.sh scripts

sbatch ./scripts/extractFPCSentinelsPixelwiseResults.sh


cd  /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsKnownLoci/

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWASfollowUp/extractpixelWiseResultsKnownLociPixelwiseResults.sh scripts

sbatch ./scripts/extractpixelWiseResultsKnownLociPixelwiseResults.sh

