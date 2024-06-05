

## extract pixelwise results for all FfPC sentinels, and plot.
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsFPCsentinels/logs/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsFPCsentinels/scripts/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsFPCsentinels/results/

mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsPixOnlysentinels/logs/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsPixOnlysentinels/scripts/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsPixOnlysentinels/results/


mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsKnownLoci/logs/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsKnownLoci/scripts/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsKnownLoci/results/


mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsAllSentinels/logs/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsAllSentinels/scripts/
mkdir -p /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsAllSentinels/results/



cd  /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsFPCsentinels/

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWASfollowUp/extractFPCSentinelsPixelwiseResults.sh scripts

sbatch ./scripts/extractFPCSentinelsPixelwiseResults.sh


cd  /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsKnownLoci/

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWASfollowUp/extractpixelWiseResultsKnownLociPixelwiseResults.sh scripts

sbatch ./scripts/extractpixelWiseResultsKnownLociPixelwiseResults.sh


cd  /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsPixOnlysentinels/

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWASfollowUp/extractpixelWiseResultsPixOnly.sh scripts

sbatch ./scripts/extractpixelWiseResultsPixOnly.sh


cd  /vast/scratch/users/jackson.v/retThickness/GWASfollowup/pixelWiseResultsAllsentinels/

rsync -av /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/retinalThicknessGWAS/GWASfollowUp/extractpixelWiseResultsAllSentinels.sh scripts

sbatch ./scripts/extractpixelWiseResultsAllSentinels.sh