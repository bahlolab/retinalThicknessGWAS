#!/bin/sh

#SBATCH -J runningMTAG
#SBATCH -o /wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWAS/output/MTAGinput/mtag.log
#SBATCH -t 16:0:0
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au

# module load git
# module load python

# cd /vast/scratch/users/jackson.v/retThickness/fpcGWAS/
# git clone https://github.com/omeed-maghzian/mtag.git
# cd mtag

# mkdir output



module load python

cd /vast/scratch/users/jackson.v/retThickness/fpcGWAS/mtag

inFiles=("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWAS/output/MTAGinput/allChr_fpc"{1..25}"_MTAG.txt")

# inFiles=("/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/fpcGWAS/output/MTAGinput/chr5fpc"{1..10}"_MTAG.txt")


delim=""
joined=""
for item in "${inFiles[@]}"; do
  joined="$joined$delim$item"
  delim=","
done
inFilesList=$(echo "$joined")


python mtag.py  \
	--sumstats $inFilesList \
	--out ../mtagOutput/allChr_fpcs1to25 \
    --snp_name ID \
    --n_name OBS_CT \
    --z_name T_STAT \
    --eaf_name A1_FREQ \
    --maf_min 0 \
    --force \
    --fdr \
    --stream_stdout 