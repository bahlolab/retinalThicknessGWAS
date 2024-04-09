#!/bin/bash

#SBATCH -J scanSplit
#SBATCH -t 2:0:0
#SBATCH --mem=1G
#SBATCH --output=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/ERRORS/scanSplit_%A_%a.log
#SBATCH -a 1-994


job=$SLURM_ARRAY_TASK_ID


let "end=$job * 101"
let "start=$end-100"


for fi in $(seq $start $end)
do
read PATIENT SCAN < <(sed -n ${fi}p "/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/scanIDs.txt")

mkdir -p /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/rawPerScan/$PATIENT

FILE=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/rawPerScan/${PATIENT}/${SCAN}_scan.csv

  if [ -f $FILE ]; then
    # nLines=$(wc -l < "$FILE")
    #
    # if [ $nLines -lt 128 ]; then

      # echo "Incomplete File: $FILE"
      echo "Removing File: $FILE"
      rm -f $FILE

      echo "Creating File: $FILE"
      for i in {0..127}
      do
        grep $SCAN /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/rawData/height_${i}.csv >> $FILE
     done

    # else
    #
    #  echo "File $FILE exists."

   # fi

  else

    echo "Creating File: $FILE"
    for i in {0..127}
    do
      grep $SCAN /vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/rawData/height_${i}.csv >> $FILE
    done
  fi



done
