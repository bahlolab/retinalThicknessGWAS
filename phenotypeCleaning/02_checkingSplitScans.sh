#!/bin/bash

cd /vast/scratch/users/jackson.v/retThickness/pheno/

start=1
end=$(wc -l < "/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/scanIDs.txt")


for fi in $(seq $start $end)
do
read PATIENT SCAN < <(sed -n ${fi}p "/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/scanIDs.txt")

FILE=/vast/projects/bahlo_ukbiobank/app28541_retinal/retinalThickness/data/rawPerScan/${PATIENT}/${SCAN}_scan.csv

  if [ -f $FILE ]; then

   nLines=$(wc -l < "$FILE")
   lastEdit=$(date -r "$FILE"  +"%Y-%m-%d")


  echo "$nLines, $lastEdit, $FILE" >> filesFullSummary.csv


  else

    echo "Missing, NA, $FILE" >> filesFullSummary.csv

  fi
done
