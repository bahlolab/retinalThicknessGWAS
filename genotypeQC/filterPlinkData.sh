#!/bin/bash

while getopts p:l:k:e:o: option
do
  case "${option}"
  in
    p) plinkPath=${OPTARG};;
    l) listPlinkFiles=${OPTARG};;
    k) keepSamps=${OPTARG};;
    e) extractSNPs=${OPTARG};;
    o) outputPrefix=${OPTARG};;
  esac
done

export PATH=${PATH}:${plinkPath}

mkdir cleaningTemp

rand=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)

## merge chromosomes together
plink \
  --merge-list ${listPlinkFiles} \
  --make-bed \
  --out cleaningTemp/allChr_$rand

## extract variants and samples for GRM calculations
plink \
  --bfile cleaningTemp/allChr_$rand \
  --keep ${keepSamps} \
  --extract ${extractSNPs} \
  --make-bed \
  --out ${outputPrefix}
