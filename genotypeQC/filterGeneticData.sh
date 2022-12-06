#!/bin/bash

while getopts q:g:s:x:k:e:o:n: option
do
  case "${option}"
  in
    q) plinkPath=${OPTARG};;
    g) genDataDir=${OPTARG};;
    s) sampFile=${OPTARG};;
    x) xSampFile=${OPTARG};;
    k) keepSamps=${OPTARG};;
    e) extractSNPs=${OPTARG};;
    o) outputDir=${OPTARG};;
    n) outName=${OPTARG};;
  esac
done

## make temp directories
mkdir -p cleaningTemp/plinkScripts
mkdir -p cleaningTemp/plinkErrors


for chr in $(seq 1 22)
do

    cat <<- EOF > cleaningTemp/plinkScripts/plinkFiltering_chr${chr}_${outName}.sh
#!/bin/sh

#SBATCH -J plink-${chr}_${outName}
#SBATCH -o cleaningTemp/plinkErrors/CHR${chr}_${outName}
#SBATCH -t 48:0:0
#SBATCH --mem=60G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au


${plinkPath}/plink2 \
  --bgen ${genDataDir}/ukb_imp_chr${chr}_v3.bgen ref-first \
  --sample ${sampFile} \
  --extract ${extractSNPs} \
  --keep ${keepSamps} \
  --threads 2 \
  --make-pgen \
  --out $outputDir/plink2Bin/${outName}_chr${chr}

${plinkPath}/plink2 \
  --bgen ${genDataDir}/ukb_imp_chr${chr}_v3.bgen ref-first \
  --sample ${sampFile} \
  --extract ${extractSNPs} \
  --keep ${keepSamps} \
  --threads 2 \
  --make-bed \
  --out $outputDir/plinkBin/${outName}_chr${chr}

EOF


  sbatch cleaningTemp/plinkScripts/plinkFiltering_chr${chr}_${outName}.sh

  sleep 5
done

chr=X
cat <<- EOF > cleaningTemp/plinkScripts/plinkFiltering_chr${chr}_${outName}.sh
#!/bin/sh

#SBATCH -J plink-${chr}_${outName}
#SBATCH -o cleaningTemp/plinkErrors/CHR${chr}_${outName}
#SBATCH -t 48:0:0
#SBATCH --mem=60G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au


${plinkPath}/plink2 \
  --bgen ${genDataDir}/ukb_imp_chr${chr}_v3.bgen ref-first \
  --sample ${xSampFile} \
  --extract ${extractSNPs} \
  --keep ${keepSamps} \
  --threads 2 \
  --make-pgen \
  --out $outputDir/plink2Bin/${outName}_chr${chr}

${plinkPath}/plink2 \
  --bgen ${genDataDir}/ukb_imp_chr${chr}_v3.bgen ref-first \
  --sample ${xSampFile} \
  --extract ${extractSNPs} \
  --keep ${keepSamps} \
  --threads 2 \
  --make-bed \
  --out $outputDir/plinkBin/${outName}_chr${chr}

EOF

sbatch cleaningTemp/plinkScripts/plinkFiltering_chr${chr}_${outName}.sh
