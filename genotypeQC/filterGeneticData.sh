#!/bin/bash

while getopts q:p:g:s:x:k:e:o:n: option
do
  case "${option}"
  in
    q) qctoolPath=${OPTARG};;
    p) plinkPath=${OPTARG};;
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
mkdir -p cleaningTemp/filtScripts
mkdir -p cleaningTemp/filtErrors


for chr in $(seq 1 22)
do

    cat <<- EOF > cleaningTemp/filtScripts/plinkFiltering_chr${chr}_${outName}.sh
#!/bin/sh

#SBATCH -J chr${chr}_${outName}
#SBATCH -o cleaningTemp/plinkErrors/CHR${chr}_${outName}
#SBATCH -t 48:0:0
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au

module load gcc
export PATH=\${PATH}:${qctoolPath}

qctool \
-g ${genDataDir}/ukb_imp_chr${chr}_v3.bgen  \
-s ${sampFile} \
-incl-snpids ${extractSNPs} \
-incl-samples ${keepSamps} \
-og $outputDir/bgenFilt/${outName}_chr${chr}.bgen \
-os $outputDir/bgenFilt/${outName}_chr${chr}.sample

${plinkPath}/plink2 \
  --bgen $outputDir/bgenFilt/${outName}_chr${chr}.bgen ref-first \
  --sample $outputDir/bgenFilt/${outName}_chr${chr}.sample \
  --rm-dup exclude-all \
  --threads 2 \
  --memory 300000 \
  --make-pgen \
  --out $outputDir/plink2Bin/${outName}_chr${chr}

${plinkPath}/plink2 \
  --pfile $outputDir/plink2Bin/${outName}_chr${chr} \
  --threads 2 \
  --make-bed \
  --out $outputDir/plinkBin/${outName}_chr${chr}

EOF


  sbatch cleaningTemp/filtScripts/plinkFiltering_chr${chr}_${outName}.sh

  sleep 5
done

chr=X
cat <<- EOF > cleaningTemp/plinkScripts/plinkFiltering_chr${chr}_${outName}.sh
#!/bin/sh

#SBATCH -J plink-${chr}_${outName}
#SBATCH -o cleaningTemp/plinkErrors/CHR${chr}_${outName}
#SBATCH -t 48:0:0
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au


module load gcc
export PATH=\${PATH}:${qctoolPath}

qctool \
-g ${genDataDir}/ukb_imp_chr${chr}_v3.bgen  \
-s ${xSampFile} \
-incl-snpids ${extractSNPs} \
-incl-samples ${keepSamps} \
-og $outputDir/bgenFilt/${outName}_chr${chr}.bgen \
-os $outputDir/bgenFilt/${outName}_chr${chr}.sample

${plinkPath}/plink2 \
  --bgen $outputDir/bgenFilt/${outName}_chr${chr}.bgen ref-first \
  --sample $outputDir/bgenFilt/${outName}_chr${chr}.sample \
  --rm-dup exclude-all \
  --threads 2 \
  --memory 150000 \
  --make-pgen \
  --out $outputDir/plink2Bin/${outName}_chr${chr}

${plinkPath}/plink2 \
  --pfile $outputDir/plink2Bin/${outName}_chr${chr} \
  --threads 2 \
  --make-bed \
  --out $outputDir/plinkBin/${outName}_chr${chr}

EOF

sbatch cleaningTemp/plinkScripts/plinkFiltering_chr${chr}_${outName}.sh







