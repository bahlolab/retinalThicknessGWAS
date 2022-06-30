#!/bin/bash

while getopts q:g:s:x:k:e:o:n: option
do
  case "${option}"
  in
    q) qctoolPath=${OPTARG};;
    g) genDataDir=${OPTARG};;
    s) sampFile=${OPTARG};;
    x) xSampFile=${OPTARG};;
    k) keepSamps=${OPTARG};;
    e) extractSNPs=${OPTARG};;
    o) outputPrefix=${OPTARG};;
    n) jobName=${OPTARG};;
  esac
done

## make temp directories
mkdir -p cleaningTemp/qctoolScripts
mkdir -p cleaningTemp/qctoolErrors

## Create cleaned,  data

for chr in $(seq 1 22)
do

    cat <<- EOF > cleaningTemp/qctoolScripts/qctoolFiltering_chr${chr}_${jobName}.sh
#!/bin/sh

#SBATCH -J qctool-${chr}
#SBATCH -o cleaningTemp/qctoolErrors/CHR${chr}_${jobName}
#SBATCH -t 48:0:0
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au

module load gcc
export PATH=\${PATH}:${qctoolPath}

qctool \
-g ${genDataDir}/ukb_imp_chr${chr}_v3.bgen  \
-s ${sampFile} \
-incl-snpids ${extractSNPs} \
-incl-samples ${keepSamps} \
-og ${outputPrefix}_chr${chr}.bgen \
-os ${outputPrefix}_chr${chr}.sample

EOF


  sbatch cleaningTemp/qctoolScripts/qctoolFiltering_chr${chr}_${jobName}.sh

  sleep 5
done

chr=X
cat <<- EOF > cleaningTemp/qctoolScripts/qctoolFiltering_chr${chr}_${jobName}.sh
#!/bin/sh

#SBATCH -J qctool-${chr}
#SBATCH -o cleaningTemp/qctoolErrors/CHR${chr}_${jobName}
#SBATCH -t 48:0:0
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au

module load gcc
export PATH=\${PATH}:${qctoolPath}

qctool \
-g ${genDataDir}/ukb_imp_chr${chr}_v3.bgen  \
-s ${xSampFile} \
-incl-snpids ${extractSNPs} \
-incl-samples ${keepSamps} \
-og ${outputPrefix}_chr${chr}.bgen \
-os ${outputPrefix}_chr${chr}.sample

EOF

sbatch cleaningTemp/qctoolScripts/qctoolFiltering_chr${chr}_${jobName}.sh
