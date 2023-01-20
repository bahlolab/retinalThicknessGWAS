

chr=22
sentinelsList=/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/chr${chr}sentinelsIDonly.txt

nSNPs=$(wc -l  $sentinelsList | cut -d " " -f 1)


  cat <<- EOF > /vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/scripts/regionPlots_chr${chr}.sh
#!/bin/bash

#SBATCH -J extractGWsig
#SBATCH -o /vast/scratch/users/jackson.v/retThickness/GWAS/logs/extractResultGWsig_%A_%a.log
#SBATCH -t 1:0:0
#SBATCH --mem=10GB
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jackson.v@wehi.edu.au
#SBATCH -a 1-$nSNPs

cd /vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/chr${chr}

export PATH=\$PATH:/wehisan/bioinf/lab_bahlo/users/jackson.v/resources/plink_linux_x86_64_20190617/

sentinelsList=/wehisan/bioinf/lab_bahlo/projects/misc/retinalThickness/GWAS/output/chr${chr}sentinelsIDonly.txt

snp=\$SLURM_ARRAY_TASK_ID


  read sentinel < <(sed -n \${snp}p \$sentinelsList)

  file=/vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/chr${chr}/\${sentinel}_METAL.txt


/wehisan/bioinf/lab_bahlo/users/jackson.v/resources/locuszoom/bin/locuszoom \
  --metal \$file \
  --refsnp \$sentinel \
  --flank 500kb \
  --markercol ID \
  --pvalcol P \
  --pop EUR \
  --build hg19 \
  --source 1000G_March2012 \
  --prefix \$sentinel \
  showAnnot=T \
  showRecomb=T | tee \${sentinel}_regionplot.txt

done

rsync -av ./\${sentinel}*/*.pdf /vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/chr${chr}

EOF

sbatch -a 1 /vast/scratch/users/jackson.v/retThickness/GWAS/regionPlots/scripts/regionPlots_chr${chr}.sh

