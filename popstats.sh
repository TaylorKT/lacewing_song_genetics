#!/bin/bash
#SBATCH --job-name=popstats
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=6G
#SBATCH --mail-type=END
#SBATCH --mail-user=katherine.l.taylor@uconn.edu
#SBATCH -o ./errandout/popstats%j.out
#SBATCH -e ./errandout/popstats%j.err

module load vcftools/0.1.16
cd ~/song_genes_reanalysis/popstats

## Calculate FST
vcftools --vcf ../VCF/song_selection_filter.vcf \
--weir-fst-pop plora.txt --weir-fst-pop adam.txt --fst-window-size 50000 \
--fst-window-step 5000 --out weirfst_plora_vs_adam

# Calculate pi

vcftools --vcf ../VCF/song_selection_filter.vcf \
--keep plora_adam.txt --recode --recode-INFO-all --out song_selection_filter_pa

vcftools --vcf song_selection_filter_pa.recode.vcf \
--window-pi 50000 --window-pi-step 5000 --out pi_plora_vs_adam