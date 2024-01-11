#!/bin/bash
#SBATCH --job-name=gemma
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=10G
#SBATCH --mail-type=END
#SBATCH --mail-user=katherine.l.taylor@uconn.edu
#SBATCH -o ./errandout/gemma%j.out
#SBATCH -e ./errandout/gemma%j.err

module load plink/1.90.beta.4.4
cd ~/song_genes_reanalysis/gemma

### Format file for gemma
plink --ped ../VCF/song_family_filter_fixed_impute.ped --map ../VCF/song_family_filter_fixed.map \
--make-bed --out song_family_filter_fixed_impute_gemma --allow-extra-chr \
--hwe 0.01 --maf 0.05 --set-missing-var-ids @:#

##### Gemma LMM

~/tools/gemma-0.98.4-linux-static-AMD64 -bfile song_family_filter_fixed_impute_gemma -lm 4 -o song_family

##### Gemma BSLMM

~/tools/gemma-0.98.4-linux-static-AMD64 -bfile song_family_filter_fixed_impute_gemma -bslmm 1 -o song_family_BSLMM_1 -w 500000 -s 5000000
~/tools/gemma-0.98.4-linux-static-AMD64 -bfile song_family_filter_fixed_impute_gemma -bslmm 1 -o song_family_BSLMM_2 -w 500000 -s 5000000
~/tools/gemma-0.98.4-linux-static-AMD64 -bfile song_family_filter_fixed_impute_gemma -bslmm 1 -o song_family_BSLMM_3 -w 500000 -s 5000000
~/tools/gemma-0.98.4-linux-static-AMD64 -bfile song_family_filter_fixed_impute_gemma -bslmm 1 -o song_family_BSLMM_4 -w 500000 -s 5000000
~/tools/gemma-0.98.4-linux-static-AMD64 -bfile song_family_filter_fixed_impute_gemma -bslmm 1 -o song_family_BSLMM_5 -w 500000 -s 5000000
