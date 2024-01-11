#!/bin/bash
#SBATCH --job-name=linkimpute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=30G
#SBATCH --mail-type=END
#SBATCH --mail-user=katherine.l.taylor@uconn.edu
#SBATCH -o ./errandout/linkimpute%j.out
#SBATCH -e ./errandout/linkimpute%j.err

cd ~/song_genes_reanalysis/VCF

#java -jar ~/tools/LinkImpute.jar -q song_family_filter.ped song_family_filter_imputeC.ped

java -jar ~/tools/LinkImpute.jar -q song_family_filter_fixed.ped song_family_filter_fixed_impute.ped