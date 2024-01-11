#!/bin/bash
#SBATCH --job-name=entap
#SBATCH --mail-user=katherine.l.taylor@uconn.edu
#SBATCH --mail-type=all
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH -o ./errandout/entap_%j.out
#SBATCH -e ./errandout/entap_%j.err
#SBATCH --partition=general
#SBATCH --qos=general

module load EnTAP/0.10.8
module load anaconda/2.4.0
module load perl/5.30.1
module load diamond/2.0.6
module load interproscan/5.25-64.0
module load TransDecoder/5.3.0

EnTAP --runP --ini  ~/song_genes_reanalysis/scripts/entap_config.ini  \
-i ~/song_genes_reanalysis/reference_genome/ChrCarn1.1_protein.faa \
-d /isg/shared/databases/Diamond/Uniprot/uniprot_sprot.fa.2.0.6.dmnd \
--threads 8 --out-dir ~/song_genes_reanalysis/reference_genome/entap/ChrCarn1.1