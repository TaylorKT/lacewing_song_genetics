#!/bin/bash
#SBATCH --job-name=mpileup_family
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=6G
#SBATCH --mail-type=END
#SBATCH --mail-user=katherine.l.taylor@uconn.edu
#SBATCH -o ./errandout/mpileup_family%j.out
#SBATCH -e ./errandout/mpileup_family%j.err

cd ~/song_genes_reanalysis/aligned

module load bcftools/1.9

bcftools mpileup -f ~/song_genes_reanalysis/reference_genome/ChrCarn1.1_genomic.fasta \
-b ./bam_list_family.txt --threads 1 -O u -o ../VCF/song_family.bcf

bcftools call -vmO v ../VCF/song_family.bcf -o ../VCF/song_family.vcf