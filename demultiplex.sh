#!/bin/bash
#SBATCH --job-name=demultiplex
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=10G
#SBATCH --mail-type=END
#SBATCH --mail-user=katie.leah.taylor@gmail.com
#SBATCH -o ./errandout/demultiplex_%j.out
#SBATCH -e ./errandout/demultiplex_%j.err

cd ~/song_genes/

module load stacks/2.2

process_radtags -f ./raw_rad/pooled.fastq -o ./demultiplexed -b ./barcodes.txt -e pstI -c -q -t 140 -s 30
