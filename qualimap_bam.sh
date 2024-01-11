#!/bin/bash
#SBATCH --job-name=qualimap
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mem=1G
#SBATCH --mail-type=END
#SBATCH --mail-user=katherine.l.taylor@uconn.edu
#SBATCH -o ./errandout/qualimap/qualimap_%j.out
#SBATCH -e ./errandout/qualimap/qualimap_%j.err
#SBATCH --array=[0-121]%20

cd /home/FCAM/ktaylor/song_genes_reanalysis/aligned

BAMS=($(ls -1 /home/FCAM/ktaylor/song_genes_reanalysis/aligned/*.bam))
INFILE=${BAMS[$SLURM_ARRAY_TASK_ID]}

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID
echo $INFILE

module load qualimap

qualimap bamqc -bam $INFILE