#!/bin/bash
#SBATCH --job-name=song_bwa
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mem=8G
#SBATCH --mail-type=END
#SBATCH --mail-user=katherine.l.taylor@uconn.edu
#SBATCH -o ./errandout/song_bwa_%j.out
#SBATCH -e ./errandout/song_bwa_%j.err


files="F2M13
		F2F37
		F2M5
		F2M36
		F2M27
		F2M19
		F2M10
		F2F20
		F2M21
		F2M34
		F2F6
		F2M40
		F2M8
		F2M22
		F2M20
		F2M39
		F2F8
		F2M17
		F2M4
		F2M29
		F2F36
		F2M9
		F2F1
		F2M28
		F2M11
		F2M30
		F2M35
		F2M32
		F2M26
		F2M3
		F2F31
		F2F23
		F2F7
		F2F34
		F2F27
		F2F29
		F2F2
		F2M6
		F2F24
		F2F30
		F2F12
		F2F5
		F2F28
		F2M7
		F2M24
		F2F11
		F2F33
		F2F3
		F2M42
		F2F18
		F2M18
		F2F17
		F2M16
		F2M25
		F2F26
		F2F16
		F2F25
		F2F15
		F2F14
		F2F4
		F2F32
		F2F22
		F2F35
		F2F10
		F2F21
		F2F13
		F2M38
		F2M31
		F2M12
		F2M41
		F2F9
		F2M14
		F2M15
		F2M23
		F2M44
		F2M37
		F2M43
		F2M33
		F2M45
		F2M46
		CF78
		DL20
		F1F2
		F1M2
		F1M1
		F1M5
		F1F4
		F1F5
		F1F3
		F1M3
		F1F1
		F1M6
		F1M4
		F1F6
		F1112
		FGXCONTROL"

cd ~/song_genes_reanalysis/aligned

module load bwa/0.7.17

for file in $files
do
bwa mem -t 8 ../reference_genome/Ccarn1 ~/song_genes/demultiplexed/${file}.fq.gz > ${file}.sam
done  

module load samtools/1.7

for file in $files
do 
samtools view -b ${file}.sam | samtools sort --threads 8 > ${file}.bam;
done

for file in $files
do 
    echo ${file};
    samtools flagstat ${file}.sam;
done

