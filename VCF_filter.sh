#!/bin/bash
#SBATCH --job-name=VCF_filter
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=1G
#SBATCH --mail-type=END
#SBATCH --mail-user=katherine.l.taylor@uconn.edu
#SBATCH -o ./errandout/VCF_filter%j.out
#SBATCH -e ./errandout/vcf_filter%j.err

module load bcftools

cd ~/song_genes_reanalysis/VCF

# Filter VCF for selection analysis
bcftools view --samples "^136_adamsi.bam,137_adamsi.bam,162_adamsi.bam,163_adamsi.bam,171_adamsi.bam,172_adamsi.bam,DL20.bam,132_plorabunda.bam,133_plorabunda.bam,181_plorabunda.bam,182_plorabunda.bam,184_plorabunda.bam,CF78.bam" \
song_selection.vcf > song_selection_pa.vcf

bcftools view --min-alleles 2 --max-alleles 2 \
--include 'TYPE="snp" & QUAL>30 & F_MISSING<0.8 & MAF>0.05' \
song_selection_pa.vcf > song_selection_pa_filter.vcf

bcftools stats song_selection.vcf > song_selection_stats.txt
bcftools stats song_selection_pa_filter.vcf > song_selection_pa_filter_stats.txt

# Filter VCF for family

bcftools view --min-alleles 2 --max-alleles 2 \
--include 'TYPE="snp" & QUAL>30 & F_MISSING<0.8 & MAF>0.05' \
song_family.vcf > song_family_filter.vcf


bcftools stats song_family.vcf > song_family_stats.txt
bcftools stats song_family_filter.vcf > song_family_filter_stats.txt

## Extract a VCF with just F0 
bcftools view ../VCF/song_family_filter.vcf -s CF78.bam,DL20.bam -o ../VCF/family_filter_f0.vcf

## Extract VCF with no het calls for parents
bcftools view family_filter_f0.vcf --genotype ^het -o family_filter_f0_fixed.vcf

## Extract VCF with no het or missing calls for parents

bcftools view family_filter_f0_fixed.vcf --genotype ^miss -o family_filter_f0_fixed_nomiss.vcf

# Extract VCF with just F1 females

bcftools view ../VCF/song_family_filter.vcf -s F1F1.bam,F1F2.bam,F1F3.bam,F1F4.bam,F1F5.bam,F1F6.bam -o ../VCF/family_filter_f1.vcf

## Extract VCF with no hom calls for F1 females

bcftools view family_filter_f1.vcf --genotype ^hom -o family_filter_f1_het.vcf

## Extract with no hom or missing calls for F1 females

bcftools view family_filter_f1_het.vcf --genotype ^miss -o family_filter_f1_het_nomiss.vcf

## Make gzipped files

bcftools view -I song_family_filter.vcf -O z -o song_family_filter.vcf.gz
bcftools view -I family_filter_f0_fixed_nomiss.vcf -O z -o family_filter_f0_fixed_nomiss.vcf.gz
bcftools view -I family_filter_f1_het_nomiss.vcf -O z -o family_filter_f1_het_nomiss.vcf.gz

bcftools index song_family_filter.vcf.gz
bcftools index family_filter_f0_fixed_nomiss.vcf.gz
bcftools index family_filter_f1_het_nomiss.vcf.gz

##Intersect the three files

bcftools isec -n=3 -w1 song_family_filter.vcf.gz family_filter_f0_fixed_nomiss.vcf.gz family_filter_f1_het_nomiss.vcf.gz -o song_family_filter_parfixed_f1het.vcf 

## Do hwe filtering for family cross
bcftools +fill-tags song_family_filter_parfixed_f1het.vcf | bcftools view -e'HWE<=0.01' > song_family_filter_parfixed_f1het_hwe.vcf

## Calculate stats on that dataset

bcftools stats song_family_filter_parfixed_f1het_hwe.vcf > song_family_filter_parfixed_f1het_stats.txt

bcftools query -f '%DP\n' song_family_filter_parfixed_f1het_hwe.vcf > DP.txt

## Do file convert for linkimpute

module load plink/1.90.beta.4.4 
plink --vcf song_family_filter_parfixed_f1het_hwe.vcf --make-bed \
--const-fid "progeny" --out song_family_filter_fixed --allow-extra-chr 

plink --bfile song_family_filter_fixed --allow-extra-chr \
--pheno pxa_phen.txt --update-sex sex.txt --recode --make-bed \
--out song_family_filter_fixed --const-fid "progeny"


