#!/bin/bash
#SBATCH --job-name=ld
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=75G
#SBATCH --mail-type=END
#SBATCH --mail-user=katherine.l.taylor@uconn.edu
#SBATCH -o ./errandout/ld%j.out
#SBATCH -e ./errandout/ld%j.err

module load plink/1.90.beta.4.4
module load bcftools/1.9 

cd ../ld


bcftools view -I ../VCF/song_selection_cp_filter.vcf -O z -o ../VCF/song_selection_pa_filter.vcf.gz
bcftools index ../VCF/song_selection_pa_filter.vcf.gz


bcftools view ../VCF/song_selection_pa_filter.vcf.gz --regions NC_058337.1 \
-s 132_plorabunda.bam,133_plorabunda.bam,134_plorabunda.bam,181_plorabunda.bam,182_plorabunda.bam,184_plorabunda.bam,CF78.bam \
-o song_selection_p_filter_1.vcf 

bcftools view ../VCF/song_selection_pa_filter.vcf.gz --regions NC_058338.1 \
-s 132_plorabunda.bam,133_plorabunda.bam,134_plorabunda.bam,181_plorabunda.bam,182_plorabunda.bam,184_plorabunda.bam,CF78.bam \
-o song_selection_p_filter_2.vcf 

bcftools view ../VCF/song_selection_pa_filter.vcf.gz --regions NC_058339.1 \
-s 132_plorabunda.bam,133_plorabunda.bam,134_plorabunda.bam,181_plorabunda.bam,182_plorabunda.bam,184_plorabunda.bam,CF78.bam \
-o song_selection_p_filter_3.vcf 

bcftools view ../VCF/song_selection_pa_filter.vcf.gz --regions NC_058340.1 \
-s 132_plorabunda.bam,133_plorabunda.bam,134_plorabunda.bam,181_plorabunda.bam,182_plorabunda.bam,184_plorabunda.bam,CF78.bam \
-o song_selection_p_filter_4.vcf 

bcftools view ../VCF/song_selection_pa_filter.vcf.gz --regions NC_058341.1 \
-s 132_plorabunda.bam,133_plorabunda.bam,134_plorabunda.bam,181_plorabunda.bam,182_plorabunda.bam,184_plorabunda.bam,CF78.bam \
-o song_selection_p_filter_5.vcf 

bcftools view ../VCF/song_selection_pa_filter.vcf.gz --regions NC_058342.1 \
-s 132_plorabunda.bam,133_plorabunda.bam,134_plorabunda.bam,181_plorabunda.bam,182_plorabunda.bam,184_plorabunda.bam,CF78.bam \
-o song_selection_p_filter_6.vcf 

bcftools view ../VCF/song_selection_pa_filter.vcf.gz --regions NC_058337.1 \
-s DL20.bam,136_adamsi.bam,137_adamsi.bam,162_adamsi.bam,163_adamsi.bam,171_adamsi.bam,172_adamsi.bam  \
-o song_selection_a_filter_1.vcf 

bcftools view ../VCF/song_selection_pa_filter.vcf.gz --regions NC_058338.1 \
-s DL20.bam,136_adamsi.bam,137_adamsi.bam,162_adamsi.bam,163_adamsi.bam,171_adamsi.bam,172_adamsi.bam  \
-o song_selection_a_filter_2.vcf 

bcftools view ../VCF/song_selection_pa_filter.vcf.gz --regions NC_058339.1 \
-s DL20.bam,136_adamsi.bam,137_adamsi.bam,162_adamsi.bam,163_adamsi.bam,171_adamsi.bam,172_adamsi.bam  \
-o song_selection_a_filter_3.vcf 

bcftools view ../VCF/song_selection_pa_filter.vcf.gz --regions NC_058340.1 \
-s DL20.bam,136_adamsi.bam,137_adamsi.bam,162_adamsi.bam,163_adamsi.bam,171_adamsi.bam,172_adamsi.bam  \
-o song_selection_a_filter_4.vcf 

bcftools view ../VCF/song_selection_pa_filter.vcf.gz --regions NC_058341.1 \
-s DL20.bam,136_adamsi.bam,137_adamsi.bam,162_adamsi.bam,163_adamsi.bam,171_adamsi.bam,172_adamsi.bam  \
-o song_selection_a_filter_5.vcf 

bcftools view ../VCF/song_selection_pa_filter.vcf.gz --regions NC_058342.1 \
-s DL20.bam,136_adamsi.bam,137_adamsi.bam,162_adamsi.bam,163_adamsi.bam,171_adamsi.bam,172_adamsi.bam  \
-o song_selection_a_filter_6.vcf 

cd ./chr1/

plink --vcf ../song_selection_p_filter_1.vcf \
--r2 --double-id --allow-extra-chr \
--ld-window-r2 0 --ld-window-kb 200000  --ld-window 999999 

cd ../chr2/

plink --vcf ../song_selection_p_filter_2.vcf \
--r2 --double-id --allow-extra-chr \
--ld-window-r2 0 --ld-window-kb 200000  --ld-window 999999 

cd ../chr3/

plink --vcf ../song_selection_p_filter_3.vcf \
--r2 --double-id --allow-extra-chr \
--ld-window-r2 0 --ld-window-kb 200000  --ld-window 999999 

cd ../chr4/

plink --vcf ../song_selection_p_filter_4.vcf \
--r2 --double-id --allow-extra-chr \
--ld-window-r2 0 --ld-window-kb 200000  --ld-window 999999 

cd ../chr5/

plink --vcf ../song_selection_p_filter_5.vcf \
--r2 --double-id --allow-extra-chr \
--ld-window-r2 0 --ld-window-kb 200000  --ld-window 999999 

cd ../chr6/

plink --vcf ../song_selection_p_filter_6.vcf \
--r2 --double-id --allow-extra-chr \
--ld-window-r2 0 --ld-window-kb 200000  --ld-window 999999 

cd ..

bcftools view ../VCF/song_selection_pa_filter.vcf.gz \
-s 132_plorabunda.bam,133_plorabunda.bam,134_plorabunda.bam,181_plorabunda.bam,182_plorabunda.bam,184_plorabunda.bam,CF78.bam \
-o song_selection_p_filter.vcf 

plink --vcf song_selection_p_filter.vcf \
--r2 --double-id --allow-extra-chr \
--ld-window-r2 0 --ld-window-kb 200000  --ld-window 999999 

bcftools view ../VCF/song_selection_pa_filter.vcf.gz \
-s DL20.bam,136_adamsi.bam,137_adamsi.bam,162_adamsi.bam,163_adamsi.bam,171_adamsi.bam,172_adamsi.bam \
-o song_selection_a_filter.vcf 

plink --vcf song_selection_a_filter.vcf \
--r2 --double-id --allow-extra-chr \
--ld-window-r2 0 --ld-window-kb 200000  --ld-window 999999 

#cat plink.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > plink.ld.summary
