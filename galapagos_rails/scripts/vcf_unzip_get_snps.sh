#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=10G
#$ -N SLiMInference
#$ -o /u/scratch/p/pkalhori
#$ -e /u/scratch/p/pkalhori
#$ -m abe
#$ -M pkalhori
#$ -t 10-15:1

Daniel_vcfdir=/u/project/rwayne/software/rails/VCF_FILES_MISSING_SITES

vcfdir=/u/scratch/p/pkalhori/rails/VCFs_Missing_sites

source /u/local/Modules/default/init/modules.sh
module load samtools
module load bcftools

allVCF=Neutral_sites_SFS_ALL_${SGE_TASK_ID}.vcf.gz

cp $Daniel_vcfdir/$allVCF $vcfdir

gunzip $vcfdir/$allVCF ##it must be unzipped

allVCF_unzipped=Neutral_sites_SFS_ALL_${SGE_TASK_ID}.vcf

#extract only SNPs
snpVCF=Neutral_sites_SNPs_only_${SGE_TASK_ID}.vcf
bcftools view -c 1:minor ${vcfdir}/${allVCF_unzipped} > ${vcfdir}/${snpVCF}
