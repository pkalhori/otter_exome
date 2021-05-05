Daniel_vcfdir=/u/project/rwayne/software/rails/VCF_FILES_MISSING_SITES

vcfdir=/u/scratch/p/pkalhori/rails/VCFs_Missing_sites

source /u/local/Modules/default/init/modules.sh
module load samtools
module load bcftools

allVCF=Neutral_sites_SFS_ALL_1.vcf.gz

cp $Daniel_vcfdir/$allVCF $vcfdir

gunzip $vcfdir/$allVCF ##it must be unzipped
#extract only SNPs
snpVCF=Neutral_sites_SNPs_only_1.vcf
bcftools view -c 1:minor ${vcfdir}/${allVCF} > ${vcfdir}/${snpVCF}