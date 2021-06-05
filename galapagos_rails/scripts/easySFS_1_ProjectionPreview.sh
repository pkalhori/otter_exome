#! /bin/bash
#$ -cwd
#$ -l rh7,h_data=50G
#$ -N easySFSPreview
#$ -o /u/scratch/p/pkalhori/rails/logs
#$ -e /u/scratch/p/pkalhori/rails/logs
#$ -m abe
#$ -M pkalhori
#$ -t 1-35:1
#$ -V
####### Easy SFS
# https://github.com/isaacovercast/easySFS
# install:
# git clone git@github.com:isaacovercast/easySFS.git
# cd easySFS
# chmod +x *.py
# easySFS.py


source /u/local/Modules/default/init/modules.sh
#module load python/2.7
module load samtools
module load anaconda3/2020.11
. /u/local/apps/anaconda3/2020.11/etc/profile.d/conda.sh

conda activate dadi

#bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip
maxHetFilter=0.75 # het filter used across all samples (per population het filter occurs during easy sfs)
todaysdate=`date +%Y%m%d`
#genotypeDate=20181119
vcfdir=/u/scratch/p/pkalhori/rails/VCFs_Missing_sites_Test
scriptdir=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/galapagos_rails/scripts
popFile=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/galapagos_rails/pinta.PCA.txt # this doesn't have baja on it; doesn't have any admixed/bad inds on it. 
# this has admixed in it , but they aren't in pop file
easySFS=$scriptdir/easySFS.abModified.3.noInteract.Exclude01Sites.HetFiltering.20181121.py  # this is my modification
outdir=/u/scratch/p/pkalhori/rails/easySFS/projection_preview/$todaysdate
mkdir -p $outdir
# had to modify easySFS so that it wouldn't prompt a "yes/no" response about samples that are missing from VCF file
# want it to just output that info and continue , not prompt yes/no.
# this vcf has all snps across all categories (cds, neutral, etc.) with 0.9 max no call frac (v. liberal)
# and has had all individuals removed that won't go into the SFS
# going to do the actual projection for each category of site
#vcf=neutral.snp_8b_forEasySFS_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf
#vcf=neutral.snp_9b_maxHetFilter_${maxHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf
vcf=Neutral_sites_SNPs_only_${SGE_TASK_ID}.vcf
#gunzip $vcfdir/populationVCFs/$vcf it must be unzipped
#( you are here )
$easySFS -i $vcfdir/${vcf} -p $popFile --preview -a -v > $outdir/neutral.snp_${SGE_TASK_ID}.easySFS.projPreview.txt

### now plot projections in R and decide on your levels. Actually DO the projections on 

# need to add -a otherwise it selects one snp per locus (was built for RAD data)
# eventually make popFile list consistent with final dataset individuals
# otherwise have to say "yes" on screen
    #Running preview mode. We will print out the results for # of segregating sites
    #for multiple values of projecting down for each population. The dadi
    #manual recommends maximizing the # of seg sites for projections, but also
   	#a balance must be struck between # of seg sites and sample size.

    #For each population you should choose the value of the projection that looks
   # best and then rerun easySFS with the `--proj` flag.
