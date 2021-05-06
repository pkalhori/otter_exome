#! /bin/bash
#$ -cwd
#$ -l h_rt=20:00:00,h_data=16G
#$ -N easySFSProjection2
#$ -o /u/scratch/pkalhori/rails
#$ -e /u/scratch/pkalhori/rails
#$ -m abe
#$ -M pkalhori

####### Easy SFS
#https://github.com/isaacovercast/easySFS
#install:
#git clone git@github.com:isaacovercast/easySFS.git
#cd easySFS
#chmod +x *.py
#easySFS.py
source /u/local/Modules/default/init/modules.sh
module load python/2.7
module load samtools
module load bcftools
#bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip
#todaysdate=20181212
todaysdate=`date +%Y%m%d`
genotypeDate=20181119
noCallFrac=1.0
vcfdir=/u/project/rwayne/software/rails/VCF_FILES
pops=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/galapagos_rails/pinta.PCA.txt
popfile=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/galapagos_rails/pinta.PCA.txt
# this has admixed in it , but they aren't in pop file
easySFS=/u/home/a/ab08028/klohmueldata/annabel_data/bin/easySFS/easySFS.abContinueMod.py 

gitdir=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/galapagos_rails
scriptdir=${gitdir}/scripts


#easySFS=$scriptdir/easySFS.abModified.3.noInteract.Exclude01Sites.HetFiltering.20181121.py  # this is my modification
# this version of script excludes sites that are 0-1 across all populations (maybe) -- not sure if it does yet. 

## choose your projections: choosing for now: 
# CA,AK,AL,COM,KUR 
#projections="12,14,16,16,12" # may change these after lab meeting
projections="18"
### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 

outdir=/u/scratch/pkalhori/rails/easySFS/projection-${todaysdate}
snpVCFdir=/u/scratch/pkalhori/rails/snpVCFs
mkdir -p $outdir
mkdir -p $snpVCFdir
# had to modify easySFS so that it wouldn't prompt a "yes/no" response about samples that are missing from VCF file
# write projection choices into a readme

echo "PIN : $projections " > $outdir/projectionChoices.${todaysdate}.txt
# make sure vcf isn't zipped

allVCF=Neutral_sites_SFS_ALL_1.vcf.gz
#extract only SNPs
snpVCF=Neutral_sites_SNPs_only_1.vcf
bcftools view -c 1:minor ${vcfdir}/${allVCF} > ${snpVCFdir}/${snpVCF}


### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 
$easySFS -i $snpVCFdir/${snpVCF} -p $popfile -a -v --proj $projections -f -o $outdir
# test run only had --proj 20 (maybe just projects 1 population? will have to see)
# -f forces overwrite of outdir
# $bgzip ${vcf}
# then do for SYN and MIS (eventually)
########## get counts of monomorphic sites to add to the SFSes ############
python $scriptdir/getMonomorphicProjectionCounts.1D.2DSFS.py --vcf $vcfdir/${allSitesvcf} --popMap $pops --proj $projections --popIDs PIN --outdir $outdir

