#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=6G
#$ -m abe
#$ -e /u/scratch/p/pkalhori/slim/R_load_calc/logs
#$ -o /u/scratch/p/pkalhori/slim/R_load_calc/logs
#$ -M pkalhori


pop=CA_AK
model=2D.3Epoch.Translocation.25for2Gen
wd=/u/scratch/p/pkalhori/slim/R_load_calc
ld=$wd/$pop/$model
mkdir -p $ld

pathToScript=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/SLIM/R_scripts/Linux_Format
myscript=load_calcs_2D.R

source /u/local/Modules/default/init/modules.sh
module load R
Rscript $pathToScript/$myscript

cp load_data.txt $ld

 
sleep 2m
