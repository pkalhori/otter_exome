#! /bin/bash
#$ -cwd
#$ -l h_rt=20:00:00,h_data=28G
#$ -N SLiMInference
#$ -o /u/scratch/p/pkalhori/fastsimcoal/reports
#$ -e /u/scratch/p/pkalhori/fastsimcoal/reports
#$ -m abe
#$ -M pkalhori
#$ -t 1-50:1

rundate=`date +%Y%m%d`
models="1D.1Epoch_CA 1D.2Epoch_CA 1D.3Epoch_CA"
for model in $models
do
for num in {1..10}
do
pops=CA.simulation.rep.${num}.CA.1D.2Epoch.35Gen.200Inds

#this is the string of populations to loop through
# wd stands for "working directory"
wd=/u/scratch/p/pkalhori/fastsimcoal/slim_simulated_inference
SFSdir=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/fastsimcoal/SFSes/fastsimcoal2Format_folded
md=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/fastsimcoal/github_files_CA
fsc=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/software/fsc26_linux64/fsc26
for model in $models
#iterate through each model
do
for pop in $pops
#iterate through each population
do
cd $wd
header=${model}_${pop}_$rundate
mkdir $header
cd $header
cp $md/$model/$model.tpl ${model}_${pop}.tpl
cp $md/$model/$model.est ${model}_${pop}.est
cp $wd/SFSdir/${pop}_MAFpop0.obs ${model}_${pop}_MAFpop0.obs
mkdir $wd/$header/run_${SGE_TASK_ID}
cd $wd/$header/run_${SGE_TASK_ID}
cp $wd/$header/${model}_${pop}.tpl $wd/$header/${model}_${pop}.est $wd/$header/${model}_${pop}_MAFpop0.obs ./
#ss=`grep -w $pop $wd/projectionValues.txt|awk '{print$2}'`
#sed -i "s/SAMPLE_SIZE/$ss/g" ${model}_${pop}.tpl
$fsc -t ${model}_${pop}.tpl -n100000 -m -e ${model}_${pop}.est -M -L 50 -q
done
done
done
cd $wd
sleep 5m 


