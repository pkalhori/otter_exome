\#! /bin/bash
#$ -cwd
#$ -l h_rt=20:00:00,h_data=10G
#$ -N {name of job}
#$ -o {path to output files}
#$ -e {path to error files"
#$ -m abe
#$ -M {username}
#$ -t {number of runs in an array}

deme0={population 1}
deme1={population 2}
pops="neutral.${deme1}.${deme0}"
models={models to loop through}
rundate=`date +%Y%m%d`
#this is the string of populations to loop through
# wd stands for "working directory"
wd=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/fastsimcoal
md=$wd/modeldir
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
cp $md/$model.tpl ${model}_${pop}.tpl
cp $md/$model.est ${model}_${pop}.est
cp $wd/SFSdir/${pop}.NewNames_jointMAFpop1_0.obs ${model}_${pop}_jointMAFpop1_0.obs
mkdir $wd/$header/run_${SGE_TASK_ID}
cd $wd/$header/run_${SGE_TASK_ID}
cp $wd/$header/${model}_${pop}.tpl $wd/$header/${model}_${pop}.est $wd/$header/${model}_${pop}_jointMAFpop1_0.obs ./
ss0=`grep -w $deme0 $wd/projectionValues.txt|awk '{print$2}'`
sed -i "s/SAMPLE_SIZE_0/$ss0/g" ${model}_${pop}.tpl
ss1=`grep -w $deme1 $wd/projectionValues.txt|awk '{print$2}'`
sed -i "s/SAMPLE_SIZE_1/$ss1/g" ${model}_${pop}.tpl
$fsc -t ${model}_${pop}.tpl -n100000 -m -e ${model}_${pop}.est -M -L 50 
done
done
cd $wd
sleep 10m 


