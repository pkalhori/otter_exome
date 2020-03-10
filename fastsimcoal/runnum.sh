#wd=/u/scratch/p/pkalhori/fastsimcoal/bird_data
wd=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/fastsimcoal
models="2D.2Epoch 2D.2Epoch.Mig.Symmetric 2D.3Epoch.FixedContraction.Time.Sizes 2D.3Epoch.Mig.Symmetric.FixedContraction.Time.Sizes"
pops=neutral.CA.AK
#muts="4.6e-9"
rundate=20200305
for model in $models
do
for pop in $pops
do
#for mut in $muts
#do
outfile=$wd/${model}_${pop}_${rundate}.all.output.concatted.txt
#get header:
header=`head -n1 $wd/${model}_${pop}_${rundate}/run_1/${model}_${pop}/*bestlhoods`
echo -e "runNum\t$header" > $outfile
for i in {1..50}
do 
outdir=$wd/${model}_${pop}_${rundate}/run_${i}/${model}_${pop}
results=`grep -v [A-Z] $outdir/*.bestlhoods`
echo -e "${i}\t$results" >> $outfile
#done
done
done
done
