#wd=/u/scratch/p/pkalhori/fastsimcoal/bird_data
#wd=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/fastsimcoal
wd=/u/scratch/p/pkalhori/fastsimcoal/slim_simulated_inference
models="1D.1Epoch_CA 1D.2Epoch_CA 1D.3Epoch_CA"
for rep in {1..10}
do
pop=CA.simulation.rep.${rep}.CA.1D.2Epoch.35Gen.200Inds
#muts="4.6e-9"
rundate=20200708
for model in $models
do
outfile=$wd/${model}_${pop}_${rundate}.all.output.concatted.txt
#get header:
header=`head -n1 $wd/${model}_${pop}_${rundate}/run_1/${model}_${pop}/*bestlhoods`
echo -e "runNum\t$header" > $outfile
for i in {1..50}
do 
outdir=$wd/${model}_${pop}_${rundate}/run_${i}/${model}_${pop}
results=`grep -v [A-Z] $outdir/*.bestlhoods`
echo -e "${i}\t$results" >> $outfile
done
done
done
