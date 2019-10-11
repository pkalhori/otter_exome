states="PreContraction PostContraction PostRecovery"                                                                                                                        for k in $(seq 50002 2 50034)
for k in $(seq 50002 2 50034)
do
states="$states Contraction.${k}gen"
done
for l in $(seq 50038 2 50052)
do
states="$states Recovery.${l}gen"
done
for m in $(seq 50056 2 50104)
do
states="$states Future.${m}gen"
done
#states="PreContraction PostContraction PostRecovery SpillRecovery.${j}gen PostSpill"
