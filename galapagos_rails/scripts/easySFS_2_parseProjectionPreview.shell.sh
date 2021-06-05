####### Parse EasySFS output
projection_date=20210530
wd=/u/scratch/p/pkalhori/rails/easySFS/projection_preview/$projection_date

for i in {1..35}
do

easyOut=neutral.snp_$i.easySFS.projPreview.txt
pops="PIN"

for pop in $pops
do
# put a "$" at the end of pop to signify that it should be at the end of the line
# get the line from the easy SFS output; get rid of of the parentheses and make a column
echo "projection,snps" > $wd/$pop.${easyOut%.txt}.R.format.txt
grep -A1 "$pop$" $wd/$easyOut | tail -1 | sed 's/(//g' | sed 's/)//g' | sed 's/, /,/g' |  tr '\t' '\n' >> $wd/$pop.${easyOut%.txt}.R.format.txt
done
done
