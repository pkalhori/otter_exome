//Number of population samples (demes)
2
//Population effective sizes (number of genes)
NCUR0
NCUR1
//Sample sizes
14
12
//Growth rates : negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 0
0 0
//Migration matrix 1
0 Mig_Back_0_1
Mig_Back_1_0 0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix
3 historical events
56 0 0 0 RESIZE_0 0 1
73 1 1 0 RESIZE_1 0 1
TDIV 1 0 1 RESIZE_TOTAL 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 8.6e-9 OUTEXP
