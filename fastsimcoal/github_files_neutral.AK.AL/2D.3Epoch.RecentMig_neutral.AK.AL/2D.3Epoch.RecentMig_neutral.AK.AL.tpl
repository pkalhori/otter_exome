//Number of population samples (demes)
2
//Population effective sizes (number of genes)
3600
780
//Sample sizes
SAMPLE_SIZE_0
SAMPLE_SIZE_1
//Growth rates : negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0
0 R_Mig_0_1
R_Mig_1_0 0
//Migration matrix 1
0 Mig_Back_0_1
Mig_Back_1_0 0
//Migration matrix 2
0 0
0 0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix
3 historical events
56 1 1 0 11.282 0 1
306 0 0 0 2.694 0 1
TDIV 1 0 1 RESIZE_TOTAL 0 2
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 8.6e-9 OUTEXP
