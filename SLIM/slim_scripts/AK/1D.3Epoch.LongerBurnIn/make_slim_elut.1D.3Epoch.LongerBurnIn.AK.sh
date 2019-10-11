# Script to make SLIM job script
# USAGE: ./make_slim_wolf_101517.job.sh [h]
pop=AK
model=1D.3Epoch.LongerBurnIn
#model=1D.2Epoch.1.5Mb.cds.LongerContract
gitdir=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/SLIM
scriptdir=$gitdir/slim_scripts/$pop/$model
mkdir -p $scriptdir # set script dir
todaysdate=`date +%Y%m%d`
########## population specific parameters ########
# set sample size for vcf (should match empirical for AK)
# can set manually or pull from a table
ss=7 # in diploids # will depend on population, note you can find these in projectionValues.txt file
#population	ProjectionValueHaploids Diploids
#KUR	12	6
#CA	12	6
#COM	34	17
#AK	14	7
#AL	20	10
# variables:
#stats_sample=30
nanc=4500 # a ball park Nanc for AK
nu=250 # contraction size  -- a medium contraction size
tcontract=35 # contraction duration before you sample -- a longer contraction time than dadi ifnerence (more like FSC)
trecovery=17 #1911-2019 is around 17 generations (6y/gen)
tfuture=50 #allow simulation to continue running at recovered pop size  
nrec=2500 #not sure what population size to use; best I could find for WPWS
# this nu / T combo is within the MLE estimate for Alaska from my grid search. also want to try the California version (100 for 30 gens with smaller Nanc of ~3500)
######### general parameters ; can set here or in command line ##########
# Set g, number of genes (exons)
g=1000
# Set gLen, the length of exons
gLen=1500
# Set t, number of burn-in generations
t=150000
# set mutation rate
mu=8.64e-9 # mutation rate
# Set h, dominance coefficient
h=$1  # loop through hs
# Set j, the chunk number (for 14 chunks?)

#total number of chromosomes, since all "chunks" will be simulated together 


# Make script
# chunk gets set when you run slim (based on SGE task id) # or something?? how to do this part? I don't really want to make a separte slim script each time? have it be a -d thing maybe?
# have to figure out the chunks/replicates situation.

cat > $scriptdir/slim_elut_${model}_${pop}_h${h}.job << EOM

// changes to make: apparently 1e-03 is reasonable between-gene recomb rate
// and then want to make separate chromosomes with 0.5 between them (or just simulate them separately)
// okay I think I am going to model 1.5Mb stretches of seqeunce, each containing 1000 genes of size 1500bp
// recomb w/in each gene will be 1e-08, between genes will be 1e-3
// will either do 14x1.5mb so they are independent, or will simulate them all together if want to think about overall genetic load (then would have to set up chromosomes within the simulation)
initialize() {
	defineConstant("g",$g); //number of genes; starting with 1000 (AB)
	defineConstant("geneLength", $gLen); // length of each gene
	defineConstant("seqLength", g*geneLength); // total sequence length starting with 1.5Mb (AB)
	//defineConstant("outdir",\"$outdir\"); -- set in command line
	//defineConstant("v_CHUNK",$chunk); // portion of genome  -- set in command line
	//defineConstant("v_REP",$rep); // overall replicate number -- set in command line
	defineConstant("v_h",$h); // dominance coefficient
	defineConstant("v_SS",$ss); // sample size
	defineConstant("v_MU",$mu);
	defineConstant("v_NANC",$nanc); // ancestral size
	defineConstant("v_NU",$nu); // contraction size
	defineConstant("v_NREC",$nrec); // contraction size
	defineConstant("sampleSize", 30);

	//cat("Exome portion length:"+seqLength+"\n");
	initializeMutationRate(v_MU);
	// m1 mutation type: neutral *synonymous*
	initializeMutationType("m1", 0.5, "f", 0.0);
	// m2 mutation type: missense(?) -- This is from Chris; ask where he got params (Kim et al?)
	initializeMutationType("m2", v_h, "g",-0.01314833, 0.186); // set H in array (0 for recessive, 0.5 for additive) 
	m2.convertToSubstitution = F; // okay cool if you have this, then fixed sites are in the vcf file
	m1.convertToSubstitution = F; // keeps fixed sites in vcf file (do I want?)
	// initialize exon: g1, has both neutral (m1) and misense (m2) mutations
	// syn happens at rate 1, mis at rate 2.31:1 since there are more missense sites (Christian/Chris); note chris names neutral as m2 and missense at m1 so I reversed things here 
	initializeGenomicElementType("g1", c(m1,m2), c(1.0,2.31)); // 2.31 is the NS:S ratio
	
	for (i in 0:(g-1)){
		initializeGenomicElement(g1, ((i)*geneLength+1), ((i+1)*geneLength));
	}
	// figured out a way for recomb to be 0.5 between blocks, but 1e-08 within blocks
	// so that they aren't actually linked; simulating  independent genes/exons
	//initializeRecombinationRate(1e-8);
	// this sets up rates that alternate between 1e-08 and 1e-03
	// and ends that are in pattern 999 1000 1999 2000 2999 ...
	// so that the pattern is that r is 1e-8 for 0-999, then 1e-03 (a reasonable between-gene recomb rate) between 999 and 1000 (each gene), then is 1e-8 through next gene, and so on.
	// making 5000 independent blocks.
	rates=c(rep(c(1e-08,1e-3),g));
	ends=NULL;
	for (index in 0:(g-1))
	{
		ends=c(ends,index*geneLength+(geneLength-1),index*geneLength+geneLength);
	}
	initializeRecombinationRate(rates,ends);


}

//Chris's function for getting stats 
function (s) getStats(o pop, i sampSize)
{
	i = sample(pop.individuals, sampSize, F); //do i want to try doinng this for the whole popultion or will it take too long?
	
	m = sortBy(i.genomes.mutations, "position"); //get all mutations in sample
	m_uniq = unique(m); // get rid of redundant muts
	DAF = sapply(m_uniq, "sum(m == applyValue);"); // count number of each mut in pop
	m_uniq_polym = m_uniq[DAF != i.genomes.size()]; //remove fixed mutations??
	
	//initialize vectors
	//ROH_length_sumPerInd = c();
	Num_VstrDel_muts = c();
	Num_strDel_muts = c();
	Num_modDel_muts = c();
	Num_wkDel_muts = c();
	ind_het = c();
	//fitness_population = c();
	
	for (individual in i) {
		
		indm = sortBy(individual.genomes.mutations, "position");
		indm = indm[match(indm, m_uniq_polym) >= 0];   // Check that individual mutations are not fixed 
		indm_uniq = unique(indm);
		
		genotype = sapply(indm_uniq, "sum(indm == applyValue);");
		
		// tally number of mutations for different classes of selection coefficient per individual
		s = individual.genomes.mutations.selectionCoeff;
		
		Num_VstrDel_muts = c(Num_VstrDel_muts, sum(s<=-0.05));
		Num_strDel_muts = c(Num_strDel_muts, sum(s<=-0.01));
		Num_modDel_muts = c(Num_modDel_muts, sum(s<=-0.001 & s > -0.01));
		Num_wkDel_muts = c(Num_wkDel_muts, sum(s<=-0.00001 & s > -0.001));
		
		if (isNULL(genotype)) {
			ind_het = c(ind_het, 0); //putting this here to avoid error when trying to sum null vector
			next;
		}
		
		ind_het = c(ind_het, sum(genotype==1)/(seqLength));
		
	}
	
	return(pop.individuals.size() + "," + "," + mean(ind_het) + "," + "," + mean(Num_VstrDel_muts) + "," + mean(Num_strDel_muts)+ "," + mean(Num_modDel_muts) + "," + mean(Num_wkDel_muts));
}


// create a population of variable v_NANC individuals
1 {
	sim.addSubpop("p1", v_NANC);
			// set up outfile: 
	writeFile(paste(c(outdir,"/slim.output.BurnIn.",v_CHUNK,".summary.txt"),sep=""),"replicate,chunk,generation,popSize,meanHet,avgVstrDel,avgStrDel,avgModDel,avgWkDel\n",append=F); // open fresh file
}

// output generation number so I can track progress and stats during burn in 

1:${t} late() {
	if (sim.generation % 1000 == 0){
		cat(sim.generation+"\n");
		stats = getStats(p1, sampleSize);
		//add data to outfile
		writeFile(paste(c(outdir,"/slim.output.BurnIn.",v_CHUNK,".summary.txt"),sep=""),paste(c(v_REP,v_CHUNK,sim.generation,stats),sep=","),append=T);
	}
}
// after t generation burn in, sample individuals and output counts across whole population as well (from JAR script) ;
// then do this again after the contraction -- then only need to simulate once instead of doing 1 and 2 epoch separately. 
${t} late() {
	p1.outputVCFSample(v_SS, F,filePath=paste(c(outdir,"/slim.output.PreContraction.",v_CHUNK,".vcf"),sep=""));
	}
	
// set up outfile: 
${t} late() {
	writeFile(paste(c(outdir,"/slim.output.",v_CHUNK,".summary.txt"),sep=""),"replicate,chunk,generation,mutid,type,s,age,subpop,p1numhet,p1numhom,popsizeDIP\n",append=F); // open fresh file
	}
//Calculate load every two generations after burn in until end of simulation
${t}: late() {
	if (sim.generation % 2 == 0){
	//file header
	//mutation id
	//mutation type
	//selection coefficient
	//age of mutation in generations
	//subpopulation it arose in
	//number of heterozygote derived in p1
	//number of homozygote derived in p1
	//number of heterozygote derived in p2
	//number of homozygote derived in p2
	//these are genotype counts not allele counts

	//for every mutation in the simulation
	for (mut in sim.mutations){
		id = mut.id;
		s = mut.selectionCoeff;
		generation= sim.generation - 50000;
		subpop = mut.subpopID;
		age = sim.generation - mut.originGeneration;
		type = mut.mutationType;
		popsize=size(p1.individuals);
		//initialize genotype counts
		p1numhet = 0;
		p1numhom = 0;
		
		//count hom and het derived in p1
		for (p1i in p1.individuals){
			gt = sum(c(p1i.genomes[1].containsMutations(mut), p1i.genomes[0].containsMutations(mut)));
			if (gt == 1){
				p1numhet = p1numhet + 1;
			} else if (gt == 2){
				p1numhom = p1numhom + 1;
			}
		}
		
	
		// string for mutation type. add m3, m4, etc. if you have multiple types
		if (type == m1){
			type = "m1";
		} else if (type == m2){
			type = "m2";
		}
		
		//print results
		writeFile(paste(c(outdir,"/slim.output.",v_CHUNK,".summary.txt"),sep=""),paste(c(v_REP,v_CHUNK,generation,id,type,s,age,subpop,p1numhet,p1numhom,popsize),sep=","),append=T);
	}
}
}

// contract the population 1 gen after burn in:
$((${t} + 1)) {
	p1.setSubpopulationSize(v_NU);
	}

// Tcollect VCF file at the end of contraction 
$((${t} + 1+ ${tcontract})) late() {
	p1.outputVCFSample(v_SS, F,filePath=paste(c(outdir,"/slim.output.PostContraction.",v_CHUNK,".vcf"),sep=""));
}
// expand the population 1 generation after sampling (recovery):
$((${t} + 2 + ${tcontract})) {
	p1.setSubpopulationSize(v_NREC);
	}

// collect VCF and summary at present day
$((${t} + 2+ ${tcontract}+${trecovery})) late() {
	p1.outputVCFSample(v_SS, F,filePath=paste(c(outdir,"/slim.output.PostRecovery.",v_CHUNK,".vcf"),sep=""));
}

//Collect VCF at the end of Simulatoin 
$((${t} + 3+ ${tcontract}+${trecovery}+${tfuture})) late() {
	p1.outputVCFSample(v_SS, F,filePath=paste(c(outdir,"/slim.output.EndSimulation.",v_CHUNK,".vcf"),sep=""));
	//file header

}
$((${t} + 3+ ${tcontract}+${trecovery}+${tfuture})) late() {
sim.simulationFinished();
} 
EOM
