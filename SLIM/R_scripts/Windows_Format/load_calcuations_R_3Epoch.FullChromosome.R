require(ggplot2)
todaysdate=format(Sys.Date(),format="%Y%m%d")


data.dir="C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\"

plot.dir="C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\genetic_load_calcs\\"
dir.create(plot.dir)
#pops=c("AK","AL","genericPop.LongerContract")
#models=c("1D.2Epoch.1.5Mb.cds")
#simdates=c(20190424,20190607)
# skipping AL "AL/1D.2Epoch.1.5Mb.cds/20190424/" and CA etc -- add those in next.

popModDates=c("AK\\1D.3Epoch.LongerRecovery.FullChromosome\\20190828")

#popModDates=c("AK/1D.2Epoch.1.5Mb.cds/20190424/","AK/1D.2Epoch.1.5Mb.cds.LongerContract/20190607","genericPop/1D.2Epoch.1.5Mb.cds.20KAncSize/20190611") # AK and AL have dadi parameters, genericPop has parameters based on AK MLE grid that is fur-trade relevant. ### need to come up with better classification system for this. 
reps=c(seq(1,10))
#reps=c(seq(1,25)) # some reps don't make it through Hoffman; so I have a file.exists() test in the loop to skip 
#reps that didn't yield output



hset=c(0,0.5)

states=c("PreContraction.50000gen","PostContraction.50036gen","PostRecovery.50054gen")
for (gen in seq(50002,50034,2)){
  x<- paste("Contraction.",gen,"gen",sep = "")
  states<- append(states,x)
}
for (gen in seq(50038,50052,2)){
  x<- paste("Recovery.",gen,"gen",sep = "")
  states<- append(states,x)
}
for (gen in seq(50056,50104,2)){
  x<- paste("Future.",gen,"gen",sep = "")
  states<- append(states,x)
}

allLoads=data.frame()
states[3]

for(popModDate in popModDates){
  #for(model in models){
  #for(simdate in simdates){
  #  for(pop in pops){
  for(rep in reps){
    for(state in states){
      for(h in hset){
        # check if rep exists (some have random hoffman failures)
        infile=paste(data.dir,popModDate,"\\h_",h,"\\replicate_",rep,"\\slim.output.",state,".1.summary.txt",sep="")
        #infile=paste(data.dir,popModDate,"\\h_",h,"\\replicate_",rep,".slim.output.",state,".allConcatted.summary.txt.gz",sep="")
        if(file.exists(infile)){
          input = read.table(infile,sep=",",header=T)
          # calculate q (alt allele frequency) and p (ref allele frequency) per site
          # be careful about which you use in equation! s*q^2 means that q is frequency of ALT allele with associated "s". so p is freq of ref allele
          input$qFreq <- (input$p1numhet + (2*input$p1numhom)) / (2*input$popsizeDIP)
          input$pFreq <- 1 - input$qFreq
          # equation is, per site: 2hspq + sq^2 is the contribution to load; p = 1-qFreq
          # my "s" is negative, so I want to absolute value s --> |s|
          input$loadComponent <- (2*h*abs(input$s)*input$qFreq*input$pFreq) + (abs(input$s)*((input$qFreq)^2))
          # total sites:
          generation= -1*input$generation[1]
          S = sum(input$loadComponent)
          W = exp(-S) # mean fitness e^-S
          L  = 1 - W # mutation load 
          #### add to dataframe: #####
          # pull out population etc from popModDate
          pop= unlist(lapply(strsplit(popModDate, "[\\\\]|[^[:print:]]"),"[",1))
          model= unlist(lapply(strsplit(popModDate,"[\\\\]|[^[:print:]]"),"[",2))
          date= unlist(lapply(strsplit(popModDate,"[\\\\]|[^[:print:]]"),"[",3))
          loadDF <- data.frame(population=pop)
          loadDF$model <- model
          loadDF$date <- date
          loadDF$generation <- generation
          loadDF$rep <- rep
          loadDF$state <- state
          loadDF$h <- h
          loadDF$S_allsites <- S
          loadDF$W_meanFitness <- W
          loadDF$L_mutationLoad <- L
          #### combine with other reps: #####
          allLoads = rbind(allLoads,loadDF)
        }
      }
    }
  }
}

allLoads$generation <- (allLoads$generation)*-1


# change order of factors:
allLoads$state <- factor(allLoads$state,levels=c("PreContraction","PostContraction","PostRecovery"))
# label H:
allLoads$hLabel <- paste("h = ",allLoads$h)
allLoads$model<- "1D.3Epoch.LongerRecovery.FullChromosome"
allLoads$model <- factor(allLoads$model,levels=c("1D.5Epoch"))
p1 <- 
  ggplot(allLoads,aes(x=generation,y=L_mutationLoad,fill=state))+
  geom_violin(position=position_dodge(.5))+
  geom_point(position=position_dodge(.5),size = 1,alpha=0.5)+
  theme_bw()+
  facet_grid(hLabel~interaction(population,model))+
  ylab("Genetic Load")+
  xlab("") +
  theme(legend.position = "none")
p1

##try error bars
require(ggplot2)
p2 <- 
  ggplot(allLoads,aes(x=generation,y=L_mutationLoad))+
  geom_point(position=position_dodge(.5),size = .1,alpha=0.5)+
  stat_summary(fun.y = "mean", geom = "point", color = "red", size = 1)+
  stat_summary(fun.y = "mean", geom = "line", color = "black", size = 0.5)+
  theme_bw()+
  facet_grid(hLabel~interaction(allLoads$population,model),scales="free")+
  ylab("Genetic Load")+
  xlab("Generation") +
  theme(legend.position = "none")+
  ggtitle("Genetic Load for 3 Epoch Model, Full Chromosome")+
  geom_vline(xintercept=35,color="blue",size=0.5)

p2
ggsave(paste(plot.dir,"Load.perGen.3Epoch.FullChr",todaysdate,".pdf",sep=""),p2,height=6,width=8)

allLoads<-allLoads[order(allLoads$h,allLoads$generation),]
