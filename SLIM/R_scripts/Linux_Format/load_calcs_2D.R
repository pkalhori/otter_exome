#require(ggplot2)
todaysdate=format(Sys.Date(),format="%Y%m%d")


data.dir="/u/scratch/p/pkalhori/slim/concattedSummaries/"


#pops=c("AK","AL","genericPop.LongerContract")
#models=c("1D.2Epoch.1.5Mb.cds")
#simdates=c(20190424,20190607)
# skipping AL "AL/1D.2Epoch.1.5Mb.cds/20190424/" and CA etc -- add those in next 

popModDates=c("CA_AK/2D.3Epoch.Translocation/20191012")

#popModDates=c("AK/1D.2Epoch.1.5Mb.cds/20190424/","AK/1D.2Epoch.1.5Mb.cds.LongerContract/20190607","genericPop/1D.2Epoch.1.5Mb.cds.20KAncSize/20190611") # AK and AL have dadi parameters, genericPop has parameters based on AK MLE grid that is fur-trade relevant. ### need to come up with better classification system for this. 
#reps=c(seq(1,23))

reps=c(seq(1,25)) # some reps don't make it through Hoffman; so I have a file.exists() test in the loop to skip reps that didn't yield output
hset=c(0,0.5)

allLoads=data.frame()


for (popModDate in popModDates){
print(popModDate)
  for(rep in reps){
	print(rep)
    for(h in hset){
	print(h)
      # check if rep exists (some have random hoffman failures)
      infile=paste(data.dir,popModDate,"/h_",h,"/replicate_",rep,".slim.output.allConcatted.summary.txt.gz",sep="")
print(infile)
      if(file.exists(infile)){
        input <- read.csv(infile,sep=",")
        # equation is, per site: 2hspq + sq^2 is the contribution to load; p = 1-qFreq
        # my "s" is negative, so I want to absolute value s --> |s|
        #### add to dataframe: #####
        # pull out population etc from popModDate
       #loadDF <- data.frame(population=pop)  
        generations<- unique(input$generation)
        subpops=c(1,2)
        
        for (gen in generations){
print(gen)
          for (spop in subpops){
            input$qFreq[input$generation==gen & input$subpop==spop] <- (input$numhet[input$generation==gen & input$subpop==spop] + (2*input$numhom[input$generation==gen & input$subpop==spop])) / (2*input$popsizeDIP[input$generation==gen & input$subpop==spop])
            input$pFreq[input$generation==gen & input$subpop==spop] <- 1 - input$qFreq[input$generation==gen & input$subpop==spop]
            # equation is, per site: 2hspq + sq^2 is the contribution to load; p = 1-qFreq
            # my "s" is negative, so I want to absolute value s --> |s|
            input$loadComponent[input$generation==gen & input$subpop==spop] <- (2*h*abs(input$s[input$generation==gen & input$subpop==spop])*input$qFreq[input$generation==gen & input$subpop==spop]*input$pFreq[input$generation==gen & input$subpop==spop]) + (abs(input$s[input$generation==gen & input$subpop==spop])*((input$qFreq[input$generation==gen & input$subpop==spop])^2))
            # total sites:
            
            S = sum(input$loadComponent[input$generation==gen & input$subpop==spop])
            W = exp(-S) # mean fitness e^-S
            L  = 1 - W # mutation load 
            #### add to dataframe: #####
            # pull out population etc from popModDate
            pop= unlist(lapply(strsplit(popModDate,"/"),"[",1))
            model= unlist(lapply(strsplit(popModDate,"/"),"[",2))
            date= unlist(lapply(strsplit(popModDate,"/"),"[",3))
            loadDF <- data.frame(population=pop)
            loadDF$model <- model
            loadDF$date <- date
            loadDF$rep <- rep
            loadDF$h <- h
            loadDF$S_allsites <- S
            loadDF$W_meanFitness <- W
            loadDF$L_mutationLoad <- L
            loadDF$generation <- gen
            loadDF$subpopulation <- spop
            #### combine with other reps: #####
            allLoads = rbind(allLoads,loadDF)
          }
          
          
        }
      }
    }
  }
}
write.table(allLoads, file="BurnInData2D.txt",quote = F, row.names = F, sep="\t")


