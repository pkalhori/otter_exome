require(ggplot2)
todaysdate=format(Sys.Date(),format="%Y%m%d")


data.dir="C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\"


plot.dir="C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\genetic_load_calcs\\"
dir.create(plot.dir)
#pops=c("AK","AL","genericPop.LongerContract")
#models=c("1D.2Epoch.1.5Mb.cds")
#simdates=c(20190424,20190607)
# skipping AL "AL/1D.2Epoch.1.5Mb.cds/20190424/" and CA etc -- add those in next 

popModDates=c("AK\\1D.3Epoch.LongerBurnIn\\20190918")

#popModDates=c("AK/1D.2Epoch.1.5Mb.cds/20190424/","AK/1D.2Epoch.1.5Mb.cds.LongerContract/20190607","genericPop/1D.2Epoch.1.5Mb.cds.20KAncSize/20190611") # AK and AL have dadi parameters, genericPop has parameters based on AK MLE grid that is fur-trade relevant. ### need to come up with better classification system for this. 
#reps=c(seq(1,23))
reps=c(seq(1,10)) # some reps don't make it through Hoffman; so I have a file.exists() test in the loop to skip reps that didn't yield output
hset=c(0,0.5)
states="BurnIn"


allLoads=data.frame()
popModDates=c("AK\\1D.3Epoch.LongerBurnIn\\20190918")
for (h in hset){
  for (rep in reps){
    infile=paste(data.dir,popModDate,"\\h_",h,"\\replicate_",rep,".slim.output.",state,".allConcatted.summary.txt",sep="")
    input<-read.csv(infile,sep = ",")
    allData=rbind(allData,input)
  }
}
generations<- seq(1000,150000,1000)
for (popModDate in popModDates){
  for(rep in reps){
    for(state in states){
      for(h in hset){
          # check if rep exists (some have random hoffman failures)
        infile=paste(data.dir,popModDate,"\\h_",h,"\\replicate_",rep,".slim.output.",state,".allConcatted.summary.txt",sep="")
        if(file.exists(infile)){
          input <- read.csv(infile,sep=",")
          pop= unlist(lapply(strsplit(popModDate, "[\\\\]|[^[:print:]]"),"[",1))
          model= unlist(lapply(strsplit(popModDate,"[\\\\]|[^[:print:]]"),"[",2))
          date= unlist(lapply(strsplit(popModDate,"[\\\\]|[^[:print:]]"),"[",3))
            # equation is, per site: 2hspq + sq^2 is the contribution to load; p = 1-qFreq
            # my "s" is negative, so I want to absolute value s --> |s|
            #### add to dataframe: #####
            # pull out population etc from popModDate
          loadDF <- data.frame(population=pop)  
          for (gen in generations){
              loadDF$vstr<- sum(input$avgVstrDel[input$generation==gen])
              loadDF$str <- sum(input$avgStrDel[input$generation==gen])
              loadDF$mstr <- sum(input$avgModDel[input$generation==gen])
              loadDF$wstr <- sum(input$avgWkDel[input$generation==gen])
              loadDF$het <- mean(input$meanHet[input$generation==gen])

              
              loadDF$model <- model
              loadDF$date <- date
              loadDF$generation <- gen
              loadDF$rep <- rep
              loadDF$state <- state
              loadDF$h <- h
  
            #### combine with other reps: #####
            allLoads = rbind(allLoads,loadDF)
          }
        }
      }
    }
  }
}


par(mfrow=c(2,2))
##h=0
plot(y=rep1_data_h0$L_mutationLoad,x=rep1_data_h0$generation,type="l")
#####h=0.5
plot(y=rep1_data_h05$L_mutationLoad,x=rep1_data_h05$generation,type="l")

# change order of factors:
allLoads$state <- factor(allLoads$state,levels=c("PreContraction","PostContraction","PostRecovery"))
# label H:

allLoads$hLabel <- paste("h = ",allLoads$h)
allLoads$model<- "1D.3Epoch.LongerBurnIn"
allLoads$population <- factor(allLoads$population,levels=c("AK"))
p1 <- 
  ggplot(allLoads,aes(x=generation,y=L_mutationLoad,fill=state))+
  geom_violin(position=position_dodge(.5))+
  geom_point(position=position_dodge(.5),size = 1,alpha=0.5)+
  theme_bw()+
  facet_grid(hLabel~interaction(date,model))+
  ylab("Genetic Load")+
  xlab("Generation") +
  theme(legend.position = "none")
p1
ggsave(paste(plot.dir,"Load.perPop.PrePostContract.",todaysdate,".pdf",sep=""),p1,height=6,width=8)

p2 <- 
  ggplot(allLoads,aes(x=generation,y=het))+
  geom_point(position=position_dodge(.5),size = .1,alpha=0.5)+
  stat_summary(fun.y = "mean", geom = "point", color = "red", size = 1)+
  stat_summary(fun.y = "mean", geom = "line", color = "black", size = 0.5)+
  theme_bw()+
  facet_grid(hLabel~interaction(allLoads$population,model),scales="free")+
  ylab("Per Site Heterozygosity")+
  xlab("Generation") +
  theme(legend.position = "none")+
  ggtitle("Mean Heterozygosity for AK Burn In")

p2

write.table(allLoads, file="BurnInData3Epoch.txt",quote = F, row.names = F, sep="\t")

p3 <- 
  ggplot(allLoads,aes(x=generation,y=vstr))+
  #geom_point(position=position_dodge(.5),size = .1,alpha=0.5)+
  #stat_summary(fun.y = "mean", geom = "point", color = "black", size = 1)+
  stat_summary(fun.y = "mean", geom = "line", color = "black", size = 0.5)+
  stat_summary(data=allLoads, aes(x=generation,y=str),fun.y = "mean", geom = "line",color = "red", size = 0.5)+
  #stat_summary(data=allLoads, aes(x=generation,y=wstr),fun.y = "mean", geom = "line",color = "blue", size = 0.5)+
  stat_summary(data=allLoads, aes(x=generation,y=mstr),fun.y = "mean", geom = "line",color = "green", size = 0.5)+
  theme_bw()+
  facet_grid(hLabel~interaction(allLoads$population,model),scales="free")+
  ylab("Number of Variants")+
  xlab("Generation") +
  theme(legend.position = "right")+
  ggtitle("Del Muts for AK Burn In")

p3
ggsave(paste(plot.dir,"BurnInHet.",todaysdate,".pdf",sep=""),p2,height=6,width=8)
length(unique(allLoads$rep[allLoads$h==0]))
?tapply()
summary(allLoads)

allLoads$L_mutationLoad[allLoads$rep==1,allLoads$h==0]

write.table(allLoads, file = "Genetic_Load_AK_3Epoch.txt", quote=F, sep="\t", row.names = F)
ggsave(paste(plot.dir,"Load.perGen.3EpochLongerRecovery",todaysdate,".pdf",sep=""),p2,height=6,width=8)


sd_load<-tapply(allLoads$L_mutationLoad, list(allLoads$generation,allLoads$h), FUN = sd)

avg_load<-tapply(allLoads$L_mutationLoad, list(allLoads$generation,allLoads$h), FUN = mean)
aggregate(avg_load_melted,sd_load_melted)
avg_load_melted <- melt(avg_load)
sd_load_melted <- melt(sd_load)

library(dplyr)
avg_load_data<- avg_load_melted %>% rename(Generation=Var1, h_label=Var2, L_mutationLoad=value)
write.table(avg_load_data,file="Average_Load_3Epoch_AK.txt", quote = F, row.names = F, sep = "\t")
load_min<- min(avg_load_data$L_mutationLoad[avg_load_data$h_label==0])

(load_max-load_min)/load_min
