require(ggplot2)
todaysdate=format(Sys.Date(),format="%Y%m%d")


data.dir="C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\"


plot.dir="C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\genetic_load_calcs\\"
dir.create(plot.dir)
#pops=c("AK","AL","genericPop.LongerContract")
#models=c("1D.2Epoch.1.5Mb.cds")
#simdates=c(20190424,20190607)
# skipping AL "AL/1D.2Epoch.1.5Mb.cds/20190424/" and CA etc -- add those in next 

popModDates=c("CA_AK\\2D.3Epoch.Translocation\\20191003")

#popModDates=c("AK/1D.2Epoch.1.5Mb.cds/20190424/","AK/1D.2Epoch.1.5Mb.cds.LongerContract/20190607","genericPop/1D.2Epoch.1.5Mb.cds.20KAncSize/20190611") # AK and AL have dadi parameters, genericPop has parameters based on AK MLE grid that is fur-trade relevant. ### need to come up with better classification system for this. 
#reps=c(seq(1,23))
reps=1
reps=c(seq(1,10)) # some reps don't make it through Hoffman; so I have a file.exists() test in the loop to skip reps that didn't yield output
hset=c(0,0.5)
unique(input$chunk)
tail(input)
allLoads=data.frame()
popModDates=c("CA_AK\\2D.3Epoch.Translocation\\20191003")

for (h in hset){
  for (rep in reps){
    infile=paste(data.dir,popModDate,"\\h_",h,"\\replicate_",rep,".slim.output.",".allConcatted.summary.txt.gz",sep="")
    input<-read.csv(infile,sep = ",")
    allData=rbind(allData,input)
  }
}
#generations<- c(seq(100,8000,100),seq(8000,8000)
for (popModDate in popModDates){
  for(rep in reps){
    for(h in hset){
        # check if rep exists (some have random hoffman failures)
        infile=paste(data.dir,popModDate,"\\h_",h,"\\replicate_",rep,".slim.output.allConcatted.summary.txt.gz",sep="")
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
          generations<- unique(input$generation)
          subpops=c(1,2)
          
          for (gen in generations){
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
            pop= unlist(lapply(strsplit(popModDate, "[\\\\]|[^[:print:]]"),"[",1))
            model= unlist(lapply(strsplit(popModDate,"[\\\\]|[^[:print:]]"),"[",2))
            date= unlist(lapply(strsplit(popModDate,"[\\\\]|[^[:print:]]"),"[",3))
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

read.table(infile)
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
allLoads$subpopulation <- factor(allLoads$subpopulation,levels=c(1,2))
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

allLoads <- allLoads[allLoads$h==0,]
allLoads[allLoads$subpopulation==1,]
subpop1<- allLoads[allLoads$subpopulation==1,]
subpop1 <- subpop1[subpop1$generation>=8000,]
subpop2<- allLoads[allLoads$subpopulation==1,]
subpop2 <- subpop2[subpop2$generation>=8000,]
p2 <- 
  ggplot(subpop1,aes(x=generation,y=L_mutationLoad))+
  geom_point(position=position_dodge(.5),size = .1,alpha=0.5)+
  stat_summary(fun.y = "mean", geom = "point", color = "red", size = 1)+
  stat_summary(fun.y = "mean", geom = "line", color = "black", size = 0.5)+
  theme_bw()+
  facet_grid(hLabel~interaction(subpop1$population,model),scales="free")+
  ylab("Genetic Load")+
  xlab("Generation") +

  theme(legend.position = "none")+
  ggtitle("2D Load AK")

p2

write.table(allLoads, file="BurnInData2D.txt",quote = F, row.names = F, sep="\t")

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
