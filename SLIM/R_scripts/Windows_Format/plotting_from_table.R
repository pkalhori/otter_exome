require(ggplot2)
todaysdate=format(Sys.Date(),format="%Y%m%d")


data.dir="C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\"
popModDate=c("CA_AK\\2D.3Epoch.Translocation\\20191103\\")
popModDate2=c("CA_AK\\2D.3Epoch.Translocation\\20191023\\")
popModDate3=c("CA_AK\\2D.3Epoch.Translocation\\20191113\\")

##mig burst
allLoads1<- read.table(paste(data.dir,popModDate,"25for2Gen.txt",sep = ""), header = T)

#cntinual mig
allLoads2<- read.table(paste(data.dir,popModDate2,"migs_data.txt",sep = ""), header = T) 

##no mig
allLoads3 <- read.table(paste(data.dir,popModDate3,"nomig_data.txt",sep = ""), header = T) 

allLoads<- rbind.data.frame(allLoads1,allLoads2,allLoads3)
?write.csv
write.csv(allLoads,file="translocation_data.csv",row.names = FALSE)
getwd()
# label H:

allLoads$hLabel <- paste("h = ",allLoads$h)
allLoads <- allLoads[allLoads$generation>=8000,]
allLoads$subpopulation[allLoads$subpopulation==1] <- "AK"
allLoads$subpopulation[allLoads$subpopulation==2] <- "CA"

allLoads$migLabel <- NA
allLoads[allLoads$model=="2D.3Epoch.Translocation.1perGen",]$migLabel <- "1 ind/gen"
allLoads[allLoads$model=="2D.3Epoch.Translocation.10perGen",]$migLabel <- "10 ind/gen"
allLoads[allLoads$model=="2D.3Epoch.Translocation.5perGen",]$migLabel <- "5 ind/gen"
allLoads[allLoads$model=="2D.3Epoch.Translocation.25perGen",]$migLabel <- "25 ind/gen"
allLoads[allLoads$model=="2D.3Epoch.Translocation.25for2Gen",]$migLabel <- "25 ind for 2 gen"
allLoads[allLoads$model=="2D.3Epoch.NoTranslocation",]$migLabel <- "No Migrants"

AK <- c(8002,8038,8056)
CA <- c(8012,8038,8056)
library(reshape2)
dates <- melt(rbind(AK,CA))
library(tidyverse)
dates <- dates %>% group_by(Var1)

library(tidyverse)

dates_df <- data.frame(
  subpopulation = c("AK","CA"),
  dates = c(8002,8038,8056,8012,8038,8056)
)
allLoads_means_se <- allLoads %>% 
  group_by(generation,hLabel,subpopulation,migLabel) %>% # Group the data by manufacturer
  summarize(mean_load=mean(L_mutationLoad), # Create variable with mean of cty per group
            sd_load=sd(L_mutationLoad), # Create variable with sd of cty per group
            N_load=n(), # Create new variable N of cty per group
            se=sd_load/sqrt(N_load), # Create variable with se of cty per group
            upper_limit=mean_load+sd_load, # Upper limit
            lower_limit=mean_load-sd_load # Lower limit
  ) 


allLoads_means_se<- as.data.frame(allLoads_means_se)
class(allLoads_means_se)

levels(as.factor(allLoads_means_se$migLabel))
allLoads_means_se$migLabel <- factor(allLoads_means_se$migLabel,levels=c("No Migrants", "1 ind/gen" ,"5 ind/gen", "10 ind/gen", "25 ind/gen", "25 ind for 2 gen"))



p2 <- 
  ggplot(allLoads_means_se,aes(x=generation,y=mean_load))+
  geom_line(position=position_dodge(.5),size = 1,alpha=0.5,color="blue")+
  #stat_summary(fun.y = "mean", geom = "point", size = 1, color=model)+
  #stat_summary(fun.y = "mean", geom = "line", size = 0.5)+
  theme_bw()+
  
  ylab("Genetic Load")+
  xlab("Generation") +
  theme(legend.position = "left")+
  ggtitle("Genetic Load for 3 Epoch Model")+
  facet_grid(hLabel~subpopulation,scales="free")+

  #stat_summary(fun.data = allLoads_means_se, geom = "errorbar")+
  geom_errorbar(data=allLoads_means_se, mapping=aes(ymin=lower_limit,ymax=upper_limit),color="green",size=0.1)+
  geom_vline(data=dates_df, aes(xintercept=dates))
  
  #theme(legend.title=element_blank())


p3 <- 
  ggplot(allLoads,aes(x=generation,y=L_mutationLoad))+
  #geom_line(position=position_dodge(.5),size = .1,alpha=0.5)+
  #stat_summary(fun.y = "mean", geom = "point", size = 1, color=model)+
  stat_summary(fun.y = "mean", geom = "line", size = 0.5)+
  theme_bw()+
  geom_vline(data=dates, aes(xintercept=value))+
  facet_grid(hLabel~interaction(allLoads$subpopulation),scales="free")+
  ylab("Genetic Load")+
  xlab("Generation") +
  theme(legend.position = "left")+
  ggtitle("Genetic Load for 3 Epoch Model")
  
  
  #stat_summary(fun.data = allLoads_means_se, geom = "errorbar")+
  #geom_errorbar(data=allLoads_means_se, mapping=aes(ymin=lower_limit,ymax=upper_limit))+
  #theme(legend.title=element_blank())
              
class(allLoads$subpopulation)
dates<-as.data.frame(dates)
as.tibble(dates)



  




plot.dir="C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\genetic_load_calcs\\"
ggsave(paste(plot.dir,"3Epoch_BothPops",todaysdate,".pdf",sep=""),p2,height=6,width=8)





library(dplyr)
####Pre and Post contraction boxplots withhout translocation
allLoads$hLabel <- paste("h = ",allLoads$h)
allLoads <- allLoads[allLoads$generation>=8000,]
allLoads$subpopulation[allLoads$subpopulation==1] <- "AK"
allLoads$subpopulation[allLoads$subpopulation==2] <- "CA"
boxplot_loads <- allLoads[allLoads$generation%in%c(8002,8012,8038),]


boxplot_loads$state<-NA
boxplot_loads<-filter(boxplot_loads,(generation==8002 & subpopulation=="AK")|(generation==8012 & subpopulation=="CA")|(generation==8038))
boxplot_loads[(boxplot_loads$generation==8002 & boxplot_loads$subpopulation=="AK"),]$state <- "PreContraction"
boxplot_loads[(boxplot_loads$generation==8012 & boxplot_loads$subpopulation=="CA"),]$state <- "PreContraction"
boxplot_loads[boxplot_loads$generation==8038,]$state <- "PostContraction"
boxplot_loads$miglabel<- NA
boxplot_loads$state <- factor(boxplot_loads$state,levels=c("PreContraction","PostContraction"))
p1 <- 
  ggplot(boxplot_loads,aes(x=state,y=L_mutationLoad))+
  facet_grid(hLabel~interaction(boxplot_loads$subpopulation))+
  geom_boxplot()+

  #theme_bw()+


  ylab("Genetic Load")+
  xlab("State") +
  #theme(legend.position = "left")+
  ggtitle("Genetic Load Before and After Contraction")
ggsave(paste(plot.dir,"Boxplot_BothPops",todaysdate,".pdf",sep=""),p1,height=6,width=8)

AK_means <- allLoads_means_se[allLoads_means_se$subpopulation=="AK",]
max(CA_means$mean_load)