require(ggplot2)
todaysdate=format(Sys.Date(),format="%Y%m%d")


data.dir="C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\"
popModDate2=c("CA\\1D.3Epoch.LongerRecovery\\20191012\\")
popModDate=c("AK\\1D.5Epoch\\20200814\\")

##AK
allLoads_old<- read.table(paste(data.dir,popModDate2,"CA_data.txt",sep = ""), header = T)

allLoads<- read.table(paste(data.dir,popModDate,"20200814LoadPerGeneration.ThroughTime.AllReps.RemovedBurninFixedVar.txt",sep = ""), header = T)



allLoads<- rbind.data.frame(allLoads_AK)


# label H:

allLoads$hLabel[allLoads$h==0] <- paste("Recessive")
allLoads$hLabel[allLoads$h==0.5] <- paste("Additive")


AK <- c(36)
CA <- c(36)
library(reshape2)
dates <- melt(rbind(AK,CA))
library(tidyverse)
dates <- dates %>% group_by(Var1)
install.packages("tidyselect")
install.packages("tidyverse")
library(tidyverse)

dates_df <- data.frame(
  population = c("AK","CA"),
  years = c(36,36)
)
allLoads_means_se <- allLoads %>% 
  group_by(generation) %>% # Group the data by manufacturer
  summarize(mean_load=mean(L), # Create variable with mean of cty per group
            sd_load=sd(L), # Create variable with sd of cty per group
            N_load=n(), # Create new variable N of cty per group
            se=sd_load/sqrt(N_load), # Create variable with se of cty per group
            upper_limit=mean_load+sd_load, # Upper limit
            lower_limit=mean_load-sd_load # Lower limit
  ) 
allLoads_means_se_old <- allLoads_old %>% 
  group_by(generation,h,population) %>% # Group the data by manufacturer
  summarize(mean_load=mean(L_mutationLoad), # Create variable with mean of cty per group
            sd_load=sd(L_mutationLoad), # Create variable with sd of cty per group
            N_load=n(), # Create new variable N of cty per group
            se=sd_load/sqrt(N_load), # Create variable with se of cty per group
            upper_limit=mean_load+sd_load, # Upper limit
            lower_limit=mean_load-sd_load # Lower limit
  ) 


allLoads_means_se<- as.data.frame(allLoads_means_se)
class(allLoads_means_se)

##Recessive

Pre_Contraction <-allLoads_means_se$mean_load[allLoads_means_se$generation==0 & allLoads_means_se$hLabel=="Recessive"]
Post_Contraction <- allLoads_means_se$mean_load[allLoads_means_se$generation==36 & allLoads_means_se$hLabel=="Recessive"]
(Post_Contraction-Pre_Contraction)/Pre_Contraction


##Additive
Pre_Contraction <-allLoads_means_se$mean_load[allLoads_means_se$generation==0 & allLoads_means_se$hLabel=="Additive"]
Post_Contraction <- allLoads_means_se$mean_load[allLoads_means_se$generation==36 & allLoads_means_se$hLabel=="Additive"]
(Post_Contraction-Pre_Contraction)/Pre_Contraction

####not removing fixed burn in load
##Recessive

Pre_Contraction <-allLoads_means_se$mean_load[allLoads_means_se$generation==0 & allLoads_means_se$hLabel=="Recessive"]
Post_Contraction <- allLoads_means_se$mean_load[allLoads_means_se$generation==36 & allLoads_means_se$hLabel=="Recessive"]
(Post_Contraction-Pre_Contraction)/Pre_Contraction


##Additive
Pre_Contraction <-allLoads_means_se_old$mean_load[allLoads_means_se_old$generation==0 & allLoads_means_se_old$h==0]
Post_Contraction <- allLoads_means_se_old$mean_load[allLoads_means_se_old$generation==26 & allLoads_means_se_old$h==0]
max(allLoads_means_se_old$mean_load[allLoads_means_se_old$h==0])
(Post_Contraction-Pre_Contraction)/Pre_Contraction

Pre_Contraction/Post_Contraction



p2 <- 
  ggplot(allLoads_means_se,aes(x=generation,y=mean_load))+
  geom_line(position=position_dodge(.5),size = 1,alpha=0.5,color="blue")+
  #stat_summary(fun.y = "mean", geom = "point", size = 1, color=model)+
  #stat_summary(fun.y = "mean", geom = "line", size = 0.5)+
  theme_bw()+
  
  ylab("Genetic Load")+
  xlab("Generation") +
  theme(legend.position = "left")+
  ggtitle("Genetic Load for AK 5Epoch")+
  #facet_grid(hLabel~population,scales="free")+
  
  #stat_summary(fun.data = allLoads_means_se, geom = "errorbar")+
  geom_errorbar(data=allLoads_means_se, mapping=aes(ymin=lower_limit,ymax=upper_limit),color="green",size=0.1)+
  geom_vline(xintercept=36)
  #geom_vline(xintercept=50)+
  #geom_vline(xintercept=56)

plot.dir<-"C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\AK\\1D.5Epoch\\20200814\\"
ggsave(paste(plot.dir,"5Epoch_AK_hs_",todaysdate,".pdf",sep=""),p2,height=6,width=8)
#theme(legend.title=element_blank())

head(allLoads)
p3 <- 
  ggplot(allLoads,aes(x=generation,y=L))+
  #geom_line(position=position_dodge(.5),size = .1,alpha=0.5)+
  #stat_summary(fun.y = "mean", geom = "point", size = 1, color=model)+
  stat_summary(fun = "mean", geom = "line", size = 0.5)+
  theme_bw()+
  #geom_vline(data=dates, aes(xintercept=value))+
  #facet_grid(hLabel~interaction(allLoads$subpopulation),scales="free")+
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
ggsave(paste(plot.dir,"5Epoch_AK",todaysdate,".pdf",sep=""),p2,height=6,width=8)





library(dplyr)
####Pre and Post contraction boxplots withhout translocation


#boxplot_loads <- allLoads[allLoads$generation%in%c(8002,8012,8038),]



boxplot_loads<-filter(allLoads,(generation==36 & population=="AK")|(generation==26 & population=="CA")|(generation==2))

boxplot_loads$state<-NA
boxplot_loads[(boxplot_loads$generation==36 & boxplot_loads$population=="AK"),]$state <- "PostContraction"
boxplot_loads[(boxplot_loads$generation==26 & boxplot_loads$population=="CA"),]$state <- "PostContraction"
boxplot_loads[boxplot_loads$generation==2,]$state <- "PreContraction"

boxplot_loads$state <- factor(boxplot_loads$state,levels=c("PreContraction","PostContraction"))
p1 <- 
  ggplot(boxplot_loads,aes(x=state,y=L_mutationLoad))+
  facet_grid(hLabel~interaction(boxplot_loads$population),scales="free")+
  geom_boxplot()+
  
  theme_bw()+
  
  
  ylab("Genetic Load")+
  xlab("State") +
  #theme(legend.position = "left")+
  ggtitle("Genetic Load Before and After Contraction")
ggsave(paste(plot.dir,"Boxplot_SeparatePops",todaysdate,".pdf",sep=""),p1,height=6,width=8)

AK_means <- allLoads_means_se[allLoads_means_se$subpopulation=="AK",]
max(CA_means$mean_load)
