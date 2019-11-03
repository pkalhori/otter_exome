require(ggplot2)
todaysdate=format(Sys.Date(),format="%Y%m%d")


data.dir="C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\"
popModDate=c("CA_AK\\2D.3Epoch.Translocation\\20191023\\")

allLoads<- read.table(paste(data.dir,popModDate,"migs_data.txt",sep = ""), header = T)


# label H:

allLoads$hLabel <- paste("h = ",allLoads$h)
allLoads <- allLoads[allLoads$generation>=8000,]
allLoads$subpopulation[allLoads$subpopulation==1] <- "AK"
allLoads$subpopulation[allLoads$subpopulation==2] <- "CA"



p2 <- 
  ggplot(allLoads,aes(x=generation,y=L_mutationLoad, colour=migLabel))+
  #geom_point(position=position_dodge(.5),size = .1,alpha=0.5)+
  #stat_summary(fun.y = "mean", geom = "point", size = 1, color=model)+
  stat_summary(fun.y = "mean", geom = "line", size = 0.5)+
  theme_bw()+
  facet_grid(hLabel~interaction(allLoads$subpopulation),scales="free")+
  ylab("Genetic Load")+
  xlab("Generation") +
  theme(legend.position = "left")+
  ggtitle("Genetic Load for 3 Epoch Model With Translocation")+
  geom_vline(data=dates, aes(xintercept=value))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar")+
  theme(legend.title=element_blank())
              
p2head(allLoads)
allLoads$migLabel <- NA
allLoads[allLoads$model=="2D.3Epoch.Translocation.1perGen",]$migLabel <- "1 ind/gen"
allLoads[allLoads$model=="2D.3Epoch.Translocation.10perGen",]$migLabel <- "10 ind/gen"
allLoads[allLoads$model=="2D.3Epoch.Translocation.5perGen",]$migLabel <- "5 ind/gen"
allLoads[allLoads$model=="2D.3Epoch.Translocation.25perGen",]$migLabel <- "25 ind/gen"
dates$
AK <- c(8002,8038,8056)
CA <- c(8012,8038,8056)
library(reshape2)
dates <- melt(rbind(AK,CA))
library(tidyverse)
dates <- dates %>% group_by(Var1)

plot.dir="C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\genetic_load_calcs\\"
ggsave(paste(plot.dir,"Multiple_Migs",todaysdate,".pdf",sep=""),p2,height=6,width=8)

library(reshape2)




max(allLoads$L_mutationLoad[allLoads$h==0])
min(allLoads$L_mutationLoad[allLoads$h==0])
mean(allLoads$L_mutationLoad[allLoads$h==0 & allLoads$generation==96])
mean(allLoads$L_mutationLoad[allLoads$h==0 & allLoads$generation==0])
mean(allLoads$L_mutationLoad[allLoads$h==0 & allLoads$generation==28])
