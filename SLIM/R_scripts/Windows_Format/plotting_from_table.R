require(ggplot2)
todaysdate=format(Sys.Date(),format="%Y%m%d")


data.dir="C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\"
popModDate=c("CA_AK\\2D.3Epoch.Translocation\\20191012\\")

allLoads<- read.table(paste(data.dir,popModDate,"BurnInData2D.txt",sep = ""), header = T)

pop_1 <- allLoads[allLoads$subpopulation==1,]
pop_2 <- pop_2[pop_2$generation>8000,]
pop_2 <- allLoads[allLoads$subpopulation==2,]
# label H:

pop_2$hLabel <- paste("h = ",pop_2$h)
pop_2$model<- "2D.3Epoch.Translocation"
pop_2$population <- factor(pop_2$population,levels=c("CA_AK"))
p2 <- 
  ggplot(pop_2,aes(x=generation,y=L_mutationLoad))+
  #geom_point(position=position_dodge(.5),size = .1,alpha=0.5)+
  stat_summary(fun.y = "mean", geom = "point", color = "red", size = 1)+
  stat_summary(fun.y = "mean", geom = "line", color = "black", size = 0.5)+
  theme_bw()+
  facet_grid(hLabel~interaction(pop_1$model),scales="free")+
  ylab("Genetic Load")+
  xlab("Generation") +
  theme(legend.position = "none")+
  ggtitle("Genetic Load for 3 Epoch Model")+
  geom_vline(xintercept=8038,color="blue",size=0.5)+
  geom_vline(xintercept=8056,color="blue",size=0.5)+
  geom_vline(xintercept=8012,color="blue",size=0.5)

plot.dir="C:\\Users\\poone\\OneDrive\\Documents\\Otter_Exome_Project\\SLIM_results\\genetic_load_calcs\\"
ggsave(paste(plot.dir,"Load.2D.CA",todaysdate,".pdf",sep=""),p2,height=6,width=8)
max(allLoads$L_mutationLoad[allLoads$h==0])
min(allLoads$L_mutationLoad[allLoads$h==0])
mean(allLoads$L_mutationLoad[allLoads$h==0 & allLoads$generation==96])
mean(allLoads$L_mutationLoad[allLoads$h==0 & allLoads$generation==0])
mean(allLoads$L_mutationLoad[allLoads$h==0 & allLoads$generation==28])
