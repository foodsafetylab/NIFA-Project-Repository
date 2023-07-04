library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggsci)

###NIFA Poster 2023 - Figures for papers too - Must use Rdata = NIFA_IAFP_posterVisuals.Rdata - located C:\Users\<user>\Box Sync\NIFA Project\NIFA Project Repository
master<-read.csv("C:/Users/jfq/Box Sync/NIFA Project/CA 9-2021/Plate results/Master plate counts.csv")
#Colorblind Friendly palette
cbp1 <- c( "#CC79A7","#56B4E9", "#009E73",
           "#0072B2", "#D55E00")

#Preharvest
Overall_PrHarvest_comparison<-ggplot(data=PreHarvest, aes(x=Sample.Type, y=Best.Estimate, col=Test))+
  geom_boxplot(outlier.colour="black", outlier.shape = 19)+
  geom_hline(yintercept = 0.47, color="red")+
  geom_hline(yintercept = 0.6, color="steelblue")+
  scale_color_aaas()+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=0.7, size=1)+
  ggtitle("Overall PreHarvest Bacterial Recovery Comparison - Swab v Composite Produce Samples ")

Overall_PrHarvest_comparison 

#InHarvest
Harvesters<-InHarvest_noEcoli%>%filter(Location=="Leftover"|Location=="Gloves")

Overall_OnHarvester_comparison<-ggplot(data=Harvesters, aes(x=Sample.Type, y=Best.Estimate.log.CFU.g., col=Sample.Type))+
  geom_boxplot(outlier.colour="black", outlier.shape = 19)+
  #geom_hline(yintercept = 0.47, color="red")+
  #geom_hline(yintercept = 0.6, color="steelblue")+
  scale_color_aaas()+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=0.7, size=1)+
  ggtitle("Overall In-Harvest Bacterial Recovery Comparison - Aggregative Swabs, Gloves and Composite Produce Samples")

Overall_OnHarvester_comparison+facet_grid(.~Test)

Overall_BinHarvester_comparison<-ggplot(data=InHarvest_noEcoli, aes(x=Sample.Type, y=Best.Estimate.log.CFU.g., col=Test))+
  geom_boxplot(outlier.colour="black", outlier.shape = 19)+
  #geom_hline(yintercept = 0.47, color="red")+
  #geom_hline(yintercept = 0.6, color="steelblue")+
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=0.7, size=1)+
  ggtitle("Overall In-Harvest Bacterial Recovery Comparison - Aggregative Swabs, Gloves and Composite Produce Samples")

Overall_BinHarvester_comparison+facet_wrap(.~Test)

#remove the ecoli sample from master - call it, master no ecoli
master_noecoli<-master%>%filter(Test=="APC"|Test=="Coliform")

master_noecoli$Sample.Time.Point_factor<- factor(master_noecoli$Sample.Time.Point, levels=c('Preharvest', 'Harvest', 'Post Harvest'))
OvSample<-ggplot(data= master_noecoli, aes(x=Sample.Type, y=Log.CFU.g., col= Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=15),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  stat_summary(aes(x=Sample.Type, y=Log.CFU.g.), fun="mean", shape=4) +
  facet_grid(Test~Sample.Time.Point, scales = "free_x") + 
  ggtitle("Overall Bacterial Recovery Comparison")
OvSample

##Renaming of master Sample Description
master_noecoli$Sample.Description[master_noecoli$Sample.Description=="Cloth produce bottom" & master_noecoli$Sample.Time.Point == "Preharvest"] <- "Aggregative Swab"
master_noecoli$Sample.Description[master_noecoli$Sample.Description=="Cloth produce sides" & master_noecoli$Sample.Time.Point == "Preharvest"] <- "Aggregative Swab"
master_noecoli$Sample.Description[master_noecoli$Sample.Description=="Produce bottom" & master_noecoli$Sample.Time.Point == "Preharvest"] <- "Composite Produce"
master_noecoli$Sample.Description[master_noecoli$Sample.Description=="Produce sides" & master_noecoli$Sample.Time.Point == "Preharvest"] <- "Composite Produce"
master_noecoli$Sample.Description[master_noecoli$Sample.Description=="Produce sides HR" & master_noecoli$Sample.Time.Point == "Preharvest"] <- "Produce High Resolution"


master_noecoli$Sample.Description[master_noecoli$Sample.Description=="Chute" & master_noecoli$Sample.Time.Point == "Harvest"] <- "Swab on Chute"
master_noecoli$Sample.Description[master_noecoli$Sample.Description=="Cloth bin tops" & master_noecoli$Sample.Time.Point == "Harvest"] <- "Swab of Bin Tops"
master_noecoli$Sample.Description[master_noecoli$Sample.Description=="Produce bins" & master_noecoli$Sample.Time.Point == "Harvest"] <- "Produce from Bins"
master_noecoli$Sample.Description[master_noecoli$Sample.Description=="Produce leftover trims" & master_noecoli$Sample.Time.Point == "Harvest"] <- "Produce Leftover Trims"


master_noecoli$Sample.Description[master_noecoli$Sample.Description=="Produce sides" & master_noecoli$Sample.Time.Point == "Post Harvest"] <- "Produce Exterior"
master_noecoli$Sample.Description[master_noecoli$Sample.Description=="Chute" & master_noecoli$Sample.Time.Point == "Post Harvest"] <- "Swab on Chute"
master_noecoli$Sample.Description[master_noecoli$Sample.Description=="Cloth produce interior" & master_noecoli$Sample.Time.Point == "Post Harvest"] <- "Swab of Produce Interior"
master_noecoli$Sample.Description[master_noecoli$Sample.Description=="Produce interior" & master_noecoli$Sample.Time.Point == "Post Harvest"] <- "Produce Interior"

master_noecoli$Sample.Type[master_noecoli$Sample.Type=="Cloth"] <- "Aggregative Swab"


master_noecoli$Sample.Description_factor<- factor(master_noecoli$Sample.Description, levels=c('Aggregative Swab', 'Composite Produce', 'Produce High Resolution', 
                                                                                              'Glove', 'Produce Leftover Trims', 'Swab on Chute', 'Swab of Bin Tops',  "Produce from Bins", 
                                                                                              'Produce Exterior', 'Swab of Produce Interior', 'Produce Interior'))

all_plot<-ggplot(data= master_noecoli, aes(x=Sample.Description_factor, y=Log.CFU.g., col=Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1, size = 10), legend.position = "top")+
  expand_limits(y=0)+
  scale_color_brewer(palette="Dark2")+
  labs(col="Sample Type")+
  #scale_color_manual(values=cbp1)+
  scale_y_continuous("Bacterial Counts Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1, size=0.5)+
  stat_summary(aes(x=Sample.Description_factor, y=Log.CFU.g.), fun="mean", shape=4) +
  ggtitle("Bacterial recovery of Aggregative and Composite Produce Samples Collected at Different Stages")+
  xlab("Sample Description") +
  facet_grid(Test~Sample.Time.Point_factor, scales="free", switch="y")
all_plot
ggsave("Samples 2021_Master plot.jpeg",all_plot, "jpeg",width = 10, height = 6, units = "in", dpi=300, path = "C:/Users/jorge/Box Sync/NIFA Project/NIFA Project Repository/Outputs 2021")
####END OF RESULTS FROM 2021 FOR NIFA POSTER 

#Sampling technique characterization - Prehydrated and Dry cloths collected by hand or using MSD
OvPH_MT_pressure_hyd<-ggplot(data= Only_Swabs_noEcoli, aes(x=Hand.msd, y=Best..Estimate, col= Wet.dry))+
  geom_boxplot(outlier.alpha = 0, position = position_dodge2(preserve="single"))+
  theme(plot.title = element_text(size = rel(1.2), face = "bold"), axis.text = element_text(size= 10, angle=0, hjust=0.5), 
  axis.title = element_text(size=12))+
  geom_point(position=position_jitterdodge(jitter.width = 0.05, jitter.height = 0.05))+
  #stat_compare_means(method = "anova", label.y = 40)+
  expand_limits(y=0)+
  scale_y_continuous("Bacterial Counts Estimate (Log CFU/g of swab)", breaks= c(1,2,3,4,5,6,7,8,9))+
  scale_x_discrete("Collection method (Pressure)", labels=c("Hand", "Manual Sampling Device"))+
  scale_color_manual(values=c("#1B9E77", "#d95f02"), labels=c("Dry Cloth", "Prehydrated Cloth")) +
  #scale_color_discrete(labels=c("Dry Cloth", "Prehydrated Cloth"))+
  labs(col="Hydration Factor")+
  ggtitle("Pressure and Hydration Effects on Bacterial Recovery in Aggregative Swab Sampling Technique Optimization")

OvPH_MT_pressure_hyd+facet_grid(.~ Test)+stat_compare_means(aes(group=Wet.dry),label="p.format")+stat_compare_means(aes(group=Hand.msd),label="p.format")

ggsave("Swab-characterization-overall2.tiff",device = "tiff", dpi=300, path = "C:/Users/jfq/Box Sync/NIFA Project/NIFA Project Repository/Output 2022", height= 5, width= 10, units = "in")
ggsave("Swab-characterization-overall2_3.jpeg", device = "jpeg", dpi=300, path = "C:/Users/jfq/Box Sync/NIFA Project/NIFA Project Repository/Output 2022", height= 5, width= 10, units = "in")

#Developed sampling technique compared to composite produce samples
OvPH_MTvPS<-ggplot(data= Everything_wet, aes(x=Sample.Type, y=Best..Estimate, col=Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(plot.title = element_text(size = rel(1.2), face = "bold"), legend.position = "none", axis.text = element_text(size= 10, angle=0, hjust=0.5), axis.title = element_text(size=12))+
  geom_point(alpha=0.5)+
  expand_limits(y=0)+
  scale_color_manual(values=c("#1B9E77", "#7570B3")) +
  scale_y_continuous("Bacterial Counts Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7,8,9))+
  scale_x_discrete("Sample Type")+
  labs(fill= "Sample Type")+
  ggtitle("Bacterial Recovery in Aggregative Swab Cloths and Composite Produce Samples")
OvPH_MTvPS+facet_grid(.~ Test)

ggsave("Swab-v-Produce3.tiff",device = "tiff", dpi=300, path = "C:/Users/jfq/Box Sync/NIFA Project/NIFA Project Repository/Output 2022", height= 4, width= 6, units = "in")
ggsave("Swab-v-Produce3_2.jpeg",device = "jpeg", dpi=300, path = "C:/Users/jorge/Box Sync/NIFA Project/NIFA Project Repository/Output 2022", height= 5, width= 10, units = "in")


   df <- read.table(text = 
                   "Sample-type  Generic-E.coli-Positive  Generic-E.coli-Negative  Total
Aggregative-Swab  19   6 25
Produce-Sample  15   11 26
Total  34 17 51"
)
names(df) <- c("Sample Type","Generic E.coli Positive","Generic E. coli Negative", "Total")

EqDatedf <- as.data.frame(df[1,])
#EmptyLine <- data.frame(Sample Type = "", Generic E.coli Positive = "", Generic E. coli Negative = "", Total = "")

pdf(file = "q.pdf")

for (i in 2:nrow(df)) 
{
  if (as.vector(df$Date[i])  ==  as.vector(df$Date[i-1])) 
  {EqDatedf <- rbind(EqDatedf, df[i,])}
  
  else {
    EqDatedf <- rbind(EqDatedf, EmptyLine)
    EqDatedf <- rbind(EqDatedf, df[i,]) 
  }
}

grid.table(EqDatedf, show.rownames = FALSE)
dev.off()


