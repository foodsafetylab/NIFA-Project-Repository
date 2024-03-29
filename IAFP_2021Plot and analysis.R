library(tidyverse)
library(ggplot2)


PreHarvest<-read.csv("PreHarvest Trimmed.csv")
View(PreHarvest)

PreHarvest$Sample.Type<-as.factor(PreHarvest$Sample.Type)
PreHarvest$Area<-as.factor(PreHarvest$Area)
PreHarvest$Test<-as.factor(PreHarvest$Test)

PreHarvest$Sample_Type<- ifelse(PreHarvest$Sample.Type == "Swab", "Aggregative Swab", 
                                ifelse(PreHarvest$Sample.Type == "Grab", "Composite Produce Sample", "High Resolution Produce Sample"))
PreHarvest$Areas<-ifelse(PreHarvest$Area == "Edge", "Road Edge",
                         ifelse(PreHarvest$Area == "Bed 1", "Outer Edge", "Center"))
PreHarvest$Tests<-ifelse(PreHarvest$Test == "APC", "Aerobic Plate Counts", "Coliform Counts")

PreHarvest$Sample_Type<-as.factor(PreHarvest$Sample_Type)
PreHarvest$Areas<-as.factor(PreHarvest$Areas)
PreHarvest$Tests<-as.factor(PreHarvest$Tests)

Pr_ASvCPSvHR_Loc<-ggplot(data= PreHarvest, aes(x=Areas, y=Best.Estimate , col= Sample_Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=10),axis.title=element_text(size=10,face="bold"),axis.text.x = element_text(size=8,angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_colour_discrete(labels = function(x) str_wrap(x, width = 15))+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15))+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("PreHarvest Bacterial Recovery Comparison between Field Locations")
Pr_ASvCPSvHR_Loc+ facet_grid(.~Test) +labs(col="Sample Type") 

ggsave("Swab-Produce-Loc_2021.tiff",device = "tiff", dpi=300, path = "C:/Users/jorge/Box Sync/NIFA Project/NIFA Project Repository/Outputs 2021", height= 4, width= 8, units = "in")
ggsave("Swab-Produce-Loc_2021.jpeg", device = "jpeg", dpi=300, path = "C:/Users/jorge/Box Sync/NIFA Project/NIFA Project Repository/Outputs 2021", height= 4, width= 8, units = "in")


Pr_ASvCPSvHR<-ggplot(data= PreHarvest, aes(x=Sample_Type, y=Best.Estimate , col= Test))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=10),axis.title=element_text(size=10,face="bold"),axis.text.x = element_text(size=8,angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15))+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("Overall PreHarvest Bacterial Recovery Comparison between Sample Types")
Pr_ASvCPSvHR+ facet_grid(.~Test) +labs(x= "Sample Type", col="Microbiological Test")

ggsave("Swab-Produce_2021.tiff",device = "tiff", dpi=300, path = "C:/Users/jorge/Box Sync/NIFA Project/NIFA Project Repository/Outputs 2021", height= 4, width= 8, units = "in")
ggsave("Swab-Produce_2021.jpeg", device = "jpeg", dpi=300, path = "C:/Users/jorge/Box Sync/NIFA Project/NIFA Project Repository/Outputs 2021", height= 4, width= 8, units = "in")

Preharvest_noHR<-PreHarvest%>%filter(Sample_Type!="High Resolution Produce Sample")

Pr_ASvCPS<-ggplot(data= Preharvest_noHR, aes(x=Sample_Type, y=Best.Estimate , col= Test))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=10),axis.title=element_text(size=10,face="bold"),axis.text.x = element_text(size=8,angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15))+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("Overall PreHarvest Bacterial Recovery Comparison between Sample Types")
Pr_ASvCPS+ facet_grid(.~Test) +labs(x= "Sample Type", col="Microbiological Test")

ggsave("Swab-Produce-noHR_2021.tiff",device = "tiff", dpi=300, path = "C:/Users/jorge/Box Sync/NIFA Project/NIFA Project Repository/Outputs 2021", height= 4, width= 8, units = "in")
ggsave("Swab-Produce-noHR_2021.jpeg", device = "jpeg", dpi=300, path = "C:/Users/jorge/Box Sync/NIFA Project/NIFA Project Repository/Outputs 2021", height= 4, width= 8, units = "in")

PreHarvestAPC<-PreHarvest%>%filter(Test=="APC")

PreHarvest%>%
  count(Sample.Type, Test, sort=T)


