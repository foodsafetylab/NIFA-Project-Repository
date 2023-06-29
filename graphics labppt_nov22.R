library(tidyverse)

LabPPT_AdjustedFanoe_Exp_1 <- read.csv("LabPPT_AdjustedFanoe Exp_1.csv", header = T)
Only_Swabs<-LabPPT_AdjustedFanoe_Exp_1%>% filter(Sample.Type=="Aggregative Swab")
Only_Swabs_noEcoli<-Only_Swabs%>%filter((Test != "E.coli"))
Everything_but_Ecoli<-LabPPT_AdjustedFanoe_Exp_1%>% filter(Test != "E.coli")
Everything_wet<-Everything_but_Ecoli%>%filter(Sample.Type=="Produce Sample"|Sample.Type=="Aggregative Swab" & Wet.dry!="DRY")

labels<-c("Dry", "PreHydrated") 

#summary statistics
Only_Swabs%>%group_by(Test, Wet.dry)%>%summarise(mean=mean(Best..Estimate), sd=sd(Best..Estimate))
Everything_wet%>%group_by(Test, Sample.Type)%>%summarise(mean=mean(Best..Estimate), sd=sd(Best..Estimate))

#only hydrated v dry
OvPH_MT<-ggplot(data= Only_Swabs_noEcoli, aes(x=Wet.dry, y=Best..Estimate, col=Wet.dry))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=10),axis.text.x = element_text(angle=0, hjust=0.5))+
  geom_point(position=position_jitter(width = 0.1, height = 0.1))+
  expand_limits(y=0)+
  scale_color_brewer(palette="Dark2")+
  scale_y_continuous("Best Estimate (Log CFU/g of swab)", breaks= c(1,2,3,4,5,6,7,8,9))+
  scale_x_discrete("Sample Description")+
  labs(col="Hydration Factor", labels=c("Dry", "PreHydrated"))+
  ggtitle("Bacterial Recovery of Dry and Prehydrated Aggregative Swabs")
OvPH_MT+facet_grid(.~ Test)

ggsave("Swab-characterization_2.tiff",device = "tiff", dpi=300, path = "C:/Users/jfq/Box Sync/NIFA Project/NIFA Project Repository/Output 2022", height= 4, width= 8, units = "in")
ggsave("Swab-characterization_2.jpeg", device = "jpeg", dpi=300, path = "C:/Users/jfq/Box Sync/NIFA Project/NIFA Project Repository/Output 2022", height= 4, width= 8, units = "in")


OvPH_MTvPS<-ggplot(data= Everything_wet, aes(x=Sample.Type, y=Best..Estimate, fill=Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=8),axis.text.x = element_text(size= 8, angle=0, hjust=0.5), title = element_text(size=8))+
  geom_point(alpha=0.5)+
  expand_limits(y=0)+
  scale_fill_manual(values=c("#69b3a2", "white")) +
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7,8,9))+
  scale_x_discrete("Sample Type")+
  labs(fill= "Sample Type")+
  ggtitle("Bacterial Recovery in Aggregative Swab Cloths and Composite Produce Samples")
OvPH_MTvPS+facet_grid(.~ Test)

ggsave("Swab-v-Produce3.tiff",device = "tiff", dpi=300, path = "C:/Users/jfq/Box Sync/NIFA Project/NIFA Project Repository/Output 2022", height= 4, width= 6, units = "in")
ggsave("Swab-v-Produce3.jpeg",device = "jpeg", dpi=300, path = "C:/Users/jfq/Box Sync/NIFA Project/NIFA Project Repository/Output 2022", height= 4, width= 6, units = "in")


#NIFA and IAFP Poster
OvPH_MT_pressure_hyd<-ggplot(data= Only_Swabs_noEcoli, aes(x=Hand.msd, y=Best..Estimate, col= Wet.dry))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=10),axis.text.x = element_text(angle=0, hjust=0.5))+
  geom_point(position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g of swab)", breaks= c(1,2,3,4,5,6,7,8,9))+
  scale_x_discrete("Collection method (Pressure)")+
  scale_color_brewer(palette="Dark2")+
  labs(col="Hydration Factor")+
  ggtitle("Pressure and Hydration Effects on Bacterial Recovery in Aggregative Swab Sampling Technique Optimization")

OvPH_MT_pressure_hyd+facet_grid(.~ Test)

ggsave("Swab-characterization-overall2.tiff",device = "tiff", dpi=300, path = "C:/Users/jfq/Box Sync/NIFA Project/NIFA Project Repository/Output 2022", height= 5, width= 10, units = "in")
ggsave("Swab-characterization-overall2.jpeg", device = "jpeg", dpi=300, path = "C:/Users/jfq/Box Sync/NIFA Project/NIFA Project Repository/Output 2022", height= 5, width= 10, units = "in")

