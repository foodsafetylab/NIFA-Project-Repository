#C:/Users/jorge/Box Sync/NIFA Project/CA 9-2021/Plate results for Laptop
#C:/Users/jfq/Box Sync/NIFA Project/CA 9-2021/Plate results for lab Pc
# C:/Users/rgathm2/Box Sync/NIFA Project/CA 9-2021/Plate results
Pr<-read.csv("R_Colif Swabs and Tissue.csv")
InPst<-read.csv("R_In-PostH_Swabs and Tissues.csv")

#####Run from here
library(tidyverse)
library(ggplot2)

PreHarvest<-read.csv("PreHarvest Trimmed.csv")
Pr0.1<-PreHarvest[,-1]
# colnames(Pr0.1)<-c("Sampling Day", "Sample Type", "Area", "Field Location", "Class", "Strata", "Test", "Best Estimate Log(CFU/g)")
# write.csv(Pr0.1, "PreHarvest Trimmed.csv")
InPostharvest<-read.csv("In-PostHarvest Trimmed.csv")
InPst0.1<-InPostharvest[,-1]
# colnames(InPst0.1)<-c("Sampling Day", "Sample Type", "Location", "Class", "Strata", "Test", "Best Estimate Log(CFU/g)")
# write.csv(InPst0.1, "In-PostHarvest Trimmed.csv")

#Data Filtering and Tidying 
##PreHarvest - (1)Swab, Tissue and Tissue HR in Perimeter 
##PreHarvest - (2) Swab and Tissue in Perimeter v Center
Pr1<- Pr0.1%>% filter(Sample.Type == "Grab"|Sample.Type == "Swab"|Field.Loc=="Bed 1" | Field.Loc== "Edge")
Pr2<- Pr0.1%>% filter(Sample.Type == "Grab"|Sample.Type == "Swab")
#PreHarvest (3) Swab and Tissue Center beds
Pr3<- Pr0.1%>% filter(Area == "Center")
Pr3.1<-Pr3%>% filter(Sample.Type == "Grab"|Sample.Type == "Swab")
#PreHarvest HR v Swab permiter
Pr4<- Pr0.1%>% filter(Sample.Type == "Grab HR"|Sample.Type == "Swab"|Sample.Type=="Grab")
Pr4.1<-Pr4%>% filter(Field.Loc=="Bed 1" | Field.Loc== "Edge"|Field.Loc=="Bed 8")
#PreHarvest HR v Grab perimeter
Pr5<- Pr0.1%>% filter(Sample.Type == "Grab HR"|Sample.Type == "Grab")
Pr5.1<-Pr5%>% filter(Field.Loc=="Bed 1" | Field.Loc== "Edge"|Field.Loc=="Bed 8")

####InHarvest and Post Harvest
IH<- InPst0.1 %>% filter(Sampling.Day == "InHarvest D1"|Sampling.Day == "InHarvest D2")
PH<- InPst0.1 %>% filter(Sampling.Day == "PostHarvest D1"|Sampling.Day == "PostHarvest D2")

#InHarvest Swabs Bins v Chute  v Tissue Bins 
IH_MTvTissue<-IH %>% filter(Sample.Type == "MT"| Sample.Type == "Tissue")
IH_MTTiChuteBin<-IH_MTvTissue %>% filter(Location == "Bins" | Location == "Chute")
#InHarvest Gloves v Tissue 
IH2<-IH %>% filter(Sample.Type == "Gloves"| Sample.Type == "Tissue")
IH2.1<-IH2 %>% filter(Sample.Type == "Gloves" | Location == "Leftover")

#PostHarvest Gloves v MT Chute v Tissue Sides
PH1<-PH %>% filter(Sample.Type == "Gloves"| Location == "Sides" | Location == "Chute")

#PostHarvest Interior Swab v Tissue 
PH2<-PH %>% filter(Location == "Interior")


##Previous Trials
# Pr4<-InPst%>% filter(Sampling.Day == "InHarvest D2")
# Pr3.1<-Pr3%>% filter(Location == "Bins" | Location == "Chute" | Location == "Interior") 
# Pr3.2<-Pr3%>% filter(Location == "Chute" ) 
# Pr4<- InPst %>% filter(Sampling.Day == "InHarvest D2" | Strata == "Bed 1" | Strata == "Bed 5")
# Pr5<- InPst %>% filter(Sampling.Day == "InHarvest D2" | Sampling.Day== "InHarvest D1")
# Pr5<-Pr5[-115,]
# 
# PP<-Pr%>% filter(Sampling.Day == "Preharvest D2")
# PPP<-PPP1%>%filter(Sample.Type == "Grab"|Sample.Type == "Swab")
# PPP1<-PP%>% filter(Field.Loc=="Bed 1" | Field.Loc=="Bed 5")
# PPP<-PPP[,-3]
# 
# II<-InPst%>%filter(Sampling.Day=="InHarvest D2")
# II2<-II%>%filter(Strata=="Bed 1"| Strata == "Bed 5")
# II3<-II2[-26,]
# 
# ant<-ggplot(data= PPP, aes(x=Sample.Type, y=Best.Estimate , col= Field.Loc))+
#   geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=25),axis.text.x = element_text(angle=0, hjust=0.5))+
#   expand_limits(y=0)+
#   scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
#   geom_point(position=position_jitterdodge(),alpha=1)+
#   ggtitle("Matching PreHarvest Bacterial Recovery Comparison - Swab v Composite Tissue ")
# ant+ facet_grid(.~Test)
# 
# ant1<-ggplot(data= II3, aes(x=Sample.Type, y=Best.Estimate.log.CFU.g. , col= Strata))+
#   geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=25),axis.text.x = element_text(angle=0, hjust=0.5))+
#   expand_limits(y=0)+
#   scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
#   geom_point(position=position_jitterdodge(),alpha=1)+
#   ylim(0, 6)+
#   ggtitle("Matching InHarvest Bacterial Recovery Comparison - Swab v Composite Tissue ")
# ant1+ facet_grid(.~Test)

##
#Swab v Tissue Overall Perimeter
OvPH_MTvTissue<-ggplot(data= Pr2, aes(x=Area, y=Best.Estimate , col= Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=25),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("Overall PreHarvest Bacterial Recovery Comparison - Swab v Composite Tissue ")
OvPH_MTvTissue+ facet_grid(.~Test)

#Swab v Tissue v Tissue HR Center v Perimeter
OvPH_MTvCTvHR<-ggplot(data= Pr1, aes(x=Area, y=Best.Estimate, col= Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=25),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("PreHarvest Bacterial Recovery Comparison - Center and Perimeter")
OvPH_MTvCTvHR+ facet_grid(.~Test)
#ggsave(filename = "MT_Tissue&HR Center v Peri.pdf",plot = Papi1, device = pdf, path ="C:/Users/jorge/Box Sync/NIFA Project/CA 9-2021/Plate results", width = 50, height = 50, limitsize = F) 

# #Swab Tissue Center
# Papi2<-ggplot(data= Pr3.1, aes(x=Sample.Type, y=Best.Estimate , col= Area))+
#   geom_boxplot(outlier.colour = NULL)+
#   expand_limits(y=0)+
#   scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
#   geom_point(position=position_jitterdodge(),alpha=1)+
#   ggtitle("PreHarvest Bacterial Recovery Comparison - Samples from the Center")
# Papi2+ facet_grid(.~Test)

#PreHarvest Composite grab and HR v Swab perimeter
PH_MTvCTvHR<-ggplot(data= Pr4.1, aes(x=Sample.Type, y=Best.Estimate , col= Field.Loc))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=25),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("PreHarvest Bacterial Recovery Comparison - Samples in the Perimeter")
PH_MTvCTvHR+ facet_grid(.~Test)

# #PreHarvest HR v Grab perimeter
# Papi4<-ggplot(data= Pr5.1, aes(x=Sample.Type, y=Best.Estimate , col= Field.Loc))+
#   geom_boxplot()+
#   expand_limits(y=0)+
#   scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
#   geom_point(position=position_jitterdodge(),alpha=1)+
#   ggtitle("PreHarvest Bacterial Recovery Comparison - Tissue in the Perimeter")
# Papi4+ facet_grid(.~Test)

#InHarvest Swabs Bins v Chute  v Tissue Bins 
Ltc<-ggplot(data= IH_MTTiChuteBin, aes(x=Location, y=Best.Estimate.log.CFU.g. , col= Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=25),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("In-Harvest Bacterial Recovery Comparison - On top of the trailer")
Ltc + facet_grid(.~Test)

#InHarvest Gloves v Tissue 
Ltc1<-ggplot(data= IH2.1, aes(x=Sample.Type, y=Best.Estimate.log.CFU.g. , col= Location))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=25),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("In Harvest Bacterial Recovery Comparison - Samples taken from the ground")
Ltc1 + facet_grid(.~Test)

#PostHarvest Gloves v MT Chute v Tissue Sides
Ltc2<-ggplot(data= PH1, aes(x=Sample.Type, y=Best.Estimate.log.CFU.g. , col= Location))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=25),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("PostHarvest Bacterial Recovery Comparison - Exterior Sampling")
Ltc2 + facet_grid(.~Test)

#PostHarvest Interior Swab v Tissue 
Ltc3<-ggplot(data= PH2, aes(x=Sample.Type, y=Best.Estimate.log.CFU.g., col= Location))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=25),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("PostHarvest Bacterial Recovery Comparison - Interior Sampling")
Ltc3 + facet_grid(.~Test)

#Overall IH
Ltc5<-ggplot(data= IH, aes(x=Sample.Type, y=Best.Estimate.log.CFU.g. , col= Location))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=25),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("In-Harvest Bacterial Recovery Comparison - Overall Sampling")
Ltc5 + facet_grid(.~Test)

#Overall PH
Ltc4<-ggplot(data= PH, aes(x=Sample.Type, y=Best.Estimate.log.CFU.g., col= Location))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=25),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("PostHarvest Bacterial Recovery Comparison - Overall Sampling")
Ltc4 + facet_grid(.~Test)


##Overall Sample Type Plots Pre/ In/ Post 
Qtzy<-ggplot(Pr)+
  geom_boxplot(aes(x=Sample.Type, y=Best.Estimate , col=Test))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  ggtitle("Overall Bacterial Recovery Comparison")+
  theme_bw(base_size = 10)+
  geom_jitter(aes(x=Sample.Type, y=Best.Estimate))
Qtzy + facet_grid(.~Sampling.Day)

Qtzy2<-ggplot(InPst)+
  geom_boxplot(aes(x=Sample.Type, y=Best.Estimate.log.CFU.g. , col=Test))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_jitter(aes(x=Sample.Type, y=Best.Estimate.log.CFU.g.))
Qtzy2 + facet_grid(.~Sampling.Day)

#Filtered PreHarvest based on sample location - no distinguish between days
Pr9<- Pr %>% filter(Field.Loc=="Bed 1" | Field.Loc == "Bed 8" | Field.Loc=="Edge" | Field.Loc == "Bed 3" | Field.Loc == "Bed 7" | Field.Loc == "Bed 5")

Qtzy3<-ggplot(Pr2)+
  geom_boxplot(aes(x=Sample.Type, y=Best.Estimate , col= Field.Loc))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))

Qtzy3+ facet_grid(.~Test)

#Filtered PreHarvest based on sample location - distinguishing between days

Pr6<- Pr %>% filter(Field.Loc=="Bed 1" | Field.Loc == "Bed 8" | Field.Loc=="Edge" | Field.Loc == "Bed 3" | Field.Loc == "Bed 7" | Field.Loc == "Bed 5")

Qtzy7<-ggplot(Pr2)+
  geom_boxplot(aes(x=Sample.Type, y=Best.Estimate , col= Field.Loc))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))

Qtzy7+ facet_grid(.~Test)


#Filtered In-Harvest and Post Harvest
Pr3<- InPst %>% filter(Sampling.Day == "PostHarvest D1" | Sampling.Day == "PostHarvest D2")
Pr4<-InPst%>% filter(Sampling.Day == "InHarvest D2")
Pr3.1<-Pr3%>% filter(Location == "Bins" | Location == "Chute" | Location == "Interior") 
Pr3.2<-Pr3%>% filter(Location == "Chute" ) 
Pr4<- InPst %>% filter(Sampling.Day == "InHarvest D2" | Strata == "Bed 1" | Strata == "Bed 5")
Pr5<- InPst %>% filter(Sampling.Day == "InHarvest D2" | Sampling.Day== "InHarvest D1")
Pr5<-Pr5[-115,]

#Filter for round 1
Pr5.1<-Pr5 %>% filter(Class == "Round 1")

#Beds v Edge; Pr MT Bed 1&5 v In Gloves Bed 1&5; 

#Pst Gloves, MT CHUTE, MT Interior
# Qtzy4<-ggplot(Pr3)+
#   geom_boxplot(aes(x=Sample.Type, y=Best.Estimate.log.CFU.g. , col= Location))+
#   expand_limits(y=0)+
#   scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
#   geom_jitter(aes(x=Sample.Type, y=Best.Estimate.log.CFU.g.), width = 0.20)
# 
# Qtzy4+ facet_grid(.~Test)

Qtzy4.1<-ggplot(data= Pr3, aes(x=Sample.Type, y=Best.Estimate.log.CFU.g. , col= Location))+
  geom_boxplot()+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)
Qtzy4.1+ facet_grid(.~Test)


Qtzy5<-ggplot(Pr5)+
  geom_boxplot(aes(x=Sample.Type, y=Best.Estimate.log.CFU.g. , col= Class))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))

Qtzy5+ facet_grid(.~Test)

Qtzy6<-ggplot(Pr5.1)+
  geom_boxplot(aes(x=Sample.Type, y=Best.Estimate.log.CFU.g. , col= Class))+ #col= Location))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))

Qtzy6+ facet_grid(.~Test)


## FOR LEE

fLee<- Pr%>% filter(Sample.Type == "Grab" | Sample.Type == "Swab")


Lee<-ggplot(fLee)+
  geom_boxplot(aes(x=Sample.Type, y=Best.Estimate , col=Test))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  ggtitle("Overall Bacterial Recovery Comparison")+
  theme_bw(base_size = 10)
Lee 


## END FOR LEE

str(unique(Pr$Sample.Type))

ssdaf<-Pr[Pr$Sample.Type != "Tissue",]
tissue11<-ssdaf[59,1]

Pr1[59,1] = "Tissue"

###### Beginning RJG
#preharvest swab vs tissue
PreSwabVsTissue<- ggplot(Pr) + geom_boxplot(aes(x= Sample.Type, y= Best.Estimate, col= Test), outlier.shape = NA) +
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7)) +
  geom_point(aes(x= Sample.Type, y= Best.Estimate, col= Test), position=position_jitterdodge(jitter.width = 0.2),alpha=0.5) +
  facet_wrap("Sample.Type", scales="free_x") +
  ggtitle("Preharvest counts of swabs vs tissues")
PreSwabVsTissue

#preharvest bottoms vs sides
BottomvsSide<- ggplot(Pr2) +
  geom_boxplot(aes(x= Sample.Type, y= Best.Estimate, col= Test), outlier.shape = NA) +
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7)) +
  geom_point(aes(x= Sample.Type, y= Best.Estimate, col= Test), position=position_jitterdodge(),alpha=0.5) +
  facet_wrap("Class") +
  xlab("Sample Type") +
  ggtitle("Preharvest counts of bottoms vs sides of romaine")
BottomvsSide

#preharvest day 1, HR strata vs comp
PrD1StrataVsComp<- Pr %>% filter(Sampling.Day == "Preharvest D1")
PrD1StrataVsComp<- PrD1StrataVsComp %>% filter(Field.Loc == "Bed 1" | Field.Loc == "Bed 8" | Field.Loc == "Edge")
PrD1StrataVsCompAPC<- PrD1StrataVsComp %>% filter(Test == "APC")
PrD1StrataVsCompColi<- PrD1StrataVsComp %>% filter(Test == "Coliforms")

#preharvest day 2, HR strata vs comp, UI379 and UI380 strata is listed inner, should be strata 1
PrD2StrataVsComp<- Pr %>% filter(Sampling.Day == "Preharvest D2")
PrD2StrataVsComp<- PrD2StrataVsComp %>% filter(Field.Loc == "Bed 1" | Field.Loc == "Bed 8" | Field.Loc == "Edge")
PrD2StrataVsCompAPC<- PrD2StrataVsComp %>% filter(Test == "APC")
PrD2StrataVsCompColi<- PrD2StrataVsComp %>% filter(Test == "Coliforms")

StrataVsCompAPC<- ggplot(PrD1StrataVsCompAPC, aes(x= Strata, y= Best.Estimate), col= Sample.Type) +
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7)) +
  facet_wrap("Field.Loc", scales = "free_x") +
  labs(x= "Bed and Strata Location") +
  geom_point(aes(x= Strata, y= Best.Estimate, col= Sample.Type)) +
  theme(axis.text.x = element_text(angle=90)) +
  ggtitle("Preharvest Day 1 HR Strata vs. Composite APCs")
StrataVsCompAPC

StrataVsCompColi<- ggplot(PrD1StrataVsCompColi, aes(x= Strata, y= Best.Estimate), col= Sample.Type) +
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7)) +
  facet_wrap("Field.Loc", scales = "free_x") +
  labs(x= "Bed and Strata Location") +
  geom_point(aes(x= Strata, y= Best.Estimate, col= Sample.Type)) +
  theme(axis.text.x = element_text(angle=90)) +
  ggtitle("Preharvest Day 1 HR Strata vs. Composite Coliforms")
StrataVsCompColi
