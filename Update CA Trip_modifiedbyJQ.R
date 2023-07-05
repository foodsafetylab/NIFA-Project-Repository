library(tidyverse)
library(ggplot2)
library(rstatix)

#load workspace, global environment PreIAFP_NIFA_analysis.RData
master<- read.csv("C:/Users/jfq/Box Sync/NIFA Project/CA 9-2021/Plate results/Master plate counts.csv")
count_master<-master%>%group_by(Sample.Time.Point, Sample.Description, Test)%>%summarise(count=n())

##File Reading
#PreHarvest<-read.csv("PreHarvest Trimmed.csv")
#Pr0.1<-PreHarvest[,-1]
#InPostharvest<-read.csv("In-PostHarvest Trimmed.csv")
#InPst0.1<-InPostharvest[,-1]

#Data Filtering and Tidying 

##PreHarvest - (1)Swab, Tissue and Tissue HR in Perimeter
##PreHarvest - (2) Swab and Tissue in Perimeter v Center - Overall no HR

Overall_PH<- PreHarvest
Perimeter1<- PreHarvest%>% filter(Field.Loc=="Bed 1" | Field.Loc== "Edge" | Field.Loc == "Bed 8")

PreHarvest<- master %>% filter(Sample.Time.Point == "Preharvest")

#Visualization of Overall Sample Type Performance

OvSample<-ggplot(data= PreHarvest, aes(x=Sample.Type, y=Log.CFU.g., col= Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=15),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  stat_summary(aes(x=Sample.Type, y=Log.CFU.g.), fun="mean", shape=4) +
  facet_wrap(~Test, scales = "free_x") + 
  ggtitle("Overall PreHarvest Bacterial Recovery Comparison")
OvSample

#APC
PreHarvest_swab_APC<- subset(PreHarvest, Sample.Type== "Cloth" & Test== "APC")
range(PreHarvest_swab_APC$Log.CFU.g.)
mean(PreHarvest_swab_APC$Log.CFU.g.)
median(PreHarvest_swab_APC$Log.CFU.g.)
sd(PreHarvest_swab_APC$Log.CFU.g.)

PreHarvest_tissue_APC<- subset(PreHarvest, Sample.Type== "Produce" & Test== "APC")
range(PreHarvest_tissue_APC$Log.CFU.g.)
mean(PreHarvest_tissue_APC$Log.CFU.g.)
median(PreHarvest_tissue_APC$Log.CFU.g.)
sd(PreHarvest_tissue_APC$Log.CFU.g.)

PreHarvest_tissueHR_APC<- subset(PreHarvest, Sample.Type== "Produce High Resolution" & Test== "APC")
range(PreHarvest_tissueHR_APC$Log.CFU.g.)
mean(PreHarvest_tissueHR_APC$Log.CFU.g.)
median(PreHarvest_tissueHR_APC$Log.CFU.g.)
sd(PreHarvest_tissueHR_APC$Log.CFU.g.)

#coliforms
PreHarvest_swab_coli<- subset(PreHarvest, Sample.Type== "Cloth" & Test== "Coliform")
range(PreHarvest_swab_coli$Log.CFU.g.)
mean(PreHarvest_swab_coli$Log.CFU.g.)
median(PreHarvest_swab_coli$Log.CFU.g.)
sd(PreHarvest_swab_coli$Log.CFU.g.)

PreHarvest_tissue_coli<- subset(PreHarvest, Sample.Type== "Produce" & Test== "Coliform")
range(PreHarvest_tissue_coli$Log.CFU.g.)
mean(PreHarvest_tissue_coli$Log.CFU.g.)
median(PreHarvest_tissue_coli$Log.CFU.g.)
sd(PreHarvest_tissue_coli$Log.CFU.g.)

PreHarvest_tissueHR_coli<- subset(PreHarvest, Sample.Type== "Produce High Resolution" & Test== "Coliform")
range(PreHarvest_tissueHR_coli$Log.CFU.g.)
mean(PreHarvest_tissueHR_coli$Log.CFU.g.)
median(PreHarvest_tissueHR_coli$Log.CFU.g.)
sd(PreHarvest_tissueHR_coli$Log.CFU.g.)

PreHarvestAPC<- subset(PreHarvest, Test== "APC")
PreHarvestColiform<- subset(PreHarvest, Test== "Coliform")

#anova on sample type and top vs sides and sample day
#normality assumption, p= 0.0005
PreHarvest_anova_APC<- lm(Log.CFU.g. ~ Sample.Type*Location*Sample.Day, data=PreHarvestAPC)
PreHarvest_lm_APC<- lm(Log.CFU.g. ~ Sample.Type+Field.Location+Location+Sample.Type*Location+Sample.Type*Field.Location+Sample.Type*Location*Field.Location, data=PreHarvestAPC)

shapiro_test(residuals(PreHarvest_anova_APC))
#homogeneity of variance, p=0.04
PreHarvestAPC %>% levene_test(Log.CFU.g. ~ Sample.Type*Location*Sample.Day)
#anova
anova_test(PreHarvest_anova_APC)
PreHarvest_APC_ANOVA<-anova(PreHarvest_lm_APC)
#tukey
tukey_hsd(PreHarvestAPC, Log.CFU.g. ~ Sample.Type)
tukey_hsd(PreHarvestAPC, Log.CFU.g. ~ Sample.Day)
tukey_hsd(PreHarvestAPC, Log.CFU.g. ~ Sample.Type*Location)

#normality assumption, p= 
PreHarvest_anova_coli<- lm(Log.CFU.g. ~ Sample.Type*Location*Sample.Day, data=PreHarvestColiform)
shapiro_test(residuals(PreHarvest_anova_coli))
#homogeneity of variance, p=
PreHarvestColiform %>% levene_test(Log.CFU.g. ~ Sample.Type*Location*Sample.Day)
#anova
anova_test(PreHarvest_anova_coli)
#tukey
tukey_hsd(PreHarvestColiform, Log.CFU.g. ~ Sample.Type)
tukey_hsd(PreHarvestColiform, Log.CFU.g. ~ Location)
tukey_hsd(PreHarvestColiform, Log.CFU.g. ~ Sample.Type*Location)

#normality assumption, p= 
PreHarvest_anova<- lm(Log.CFU.g. ~ Sample.Type*Location*Sample.Day, data=PreHarvest)
shapiro_test(residuals(PreHarvest_anova))
#homogeneity of variance, p=
PreHarvest %>% levene_test(Log.CFU.g. ~ Sample.Type*Location)
#anova
anova_test(PreHarvest_anova)

#variance test
PreHarvestAPC_tis_swab<- subset(PreHarvestAPC, Sample.Type== "Produce" | Sample.Type== "Cloth")
var.test(Log.CFU.g. ~ Sample.Type, data = PreHarvestAPC_tis_swab)#Sig dif, cloth higher than produce

PreHarvestAPC_tisHR_swab<- subset(PreHarvestAPC, Sample.Type== "Produce High Resolution" | Sample.Type== "Cloth")
var.test(Log.CFU.g. ~ Sample.Type, data = PreHarvestAPC_tisHR_swab)#Sig dif, cloth higher than produce

PreHarvestAPC_tisHR_tis<- subset(PreHarvestAPC, Sample.Type== "Produce High Resolution" | Sample.Type== "Produce")
var.test(Log.CFU.g. ~ Sample.Type, data = PreHarvestAPC_tisHR_tis)#Not Sig dif

PreHarvestcoli_tis_swab<- subset(PreHarvestColiform, Sample.Type== "Produce" | Sample.Type== "Cloth")
var.test(Log.CFU.g. ~ Sample.Type, data = PreHarvestcoli_tis_swab)#Sig dif, cloth lower than produce

PreHarvestcoli_tisHR_swab<- subset(PreHarvestColiform, Sample.Type== "Produce High Resolution" | Sample.Type== "Cloth")
var.test(Log.CFU.g. ~ Sample.Type, data = PreHarvestcoli_tisHR_swab)#Sig dif, cloth lower than produce

PreHarvestcoli_tisHR_tis<- subset(PreHarvestColiform, Sample.Type== "Produce High Resolution" | Sample.Type== "Produce")
var.test(Log.CFU.g. ~ Sample.Type, data = PreHarvestcoli_tisHR_tis)#Not Sig dif

# Visualization of Overall Pre Harvest Sampling - No distiction bw days

OvPH<-ggplot(data= PreHarvest, aes(x=Field.Location, y=Log.CFU.g., col= Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=15),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  stat_summary(aes(x=Field.Location, y=Log.CFU.g.), fun="mean", shape=4) +
  ggtitle("Overall PreHarvest Bacterial Recovery Comparison")
OvPH+ facet_grid(.~Test, scales = "free_x")

#anova on sample type and top vs sides and sample day
#normality assumption, p= 
PreHarvest_anova_APC_edge<- lm(Log.CFU.g. ~ Sample.Type*Field.Location, data=PreHarvestAPC)
shapiro_test(residuals(PreHarvest_anova_APC_edge))
#homogeneity of variance, p=
PreHarvestAPC %>% levene_test(Log.CFU.g. ~ Sample.Type*Field.Location)
#anova
anova_test(PreHarvest_anova_APC_edge)
#tukey
tukey_hsd(PreHarvestAPC, Log.CFU.g. ~ Sample.Type)

#normality assumption, p= 
PreHarvest_anova_coli_edge<- lm(Log.CFU.g. ~ Sample.Type*Field.Location, data=PreHarvestColiform)
shapiro_test(residuals(PreHarvest_anova_coli_edge))
#homogeneity of variance, p=
PreHarvestColiform %>% levene_test(Log.CFU.g. ~ Sample.Type*Field.Location)
#anova
anova_test(PreHarvest_anova_coli_edge)
#tukey
tukey_hsd(PreHarvestColiform, Log.CFU.g. ~ Sample.Type)
interaction_type_location<- tukey_hsd(PreHarvestColiform, Log.CFU.g. ~ Sample.Type*Field.Location)



#PreHarvest HR v Swab Perimeter and HR v Composite Tissue

#PreHarvest HR v Swab permiter
Pr4<-PreHarvest%>% filter(Sample.Type == "Tissue High Resolution"|Sample.Type == "Swab")
Pr4.1<-Pr4%>% filter(Bed.Number=="Bed 1" | Bed.Number== "Edge"|Bed.Number=="Bed 8")
#PreHarvest HR v Grab perimeter
Pr5<-PreHarvest%>% filter(Sample.Type == "Tissue High Resolution"|Sample.Type == "Tissue")
Pr5.1<-Pr5%>% filter(Bed.Number=="Bed 1" | Bed.Number== "Edge"| Bed.Number=="Bed 8")

##Visualization of MT v HR

MTvHR<-ggplot(data= Pr4.1, aes(x=Sample.Type, y= Log.CFU.g. , col= Bed.Number))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=15),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("PreHarvest Bacterial Recovery Comparison - Samples in the Perimeter")
MTvHR+ facet_grid(.~Test)

#Visualization of Composite Tissue (CT) and HR

CTvHR<-ggplot(data= Pr5.1, aes(x=Sample.Type, y= Log.CFU.g. , col= Bed.Number))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=15),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("PreHarvest Bacterial Recovery Comparison - Tissue Samples in the Perimeter")
CTvHR+ facet_grid(.~Test)

#Filtering for Grab and Swab for Sides v Bottom comparison

SidesBottom<-PreHarvest%>% filter(Sample.Type == "Tissue"|Sample.Type == "Swab")

#Visualization of Side v Bottom in MT v CT

SideBottom_CTvMT<-ggplot(data= SidesBottom, aes(x=Sample.Type, y=Log.CFU.g., col= Location))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=15),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("PreHarvest Bacterial Recovery Comparison - Specimen Location in Lettuce Head")
SideBottom_CTvMT+facet_grid(.~Test)

####Overall Look at In-Harvest Sample Type Performance 

#IH<- InPst0.1 %>% filter(Sampling.Day == "InHarvest D1"|Sampling.Day == "InHarvest D2")
#IH_NoGloves<-IH %>% filter(Sample.Type == "Swab"| Sample.Type == "Grab")
#IH_OnTrailer<-IH_NoGloves %>% filter(Location == "Bins" | Location == "Chute")
IH<- subset(master, Sample.Time.Point== "Harvest")
In_OnTrailer<- subset(IH, Location=="Chute" | Location == "Bins")

##Visualization Of Comparison between sampling performance of MT and Composite Tissue taken from the Chute and Bins

Chute_Bins_OnTrailer<-ggplot(data= In_OnTrailer, aes(x=Location, y=Log.CFU.g. , col= Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=15),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  stat_summary(aes(x=Location, y=Log.CFU.g.), fun="mean", shape=4) +
  ggtitle("In-Harvest Bacterial Recovery Comparison - On top of the trailer")
Chute_Bins_OnTrailer + facet_grid(.~Test)

HarvestAPC<- subset(IH, Test== "APC")
HarvestColiform<- subset(IH, Test== "Coliform")
In_OnTrailerAPC<- subset(In_OnTrailer, Test== "APC")
In_OnTrailerColi<- subset(In_OnTrailer, Test== "Coliform")

###Not needed if analysis is only done overall HARVEST
#APC
Harvest_chute_APC<- subset(HarvestAPC, Location== "Chute")
range(Harvest_chute_APC$Log.CFU.g.)
mean(Harvest_chute_APC$Log.CFU.g.)
median(Harvest_chute_APC$Log.CFU.g.)
sd(Harvest_chute_APC$Log.CFU.g.)

Harvest_bin_swab_APC<- subset(HarvestAPC, Sample.Type== "Cloth" & Location== "Bins")
range(Harvest_bin_swab_APC$Log.CFU.g.)
mean(Harvest_bin_swab_APC$Log.CFU.g.)
median(Harvest_bin_swab_APC$Log.CFU.g.)
sd(Harvest_bin_swab_APC$Log.CFU.g.)

Harvest_bin_tissue_APC<- subset(HarvestAPC, Sample.Type== "Produce" & Location== "Bins")
range(Harvest_bin_tissue_APC$Log.CFU.g.)
mean(Harvest_bin_tissue_APC$Log.CFU.g.)
median(Harvest_bin_tissue_APC$Log.CFU.g.)
sd(Harvest_bin_tissue_APC$Log.CFU.g.)

#coliform
Harvest_chute_coli<- subset(HarvestColiform, Location== "Chute")
range(Harvest_chute_coli$Log.CFU.g.)
mean(Harvest_chute_coli$Log.CFU.g.)
median(Harvest_chute_coli$Log.CFU.g.)
sd(Harvest_chute_coli$Log.CFU.g.)

Harvest_bin_swab_coli<- subset(HarvestColiform, Sample.Type== "Cloth" & Location== "Bins")
range(Harvest_bin_swab_coli$Log.CFU.g.)
mean(Harvest_bin_swab_coli$Log.CFU.g.)
median(Harvest_bin_swab_coli$Log.CFU.g.)
sd(Harvest_bin_swab_coli$Log.CFU.g.)

Harvest_bin_tissue_coli<- subset(HarvestColiform, Sample.Type== "Produce" & Location== "Bins")
range(Harvest_bin_tissue_coli$Log.CFU.g.)
mean(Harvest_bin_tissue_coli$Log.CFU.g.)
median(Harvest_bin_tissue_coli$Log.CFU.g.)
sd(Harvest_bin_tissue_coli$Log.CFU.g.)
####


#anova on sample type and location and sample day for overall harvest 
#normality assumption, p= 
Harvest_lm_APC<- lm(Log.CFU.g. ~ Sample.Description, data=HarvestAPC)
shapiro_test(residuals(Harvest_trailer_lm_APC))
#homogeneity of variance, p=
#In_OnTrailerAPC %>% levene_test(Log.CFU.g. ~ Sample.Type*Location*Sample.Day)
#anova
Harvest_anova_APC<-anova(Harvest_lm_APC)
#tukey
Dif_Harvest_anova_APC<-tukey_hsd(HarvestAPC, Log.CFU.g. ~ Sample.Description)
#tukey_hsd(In_OnTrailerAPC, Log.CFU.g. ~ Sample.Day)

#variance test

HarvestAPC_trailer_bin_swab_chute<- subset(In_OnTrailerAPC, Location== "Chute" | Sample.Type== "Swab" & Location== "Bins")
var.test(Log.CFU.g. ~ Location, data = HarvestAPC_trailer_bin_swab_chute)

HarvestAPC_trailer_bin_tis_chute<- subset(In_OnTrailerAPC, Location== "Chute" | Sample.Type== "Tissue" & Location== "Bins")
var.test(Log.CFU.g. ~ Location, data = HarvestAPC_trailer_bin_tis_chute)

HarvestAPC_trailer_bin_tis_swab<- subset(In_OnTrailerAPC, Location== "Bins")
var.test(Log.CFU.g. ~ Sample.Type, data = HarvestAPC_trailer_bin_tis_swab)

#anova on sample type and location and sample day for harvest on trailer
#normality assumption, p= 
Harvest_trailer_anova_coli<- lm(Log.CFU.g. ~ Sample.Type*Location*Sample.Day, data=In_OnTrailerColi)
shapiro_test(residuals(Harvest_trailer_anova_coli))
#homogeneity of variance, p=
In_OnTrailerColi %>% levene_test(Log.CFU.g. ~ Sample.Type*Location*Sample.Day)
#anova
anova_test(Harvest_trailer_anova_coli)
#tukey
tukey_hsd(In_OnTrailerColi, Log.CFU.g. ~ Sample.Type)
tukey_hsd(In_OnTrailerColi, Log.CFU.g. ~ Sample.Day)
tukey_hsd(In_OnTrailerColi, Log.CFU.g. ~ Sample.Type*Sample.Day)
#variance test
Harvestcoli_trailer_bin_swab_chute<- subset(In_OnTrailerColi, Location== "Chute" | Sample.Type== "Swab" & Location== "Bins")
var.test(Log.CFU.g. ~ Location, data = Harvestcoli_trailer_bin_swab_chute)

Harvestcoli_trailer_bin_tis_chute<- subset(In_OnTrailerColi, Location== "Chute" | Sample.Type== "Tissue" & Location== "Bins")
var.test(Log.CFU.g. ~ Location, data = Harvestcoli_trailer_bin_tis_chute)

Harvestcoli_trailer_bin_tis_swab<- subset(In_OnTrailerColi, Location== "Bins")
var.test(Log.CFU.g. ~ Sample.Type, data = Harvestcoli_trailer_bin_tis_swab)


#IH_NoMT<-IH %>% filter(Sample.Type == "Gloves"| Sample.Type == "Grab")
IH_Ground<- subset(IH, Location == "Leftover" | Location == "Glove")
IH_Ground_APC<- subset(IH_Ground, Test== "APC")
IH_Ground_coli<- subset(IH_Ground, Test== "Coliform")

#Visualization Of Comparison between sampling performance of Gloves and Tissue Leftover Trims taken on the ground
GlovesLeftovers_Ground<-ggplot(data= IH_Ground, aes(x=Location, y=Log.CFU.g. , col= Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=15),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("In-Harvest Bacterial Recovery Comparison - On the ground") +
  stat_summary(aes(x=Location, y=Log.CFU.g.), fun="mean", shape=4)
GlovesLeftovers_Ground + facet_grid(.~Test)

#APC
Harvest_glove_APC<- subset(HarvestAPC, Location== "Glove")
range(Harvest_glove_APC$Log.CFU.g.)
mean(Harvest_glove_APC$Log.CFU.g.)
median(Harvest_glove_APC$Log.CFU.g.)
sd(Harvest_glove_APC$Log.CFU.g.)

Harvest_leftover_APC<- subset(HarvestAPC, Location== "Leftover")
range(Harvest_leftover_APC$Log.CFU.g.)
mean(Harvest_leftover_APC$Log.CFU.g.)
median(Harvest_leftover_APC$Log.CFU.g.)
sd(Harvest_leftover_APC$Log.CFU.g.)

#coliforms
Harvest_glove_coli<- subset(HarvestColiform, Location== "Glove")
range(Harvest_glove_coli$Log.CFU.g.)
mean(Harvest_glove_coli$Log.CFU.g.)
median(Harvest_glove_coli$Log.CFU.g.)
sd(Harvest_glove_coli$Log.CFU.g.)

Harvest_leftover_coli<- subset(HarvestColiform, Location== "Leftover")
range(Harvest_leftover_coli$Log.CFU.g.)
mean(Harvest_leftover_coli$Log.CFU.g.)
median(Harvest_leftover_coli$Log.CFU.g.)
sd(Harvest_leftover_coli$Log.CFU.g.)

#anova on sample type and sample day for harvest on ground
#normality assumption, p= 
Harvest_ground_anova_APC<- lm(Log.CFU.g. ~ Sample.Type*Sample.Day, data=IH_Ground_APC)
shapiro_test(residuals(Harvest_ground_anova_APC))
#homogeneity of variance, p=
IH_Ground_APC %>% levene_test(Log.CFU.g. ~ Sample.Type*Sample.Day)
#anova
anova_test(Harvest_ground_anova_APC)
#tukey
tukey_hsd(IH_Ground_APC, Log.CFU.g. ~ Sample.Type)
#variance test
var.test(Log.CFU.g. ~ Sample.Type, data = IH_Ground_APC)

#anova on sample type and sample day for harvest on ground
#normality assumption, p= 
Harvest_ground_anova_coli<- lm(Log.CFU.g. ~ Sample.Type*Sample.Day, data=IH_Ground_coli)
shapiro_test(residuals(Harvest_ground_anova_coli))
#homogeneity of variance, p=
IH_Ground_coli %>% levene_test(Log.CFU.g. ~ Sample.Type*Sample.Day)
#anova
anova_test(Harvest_ground_anova_coli)
#tukey
tukey_hsd(IH_Ground_coli, Log.CFU.g. ~ Sample.Type*Sample.Day)
#variance test
var.test(Log.CFU.g. ~ Sample.Type, data = IH_Ground_coli)


####Overall Look at Post-Harvest Sample Type Performance 

#PostH<- InPst0.1 %>% filter(Sampling.Day == "PostHarvest D1"|Sampling.Day == "PostHarvest D2")
PostH<- subset(all_data, Sample.Time.Point == "Post Harvest")
Post_Exterior<-PostH %>% filter(Location == "Glove"| Location == "Produce sides" | Location == "Chute")
Post_Interior<-PostH %>% filter(Location == "Produce interior")

#External Sampling only 

ExtSampling<-ggplot(data= Post_Exterior, aes(x=Location, y=Log.CFU.g. , col= Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=15),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  stat_summary(aes(x=Location, y=Log.CFU.g.), fun="mean", shape=4) +
  ggtitle("In-Harvest Bacterial Recovery Comparison - Exterior Sampling")
ExtSampling + facet_grid(.~Test)

#Internal Sampling only

IntSampling<-ggplot(data= Post_Interior, aes(x=Location, y=Log.CFU.g., col= Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=15),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  stat_summary(aes(x=Location, y=Log.CFU.g.), fun="mean", shape=4) +
  ggtitle("In-Harvest Bacterial Recovery Comparison - Interior Sampling")
IntSampling + facet_grid(.~Test)


#Internal v External Sampling performance

IntVsExt<-ggplot(data= PostH, aes(x=Sample.Type, y=Log.CFU.g., col=Location))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=15),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  stat_summary(aes(x=Sample.Type, y=Log.CFU.g.), fun="mean", shape=4) +
  ggtitle("Post-Harvest Bacterial Recovery Comparison - Interior vs Exterior Sampling")+
  facet_grid(Test~Strata, scales="free")
IntVsExt 

PostHAPC<- subset(PostH, Test == "APC")
PostHcoli<- subset(PostH, Test == "Coliform")

#APC
Post_glove_APC<- subset(PostHAPC, Location== "Glove")
range(Post_glove_APC$Log.CFU.g.)
mean(Post_glove_APC$Log.CFU.g.)
median(Post_glove_APC$Log.CFU.g.)
sd(Post_glove_APC$Log.CFU.g.)

Post_chute_APC<- subset(PostHAPC, Location== "Chute")
range(Post_chute_APC$Log.CFU.g.)
mean(Post_chute_APC$Log.CFU.g.)
median(Post_chute_APC$Log.CFU.g.)
sd(Post_chute_APC$Log.CFU.g.)

Post_tisex_APC<- subset(PostHAPC, Location== "Tissue sides")
range(Post_tisex_APC$Log.CFU.g.)
mean(Post_tisex_APC$Log.CFU.g.)
median(Post_tisex_APC$Log.CFU.g.)
sd(Post_tisex_APC$Log.CFU.g.)

Post_tisin_APC<- subset(PostHAPC, Location== "Tissue interior" & Sample.Type== "Tissue")
range(Post_tisin_APC$Log.CFU.g.)
mean(Post_tisin_APC$Log.CFU.g.)
median(Post_tisin_APC$Log.CFU.g.)
sd(Post_tisin_APC$Log.CFU.g.)

Post_swabin_APC<- subset(PostHAPC, Location== "Tissue interior" & Sample.Type== "Swab")
range(Post_swabin_APC$Log.CFU.g.)
mean(Post_swabin_APC$Log.CFU.g.)
median(Post_swabin_APC$Log.CFU.g.)
sd(Post_swabin_APC$Log.CFU.g.)

#coliform
Post_glove_coli<- subset(PostHcoli, Location== "Glove")
range(Post_glove_coli$Log.CFU.g.)
mean(Post_glove_coli$Log.CFU.g.)
median(Post_glove_coli$Log.CFU.g.)
sd(Post_glove_coli$Log.CFU.g.)

Post_chute_coli<- subset(PostHcoli, Location== "Chute")
range(Post_chute_coli$Log.CFU.g.)
mean(Post_chute_coli$Log.CFU.g.)
median(Post_chute_coli$Log.CFU.g.)
sd(Post_chute_coli$Log.CFU.g.)

Post_tisex_coli<- subset(PostHcoli, Location== "Tissue sides")
range(Post_tisex_coli$Log.CFU.g.)
mean(Post_tisex_coli$Log.CFU.g.)
median(Post_tisex_coli$Log.CFU.g.)
sd(Post_tisex_coli$Log.CFU.g.)

Post_tisin_coli<- subset(PostHcoli, Location== "Tissue interior" & Sample.Type== "Tissue")
range(Post_tisin_coli$Log.CFU.g.)
mean(Post_tisin_coli$Log.CFU.g.)
median(Post_tisin_coli$Log.CFU.g.)
sd(Post_tisin_coli$Log.CFU.g.)

Post_swabin_coli<- subset(PostHcoli, Location== "Tissue interior" & Sample.Type== "Swab")
range(Post_swabin_coli$Log.CFU.g.)
mean(Post_swabin_coli$Log.CFU.g.)
median(Post_swabin_coli$Log.CFU.g.)
sd(Post_swabin_coli$Log.CFU.g.)

#anova on sample type, sample day, head location for post harvest
#normality assumption, p= 
Post_anova_APC<- lm(Log.CFU.g. ~ Sample.Type*Sample.Day*Strata, data=PostHAPC)
shapiro_test(residuals(Post_anova_APC))
#homogeneity of variance, p=
PostHAPC %>% levene_test(Log.CFU.g. ~ Sample.Type*Sample.Day*Strata)
#anova
anova_test(Post_anova_APC)
#tukey
tukey_hsd(PostHAPC, Log.CFU.g. ~ Strata)

#normality assumption, p= 
Post_anova_coli<- lm(Log.CFU.g. ~ Sample.Type*Sample.Day*Strata, data=PostHcoli)
shapiro_test(residuals(Post_anova_coli))
#homogeneity of variance, p=
PostHcoli %>% levene_test(Log.CFU.g. ~ Sample.Type*Sample.Day*Strata)
#anova
anova_test(Post_anova_coli)
#variance test
GloveSwabEx<- subset(PostHAPC, Sample.Type== "Glove" | Sample.Type== "Swab" & Location== "Chute")
var.test(Log.CFU.g. ~ Sample.Type, data = GloveSwabEx)
GloveTissueEx<- subset(PostHAPC, Sample.Type== "Glove" | Sample.Type== "Tissue" & Location== "Tissue sides")
var.test(Log.CFU.g. ~ Sample.Type, data = GloveTissueEx)
GloveSwabIn<- subset(PostHAPC, Sample.Type== "Glove" | Sample.Type== "Swab" & Location== "Tissue interior")
var.test(Log.CFU.g. ~ Sample.Type, data = GloveSwabIn)
GloveTissueIn<- subset(PostHAPC, Sample.Type== "Glove" | Sample.Type== "Tissue" & Location== "Tissue interior")
var.test(Log.CFU.g. ~ Sample.Type, data = GloveTissueIn)
ChuteTissueEx<- subset(PostHAPC, Sample.Type== "Swab" & Location== "Chute" | Sample.Type== "Tissue" & Location== "Tissue sides")
var.test(Log.CFU.g. ~ Sample.Type, data = ChuteTissueEx)
ChuteSwabIn<- subset(PostHAPC, Sample.Type== "Swab")
var.test(Log.CFU.g. ~ Location, data = ChuteSwabIn)
ChuteTissueIn<- subset(PostHAPC, Sample.Type== "Swab" & Location== "Chute" | Sample.Type== "Tissue" & Location== "Tissue interior")
var.test(Log.CFU.g. ~ Sample.Type, data = ChuteTissueIn)
TissueExSwabIn<- subset(PostHAPC, Sample.Type== "Swab" & Location== "Tissue interior" | Sample.Type== "Tissue" & Location== "Tissue sides")
var.test(Log.CFU.g. ~ Sample.Type, data = TissueExSwabIn)
TissueExTissueIn<- subset(PostHAPC, Sample.Type== "Tissue" & Location== "Tissue interior" | Sample.Type== "Tissue" & Location== "Tissue sides")
var.test(Log.CFU.g. ~ Location, data = TissueExTissueIn)
SwabInTissueIn<- subset(PostHAPC, Sample.Type== "Swab" & Location== "Tissue interior" | Sample.Type== "Tissue" & Location== "Tissue interior")
var.test(Log.CFU.g. ~ Sample.Type, data = SwabInTissueIn)

GloveSwabExcoli<- subset(PostHcoli, Sample.Type== "Glove" | Sample.Type== "Swab" & Location== "Chute")
var.test(Log.CFU.g. ~ Sample.Type, data = GloveSwabExcoli)
GloveTissueExcoli<- subset(PostHcoli, Sample.Type== "Glove" | Sample.Type== "Tissue" & Location== "Tissue sides")
var.test(Log.CFU.g. ~ Sample.Type, data = GloveTissueExcoli)
GloveSwabIncoli<- subset(PostHcoli, Sample.Type== "Glove" | Sample.Type== "Swab" & Location== "Tissue interior")
var.test(Log.CFU.g. ~ Sample.Type, data = GloveSwabIncoli)
GloveTissueIncoli<- subset(PostHcoli, Sample.Type== "Glove" | Sample.Type== "Tissue" & Location== "Tissue interior")
var.test(Log.CFU.g. ~ Sample.Type, data = GloveTissueIncoli)
ChuteTissueExcoli<- subset(PostHcoli, Sample.Type== "Swab" & Location== "Chute" | Sample.Type== "Tissue" & Location== "Tissue sides")
var.test(Log.CFU.g. ~ Sample.Type, data = ChuteTissueExcoli)
ChuteSwabIncoli<- subset(PostHcoli, Sample.Type== "Swab")
var.test(Log.CFU.g. ~ Location, data = ChuteSwabIncoli)
ChuteTissueIncoli<- subset(PostHcoli, Sample.Type== "Swab" & Location== "Chute" | Sample.Type== "Tissue" & Location== "Tissue interior")
var.test(Log.CFU.g. ~ Sample.Type, data = ChuteTissueIncoli)
TissueExSwabIncoli<- subset(PostHcoli, Sample.Type== "Swab" & Location== "Tissue interior" | Sample.Type== "Tissue" & Location== "Tissue sides")
var.test(Log.CFU.g. ~ Sample.Type, data = TissueExSwabIncoli)
TissueExTissueIncoli<- subset(PostHcoli, Sample.Type== "Tissue" & Location== "Tissue interior" | Sample.Type== "Tissue" & Location== "Tissue sides")
var.test(Log.CFU.g. ~ Location, data = TissueExTissueIncoli)
SwabInTissueIncoli<- subset(PostHcoli, Sample.Type== "Swab" & Location== "Tissue interior" | Sample.Type== "Tissue" & Location== "Tissue interior")
var.test(Log.CFU.g. ~ Sample.Type, data = SwabInTissueIncoli)

##All data plot
all_data$Sample.Time.Point_factor<- factor(all_data$Sample.Time.Point, levels=c('Preharvest', 'Harvest', 'Post Harvest'))
all_data$Sample.Description_factor<- factor(all_data$Sample.Description, levels=c('Cloth produce bottom', 'Cloth produce sides', 'Produce bottom', 'Produce sides', 'Produce sides HR', 'Produce leftover trims', 'Glove', 'Chute', 'Produce bins', 'Cloth bin tops', 'Cloth produce interior', 'Produce interior'))
all_plot<-ggplot(data= all_data, aes(x=Sample.Description_factor, y=Log.CFU.g., col=Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=15),axis.text.x = element_text(angle=90, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  stat_summary(aes(x=Sample.Description_factor, y=Log.CFU.g.), fun="mean", shape=4) +
  ggtitle("Bacterial Recovery Overview")+
  xlab("Sample Description") +
  facet_grid(Test~Sample.Time.Point_factor, scales="free")
all_plot

###Appendix

# Visualization of Perimeter Effect - Appendix

OvPerimeter<-ggplot(data= Perimeter1, aes(x=Field.Loc, y=Best.Estimate , col= Sample.Type))+
  geom_boxplot(outlier.alpha = 0)+theme(text = element_text(size=15),axis.text.x = element_text(angle=0, hjust=0.5))+
  expand_limits(y=0)+
  scale_y_continuous("Best Estimate (Log CFU/g)", breaks= c(1,2,3,4,5,6,7))+
  geom_point(position=position_jitterdodge(),alpha=1)+
  ggtitle("Overall PreHarvest Bacterial Recovery Comparison - Swab v Composite Tissue ")+
  facet_grid(Field.Loc~Test, scales="free")
OvPerimeter

#PreHarvest Swab and Tissue only in Center - Appendix

#PreHarvest (3) Swab and Tissue Center beds
Pr3<- Overall_PH%>% filter(Area == "Center")
Pr3.1<-Pr3%>% filter(Sample.Type == "Grab"|Sample.Type == "Swab")
