##C:/Users/jorge/Box Sync/NIFA Project/CA 9-2021/Plate results/
##"C:/Users/rgathm2/Box Sync/NIFA Project/CA 9-2021/Plate results"
##No need if you open the project on C:\Users\jfq\Box Sync\NIFA Project\NIFA Project Repository
library(tidyverse)
library(MASS)
library(car)
library(FSA)
#Subsetting data for analysis
PreHarvest<-read.csv("PreHarvest Trimmed.csv")

PreHarvestAPC<-PreHarvest%>%filter(Test=="APC")
PreHarvestColif<-PreHarvest%>%filter(Test=="Coliforms")
ASHR<-PreHarvestAPC%>%filter(Sample.Type=="Grab HR"|Sample.Type=="Swab")
CTHR<-PreHarvestAPC%>%filter(Sample.Type=="Grab HR"| Sample.Type=="Grab")
ASCT<-PreHarvestAPC%>%filter(Sample.Type=="Swab"| Sample.Type=="Grab")
AS<-PreHarvestAPC%>%filter(Sample.Type=="Swab")
CT<-PreHarvestAPC%>%filter(Sample.Type=="Grab HR")
HR<-PreHarvestAPC%>%filter(Sample.Type=="Grab")
cASHR<-PreHarvestColif%>%filter(Sample.Type=="Grab HR"|Sample.Type=="Swab")
cCTHR<-PreHarvestColif%>%filter(Sample.Type=="Grab HR"| Sample.Type=="Grab")
cASCT<-PreHarvestColif%>%filter(Sample.Type=="Swab"| Sample.Type=="Grab")
cAS<-PreHarvestColif%>%filter(Sample.Type=="Swab")
cCT<-PreHarvestColif%>%filter(Sample.Type=="Grab HR")
cHR<-PreHarvestColif%>%filter(Sample.Type=="Grab")

#Check Assumptions
APClm<-lm(Best.Estimate ~ Sample.Type, data=PreHarvestAPC)
locAPClm<-lm(Best.Estimate ~ Sample.Type*Area, data=PreHarvestAPC)
shapiro.test(residuals(APClm))##normatility p=0.11
leveneTest(Best.Estimate~Sample.Type, data=PreHarvestAPC)#heterogenous p=0.00012
#leveneTest(APClm)
anova(APClm) #effect in sample type p<2.2e-16

#Variance test
#var.test(Best.Estimate ~ Sample.Type, data=PreHarvestAPC)##just for two groups works
bartlett.test(Best.Estimate~Sample.Type, data=PreHarvestAPC)
var.test(Best.Estimate~Sample.Type, data=ASHR)#var of AS:0.2899263; CT:0.1031136; HR:0.0706667
var.test(Best.Estimate~Sample.Type, data=CTHR)
var.test(Best.Estimate~Sample.Type, data=ASCT)


#Boxcox Transf = lambda: -1.434343 - for both anovas
tAPClm<-boxcox(Best.Estimate~Sample.Type, data=PreHarvestAPC)
lambda<-tAPClm$x[which.max(tAPClm$y)]

tLocAPClm<-boxcox(Best.Estimate~Sample.Type*Area, data=PreHarvestAPC)
lambdaloc<-tLocAPClm$x[which.max(tLocAPClm$y)]

newAPClm<-lm(((Best.Estimate^lambda-1)/lambda)~Sample.Type, data= PreHarvestAPC)
shapiro.test(residuals(newAPClm))
leveneTest(newAPClm)
newanova<-anova(newAPClm)#significant effect of sample type p<2.2e-16

#transformation for twoway anova - For location
newlocAPClm<-lm(((Best.Estimate^lambda-1)/lambda)~Sample.Type*Area, data= PreHarvestAPC)
shapiro.test(residuals(newlocAPClm))
leveneTest(newlocAPClm)
newLOCanova<-anova(newlocAPClm)
APCaov<-aov(((Best.Estimate^lambda-1)/lambda)~Sample.Type, data= PreHarvestAPC)
APCdiff<-TukeyHSD(APCaov)
APCdiff2<-TukeyHSD(testanova)
plot(APCdiff) #Check for intervals in Tukey HSD

#SimTestDiff(data=PreHarvestAPC, grp= "Sample.Type", resp = "Best.Estimate", type= "Tukey", covar.equal = TRUE)

Coliflm<-lm(Best.Estimate ~ Sample.Type, data=PreHarvestColif)
shapiro.test(residuals(Coliflm))##normatility p=0.054
leveneTest(Best.Estimate~Sample.Type, data=PreHarvestColif) #heterogenous p=1.59e-9
anova(Coliflm)#effect on sample type p<1.96e-5
#BOxcox Transformation 
tColiflm<-boxcox(Best.Estimate~Sample.Type, data=PreHarvestColif)
lambda2<-tColiflm$x[which.max(tColiflm$y)]#lambda for coliforms: 0.3838383838...4
newColiflm<-lm(((Best.Estimate^lambda2-1)/lambda2)~Sample.Type, data= PreHarvestColif)
shapiro.test(residuals(newColiflm))
leveneTest(newColiflm)
newanovaColif<-anova(newColiflm)

# ColifWilx<- wilcox.test(Best.Estimate~Sample.Type, data= PreHarvestColif)
# ColifWilx

ColifKrusk<-kruskal.test(Best.Estimate~Sample.Type, data= PreHarvestColif)
ColifKrusk
posthocKrusk<-dunnTest(Best.Estimate~Sample.Type, data= PreHarvestColif)


#Variance
var.test(Best.Estimate~Sample.Type, data=cASHR)#var of cAS:0.1310509; cCT:0.8742782; cHR:1.667079
var.test(Best.Estimate~Sample.Type, data=cCTHR)
var.test(Best.Estimate~Sample.Type, data=cASCT)

# Edge<-PreHarvest%>%filter(Field.Loc=="Edge")
# EdgeAPC<-Edge%>%filter(Test=="APC")


# Bed8<-
# 
# Bed1<-
#########################
##In- and post-harvest
InPostHarvest<- read.csv("In-PostHarvest Trimmed.csv")

##Subsetting for in harvest
InHarvest<- InPostHarvest%>%filter(Sampling.Day=="InHarvest D1"|Sampling.Day=="InHarvest D2")

InHarvestAPC<- InHarvest%>%filter(Test=="APC")
InHarvestColi<- InHarvest%>%filter(Test=="Coliforms")

#On the ground (glove vs leftover trims)
InHarvestGroundAPC<- InHarvestAPC%>%filter(Location=="Leftover"|Location=="Gloves")
InHarvestGroundColi<- InHarvestColi%>%filter(Location=="Leftover"|Location=="Gloves")

#Stats on the ground (APC)
#Check Assumptions
InHarvestGroundAPC_lm<-lm(Best.Estimate.log.CFU.g. ~ Sample.Type, data=InHarvestGroundAPC)
shapiro.test(residuals(InHarvestGroundAPC_lm))##normatility p=0.34
leveneTest(InHarvestGroundAPC_lm)#heterogenous p=0.94
#Anova
anova(InHarvestGroundAPC_lm) #p=3.47e-05
#Variance
var.test(Best.Estimate.log.CFU.g. ~ Sample.Type, data=InHarvestGroundAPC) #p=0.91

#Stats on the ground (Coliforms)
#Check Assumptions
InHarvestGroundColi_lm<-lm(Best.Estimate.log.CFU.g. ~ Sample.Type, data=InHarvestGroundColi)
shapiro.test(residuals(InHarvestGroundColi_lm))##normatility p=0.05
leveneTest(InHarvestGroundColi_lm)#heterogenous p=0.04
#Anova
anova(InHarvestGroundColi_lm) #p=0.94
#Variance
var.test(Best.Estimate.log.CFU.g. ~ Sample.Type, data=InHarvestGroundColi) #p=0.0023

#On the trailer
TrailerBinsAPC<- InHarvestAPC%>%filter(Location=="Bins")

#######################################################
#For prelim analysis

#Removal of Grab HR
ASCT_APClm<-lm(Best.Estimate ~ Sample.Type, data=ASCT)
ASCT_locAPClm<-lm(Best.Estimate ~ Sample.Type*Area, data=ASCT)
shapiro.test(residuals(ASCT_APClm))##normatility p=0.9038
leveneTest(Best.Estimate~Sample.Type, data=ASCT)#heterogeneous 0.00058

#Boxcox Transf = ASCT_lambda: -1.7979798/ ASCT_lambdaloc= -1.7171717172
t_ASCT_APClm<-boxcox(Best.Estimate~Sample.Type, data=ASCT)
ASCT_lambda<-t_ASCT_APClm$x[which.max(t_ASCT_APClm$y)]

t_ASCT_LocAPClm<-boxcox(Best.Estimate~Sample.Type*Area, data=ASCT)
ASCT_lambdaloc<-t_ASCT_LocAPClm$x[which.max(t_ASCT_LocAPClm$y)]

ASCT_newAPClm<-lm(((Best.Estimate^ASCT_lambda-1)/ASCT_lambda)~Sample.Type, data= ASCT)
shapiro.test(residuals(ASCT_newAPClm))
leveneTest(ASCT_newAPClm)#homogeneous variance p=0.51
ASCT_newanova<-anova(ASCT_newAPClm)#significant effect of sample type p<1.447e-15

ASCT_locnewAPClm<-lm(((Best.Estimate^ASCT_lambdaloc-1)/ASCT_lambdaloc)~Sample.Type+Area+Class+Sample.Type*Area+Sample.Type*Class+Area*Class, data= ASCT)
shapiro.test(residuals(ASCT_locnewAPClm))
leveneTest(ASCT_locnewAPClm)#homogeneous variance p=0.51
ASCT_locnewanova<-anova(ASCT_locnewAPClm)
ASCT_locAOV<-aov(((Best.Estimate^ASCT_lambdaloc-1)/ASCT_lambdaloc)~Sample.Type+Area+Class+Sample.Type*Area+Sample.Type*Class+Area*Class, data= ASCT)
APCdiff3<-TukeyHSD(ASCT_locAOV)
##

##Coliforms
cASCT_CClm<-lm(Best.Estimate ~ Sample.Type, data=cASCT)
cASCT_locCClm<-lm(Best.Estimate ~ Sample.Type*Area, data=cASCT)
shapiro.test(residuals(cASCT_CClm))##normatility p=0.9038
leveneTest(Best.Estimate~Sample.Type, data=cASCT)#heterogeneous 0.00058

#Boxcox Transf = cASCT_lambda: 0.10101010/ cASCT_lambdaloc= 0.181818182
t_cASCT_CClm<-boxcox(Best.Estimate~Sample.Type, data=cASCT)
cASCT_lambda<-t_cASCT_CClm$x[which.max(t_cASCT_CClm$y)]

t_cASCT_LocCClm<-boxcox(Best.Estimate~Sample.Type*Area, data=cASCT)
cASCT_lambdaloc<-t_cASCT_LocCClm$x[which.max(t_cASCT_LocCClm$y)]

cASCT_newCClm<-lm(((Best.Estimate^cASCT_lambda-1)/cASCT_lambda)~Sample.Type, data= cASCT)
shapiro.test(residuals(cASCT_newCClm))
leveneTest(cASCT_newCClm)#homogeneous variance p=0.51
cASCT_newanova<-anova(cASCT_newCClm)#significant effect of sample type p<1.447e-15

cASCT_locnewCClm<-lm(((Best.Estimate^cASCT_lambdaloc-1)/cASCT_lambdaloc)~Sample.Type+Area+Class+Sample.Type*Area+Sample.Type*Class+Area*Class, data= cASCT)
shapiro.test(residuals(cASCT_locnewCClm))#non normal p=0.0034
leveneTest(cASCT_locnewCClm)#Model must be completely crossed formula only
ASCT_locnewanova<-anova(ASCT_locnewAPClm)
cASCT_locAOV<-aov(((Best.Estimate^cASCT_lambdaloc-1)/cASCT_lambdaloc)~Sample.Type+Area+Class+Sample.Type*Area+Sample.Type*Class+Area*Class, data= cASCT)
CCdiff4<-TukeyHSD(cASCT_locAOV)

 


