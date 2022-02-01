##C:/Users/jorge/Box Sync/NIFA Project/CA 9-2021/Plate results/
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