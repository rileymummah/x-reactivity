library(readr)
library(tidyverse)
library(ggplot2)

##Confirmed Shedders Titer Decline Plots
setwd()
DF <- read.csv("SheddersTiterDecline.csv")

new_Labels <- c("02107"="1","07478"="2","14A0D"="3","23727"="4","23850"="5","24269"="6","37531"="7","37539"="8","39031"="9","39037"="10","59843"="11","61904"="12","78624"="13","86076"="14","90599"="15","90697"="16","98205"="17","B0F3B"="18","C0D60"="19","F4A38"="20")

plt1 <- ggplot(DF[!is.na(DF$PomonaLevel),],aes(DaysSinceKnownPos)) + 
  geom_line(aes(y=PomonaLevel,colour="Pomona")) + 
  geom_point(aes(y=PomonaLevel,colour="Pomona")) + 
  geom_line(aes(y=AutumnalisLevel,colour="Autumnalis")) + 
  geom_point(aes(y=AutumnalisLevel,colour="Autumnalis")) + 
  geom_text(aes(label=CurrentShedPN,x=DaysSinceKnownPos,y=-2)) + 
  geom_text(aes(label= "PCR",x=-500,y=-2))  + 
  geom_text(aes(label= "MAT",x=-500,y=10)) + 
  geom_hline(aes(yintercept=-0.8),linewidth=0.5) + 
  facet_wrap(vars(Pittag), labeller = labeller(Pittag=new_Labels)) +
  scale_y_continuous(limits=c(-3,12),breaks=c(0,2,4,6,8,10))+ 
  scale_color_manual(values=c("purple","black"),labels=c("Autumnalis Titer","Pomona Titer"),name="") + 
  theme_light() + labs(y="PCR Status | Titer Level", x="Days Since Known Positive", title = "Confirmed Shedders Titer Decline Plots") + 
  theme(plot.title = element_text(hjust = 0.5)  )

ggsave(plt1)


#### Figure 4 B
LabComparePom3 <- read.csv(LabComparePom3)

ggplot(data=LabComparePom3,aes(Lab,Pomona)) + geom_boxplot() + scale_y_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + labs(y="log MAT",x=NULL) + theme_classic(base_size = 14)

ggsave()

## Figure 5

LabCompareBoth3 <- read.csv(LabCompareBoth3)

ggplot(data=LabCompareBoth3,aes(Lab,Titer,color=Serovar)) + geom_boxplot()+ theme_classic(base_size = 14) + scale_color_manual(values=c("purple","black"))  + scale_y_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + labs(y="log MAT",x=NULL)

ggsave()


#### Figure 4? A
CDCCSUCornellCompare <- read.csv(CDCCSUCornellCompare)

library(viridis)
labcolors <- viridis(6,option="D")

##3 lab compare sep graphs pom vs aut
plt6 <- ggplot(CDCCSUCornellCompare,aes(x=CDCAutumnalis,y=CDCPomona,shape=factor(CDCAutPCH))) + geom_jitter(width = 0.2, height = 0.2, color=labcolors[4],size=2) + theme_minimal(base_size = 14) + theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + scale_y_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + geom_abline(intercpet=0, slope=1,size =.8) + labs(x="Aut", y="Pom",title="Lab A")+scale_shape_manual(values=c(0,16))+theme(legend.position = "none")

plt7 <- ggplot(CDCCSUCornellCompare,aes(x=CornellAutumnalis,y=CornellPomona)) + geom_jitter(width = 0.2, height = 0.2, color=labcolors[4],size=2) + theme_minimal(base_size = 14) + theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + scale_y_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + geom_abline(intercpet=0, slope=1,size=.8) + labs(x="Aut", y="Pom",title="Lab B")

plt8 <- ggplot(CDCCSUCornellCompare,aes(x=CSUAutumnalis,y=CSUPomona,shape=factor(CSUAutPCH))) + geom_jitter(width = 0.2, height = 0.2, color=labcolors[4],size=2) + theme_minimal(base_size = 14) + theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + scale_y_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + geom_abline(intercpet=0, slope=1,size=.8) + labs(x="Aut", y="Pom",title="Lab C")+scale_shape_manual(values=c(1,16))+theme(legend.position = "none")

setwd("~/Dropbox/R/Leptospirosis/WorkingCombinedTestingFoxes/LabCompare")
jpeg("LabABCcomparePA_3plots.jpeg",res=300,height = 3, width=9, units = "in")

grid.arrange(plt6,plt7,plt8,nrow=1,ncol=3)

dev.off()