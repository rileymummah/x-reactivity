## ---------------------------
## This code was written by: acr gomez & ro mummah
## For questions: rom5173@ucla.edu
## ---------------------------

## ---------------------------
## Objective: 
## Host-specific patterns of relative MAT antibody titers detected against 
## five Leptospira serovars (Pomona, Djasiman, Autumnalis, Bratislava, and 
## Icterohaemorrhagiae) when the infecting serovar is L. interrogans serovar 
## Pomona
## 
## Input:
##   Foxes_shedding_or_pfge_without_duplicates.csv
##   Skunks_pfge_or_shedding.csv
##   CSL_shedding_or_pfge.csv
##
## Output: 
##   Fig1
##
## ---------------------------

## load packages ---------------------------
library(tidyverse)
library(fmsb)
library(viridis)
library(lubridate)


## load data ---------------------------
options(stringsAsFactors = FALSE)
foxes <- read.csv("data/Foxes_shedding_or_pfge_without_duplicates.csv")
skunks <- read.csv("data/Skunks_pfge_or_shedding.csv")
csl <- read.csv("data/CSL_shedding_or_pfge.csv")


## Set colors for plotting ---------------------------
colorFoxline  <- viridis(6,0.05,option="D")[4]
colorFox  <- viridis(9,0.025,option='D')[5]

colorCSLline <-  viridis(1,0.05,option="D")
colorCSL <-  viridis(9,0.03,option='D')[2]

colorSkunk <-  viridis(9,0.2,option='D')[8]



fullnames <- c("Pom","Ict","Bra","Aut","Dja")


csl <- transmute_at(csl, 2:6, as.numeric)
foxes <- transmute_at(foxes, 2:6, as.numeric)
skunks <- transmute_at(skunks, 2:6, as.numeric)

numcsl <- dim(csl)[1]
numfoxes <- dim(foxes)[1]
numskunks <- dim(skunks)[1]

relcsl <- csl/apply(csl, FUN=max, MARGIN = 1)
relfoxes <- foxes/apply(foxes, FUN=max, MARGIN = 1)
relskunks <- skunks/apply(skunks, FUN=max, MARGIN = 1)


####
meancolcsl <- colSums(relcsl)/dim(relcsl)[1]
cslsample <- relcsl[sample(seq_along(relcsl$LogPom), numfoxes),] #subsampling with same number of samples as foxes
meancolfox <- colSums(relfoxes)/dim(relfoxes)[1]
meansku <- colSums(relskunks)/dim(relskunks)[1]

# percent max -------------------------------------------------------------

cslmaxpercent <- csl*NA
for (ii in seq_along(csl[,1])){
  cslmaxpercent[ii,] <- csl[ii,] == max(csl[ii,])
}
cslmaxp <- colSums(cslmaxpercent)/length(csl[,1])

foxesmaxpercent <- foxes*NA
for (ii in seq_along(foxes[,1])){
  foxesmaxpercent[ii,] <- foxes[ii,] == max(foxes[ii,])
}
foxesmaxp <- colSums(foxesmaxpercent)/length(foxes[,1])


skunksmaxpercent <- skunks*NA
for (ii in seq_along(skunks[,1])){
  skunksmaxpercent[ii,] <- skunks[ii,] == max(skunks[ii,])
}
skunksmaxp <- colSums(skunksmaxpercent)/length(skunks[,1])


######
minicsl <- rbind(1, 0, cslsample, meancolcsl, cslmaxp) #first row is max, second row is min,last 2 rows is average
minifox <- rbind(1, 0, relfoxes, meancolfox, foxesmaxp) #average and total
minisku <- rbind(1, 0, relskunks, meansku, skunksmaxp)

foxtitle <- paste(fullnames," \n",round(foxesmaxp*100,digits = 1),"%",sep="")
csltitle <- paste(fullnames," \n",round(cslmaxp*100,digits = 1),"%",sep="")
skunktitle <- paste(fullnames," \n",round(skunksmaxp*100,digits = 1),"%",sep="")

# dashed: proportion of samples for which this is the max titer"
# continuous: average titer/(max titer) across all samples for this serovar"


png("figures/fig1.png",width = 1740,height = 650)
par(mfrow=c(1,3))
radarchart( minicsl,title="Sea Lions", axistype=1 , pty=32, seg=5,#centerzero = T,
            vlabels=csltitle, vlcex=3, cex.main=5,
            plty=c(rep(1,numfoxes+1),2), 
            plwd=2,
            pcol=c(rep(colorCSLline,numfoxes),"black","black"),
            pfcol=c(rep(colorCSL,numfoxes),"#00000000","#00000000"),
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="transparent", 
            calcex=2,
            caxislabels=seq(0,1,.2))

radarchart( minifox, title= "Foxes",axistype=1 , pty=32,seg=5,#centerzero = T,
            vlabels=foxtitle,cex.main=5,
            plty=c(rep(1,numfoxes+1),2),plwd=2,
            pcol=c(rep(colorFoxline,numfoxes),"black","black"),
            pfcol=c(rep(colorFox,numfoxes),"#00000000","#00000000"),
            cglcol="grey", cglty=1, axislabcol="transparent", 
            caxislabels=seq(0,1,.2), #cglwd=0.8,
            #custom labels
            calcex=2,
            vlcex=3)
            

radarchart( minisku, title="Skunks",axistype=1 , pty=32, seg=5,#centerzero = T,
            vlabels=skunktitle, cex.main=5,
            plty=c(rep(1,numskunks+1),2), plwd=2,
            pcol=c(rep("#00000000",numskunks),"black","black"),
            pfcol=c(rep(colorSkunk,numskunks),"#00000000","#00000000"),
            caxislabels=seq(0,1,.2), #cglwd=0.8,
            cglcol="grey", cglty=1, axislabcol="transparent", 
            #custom labels
            calcex=2,
            vlcex=3)

dev.off()

# End script