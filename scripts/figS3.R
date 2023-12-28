## ---------------------------
## This code was written by: r.o. mummah
## For questions: rmummah@usgs.gov
## Date Created: 2023-08-08
## ---------------------------

## ---------------------------
## Objective: Extract necessary data for plotting & analysis
## Patterns of relative MAT antibody titers detected against five Leptospira 
## serovars when the infecting serovar is L. interrogans serovar Pomona for 
## MAT-positive skunks and PCR-positive skunks.
## 
## Input:
##   Skunks_matpos_nodup.csv
##   Skunks_pfge_or_shedding.csv
## Output: 
##   FigS3
##
## ---------------------------

## load packages ---------------------------
library(tidyverse)
library(fmsb)
library(viridis)
library(lubridate)


## load data ---------------------------
options(stringsAsFactors = FALSE)
skunk.mat <- read.csv("data/Skunks_matpos_nodup.csv")
skunk.pcr <- read.csv("data/Skunks_pfge_or_shedding.csv")


## Set colors for plotting ---------------------------
# colorFoxline  <- viridis(6,0.05,option="D")[4]
# colorFox  <- viridis(9,0.025,option='D')[5]#viridis(6,0.025,option="D")[4]
# # colorFoxopaque  <- viridis(6,option="D")[4]
# colorCSLline <-  viridis(1,0.05,option="D")
# colorCSL <-  viridis(9,0.03,option='D')[2]#viridis(1,0.03,option="D")
# colorCSLopaque <- viridis(1,option="D")
colorSkunk <-  viridis(9,0.2,option='D')[8]#viridis(6,0.2,option="D")[6]


fullnames <- c("Pom","Ict","Bra","Aut","Dja")

skunk.mat <- transmute_at(skunk.mat, 2:6, as.numeric)
skunk.pcr <- transmute_at(skunk.pcr, 2:6, as.numeric)

numskunk.mat <- dim(skunk.mat)[1]
numskunk.pcr <- dim(skunk.pcr)[1]

relskunk.mat <- skunk.mat/apply(skunk.mat, FUN=max, MARGIN = 1)
relskunk.pcr <- skunk.pcr/apply(skunk.pcr, FUN=max, MARGIN = 1)


####
meansku.mat <- colSums(relskunk.mat)/dim(relskunk.mat)[1]
meansku.pcr <- colSums(relskunk.pcr)/dim(relskunk.pcr)[1]

# percent max -------------------------------------------------------------
# START HERE WITH REPLACEMENT TO TWO SKUNK FILES

maxpercentsku.mat <- skunk.mat*NA
for (ii in seq_along(skunk.mat[,1])){
  maxpercentsku.mat[ii,] <- skunk.mat[ii,] == max(skunk.mat[ii,])
}
maxp.skunk.mat <- colSums(maxpercentsku.mat)/length(skunk.mat[,1])


maxpercentsku.pcr <- skunk.pcr*NA
for (ii in seq_along(skunk.pcr[,1])){
  maxpercentsku.pcr[ii,] <- skunk.pcr[ii,] == max(skunk.pcr[ii,])
}
maxp.skunk.pcr <- colSums(maxpercentsku.pcr)/length(skunk.pcr[,1])


######
minisku.mat <- rbind(1, 0, relskunk.mat, meansku.mat, maxp.skunk.mat)
minisku.pcr <- rbind(1, 0, relskunk.pcr, meansku.pcr, maxp.skunk.pcr)

title.sku.mat <- paste(fullnames," \n",round(maxp.skunk.mat*100,digits = 1),"%",sep="")
title.sku.pcr <- paste(fullnames," \n",round(maxp.skunk.pcr*100,digits = 1),"%",sep="")

# dashed: proportion of samples for which this is the max titer"
# continuous: average titer/(max titer) across all samples for this serovar"


png("figures/figS3.png",width = 1160,height = 650)
par(mfrow=c(1,2), mai = c(0.1, 0.1, 1, 0.1))
radarchart( minisku.mat, title="Skunk MAT+",
            axistype=1, pty=32, seg=5,#centerzero = T,
            vlabels=title.sku.mat, cex.main=4,
            plty=c(rep(1, numskunk.mat+1),2), plwd=2,
            pcol=c(rep("#00000000", numskunk.mat),"black","black"),
            pfcol=c(rep(colorSkunk, numskunk.mat),"#00000000","#00000000"),
            caxislabels=seq(0,1,.2), #cglwd=0.8,
            cglcol="grey", cglty=1, axislabcol="transparent", 
            #custom labels
            calcex=2,
            vlcex=2)


radarchart( minisku.pcr, title="Skunk PCR+",
            axistype=1 , pty=32, seg=5,#centerzero = T,
            vlabels=title.sku.pcr, cex.main=4,
            plty=c(rep(1,numskunk.pcr+1),2), plwd=2,
            pcol=c(rep("#00000000",numskunk.pcr),"black","black"),
            pfcol=c(rep(colorSkunk,numskunk.pcr),"#00000000","#00000000"),
            caxislabels=seq(0,1,.2), #cglwd=0.8,
            cglcol="grey", cglty=1, axislabcol="transparent", 
            #custom labels
            calcex=2,
            vlcex=2)

dev.off()

# End script
