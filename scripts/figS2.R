## ---------------------------
## This code was written by: r.o. mummah
## For questions: rom5173@ucla.edu
## ---------------------------

## ---------------------------
## Objective: Create FigS1
## Host-specific patterns of relative MAT antibody titers detected against five 
## Leptospira serovars (Pomona, Djasiman, Autumnalis, Bratislava, and 
## Icterohaemorrhagiae) when the infecting serovar is L. interrogans serovar 
## Pomona for individuals with positive PFGE only
## 
## Input:
##   CSL_pfge_confirmed_only.csv
##   Foxes_pfge_only_noduplicates.csv
##   Skunks_pfge_or_shedding.csv
##
## Output: 
##   FigS2
##
## ---------------------------

## load packages ---------------------------
library(tidyverse)
library(magrittr)
library(fmsb)
library(viridis)
library(lubridate)

## load functions ---------------------------
convertMAT <- function(x){
  x <- as.numeric(x)
  titers <- c(0,2^(0:15)*100)
  logtiters <- 0:16
  if (x %in% titers) {
    return(logtiters[titers==x])
  } else {
    return(NA)
  }
}

## load data ---------------------------

# CSL
csl.pfge <- read.table("data/CSL_pfge_confirmed_only.csv", 
                       header=T, sep=',') %>% # PFGE
              mutate(Species = 'CSL')

# Foxes
fox.pfge <- read.table("data/Foxes_pfge_only_noduplicates.csv", 
                       header=T, sep=',') %>% # PFGE only
            mutate(Species = 'fox')


# Skunks
# Only 1 PFGE skunk, so we use all of them for both datasets
skunks <- read.csv("data/Skunks_pfge_or_shedding.csv", 
                   header=T, sep=',') %>%
            mutate(Species = 'skunk')


colorFoxline  <- viridis(6,0.05,option="D")[4]
colorFox  <- viridis(9,0.025,option='D')[5]

colorCSLline <-  viridis(1,0.05,option="D")
colorCSL <-  viridis(9,0.03,option='D')[2]

colorSkunk <-  viridis(9,0.2,option='D')[8]



fullnames <- c("Pom","Ict","Bra","Aut","Dja")


csl.pfge <- transmute_at(csl.pfge, 2:6, as.numeric)
fox.pfge <- transmute_at(fox.pfge, 2:6, as.numeric)
skunks <- transmute_at(skunks, 2:6, as.numeric)

numcsl <- dim(csl.pfge)[1]
numfoxes <- dim(fox.pfge)[1]
numskunks <- dim(skunks)[1]

relcsl <- csl.pfge/apply(csl.pfge, FUN=max, MARGIN = 1)
relfoxes <- fox.pfge/apply(fox.pfge, FUN=max, MARGIN = 1)
relskunks <- skunks/apply(skunks, FUN=max, MARGIN = 1)


####
meancolcsl=colSums(relcsl)/dim(relcsl)[1]
# cslsample <- relcsl[sample(seq_along(relcsl$LogPom),numfoxes),] #subsampling with same number of samples as foxes
meancolfox=colSums(relfoxes)/dim(relfoxes)[1]
meansku = colSums(relskunks)/dim(relskunks)[1]

#####percent max
cslmaxpercent <- csl.pfge*NA
for (ii in seq_along(csl.pfge[,1])){
  cslmaxpercent[ii,] <- csl.pfge[ii,] == max(csl.pfge[ii,])
}
cslmaxp <- colSums(cslmaxpercent)/length(csl.pfge[,1])

foxesmaxpercent <- fox.pfge*NA
for (ii in seq_along(fox.pfge[,1])){
  foxesmaxpercent[ii,] <- fox.pfge[ii,] == max(fox.pfge[ii,])
}
foxesmaxp <- colSums(foxesmaxpercent)/length(fox.pfge[,1])


skunksmaxpercent <- skunks*NA
for (ii in seq_along(skunks[,1])){
  skunksmaxpercent[ii,] <- skunks[ii,] == max(skunks[ii,])
}
skunksmaxp <- colSums(skunksmaxpercent)/length(skunks[,1])


######
minicsl <- rbind(1, 0, relcsl,meancolcsl,cslmaxp) #first row is max, second row is min,last 2 rows is average
minifox <- rbind(1 , 0 , relfoxes, meancolfox,foxesmaxp) #average and total
minisku <- rbind(1,0,relskunks,meansku,skunksmaxp)

foxtitle <- paste(fullnames," \n",round(foxesmaxp*100,digits = 2),"%",sep="")
csltitle <- paste(fullnames," \n",round(cslmaxp*100,digits = 2),"%",sep="")
skunktitle <- paste(fullnames," \n",round(skunksmaxp*100,digits = 2),"%",sep="")




png("figures/figS2.png", width = 1740, height = 650)
par(mfrow=c(1,3))

radarchart( minicsl, title="Sea Lions", axistype=1, pty=32, seg=5,#centerzero = T,
            vlabels=csltitle, vlcex=2.5, cex.main=5,
            plty=c(rep(1,numcsl+1),2), 
            plwd=2,
            pcol=c(rep(colorCSLline,numcsl),"black","black"),
            pfcol=c(rep(colorCSL,numcsl),"#00000000","#00000000"),
            cglcol="grey", cglty=1, axislabcol="transparent", 
            calcex = 2,
            caxislabels=seq(0,1,.2))

radarchart( minifox, title="Foxes", axistype=1, pty=32, seg=5,
            vlabels=foxtitle, vlcex=2.5, cex.main=5,
            plty=c(rep(1,numfoxes+1),2),
            plwd=2,
            pcol=c(rep(colorFoxline,numfoxes),"black","black"),
            pfcol=c(rep(colorFox,numfoxes),"#00000000","#00000000"),
            cglcol="grey", cglty=1, axislabcol="transparent", 
            caxislabels=seq(0,1,.2))

radarchart( minisku, title="Skunks", axistype=1, pty=32, seg=5,
            vlabels=skunktitle, vlcex=2.5, cex.main=5,
            plty=c(rep(1,numskunks+1),2),
            plwd=2,
            pcol=c(rep("#00000000",numskunks),"black","black"),
            pfcol=c(rep(colorSkunk,numskunks),"#00000000","#00000000"),
            caxislabels=seq(0,1,.2), 
            cglcol="grey", cglty=1, axislabcol="transparent") 

dev.off()


# End script

