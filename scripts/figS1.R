## ---------------------------
## This code was written by: r.o. mummah
## For questions: rmummah@umass.edu
## Date Created: 2023-07-07
## ---------------------------

## ---------------------------
## Objective: 
##
## 
## Input:
##   
##
## Output: 
##
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
csl.pfge <- read.table("selected samples/CSL_pfge_confirmed_only.csv", 
                       header=T, sep=',') %>% # PFGE
  mutate(Species = 'CSL')

# Foxes
fox.pfge <- read.table("selected samples/Foxes_pfge_only_noduplicates.csv", 
                       header=T, sep=',') %>% # PFGE only
  # read.table("selected samples/Foxes_shedding_or_pfge_without_duplicates.csv", 
  #                      header=T, sep=',') %>% # PFGE + PCR+
  mutate(Species = 'fox') #%>%
  # # Found in PFGE_shedding_IDS_foxes_skunks.csv
  # filter(ID %in% c('07180','23850','26195','32256','36401','46621',
  #                  '84D5D','86076','B0F3B','C0D60','E6D47'))


# Skunks
# Only 1 PFGE skunk, so we use all of them for both datasets
skunks <- read.csv("selected samples/Skunks_pfge_or_shedding.csv", 
                   header=T, sep=',')[-5,] %>%
  mutate(Species = 'skunk')


colorFoxline  <- viridis(6,0.05,option="D")[4]
colorFox  <- viridis(6,0.03,option="D")[4]
colorCSLline <-  viridis(1,0.05,option="D")
colorCSL <-  viridis(1,0.03,option="D")
colorSkunk <-  viridis(6,0.5,option="D")[6]


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
cslsample <- relcsl[sample(seq_along(relcsl$LogPom),numfoxes),] #subsampling with same number of samples as foxes
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
minicsl=rbind(1 , 0 ,cslsample,meancolcsl,cslmaxp) #first row is max, second row is min,last 2 rows is average
minifox <- rbind(1 , 0 , relfoxes, meancolfox,foxesmaxp) #average and total
minisku <- rbind(1,0,relskunks,meansku,skunksmaxp)

foxtitle <- paste(fullnames," \n",round(foxesmaxp*100,digits = 2),"%",sep="")
csltitle <- paste(fullnames," \n",round(cslmaxp*100,digits = 2),"%",sep="")
skunktitle <- paste(fullnames," \n",round(skunksmaxp*100,digits = 2),"%",sep="")




png("figures/figS1.png", width = 1740, height = 650)
par(mfrow=c(1,3))

radarchart( minicsl, title="Sea Lions", axistype=1, pty=32, seg=5,#centerzero = T,
            vlabels=csltitle, vlcex=2.5, cex.main=5,
            plty=c(rep(1,numfoxes+1),2), 
            plwd=2,
            # title="dashed: proportion of samples for which this is the max titer",
            pcol=c(rep(colorCSLline,numfoxes),"black","black"),
            pfcol=c(rep(colorCSL,numfoxes),"#00000000","#00000000"),
            #pfcol=colors_in , plwd=1 ,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="transparent", 
            calcex = 2,
            caxislabels=seq(0,1,.2))

radarchart( minifox, title="Foxes", axistype=1, pty=32, seg=5,#centerzero = T,
            vlabels=foxtitle, vlcex=2.5, cex.main=5,
            plty=c(rep(1,numfoxes+1),2),
            plwd=2,
            pcol=c(rep(colorFoxline,numfoxes),"black","black"),
            pfcol=c(rep(colorFox,numfoxes),"#00000000","#00000000"),
            cglcol="grey", cglty=1, axislabcol="transparent", 
            caxislabels=seq(0,1,.2))
# title="continuous: average titer/(max titer) across all samples for this serovar")

radarchart( minisku, title="Skunks", axistype=1, pty=32, seg=5,#centerzero = T,
            vlabels=skunktitle, vlcex=2.5, cex.main=5,
            plty=c(rep(1,numskunks+1),2),
            plwd=2,
            pcol=c(rep("#00000000",numskunks),"black","black"),
            pfcol=c(rep(colorSkunk,numskunks),"#00000000","#00000000"),
            caxislabels=seq(0,1,.2), 
            cglcol="grey", cglty=1, axislabcol="transparent") 

dev.off()

