## ---------------------------
## This code was written by: acr gomez & ro mummah
## For questions: rom5173@ucla.edu
## ---------------------------

## ---------------------------
## Objective: Create Fig3
## Selected example of longitudinal antibody titer dynamics in a 
## Channel Island fox.
## 
## Input:
##   fox_02107_only.csv
##
## Output: 
##   Fig3
##
## ---------------------------

## load packages ---------------------------
library(tidyverse)
library(lubridate)
library(patchwork)

## load functions ---------------------------
convertMAT <- function(x){
  x <- as.numeric(x)
  titers <- c(0, 2^(0:15)*100)
  logtiters <- 0:16
  ii=match(x,titers)
  return(logtiters[ii])
}

## load data ---------------------------
fox02107 <- read.csv("data/fox_02107.csv") %>%
              mutate(Date = ymd(Date))


## ---------------------------

pcr <- filter(fox02107, Test=="PCR")
mat <- filter(fox02107, Test=="MAT")

leg <- data.frame(lab = c("Pomona","Autumnalis"),
                  x = c(as.Date("2015-06-01"), as.Date("2015-07-15")),
                  xdot = c(as.Date("2014-08-15"), as.Date("2014-08-15")),
                  color = c("black","cornflowerblue"),
                  lty = c('solid','dashed'),
                  y = 9:8)

top <- ggplot(mat) +
        theme_classic(base_size = 14) +
        ylab(expression(paste(log[2],' MAT titer', sep=""))) +
        xlab("") +
        scale_x_date(labels=NULL, limits = c(as.Date("2011-08-01"), as.Date("2016-01-01"))) +
        scale_y_continuous(breaks=seq(0,8,2), labels=seq(0,8,2)) +
        geom_line(aes(x = Date, y = as.numeric(PomonaLevel)), 
                      col = 'black') +
        geom_point(aes(x = Date, y = as.numeric(PomonaLevel)), 
                   col = 'black') +
        geom_line(aes(x = Date, y = as.numeric(AutLevel)), 
                  color="cornflowerblue", lty='dashed') +
        geom_point(aes(x = Date, y = as.numeric(AutLevel)), 
                   color="cornflowerblue", pch=18, cex=2.5) +
        geom_text(data = leg, aes(x=x,y=y,label=lab)) +
        geom_point(data = leg, aes(x=xdot, y=y), 
                   col=leg$color, pch=c(19,18), cex=c(1.5,2.5)) +
        geom_segment(data = leg, aes(x=xdot-80, xend=xdot+80, y=y, yend=y), 
                     col=leg$color, lty=leg$lty)
        

bot <- ggplot(data=pcr) +
        theme_classic(base_size = 14)+
        scale_x_date(limits = c(as.Date("2011-08-01"), as.Date("2016-01-01"))) +
        geom_point(aes(x=Date, y=1), pch=c("+","-","-"), cex=6) +
        scale_y_continuous(breaks=0, limits=c(0.75,1.25), labels="") +
        ylab("PCR") +
        xlab("Date Sampled")

png("figures/fig3.png", width = 4,height = 3.5,res=300, units = "in")
    top / bot + plot_layout(ncol=1, heights=c(3,1))
dev.off()


# End script
