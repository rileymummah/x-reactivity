## ---------------------------
## This code was written by: r.o. mummah
## For questions: rmummah@umass.edu
## Date Created: 2023-07-12
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
library(GGally)
library(patchwork)

## load functions ---------------------------


## load data ---------------------------
options(stringsAsFactors = FALSE)
foxes <- read.csv("selected samples/Foxes_shedding_or_pfge_without_duplicates.csv")
skunks <- read.csv("selected samples/Skunks_pfge_or_shedding.csv")
csl <- read.csv("selected samples/CSL_shedding_or_pfge.csv")

## Set up ---------------------------

fullnames <- c("Pomona","Icterohaemorrhagiae","Bratislava","Autumnalis","Djasiman")


csl <- transmute_at(csl,2:6,as.numeric)
csl$Species ="Sea Lion"
foxes <- transmute_at(foxes,2:6,as.numeric)
foxes$Species = "Island Fox"
skunks <- transmute_at(skunks,2:6,as.numeric)
skunks$Species = "Skunk"


all <- rbind(csl,foxes,skunks) %>%
        reshape::melt(all) %>%
        mutate(Species = factor(Species, 
                                levels = c("Sea Lion","Island Fox","Skunk"),
                                labels = c("Sea Lion","Fox","Skunk"),
                                ordered=TRUE),
               both = paste(variable, Species)) %>%
        mutate(both = factor(both, both, both, ordered=T))


# Boxplot -----------------------------------------------------------------

colorFox  <- viridis(6,option="D")[4]
colorCSL <-  viridis(1,option="D")
colorSkunk <-  viridis(6,option="D")[6]

boxplot <- ggplot(all) +
            geom_boxplot(aes(both, value, color=Species), show.legend = F)+
            geom_rect(aes(ymin=0, ymax=15, 
                          xmin="LogIct Sea Lion", xmax="LogBra Sea Lion"), 
                      alpha =0.5, fill="grey90",
                      position = position_nudge(x=-.5)) +
            geom_rect(aes(ymin=0, ymax=15, 
                          xmin="LogAut Sea Lion", xmax="LogDja Sea Lion"), 
                      alpha =0.5, fill="grey90",
                      position = position_nudge(x=-.5)) +
            geom_boxplot(aes(both, value, color=Species), show.legend = F) +
            theme_bw(base_size=18) +
            theme(panel.grid = element_blank()) +
            ylab("log MAT") +
            xlab("") +
            theme(plot.margin = unit(c(1,1,1,1), "cm"), aspect.ratio =.5) +
            scale_y_continuous(limits = c(-0.5, 15), 
                               breaks = seq(0,15, by=2), 
                               labels = seq(0,15, by=2)) +
            scale_x_discrete(breaks = c("LogPom Fox","LogBra Fox",
                                        "LogDja Fox","LogIct Fox",
                                        "LogAut Fox"),
                             labels=c("Pom","Bra","Dja","Ict","Aut")) +
            scale_color_manual(values = c(colorCSL,colorFox,colorSkunk)) 


png("figures/fig2-boxplot.png", width=500, height = 300)
  boxplot
dev.off()

