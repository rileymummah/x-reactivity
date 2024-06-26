## ---------------------------
## This code was written by: ro mummah & acr gomez
## For questions: rom5173@ucla.edu
## ---------------------------

## ---------------------------
## Objective: Create the boxplot inset in Fig2
## Pairwise antibody titer levels against Leptospira interrogans serovars 
## Pomona, Djasiman, Autumnalis, Bratislava, and Icterohaemorrhagiae in three 
## host species
##
## Input:
##   Foxes_PCR-PFGE.csv
##   Skunks_PCR-PFGE.csv
##   CSL_PCR-PFGE.csv
##
## Output: 
##   Boxplot inset for Fig2
##
## ---------------------------

## load packages ---------------------------
library(tidyverse)
library(magrittr)
library(reshape)
library(fmsb)
library(viridis)
library(lubridate)
library(GGally)
library(patchwork)


## load data ---------------------------
options(stringsAsFactors = FALSE)
foxes <- read.csv("data/Foxes_PCR-PFGE.csv")
skunks <- read.csv("data/Skunks_PCR-PFGE.csv")
csl <- read.csv("data/CSL_PCR-PFGE.csv")

## Set up ---------------------------

fullnames <- c("Pomona","Icterohaemorrhagiae","Bratislava","Autumnalis","Djasiman")


csl <- transmute_at(csl,2:6,as.numeric)
csl$Species ="Sea Lion"
foxes <- transmute_at(foxes,2:6,as.numeric)
foxes$Species = "Island Fox"
skunks <- transmute_at(skunks,2:6,as.numeric)
skunks$Species = "Skunk"


all <- rbind(csl,foxes,skunks) %>% 
        reshape::melt(., id='Species') %>% 
        mutate(Species = factor(Species, 
                                levels = c("Sea Lion","Island Fox","Skunk"),
                                labels = c("Sea Lion","Fox","Skunk"),
                                ordered=TRUE),
               both = paste(variable, Species)) %>% 
        # arrange(variable) %>%
        mutate(both = factor(both, both, both, ordered=T))


# Boxplot -----------------------------------------------------------------

colorFox  <- viridis(6,option="D")[4]
colorCSL <-  viridis(1,option="D")
colorSkunk <- viridis(9, option="D")[8]

boxplot <- ggplot(all) +
            theme_bw(base_size=20) +
            theme(panel.grid = element_blank(),
                  plot.margin = unit(c(1,1,1,1), "cm"),
                  plot.background = element_blank(),
                  aspect.ratio = .5) +
            ylab(expression(paste(log[2]," MAT titer", sep=""))) +
            xlab("") +
            geom_boxplot(aes(both, value, color=Species), show.legend = F) +
            geom_rect(aes(ymin=0, ymax=15, 
                          xmin="LogIct Sea Lion", xmax="LogBra Sea Lion"), 
                      alpha =0.5, fill="grey90",
                      position = position_nudge(x=-.5)) +
            geom_rect(aes(ymin=0, ymax=15, 
                          xmin="LogAut Sea Lion", xmax="LogDja Sea Lion"), 
                      alpha =0.5, fill="grey90",
                      position = position_nudge(x=-.5)) +
            geom_boxplot(aes(both, value, color=Species), show.legend = F) +
            scale_y_continuous(limits = c(-0.5, 15), 
                               breaks = seq(0,15, by=2), 
                               labels = seq(0,15, by=2)) +
            scale_x_discrete(breaks = c("LogPom Fox","LogBra Fox",
                                        "LogDja Fox","LogIct Fox",
                                        "LogAut Fox"),
                             labels=c("Pom","Bra","Dja","Ict","Aut")) +
            scale_color_manual(values = c(colorCSL,colorFox,colorSkunk)) 


# tiff("figures/fig2-boxplot.tiff", width=500, height = 300)
#   boxplot
# dev.off()


# End script
