## ---------------------------
## This code was written by: r.o. mummah
## For questions: rmummah@umass.edu
## Date Created: 2023-07-07
## ---------------------------

## ---------------------------
## Objective: Create FigS2
## Pairwise antibody titer levels against Leptospira interrogans serovars 
## Pomona, Djasiman, Autumnalis, Bratislava, and Icterohaemorrhagiae in three 
## host species for individuals which are PFGE positive only
## 
## Input:
##   Foxes_pfge_only_noduplicates.csv
##   Skunks_pfge_or_shedding.csv
##   CSL_pfge_confirmed_only.csv
##
## Output: 
##   FigS2
##
## ---------------------------

## load packages ---------------------------
library(tidyverse)
library(fmsb)
library(viridis)
library(lubridate)
library(GGally)
library(patchwork)


## load data ---------------------------
options(stringsAsFactors = FALSE)
foxes <- read.csv("data/Foxes_pfge_only_noduplicates.csv")
skunks <- read.csv("data/Skunks_pfge_or_shedding.csv")
csl <- read.csv("data/CSL_pfge_confirmed_only.csv")


## Set up ---------------------------

fullnames <- c("Pomona","Icterohaemorrhagiae","Bratislava","Autumnalis","Djasiman")


csl <- transmute_at(csl,2:6,as.numeric)
csl$Species ="Sea Lion"
foxes <- transmute_at(foxes,2:6,as.numeric)
foxes$Species = "Island Fox"
skunks <- transmute_at(skunks,2:6,as.numeric)
skunks$Species = "Skunk"

all=rbind(csl,foxes,skunks)

# Universal plotting ------------------------------------------------------

colorCSL <- viridis(9, option="D")[2]
colorFox <- viridis(9, option="D")[5]
colorSkunk <- viridis(9, option="D")[8]

uni <- ggplot(all) + 
        geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1) + 
        scale_x_continuous(limits = c(-0.5, 15), 
                           breaks = seq(0,15, by=2)) +
        scale_y_continuous(limits = c(-0.5, 15), 
                           breaks = seq(0,15, by=2)) + 
        scale_color_manual(values = c(colorFox, colorCSL, colorSkunk)) +
        theme_bw(base_size = 20)

noaxes <- theme(plot.margin = unit(c(0,0,0,0), "cm"),
                axis.title.x  = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y.left  = element_blank(),
                axis.text.y.left = element_blank(),
                axis.ticks.y.left = element_blank(),
                aspect.ratio = 1)


yonly <-  theme(plot.margin = unit(c(0,0,0,0), "cm"),
                axis.title.x  = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank())


xonly <- theme(plot.margin = unit(c(0,0,0,0), "cm"),
               axis.title.y =  element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(), aspect.ratio = 1)


both <- theme(plot.margin = unit(c(0,0,0,0), "cm"), aspect.ratio = 1)


# Y-labels only -----------------------------------------------------------

autbra <- uni +
          ylab("Bra") +
          geom_jitter(aes(LogAut,LogBra,color=Species), 
                      alpha=.5,show.legend = F,cex=3) +
          yonly


autdja <- uni +
          ylab("Dja") +
          geom_jitter(aes(LogAut,LogDja,color=Species),
                      alpha=.5,show.legend = F,cex=3) +
          yonly


autict <- uni +
          ylab("Ict")+
          geom_jitter(aes(LogAut,LogIct,color=Species),
                      alpha=.5,show.legend = F,cex=3)+
          yonly


# X- and Y-labels -------------------------------------------------------------

autpom <- uni +
          ylab("Pom") +
          xlab("Aut") +
          geom_jitter(aes(LogAut,LogPom,color=Species),
                      alpha=.5,show.legend = F,cex=3) +
          both


# X-labels only -----------------------------------------------------------

brapom <- uni +
          xlab("Bra") +
          geom_jitter(aes(LogBra,LogPom,color=Species),
                      alpha=.5,show.legend = F,cex=3) +
          xonly


djapom <- uni +
          xlab("Dja") +
          geom_jitter(aes(LogDja,LogPom,color=Species),
                      alpha=.5,show.legend = F,cex=3) +
          xonly


ictpom  <- uni +
            xlab("Ict") +
            geom_jitter(aes(LogIct,LogPom,color=Species),
                        alpha=.5,show.legend = F,cex=3) +
            xonly


# No axis labels ----------------------------------------------------------
bradja <- uni +
          geom_jitter(aes(LogBra,LogDja,color=Species),
                      alpha=.5,show.legend = F,cex=3)+
          noaxes


braict <- uni +
          geom_jitter(aes(LogBra,LogIct,color=Species),
                      alpha=.5,show.legend = F,cex=3) +
          noaxes


djaict <- uni +
          geom_jitter(aes(LogDja,LogIct,color=Species),
                      alpha=.5,show.legend = F,cex=3) +
          noaxes



# Blank plots
nothing = ggplot() + theme_void()

# Create legend
legen <- data.frame(Species=c("Sea Lion","Fox","Skunk"),
                    col=c(colorFox,colorCSL,colorSkunk),
                    xx=c(4,4,4),yy=c(10,8,6))

legen$Species <- factor(legen$Species,legen$Species,legen$Species,ordered=TRUE)
leg <- ggplot(legen) +
        theme_void(base_size = 20) +
        theme(legend.position=c(0.5,0.5),plot.margin = unit(c(0,0,0,0), "cm"),
              aspect.ratio =1) +
        geom_point(aes(x=xx,y=yy,color=Species,fill=Species),cex=3,alpha=.5) +
        geom_point(aes(x=xx,y=yy,color=Species),cex=3,pch=1) +
        geom_point(aes(x=xx,y=yy),color="white",fill="white",cex=4) +
        scale_color_manual(values = c(colorCSL,colorFox,colorSkunk)) 


# Create scatterplot ------------------------------------------------------

png("figures/figS2.png",width = 800,height = 800)

autbra+nothing+nothing+nothing+
  autdja+bradja+nothing+leg+
  autict+braict+djaict+nothing+
  autpom+brapom+djapom+ictpom+
  plot_layout(nrow=4,ncol=4)

dev.off()



# Boxplot -----------------------------------------------------------------

all2 <- all %>%
        reshape::melt(., id = 'Species') %>%
        mutate(Species = factor(Species,
                                levels = c("Sea Lion","Island Fox","Skunk"),
                                labels = c("Sea Lion","Fox","Skunk"),
                                ordered=TRUE),
               both = paste(variable, Species)) %>%
        mutate(both = factor(both, both, both, ordered=T))


boxplot <- ggplot(all2) +
            theme_bw(base_size=18) +
            theme(panel.grid = element_blank(),
                  plot.margin = unit(c(1,1,1,1), "cm"),
                  aspect.ratio =.5) +
            labs(y = expression(paste(log[2]," MAT titer", sep="")), 
                 x = "") +
            geom_boxplot(aes(both,value,color=Species),show.legend = F) +
            geom_rect(aes(ymin=0, ymax=15, 
                          xmin="LogIct Sea Lion", xmax="LogBra Sea Lion"), 
                      alpha =0.5, fill="grey90",position = position_nudge(x=-.5)) +
            geom_rect(aes(ymin=0, ymax=15, 
                          xmin="LogAut Sea Lion", xmax="LogDja Sea Lion"), 
                      alpha =0.5, fill="grey90", position = position_nudge(x=-.5)) +
            geom_boxplot(aes(both,value,color=Species), show.legend = F) +
            scale_y_continuous(limits = c(-0.5, 15), 
                               breaks = seq(0,15, by=2),
                               labels=seq(0,15, by=2)) +
            scale_x_discrete(breaks = c("LogPom Fox","LogBra Fox","LogDja Fox",
                                        "LogIct Fox","LogAut Fox"),
                             labels=c("Pom","Bra","Dja","Ict","Aut")) +
            scale_color_manual(values = c(colorCSL,colorFox,colorSkunk)) 


png("figures/figS2-boxplot.png", width=500, height = 300)
  boxplot
dev.off()


# End script
