## ---------------------------
## This code was written by: ro mummah & acr gomez
## For questions: rom5173@ucla.edu
## ---------------------------

## ---------------------------
## Objective: Create scatterplot in Fig2
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
##   Scatterplot for Fig2
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

## load data ---------------------------
options(stringsAsFactors = FALSE)
foxes <- read.csv("data/Foxes_PCR-PFGE.csv")
skunks <- read.csv("data/Skunks_PCR-PFGE.csv")
csl <- read.csv("data/CSL_PCR-PFGE.csv")


## CREATE SCATTERPLOT ---------------------------

fullnames <- c("Pomona","Icterohaemorrhagiae","Bratislava","Autumnalis","Djasiman")

csl <- transmute_at(csl,2:6,as.numeric)
csl$Species ="Sea Lion"
foxes <- transmute_at(foxes,2:6,as.numeric)
foxes$Species = "Island Fox"
skunks <- transmute_at(skunks,2:6,as.numeric)
skunks$Species = "Skunk"

all <- rbind(csl,foxes,skunks)


# Universal plotting ------------------------------------------------------

# Set colors
colorCSL <- viridis(9, option="D")[2]
colorFox <- viridis(9, option="D")[5]
colorSkunk <- viridis(9, option="D")[8]

# Universal baseplot
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

yonly <- theme(plot.margin = unit(c(0,0,0,0), "cm"),
               axis.title.x  = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               aspect.ratio = 1)

xonly <- theme(plot.margin = unit(c(0,0,0,0), "cm"),
               axis.title.y =  element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               aspect.ratio = 1)

both <- theme(plot.margin = unit(c(0,0,0,0), "cm"),
              aspect.ratio = 1)


# Y-labels only -----------------------------------------------------------

autbra <- uni +
          geom_jitter(aes(LogAut, LogBra, color=Species),
                      alpha=.5, show.legend = F, cex=3) + 
          ylab("Bra") +
          yonly


autdja <- uni +
          geom_jitter(aes(LogAut, LogDja, color=Species), 
                      alpha=.5, show.legend = F, cex=3) +
          ylab("Dja") +
          yonly

autict <- uni +
          geom_jitter(aes(LogAut,LogIct,color=Species), 
                      alpha=.5, show.legend = F, cex=3) +
          ylab("Ict") +
          yonly


# X- and Y-labels -------------------------------------------------------------

autpom <- uni +
          geom_jitter(aes(LogAut, LogPom, color=Species), 
                      alpha=.5, show.legend = F, cex=3) +
          xlab("Aut") + ylab("Pom") +
          both


# X-labels only -----------------------------------------------------------

brapom <- uni +
          geom_jitter(aes(LogBra, LogPom, color=Species), 
                      alpha=.5, show.legend = F, cex=3) +
          xlab("Bra") +
          xonly

djapom <- uni +
          geom_jitter(aes(LogDja, LogPom, color=Species), 
                      alpha=.5, show.legend = F, cex=3) +
          xlab("Dja") +
          xonly

ictpom <- uni +
          geom_jitter(aes(LogIct, LogPom, color=Species), 
                      alpha=.5, show.legend = F, cex=3) +
          xlab("Ict") + 
          xonly



# No axis labels ----------------------------------------------------------

bradja <- uni +
          geom_jitter(aes(LogBra, LogDja, color=Species), 
                      alpha=.5, show.legend = F, cex=3) +
          noaxes

braict <- uni +
          geom_jitter(aes(LogBra, LogIct, color=Species),
                      alpha=.5, show.legend = F, cex=3) +
          noaxes

djaict <- uni +
          geom_jitter(aes(LogDja, LogIct, color=Species),
                      alpha=.5, show.legend = F, cex=3) +
          noaxes



# Blank plots
nothing <- ggplot() + theme_void()

# Create legend
legen <- data.frame(Species=c("Sea Lion", "Fox", "Skunk"),
                    col=c(colorFox, colorCSL, colorSkunk),
                    xx=c(4, 4, 4),
                    yy=c(10, 8, 6))

leg <- ggplot(legen) +
        theme_void(base_size = 20) +
        theme(legend.position=c(0.5,0.5)) +
        geom_point(aes(x=xx, y=yy, color=Species, fill=Species),
                   cex=3, alpha=.5) +
        geom_point(aes(x=xx, y=yy, color=Species), cex=3,pch=1) +
        geom_point(aes(x=xx, y=yy), color="white",fill="white",cex=4) +
        scale_color_manual(values = c(colorFox, colorCSL, colorSkunk))


# Create scatterplot ------------------------------------------------------

autbra+nothing+nothing+nothing+
autdja+bradja+nothing+leg+
autict+braict+djaict+nothing+
autpom+brapom+djapom+ictpom +
plot_layout(nrow=4, ncol=4) -> scatter


## CREATE BOXPLOT ---------------------------



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


# Combine plots -----------------------------
library(cowplot)

tiff("figures/fig2.tiff", width = 800, height = 800)
ggdraw() +
  draw_plot(scatter) +
  draw_plot(boxplot, x = 0.53, y = .6, 0.45, 0.45)
dev.off()
  

# End script
