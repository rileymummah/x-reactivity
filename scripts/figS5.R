## ---------------------------
## This code was written by: r.o. mummah
## For questions: rom5173@ucla.edu
## ---------------------------

## ---------------------------
## Objective: Create Fig S4
## Longitudinal antibody titer dynamics in Channel Island foxes
## 
## Input:
##   SheddersTiterDecline.csv
##
## Output: 
##   FigS5
##
## ---------------------------

## load packages ---------------------------
library(tidyverse)
library(magrittr)
library(readr)


## load data ---------------------------

## Confirmed Shedders Titer Decline Plots
DF <- read.csv("data/SheddersTiterDecline.csv")

## ---------------------------

new_Labels <- c("02107"="1","07478"="2","14A0D"="3","23727"="4","23850"="5","24269"="6","37531"="7","37539"="8","39031"="9","39037"="10","59843"="11","61904"="12","78624"="13","86076"="14","90599"="15","90697"="16","98205"="17","B0F3B"="18","C0D60"="19","F4A38"="20")

ggplot(DF[!is.na(DF$PomonaLevel),],aes(DaysSinceKnownPos)) + 
  geom_line(aes(y = PomonaLevel, col='Pomona'), lty = "solid") + 
  geom_point(aes(y = PomonaLevel, col='Pomona')) + 
  geom_line(aes(y = AutumnalisLevel, col='Autumnalis'), lty = "dashed") + 
  geom_point(aes(y = AutumnalisLevel, col='Autumnalis'), pch=18, cex=2.5) + 
  geom_text(aes(label = CurrentShedPN, x = DaysSinceKnownPos, y = -2)) + 
  geom_text(aes(label = "PCR", x = -500, y = -2))  + 
  geom_text(aes(label = "MAT", x = -500, y = 10)) + 
  geom_hline(aes(yintercept = -0.8), linewidth = 0.5) + 
  facet_wrap(vars(Pittag), labeller = labeller(Pittag = new_Labels)) +
  scale_y_continuous(limits = c(-3,12), breaks = c(0,2,4,6,8,10))+ 
  scale_color_manual(values = c("cornflowerblue","black"),
                     labels = c("Autumnalis Titer","Pomona Titer"), name = "") +
  labs(y = expression(paste("PCR status | ",log[2]," MAT titer")), 
       x = "Days Since Known Positive") +
  theme_bw() +
  theme(strip.text = element_text(face='bold'),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        axis.title.x = element_text(vjust = -0.5))

ggsave('figures/figS5.png', 
       width = 9, height = 8, units = 'in', dpi=600)

# End script
