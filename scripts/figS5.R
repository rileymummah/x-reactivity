## ---------------------------
## This code was written by: r.o. mummah
## For questions: rom5173@ucla.edu
## Date Created: 2024-04-08
## ---------------------------

## ---------------------------
## Objective: 
##   Create a boxplot figure of MAT titers stratefied by Lab and Serovar
## 
## Input:
##   data/LabComparisons.csv
##
## Output: 
##   figures/figS5.png
##
## ---------------------------

## load packages ---------------------------
library(tidyverse)
library(magrittr)


## load data ---------------------------

data <- read.csv('data/LabComparisons.csv')


## Figure S5 ---------------------------

ggplot(data=data, aes(Lab, Titer, color = Serovar)) + 
  geom_boxplot() + 
  scale_color_manual(values = c("cornflowerblue","black"))  + 
  scale_y_continuous(limits=c(0,10), breaks=seq(0,10)) +
  labs(y = expression(paste(log[2]," MAT titer")),
       x = NULL) +
  theme_classic(base_size = 14)

ggsave('figures/figS5.png', 
       width = 6.5, height = 4, units = 'in', dpi=600)

# End script