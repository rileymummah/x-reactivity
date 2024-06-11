## ---------------------------
## This code was written by: r.o. mummah
## For questions: rom5173@ucla.edu
## ---------------------------

## ---------------------------
## Objective: Comparison of antibody titer results for fox serum samples 
## evaluated at three testing laboratories
##
## Input:
##   data/LabComparisons.csv
##   data_raw/CDCCSUCornellCompare_reformatted.csv
##
## Output: 
##   Fig4
##
## ---------------------------

## load packages ---------------------------
library(tidyverse)
library(magrittr)
library(viridis)
library(patchwork)


## load data ---------------------------

data <- read.csv('data/labcomparisonEndpoint.csv')

labdata <- read.csv('data/labcomparisonNonEndpoint.csv')



## Figure 4A ---------------------------

base <- ggplot() +
        scale_x_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14)) +
        scale_y_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14)) +
        geom_abline(intercept=0, slope=1, size =.8) +
        labs(x="Aut", y="Pom") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5),
              aspect.ratio = 1,
              legend.position = "none")


labcolors <- viridis(6, option="D")


plt1 <- base +
        geom_jitter(data = labdata %>% filter(Lab == 'Lab A'),
                    aes(x=Autumnalis, y=Pomona, shape=factor(AutEndpoint)),
                    width = 0.2, height = 0.2, color = labcolors[4], size=2, alpha = 0.75) +
        scale_shape_manual(values=c(0,16)) +
        labs(title="Lab A")


plt2 <- base +
        geom_jitter(data = labdata %>% filter(Lab == 'Lab B'),
                    aes(x=Autumnalis, y=Pomona),
                    width = 0.2, height = 0.2, color = labcolors[4], size=2, alpha = 0.75) +
        labs(title="Lab B")


plt3 <- base +
        geom_jitter(data = labdata %>% filter(Lab == 'Lab C'),
                    aes(x=Autumnalis, y=Pomona, shape=factor(AutEndpoint)),
                    width = 0.2, height = 0.2, color = labcolors[4], size=2, alpha = 0.75) +
        scale_shape_manual(values=c(1,16)) +
        labs(title="Lab C")


# Assemble plot
top <- plt1 + plt2 + plt3 + plot_layout(nrow = 1)



## Figure 4B ---------------------------

ggplot(data=data, aes(Lab, Titer, color = Serovar)) + 
  geom_boxplot() + 
  scale_color_manual(values = c("cornflowerblue","black"))  + 
  scale_y_continuous(limits=c(0,10), breaks=seq(0,10)) +
  labs(y = expression(paste(log[2]," MAT titer")),
       x = NULL) +
  theme_classic() +
  theme(legend.position = 'bottom') -> bottom


# Assemble plot
top / bottom + plot_layout(nrow = 2) + plot_annotation(tag_levels = list(c('A','','','B')))


ggsave('figures/fig4.png', 
       width = 6.5, height = 5.5, units = 'in', dpi=600)

# End script