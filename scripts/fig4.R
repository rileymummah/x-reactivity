## ---------------------------
## This code was written by: r.o. mummah
## For questions: rmummah@usgs.gov
## Date Created: 2024-04-08
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
library(viridis)
library(patchwork)


## load data ---------------------------

data <- read.csv('data/LabComparisons.csv')

labdata <- read.csv('data_raw/CDCCSUCornellCompare.csv')

# labdata <- read.csv('data_raw/CDCCSUCornellCompare_reformatted.csv') %>%
#   mutate(Endpoint = ifelse(Endpoint == 16, 'Y', 'N'),
#          Titer = as.numeric(Titer)) %>%
#   select(-Date, -Pittag) %>%
#   pivot_wider(id_cols = c(RowID, Lab, Endpoint),
#               names_from = 'Serovar', values_from = 'Titer')



## Figure 4B ---------------------------

labcolors <- viridis(6, option="D")


plt6 <- ggplot(labdata,
               aes(x=CDCAutumnalis, y=CDCPomona, shape=factor(CDCAutPCH))) + 
  geom_jitter(width = 0.2, height = 0.2, color = labcolors[4], size=2, alpha = 0.75) + 
  scale_x_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14)) + 
  scale_y_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14)) + 
  scale_shape_manual(values=c(0,16)) +
  geom_abline(intercept=0, slope=1, size =.8) + 
  labs(x="Aut", y="Pom", title="Lab A") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5),
        aspect.ratio = 1, 
        legend.position = "none")

plt7 <- ggplot(CDCCSUCornellCompare,
               aes(x=CornellAutumnalis, y=CornellPomona)) + 
  geom_jitter(width = 0.2, height = 0.2, color = labcolors[4], size=2, alpha = 0.75) + 
  scale_x_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + 
  scale_y_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + 
  scale_shape_manual(values=c(0,16)) +
  geom_abline(intercept=0, slope=1, size =.8) + 
  labs(x="Aut", y="Pom", title="Lab B") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5),
        aspect.ratio = 1, 
        legend.position = "none")

plt8 <- ggplot(CDCCSUCornellCompare,
               aes(x=CSUAutumnalis, y=CSUPomona, shape=factor(CSUAutPCH))) + 
  geom_jitter(width = 0.2, height = 0.2, color = labcolors[4], size=2, alpha = 0.75) + 
  scale_x_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + 
  scale_y_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + 
  scale_shape_manual(values=c(0,16)) +
  geom_abline(intercept=0, slope=1, size =.8) + 
  labs(x="Aut", y="Pom", title="Lab C") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5),
        aspect.ratio = 1, 
        legend.position = "none")


# Assemble plot
top <- plt6 + plt7 + plt8 + plot_layout(nrow = 1) 


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