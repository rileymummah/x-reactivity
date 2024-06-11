## ---------------------------
## This code was written by: acr gomez & ro mummah
## For questions: rom5173@ucla.edu
## ---------------------------

## ---------------------------
## Objective: 
## Analysis of host x serovar patterns of absolute and relative MAT antibody 
## titers detected againstfive Leptospira serovars (Pomona, Djasiman, Autumnalis, 
## Bratislava, and Icterohaemorrhagiae) when the infecting serovar is 
## L. interrogans serovar Pomona
## 
## Input:
##   Foxes_shedding_or_pfge_without_duplicates.csv
##   Skunks_pfge_or_shedding.csv
##   CSL_shedding_or_pfge.csv
##
## Output: 
##   Statistical tests for supplemental tables
##
## ---------------------------

## load packages ---------------------------
library(tidyverse)
library(fmsb)
library(viridis)
library(lubridate)
library(npmv)


## main text analysis ---------------------------
options(stringsAsFactors = FALSE)
foxes <- read.csv("data/Foxes_shedding_or_pfge_without_duplicates.csv") %>%
          transmute_at(2:6, as.numeric) %>%
          mutate(spp = 'fox')

skunks <- read.csv("data/Skunks_pfge_or_shedding.csv") %>%
          transmute_at(2:6, as.numeric) %>%
          mutate(spp = 'skunk')

csl <- read.csv("data/CSL_shedding_or_pfge.csv") %>%
        transmute_at(2:6, as.numeric) %>%
        mutate(spp = 'csl')



### CSL ---------------------------------------------------------------------

# Do titer values differ within csl? YES
csl %>%
  select(-spp) %>%
  pivot_longer(cols = 1:5, 
               names_to = 'serovar',
               values_to = 'titer') -> csl2
kruskal.test(titer ~ serovar, data = csl2)

# Kruskal-Wallis rank sum test
# data:  titer by serovar
# Kruskal-Wallis chi-squared = 152.23, df = 4, p-value < 2.2e-16

# Which serovars are different?
pairwise.wilcox.test(csl2$titer, csl2$serovar,
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  csl2$titer and csl2$serovar 
# 
#        LogAut  LogBra  LogDja  LogIct 
# LogBra 0.17143 -       -       -      
# LogDja 0.57523 0.39500 -       -      
# LogIct 3.6e-05 0.01569 0.00095 -      
# LogPom 4.6e-16 < 2e-16 < 2e-16 < 2e-16
# 
# P value adjustment method: BH 



### Fox ---------------------------------------------------------------------

# Do titer values differ within foxes? YES
foxes %>%
  select(-spp) %>%
  pivot_longer(cols = 1:5, 
               names_to = 'serovar',
               values_to = 'titer') -> foxes2
kruskal.test(titer ~ serovar, data = foxes2)

# Kruskal-Wallis rank sum test
# data:  titer by serovar
# Kruskal-Wallis chi-squared = 166.41, df = 4, p-value < 2.2e-16

# Which serovars are different?
pairwise.wilcox.test(foxes2$titer, foxes2$serovar,
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# data:  foxes2$titer and foxes2$serovar 
# 
#        LogAut  LogBra  LogDja  LogIct 
# LogBra 5.6e-15 -       -       -      
# LogDja 0.0120  2.5e-11 -       -      
# LogIct < 2e-16 1.8e-05 < 2e-16 -      
# LogPom 0.0017  2.0e-12 0.5711  < 2e-16
# 
# P value adjustment method: BH 



### Skunk -------------------------------------------------------------------

# Do titer values differ within skunks? YES
skunks %>%
  select(-spp) %>%
  pivot_longer(cols = 1:5, 
               names_to = 'serovar',
               values_to = 'titer') -> skunks2
kruskal.test(titer ~ serovar, data = skunks2)

# Kruskal-Wallis rank sum test
# data:  titer by serovar
# Kruskal-Wallis chi-squared = 11.728, df = 4, p-value = 0.01949

# Which serovars are different?
pairwise.wilcox.test(skunks2$titer, skunks2$serovar,
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  skunks2$titer and skunks2$serovar 
# 
#        LogAut LogBra LogDja LogIct
# LogBra 0.208  -      -      -     
# LogDja 0.330  0.601  -      -     
# LogIct 0.208  1.000  0.601  -     
# LogPom 0.301  0.095  0.095  0.095 
# 
# P value adjustment method: BH 



### Are the titer profiles different across species? ------------------------
# ANOSIM (ANalysis Of Similarities) is a non-parametric test of significant difference between two or more groups, based on any distance measure (Clarke 1993)
data <- rbind(foxes, skunks, csl)

# ANOSIM test for dissimilarity
mat <- as.matrix(data[,1:5])

mod <- vegan::anosim(mat, data$spp, distance = "bray", permutations = 999)

# Strong statistical evidence that the titer profiles are different across host species

# vegan::anosim(x = mat, grouping = data$spp, permutations = 999, distance = "bray") 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: 0.2854 
# Significance: 0.001 
# 
# Permutation: free
# Number of permutations: 999

veganEx::anosim.pairwise(mat, data$spp)

#        pairs   anosimR p.value p.adj
# fox.vs.skunk 0.7855394   0.001 0.003
#   fox.vs.csl 0.2573589   0.001 0.003
# skunk.vs.csl 0.5804374   0.002 0.006

### Does Pomona differ across host species? ----
data %>%
  pivot_longer(cols = 1:5,
               names_to = 'serovar',
               values_to = 'titer') -> data2

csl.pom <- data2 %>% filter(spp == 'csl', serovar == 'LogPom')
fox.pom <- data2 %>% filter(spp == 'fox', serovar == 'LogPom')
skunk.pom <- data2 %>% filter(spp == 'skunk', serovar == 'LogPom')

kruskal.test(list(csl.pom$titer, fox.pom$titer, skunk.pom$titer))

# Kruskal-Wallis rank sum test
# data:  list(csl.pom$titer, fox.pom$titer, skunk.pom$titer)
# Kruskal-Wallis chi-squared = 70.444, df = 2, p-value = 5.049e-16

pairwise.wilcox.test(c(csl.pom$titer, fox.pom$titer, skunk.pom$titer),
                     c(csl.pom$spp, fox.pom$spp, skunk.pom$spp),
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# data:  c(csl.pom$titer, fox.pom$titer, skunk.pom$titer) and c(csl.pom$spp, fox.pom$spp, skunk.pom$spp) 
# 
#       csl     fox   
# fox   1.3e-15 -     
# skunk 0.0079  0.1662
# 
# P value adjustment method: BH 

### Does Djasiman differ across host species? ----

csl.dja <- data2 %>% filter(spp == 'csl', serovar == 'LogDja')
fox.dja <- data2 %>% filter(spp == 'fox', serovar == 'LogDja')
skunk.dja <- data2 %>% filter(spp == 'skunk', serovar == 'LogDja')

kruskal.test(list(csl.dja$titer, fox.dja$titer, skunk.dja$titer))

# Kruskal-Wallis rank sum test
# data:  list(csl.dja$titer, fox.dja$titer, skunk.dja$titer)
# Kruskal-Wallis chi-squared = 8.0699, df = 2, p-value = 0.01769

pairwise.wilcox.test(c(csl.dja$titer, fox.dja$titer, skunk.dja$titer),
                     c(csl.dja$spp, fox.dja$spp, skunk.dja$spp),
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# data:  c(csl.dja$titer, fox.dja$titer, skunk.dja$titer) and c(csl.dja$spp, fox.dja$spp, skunk.dja$spp) 
# 
#       csl    fox   
# fox   0.2270 -     
# skunk 0.0305 0.0049
# 
# P value adjustment method: BH 


### Does Autumnalis differ across host species? ----

csl.aut <- data2 %>% filter(spp == 'csl', serovar == 'LogAut')
fox.aut <- data2 %>% filter(spp == 'fox', serovar == 'LogAut')
skunk.aut <- data2 %>% filter(spp == 'skunk', serovar == 'LogAut')

kruskal.test(list(csl.aut$titer, fox.aut$titer, skunk.aut$titer))

# Kruskal-Wallis rank sum test
# data:  list(csl.aut$titer, fox.aut$titer, skunk.aut$titer)
# Kruskal-Wallis chi-squared = 5.7163, df = 2, p-value = 0.05737

pairwise.wilcox.test(c(csl.aut$titer, fox.aut$titer, skunk.aut$titer),
                     c(csl.aut$spp, fox.aut$spp, skunk.aut$spp),
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# data:  c(csl.aut$titer, fox.aut$titer, skunk.aut$titer) and c(csl.aut$spp, fox.aut$spp, skunk.aut$spp) 
# 
#       csl    fox   
# fox   0.9911 -     
# skunk 0.0617 0.0064
# 
# P value adjustment method: BH

### Does Ictero differ across host species? ----

csl.ict <- data2 %>% filter(spp == 'csl', serovar == 'LogIct')
fox.ict <- data2 %>% filter(spp == 'fox', serovar == 'LogIct')
skunk.ict <- data2 %>% filter(spp == 'skunk', serovar == 'LogIct')

kruskal.test(list(csl.ict$titer, fox.ict$titer, skunk.ict$titer))

# Kruskal-Wallis rank sum test
# data:  list(csl.ict$titer, fox.ict$titer, skunk.ict$titer)
# Kruskal-Wallis chi-squared = 81.767, df = 2, p-value < 2.2e-16

pairwise.wilcox.test(c(csl.ict$titer, fox.ict$titer, skunk.ict$titer),
                     c(csl.ict$spp, fox.ict$spp, skunk.ict$spp),
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  c(csl.ict$titer, fox.ict$titer, skunk.ict$titer) and c(csl.ict$spp, fox.ict$spp, skunk.ict$spp) 
# 
#       csl    fox   
# fox   <2e-16 -     
# skunk 0.0027 0.1245
# 
# P value adjustment method: BH 


### Does Bratislava differ across host species? ----

csl.bra <- data2 %>% filter(spp == 'csl', serovar == 'LogBra')
fox.bra <- data2 %>% filter(spp == 'fox', serovar == 'LogBra')
skunk.bra <- data2 %>% filter(spp == 'skunk', serovar == 'LogBra')

kruskal.test(list(csl.bra$titer, fox.bra$titer, skunk.bra$titer))

# Kruskal-Wallis rank sum test
# data:  list(csl.bra$titer, fox.bra$titer, skunk.bra$titer)
# Kruskal-Wallis chi-squared = 36.867, df = 2, p-value = 9.872e-09

pairwise.wilcox.test(c(csl.bra$titer, fox.bra$titer, skunk.bra$titer),
                     c(csl.bra$spp, fox.bra$spp, skunk.bra$spp),
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# data:  c(csl.bra$titer, fox.bra$titer, skunk.bra$titer) and c(csl.bra$spp, fox.bra$spp, skunk.bra$spp) 
# 
#       csl     fox  
# fox   6.1e-08 -    
# skunk 0.011   0.011
# 
# P value adjustment method: BH 





# supplemental text analysis (PFGE only) ----------------------------------

# CSL
csl.pfge <- read.table("data/CSL_pfge_confirmed_only.csv", 
                       header=T, sep=',') %>% # PFGE
            mutate(Species = 'CSL')

# Foxes
fox.pfge <- read.table("data/Foxes_pfge_only_noduplicates.csv", 
                       header=T, sep=',') %>% # PFGE only
            mutate(Species = 'fox')


# Skunks
# Only 1 PFGE skunk, so we use all of them for both datasets
skunks <- read.csv("data/Skunks_pfge_or_shedding.csv", 
                   header=T, sep=',') %>%
          mutate(Species = 'skunk')


### CSL ---------------------------------------------------------------------

# Do titer values differ within csl? YES
csl.pfge %>%
  select(-Species, -ID) %>%
  pivot_longer(cols = 1:5, 
               names_to = 'serovar',
               values_to = 'titer') -> csl2
kruskal.test(titer ~ serovar, data = csl2)

# Kruskal-Wallis rank sum test
# data:  titer by serovar
# Kruskal-Wallis chi-squared = 23.933, df = 4, p-value = 8.239e-05

# Which serovars are different?
pairwise.wilcox.test(csl2$titer, csl2$serovar,
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  csl2$titer and csl2$serovar 
# 
#        LogAut  LogBra  LogDja  LogIct 
# LogBra 0.85319 -       -       -      
# LogDja 1.00000 0.85319 -       -      
# LogIct 0.22994 0.22994 0.22994 -      
# LogPom 0.00399 0.00169 0.00304 0.00016
# 
# P value adjustment method: BH 



### Fox ---------------------------------------------------------------------

# Do titer values differ within foxes? YES
fox.pfge %>%
  select(-ID, -Species) %>%
  pivot_longer(cols = 1:5, 
               names_to = 'serovar',
               values_to = 'titer') -> foxes2
kruskal.test(titer ~ serovar, data = foxes2)

# Kruskal-Wallis rank sum test
# data:  titer by serovar
# Kruskal-Wallis chi-squared = 30.144, df = 4, p-value = 4.575e-06

# Which serovars are different?
pairwise.wilcox.test(foxes2$titer, foxes2$serovar,
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction
# 
# data:  foxes2$titer and foxes2$serovar
# 
#        LogAut  LogBra  LogDja  LogIct 
# LogBra 0.00223 -       -       -      
# LogDja 0.40709 0.03007 -       -      
# LogIct 0.00047 0.30984 0.00223 -      
# LogPom 0.45069 0.00223 0.61592 0.00047
# 
# P value adjustment method: BH



### Skunk -------------------------------------------------------------------
## Analysis is the same as the main text because only 2 skunks were PFGE confirmed
# Do titer values differ within skunks? YES
skunks %>%
  select(-spp) %>%
  pivot_longer(cols = 1:5, 
               names_to = 'serovar',
               values_to = 'titer') -> skunks2
kruskal.test(titer ~ serovar, data = skunks2)

# Kruskal-Wallis rank sum test
# data:  titer by serovar
# Kruskal-Wallis chi-squared = 11.728, df = 4, p-value = 0.01949

# Which serovars are different?
pairwise.wilcox.test(skunks2$titer, skunks2$serovar,
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  skunks2$titer and skunks2$serovar 
# 
#        LogAut LogBra LogDja LogIct
# LogBra 0.208  -      -      -     
# LogDja 0.330  0.601  -      -     
# LogIct 0.208  1.000  0.601  -     
# LogPom 0.301  0.095  0.095  0.095 
# 
# P value adjustment method: BH 



### Are the titer profiles different across species? ------------------------

data <- rbind(fox.pfge, skunks, csl.pfge)

# npmv::nonpartest(LogPom|LogIct|LogAut|LogBra|LogDja ~ spp, data = data, permtest=TRUE, 
#            permreps=10000, plots=T)

# ANOSIM test for dissimilarity
mat <- as.matrix(data[,2:6])

mod <- vegan::anosim(mat, data$Species, distance = "bray", permutations = 999)

# Strong statistical evidence that the titer profiles are different across host species

# vegan::anosim(x = mat, grouping = data$spp, permutations = 999, distance = "bray") 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: 0.4645 
# Significance: 0.001 
# 
# Permutation: free
# Number of permutations: 999


### Does Pomona differ across host species? ----
data %>%
  select(-ID) %>%
  pivot_longer(cols = 1:5,
               names_to = 'serovar',
               values_to = 'titer') %>%
  rename('spp' = 'Species') -> data2

csl.pom <- data2 %>% filter(spp == 'CSL', serovar == 'LogPom')
fox.pom <- data2 %>% filter(spp == 'fox', serovar == 'LogPom')
skunk.pom <- data2 %>% filter(spp == 'skunk', serovar == 'LogPom')

kruskal.test(list(csl.pom$titer, fox.pom$titer, skunk.pom$titer))

# Kruskal-Wallis rank sum test
# data:  list(csl.pom$titer, fox.pom$titer, skunk.pom$titer)
# Kruskal-Wallis chi-squared = 17.4, df = 2,
# p-value = 0.0001666

pairwise.wilcox.test(c(csl.pom$titer, fox.pom$titer, skunk.pom$titer),
                     c(csl.pom$spp, fox.pom$spp, skunk.pom$spp),
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# data:  c(csl.pom$titer, fox.pom$titer, skunk.pom$titer) and c(csl.pom$spp, fox.pom$spp, skunk.pom$spp) 
# 
#       CSL    fox   
# fox   0.0014 -     
# skunk 0.0070 0.1066
# 
# P value adjustment method: BH 

### Does Djasiman differ across host species? ----

csl.dja <- data2 %>% filter(spp == 'CSL', serovar == 'LogDja')
fox.dja <- data2 %>% filter(spp == 'fox', serovar == 'LogDja')
skunk.dja <- data2 %>% filter(spp == 'skunk', serovar == 'LogDja')

kruskal.test(list(csl.dja$titer, fox.dja$titer, skunk.dja$titer))

# Kruskal-Wallis rank sum test
# data:  list(csl.dja$titer, fox.dja$titer, skunk.dja$titer)
# Kruskal-Wallis chi-squared = 9.3394, df = 2, p-value = 0.009375

pairwise.wilcox.test(c(csl.dja$titer, fox.dja$titer, skunk.dja$titer),
                     c(csl.dja$spp, fox.dja$spp, skunk.dja$spp),
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# data:  c(csl.dja$titer, fox.dja$titer, skunk.dja$titer) and c(csl.dja$spp, fox.dja$spp, skunk.dja$spp) 
# 
#       CSL   fox  
# fox   0.171 -    
# skunk 0.017 0.017
# 
# P value adjustment method: BH 


### Does Autumnalis differ across host species? ----

csl.aut <- data2 %>% filter(spp == 'CSL', serovar == 'LogAut')
fox.aut <- data2 %>% filter(spp == 'fox', serovar == 'LogAut')
skunk.aut <- data2 %>% filter(spp == 'skunk', serovar == 'LogAut')

kruskal.test(list(csl.aut$titer, fox.aut$titer, skunk.aut$titer))

# Kruskal-Wallis rank sum test
# data:  list(csl.aut$titer, fox.aut$titer, skunk.aut$titer)
# Kruskal-Wallis chi-squared = 6.0394, df = 2, p-value = 0.04882

pairwise.wilcox.test(c(csl.aut$titer, fox.aut$titer, skunk.aut$titer),
                     c(csl.aut$spp, fox.aut$spp, skunk.aut$spp),
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# data:  c(csl.aut$titer, fox.aut$titer, skunk.aut$titer) and c(csl.aut$spp, fox.aut$spp, skunk.aut$spp) 
# 
#       CSL   fox  
# fox   0.528 -    
# skunk 0.082 0.016
# 
# P value adjustment method: BH 

### Does Ictero differ across host species? ----

csl.ict <- data2 %>% filter(spp == 'CSL', serovar == 'LogIct')
fox.ict <- data2 %>% filter(spp == 'fox', serovar == 'LogIct')
skunk.ict <- data2 %>% filter(spp == 'skunk', serovar == 'LogIct')

kruskal.test(list(csl.ict$titer, fox.ict$titer, skunk.ict$titer))

# Kruskal-Wallis rank sum test
# data:  list(csl.ict$titer, fox.ict$titer, skunk.ict$titer)
# Kruskal-Wallis chi-squared = 21.308, df = 2, p-value = 2.36e-05

pairwise.wilcox.test(c(csl.ict$titer, fox.ict$titer, skunk.ict$titer),
                     c(csl.ict$spp, fox.ict$spp, skunk.ict$spp),
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# data:  c(csl.ict$titer, fox.ict$titer, skunk.ict$titer) and c(csl.ict$spp, fox.ict$spp, skunk.ict$spp) 
# 
#       CSL     fox    
# fox   0.00023 -      
# skunk 0.00367 0.08191
# 
# P value adjustment method: BH 


### Does Bratislava differ across host species? ----

csl.bra <- data2 %>% filter(spp == 'CSL', serovar == 'LogBra')
fox.bra <- data2 %>% filter(spp == 'fox', serovar == 'LogBra')
skunk.bra <- data2 %>% filter(spp == 'skunk', serovar == 'LogBra')

kruskal.test(list(csl.bra$titer, fox.bra$titer, skunk.bra$titer))

# Kruskal-Wallis rank sum test
# data:  list(csl.bra$titer, fox.bra$titer, skunk.bra$titer)
# Kruskal-Wallis chi-squared = 14.016, df = 2, p-value = 0.0009044

pairwise.wilcox.test(c(csl.bra$titer, fox.bra$titer, skunk.bra$titer),
                     c(csl.bra$spp, fox.bra$spp, skunk.bra$spp),
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# data:  c(csl.bra$titer, fox.bra$titer, skunk.bra$titer) and c(csl.bra$spp, fox.bra$spp, skunk.bra$spp) 
# 
#       CSL    fox   
# fox   0.0091 -     
# skunk 0.0148 0.0336
# 
# P value adjustment method: BH 


# Are MAT+ skunks comparable to PFGE/culture-confirmed skunks? ------------
skunk.mat <- read.csv("data/Skunks_matpos_nodup.csv") %>%
              mutate(test = 'MAT')
skunk.pcr <- read.csv("data/Skunks_pfge_or_shedding.csv") %>%
              mutate(test = 'PFGE/culture')

data <- rbind(skunk.mat, skunk.pcr)

# ANOSIM test for dissimilarity
mat <- as.matrix(data[,2:6])

mod <- vegan::anosim(mat, data$test, distance = "bray", permutations = 999)

# Strong statistical evidence that the titer profiles are different across host species

data <- rbind(foxes, skunks, csl)

# ANOSIM test for dissimilarity
mat <- as.matrix(data[,1:5])

mod <- vegan::anosim(mat, data$spp, distance = "bray", permutations = 999)

# No statistical evidence that the titer profiles are different between MAT+ and PFGE/culture+ skunks.

# vegan::anosim(x = mat, grouping = data$test, permutations = 999,      distance = "bray") 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: 0.04494 
# Significance: 0.356 
# 
# Permutation: free
# Number of permutations: 999


# Do Pom and Aut titers differ across laboratories? -----------------------

data <- read.csv('data_raw/LabCompareBoth3.csv') %>%
  pivot_wider(id_cols = c('Pittag','Lab'), 
              names_from = 'Serovar', 
              values_from = 'Titer') %>%
  select(-Pittag) %>%
  pivot_longer(cols = 2:3, 
               names_to = 'serovar',
               values_to = 'titer')

# Pomona
data %>%
  filter(serovar == "Pomona") %>%
  kruskal.test(titer ~ Lab, data = .)

# Kruskal-Wallis rank sum test
# data:  titer by Lab
# Kruskal-Wallis chi-squared = 15.927, df = 2, p-value = 0.000348

# Which labs are different?
pairwise.wilcox.test(data$titer[data$serovar == "Pomona"], data$Lab[data$serovar == "Pomona"],
                     p.adjust.method = "BH")

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  data$titer[data$serovar == "Pomona"] and data$Lab[data$serovar == "Pomona"] 
# 
#       Lab A   Lab B  
# Lab B 0.02793 -      
# Lab C 0.00073 0.02356



# Autumnalis
data %>%
  filter(serovar == "Autumnalis") %>%
  kruskal.test(titer ~ Lab, data = .)

# Kruskal-Wallis rank sum test
# data:  titer by Lab
# Kruskal-Wallis chi-squared = 2.8709, df = 2, p-value = 0.238


# Calculate statistic for 4-fold increase

read.csv("data/Foxes_shedding_or_pfge_without_duplicates.csv") %>%
  select(ID, LogPom, LogAut) %>%
  mutate(`4fold` = ifelse(LogAut >= (LogPom + 2), 1, 0)) %>%
  summarize(n = sum(`4fold`),
            total = n(),
            percent = n/total*100)





# End script