library(viridis)
labcolors <- viridis(6, option="D")

#### Figure 4? A
# labdata <- read.csv('data_raw/CDCCSUCornellCompare_reformatted.csv') %>%
#   mutate(Endpoint = ifelse(Endpoint == 16, 'Y', 'N'),
#          Titer = as.numeric(Titer)) %>%
#   select(-Date, -Pittag) %>%
#   pivot_wider(id_cols = c(RowID, Lab, Endpoint),
#               names_from = 'Serovar', values_from = 'Titer')


labdata <- read.csv('data_raw/CDCCSUCornellCompare.csv')


##3 lab compare sep graphs pom vs aut
plt6 <- ggplot(labdata,
               aes(x=CDCAutumnalis,y=CDCPomona,shape=factor(CDCAutPCH))) + 
  geom_jitter(width = 0.2, height = 0.2, color = labcolors[4], size=2) + 
  scale_x_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14)) + 
  scale_y_continuous(limits=c(0,14), breaks=c(0,2,4,6,8,10,12,14)) + 
  scale_shape_manual(values=c(0,16)) +
  geom_abline(intercept=0, slope=1,size =.8) + 
  labs(x="Aut", y="Pom",title="Lab A") +
  theme_minimal(base_size = 14) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

plt7 <- ggplot(CDCCSUCornellCompare,
               aes(x=CornellAutumnalis,y=CornellPomona)) + 
  geom_jitter(width = 0.2, height = 0.2, color = labcolors[4], size=2) + 
  scale_x_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + 
  scale_y_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + 
  scale_shape_manual(values=c(0,16)) +
  geom_abline(intercept=0, slope=1,size =.8) + 
  labs(x="Aut", y="Pom",title="Lab B") +
  theme_minimal(base_size = 14) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

plt8 <- ggplot(CDCCSUCornellCompare,
               aes(x=CSUAutumnalis,y=CSUPomona,shape=factor(CSUAutPCH))) + 
  geom_jitter(width = 0.2, height = 0.2, color = labcolors[4], size=2) + 
  scale_x_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + 
  scale_y_continuous(limits=c(0,14),breaks=c(0,2,4,6,8,10,12,14)) + 
  scale_shape_manual(values=c(0,16)) +
  geom_abline(intercept=0, slope=1,size =.8) + 
  labs(x="Aut", y="Pom",title="Lab C") +
  theme_minimal(base_size = 14) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")


library(gridExtra)


jpeg("figures/LabABCcomparePA_3plots.jpeg",res=300,height = 3, width=9, units = "in")

grid.arrange(plt6,plt7,plt8,nrow=1,ncol=3)

dev.off()