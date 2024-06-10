foxes <- read.csv("data/Foxes_shedding_or_pfge_without_duplicates.csv")
skunks <- read.csv("data/Skunks_pfge_or_shedding.csv")
csl <- read.csv("data/CSL_shedding_or_pfge.csv")

dplyr::bind_cols(csl$ID, foxes$ID, skunks$ID)
data.frame(csl = csl$ID,
           fox = c(foxes$ID, rep(NA,78)),
           skunk = c(skunks$ID, rep(NA,133))) %>%
  write.table(file = 'sppID.txt')
