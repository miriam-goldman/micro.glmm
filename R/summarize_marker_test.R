## combine marker files and output a summary

library("optparse")
library("tidyverse")
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="all marker files", metavar="character")
)
file_list<-option_list$i
all_marker_files<-lapply(file_list,read.csv)
all_marker_files %>% group_by(species_id,tau) %>% summarize(n()) %>% write.csv("summary_file.csv")