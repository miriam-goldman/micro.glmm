suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(tidyverse))
source("../R/helper_functions.R")
suppressPackageStartupMessages(library(logr))
suppressPackageStartupMessages(library(data.table))
source("../R/glmm_functions.R")
option_list = list(
  make_option(c("-s", "--species_id"), type="character", default="s_id", 
              help="species id for labeling output file names (optional)", metavar="character"),
  make_option(c("-v", "--verbose"), type="logical", default=TRUE, 
              help="verbose, TRUE or FALSE", metavar="logical"),
  make_option(c("-o", "--out_folder"), type="character", default=NULL,
              help="output folder name", metavar="character"),
  make_option(c("--Rdata"), type="character", default=NULL,
              help="file location of Rdata made by pop_structure_test_command_line.R",metavar="character"),
  make_option(c("--copy_number"), type="character", default=NULL,
              help=".tsv file with gene_id, sample_name, and copynumber for CNV data",metavar="character"),
  make_option(c("--SPA"), type="logical", default=FALSE,
              help="whether to use SPA for pvalues, will slow down output (optional)", metavar="logical"),
  make_option(c("--scale_copynumber"), type="logical", default=FALSE,
              help="whether to scale the copynumber data for each gene (optional)", metavar="logical"),
  make_option(c("--log_copynumber"), type="logical", default=FALSE,
              help="whether to log the copynumber data for each gene (optional)", metavar="logical"),
  make_option(c("--q_value"), type="numeric", default=.05,
              help="q value for marker test plot", metavar="numeric"),
  make_option(c("--alpha_value"), type="numeric", default=.05,
              help="alpha value for marker test plot", metavar="numeric")
  
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
verbose<-opt$verbose
output_dir<-test_dir(opt$out_folder,verbose)
s_id<-opt$species_id
tmp <- file.path(output_dir, paste0(s_id,".marker_test.log"))

lf <- log_open(tmp,show_notes=FALSE)
put("running validate genes command is",console = verbose)
command_arg_list<-commandArgs(trailingOnly = TRUE)
commad_message<-paste(unlist(command_arg_list),collapse = " ")

put(paste("Rscript marker_test_command_line.R",commad_message),console = verbose)

validate_marker_test(opt)



output_dir<-test_dir(opt$out_folder,verbose)
marker_test_df<-micro_glmm(glmm_fit,glm_fit0,GRM,copy_number_df,SPA=spa_opt,scale_g=scale_copynumber_opt,log_g=log_copynumber_opt)
pdf(file= file.path(output_dir,paste0(s_id,".marker_test_output.pdf")) ) 

# create a 2X2 grid 
par( mfrow= c(4,1) )
m<-nrow(marker_test_df)
bonferroni_cutoff<-opt$alpha_value/m
marker_test_df<-marker_test_df %>% arrange(pvalue)
q=opt$q_value
marker_test_df<-marker_test_df%>% mutate(rank=seq(1,m)) %>% mutate(bh_pvalue=(rank/m)*q)
write.table(marker_test_df,file.path(output_dir, paste0(s_id,".marker_test.tsv")),sep="\t",row.names=FALSE)
filtered_test_df<-marker_test_df %>% filter(pvalue<=bh_pvalue)
bh_cutoff<-max(filtered_test_df$pvalue)
print(marker_test_df[1:50,])

marker_test_df %>% ggplot(aes(x=pvalue)) +geom_histogram()+ggtitle(paste("Pvalue histogram for species ",s_id))
marker_test_df %>% ggplot(aes(y=-log10(pvalue),x=beta))+geom_hline(color="green",yintercept =-log10(bh_cutoff))+geom_hline(color="red",yintercept =-log10(bonferroni_cutoff))+geom_point()+ggtitle(paste("volcano plot for species ",s_id,"tau value is",glmm_fit$tau[2]))





log_close()