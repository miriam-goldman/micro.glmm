library("optparse")
suppressPackageStartupMessages(library(tidyverse))
source("helper_functions.R")
library(logr)
library(data.table)
source("glmm_functions.R")
option_list = list(
  make_option(c("-s", "--species_id"), type="character", default="s_id", 
              help="species id for labeling output file names (optional)", metavar="character"),
  make_option(c("-v", "--verbose"), type="logical", default=TRUE, 
              help="verbose, TRUE or FALSE"),
  make_option(c("-o", "--out_folder"), type="character", default=NULL,
              help="output folder name", metavar="character"),
  make_option(c("--Rdata"), type="character", default=NULL,
              help="Rdata made by pop_structure_test_command_line.R" ),
  make_option(c("--copy_number"), type="character", default=NULL,
              help=".tsv file with gene_id, sample_name, and copynumber for CNV data"),
  make_option(c("--SPA"), type="logical", default=FALSE,
              help="whether to use SPA for pvalues, will slow down output (optional)"),
  make_option(c("--scale_copynumber"), type="logical", default=FALSE,
              help="whether to scale the copynumber data for each gene (optional)"),
  make_option(c("--log_copynumber"), type="logical", default=FALSE,
              help="whether to log the copynumber data for each gene (optional)")
  
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
write.csv2(marker_test_df,file.path(output_dir, paste0(s_id,".marker_test.tsv")),sep="\t")

log_close()