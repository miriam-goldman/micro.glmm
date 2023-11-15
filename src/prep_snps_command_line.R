library("optparse",quietly=TRUE)
#source("../R/prep_snps.R")
library(devtools)
suppressPackageStartupMessages(library(micro.glmm))
suppressPackageStartupMessages(library(tidyverse,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE))
suppressPackageStartupMessages(library(logr))
suppressPackageStartupMessages(library(data.table))
library(RcppXPtrUtils,quietly=TRUE,warn.conflicts=FALSE)
option_list = list(
  make_option(c("-s", "--species_id"), type="character", default="s_id", 
              help="species id for labeling output file names (optional)", metavar="character"),
  make_option(c("-v", "--verbose"), type="logical", default=TRUE, 
              help="verbose, TRUE or FALSE", metavar="logical"),
  make_option(c("--make_plots"), type="logical", default=TRUE, 
              help="make plots, TRUE or FALSE", metavar="logical"),
  make_option(c("-o", "--out_folder"), type="character", default=NULL,
              help="output folder name", metavar="character"),
  make_option(c("-i", "--snps_info_file"), type="character", default=NULL,
              help="snps info file name and path", metavar="character"),
  make_option(c("-f", "--snps_freq_file"), type="character", default=NULL,
              help="snps freq file name and path", metavar="character"),
  make_option(c("-d", "--snps_depth_file"), type="character", default=NULL,
              help="snps depth file name and path", metavar="character"),
  make_option(c("-c", "--centroid_to_repgenes_file"), type="character", default=NULL,
              help="pangenome centroid to repgenes file (optional)", metavar="character"),
  make_option(c("-p", "--centroid_prevalence_file"), type="character", default=NULL,
              help="centriod prevalence file (optional)" , metavar="character"),
  make_option(c("--centroid_prevalence_cutoff"), type="numeric", default=.7,
              help="centriod prevalence file (optional)" ,metavar="numeric"),
  make_option(c("--genes_summary"), type="character", default=NULL,
              help="genes summary from MIDAS2 output used to filter sample list (optional)",metavar="character" ),
  make_option(c("-q", "--run_qp"), type="logical", default=FALSE, 
              help="option to run qp or not (optional)",metavar="logical"),
  make_option(c("-u", "--median_upper_filter"), type="integer", default=3,
              help="upper limit on median filter num of times median value of depth to let pass for each sample for QP (optional)",metavar="numeric"),
  make_option(c("-l", "--median_lower_filter"), type="numeric", default=.3,
              help="lower limit on median filter num of times median value of depth to let pass for each sample for QP (optional)",metavar="numeric"),
  make_option(c("-a", "--abosulte_filter"), type="numeric", default=3,
              help="lowest depth to all to pass filters",metavar="integer"),
  make_option(c("-m", "--sample_median_depth_filter"), type="integer", default=10,
              help="median depth needed for sample to be considered to have the species",metavar="integer"),
  make_option(c("-n","--number_of_samples_for_sites"), type="integer", default=5,
              help="number of samples needed for site to be considered to have the species",metavar="integer")
  
)

## filter SNPs make plots to explain SNPs and make distances

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
verbose<-opt$verbose
output_dir<-test_dir(opt$out_folder,verbose)
s_id<-opt$species_id
tmp <- file.path(output_dir, paste0(s_id,".prep_snps.log"))

lf <- log_open(tmp,show_notes=FALSE)
put("running validate genes command is",console = verbose)
command_arg_list<-commandArgs(trailingOnly = TRUE)
commad_message<-paste(unlist(command_arg_list),collapse = " ")

put(paste("Rscript prep_snps_command_line.R",commad_message),console = verbose)
validate_snps_output=validate_snps_input(opt)#parse_genes list

do.call(prep_snps_function_R,validate_snps_output)
# Close log
log_close()