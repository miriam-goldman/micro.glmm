library("optparse")
library(logr)
option_list = list(
  make_option(c("-s", "--species_id"), type="character", default="s_id", 
              help="species id for labeling output file names (optional)", metavar="character"),
  make_option(c("-v", "--verbose"), type="logical", default=TRUE, 
              help="verbose, TRUE or FALSE", metavar="logical"),
  make_option(c("--make_plots"), type="logical", default=FALSE, 
              help="make_plots, TRUE or FALSE", metavar="logical"),
  make_option(c("--write_csv"), type="logical", default=FALSE, 
              help="write_csv, TRUE or FALSE", metavar="logical"),
  make_option(c("-o", "--out_folder"), type="character", default=NULL,
              help="output folder name", metavar="character"),
  make_option(c("-c", "--genes_copynumber_file"), type="character", default=NULL,
              help="snps freq file name and path", metavar="character"),
  make_option(c("-d", "--genes_depth_file"), type="character", default=NULL,
              help="snps depth file name and path", metavar="character"),
  make_option(c("-p", "--centroid_prevalence_file"), type="character", default=NULL,
              help="centriod prevalence file (optional)" , metavar="character"),
  make_option(c("--centroid_prevalence_cutoff"), type="numeric", default=.7,
              help="centriod prevalence file (optional)",metavar="numeric"),
  make_option(c("--genes_summary"), type="character", default=NULL,
              help="genes summary from MIDAS2 output used to filter sample list (optional)",metavar="character" ),
  make_option(c("--metadata"), type="character", default=NULL,
              help="metadata for samples; controls coded as 0s if using to filter (optional)",metavar="character" ),
  make_option(c("--min_num_control"), type="numeric", default=5,
              help="min number of samples that are controls to use gene (optional)",metavar="numeric" ),
  make_option(c("--min_num_case"), type="numeric", default=5,
              help="min number of samples that are case to use gene (optional)",metavar="numeric" ),
  make_option(c("--GRM"), type="character", default=NULL,
              help="GRM used to check sample list location of tsv(optional)",metavar="character" ),
  make_option(c("-n","--number_of_samples_for_copynumber"), type="integer", default=5,
              help="number of samples needed for gene to be considered",metavar="integer"),
  make_option(c("-m","--depth_cutoff"), type="integer", default=5,
              help="depth needed for gene to be included",metavar="integer"),
  make_option(c("--log_scale"),type="logical",default = TRUE,
              help="whether to log scale copynumber",metavar="logical"),
  make_option(c("--mean_center"),type="logical",default = TRUE,
              help="whether to mean scale copynumber",metavar="logical"),
  make_option(c("--var_filter"),type="numeric",default = 0,
              help="whether to filter copynumbers on variability (optional)",metavar="numeric")
)
suppressPackageStartupMessages(library(micro.glmm))
#' Rscript prep_genes.R 
#' --species_id 104158 
#' -v TRUE 
#' -o "~/Projects/micro-glmm/test_output" 
#' -c "/Users/miriamgoldman/Projects/IBD_folder/midas2_data/merge/genes/104158/104158.snps_info.tsv" 
#' -d "/Users/miriamgoldman/Projects/IBD_folder/midas2_data/merge/genes/104158/104158.snps_depth.tsv" 
#' -m 5 -n 5 
###check inputs given and read in data
## filter Genes make plots to explain Genes and make a table of copynumber and sample name

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
verbose<-opt$verbose
output_dir<-test_dir(opt$out_folder,verbose)
s_id<-opt$species_id
tmp <- file.path(output_dir, paste0(s_id,".prep_genes.log"))

lf <- log_open(tmp,show_notes=FALSE)
put("running validate genes command is",console = verbose)
command_arg_list<-commandArgs(trailingOnly = TRUE)
commad_message<-paste(unlist(command_arg_list),collapse = " ")

put(paste("Rscript prep_genes_command_line.R",commad_message),console = verbose)
validate_genes_output=validate_genes_input(opt)#parse_genes list
do.call(prep_genes_function_R,validate_genes_output)

# Close log
log_close()
