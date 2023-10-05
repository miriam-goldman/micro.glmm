###Run structue test
#after prep_snps

source("glmm_functions.R")
library("optparse")
library(tidyverse)
source("helper_functions.R")
library(logr)
library(data.table)
library(Matrix)
library(generics)
option_list = list(
  make_option(c("-s", "--species_id"), type="character", default="s_id", 
              help="species id for labeling output file names (optional)", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default="s_id", 
              help="file path to metadata", metavar="character"),
  make_option(c("-v", "--verbose"), type="logical", default=TRUE, 
              help="verbose, TRUE or FALSE"),
  make_option(c("-o", "--out_folder"), type="character", default=NULL,
              help="output folder name", metavar="character"),
  make_option(c("--GRM"), type="character", default=NULL,
              help="GRM used to check sample list (optional)" ),
  make_option(c("--tau"), type="numeric", default=0,
              help="initail values of tau (optional)" ),
  make_option(c("--maxiter"), type="numeric", default=100,
              help="max iteration of fitting parameters (optional)" ),
  make_option(c("--tol"), type="numeric", default=.0001,
              help="tolerance for parameters (optional)" ),
  make_option(c("--formula_to_fit"), type="character", default="y~1",
              help="formula used to fit population structure; the phenotype and any covariates in the metadata file that are to be used" )
 
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
verbose<-opt$verbose
output_dir<-test_dir(opt$out_folder,verbose)
s_id<-opt$species_id
tmp <- file.path(output_dir, paste0(s_id,".pop_structure_test.log"))

lf <- log_open(tmp,show_notes=FALSE)
put("running validate genes command is",console = verbose)
command_arg_list<-commandArgs(trailingOnly = TRUE)
commad_message<-paste(unlist(command_arg_list),collapse = " ")

put(paste("Rscript pop_structure_test_command_line.R",commad_message),console = verbose)

if(isTRUE(!is.na(opt$GRM))){
  if(file_test("-f",opt$GRM)){
    GRM <- fread(opt$GRM,sep="\t",header=FALSE) 
    colnames(GRM)<-c("sample_name",GRM$V1)
    setindexv(GRM,'sample_name')
    GRM<-GRM %>% select(-sample_name)
    stopifnot(ncol(GRM)==(nrow(GRM)))
    GRM<-Matrix(as.matrix(GRM))
    dimnames(GRM)<-c(dimnames(GRM)[2],dimnames(GRM)[2])
  }else{
    put("GRM invalid",console = verbose)
    stop()
  }
}else{
  put("GRM invalid",console = verbose)
  stop()
}

if(isTRUE(!is.na(opt$metadata))){
  if(file_test("-f",opt$metadata)){
    metadata <- fread(opt$metadata) 
    metadata_overlaps<-intersect(metadata$sample_name,colnames(GRM))
    if(isFALSE(length(metadata_overlaps)>0)){
      put("samples in metadata do not match MIDAS output please make sure you metadata has sample names
            in a colmun labeled sample_name and binary phenotypes in a column labeled disease_status",console = verbose)
      stop()
    }else{
      metadata<-metadata %>% filter(sample_name %in% metadata_overlaps)
      GRM<-GRM[metadata_overlaps,metadata_overlaps]
      print(metadata)
    }
  }else{
    put("metadata invalid",console = verbose)
    stop()
  }
}else{
  put("metadata invalid",console = verbose)
  stop()
  }

glm_fit0 = glm(opt$formula_to_fit, data = metadata, family = "binomial")

pop_structure_test(glm_fit0, GRM,species_id=s_id,tau=c(0,opt$tau),maxiter =opt$maxiter, verbose = verbose,tol=opt$tol)


#validate_genes_output=validate_genes_input(opt)#parse_genes list
#do.call(prep_genes_function_R,validate_genes_output)
# Close log
log_close()



