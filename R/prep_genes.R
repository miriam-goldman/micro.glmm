library("optparse")
suppressPackageStartupMessages(library(tidyverse))
library(data.table,quietly=TRUE,warn.conflicts=FALSE)
source("helper_functions.R")
library(magrittr,quietly=TRUE,warn.conflicts=FALSE)
option_list = list(
  make_option(c("-s", "--species_id"), type="character", default="s_id", 
              help="species id for labeling output file names (optional)", metavar="character"),
  make_option(c("-v", "--verbose"), type="logical", default=TRUE, 
              help="verbose, TRUE or FALSE"),
  make_option(c("-o", "--out_folder"), type="character", default=NULL,
              help="output folder name", metavar="character"),
  make_option(c("-c", "--genes_copynumber_file"), type="character", default=NULL,
              help="snps freq file name and path", metavar="character"),
  make_option(c("-d", "--genes_depth_file"), type="character", default=NULL,
              help="snps depth file name and path", metavar="character"),
  make_option(c("-p", "--centroid_prevalence_file"), type="character", default=NULL,
              help="centriod prevalence file (optional)" , metavar="character"),
  make_option(c("--centroid_prevalence_cutoff"), type="numeric", default=.7,
              help="centriod prevalence file (optional)" ),
  make_option(c("--genes_summary"), type="character", default=NULL,
              help="genes summary from MIDAS2 output used to filter sample list (optional)" ),
  make_option(c("--metadata"), type="character", default=NULL,
              help="metadata for samples; controls coded as 0s if using to filter (optional)" ),
  make_option(c("--min_num_control"), type="numeric", default=5,
              help="min number of samples that are controls to use gene (optional)" ),
  make_option(c("--min_num_case"), type="numeric", default=5,
              help="min number of samples that are case to use gene (optional)" ),
  make_option(c("--GRM"), type="character", default=NULL,
              help="GRM used to check sample list (optional)" ),
  make_option(c("-n","--number_of_samples_for_copynumber"), type="integer", default=5,
              help="number of samples needed for gene to be considered"),
  make_option(c("-m","--depth_cutoff"), type="integer", default=5,
              help="depth needed for gene to be included")
)

#' Rscript prep_genes.R 
#' --species_id 104158 
#' -v TRUE 
#' -o "~/Projects/micro-glmm/test_output" 
#' -c "/Users/miriamgoldman/Projects/IBD_folder/midas2_data/merge/genes/104158/104158.snps_info.tsv" 
#' -d "/Users/miriamgoldman/Projects/IBD_folder/midas2_data/merge/genes/104158/104158.snps_depth.tsv" 
#' -m 5 -n 5 




## filter Genes make plots to explain Genes and make a table of copynumber and sample name

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
verbose=opt$verbose
###check inputs given and read in data
s_id=opt$species_id
if(verbose){
  message(paste("species id for out files is:",s_id))
}
output_dir<-test_dir(opt$out_folder,verbose)

##read in gene data
if(isTRUE(!is.na(opt$genes_copynumber_file))){
  if(file_test("-f",opt$genes_copynumber_file)){
    gcopynumber<-fread(opt$genes_copynumber_file)
  }else{
    message("copy number file doesnt exist")
    stop()
  }
}else{
  message("copy number file doesnt exist")
  stop()
}

if(isTRUE(!is.na(opt$genes_copynumber_file))){
  if(file_test("-f",opt$genes_depth_file)){
    gdepth<-fread(opt$genes_depth_file)
  }else{
    message("genes depth file doesnt exist")
  }
}else{
  message("genes depth file doesnt exist")
}

if(verbose){
  message("gene files read in")
}


if(isTRUE(!is.na(opt$genes_summary))){
  if(file_test("-f",opt$genes_summary)){
    if(verbose){
      message("reading in genes summary")
    }
    genes_summary<-fread(opt$genes_summary) %>% filter(marker_coverage > 0)
    genes_summary_by_sample <- genes_summary %>% group_by(sample_name) %>% summarize(num_species=n(),avg_mapped_reads=mean(mapped_reads))
    genes_summary_used=TRUE
  }else{
    genes_summary_used=FALSE
    message("genes depth file doesnt exist")
  }
}else{
  genes_summary_used=FALSE
  message("genes depth file doesnt exist")
  }

if(isTRUE(!is.na(opt$centroid_prevalence_file))){
  if(file_test("-f",opt$centroid_prevalence_file)){
    centroid_prevalence_file<-fread(opt$centroid_prevalence_file)
    pangenome_used=TRUE
    if(verbose){
      message(paste("pangenome used and read in centriod cutoff at",centroid_prevalence_cutoff))
    }
  }else{
    pangenome_used=FALSE
  }
  }else{
  pangenome_used=FALSE
}
if(isTRUE(!is.na(opt$metadata))){
  if(file_test("-f",opt$metadata)){
    metadata <- fread(opt$metadata) 
    print(metadata)
    metadata_overlaps<-any(metadata$sample_name %in% colnames(gcopynumber))
    if(!metadata_overlaps){
      message("samples in metadata do not match MIDAS output please make sure you metadata has sample names in a colmun labeled sample_name and binary phenotypes in a column labeled disease_status")
      }else{
        metadata_used=TRUE
      if(verbose){
        message(paste("using metadata to filter, filtering controls at",opt$min_num_control,"per gene"))
      }
    }
  }else{
    metadata_used=FALSE
  }
}else{
  metadata_used=FALSE
  if(verbose){
    message("not using metadata")
  }
  
}
if(isTRUE(!is.na(opt$GRM))){
  if(file_test("-f",opt$GRM)){
    GRM <- fread(opt$GRM,sep="\t") 
    GRM_overlaps<-any(colnames(GRM) %in% colnames(gcopynumber))
    if(!GRM_overlaps){
      message("samples in GRM do not match MIDAS output please make sure your GRM has sample names as the colmun names")
    }else{
    GRM_used=TRUE
    if(verbose){
      message(paste("using GRM to filter sample list from GRM will be used"))
    }
    }
  }
}else{
  GRM_used=FALSE
  if(verbose){
    message("not using GRM to filter")
  }
  
}
samples_per_copynumber<-opt$number_of_samples_for_copynumber
if(!is.numeric(samples_per_copynumber)){
  message("samples per copynumber invlaid use number or 0")
}
depth_cutoff<-opt$depth_cutoff
if(!is.numeric(depth_cutoff)){
  message("depth cutoff invlaid use number or 0")
}
start_genes<-ncol(gcopynumber)
if(verbose){
  message("optional files read in")
  message(paste("for Species ID",s_id,
    "filtering with samples per copynumber",samples_per_copynumber,
                "depth at",depth_cutoff,
                "GRM used:",GRM_used,
                "metadata used:",metadata_used,
                "genes summary used:",genes_summary_used,
                  "starting with:", start_genes,"samples", 
                  "and with:", length(unique(gcopynumber$gene_id)),"genes"))
}


##filter by samples first order of filtering... 1) GRM, 2) genes_summary, 3) pangeome, 4) depth and number of samples per, 5) metadata

list_of_samples <- colnames(gcopynumber)[-1]

if(GRM_used){
  list_of_samples<-list_of_samples[which(list_of_samples %in% colnames(GRM))]
  if(verbose){
    message(paste("number of samples after GRM filter:",length(list_of_samples)))
  }
}
if(genes_summary_used){
  list_of_samples_gs <- genes_summary %>% filter(species_id == s_id) %>% .$sample_name %>% unique()
  list_of_samples<-list_of_samples[which(list_of_samples %in% list_of_samples_gs)]
  if(verbose){
    message(paste("number of samples after genes summary filter:",length(list_of_samples)))
  }
}
stopifnot(length(list_of_samples)>0)
gdepth %<>% gather(sample_name, gene_depth, setdiff(colnames(gdepth), "gene_id"))
gdepth %<>% filter(sample_name %in% list_of_samples)
gcopynumber %<>% gather(sample_name, copy_number, setdiff(colnames(gcopynumber), "gene_id"))
gcopynumber %<>% filter(sample_name %in% list_of_samples)
gdepth %<>% filter(gene_depth >= depth_cutoff)
if(verbose){
  print(paste("Gene-level first filter: average gene depth >=",depth_cutoff))
}

# compute gene occurrence frequency
total_sample_counts<-length(list_of_samples)
byGene <- gdepth %>% group_by(gene_id) %>% summarize(sample_counts = n()) %>% ungroup() %>%
  mutate(sample_freq = round(sample_counts / total_sample_counts, 2))
if(verbose){
  print(paste("Filter out gene is present in less than", samples_per_copynumber, "samples (for modeling purpose)"))
}

byGene %<>% filter(sample_counts >= samples_per_copynumber)



if(pangenome_used){
  byGene<-left_join(byGene,centroid_prevalence_file,by=c("gene_id"="centroid_95"))
  byGene %<>% mutate(isCore = ifelse(centroid_prevalence >= centroid_prevalence_cutoff, "core", "accessory"))
  byGene %>% ggplot(aes(x=sample_freq,y=centroid_prevalence))+  stat_bin_hex(bins=50)+labs(title=paste("For Species",s_id,"sample_freq verse centriod prevalence"))
  ggsave(file.path(output_dir, paste0(s_id,".gene_2d_histogram.pdf")), width = 7, height = 6)
  byGene %>% ggplot(aes(x = sample_freq,color=isCore)) + geom_histogram(bins = 30) + ggtitle(paste("Gene occurrence for species:", s_id))
  ggsave(file.path(output_dir, paste0(s_id,".gene_histogram.pdf")), width = 7, height = 6)
  byGene<-byGene %>% filter(!isCore)
  if(verbose){
    print(paste("labeled genes core if >=", centroid_prevalence_cutoff, "of examples for species had the gene"))
    message(paste("number of genes length after pangenome filter:",length(unique(byGene$gene_id))))
  }
}else{
  byGene %>% ggplot(aes(x = sample_freq)) + geom_histogram(bins = 30) + ggtitle(paste("Gene occurrence for species:", s_id))
  ggsave(file.path(output_dir, paste0(s_id,".gene_histogram.pdf")), width = 7, height = 6)
  
}

## depth and copy number togetehr
df <- left_join(gdepth, gcopynumber, by=c("gene_id", "sample_name"))
## filter to genes with enough samples and enough depth
df %<>% filter(gene_id %in% unique(byGene$gene_id))
## filter to genes with enough samples and enough depth
df %<>% left_join(byGene, by=c("gene_id"))
list_of_genes<-unique(byGene$gene_id)

if(metadata_used){
  model_df_input<-df %>% left_join(metadata, by=c("sample_name"))
  keep_genes<-model_df_input %>% group_by(gene_id,y) %>% 
    summarise(num_samples=n()) %>% 
    ungroup() %>% 
    spread(y,num_samples, fill=NA) %>% 
    filter(`0`>opt$min_num_control,`1`>opt$min_num_case) #filter for enough samples in control and not control
  list_of_genes<-list_of_genes[which(list_of_genes %in% keep_genes$gene_id)]
  if(verbose){
    message(paste("number of genes length after metadata filter:",length(list_of_genes)))
  }
}
copy_number_for_model<-df %>% filter(gene_id %in% list_of_genes)
df$s_id<-s_id
write.csv(copy_number_for_model,file.path(output_dir,paste0(s_id,".modeling_data.csv")))

if(verbose){
  message("gene filtering complete")
  message(paste("for Species ID",s_id,
                "filtering with samples per copynumber",samples_per_copynumber,
                "depth at",depth_cutoff,
                "GRM used:",GRM_used,
                "pangenome used:",pangenome_used,
                "metadata used:",metadata_used,
                "genes summary used:",genes_summary_used,
                "starting with:", start_genes,"samples", 
                "and with:", length(unique(gcopynumber$gene_id)),"genes",
                "ending with:", length(unique(copy_number_for_model$sample_name)),"samples", 
                "and with:", length(unique(copy_number_for_model$gene_id)),"genes"))
}


