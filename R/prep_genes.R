library("optparse")
library(tidyverse)
library(readr)
library(stringr)
library(reshape2)
library(readxl)
library(data.table)
library(Matrix)
library(pander)
library(pheatmap)
library(ggExtra)
library(pROC)
library(pander)
library(forcats)
library(grid)
library(RColorBrewer)
library(ggbeeswarm)
source("~/MGWAS/2023-05-03-run_model/functions_for_MWAS_LMM.R", echo=TRUE)
library(magrittr)
option_list = list(
  make_option(c("-s", "--species_id"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output folder name", metavar="character"),
  make_option(c("-d", "--data_dir"), type="character",default=NULL,
              help="input folder name", metavar="character"),
  make_option(c("-p","--pangenome"), type="character", default=NULL,
              help="pangenome folder", metavar="character"),
  make_option(c("-m","-- "), type="character", default=NULL,
              help="metadata file", metavar="character")
)
## filter SNPs make plots to explain SNPs and make distances

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
s_id=opt$s
midas_dir<-"/pollard/data/projects/mgwas/dbs/midasdb_uhgg/"
data_dir="/pollard/data/projects/mgwas/midas2_output/merge/genes/"
snps_dir<-"/pollard/data/projects/mgwas/midas2_output/merge/snps/"
output_dir="/pollard/data/projects/mgoldman/glmm_output/"
qp_dir="/pollard/data/projects/mgoldman/distance_output/"
starting_dir="~/MGWAS/Franza_metadata.csv"
scripts_df<-"~/MGWAS/2023-02-02-gut-genome-comparisons/scripts/"
## filter marker coverage more than 0 and find avg mapped reads
genes_summary <- read_delim(paste0(data_dir,"genes_summary.tsv"), delim = "\t", show_col_types = F)
genes_summary<-genes_summary %>% filter(marker_coverage > 0)
genes_summary_by_sample <- genes_summary %>% group_by(sample_name) %>% summarize(num_species=n(),avg_mapped_reads=mean(mapped_reads))

snps_summary<-read_delim(paste0(snps_dir,"snps_summary.tsv"), delim = "\t", show_col_types = F)
#read in metadata
metadata <- read_delim(starting_dir, delim = ",", show_col_types = F) %>%
  select(`NCBI.Accession.Number`, Diagnosis,Age,antibiotic) %>%
  mutate(Diagnosis = as.character(Diagnosis),Age=as.numeric(Age),antibiotic=as.factor(antibiotic))
colnames(metadata) <- c("sample_name", "disease_status","Age","antibiotic")

print(metadata)
species_list_filename<-paste0(scripts_df,"all_species_list.tsv")
species_list=read.csv(species_list_filename)

list_of_samples <- genes_summary %>% filter(species_id == s_id) %>% .$sample_name %>% unique()
snp_summary_by_species<-snps_summary %>% filter(species_id==s_id) %>% mutate(species_coverage=mean_coverage,s_id=species_id) %>% select(s_id,sample_name,species_coverage)
smetadata<-metadata %>% left_join(snp_summary_by_species) 
total_sample_counts <- length(list_of_samples) 
system(paste("echo For species", s_id, "there are", total_sample_counts, "samples"))
if(total_sample_counts>30){
  data_dir_s<-paste0(data_dir,s_id,"/")
  print(data_dir_s)
  data_file_list<-list.files(data_dir_s)
  if(all(grepl("lz4",data_file_list))){
    setwd(data_dir_s)
    unzip="lz4 -md *"
    system(unzip)
  }
  #add core label from UHGG
  pangeome_folder<-paste0(midas_dir,"pangenomes/")
  species_pangenome_folder<-paste0(pangeome_folder,s_id)
  core_label<-fread(paste0(species_pangenome_folder,"/centroid_prevalence.tsv"))
  data_file_list<-list.files(data_dir_s)
  #read in depth and copynumber
  gdepth <- fread(paste0(data_dir_s, s_id, ".genes_depth.tsv"), sep = "\t")
  gdepth %<>% gather(sample_name, gene_depth, setdiff(colnames(gdepth), "gene_id"))
  gdepth %<>% filter(sample_name %in% list_of_samples)
  gcopynum <- fread(paste0(data_dir_s, s_id, ".genes_copynum.tsv"), sep = "\t")
  gcopynum %<>% gather(sample_name, copy_number, setdiff(colnames(gcopynum), "gene_id"))
  gcopynum %<>% filter(sample_name %in% list_of_samples)
  
  
  print("Gene-level first filter: average gene depth >= 4")
  gdepth %<>% filter(gene_depth >= 4)
  
  
  # compute gene occurrence frequency
  byGene <- gdepth %>% group_by(gene_id) %>% summarize(sample_counts = n()) %>% ungroup() %>%
    mutate(sample_freq = round(sample_counts / total_sample_counts, 2))
  
  print("Filter out gene is present in less than 10 samples (for modeling purpose)")
  byGene %<>% filter(sample_counts >= 10)
  byGene<-left_join(byGene,core_label,by=c("gene_id"="centroid_95"))
  byGene %<>% mutate(isCore = ifelse(centroid_prevalence >= 0.7, "core", "accessory"))
  print("Let's look at the average gene depth histogram")
  byGene %>% ggplot(aes(x = sample_freq,color=isCore)) + geom_histogram() + ggtitle(paste("Gene occurrence for species:", s_id))
  ggsave(file.path(output_dir, paste0("/",s_id,"/","gene_histogram", s_id, ".pdf")), width = 7, height = 6)
  byGene %>% ggplot(aes(x=sample_freq,y=centroid_prevalence))+  stat_bin_hex(bins=50)+labs(title=paste("For Species",s_id,"sample_freq verse centriod prevalence"))
  ggsave(file.path(output_dir, paste0("/",s_id,"/","gene_2d_histogram", s_id, ".pdf")), width = 7, height = 6)
  ## depth and copy number togetehr
  df <- left_join(gdepth, gcopynum, by=c("gene_id", "sample_name"))
  ## filter to genes with enough samples and enough depth
  df %<>% filter(gene_id %in% unique(byGene$gene_id))
  ## filter to genes with enough samples and enough depth
  df %<>% left_join(byGene, by=c("gene_id"))
  
  ## add in metadata
  df %<>% left_join(smetadata, by=c("sample_name"))
  
  # filter out antiboitic
  model_df_input<-df %>% filter(antibiotic!=1) %>%ungroup()  %>% mutate(s_id=s_id) 
  total_samples=sum(!duplicated(model_df_input$sample_name))
  num_control=sum(smetadata$disease_status==0)
  #filter on number control for each gene
  if(num_control>5 & num_control!=total_samples){
    keep_genes_2<-df %>% group_by(gene_id,disease_status) %>% summarise(num_samples=n()) %>% ungroup() %>% spread(disease_status,num_samples, fill=NA) %>% mutate(not_control=`1`+`2`)%>%  filter(`0`>7,not_control>7) #filter for enough samples in control and not control
    df<-df %>% right_join(keep_genes_2)
    model_df_input<-model_df_input %>% filter(gene_id %in% keep_genes_2$gene_id)
    #filter to only accessory genes
    model_df_input<-model_df_input %>% ungroup() %>% filter(isCore=="accessory") %>%  mutate(disease_status=as.numeric(ifelse(disease_status>=1,1,0)))#<==========
    print(nrow(model_df_input))
    
    model_df_input<-model_df_input %>% filter(!is.na(s_id))
    
    print(nrow(model_df_input))
    
    if(nrow(model_df_input)>100){
      #s_id=101338
      s_info<-species_list %>% filter(species_id==s_id) %>% select(species_id,D,P,C,O,F,G,S)
      print(paste("species",s_id))
      GRM = data.table:::fread(paste0(qp_dir,s_id,"/by_sample_distances", s_id, ".tsv"),sep="\t")
      
      if(ncol(GRM)>30){
        rownames(GRM)<-GRM$V1
        GRM<-GRM[-1,]
        GRM<-GRM %>% select(-ref,-V1)
        GRM<-1-GRM
        GRM<-Matrix(as.matrix(GRM))
        dimnames(GRM)<-c(dimnames(GRM)[2],dimnames(GRM)[2])
        data.for.model.final=smetadata %>% filter(sample_name %in% colnames(GRM)) %>% select(sample_name,disease_status,Age,species_coverage) %>% mutate(disease_status=ifelse(as.numeric(disease_status)>0,1,0)) # n (for species from CNV) x 3
        
        rownames(data.for.model.final)<-data.for.model.final$sample_name
        data.for.model.final<-data.for.model.final[order(colnames(GRM)),]
        #copy number matches GRM
        
        copy_number_for_model<-model_df_input %>% filter(sample_name %in% data.for.model.final$sample_name)
        write.csv(copy_number_for_model,paste0(output_dir,s_id,"/",s_id,"_modeling_data.csv"))
        
      }
    }
  }
}