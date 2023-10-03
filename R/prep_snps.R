library("optparse",quietly=TRUE)
suppressPackageStartupMessages(library(tidyverse,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE))
library(data.table,quietly=TRUE,warn.conflicts=FALSE)
library(magrittr,quietly=TRUE,warn.conflicts=FALSE)
library(pander,quietly=TRUE,warn.conflicts=FALSE)
library(purrr,quietly=TRUE,warn.conflicts=FALSE)
library(ggExtra,quietly=TRUE,warn.conflicts=FALSE)
library(parallel,quietly=TRUE,warn.conflicts=FALSE)
library(paws,quietly=TRUE,warn.conflicts=FALSE)
library(logr,quietly=TRUE,warn.conflicts=FALSE)
library(RcppXPtrUtils,quietly=TRUE,warn.conflicts=FALSE)
library(pheatmap,quietly=TRUE,warn.conflicts=FALSE)
library(parallelDist,quietly=TRUE,warn.conflicts=FALSE)
source("helper_functions.R")
library(ape,quietly=TRUE,warn.conflicts=FALSE)

manhattanFuncPtr <- cppXPtr(
  "double customDist(const arma::mat &A, const arma::mat &B) {
    float dist_com=0;
    int n_sites=0;
    for (int i = 0; i < A.size(); i++){
      if(A[i]>=0 & B[i] >=0){
        dist_com += fabs(A[i]-B[i]);
        n_sites++;
      }
	          }
	          return dist_com/n_sites;
  }", depends = c("RcppArmadillo"))

alleleFuncPtr <- cppXPtr(
  "double customDist(const arma::mat &A, const arma::mat &B) {
    double allele_match=0;
    int n_sites=0;
    for (int i = 0; i < A.size(); i++){
      if(A[i]>=0 & B[i] >=0){
          n_sites++;
          if(fabs(A[i]-B[i])>.6){
            allele_match++;
          }
      }
	   }
	          return allele_match/n_sites;
  }", depends = c("RcppArmadillo"))




jaccardFuncPtr <- cppXPtr(
  "double customDist(const arma::mat &A, const arma::mat &B) {
    double num_a=arma::accu(A>=0);
    double num_b=arma::accu(B>=0);
    double num_both=arma::accu((B>=0 && A>=0));
    double dist;
    dist=1-(num_both/(num_a+num_b-num_both));
    return dist;
  }", depends = c("RcppArmadillo"))

##Input needed
#output folder
#freq file
#info file
#depth file
#pangenome info (optional)
##centroid_to_repgenes.tsv
##centroid_prevalence.tsv
#calc qp (optional)
#species id (optional)
#centriod upper filter
#centriod lower file
#median depth filter
#absoulite depth filter


option_list = list(
  make_option(c("-s", "--species_id"), type="character", default="s_id", 
              help="species id for labeling output file names (optional)", metavar="character"),
  make_option(c("-v", "--verbose"), type="logical", default=TRUE, 
              help="verbose, TRUE or FALSE"),
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
              help="centriod prevalence file (optional)" ),
  make_option(c("-q", "--run_qp"), type="logical", default=FALSE, 
              help="option to run qp or not (optional)"),
  make_option(c("-u", "--median_upper_filter"), type="integer", default=3,
              help="upper limit on median filter num of times median value of depth to let pass for each sample for QP (optional)"),
  make_option(c("-l", "--median_lower_filter"), type="numeric", default=.3,
              help="lower limit on median filter num of times median value of depth to let pass for each sample for QP (optional)"),
  make_option(c("-a", "--abosulte_filter"), type="numeric", default=3,
              help="lowest depth to all to pass filters"),
  make_option(c("-m", "--sample_median_depth_filter"), type="integer", default=10,
              help="median depth needed for sample to be considered to have the species"),
  make_option(c("-n","--number_of_samples_for_sites"), type="integer", default=5,
              help="number of samples needed for site to be considered to have the species")
  
)

## filter SNPs make plots to explain SNPs and make distances

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

verbose=opt$verbose
###check inputs given and read in data
s_id=opt$species_id
if(verbose){
  message(paste("species id for out files is:",s_id))
}





output_dir<-test_dir(opt$out_folder,verbose)

if(file_test("-f",opt$snps_info_file)){
  snp_info<-fread(opt$snps_info_file)
}else{
  message("snp info file doesnt exist")
}

if(file_test("-f",opt$snps_freq_file)){
  snp_freq<-fread(opt$snps_freq_file)
}else{
  message("snp freq file doesnt exist")
}

if(file_test("-f",opt$snps_depth_file)){
  snp_depth<-fread(opt$snps_depth_file)
}else{
  message("snp depth file doesnt exist")
}
if(verbose){
  message("SNP files read in")
}

run_pangenome<-isTRUE(!is.na(opt$centroid_to_repgenes) & !is.na(opt$centroid_prevalence))
if(run_pangenome){
    if(verbose){
      message("running pangeome step")
    }
    if(file_test("-f",opt$centroid_to_repgenes)){
      centroid_to_repgenes<-fread(opt$centroid_to_repgenes)
    }else{
      message("centriod to rep gene doesnt exist")
    }
    if(file_test("-f",opt$centroid_prevalence)){
      core_label<-fread(opt$centroid_prevalence)
    }else{
      message("centriod prevalence doesnt exist")
    }
    ref_core_labeled<-left_join(centroid_to_repgenes,core_label)
    snp_info<-snp_info %>% left_join(ref_core_labeled)
    snp_info$core<-snp_info$centroid_prevalence>opt$centroid_prevalence_cutoff
    if(verbose){
      message("core labled with centriod prevalnce")
    }
  }else{
    print("not running pangenome")
    snp_info$core<-snp_info$locus_type=="CDS"
    message("using snps in coding region")
  }
  

if(verbose){
  snp_info %>% dplyr::count(snp_type, locus_type, site_type) %>% 
    spread(site_type, n, fill = 0) %>%
    pander() %>% print()
}


######## Filter Samples based on D_median_cds 
# Compute median site depth for all protein coding genes
snp_depth <- snp_depth %>%
  filter(site_id %in% unique(snp_info$site_id))
nsamples = ncol(snp_depth)-1
if(verbose){
  print(paste("number of samples",nsamples))
}
# Site not covered in a certain <sample, species> pair is site-depth = 0, and it should not be included in the computation.
snp_depth[snp_depth == 0] <- NA
per_sample_median_depth <- apply(snp_depth[,-1], 2, function(x) median(x, na.rm =T))
samples_pass_depth <- names(per_sample_median_depth[per_sample_median_depth >= opt$sample_median_depth_filter]) #<-----
snp_depth %<>% select(site_id, all_of(samples_pass_depth))
nsamples2 = ncol(snp_depth)-1
if(verbose){
  print(paste("number of samples after filter",nsamples2))
}
D <- data.frame(sample_name = names(per_sample_median_depth), median_site_depth=per_sample_median_depth, row.names=NULL)

g2 <- D %>% 
  ggplot(aes(x = median_site_depth)) + geom_histogram(bins=30) + geom_vline(xintercept = opt$sample_median_depth_filter, color = "red") + 
  scale_x_log10() + 
  ggtitle(paste(s_id, ":", nsamples2, "out of", nsamples, "samples passing the filter at:",opt$sample_median_depth_filter))
if(verbose){
  print("number of samples passing filter plot")
}

ggsave(file.path(output_dir, paste0(s_id,".depth_histogram.pdf")),g2, width = 7, height = 6)


######### Keep sites [0.3, 3] * D_median_cds

if(is.numeric(opt$abosulte_filter)){
  a=opt$abosulte_filter
}
if(is.numeric(opt$number_of_samples_for_sites)){
  number_of_samples_for_sites=opt$number_of_samples_for_sites
}

if(opt$run_qp){
  
  if(is.numeric(opt$median_lower_filter)){
    l=opt$median_lower_filter
  }
  if(is.numeric(opt$median_upper_filter)){
    u=opt$median_upper_filter
  }
  
  info_for_qp <- snp_info %>% filter(site_type == "4D")
  depth_for_qp<-snp_depth %>% filter(site_id %in% unique(info_for_qp$site_id)) 
  depth_for_qp %<>% gather(sample_name, site_depth, all_of(samples_pass_depth))
  depth_for_qp %<>% 
    left_join(D %>% select(sample_name, median_site_depth)) %>%
    mutate(min_bound =  l* median_site_depth, max_bound = u * median_site_depth)
  depth_for_qp %<>%
    filter(site_depth >= min_bound & site_depth <= max_bound) %>%
    filter(site_depth >= a) %>% #<-------------
  select(site_id, sample_name, site_depth, median_site_depth)
  sc_df <- depth_for_qp %>%
    group_by(site_id) %>% 
    summarise(sample_counts = n(),total_reads=sum(site_depth)) %>%
    ungroup()
  sc_df$pos<-str_split_i(sc_df$site_id,"\\|",4)
  sc_df$contig<-str_split_i(sc_df$site_id,"\\|",3)
  sites_no <- unique(sc_df %>% filter(sample_counts <= number_of_samples_for_sites) %>% .$site_id)
  stopifnot(nrow(sites_no) == 0)
  if (length(sites_no) > 0) {
    depth_for_qp %<>% filter(!site_id %in% sites_no)
  }
  freq_for_qp <- snp_freq%>%
    filter(site_id %in% unique(depth_for_qp$site_id)) %>%
    select(site_id, all_of(samples_pass_depth))
  ######### Read in population minor allele frequency
  
  
  freq_for_qp %<>%
    gather(sample_name, allele_freq, all_of(samples_pass_depth)) %>%
    filter(allele_freq != -1)
  df <- left_join(depth_for_qp, freq_for_qp, by=c("site_id", "sample_name"))
  df %<>% 
    mutate(allele_direction = ifelse(allele_freq <= 0.2, "low", "int")) %>%
    mutate(allele_direction = ifelse(allele_freq >= 0.8, "high", allele_direction))
  
  
  f2 <- file.path(output_dir, paste0(s_id,".site_depth_with_allele_freq.pdf"))
  t1 <- df %>% select(site_depth, allele_freq) %>%
    group_by(site_depth, allele_freq) %>%
    dplyr::count() %>%
    dplyr::rename(occurence = n) %>%
    ungroup()
  ######### Compute ND and QP
  by_site <- df %>% 
    group_by(site_id) %>%
    summarise(n_low = sum(allele_direction == "low"), n_high = sum(allele_direction == "high")) %>%
    ungroup() %>%
    mutate(voted_direction = ifelse(n_low > n_high, "low", "high")) %>%
    mutate(voted_counted = ifelse(n_low > n_high, n_low, n_high))
  by_site %<>% left_join(sc_df, by=c("site_id")) %>%
    mutate(voted_freq = round(voted_counted / sample_counts, 4)) %>%
    select(-voted_counted) %>%
    select(site_id, sample_counts, everything())
  gene_info_1<-inner_join(by_site,info_for_qp,by=c("site_id"))
  if(verbose){
    print(head(info_for_qp))
    print(head(gene_info_1))
  }
  gene_info_1<-gene_info_1 %>% group_by(gene_id) %>% summarize(mean_sample_count=mean(sample_counts.x),sum_total_reads=sum(total_reads),mean_total_reads=mean(total_reads), centroid_prevalence=mean(centroid_prevalence))
  gene_plot<-gene_info_1 %>% ggplot(aes(mean_sample_count,centroid_prevalence,color=mean_total_reads))+geom_point()+labs(title=paste("mean sample count for gene verse centroid prevalence for species",s_id))
  gene_plot_1<-ggExtra::ggMarginal(gene_plot, type = "histogram")
  print(gene_plot_1,newpage = TRUE)
  f5 <- file.path(output_dir, paste(s_id,"gene_hist.pdf", sep=""))
  ggsave(f5, gene_plot_1, width = 7, height = 5)
  by_sample <- df %>%
    group_by(sample_name) %>%
    summarise(n_int = sum(allele_direction == "int"), n_total = n()) %>%
    ungroup()
  
  denominator  <- df %>%
    filter(allele_direction != "int") %>%
    left_join(by_site %>% select(site_id, voted_direction, voted_freq), by=c("site_id")) %>%
    mutate(panel_freq = ifelse(allele_direction == voted_direction, voted_freq, 1-voted_freq))
  
  nd <- denominator %>% 
    group_by(sample_name) %>%
    summarise(Nd = sum(panel_freq)) %>%
    ungroup()
  
  by_sample %<>% left_join(nd, by=c("sample_name")) %>%
    mutate(QP = round(n_int / Nd, 3))
  
  f3 <- file.path(output_dir, paste(s_id,"QP_hist.pdf", sep=""))
  g3 <- by_sample %>%
    ggplot(aes(x = QP)) + geom_histogram(bins=30) +
    theme_bw() + geom_vline(xintercept = 0.1, color = "red")
  print(g3)
  ggsave(f3, g3, width = 7, height = 5)
  if(verbose){
    print(file.path(output_dir, paste0(s_id,"by_sample.tsv")))
  }
  by_sample %>%
    write.table(file.path(output_dir, paste0(s_id,"by_sample.tsv")), sep = "\t", quote = F, row.names = F)
  if(verbose){
    message("QP completed")
    print(head(by_sample))
  }
}

#### re-filter everything so we arent just looking at 4D sites

info_for_distance<-snp_info %>% filter(core==TRUE) %>% filter(snp_type=="bi")
depth_for_distance<-snp_depth %>% filter(site_id %in% unique(info_for_distance$site_id)) 
depth_for_distance %<>% gather(sample_name, site_depth, all_of(samples_pass_depth))

depth_for_distance %<>% 
  left_join(D %>% select(sample_name, median_site_depth))

depth_for_distance %<>%
  filter(site_depth >= a) %>% #<-------------
select(site_id, sample_name, site_depth, median_site_depth)


sc_df_2 <- depth_for_distance %>%
  group_by(site_id) %>% 
  summarise(sample_counts = n(),total_reads=sum(site_depth)) %>%
  ungroup()
sc_df_2$pos<-str_split_i(sc_df_2$site_id,"\\|",4)
sc_df_2$contig<-str_split_i(sc_df_2$site_id,"\\|",3)
if(verbose){
  print("sample depth historgram")
}
sample_counts_per_site<-sc_df_2 %>% ggplot(aes(x = sample_counts)) + geom_histogram(bins=30)+labs(title="sample depth histogram for sites")

f4 <- file.path(output_dir, paste(s_id,".site_hist.pdf", sep=""))
ggsave(f4, sample_counts_per_site, width = 7, height = 5)

sites_no <- unique(sc_df_2 %>% filter(sample_counts <= number_of_samples_for_sites) %>% .$site_id)
stopifnot(nrow(sites_no) == 0)
if (length(sites_no) > 0) {
  depth_for_distance %<>% filter(!site_id %in% sites_no)
}


######### Read in population minor allele frequency
freq_for_distance <- snp_freq%>%
  filter(site_id %in% unique(depth_for_distance$site_id)) %>%
  select(site_id, all_of(samples_pass_depth))

freq_for_distance %<>%
  gather(sample_name, allele_freq, all_of(samples_pass_depth)) %>%
  filter(allele_freq != -1)


df_for_distance <- left_join(depth_for_distance, freq_for_distance, by=c("site_id", "sample_name"))
site_df<-df_for_distance %>% select(site_id,sample_name,allele_freq) %>% pivot_wider(names_from = site_id,values_from=allele_freq)
if(verbose){
  print(head(site_df))
  print(paste("num of samples",nrow(site_df)))
  print(paste("num of snps",ncol(site_df)-1))
}

freq_mat_dist_man<-parDist(as.matrix(site_df[,-1]), method="custom", func = manhattanFuncPtr)
freq_mat_dist_man<-as.matrix(freq_mat_dist_man)

dimnames(freq_mat_dist_man)<-c(site_df[,1],site_df[,1])
freq_mat_GRM_man<-1-freq_mat_dist_man
freq_mat_GRM_man  %>%  write.table(file.path(output_dir, paste0(s_id,".GRM.tsv")), sep = "\t", quote = F)
freq_mat_GRM_man  %>%  write.table(file.path(output_dir, paste0(s_id,".distance.tsv")), sep = "\t", quote = F)

heatmap_file_name<-file.path(output_dir, paste0(s_id,".heatmap.pdf"))
if(opt$run_qp){
  annotation_col = data.frame(qp=as.character(by_sample$QP<.1))
  rownames(annotation_col) = by_sample$sample_name
  ann_colors = list(
    qp = c("TRUE" = "#1B9E77", "FALSE" = "#D95F02")
  )
  
  pheatmap(freq_mat_dist_man,annotation_row = annotation_col,annotation_col=annotation_col,annotation_colors = ann_colors,filename=heatmap_file_name)
  
}else{
  pheatmap(freq_mat_dist_man,filename=heatmap_file_name)
  
}