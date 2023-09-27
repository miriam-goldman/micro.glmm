library("optparse")
library(tidyverse)
library(data.table)
library(magrittr)
library(knitr)
library(pander)
library(purrr)
library(ggExtra)
library(parallel)
library(paws)
library(logr)
library(RcppXPtrUtils)
library(pheatmap)
library(parallelDist)
#source("~/R/filtered_fast.R")
library(ape)

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


option_list = list(
  make_option(c("-s", "--species_id"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output folder name", metavar="character"),
  make_option(c("-d", "--data_dir"), type="character",default=NULL,
              help="input folder name", metavar="character"),
  make_option(c("-p","--pangenome"), type="character", default=NULL,
              help="pangenome folder", metavar="character")
)

## filter SNPs make plots to explain SNPs and make distances

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

s_id=opt$s
print(s_id)
unzip='lz4 -m -d *.lz4'
zip='rm *.tsv'
print("read in files")
data_dir=opt$d
out_dir=opt$o
setwd(paste0(data_dir,"/",s_id))
system(unzip)
midas_folder=opt$p
print("unzip lz4")
print("read info")
info<-fread(paste0(data_dir,"/",s_id,"/",s_id,".snps_info.tsv"))
print("read depth")
depth<-fread(paste0(data_dir,"/",s_id,"/",s_id,".snps_depth.tsv"))
print("read freq")
freq<-fread(paste0(data_dir,"/",s_id,"/",s_id,".snps_freqs.tsv"))
total_sites=nrow(freq)
system(zip)
print("removed tsv")
print("data read in")
info %>% dplyr::count(snp_type, locus_type, site_type) %>% 
  spread(site_type, n, fill = 0) %>%
  pander() %>% print()


centriod_to_ref<-fread(paste0(midas_folder,s_id,"/centroid_to_repgenes.tsv"))
core_label<-fread(paste0(midas_folder,s_id,"/centroid_prevalence.tsv"))
ref_core_labeled<-left_join(centriod_to_ref,core_label)
ref_core_labeled<-ref_core_labeled %>% mutate(gene_id=str_replace(rep_gene_id,"UHGG","GUT_GENOME"))
info<-info %>% left_join(ref_core_labeled)
info$core<-info$centroid_prevalence>.7 #<----------------------------------------------------------
######## Filter Samples based on D_median_cds >= 10
# Compute median site depth for all protein coding genes
depth <- depth %>%
  filter(site_id %in% unique(info$site_id))
nsamples = ncol(depth)-1
print(paste("number of samples",nsamples))
# Site not covered in a certain <sample, species> pair is site-depth = 0, and it should not be included in the computation.
depth[depth == 0] <- NA
per_sample_median_depth <- apply(depth[,-1], 2, function(x) median(x, na.rm =T))
samples_pass_depth <- names(per_sample_median_depth[per_sample_median_depth >= 10]) #<-----
depth %<>% select(site_id, all_of(samples_pass_depth))
nsamples2 = ncol(depth)-1
print(nsamples2)
D <- data.frame(sample_name = names(per_sample_median_depth), median_site_depth=per_sample_median_depth, row.names=NULL)

g2 <- D %>% 
  ggplot(aes(x = median_site_depth)) + geom_histogram() + geom_vline(xintercept = 10, color = "red") + 
  scale_x_log10() + 
  ggtitle(paste(s_id, ": ", nsamples2, " out of ", nsamples, " samples passing the 10X filter", sep=""))
print("number of samples passing filter plot")
print(g2)
ggsave(file.path(out_dir, paste0("D_hist_", s_id, ".pdf")), width = 7, height = 6)


######### Keep sites [0.3, 3] * D_median_cds
info_for_qp <- info %>% filter(site_type == "4D")

depth_for_qp<-depth %>% filter(site_id %in% unique(info_for_qp$site_id)) 
depth_for_qp %<>% gather(sample_name, site_depth, all_of(samples_pass_depth))

depth_for_qp %<>% 
  left_join(D %>% select(sample_name, median_site_depth)) %>%
  mutate(min_bound = 0.3 * median_site_depth, max_bound = 3 * median_site_depth)

depth_for_qp %<>%
  filter(site_depth >= min_bound & site_depth <= max_bound) %>%
  filter(site_depth >= 5) %>% #<-------------
select(site_id, sample_name, site_depth, median_site_depth)


sc_df <- depth_for_qp %>%
  group_by(site_id) %>% 
  summarise(sample_counts = n(),total_reads=sum(site_depth)) %>%
  ungroup()
sc_df$pos<-str_split_i(sc_df$site_id,"\\|",4)
sc_df$contig<-str_split_i(sc_df$site_id,"\\|",3)
sites_no <- unique(sc_df %>% filter(sample_counts <= 5) %>% .$site_id)
stopifnot(nrow(sites_no) == 0)
if (length(sites_no) > 0) {
  depth_for_qp %<>% filter(!site_id %in% sites_no)
}


######### Read in population minor allele frequency
freq_for_qp <- freq%>%
  filter(site_id %in% unique(depth_for_qp$site_id)) %>%
  select(site_id, all_of(samples_pass_depth))

freq_for_qp %<>%
  gather(sample_name, allele_freq, all_of(samples_pass_depth)) %>%
  filter(allele_freq != -1)


df <- left_join(depth_for_qp, freq_for_qp, by=c("site_id", "sample_name"))
df %<>% 
  mutate(allele_direction = ifelse(allele_freq <= 0.2, "low", "int")) %>%
  mutate(allele_direction = ifelse(allele_freq >= 0.8, "high", allele_direction))


f2 <- file.path(out_dir, paste("/","/site_depth_with_allele_freq_", s_id, ".pdf", sep=""))
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
print(head(info_for_qp))
print(head(gene_info_1))
gene_info_1<-gene_info_1 %>% group_by(gene_id) %>% summarize(mean_sample_count=mean(sample_counts.x),sum_total_reads=sum(total_reads),mean_total_reads=mean(total_reads), centroid_prevalence=mean(centroid_prevalence))

gene_plot<-gene_info_1 %>% ggplot(aes(mean_sample_count,centroid_prevalence,color=mean_total_reads))+geom_point()+labs(title=paste("mean sample count for gene verse centroid prevalence for species",s_id))
gene_plot_1<-ggExtra::ggMarginal(gene_plot, type = "histogram")
print(gene_plot_1,newpage = TRUE)
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

f3 <- file.path(out_dir, paste("/","QP_hist_", s_id, ".pdf", sep=""))
print("QP histogram")
g3 <- by_sample %>%
  ggplot(aes(x = QP)) + geom_histogram() +
  theme_bw() + geom_vline(xintercept = 0.1, color = "red")


print(g3)


ggsave(f3, g3, width = 7, height = 5)
system(paste("echo ",file.path(out_dir, paste("by_sample_", s_id, ".tsv", sep = ""))))
by_sample %>%
  write.table(file.path(out_dir, paste("/","by_sample_", s_id, ".tsv", sep = "")), sep = "\t", quote = F, row.names = F)
print(head(by_sample))
#### re-filter everything so we arent just looking at 4D sites

info_for_distance<-info %>% filter(core==TRUE) %>% filter(snp_type=="bi")
system(paste("echo",colnames(info)))
depth_for_distance<-depth %>% filter(site_id %in% unique(info_for_distance$site_id)) 
depth_for_distance %<>% gather(sample_name, site_depth, all_of(samples_pass_depth))

depth_for_distance %<>% 
  left_join(D %>% select(sample_name, median_site_depth)) %>%
  mutate(min_bound = 0.3 * median_site_depth, max_bound = 3 * median_site_depth)

depth_for_distance %<>%
  #filter(site_depth >= min_bound & site_depth <= max_bound) %>%
  filter(site_depth >= 3) %>% #<-------------
select(site_id, sample_name, site_depth, median_site_depth)


sc_df_2 <- depth_for_distance %>%
  group_by(site_id) %>% 
  summarise(sample_counts = n(),total_reads=sum(site_depth)) %>%
  ungroup()
sc_df_2$pos<-str_split_i(sc_df_2$site_id,"\\|",4)
sc_df_2$contig<-str_split_i(sc_df_2$site_id,"\\|",3)
print("sample depth historgram")
sample_counts_per_site<-sc_df_2 %>% ggplot(aes(x = sample_counts)) + geom_histogram()

sites_no <- unique(sc_df_2 %>% filter(sample_counts <= 5) %>% .$site_id)
stopifnot(nrow(sites_no) == 0)
if (length(sites_no) > 0) {
  depth_for_distance %<>% filter(!site_id %in% sites_no)
}


######### Read in population minor allele frequency
freq_for_distance <- freq%>%
  filter(site_id %in% unique(depth_for_distance$site_id)) %>%
  select(site_id, all_of(samples_pass_depth))

freq_for_distance %<>%
  gather(sample_name, allele_freq, all_of(samples_pass_depth)) %>%
  filter(allele_freq != -1)


df_for_distance <- left_join(depth_for_distance, freq_for_distance, by=c("site_id", "sample_name"))
site_df<-df_for_distance %>% select(site_id,sample_name,allele_freq) %>% pivot_wider(names_from = site_id,values_from=allele_freq)
print(head(site_df))
print(paste("num of samples",nrow(site_df)))
print(paste("num of snps",ncol(site_df)-1))
freq_mat_dist_man<-parDist(as.matrix(site_df[,-1]), method="custom", func = manhattanFuncPtr)
freq_mat_dist_man<-as.matrix(freq_mat_dist_man)
dimnames(freq_mat_dist_man)<-c(site_df[,1],site_df[,1])
annotation_col = data.frame(qp=as.character(by_sample$QP<.1))
rownames(annotation_col) = by_sample$sample_name


ann_colors = list(
  qp = c("TRUE" = "#1B9E77", "FALSE" = "#D95F02")
)

pheatmap(freq_mat_dist_man,annotation_row = annotation_col,annotation_col=annotation_col,annotation_colors = ann_colors)
ref_df<-rep(0,ncol(site_df))
names(ref_df)<-colnames(site_df)
ref_df[1]<-"ref"
site_df_with_ref<-rbind(ref_df,site_df)
site_df_with_ref[,-1]<-apply(site_df_with_ref[,-1],2,as.numeric)
freq_mat_dist_man_with_ref<-parDist(as.matrix(site_df_with_ref[,-1]), method="custom", func = manhattanFuncPtr)
freq_mat_dist_man_with_ref<-as.data.frame(as.matrix(freq_mat_dist_man_with_ref))
dimnames(freq_mat_dist_man_with_ref)<-c(site_df_with_ref[,1],site_df_with_ref[,1])
freq_mat_dist_man_with_ref  %>%  write.table(file.path(out_dir, paste("/","distance.", s_id, ".tsv", sep = "")), sep = "\t", quote = F)
f4 <- file.path(out_dir, paste("/","distance_heatmap_", s_id, ".pdf", sep=""))
pheatmap(freq_mat_dist_man,annotation_row = annotation_col,annotation_col=annotation_col,filename=f4)

pcoa_data=pcoa(freq_mat_dist_man)
biplot(pcoa_data)
projection_data<-as.data.frame(pcoa_data$vectors)
projection_data$sample_name<-rownames(projection_data)
projection_data<-left_join(by_sample,projection_data)
f5 <- file.path(out_dir, paste("/","PcoA_plot_", s_id, ".pdf", sep=""))
print("PcoA plot")
g5 <-projection_data %>%
  ggplot(aes(Axis.1,Axis.2,shape=as.factor(QP<.1),size=2))+
  geom_point()+labs(x=paste("PC 1 variance explained",round(pcoa_data$values$Relative_eig[1],2)),y=paste("PC 2 variance explained",round(pcoa_data$values$Relative_eig[2],2)),title=paste("PcoA for species",s_id)) + 
  theme_bw() 
ggsave(f5,g5, width = 7, height = 5)
print("PC plot")
projection_data<-princomp(freq_mat_dist_man)
projection_data_ouput<-as.data.frame(as.matrix(projection_data[["loadings"]][,]))
projection_data_ouput$sample_name<-rownames(projection_data_ouput)
projection_data_ouput<-left_join(by_sample,projection_data_ouput)
proportion_of_var<-projection_data$sdev^2/sum(projection_data$sdev^2)
projection_plot<-projection_data_ouput %>%
  ggplot(aes(Comp.1,Comp.2,shape=as.factor(QP<.1),size=2))+
  geom_point()+labs(x=paste("PC 1 variance explained",round(proportion_of_var[1],2)),y=paste("PC 2 variance explained",round(proportion_of_var[2],2)),title=paste("PCA for species",s_id)) + 
  theme_bw() 
print(projection_plot)

print("MDS plot")
projection_data<-cmdscale(freq_mat_dist_man,add=TRUE,eig=TRUE)
projection_data_ouput<-as.data.frame(as.matrix(projection_data[["points"]]))
colnames(projection_data_ouput)<-c("comp1","comp2")
projection_data_ouput$sample_name<-rownames(projection_data_ouput)
projection_data_ouput<-left_join(by_sample,projection_data_ouput)
proportion_of_var<-projection_data$eig^2/sum(projection_data$eig^2)
projection_plot<-projection_data_ouput %>%
  ggplot(aes(comp1,comp2,shape=as.factor(QP<.1),size=2))+
  geom_point()+labs(x=paste("PC 1 variance explained",round(proportion_of_var[1],2)),y=paste("PC 2 variance explained",round(proportion_of_var[2],2)),title=paste("MDS for species",s_id)) + 
  theme_bw() 
print(projection_plot)
