suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(tidyverse))
source("../R/helper_functions.R")
suppressPackageStartupMessages(library(logr))
suppressPackageStartupMessages(library(data.table))
library("ggExtra")
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
              help="alpha value for marker test plot", metavar="numeric"),
  make_option(c("--compare_to_glm"), type="logical", default=FALSE,
              help="whether to compare to glm", metavar="logical")
  
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

if(spa_opt){
  marker_test_df %>% ggplot(aes(x=SPA_pvalue)) +geom_histogram()+ggtitle(paste("SPA pvalue histogram for species ",s_id))
  marker_test_df %>% ggplot(aes(y=-log10(SPA_pvalue),x=beta))+geom_point()+ggtitle(paste("volcano plot for species ",s_id,"tau value is",glmm_fit$tau[2]))
  marker_test_df %>% ggplot(aes(y=-log10(SPA_pvalue),x=beta))+geom_point()+geom_point(y=-log10(pvalue))+ggtitle(paste("volcano plot for species ",s_id,"tau value is",glmm_fit$tau[2]))
  
}
if(opt$compare_to_glm){
  pdf(file= file.path(output_dir,paste0(s_id,".compare_to_glm.pdf")) ) 
  par( mfrow= c(4,1) )
  copy_number_df_with_y<-left_join(copy_number_df,glm_fit0$data) %>% filter(!is.na(y)) 
  
  glm_model<-copy_number_df_with_y %>% group_by(gene_id) %>% group_modify(~broom::tidy(glm(paste0(glm_fit0$formula,"+ copy_number"),data = .x,family=glm_fit0$family)))
  
  both_marker_test<-glm_model %>% filter(term=="copy_number") %>% right_join(marker_test_df)
  write.table(both_marker_test,file.path(output_dir, paste0(s_id,".both_marker_test.tsv")),sep="\t",row.names=FALSE)
  
  both_marker_test %>% ggplot(aes(y=-log10(p.value),x=estimate))+geom_point()+ggtitle(paste("glm volcano plot for species ",s_id,"tau value is",glmm_fit$tau[2]))
  
  mean_cov_glm<- copy_number_df_with_y %>% group_by(gene_id)  %>% mutate(offset=0) %>% group_map(~glm(glm_fit0$formula,data = .x,family=binomial(link = "logit"),offset=offset))
  mean_cov_df<-copy_number_df %>% group_by(gene_id) %>% summarize(num_samples=n())
  mean_cov_df$model<-as.vector(mean_cov_glm)
  mean_cov_df$Dev<-as.numeric(lapply(mean_cov_glm,function(x) x$deviance))
  mean_cov_df$ND<-as.numeric(lapply(mean_cov_glm,function(x) x$null.deviance))
  mean_cov_df$r2<-1-mean_cov_df$Dev/mean_cov_df$ND
  mean_cov_df$aic<-as.numeric(lapply(mean_cov_glm,function(x) x$aic))
  if(spa_opt){
    pvalue_choice="SPA_pvalue"
  }else{
    pvalue_choice="pvalue"
  }
  p_value_ver_3<-both_marker_test %>% ggplot(aes(x=-log10(p.value),y=-log10(paste0(pvalue_choice))))+geom_point()+ggtitle(paste("no adj pvalue for genes of species with Age:", s_id))+labs(y=c("adjusted_model"),x="glm model")+geom_abline(color="red")
  
  p_value_ver_4<-ggExtra::ggMarginal(p_value_ver_3, type = "histogram")
  print(p_value_ver_4,newpage = TRUE)
  
  beta_plot<-both_marker_test %>% ggplot(aes(x=estimate,y=beta))+geom_point()+ggtitle(paste("beta for genes of species with Age:", s_id))+labs(y=c("adjusted_model"),x="glm model")+geom_abline(color="red")
  print(beta_plot)
}


log_close()