suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(logr))
suppressPackageStartupMessages(library(data.table))
library("ggExtra")
suppressPackageStartupMessages(library(micro.glmm))
color_pal<<-c("#E69F00",
              "#CC79A7","#56B4E9","#009E73",
              "#F0E442","#0072B2","#D55E00",
              "#E1BE6A","#0231A9")
option_list = list(
  make_option(c("-s", "--species_id"), type="character", default="s_id", 
              help="species id for labeling output file names (optional)", metavar="character"),
  make_option(c("-v", "--verbose"), type="logical", default=TRUE, 
              help="verbose, TRUE or FALSE", metavar="logical"),
  make_option(c("-o", "--out_folder"), type="character", default=NULL,
              help="output folder name", metavar="character"),
  make_option(c("--pres_abs"), type="character", default=NULL,
              help=".tsv file with gene_id, sample_name, and copynumber for CNV data",metavar="character"),
  make_option(c("--genes_summary"), type="character", default=NULL,
              help=".tsv file with gene_id, sample_name, and copynumber for CNV data",metavar="character"),
  make_option(c("--metadata"), type="character", default=NULL,
              help=".tsv file with gene_id, sample_name, and copynumber for CNV data",metavar="character"),
  make_option(c("--SPA"), type="logical", default=FALSE,make_option(c("--pres_abs"), type="character", default=NULL,
              help=".tsv file with gene_id, sample_name, and copynumber for CNV data",metavar="character"),
              help="whether to use SPA for pvalues, will slow down output (optional)", metavar="logical"),
  make_option(c("--q_value"), type="numeric", default=.05,
              help="q value for gene test plot", metavar="numeric"),
  make_option(c("--alpha_value"), type="numeric", default=.05,
              help="alpha value for gene test plot", metavar="numeric"),
  make_option(c("--bonferroni"),type="logical",default = FALSE,
              help="whether to add bonferroni cutoff",metavar="logical")
  
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
verbose<-opt$verbose
output_dir<-test_dir(opt$out_folder,verbose)
s_id<-opt$species_id
tmp <- file.path(output_dir, paste0(s_id,".gene_test.log"))

lf <- log_open(tmp,show_notes=FALSE)
put("running validate genes command is",console = verbose)
command_arg_list<-commandArgs(trailingOnly = TRUE)
commad_message<-paste(unlist(command_arg_list),collapse = " ")

put(paste("Rscript gene_test_command_line_pres_abs.R",commad_message),console = verbose)

validate_gene_test_pres_abs(opt)
output_dir<-test_dir(opt$out_folder,verbose)






gene_test_df<-micro_glmm(glmm_fit,glm_fit0,GRM,copy_number_df,SPA=spa_opt)
pdf(file= file.path(output_dir,paste0(s_id,".gene_test_output.pdf")) ) 
# create a 2X2 grid 
par( mfrow= c(6,1) )
m<-nrow(gene_test_df)
bonferroni_cutoff<-opt$alpha_value/m

q=opt$q_value
if(spa_opt){
  gene_test_df<-gene_test_df %>% arrange(SPA_pvalue)
  gene_test_df<-gene_test_df%>% mutate(rank=seq(1,m)) %>% mutate(bh_pvalue=(rank/m)*q) %>% mutate(SPA_pvalue=ifelse(is.na(SPA_pvalue),1,SPA_pvalue))
  filtered_test_df<-gene_test_df %>% filter(SPA_pvalue<=bh_pvalue & SPA_pvalue<=.01)
  if(any(gene_test_df$SPA_pvalue<=gene_test_df$bh_pvalue & gene_test_df$SPA_pvalue<=.01)){
    bh_cutoff<-max(filtered_test_df$SPA_pvalue)
  }else{
    bh_cutoff<-bonferroni_cutoff
  }
  
}else{
  gene_test_df<-gene_test_df %>% arrange(pvalue)%>% mutate(pvalue=ifelse(is.na(pvalue),1,pvalue))
  gene_test_df<-gene_test_df%>% mutate(rank=seq(1,m)) %>% mutate(bh_pvalue=(rank/m)*q)
  filtered_test_df<-gene_test_df %>% filter(pvalue<=bh_pvalue)
  if(any(gene_test_df$pvalue<=gene_test_df$bh_pvalue)){
    bh_cutoff<-max(filtered_test_df$pvalue)
  }else{
    bh_cutoff<-bonferroni_cutoff
  }
  
}
write.table(gene_test_df,file.path(output_dir, paste0(s_id,".gene_test.tsv")),sep="\t",row.names=FALSE)


plot1<-gene_test_df %>% ggplot(aes(x=pvalue)) +geom_histogram(bins=30)+ggtitle(paste("Pvalue histogram for species ",s_id))
print(plot1)
if(opt$bonferroni==FALSE){
  bonferroni_cutoff=opt$alpha
  bh_cutoff=opt$alpha
}
plot2<-gene_test_df %>% ggplot(aes(y=-log10(pvalue),x=beta))+geom_hline(color="green",yintercept =-log10(bh_cutoff))+geom_hline(color="red",yintercept =-log10(bonferroni_cutoff))+geom_point()+ggtitle(paste("volcano plot for species ",s_id,"tau value is",glmm_fit$tau[2]))
print(plot2)
if(spa_opt){
  plot3<-gene_test_df %>% ggplot(aes(x=SPA_pvalue)) +geom_histogram(bins=30)+ggtitle(paste("SPA pvalue histogram for species ",s_id))
  print(plot3)
  plot4<-gene_test_df %>% ggplot(aes(y=-log10(SPA_pvalue),x=beta))+geom_hline(color="green",yintercept =-log10(bh_cutoff))+geom_hline(color="red",yintercept =-log10(bonferroni_cutoff))+geom_point()+ggtitle(paste("volcano plot for species with SPA",s_id,"tau value is",glmm_fit$tau[2]))
  print(plot4)
  plot5<-gene_test_df %>% ggplot(aes(x=beta))+geom_hline(color="green",yintercept =-log10(bh_cutoff))+geom_hline(color="red",yintercept =-log10(bonferroni_cutoff))+geom_point(aes(y=-log10(SPA_pvalue),color="red"))+geom_point(aes(y=-log10(pvalue),color="blue"))+ggtitle(paste("volcano plot for species ",s_id,"tau value is",glmm_fit$tau[2]))+labs(y="-log10(pvalues) blue pvalue red SPA pvalue")
  print(plot5)
}
end<-dev.off()


log_close()
