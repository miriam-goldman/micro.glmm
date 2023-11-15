###Run structue test
#after prep_snps

library("optparse")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(logr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(generics))
suppressPackageStartupMessages(library(plotROC))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(micro.glmm))
color_pal<<-c("#E69F00",
              "#CC79A7","#56B4E9","#009E73",
              "#F0E442","#0072B2","#D55E00",
              "#E1BE6A","#0231A9")
option_list = list(
  make_option(c("-s", "--species_id"), type="character", default="s_id", 
              help="species id for labeling output file names (optional)", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default="s_id", 
              help="file path to metadata", metavar="character"),
  make_option(c("-v", "--verbose"), type="logical", default=TRUE, 
              help="verbose, TRUE or FALSE", metavar="logical"),
  make_option(c("-o", "--out_folder"), type="character", default=NULL,
              help="output folder name", metavar="character"),
  make_option(c("--GRM"), type="character", default=NULL,
              help="GRM used to check sample list (optional)",metavar="character" ),
  make_option(c("--tau"), type="numeric", default=0,
              help="initail values of tau (optional)",metavar="numeric"),
  make_option(c("--phi"), type="numeric", default=1,
              help="initail values of phi (optional), only for quantitative data",metavar="numeric" ),
  make_option(c("--maxiter"), type="numeric", default=100,
              help="max iteration of fitting parameters (optional)",metavar="numeric" ),
  make_option(c("--tol"), type="numeric", default=.0001,
              help="tolerance for parameters (optional)",metavar="numeric" ),
  make_option(c("--family"),type="character",default="binomial",
              help="family to fit for the generalized linear model and generalized linear mixed model",metavar="character"),
  make_option(c("--formula_to_fit"), type="character", default="y~1",
              help="formula used to fit population structure; the phenotype and any covariates in the metadata file that are to be used",metavar="character" ),
  make_option(c("--n_tau"), type="numeric", default=0,
              help="number of permuations to test for significance of tau, if 0 no test run",metavar="numeric" ),
  make_option(c("--SPA"), type="logical", default=TRUE, 
              help="SPA, TRUE or FALSE for simulations", metavar="logical")
 
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

validate_GRM_metadata(opt)
glm_fit0 = glm(formula_to_fit, data = metadata, family = family_to_fit)
glmm_fit=pop_structure_test(glm_fit0, GRM,species_id=s_id,tau=c(phi0,tau0),maxiter =maxiter, verbose = verbose,tol=tol,log_file=TRUE)
b_df=data.frame(b=as.vector(glmm_fit$b),y=as.vector(glmm_fit$y),fitted.values=as.vector(glmm_fit$fitted.values))

# Open pdf file 
pdf(file= file.path(output_dir,paste0(s_id,".random_effect_output.pdf")) ) 

# create a 2X2 grid 
par( mfrow= c(11,1) )

annotation_col = data.frame(disease_status=paste0("group_",as.vector(glm_fit0$y)),"b"=as.vector(glmm_fit$b))
rownames(annotation_col) =colnames(GRM)
rownames(GRM)<-colnames(GRM)

ann_colors = list(
  disease_status = c("group_0" = color_pal[2], "group_1" = color_pal[9]),b=c(color_pal[5],"#FFFFFF",color_pal[3])
)

paletteLength <- 50
myColor <- colorRampPalette(c("#FFFFFF",color_pal[4]))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(.7, .85, length.out=ceiling(paletteLength/2) + 1), 
              seq(.86, 1, length.out=floor(paletteLength/2)))



pheatmap(GRM,show_rownames=FALSE,show_colnames=FALSE,labels_row="samples",labels_col="samples",
         main=paste("GRM for species",s_id),annotation_row = annotation_col,
         annotation_col=annotation_col,annotation_colors = ann_colors,color=myColor,breaks=myBreaks)


ggplot(b_df,aes(b,fill=as.factor(y)))+geom_histogram()+ scale_fill_manual(values=c(color_pal[2],color_pal[9]),name="diagnosis")+
  labs(x="values of random varible b",title=paste("Histogram of random varaible b for species",s_id,"\ncolored by diagnosis"),legend="diagnosis")

ggplot(b_df, aes(d = y, m = fitted.values)) + geom_roc()+labs(caption =paste(glmm_fit$summary))+ theme(
                                            plot.caption = element_text(hjust = 0,size=12))




save(glmm_fit, glm_fit0,GRM,s_id, file = file.path(output_dir,paste0(s_id,".model_obj.Rdata")))
if(opt$n_tau>0){
  simulate_tau<-run_tau_test(glm_fit0,GRM,opt$n_tau,s_id,tau0,phi0)
  num_more_ext=sum(simulate_tau$t>glmm_fit$t)
  ggplot(simulate_tau,aes(t))+geom_histogram()+geom_vline(xintercept=glmm_fit$t,color=color_pal[6])+
  labs(x="permuated t values",title=paste("Histogram of permuated t values for species",s_id,
                                           "\n actual t value for a tau of",glmm_fit$tau[2],"in orange",
                                           "\n number of permutations",opt$n_tau),caption =paste(glmm_fit$summary,"\n pvalue for t is", num_more_ext/opt$n_tau))+ theme(
                                             plot.caption = element_text(hjust = 0)
                                           )
  write.csv(data.frame(s_id=s_id,n_tau=opt$n_tau,tau=glmm_fit$tau[2],pvalue=num_more_ext/opt$n_tau,tvalue=glmm_fit$t),file.path(output_dir,paste0(s_id,".tau_file.csv")))
}


dev.off()
# Close log
log_close()



