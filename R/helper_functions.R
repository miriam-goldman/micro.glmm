color_pal<<-c("#E69F00","#CC79A7","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")

test_dir<-function(test_dir,verbose){tryCatch({
  
  if(file_test("-d",test_dir)){
    if(file_test("-x",test_dir)){
      if(verbose){
        message(paste("output directory is:",test_dir))
      }
    }else{
      simpleError("output dir doesnt exist")
      stop()
    }
    
  }else{
    simpleError("output dir doesnt exist")
    stop()
  }
  
},error=function(err) {
  err$message <- paste("output dir doesnt exist", err, sep = " ")
  # and re-raise
  stop(err)
},finally={
  output_directory=test_dir
 
  return(output_directory)})}


validate_GRM_metadata<-function(opt){
  if(isTRUE(!is.na(opt$GRM))){
    if(file_test("-f",opt$GRM)){
      GRM <- fread(opt$GRM,sep="\t",header=FALSE) 
      colnames(GRM)<-c("sample_name",GRM$V1)
      setindexv(GRM,'sample_name')
      GRM<-GRM %>% select(-sample_name)
      stopifnot(ncol(GRM)==(nrow(GRM)))
      GRM<-Matrix(as.matrix(GRM))
      dimnames(GRM)<-c(dimnames(GRM)[2],dimnames(GRM)[2])
      put("GRM read in",console = verbose)
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
      metadata <<- fread(opt$metadata) 
      metadata_overlaps<-base::intersect(metadata$sample_name,colnames(GRM))
      if(isFALSE(length(metadata_overlaps)>0)){
        put("samples in metadata do not match MIDAS output please make sure you metadata has sample names
            in a colmun labeled sample_name and binary phenotypes in a column labeled disease_status",console = verbose)
        stop()
      }else{
        metadata<<-metadata %>% filter(sample_name %in% metadata_overlaps)
        print(GRM)
        GRM<<-GRM[metadata_overlaps,metadata_overlaps]
        put("metadata read in",console = verbose)
      }
    }else{
      put("metadata invalid",console = verbose)
      stop()
    }
  }else{
    put("metadata invalid",console = verbose)
    stop()
  }
  if(is.numeric(opt$tau)){
    tau0<<-opt$tau
  }else{
    put("tau invalid",console = verbose)
    stop()
  } 
  if(is.numeric(opt$phi)){
    phi0<<-opt$phi
  }else{
    put("phi invalid",console = verbose)
    stop()
  } 
  if(is.numeric(opt$maxiter)){
    maxiter<<-opt$maxiter
  }else{
    put("maxiter invalid",console = verbose)
    stop()
  } 
  if(is.numeric(opt$tol)){
    tol<<-opt$tol
  }else{
    put("tolerance invalid",console = verbose)
    stop()
  }
  
 if(opt$family %in% c("binomial","gaussian","poisson")){
   family_to_fit<<-opt$family
   
 }else{
   put("family option invalid",console = verbose)
   stop()
 }
  if(grepl("~",opt$formula)){
    formula_to_fit<<-opt$formula
  }else{
    put("formula invalid, must be in this format y~covariates",console = verbose)
    stop()
  }
}
