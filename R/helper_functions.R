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
