library(pROC)
library(Matrix)
### functions needed to run GLMM for MWAS
getCoefficients<-function(Y, X, W, tau, GRM){
  # Y is working vector Y=alpha X + b
  # n number of samples 
  # N number of covarites
  # inputs are Y (nx1) diagnosis
  #X (nxN) covarities and diagnosis
  # W coefficient of variation 
  # tau is the variance of the residual errors
  simga=gen_sp_Sigma(W,tau,GRM)# V
  Y=as.vector(Y)
  Sigma_iY=solve(simga,Y) # V^-1 Y
  sigma_X=solve(simga,X) # V^-1 X
  cov_var=Matrix::solve(forceSymmetric(t(X) %*% sigma_X),sparse=TRUE,tol = 1e-10) # (Xt V^-1 X)^-1
  Sigma_iXt=t(sigma_X) # t(V^-1 X)
  SigmaiXtY= Sigma_iXt%*%Y # XtV^-1Y
  alpha = cov_var %*% SigmaiXtY # (Xt V X)^-1 XtVY 
  #b1= tau[2]*GRM %*%solve(simga)%*%(Y-X %*% alpha)
  #epsilon=Y-(X%*%alpha+b)
  epsilon=tau[1] * (t(Sigma_iY) - t(sigma_X %*% alpha)) / as.vector(W)#tau[1] to act on W
 
  eta = as.vector(Y - epsilon) # Y-tau \sigma (Y-X\alpha) 
  
  
  b= eta-X %*% alpha
  #eta=Y-epsilon
  #print(round(eta_2,2))
  re=list("Sigma_iY"=Sigma_iY,"Sigma_iX"=sigma_X,"cov_var"=cov_var,"alpha"=alpha,"eta"=eta,"b"=b,"epsilon"=epsilon)
}
gen_sp_Sigma<-function(W,tau,kinship){
  ### update kinship with W and tau
  ## value vector is kin
  #kinship is an (nxn) symetric matrix
  dtkin=W^-1 * (tau[1]) # inverse W 
  new_kin = kinship* tau[2]
  diag(new_kin)=diag(new_kin)+dtkin
  new_kin[new_kin<1e-4]=1e-4
  return(as.matrix(new_kin))
}

get_AI_score<-function(Y,X,GRM,W,tau,Sigma_iY,Sigma_iX,cov_var){
  ## get score for finding tau function from supplment of Saige paper
  ## Inputs Y, X, GRM, W, Tau, Sigma_Y, Sigma_X, cov_var
  Sigma=gen_sp_Sigma(W,tau,GRM)
  Sigma_iXt = t(Sigma_iX) #transpose X
  P=solve(Sigma) - Sigma_iX %*% cov_var %*% Sigma_iXt
  PY1 = P %*% Y# \hat{Y}-\hat(X) (Xt V X)^-1 PY
  APY = GRM %*% PY1 # GRM (\hat{Y}-\hat(X) (Xt V X)^-1) 
  YPAPY = t(PY1) %*% APY# dot product
  YPAPY=YPAPY[1]# dot product
  PA= P %*% GRM
  Trace_P_GRM = sum(diag(PA))
  score1=YPAPY-Trace_P_GRM
  PAPY=P%*% APY
  
  AI = (t(PAPY) %*% APY)# AI=t(Y)%*%P%*%GRM%*%P%*%GRM%*%P%*%Y
  return(list(YPAPY=YPAPY,PY=PY1,Trace_P_GRM=Trace_P_GRM,score1=score1,AI=AI[1]))
}

get_AI_score_quant<-function(Y,X,GRM,W,tau,Sigma_iY,Sigma_iX,cov_var){
  ## get score for finding tau function from supplment of Saige paper
  ## Inputs Y, X, GRM, W, Tau, Sigma_Y, Sigma_X, cov_var
  n=length(W)
  Sigma=gen_sp_Sigma(W,tau,GRM)
  Sigma_iXt = t(Sigma_iX) #transpose X
  P=solve(Sigma) - Sigma_iX %*% cov_var %*% Sigma_iXt
  diag_P=diag(P)/W
  PY1 = P %*% Y# \hat{Y}-\hat(X) (Xt V X)^-1 PY
  wPY=PY1/W
  YPwPY = t(PY1) %*% wPY
  YPwPY=YPwPY[1]
  APY=GRM %*% PY1
  YPAPY = t(PY1) %*% APY# dot product
  YPAPY=YPAPY[1]# dot product
  PA= P %*% GRM
  Trace_P_GRM =  (sum(solve(Sigma)*GRM)-sum(Sigma_iX*crossprod(GRM,t(cov_var %*% Sigma_iXt))))
  Trace_PW=sum(diag_P)
  score1=YPAPY-Trace_P_GRM#score 1
  score0=YPwPY-Trace_PW #score 0 good
  score_vector=as.matrix(c(score0[1],score1[1]))
  PwPY = P %*% wPY
  PAPY=P%*% APY
  
  
  AI_11 = (t(PAPY) %*% APY)#good
  AI_00=(t(PwPY) %*% wPY) 
  AI_01= (t(PAPY) %*% wPY)
  AI_mat=matrix(c(AI_00[1],AI_01[1],AI_01[1],AI_11[1]),2,2)
  Dtau <- solve(AI_mat, score_vector)
  
  
  
  return(list(YPAPY=YPAPY,PY=PY1,YPwPY=YPwPY,Trace_P_GRM=Trace_P_GRM,Trace_PW=Trace_PW,AI=AI_mat,score_vector=score_vector))
}




ScoreTest_NULL_Model = function(mu, y, X){
  ## score test for null model uses fitted mu, real y, and X
  mu2=mu*(1-mu)
  V = as.vector(mu2)
  res = as.vector(y - mu)
  XV = t(X * V)
  XVX = t(X) %*% (t(XV))
  XVX_inv = solve(XVX)
  XXVX_inv = X %*% XVX_inv
  XVX_inv_XV = XXVX_inv * V
  S_a =  colSums(X * res)
  re = list(XV = XV, XVX = XVX, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv, S_a = S_a, XVX_inv_XV = XVX_inv_XV, V = V)
  class(re) = "SA_NULL"
  return(re) 
}	

ScoreTest_NULL_Model_quant = function(mu,tau, y, X){
  V = rep(1/tau[1], length(y))
  res = as.vector(y - mu)
  XV = t(X * V)
  XVX = t(X) %*% (t(XV))
  XVX_inv = solve(XVX)
  XXVX_inv = X %*% XVX_inv
  XVX_inv_XV = XXVX_inv * V
  S_a =  colSums(X * res)
  re = list(XV = XV, XVX = XVX, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv, S_a = S_a, XVX_inv_XV = XVX_inv_XV, V = V)
  class(re) = "SA_NULL"
  return(re) 
}	



fitglmmaiRPCG<-function(Yvec, Xmat,GRM,wVec,  tauVec, Sigma_iY, Sigma_iX, cov_var,tol,quant=FALSE,verbose,write_log){
  if(!quant){ 
    print("is binom")
    re.AI = get_AI_score(Yvec, Xmat,GRM,wVec,  tauVec, Sigma_iY, Sigma_iX, cov_var)
    score1 = re.AI$score1# this is equation 8 from paper 
    AI1 = re.AI$AI
    Dtau = score1/AI1
    
    tau0 = tauVec
    
    tauVec[2] = tau0[2] + Dtau
    step = 1.0
    while(tauVec[2]<0){
      
      step = step*0.5
      tauVec[2] = tau0[2] + step*Dtau
      
    }
    
  }else{
    re.AI = get_AI_score_quant(Yvec, Xmat,GRM,wVec,  tauVec, Sigma_iY, Sigma_iX, cov_var)
    YPAPY = re.AI$YPAPY
    YPwPY = re.AI$YPwPY
    Trace_PW=re.AI$Trace_PW
    Trace_P_GRM=re.AI$Trace_P_GRM
    Dtau = solve(re.AI$AI, re.AI$score_vector)
    tau0 = tauVec
    tauVec = tau0 + Dtau
    step = 1.0
    while(any(tauVec<0)){
      
      step = step*0.5
      tauVec = tau0 + step * Dtau
      
    }
  }
  
  
  
  
  if(any(tauVec < tol)){
    tauVec[which(tauVec < tol)]=0
  }
  
  return(list("tau" = tauVec))
}

#' pop_structure_test
#' 
#' Fit the base model for SNP structure
#' 
#' @param glm_fit0 glm model. Model output with no sample relatedness accounted for
#' @param GRM Genetic Relatedness Matrix (from scripts or user) NxN matrix of sample relatedness
#' @param tau vector for initial values for the variance component parameter estimates 
#' @param maxiter maximum iterations to fit the glmm model
#' @param verbose whether outputting messages in the process of model fitting
#' @param log_file log file to write to
#' @return model output for the baseline structure glmm 
#' @export
pop_structure_test = function(glm_fit0, GRM,species_id,tau=c(1,1),maxiter =100, verbose = TRUE,tol=.0001,log_file=NA) {
  #Fits the null generalized linear mixed model for a binary trait
  #Args:
  #  glm_fit0: glm model. Logistic model output (with no sample relatedness accounted for) 
  #  GRM from SNP data
  #  tau: vector for iniial values for the variance component parameter estimates
  #  maxiter: maximum iterations to fit the glmm model
  #  verbose: whether outputting messages in the process of model fitting
  #Returns:
  #  model output for the null glmm
  write_log=!is.na(log_file)
  t_begin = proc.time()
  if(verbose){
    cat("begining time ")
    cat(t_begin)	
    
  }
  if(write_log){
    put("begining time ")
    put(t_begin)
  }
  y = glm_fit0$y
  n = length(y)
  X = model.matrix(glm_fit0)
  Xorig= model.matrix(glm_fit0)
  offset = glm_fit0$offset
  if(is.null(offset)){
    offset = rep(0, n)
  }
  
  family = glm_fit0$family
  eta = glm_fit0$linear.predictors
  mu = glm_fit0$fitted.values
  mu.eta = family$mu.eta(eta)
  Y = eta - offset + (y - mu)/mu.eta
  alpha0 = glm_fit0$coef
  eta0 = eta
  sample_ids<-colnames(GRM)
  if(family$family %in% c("poisson", "binomial")) {
    tau[1] = 1
    quant=FALSE
  }else{
    quant=TRUE
  }

  if(verbose) cat(" Fixed-effect coefficients: ", glm_fit0$coef,"\n")
  if(write_log) put(paste(" Fixed-effect coefficients: ", glm_fit0$coef))
  if(verbose) cat(" inital tau is ", tau,"\n")
  if(write_log) put(paste(" inital tau is ", tau))
  tau0=tau
  if(tau[1]<=0){
      stop("ERROR! The first variance component parameter estimate is 0\n")
    }
  sqrtW_0 = mu.eta/sqrt(family$variance(mu))
  W_0 = sqrtW_0^2
  re.coef = Get_Coef(y, X, tau, GRM,family, alpha0, eta0,  offset, maxiter=maxiter,verbose=verbose,tol.coef = tol,write_log=write_log)
  
  ######
  if(quant){
    re = get_AI_score_quant(re.coef$Y, X, GRM,re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov_var)
    tau[2] = max(0, as.numeric(tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace_P_GRM)/n))
    tau[1] = max(0, as.numeric(tau0[1] + tau0[1]^2 * (re$YPwPY - re$Trace_PW)/n))
  }else{
    re = get_AI_score(re.coef$Y, X, GRM,re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov_var)
    
    tau[2] = max(0, as.numeric(tau0[2] + tau0[2]^2 * ((re$YPAPY - re$Trace_P_GRM))/n)) #tau + Dtau dumb way
  }

  for(i in seq_len(maxiter)){
    if(verbose) {
      cat(paste("\ni",i))
    }
    if(verbose) cat("\nIteration ", i, "tau is ", tau, "\n")
    if(write_log) put(paste(" Iteration ", i, "tau is: ", tau))
    alpha0 = re.coef$alpha
    tau0 = tau
    #cat("tau0_v1: ", tau0, "\n")
    eta0 = eta
    rss_0=sum((y-mu)^2)
    # use Get_Coef before getAIScore       
    t_begin_Get_Coef = proc.time()
    re.coef = Get_Coef(y, X, tau, GRM,family, alpha0, eta0,  offset,verbose=verbose,maxiter=maxiter,tol.coef = tol,write_log=write_log)
    t_end_Get_Coef =  proc.time()
    if(verbose) {
    cat("\nt_end_Get_Coef - t_begin_Get_Coef\n")
    cat(t_end_Get_Coef - t_begin_Get_Coef)
    }
    if(write_log) {
      put("t_end_Get_Coef - t_begin_Get_Coef")
      put(t_end_Get_Coef - t_begin_Get_Coef)
    }
    ##update tau
   
    fit = fitglmmaiRPCG(re.coef$Y, X, GRM, re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov_var,tol=tol,verbose=verbose,write_log=write_log,quant=quant)
    t_end_fitglmmaiRPCG= proc.time()
    if(verbose) {
    cat("\nt_end_fitglmmaiRPCG - t_begin_fitglmmaiRPCG\n")
    cat(t_end_fitglmmaiRPCG - t_end_Get_Coef)
    }
    if(write_log) {
      put("t_end_fitglmmaiRPCG - t_begin_fitglmmaiRPCG")
      put(t_end_fitglmmaiRPCG - t_end_Get_Coef)
    }
    tau = as.numeric(fit$tau)
    cov_var = re.coef$cov_var
    alpha = re.coef$alpha
    eta = re.coef$eta
    Y = re.coef$Y
    mu = re.coef$mu
    res = y - mu
     if(verbose) {
       cat(paste("\nchange in tau",abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)))
    cat("\ntau: ", tau, "\n")
    cat("\ntau0: ", tau0, "\n")
     }
    if(write_log) {
      put(paste("change in tau",abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)))
      put(paste("tau: ", tau))
      put(paste("tau0: ", tau0))
    }
    if(tau[1]<=0){
      stop("\nERROR! The first variance component parameter estimate is 0\n")
    }
    
    
    if(tau[2] == 0) break
    # Use only tau for convergence evaluation, because alpha was evaluated already in Get_Coef
    tau_condition=max(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol
    
 
    rss=sum(res^2)
    rss_condition=rss_0-rss
    if(verbose){
      cat(paste("\nres",rss))
      cat(paste("\nrss change",rss_condition))
    }
    if(write_log){
      put(paste("res",rss))
      put(paste("rss change",rss_condition))
    }
    
    abs_condition=sum(res^2)
    if(tau_condition) break
    
    if(max(tau) > tol^(-2)) {
      i = maxiter
      break
    }
  }
  if(verbose) cat("\niter break at ",i)
  if(verbose) cat("\nFinal " ,tau, ":\n")
  if(write_log) put(paste("iter break at ",i))
  if(write_log) put(paste("Final " ,tau, ":"))
  if(max(tau) > tol^(-2)){
    cat("Model not converged")
    if(write_log){
      put("Model not converged")
    }
    return(glm_fit0)
  }
  
  re.coef = Get_Coef(y, X, tau, GRM,family, alpha, eta,  offset,verbose=verbose, maxiter=maxiter,tol.coef = tol,write_log=write_log)
  if(quant){
    re.final = get_AI_score_quant(re.coef$Y, X, GRM,re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov_var)
    tau[2] = max(0, tau0[2] + tau0[2]^2 * (re.final$YPAPY - re.final$Trace_P_GRM)/n)
    tau[1] = max(0, tau0[1] + tau0[1]^2 * (re.final$YPwPY - re.final$Trace_PW)/n)
  }else{
    re.final = get_AI_score(re.coef$Y, X, GRM,re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov_var)
    
    tau[2] = max(0, as.numeric(tau0[2] + tau0[2]^2 * ((re.final$YPAPY - re.final$Trace_P_GRM))/n)) #tau + Dtau dumb way
  }
  cov_var = re.coef$cov_var
  
  alpha = re.coef$alpha
  eta = re.coef$eta
  Y = re.coef$Y
  mu = re.coef$mu
  mu.eta = family$mu.eta(eta)
  sqrtW = mu.eta/sqrt(family$variance(mu))
  W=sqrtW^2
  Sigma=gen_sp_Sigma(W,tau,GRM)
  converged = ifelse(i < maxiter, TRUE, FALSE)
  res = y - mu
  rss=sum(res^2)
  if(quant){
    mu2 = rep(1/tau[1],length(y))
    obj.noK = ScoreTest_NULL_Model_quant(mu,tau,y,Xorig)
  }else{
    mu2 = mu * (1-mu)
    obj.noK = ScoreTest_NULL_Model(mu, y, Xorig)
  }
  
  #exp_tau=GetTrace_2(X,GRM,W,tau)
  #print(paste("exp tau", exp_tau))
  #print(paste("tau_test",pchisq((1+tau[2])/(1+exp_tau/2),1,lower.tail = FALSE)))
  ss=sum((y-mean(y))^2)
  #print(rss/ss)
  var_fixed=var(Xorig%*%alpha)
  var_random=var(as.vector(re.coef$b))
  var_error=var(res)
  #https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210x.2012.00261.x
  model_metrics=list("S"=sd(res),"R-sq"=(1-rss/ss),"R-sq(marginal)"=var_fixed/(var_fixed+var_random+var_error),"r-sq(conditional)"=(var_fixed+var_random)/(var_fixed+var_random+var_error))
  call=paste("pop_structure_test formula=",glm_fit0$formula,"family=",glm_fit0$family)
  Coefficients=paste("Coefficents:",alpha)
  ###make another option for qaunt
  if(!quant){
    AUC=paste("AUC",suppressMessages(suppressWarnings(auc(glm_fit0$y,as.vector(mu)))))
    tau_script=paste("tau is:",tau[2],"t value is",sum(re.coef$b^2))
    summary_text=paste(call,Coefficients,AUC,tau_script,sep="\n")
  }else{
    tau_script=paste("tau is:",tau[2],"t value is",sum(re.coef$b^2))
    summary_text=paste(call,Coefficients,sep="\n")
  }
  
  glmmSNPResult = list(summary_text=summary_text,tau=tau,
                    coefficients=alpha, b=re.coef$b,t=sum(re.coef$b^2),
                    linear.predictors=eta, linear_model=Xorig%*%alpha+re.coef$b, 
                    fitted.values=mu, var_mu=mu2,Y=Y, residuals=res, 
                    cov_var=cov_var, converged=converged,
                    sampleID = sample_ids, 
                    obj.noK=obj.noK, 
                    y = y, X = Xorig, 
                    traitType=glm_fit0$family,
                    iter_finised=i,
                    model_metrics=model_metrics,species_id=species_id)
  
  t_end_null = proc.time()
  if(verbose) {
  cat("\nt_end_null - t_begin,  fitting the structure model took\n")
  cat(t_end_null - t_begin)
  }
  
  if(write_log) {
    put("t_end_null - t_begin, fitting the structure model took")
    put(t_end_null - t_begin)
  }
  
  return(glmmSNPResult)
}


Get_Coef = function(y, X, tau, GRM,family, alpha0, eta0,  offset, verbose=FALSE,maxiter,tol.coef=tol,write_log=FALSE){
  mu = family$linkinv(eta0)
  mu.eta = family$mu.eta(eta0)
  Y = eta0 - offset + (y - mu)/mu.eta
  
  sqrtW = mu.eta/sqrt(family$variance(mu))
  
  W = sqrtW^2
  for(i in 1:maxiter){
    re.coef = getCoefficients(Y, X, W, tau, GRM)
    
    alpha = as.matrix(re.coef$alpha)
    eta = as.matrix(re.coef$eta + offset)
    
    if(verbose) {
      cat("\n Tau:\n")
      cat(tau)
      cat("\n Fixed-effect coefficients:\n")
      cat(alpha)
    }
    if(write_log) {
      put(paste(" Tau:",tau," Fixed-effect coefficients:",alpha))
    }
    mu = family$linkinv(eta)
    mu.eta = family$mu.eta(eta)
    
    Y = eta - offset + (y - mu)/mu.eta
    sqrtW = mu.eta/sqrt(family$variance(mu))
    W = sqrtW^2
    
    if( max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol.coef))< tol.coef){
      break
    }
    alpha0 = alpha
  }
  
  re = list(Y=Y, alpha=alpha, eta=eta, W=W, cov_var=re.coef$cov_var, sqrtW=sqrtW, Sigma_iY = re.coef$Sigma_iY, Sigma_iX = re.coef$Sigma_iX, mu=mu,eta_2=re.coef$eta_2,b=re.coef$b)
}

#' micro_glmm
#' 
#' Fit the CNV model with the random effects
#' 
#' @param obj.pop.strut output of pop_structure_test; GLMM of species with GRM accounted for
#' @param glm_fit0 glm model. Model output with no sample relatedness accounted for
#' @param GRM Genetic Relatedness Matrix (from scripts or user) NxN matrix of sample relatedness
#' @param copy_number_df data frame with, gene_id, sample_name, and copy_number 
#' @param SPA whether to run Saddle point approximation for pvalues (will slow down output) 
#' @param scale_g whether to scale CNVs after log transform if using
#' @param log_g whether to log transform CNVs 
#' @return model output for the baseline structure glmm 
#' @export
micro_glmm = function(obj.pop.strut,
                    glm_fit0,GRM,
                           copy_number_df,SPA=FALSE,scale_g=TRUE,log_g=TRUE){
  ## inputs: full null fitted model
  ## null model no random effect
  # GRM is kinship
  # copy_number_df is the copy number by gene 
  ## copy_number_df must have column for sample; column for gene_id; column for copy_number
  ## returns list of values for each gene examined
  list_vec<-NULL
  t_begin = proc.time()

  
  obj.noK = obj.pop.strut$obj.noK
  family = glm_fit0$family
  
  eta = obj.pop.strut$linear.predictors
  mu = obj.pop.strut$fitted.values
  mu.eta = family$mu.eta(eta)
  sqrtW = mu.eta/sqrt(glm_fit0$family$variance(mu))
  W1 = sqrtW^2   ##(mu*(1-mu) for binary) theses are the same
  tauVecNew = obj.pop.strut$tau
  X = obj.pop.strut$X
  Sigma=gen_sp_Sigma(W1,tauVecNew,GRM)
  obj.pop.strut$Sigma<-Sigma
  Sigma_iX<-solve(Sigma,X)
  obj.pop.strut$Sigma_iX<-Sigma_iX
  Sigma_iY<-solve(Sigma,obj.pop.strut$Y)
  obj.pop.strut$Sigma_iY<-Sigma_iY
  sample_lookup<-data.frame(sampleID=obj.pop.strut$sampleID,index=seq(1,length(obj.pop.strut$sampleID)))
  sample_genes<-unique(copy_number_df$gene_id)
  
  ##randomize the marker orders to be tested
  
  for(k in sample_genes){
    iter=which(sample_genes==k)
    if(iter %% 1000 == 0){
      cat(paste("number of genes done ",iter,"\n"))
      cat("time past:")
      t_now = proc.time()
      cat(t_now-t_begin)	
      cat("\n")
    }
    one_gene<-copy_number_df %>% ungroup %>% filter(gene_id==k)
    #one_gene_indexs<-sample_lookup %>% inner_join(one_gene,by=c("sampleID"="sample_name")) %>% select(sampleID,index)
    
    one_gene<-one_gene %>% inner_join(sample_lookup,by=c("sample_name"="sampleID"))
    ## filter obj to samples present in gene copy number
    filtered_obj.pop.strut<-filter_null_obj(obj.pop.strut,one_gene)
    empty_mat<-matrix(0,nrow(one_gene),nrow(one_gene))
    # log and scale copy_number
    if(log_g==TRUE){
      G0<-log(as.vector(one_gene$copy_number))
    }else{
      G0<-as.vector(one_gene$copy_number)
    }
    if(scale_g==TRUE){
      G0<-scale(G0)
    }
   
    G_tilde = G0  -  filtered_obj.pop.strut$obj.noK$XXVX_inv %*%  (filtered_obj.pop.strut$obj.noK$XV %*% G0) # G1 is X adjusted
    res=filtered_obj.pop.strut$residuals
    eta = filtered_obj.pop.strut$linear.predictors
    
    mu = filtered_obj.pop.strut$fitted.values
    mu.eta = family$mu.eta(eta)
    sqrtW = mu.eta/sqrt(glm_fit0$family$variance(mu))
    W=sqrtW^2#mu*(1-mu)
    Sigma_iG = solve(filtered_obj.pop.strut$Sigma,G_tilde)
    PG_tilde<-Sigma_iG-filtered_obj.pop.strut$Sigma_iX%*%(solve(t(filtered_obj.pop.strut$X)%*%filtered_obj.pop.strut$Sigma_iX))%*%t(filtered_obj.pop.strut$X)%*%Sigma_iG
    Y = eta  + (filtered_obj.pop.strut$y - mu)/mu.eta
    
    
    t_score=t(PG_tilde)%*%(Y)/tauVecNew[1] #t_score_2=t(G_tilde)%*%(filtered_obj.pop.strut$y - mu)
   
    m1 = t(G_tilde) %*% mu
    #qtilde=t(G_tilde)%*%filtered_obj.pop.strut$y
    var1 = t(G_tilde)%*%PG_tilde## same as  t(G)%*%Sigma_iG - t(G)%*%filtered_obj.pop.strut$Sigma_iX%*%(solve(t(filtered_obj.pop.strut$X)%*%filtered_obj.pop.strut$Sigma_iX))%*%t(filtered_obj.pop.strut$X)%*%Sigma_iG
    t_adj=t_score/sqrt(var1)
    t_adj_2=t_score^2/var1
    beta=t_score/var1
    pval=(pchisq(t_adj_2,lower.tail = FALSE,df=1,log.p=FALSE))
    z=(qnorm(pval/2, log.p=F, lower.tail = F))
    se_beta=abs(beta)/sqrt(abs(z))
    new_eta=beta[1,1]*G0+filtered_obj.pop.strut$b+filtered_obj.pop.strut$X %*% filtered_obj.pop.strut$coe
    new_mu=family$linkinv(new_eta)
    qtilde=t_score+m1
    if(SPA){
      if(var1<0){
        list_vec=rbind(list_vec,data.frame("species_id"=obj.pop.strut$species_id,tau=obj.pop.strut$tau[2],"gene_id"=k,"cor"=cor(G0,filtered_obj.pop.strut$y),"z"=z,"var1"=var1,"beta"=NA,"se beta"=NA,"pvalue"=NA,
                                           "t_adj"=NA,"num_control"=sum(filtered_obj.pop.strut$y==0),
                                           "num_total"=length(G0),
                                           SPA_pvalue=NA,spa_score=NA,pvalue_noadj=NA))
      }else{
        out1 = Saddle_Prob(q=qtilde, mu = mu, g = G_tilde, var1,Cutoff = 2,log.p=FALSE)
        list_vec=rbind(list_vec,data.frame("species_id"=obj.pop.strut$species_id,tau=obj.pop.strut$tau[2],"gene_id"=k,"cor"=cor(G0,filtered_obj.pop.strut$y),"z"=z,"var1"=var1,"beta"=beta,"se beta"=se_beta,"pvalue"=pval,
                                           "t_adj"=t_adj,"num_control"=sum(filtered_obj.pop.strut$y==0),
                                           "num_total"=length(G0),
                                           SPA_pvalue=out1$p.value,spa_score=out1$Score,pvalue_noadj=out1$p.value.NA))
      }
   
    }else{
      list_vec=rbind(list_vec,data.frame("species_id"=obj.pop.strut$species_id,tau=obj.pop.strut$tau[2],"gene_id"=k,"cor"=cor(G0,filtered_obj.pop.strut$y),"z"=z,"var1"=var1,"beta"=beta,"se beta"=se_beta,"pvalue"=pval,
                                         "t_adj"=t_adj,"num_control"=sum(filtered_obj.pop.strut$y==0),
                                         "num_total"=length(G0)))
      
    }
  }
  cat("total time past:")
  t_end = proc.time()
  cat(t_end-t_begin)	
  return(list_vec)
}

filter_null_obj<-function(obj.pop.strut,sample_indexs){
  obj.pop.strut$residuals<-obj.pop.strut$residuals[sample_indexs$index]
  obj.pop.strut$b<-obj.pop.strut$b[sample_indexs$index]
  obj.pop.strut$linear.predictors<-obj.pop.strut$linear.predictors[sample_indexs$index]
  obj.pop.strut$fitted.values<-obj.pop.strut$fitted.values[sample_indexs$index]
  obj.pop.strut$obj.noK$XXVX_inv<-obj.pop.strut$obj.noK$XXVX_inv[sample_indexs$index,]
  obj.pop.strut$obj.noK$XV<-obj.pop.strut$obj.noK$XV[,sample_indexs$index]
  obj.pop.strut$Sigma<- obj.pop.strut$Sigma[sample_indexs$index,sample_indexs$index]
  obj.pop.strut$X<-obj.pop.strut$X[sample_indexs$index,]
  obj.pop.strut$y<-obj.pop.strut$y[sample_indexs$index]
  obj.pop.strut$Sigma_iX<-obj.pop.strut$Sigma_iX[sample_indexs$index,]
  obj.pop.strut$Sigma_iY<-obj.pop.strut$Sigma_iY[sample_indexs$index]
  return(obj.pop.strut)
}

simulate_one_CNV<-function(sample_names,y,mean_sim,spread_sim,beta=1){
  nsamples=length(sample_names)
  y_delta=if_else(y==0,-1,1)
  new_y=beta*y_delta+rnorm(nsamples,mean_sim,spread_sim)
  return(new_y)
}

simulate_type1_error<-function(obj.pop.strut,glm_fit0,GRM,n_CNV=5000,alpha_value=.05,mean_sim=0,spread_sim=1,plot_qq=TRUE,SPA=FALSE){
  n_samples<-nrow(GRM)
  gene_ids=unlist(lapply(paste0("gene",seq(1,n_CNV)),function(x) rep(x,n_samples)))
  sample_names=rep(obj.pop.strut$sampleID,n_CNV)
  copynumbers=unlist(lapply(seq(1,n_CNV),function(x) simulate_one_CNV(obj.pop.strut$sampleID,obj.pop.strut$y,mean_sim,spread_sim,beta=0)))
  fake_copy_number_data<-data.frame(gene_id=gene_ids,
                                    sample_name=sample_names,
                                    copy_number=copynumbers)
  fake_data<-micro_glmm(obj.pop.strut,glm_fit0,GRM,fake_copy_number_data,SPA=SPA,scale_g=FALSE,log_g=FALSE)
  if(SPA){
    pvalues=fake_data$SPA_pvalue
    if(plot_qq){
      
      simpleQQPlot(pvalues,obj.pop.strut$tau,alpha_value,n_CNV,obj.pop.strut,SPA)
    }
  }
  pvalues=fake_data$pvalue
  if(plot_qq){
    
    simpleQQPlot(pvalues,obj.pop.strut$tau,alpha_value,n_CNV,obj.pop.strut,SPA)
  }
  
  return(fake_data)
  #return(fake_data)
}

simulate_power<-function(obj.pop.strut,glm_fit0,GRM,n_CNV=5000,alpha_value=.05,mean_sim=0,spread_sim=1,SPA=FALSE){
  n_samples<-nrow(GRM)
  beta_list=c(.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.75,1)
  gene_ids=unlist(lapply(paste0("gene",seq(1,n_CNV)),function(x) rep(x,n_samples)))
  sample_names=rep(obj.pop.strut$sampleID,n_CNV)
  beta_df=data.frame(beta=beta_list,num_recovered=0)
  for(beta in beta_list){
    prec_real<-round(n_CNV*.1,1)
    prec_fake<-round(n_CNV*.9,1)

    betas_to_sim<-c(rep(0,prec_fake),rep(beta,prec_real))
    copynumbers=unlist(lapply(betas_to_sim,function(x) simulate_one_CNV(obj.pop.strut$sampleID,obj.pop.strut$y,mean_sim,spread_sim,beta=x)))
    fake_copy_number_data<-data.frame(gene_id=gene_ids,
                                      sample_name=sample_names,
                                      copy_number=copynumbers)  
    fake_data<-micro_glmm(obj.pop.strut,glm_fit0,GRM,fake_copy_number_data,SPA=SPA,scale_g=TRUE,log_g=FALSE)
    
    fake_data$sim_beta<-betas_to_sim
    
    fake_data<-fake_data %>% mutate(num_found=(fake_data$pvalue <= alpha_value & sim_beta>0)) %>% arrange(-sim_beta)
    print(head(fake_data))
    beta_df[which(beta_df$beta==beta),2]=sum(fake_data$num_found)/sum(fake_data$sim_beta>0)
    if(SPA){
      pvalues=fake_data$SPA_pvalue
    }else{
      pvalues=fake_data$pvalue
    }
  }
 
  
  print(ggplot(beta_df,aes(beta_list,num_recovered))+geom_point(size=3)+geom_smooth(formula= y ~ log(x),se=FALSE,size=2)+geom_hline(yintercept = .9)+labs(title=paste("Power test for species",obj.pop.strut$species_id,"at alpha value",alpha_value,"number CNV",n_CNV),x="simulated beta",y="precentage recovered"))
  beta_df$case_contol<-fake_data$num_control[1]/fake_data$num_total[1]
  beta_df$num_total<-fake_data$num_total[1]
  beta_df$species_id<-fake_data$species_id[1]
  return(beta_df)
}



simulate_tau_inner<-function(glm_fit0,GRM,species_id=s_id,tau0,phi0){
  family_to_fit<-glm_fit0$family
  data.new<-glm_fit0$data
  formulate_to_fit<-glm_fit0$formula
  data.new_shuffled<-data.new[sample(1:nrow(data.new),nrow(data.new)),]

  fit_logistic = glm(formulate_to_fit, data = data.new_shuffled, family = family_to_fit)
  fit_glmm_snp<-tryCatch(pop_structure_test(fit_logistic,GRM,tau=c(phi0,tau0),verbose = FALSE,species_id=s_id),error=function(e) NULL)
  if(isTRUE(!is.na(fit_glmm_snp))){
    t=sum(fit_glmm_snp$b^2,na.rm=TRUE)
    tau=fit_glmm_snp$tau[2]
  }else{
    tau=0
    t=0
  }
  return(data.frame("tau"=tau,t))
}

simulate_tau_outer<-function(glm_fit0,GRM,n_tau,species_id=s_id,tau0,phi0){
  list_of_tau<-lapply(seq(1,n_tau),function(x) simulate_tau_inner(glm_fit0,GRM,species_id=s_id,tau0,phi0))
  df_of_tau<-do.call(rbind, list_of_tau)
  return(df_of_tau)
}

LRT_score<-function(obj.pop.strut,glm_fit0,GRM){
  obj.noK = obj.pop.strut$obj.noK
  family = glm_fit0$family
  
  eta = obj.pop.strut$linear.predictors
  mu = obj.pop.strut$fitted.values
  mu_mean=mean(obj.pop.strut$fitted.values)
  y=obj.pop.strut$y
  mu.eta = mean(family$mu.eta(eta))
  sqrtW = mu.eta/sqrt(glm_fit0$family$variance(mu))
  W = (mu*(1-mu))   ##(mu*(1-mu) for binary) theses are the same
  N=length(W)
  W_mat_a<-diag(N)
  diag(W_mat_a)=1/(mu*(1-mu))
  tauVecNew = obj.pop.strut$tau
  X = obj.pop.strut$X
  #Sigma=forceSymmetric(gen_sp_Sigma(W,tauVecNew,GRM))
  Sigma=(tauVecNew[2]*GRM+W_mat_a)
  det_sigma<-det(Sigma)
  alpha=obj.pop.strut$coefficients
  Y_hat<-eta + (y-mu)/W#
  new_cov=t(X)%*%solve(Sigma)%*%X
  Sigma_mat<-diag(N)
  diag(Sigma_mat)=diag(Sigma)
  
  alpha_hat<-solve(new_cov)%*%(t(X)%*%solve(Sigma)%*%Y_hat)
  Y_hat_X_alpha<-Y_hat-X%*%alpha_hat
  Y_hat_X_alpha_cov<-t(Y_hat_X_alpha)%*%solve(Sigma)%*%Y_hat_X_alpha
  if(tauVecNew[2]>.01){
    log_det_tau=log(det(tauVecNew[2]*GRM+W_mat_a))
    log_det_W=log(det(solve(W_mat_a)))
  }else{
    log_det_tau=log(det(tauVecNew[2]*GRM+W_mat_a))
    log_det_W=0
  }
  liklyhood_tau_a=-(1/2)*(log_det_tau+Y_hat_X_alpha_cov+log(det(new_cov)))
  mu_0 = glm_fit0$fitted.values
  mu_0_mean=mean(mu_0)
  eta_2 = glm_fit0$linear.predictors
  mu.eta_2 = family$mu.eta(eta_2)
  sqrtW_0 = mu.eta_2/sqrt(glm_fit0$family$variance(mu_0))
  tau_0_W=(mu_0*(1-mu_0))
  W_mat<-diag(N)
  diag(W_mat)=1/(mu_0*(1-mu_0))
  det_W<-det(W_mat)
  res=glm_fit0$y-glm_fit0$fitted.values
  y_hat_0=glm_fit0$linear.predictors+res/tau_0_W#1/tau_0_W#(y-glm_fit0$fitted.values)+#+res#

  cov_0=t(X)%*%solve(W_mat)%*%X
  alpha_hat_0=solve(cov_0)%*%(t(X)%*%solve(W_mat)%*%y_hat_0)
  Y_hat_X_alpha_0<-y_hat_0-X%*%alpha_hat_0
  cov_hat_0=t(Y_hat_X_alpha_0)%*%solve(W_mat)%*%Y_hat_X_alpha_0
  liklyhood_tau_0=-(1/2)*(log(det(W_mat))+cov_hat_0+log(det(cov_0)))

  T=2*(liklyhood_tau_a-liklyhood_tau_0)
  t_2=mean(obj.pop.strut$b^2)
  return(list("T"=T[1],"liklyhood_tau_a"=liklyhood_tau_a[1],"liklyhood_tau_0"=liklyhood_tau_0[1],
              "det_w"=log(det(W_mat)),det_w_a=log(det(W_mat_a)),det_sigma=log_det_tau,
              cov_det_0=log(det(cov_0)),cov_det_a=log(det(new_cov)),
              YPY_0=cov_hat_0[1],YpY=Y_hat_X_alpha_cov[1],t_2=t_2,tau=tauVecNew[2]
              #"cov_hat_0"=cov_hat_0[1],cov_hat=Y_hat_X_alpha_cov[1],
              #,)
         ))
}


simpleQQPlot = function (observedPValues,tau,alpha_value,n_CNV,obj.pop.strut,SPA) {
  error_rate=sum(observedPValues<alpha_value)/n_CNV
  expeded_pvalues=-log10(1:length(observedPValues)/length(observedPValues))
  observed_pvalues_tranformed=-log10(sort(observedPValues))
  pvalue_df<-data.frame(expeded_pvalues,observed_pvalues_tranformed)
  print(ggplot(pvalue_df,aes(expeded_pvalues,observed_pvalues_tranformed))+geom_point(size=3)+geom_abline(color="red",size=3)+labs(title=paste("qqplot for species,",obj.pop.strut$species_id, "tau:",round(tau[2],2),"error_rate < ",alpha_value,":",error_rate,"SPA:",SPA),x="-log10(expected P values)",y="-log10(observed p values)"))
}

#### taken from SPA

Saddle_Prob<-function(q, mu, g, var1,Cutoff=2,output="P",log.p=FALSE)
{
  m1<-sum(mu * g)
  p1=NULL
  p2=NULL
  
  Score<- q-m1
  #
  qinv = -sign(q-m1) * abs(q-m1) + m1
  t_adj_sq=(q - m1)^2/var1
  # Noadj
  pval.noadj<-pchisq(t_adj_sq, lower.tail = FALSE, df=1,log.p=FALSE)
  Is.converge=TRUE
  
  if(Cutoff=="BE"){
    rho<-sum(((abs(g))^3)*mu*(1-mu)*(mu^2+(1-mu)^2))
    B<-0.56*rho*var1^(-3/2)
    p<-B+alpha/2
    Cutoff=ifelse(p>=0.496,0.01,qnorm(p,lower.tail=F))
  } else if(Cutoff < 10^-1){
    Cutoff=10^-1
  } 			
  #
  if(output=="metaspline")
  {
    splfun<-Get_Saddle_Spline(mu,g,nodes)
  }
  if(is.na(abs(q - m1)) || is.na(var1)){
    pval=pval.noadj	
  }
  if(isTRUE(abs(q - m1)/sqrt(var1) < Cutoff )){
    
    pval=pval.noadj	
  } else {
    out.uni1<-getroot_K1(0, mu=mu, g=g, q=q)
    out.uni2<-getroot_K1(0, mu=mu, g=g, q=qinv)
    if(out.uni1$Is.converge==TRUE && out.uni2$Is.converge==TRUE)
    {
      p1<-tryCatch(Get_Saddle_Prob(out.uni1$root, mu, g, q,log.p=log.p),error=function(e) {
        if(log.p) return(pval.noadj-log(2))
        else return(pval.noadj/2)})
      p2<-tryCatch(Get_Saddle_Prob(out.uni2$root, mu, g, qinv,log.p=log.p),error=function(e) {
        if(log.p) return(pval.noadj-log(2))
        else return(pval.noadj/2)})
      if(log.p)
      {
        pval = add_logp(p1,p2)
      } else {

        pval = abs(p1)+abs(p2)
      }
      Is.converge=TRUE
    } else {
      cat("Error_Converge")
      pval<-pval.noadj
      Is.converge=FALSE	
    }				
  }
  
  if(pval!=0 && pval.noadj/pval>10^3)
  {
    return(Saddle_Prob(q, mu, g,var1, Cutoff=Cutoff*2,output,log.p=log.p))
  } else if(output=="metaspline")
  {
    return(list(p.value=pval, p.value.NA=pval.noadj, Is.converge=Is.converge,Score=Score,splfun=splfun,var=var1))
  } else {
    return(list(p.value=pval, p.value.NA=pval.noadj, Is.converge=Is.converge, Score=Score))
  }
}
getroot_K1<-function(init,mu,g,q,m1,tol=.Machine$double.eps^0.25,maxiter=1000)
{
  g.pos<-sum(g[which(g>0)])
  g.neg<- sum(g[which(g<0)])
  if(q>=g.pos || q<=g.neg)
  {
    return(list(root=Inf,n.iter=0,Is.converge=TRUE))
  } else {
    t<-init
    K1_eval<-K1_adj(t,mu,g,q)
    prevJump<- Inf
    rep<-1
    repeat
    {
      K2_eval<-K2(t,mu,g)
      tnew<-t-K1_eval/K2_eval
      if(is.na(tnew))
      {
        conv=FALSE
        break
      }
      if(abs(tnew-t)<tol)
      {
        conv<-TRUE
        break
      }
      if(rep==maxiter)
      {
        conv<-FALSE
        break
      }
      
      newK1<-K1_adj(tnew,mu,g,q)
      if(sign(K1_eval)!=sign(newK1))
      {
        if(abs(tnew-t)>prevJump-tol)
        {
          tnew<-t+sign(newK1-K1_eval)*prevJump/2
          newK1<-K1_adj(tnew,mu,g,q)
          prevJump<-prevJump/2
        } else {
          prevJump<-abs(tnew-t)
        }
      }
      
      rep<-rep+1
      t<-tnew
      K1_eval<-newK1
    }
    return(list(root=t,n.iter=rep,Is.converge=conv))
  }
}

Korg<-function(t, mu, g){
  n.t<-length(t)
  out<-rep(0,n.t)
  
  for(i in 1:n.t){
    t1<-t[i]
    temp<-log(1 - mu + mu * exp(g* t1))
    out[i]<-sum(temp)
  }
  return(out)
}

Get_Saddle_Prob<-function(zeta, mu, g, q,log.p=FALSE) 
{
  k1<-Korg(zeta, mu, g)
  k2<-K2(zeta, mu, g)
  
  if(is.finite(k1) && is.finite(k2))
  {
    temp1<-zeta * q - k1
    
    
    w<-sign(zeta) * (2 *temp1)^{1/2}
    v<- zeta * (k2)^{1/2}
    
    Z.test<-w + 1/w * log(v/w)	
    
    
    if(Z.test > 0){
      pval<-pnorm(Z.test, lower.tail = FALSE,log.p=log.p)
    } else {
      pval= -pnorm(Z.test, lower.tail = TRUE,log.p=log.p)
    }	
  } else {
    if(log.p)
    {
      pval<- -Inf
    } else {
      pval<-0
    }
  }
  
  return(pval)
}
K1_adj<-function(t, mu, g, q)
{
  n.t<-length(t)	
  out<-rep(0,n.t)
  
  for(i in 1:n.t){
    t1<-t[i]
    temp1<-(1 - mu)* exp(-g * t1) + mu
    temp2<-mu *g
    out[i]<-sum(temp2/temp1)-q
  }
  return(out)
}
K2<-function(t, mu, g)
{
  n.t<-length(t)
  out<-rep(0,n.t)
  
  for(i in 1:n.t){
    t1<-t[i]
    temp1<-((1 - mu)* exp(-g * t1) + mu)^2
    temp2<-(1-mu) * mu * g^2 * exp(-g*t1)
    out[i]<-sum(temp2/temp1,na.rm=TRUE)
  }
  return(out)
}

