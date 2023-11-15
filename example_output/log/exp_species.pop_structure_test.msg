Error in t.default(sigma_X %*% alpha) : argument is not a matrix Calls: pop_structure_test -> Get_Coef -> getCoefficients -> t -> t.default 
Traceback: 
[1] "5: t.default(sigma_X %*% alpha)"                                                
[2] "4: t(sigma_X %*% alpha)"                                                        
[3] "3: getCoefficients(Y, X, W, tau, GRM)"                                          
[4] "2: Get_Coef(y, X, tau, GRM, family, alpha0, eta0, offset, maxiter = maxiter, "  
[5] "       verbose = verbose, tol.coef = tol, write_log = write_log)"               
[6] "1: pop_structure_test(glm_fit0, GRM, species_id = s_id, tau = c(phi0, "         
[7] "       tau0), maxiter = maxiter, verbose = verbose, tol = tol, log_file = TRUE)"

