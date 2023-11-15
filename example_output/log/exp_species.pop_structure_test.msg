Error in diag(new_kin) : long vectors not supported yet: array.c:2192 Calls: pop_structure_test ... Get_Coef -> getCoefficients -> gen_sp_Sigma -> diag 
Traceback: 
[1] "5: diag(new_kin)"                                                               
[2] "4: gen_sp_Sigma(W, tau, GRM)"                                                   
[3] "3: getCoefficients(Y, X, W, tau, GRM)"                                          
[4] "2: Get_Coef(y, X, tau, GRM, family, alpha0, eta0, offset, maxiter = maxiter, "  
[5] "       verbose = verbose, tol.coef = tol, write_log = write_log)"               
[6] "1: pop_structure_test(glm_fit0, GRM, species_id = s_id, tau = c(phi0, "         
[7] "       tau0), maxiter = maxiter, verbose = verbose, tol = tol, log_file = TRUE)"

