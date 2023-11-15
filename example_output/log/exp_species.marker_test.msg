Error in diag(new_kin) : long vectors not supported yet: array.c:2192 Calls: micro_glmm -> gen_sp_Sigma -> diag 
Traceback: 
[1] "3: diag(new_kin)"                                                     
[2] "2: gen_sp_Sigma(W1, tauVecNew, GRM)"                                  
[3] "1: micro_glmm(glmm_fit, glm_fit0, GRM, copy_number_df, SPA = spa_opt)"

