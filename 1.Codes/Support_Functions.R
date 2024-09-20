# Bayesian P value
P_se = function(x){
 M_x = mean(x)

 if(M_x<0){
   N_x = length(x)
   P_x = length(which(x>0))
  } else{
   N_x = length(x)
   P_x = length(which(x<0))
  }
 return(P_x/N_x)
    }

# Summarize samples from MCMC
process_samples = function(z, y1, y2, x){
 res = c(y1, median(z), HPDI(z, 0.95), mean(z), sd(z), y2, P_se(z), x)
 return(res)
}

prepare_dat_ANTs_Rate_data = function(A){
  tie_strength = A$outcomes[,,1] / (1 + A$sampled_exposure)
  colnames(tie_strength) = rownames(tie_strength) = 1:ncol(tie_strength)
  df = ANTs:::df.create(tie_strength)
  df$outstrength = met.outstrength(tie_strength)
  df$instrength = met.instrength(tie_strength)
  df$Coloration = A$individual_predictors$Coloration
  df$Shyness = A$individual_predictors$Shyness
  df$exposure = A$sampled
  return(df)
}

prepare_dat_ANTs_SRI_data = function(A){
  m = A$outcomes[,,1]
  fa = fb = matrix(0, ncol = nrow(m), nrow = nrow(m))
  
  for (a in 1:nrow(fa)){
    fa[a, ] = A$sampled
    fb[, a] = A$sampled
  }
  
  ya = abs(fa - m)
  yb = abs(fb - m)
  
  sri = (m / (m + ya + yb))
  sri[is.na(sri)] = 0
  colnames(sri) = rownames(sri) = 1:ncol(sri)

  df = data.frame('id' = 1:ncol(m))
  
  df$outstrength = met.outstrength(sri)
  df$instrength = met.instrength(sri)
  df$Coloration = A$individual_predictors$Coloration
  df$Shyness = A$individual_predictors$Shyness
  df$exposure = A$sampled
  return(df)
}


fit_ANTs_Rate_Unweighted = function(model_dat, x){
  d = prepare_dat_ANTs_Rate_data(model_dat)
  # Unweigthed
  test_out = lm(outstrength ~ Coloration + Shyness, data = d)
  s = summary(test_out)
  mean = median = s$coefficients[-1,1]
  sd = s$coefficients[-1,2]
  p = s$coefficients[-1,4]
  conf = confint(test_out)
  L = conf[-1, 1]
  H = conf[-1,2]
  
  res = data.frame("Variable" = paste0("Focal, ",names(mean)),
                 "Median" = median,
                 "L" = L, 
                 "H" = H,
                 "Mean" = mean,
                 "SD" = sd,
                 "Model" = "ANTs Rates",
                 "P" = p)
  
  test_in = lm(instrength ~ Coloration + Shyness, data = d)
  s = summary(test_in)
  mean = median = s$coefficients[-1,1]
  sd = s$coefficients[-1,2]
  p = s$coefficients[-1,4]
  
  conf = confint(test_in)
  L = conf[-1, 1]
  H = conf[-1,2]
  
  
  res2 = data.frame("Variable" = paste0("Target, ",names(mean)),
                   "Median" = median,
                   "L" = L, 
                   "H" = H,
                   "Mean" = mean,
                   "SD" = sd,
                   "Model" = "ANTs Rates",
                   "P" = p)
  
  result1 = rbind(res, res2)
  result1$Setting = x
  rownames(result1) = NULL

  return(result1)
}
  
fit_ANTs_Rate_Weighted = function(model_dat, x){
  d = prepare_dat_ANTs_Rate_data(model_dat)
  test_out = lm(outstrength ~ Coloration + Shyness, data = d, weights = d$exposure)
  s = summary(test_out)
  mean = median = s$coefficients[-1,1]
  sd = s$coefficients[-1,2]
  p = s$coefficients[-1,4]
  conf = confint(test_out)
  L = conf[-1, 1]
  H = conf[-1, 2]
  
  res = data.frame("Variable" = paste0("Focal, ",names(mean)),
                   "Median" = median,
                   "L" = L, 
                   "H" = H,
                   "Mean" = mean,
                   "SD" = sd,
                   "Model" = "ANTs Rates Weighted",
                   "P" = p)
  
  test_in = lm(instrength ~ Coloration + Shyness, data = d, weights = d$exposure)
  s = summary(test_in)
  mean = median = s$coefficients[-1,1]
  sd = s$coefficients[-1,2]
  p = s$coefficients[-1,4]
  
  conf = confint(test_in)
  L = conf[-1, 1]
  H = conf[-1,2]
  
  res2 = data.frame("Variable" = paste0("Target, ",names(mean)),
                    "Median" = median,
                    "L" = L, 
                    "H" = H,
                    "Mean" = mean,
                    "SD" = sd,
                    "Model" = "ANTs Rates Weighted",
                    "P" = p)
  
  result2 = rbind(res, res2)
  rownames(result2) = NULL
  result2$Setting = x

  return(result2)
}

fit_ANTs_Rate_Perm_Unweighted = function(model_dat, x){
  d = prepare_dat_ANTs_Rate_data(model_dat)
  perm = perm.net.nl.str(d, labels = c('Coloration', 'Shyness'), nperm = 10000, progress = FALSE)
  
  R_out = R_in  = NULL
  for(a in 1:length(perm)) {
    test_out = lm(outstrength ~ Coloration + Shyness, data = perm[[a]])
    test_in = lm(instrength ~ Coloration + Shyness, data = perm[[a]])

    s = summary(test_out)
    R_out = rbind(R_out, s$coefficients[-1,1])

    s = summary(test_in)
    R_in = rbind(R_in, s$coefficients[-1,1])
  }
  
  p_out_C = ANTs:::stat.p(R_out[,1])
  p_in_C = ANTs:::stat.p(R_in[,1])
  
  p_out_S = ANTs:::stat.p(R_out[,2])
  p_in_S = ANTs:::stat.p(R_in[,2])
  
  mean_out_C = R_out[1,1]
  mean_in_C = R_in[1,1]
  
  mean_out_S = R_out[1,2]
  mean_in_S = R_in[1,2]

  
  res1 = data.frame("Variable" = paste0("Focal, ","Coloration"),
                   "Median" = mean_out_C,
                   "L" = NA, 
                   "H" = NA,
                   "Mean" = mean_out_C,
                   "SD" = NA,
                   "Model" = "ANTs Rates Permuted",
                   "P" = p_out_C[3])

  res2 = data.frame("Variable" = paste0("Focal, ","Shyness"),
                   "Median" = mean_out_S,
                   "L" = NA, 
                   "H" = NA,
                   "Mean" = mean_out_S,
                   "SD" = NA,
                   "Model" = "ANTs Rates Permuted",
                   "P" = p_out_S[3])
  
  res3 = data.frame("Variable" = paste0("Target, ","Coloration"),
                   "Median" = mean_in_C,
                   "L" = NA, 
                   "H" = NA,
                   "Mean" = mean_in_C,
                   "SD" = NA,
                   "Model" = "ANTs Rates Permuted",
                   "P" = p_in_C[3])

  res4 = data.frame("Variable" = paste0("Target, ","Shyness"),
                   "Median" = mean_in_S,
                   "L" = NA, 
                   "H" = NA,
                   "Mean" = mean_in_S,
                   "SD" = NA,
                   "Model" = "ANTs Rates Permuted",
                   "P" = p_in_S[3])
  
  res = rbind(res1, res2, res3, res4)
  row.names(res) = NULL
  res$Setting = x
  return(res)

}


fit_ANTs_Rate_Perm_Weighted = function(model_dat, x){
  d = prepare_dat_ANTs_Rate_data(model_dat)
  perm = perm.net.nl.str(d, labels = c('Coloration', 'Shyness'), nperm = 10000, progress = FALSE)
  
  R_out = R_in  = NULL
  for(a in 1:length(perm)) {
    test_out = lm(outstrength ~ Coloration + Shyness, data = perm[[a]], weights = d$exposure)
    test_in = lm(instrength ~ Coloration + Shyness, data = perm[[a]], weights = d$exposure)

    s = summary(test_out)
    R_out = rbind(R_out, s$coefficients[-1,1])

    s = summary(test_in)
    R_in = rbind(R_in, s$coefficients[-1,1])
  }
  
  p_out_C = ANTs:::stat.p(R_out[,1])
  p_in_C = ANTs:::stat.p(R_in[,1])
  
  p_out_S = ANTs:::stat.p(R_out[,2])
  p_in_S = ANTs:::stat.p(R_in[,2])
  
  mean_out_C = R_out[1,1]
  mean_in_C = R_in[1,1]
  
  mean_out_S = R_out[1,2]
  mean_in_S = R_in[1,2]

  
  res1 = data.frame("Variable" = paste0("Focal, ","Coloration"),
                   "Median" = mean_out_C,
                   "L" = NA, 
                   "H" = NA,
                   "Mean" = mean_out_C,
                   "SD" = NA,
                   "Model" = "ANTs Rates Weighted Permuted",
                   "P" = p_out_C[3])

  res2 = data.frame("Variable" = paste0("Focal, ","Shyness"),
                   "Median" = mean_out_S,
                   "L" = NA, 
                   "H" = NA,
                   "Mean" = mean_out_S,
                   "SD" = NA,
                   "Model" = "ANTs Rates Weighted Permuted",
                   "P" = p_out_S[3])
  
  res3 = data.frame("Variable" = paste0("Target, ","Coloration"),
                   "Median" = mean_in_C,
                   "L" = NA, 
                   "H" = NA,
                   "Mean" = mean_in_C,
                   "SD" = NA,
                   "Model" = "ANTs Rates Weighted Permuted",
                   "P" = p_in_C[3])

  res4 = data.frame("Variable" = paste0("Target, ","Shyness"),
                   "Median" = mean_in_S,
                   "L" = NA, 
                   "H" = NA,
                   "Mean" = mean_in_S,
                   "SD" = NA,
                   "Model" = "ANTs Rates Weighted Permuted",
                   "P" = p_in_S[3])
  
  res = rbind(res1, res2, res3, res4)
  row.names(res) = NULL
  res$Setting = x
  return(res)

}


fit_STRAND_basic = function(model_dat, x){
# Basic model
fit =  fit_block_plus_social_relations_model(data=model_dat,
                              block_regression = ~ Group,
                              focal_regression = ~ Coloration + Shyness,
                              target_regression = ~ Coloration + Shyness,
                              dyad_regression = ~  SizeDiff,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 600, iter_sampling = 600,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

res1 = summarize_strand_results(fit)
res1$summary$Model = "STRAND Basic"
res1$summary$P = NA 

res1$summary$P[which(res1$summary$Variable=="focal effects coeffs (out-degree), Coloration")] = P_se(res1$samples$srm_model_samples$focal_coeffs[,1])
res1$summary$P[which(res1$summary$Variable=="focal effects coeffs (out-degree), Shyness")] = P_se(res1$samples$srm_model_samples$focal_coeffs[,2])

res1$summary$P[which(res1$summary$Variable=="target effects coeffs (in-degree), Coloration")] = P_se(res1$samples$srm_model_samples$target_coeffs[,1])
res1$summary$P[which(res1$summary$Variable=="target effects coeffs (in-degree), Shyness")] = P_se(res1$samples$srm_model_samples$target_coeffs[,2])
 
resb = res1$summary 
resb$Setting = x
colnames(resb) = c("Variable", "Median", "L", "H", "Mean", "SD", "Model", "P","Setting")

in_set = c("focal effects coeffs (out-degree), Coloration","focal effects coeffs (out-degree), Shyness","target effects coeffs (in-degree), Coloration","target effects coeffs (in-degree), Shyness")
resb2 = resb[which(resb$Variable %in% in_set),]

resb2$Variable = c("Focal, Coloration", "Focal, Shyness", "Target, Coloration", "Target, Shyness")

return(resb2)
}


fit_STRAND_ME = function(model_dat, x){                                                               
 fit =  fit_block_plus_social_relations_model_with_measurement_bias(data=model_dat,
                              block_regression = ~ Group,
                              focal_regression = ~ Coloration + Shyness,
                              target_regression = ~ Coloration + Shyness,
                              sampling_regression = ~ Coloration + Shyness,
                              censoring_regression = ~ Coloration + Shyness,
                              dyad_regression = ~  SizeDiff,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 600, iter_sampling = 600,
                                                          max_treedepth = NULL, adapt_delta = .98)
 )

res2 = summarize_bsrm_results_with_measurement_bias(fit)
res2$summary$Model = "STRAND ME"
res2$summary$P = NA 

res2$summary$P[which(res2$summary$Variable=="focal effects coeffs (out-degree), Coloration")] = P_se(res2$samples$srm_model_samples$focal_coeffs[,1])
res2$summary$P[which(res2$summary$Variable=="focal effects coeffs (out-degree), Shyness")] = P_se(res2$samples$srm_model_samples$focal_coeffs[,2])

res2$summary$P[which(res2$summary$Variable=="target effects coeffs (in-degree), Coloration")] = P_se(res2$samples$srm_model_samples$target_coeffs[,1])
res2$summary$P[which(res2$summary$Variable=="target effects coeffs (in-degree), Shyness")] = P_se(res2$samples$srm_model_samples$target_coeffs[,2])

resb = res2$summary 
resb$Setting = x
colnames(resb) = c("Variable", "Median", "L", "H", "Mean", "SD", "Model", "P","Setting")

in_set = c("focal effects coeffs (out-degree), Coloration","focal effects coeffs (out-degree), Shyness","target effects coeffs (in-degree), Coloration","target effects coeffs (in-degree), Shyness")
resb2 = resb[which(resb$Variable %in% in_set),]

resb2$Variable = c("Focal, Coloration", "Focal, Shyness", "Target, Coloration", "Target, Shyness")

return(resb2)
}

fit_AMEN = function(model_dat, x){
  dyadic_preds = array(NA, c(model_dat$N_id, model_dat$N_id,2))
  dyadic_preds[,,1] = model_dat$dyadic_predictors$SizeDiff
  dyadic_preds[,,2] = model_dat$exposure

  row_preds = array(NA, c(model_dat$N_id,2))
  col_preds = array(NA, c(model_dat$N_id,2))

  row_preds[,1]= col_preds[,1] = model_dat$individual_predictors$Coloration 
  row_preds[,2]= col_preds[,2] = model_dat$individual_predictors$Shyness

  Y = model_dat$outcomes
  Y = Y[,,1]

  diag(Y)=NA
  res = ame(Y, Xdyad=dyadic_preds, Xrow=row_preds, Xcol=col_preds, family="nrm", plot=FALSE, print = FALSE, seed = 13, nscan = 20000, burn = 5000)

  res_df = matrix(NA, nrow=4, ncol=9)
  colnames(res_df) = c("Variable", "Median", "L", "H", "Mean", "SD", "Model", "P","Setting")

  res_df[1,] = process_samples(res$BETA[,2],"Focal, Coloration","AMEN", x)
  res_df[2,] = process_samples(res$BETA[,3],"Focal, Shyness","AMEN", x)
  res_df[3,] = process_samples(res$BETA[,4],"Target, Coloration","AMEN", x)
  res_df[4,] = process_samples(res$BETA[,5],"Target, Shyness","AMEN", x)

  res_df = data.frame(res_df)

return(res_df)
}


fit_bison = function(model_dat, x){
 # STRAND and other methods use matrices, rather than edge-lists, so we need to organize the data for bisonR

 # Create matrix form of N-vector data objects
 node_i_id = node_j_id = Shyness_j = Coloration_j = Shyness_i = Coloration_i = matrix(NA, nrow=model_dat$N_id, ncol=model_dat$N_id)
 
 # Loop over individuals, and fill matrices
 for(i in 1:model_dat$N_id){
  node_i_id[,i] = 1:model_dat$N_id
  node_j_id[i,] = 1:model_dat$N_id

  Coloration_i[,i] = model_dat$individual_predictors$Coloration
  Coloration_j[i,] = model_dat$individual_predictors$Coloration

  Shyness_i[,i] = model_dat$individual_predictors$Shyness
  Shyness_j[i,] = model_dat$individual_predictors$Shyness
  }

 # Turn all of the matrix data into long-form
 df_bison = data.frame(event=c(model_dat$outcomes), duration=c(model_dat$exposure), size_diff=c(model_dat$dyadic_predictors$SizeDiff),
                      node_i_id = c(node_i_id), node_j_id=c(node_j_id), coloration_i=c(Coloration_i), coloration_j=c(Coloration_j),
                      shyness_i=c(Shyness_i), shyness_j=c(Shyness_j)
                        )

 # Drop the rows that refer to ego-ego relationships
 df_bison = df_bison[which(df_bison$node_i_id != df_bison$node_j_id),]


#################################################### This comes from the tutorial
 # Fit Bison Model
  fit_edge = bison_model(
   (event | duration) ~ dyad(node_i_id, node_j_id),
   data=df_bison,
   model_type="binary", 
   directed = TRUE
  )

 fit_brm = bison_brm(
  bison(edge_weight(node_i_id, node_j_id)) ~ size_diff + coloration_i + coloration_j + shyness_i + shyness_j,
  fit_edge,
  df_bison,
  num_draws = 10
  )

  res_df = matrix(NA, nrow=4, ncol=9) # store results
  colnames(res_df) = c("Variable", "Median", "L", "H", "Mean", "SD", "Model", "P","Setting")

  res_df[1,] = process_samples(c(rstan::extract(fit_brm$fit,pars="b_coloration_i")$b_coloration_i),"Focal, Coloration","bison", x)
  res_df[2,] = process_samples(c(rstan::extract(fit_brm$fit,pars="b_coloration_j")$b_coloration_j),"Target, Coloration","bison", x)
  res_df[3,] = process_samples(c(rstan::extract(fit_brm$fit,pars="b_shyness_i")$b_shyness_i),"Focal, Shyness","bison", x)
  res_df[4,] = process_samples(c(rstan::extract(fit_brm$fit,pars="b_shyness_j")$b_shyness_j),"Target, Shyness","bison", x)

  res_df = data.frame(res_df)
  return(res_df)
}


fit_asnipe = function(model_dat, x){
node_1_id = node_2_id = Shyness_j = Coloration_j = Shyness_i = Coloration_i = matrix(NA, nrow=model_dat$N_id, ncol=model_dat$N_id)

for(i in 1:model_dat$N_id){
  Coloration_i[,i] = model_dat$individual_predictors$Coloration
  Coloration_j[i,] = model_dat$individual_predictors$Coloration

  Shyness_i[,i] = model_dat$individual_predictors$Shyness
  Shyness_j[i,] = model_dat$individual_predictors$Shyness
}

SizeDiff = model_dat$dyadic_predictors$SizeDiff

Y = model_dat$outcomes / (model_dat$exposure + 1)
Y = Y[,,1]

scrap_qap = mrqap.dsp(Y ~ Coloration_i + Shyness_i + Coloration_j + Shyness_j + SizeDiff, directed="directed")

res1 = data.frame("Variable" = "Focal, Coloration",
                    "Median" = scrap_qap$coefficients[2],
                    "L" = NA, 
                    "H" = NA,
                    "Mean" = scrap_qap$coefficients[2],
                    "SD" = NA,
                    "Model" = "asnipe",
                    "P" = scrap_qap$P.values[2])

res2 = data.frame("Variable" = "Focal, Shyness",
                    "Median" = scrap_qap$coefficients[3],
                    "L" = NA, 
                    "H" = NA,
                    "Mean" = scrap_qap$coefficients[3],
                    "SD" = NA,
                    "Model" = "asnipe",
                    "P" = scrap_qap$P.values[3])

res3 = data.frame("Variable" = "Target, Coloration",
                    "Median" = scrap_qap$coefficients[4],
                    "L" = NA, 
                    "H" = NA,
                    "Mean" = scrap_qap$coefficients[4],
                    "SD" = NA,
                    "Model" = "asnipe",
                    "P" = scrap_qap$P.values[4])

res4 = data.frame("Variable" = "Target, Shyness",
                    "Median" = scrap_qap$coefficients[5],
                    "L" = NA, 
                    "H" = NA,
                    "Mean" = scrap_qap$coefficients[5],
                    "SD" = NA,
                    "Model" = "asnipe",
                    "P" = scrap_qap$P.values[5])

 res = rbind(res1, res2, res3, res4)
 res$Setting = x

 rownames(res) =NULL
 return(res)
}


fit_lm_Rate_Unweighted = function(model_dat, x){
  # Create matrix form of N-vector data objects
 node_i_id = node_j_id = Shyness_j = Coloration_j = Shyness_i = Coloration_i = matrix(NA, nrow=model_dat$N_id, ncol=model_dat$N_id)
 
 # Loop over individuals, and fill matrices
 for(i in 1:model_dat$N_id){
  node_i_id[,i] = 1:model_dat$N_id
  node_j_id[i,] = 1:model_dat$N_id

  Coloration_i[,i] = model_dat$individual_predictors$Coloration
  Coloration_j[i,] = model_dat$individual_predictors$Coloration

  Shyness_i[,i] = model_dat$individual_predictors$Shyness
  Shyness_j[i,] = model_dat$individual_predictors$Shyness
  }

 # Turn all of the matrix data into long-form
 df_lm = data.frame(event=c(model_dat$outcomes), duration=c(model_dat$exposure), size_diff=c(model_dat$dyadic_predictors$SizeDiff),
                      node_i_id = c(node_i_id), node_j_id=c(node_j_id), coloration_i=c(Coloration_i), coloration_j=c(Coloration_j),
                      shyness_i=c(Shyness_i), shyness_j=c(Shyness_j)
                        )

 # Drop the rows that refer to ego-ego relationships
 df_lm = df_lm[which(df_lm$node_i_id != df_lm$node_j_id),]

 df_lm$ratio = df_lm$event/(1+df_lm$duration)

 test_out = lm(ratio ~ coloration_i + coloration_j + shyness_i + shyness_j + size_diff, data = df_lm)
 s = summary(test_out)
 mean = median = s$coefficients[-1,1]
 sd = s$coefficients[-1,2]
 p = s$coefficients[-1,4]
 conf = confint(test_out)
 L = conf[-1, 1]
 H = conf[-1,2]
  
  res = data.frame("Variable" = names(mean),
                 "Median" = median,
                 "L" = L, 
                 "H" = H,
                 "Mean" = mean,
                 "SD" = sd,
                 "Model" = "lm Rates",
                 "P" = p)
  res = res[1:4,]
  res$Variable  = c("Focal, Coloration","Target, Coloration","Focal, Shyness","Target, Shyness")
  res$Setting = x
  rownames(res) = NULL

  return(res)
}




prepare_series_plots = function(results_basic){
 res_all_basic = do.call(rbind, results_basic)

 valid_names = c("AMEN", "ANTs Rates", "ANTs Rates Permuted", "asnipe", "bison", "lm Rates", "STRAND Basic", "STRAND ME")

 res_all_basic = res_all_basic[which(res_all_basic$Model %in% valid_names),]


 res_all_basic$P2 = sign(as.numeric(res_all_basic$Mean))*as.numeric(res_all_basic$P)
 res_all_basic$Median = as.numeric(res_all_basic$Median)
 res_all_basic$L = as.numeric(res_all_basic$L)
 res_all_basic$H = as.numeric(res_all_basic$H)

 res_all_basic$Setting = as.numeric(res_all_basic$Setting)
 res_all_basic$P = as.numeric(res_all_basic$P)

 res_all_basic$D = as.numeric(res_all_basic$H)-as.numeric(res_all_basic$L)
 res_all_basic$Median2 = as.numeric(res_all_basic$Median)/res_all_basic$D
 res_all_basic$L2 = as.numeric(res_all_basic$L)/res_all_basic$D
 res_all_basic$H2 = as.numeric(res_all_basic$H)/res_all_basic$D

 res_all_basic$Variable = gsub(" ", "", res_all_basic$Variable)

 res_all_basic$Model = factor(res_all_basic$Model)
 res_all_basic$Model = factor(res_all_basic$Model, levels=c(
  "ANTs Rates", "ANTs Rates Permuted", "lm Rates", "AMEN", "asnipe", "bison", "STRAND Basic",  "STRAND ME"))

 res_all_basic$Variable = as.character(res_all_basic$Variable)
 res_all_basic$Variable[which(res_all_basic$Variable=="Focal,Coloration")] = "Coloration, out-degree"
 res_all_basic$Variable[which(res_all_basic$Variable=="Target,Coloration")] = "Coloration, in-degree"
 res_all_basic$Variable[which(res_all_basic$Variable=="Focal,Shyness")] = "Shyness, out-degree"
 res_all_basic$Variable[which(res_all_basic$Variable=="Target,Shyness")] = "Shyness, in-degree"
 res_all_basic$Variable = factor(res_all_basic$Variable)
 res_all_basic$Variable = factor(res_all_basic$Variable, levels=c("Coloration, out-degree","Coloration, in-degree","Shyness, out-degree","Shyness, in-degree"))

 return(res_all_basic)
}

prepare_robustness_plots = function(results_basic){
 res_all_basic = results_basic

 valid_names = c("AMEN", "ANTs Rates", "ANTs Rates Permuted", "asnipe", "bison", "lm Rates", "STRAND Basic", "STRAND ME")

 res_all_basic = res_all_basic[which(res_all_basic$Model %in% valid_names),]

 res_all_basic$P2 = sign(as.numeric(res_all_basic$Mean))*as.numeric(res_all_basic$P)
 res_all_basic$Median = as.numeric(res_all_basic$Median)
 res_all_basic$L = as.numeric(res_all_basic$L)
 res_all_basic$H = as.numeric(res_all_basic$H)

 res_all_basic$Setting = as.numeric(res_all_basic$Setting)
 res_all_basic$P = as.numeric(res_all_basic$P)

 res_all_basic$D = as.numeric(res_all_basic$H)-as.numeric(res_all_basic$L)
 res_all_basic$Median2 = as.numeric(res_all_basic$Median)/res_all_basic$D
 res_all_basic$L2 = as.numeric(res_all_basic$L)/res_all_basic$D
 res_all_basic$H2 = as.numeric(res_all_basic$H)/res_all_basic$D

 res_all_basic$Variable = gsub(" ", "", res_all_basic$Variable)


 res_all_basic$Variable = as.character(res_all_basic$Variable)
 res_all_basic$Variable[which(res_all_basic$Variable=="Focal,Coloration")] = "Coloration, out-degree"
 res_all_basic$Variable[which(res_all_basic$Variable=="Target,Coloration")] = "Coloration, in-degree"
 res_all_basic$Variable[which(res_all_basic$Variable=="Focal,Shyness")] = "Shyness, out-degree"
 res_all_basic$Variable[which(res_all_basic$Variable=="Target,Shyness")] = "Shyness, in-degree"
 res_all_basic$Variable = factor(res_all_basic$Variable)
 res_all_basic$Variable = factor(res_all_basic$Variable, levels=c("Coloration, out-degree","Coloration, in-degree","Shyness, out-degree","Shyness, in-degree"))

 res_all_basic$Model = factor(res_all_basic$Model)
 res_all_basic$Model = factor(res_all_basic$Model, levels=c(
  "ANTs Rates", "ANTs Rates Permuted", "lm Rates", "AMEN", "asnipe", "bison", "STRAND Basic",  "STRAND ME"))

 res_all_basic=res_all_basic[which(!is.na(res_all_basic$Setting)),]

 res_all_basic$Setting = factor(res_all_basic$Setting)
 levels(res_all_basic$Setting)=paste0("\u03ba = ",levels(res_all_basic$Setting))

 return(res_all_basic)
}

aggregate_error_rates = function(res_all_basic){
res_all_basic$Result = NA 

for(i in 1:length(res_all_basic$P)){
  if(res_all_basic$Variable[i] %in% c("Coloration, out-degree","Coloration, in-degree")){
    res_all_basic$Result[i] = ifelse(res_all_basic$P[i] < 0.05, 1, 0)
  } 
   if(res_all_basic$Variable[i] %in% c("Shyness, out-degree","Shyness, in-degree")){
    res_all_basic$Result[i] = ifelse(res_all_basic$P[i] > 0.05, 1, 0)
  }
}

res_all_basic_agg = aggregate(Result ~ Model+Variable+Setting, res_all_basic, mean)
res_all_basic_agg$ResultWrong = 1 - res_all_basic_agg$Result

return(res_all_basic_agg)
}

aggregate_p_values = function(res_all_basic){
res_all_basic_agg = aggregate(P ~ Model+Variable+Setting, res_all_basic, function(x) mean(x, na.rm=TRUE), na.action = na.pass)

return(res_all_basic_agg)
}



