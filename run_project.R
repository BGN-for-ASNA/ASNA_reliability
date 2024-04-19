source('1.Codes/1.load_packages.R')
source('1.Codes/2.data_simulation.R')
source('1.Codes/support_functions.R')
source('1.Codes/2.1.STRAND_censoring.R')
source('1.Codes/2.2.STRAND fit_social_relations_model_bias.R')

# 1. Evaluation of estimation by strand on simulation-----------------------------------------
#' @param effect A vector with effect to simulate.
#' @param exposure A vector with exposure to simulate. Vector needs to be of the same length as the argument 'effect'.
#' @param censoring A vector with censoring to simulate. Vector needs to be of the same length as the argument 'effect'.
#' @description
#' Generate the following tests: 1) Without exposure or censoring bias, 2) With exposure and no censoring bias, 3) With exposure and censoring bias.
#' For all tests, the argument 'effect' is used. So, for test 2 simulation 1, consider value 1 for both effect and exposure, similarly for test 3.
#' 
estimation.evaluation <- function(effect = seq(-5, 5, by = 1), exposure = seq(-5, 5, by = 1), censoring = seq(-5, 5, by = 1)){
  ## 1.1. Without exposure or censoring bias ----------------------------
  R = NULL
  for (a in 1:length(effect)) {
    N_id = 50
    # Sender-receiver parameter
    individual_effects_sr = c(-1, effect[a])
    sr_sigma = c(1.4, 0.8)
    sr_rho = 0.5
    # Dyadic parameter
    dr_sigma = 1.2
    dr_rho = 0.8
    # Observation bias parameter
    exposure_bias = FALSE
    exposure_mu = 1.9
    exposure_sigma = 2.9
    exposure_max = 40
    exposure_effects = 0
    # Censoring bias parameter
    simulate_censoring = TRUE
    censoring_mu = 0
    censoring_sigma = 0
    censoring_effects = 0
    
    ########################## First thing is create invidiual level and block data---------------------
    
    ## Individual-specific covariate
    Individual_Factor <- matrix(rnorm(N_id, 0, 1), nrow = N_id, ncol = 1)
    
    ## Individual-specific covariate
    Clique <- sample(1:3, N_id, replace = TRUE)
    
    V <- 1 # One blocking variable
    G <- 3 # Three categories in this variable
    
    B <- matrix(-8, nrow = G, ncol = G)
    diag(B) <- -4.5
    B[1, 3] <- -5.9
    B[3, 2] <- -6.9
    
    ########################## Then create simulated data---------------------
    A <- simulate_sbm_plus_srm_network_with_measurement_bias(
      N_id = N_id,
      B = list(B = B),
      V = V,
      groups = data.frame(Clique = factor(Clique)),
      individual_predictors = matrix(Individual_Factor, nrow = N_id, ncol = 1),
      individual_effects = matrix(individual_effects_sr, ncol = 1, nrow = 2),
      dyadic_predictors = array(rnorm(N_id * N_id, 0, 1), c(N_id, N_id, 1)),
      dyadic_effects = c(0.00),
      sr_mu = c(0, 0),
      sr_sigma = sr_sigma,
      sr_rho = sr_rho,
      dr_mu = c(0, 0),
      dr_sigma = dr_sigma,
      dr_rho = dr_rho,
      exposure_mu = exposure_mu,
      exposure_sigma = exposure_sigma,
      exposure_max = exposure_max,
      censoring_mu = censoring_mu,
      censoring_sigma = censoring_sigma,
      exposure_predictors = matrix(Individual_Factor, nrow = N_id, ncol = 1),
      censoring_predictors = matrix(Individual_Factor, nrow = N_id, ncol = 1),
      simulate_censoring = simulate_censoring,
      simulate_interactions = TRUE,
      exposure_effects = exposure_effects,
      censoring_effects = censoring_effects
    )
    
    ###################################### STRAND---------------------
    nets <- list(Grooming = A$network)
    indiv <- data.frame(Hairy = Individual_Factor)
    block <- data.frame(Clique = factor(Clique))
    exposure_nets <- list(Exposure = A$true_samps)
    
    model_dat <- make_strand_data_censoring(
      outcome = nets,
      individual_covariates = indiv,
      block_covariates = block,
      outcome_mode = "binomial",
      exposure = exposure_nets, detected = A$detected, trials = A$trials
    )
    
    fit <- fit_social_relations_model_bias(
      data = model_dat,
      block_regression = ~ Clique,
      focal_regression = ~ Hairy,
      target_regression = ~ Hairy,
      censoring_regression = ~ Hairy,
      dyad_regression = ~1,
      mode = "mcmc", return_predicted_network = TRUE,
      stan_mcmc_parameters = list(
        chains = 1, parallel_chains = 1, refresh = 1000,
        iter_warmup = 1000, iter_sampling = 1000,
        max_treedepth = NULL, adapt_delta = 0.98
      )
    )
    
    res <- summarize_strand_results(fit)
    R[[a]] = res
  }
  
  R2 = NULL
  for (a in 1:length(R)) {
    tmp = R[[a]]$summary
    r = as.numeric(tmp[tmp$Variable %in% 'target effects coeffs (in-degree), Hairy',2])
    
    R2 = rbind(R2, data.frame('simulated effect' = effect[a], 'estimated effect' = r ))
  }
  
  p1 = ggplot(R2, aes(x = simulated.effect, y = estimated.effect))+ geom_point()+xlab('Simulated effect')+ylab('Estimated effect')
  
  ## 1.2. With exposure and no censoring bias ----------------------------
  R3 = NULL
  for (a in 1:length(effect)) {
    N_id = 50
    # Sender-receiver parameter
    individual_effects_sr = c(-1, effect[a])
    sr_sigma = c(1.4, 0.8)
    sr_rho = 0.5
    # Dyadic parameter
    dr_sigma = 1.2
    dr_rho = 0.8
    # Observation bias parameter
    exposure_bias = TRUE
    exposure_mu = 1.9
    exposure_sigma = 2.9
    exposure_max = 40
    exposure_effects = exposure[a]
    # Censoring bias parameter
    simulate_censoring = TRUE
    censoring_mu = 0
    censoring_sigma = 0
    censoring_effects = 0
    
    ########################## First thing is create invidiual level and block data---------------------
    
    ## Individual-specific covariate
    Individual_Factor <- matrix(rnorm(N_id, 0, 1), nrow = N_id, ncol = 1)
    
    ## Individual-specific covariate
    Clique <- sample(1:3, N_id, replace = TRUE)
    
    V <- 1 # One blocking variable
    G <- 3 # Three categories in this variable
    
    B <- matrix(-8, nrow = G, ncol = G)
    diag(B) <- -4.5
    B[1, 3] <- -5.9
    B[3, 2] <- -6.9
    
    ########################## Then create simulated data---------------------
    A <- simulate_sbm_plus_srm_network_with_measurement_bias(
      N_id = N_id,
      B = list(B = B),
      V = V,
      groups = data.frame(Clique = factor(Clique)),
      individual_predictors = matrix(Individual_Factor, nrow = N_id, ncol = 1),
      individual_effects = matrix(individual_effects_sr, ncol = 1, nrow = 2),
      dyadic_predictors = array(rnorm(N_id * N_id, 0, 1), c(N_id, N_id, 1)),
      dyadic_effects = c(0.00),
      sr_mu = c(0, 0),
      sr_sigma = sr_sigma,
      sr_rho = sr_rho,
      dr_mu = c(0, 0),
      dr_sigma = dr_sigma,
      dr_rho = dr_rho,
      exposure_mu = exposure_mu,
      exposure_sigma = exposure_sigma,
      exposure_max = exposure_max,
      censoring_mu = censoring_mu,
      censoring_sigma = censoring_sigma,
      exposure_predictors = matrix(Individual_Factor, nrow = N_id, ncol = 1),
      censoring_predictors = matrix(Individual_Factor, nrow = N_id, ncol = 1),
      simulate_censoring = simulate_censoring,
      simulate_interactions = TRUE,
      exposure_effects = exposure_effects,
      censoring_effects = censoring_effects
    )
    
    ###################################### STRAND---------------------
    nets <- list(Grooming = A$network)
    indiv <- data.frame(Hairy = Individual_Factor)
    block <- data.frame(Clique = factor(Clique))
    exposure_nets <- list(Exposure = A$true_samps)
    
    model_dat <- make_strand_data_censoring(
      outcome = nets,
      individual_covariates = indiv,
      block_covariates = block,
      outcome_mode = "binomial",
      exposure = exposure_nets, detected = A$detected, trials = A$trials
    )
    
    fit <- fit_social_relations_model_bias(
      data = model_dat,
      block_regression = ~ Clique,
      focal_regression = ~ Hairy,
      target_regression = ~ Hairy,
      censoring_regression = ~ Hairy,
      dyad_regression = ~1,
      mode = "mcmc", return_predicted_network = TRUE,
      stan_mcmc_parameters = list(
        chains = 1, parallel_chains = 1, refresh = 1000,
        iter_warmup = 1000, iter_sampling = 1000,
        max_treedepth = NULL, adapt_delta = 0.98
      )
    )
    
    res <- summarize_strand_results(fit)
    R3[[a]] = res
  }
  
  R4 = NULL
  for (a in 1:length(R3)) {
    tmp = R3[[a]]$summary
    r = as.numeric(tmp[tmp$Variable %in% 'target effects coeffs (in-degree), Hairy',2])
    
    R4 = rbind(R4, data.frame('simulated effect' = effect[a], 'estimated effect' = r ))
  }
  
  p2 = ggplot(R4, aes(x = simulated.effect, y = estimated.effect))+ geom_point()+xlab('Simulated effect')+ylab('Estimated effect')
  
  ## 1.3. With exposure and no censoring bias ----------------------------
  R5 = NULL
  for (a in 1:length(effect)) {
    N_id = 50
    # Sender-receiver parameter
    individual_effects_sr = c(-1, effect[a])
    sr_sigma = c(1.4, 0.8)
    sr_rho = 0.5
    # Dyadic parameter
    dr_sigma = 1.2
    dr_rho = 0.8
    # Observation bias parameter
    exposure_bias = FALSE
    exposure_mu = 1.9
    exposure_sigma = 2.9
    exposure_max = 40
    exposure_effects = 0
    # Censoring bias parameter
    simulate_censoring = TRUE
    censoring_mu = 0
    censoring_sigma = 0
    censoring_effects = censoring[a]
    
    ########################## First thing is create invidiual level and block data---------------------
    
    ## Individual-specific covariate
    Individual_Factor <- matrix(rnorm(N_id, 0, 1), nrow = N_id, ncol = 1)
    
    ## Individual-specific covariate
    Clique <- sample(1:3, N_id, replace = TRUE)
    
    V <- 1 # One blocking variable
    G <- 3 # Three categories in this variable
    
    B <- matrix(-8, nrow = G, ncol = G)
    diag(B) <- -4.5
    B[1, 3] <- -5.9
    B[3, 2] <- -6.9
    
    ########################## Then create simulated data---------------------
    A <- simulate_sbm_plus_srm_network_with_measurement_bias(
      N_id = N_id,
      B = list(B = B),
      V = V,
      groups = data.frame(Clique = factor(Clique)),
      individual_predictors = matrix(Individual_Factor, nrow = N_id, ncol = 1),
      individual_effects = matrix(individual_effects_sr, ncol = 1, nrow = 2),
      dyadic_predictors = array(rnorm(N_id * N_id, 0, 1), c(N_id, N_id, 1)),
      dyadic_effects = c(0.00),
      sr_mu = c(0, 0),
      sr_sigma = sr_sigma,
      sr_rho = sr_rho,
      dr_mu = c(0, 0),
      dr_sigma = dr_sigma,
      dr_rho = dr_rho,
      exposure_mu = exposure_mu,
      exposure_sigma = exposure_sigma,
      exposure_max = exposure_max,
      censoring_mu = censoring_mu,
      censoring_sigma = censoring_sigma,
      exposure_predictors = matrix(Individual_Factor, nrow = N_id, ncol = 1),
      censoring_predictors = matrix(Individual_Factor, nrow = N_id, ncol = 1),
      simulate_censoring = simulate_censoring,
      simulate_interactions = TRUE,
      exposure_effects = exposure_effects,
      censoring_effects = censoring_effects
    )
    
    ###################################### STRAND---------------------
    nets <- list(Grooming = A$network)
    indiv <- data.frame(Hairy = Individual_Factor)
    block <- data.frame(Clique = factor(Clique))
    exposure_nets <- list(Exposure = A$true_samps)
    
    model_dat <- make_strand_data_censoring(
      outcome = nets,
      individual_covariates = indiv,
      block_covariates = block,
      outcome_mode = "binomial",
      exposure = exposure_nets, detected = A$detected, trials = A$trials
    )
    
    fit <- fit_social_relations_model_bias(
      data = model_dat,
      block_regression = ~ Clique,
      focal_regression = ~ Hairy,
      target_regression = ~ Hairy,
      censoring_regression = ~ Hairy,
      dyad_regression = ~1,
      mode = "mcmc", return_predicted_network = TRUE,
      stan_mcmc_parameters = list(
        chains = 1, parallel_chains = 1, refresh = 1000,
        iter_warmup = 1000, iter_sampling = 1000,
        max_treedepth = NULL, adapt_delta = 0.98
      )
    )
    
    res <- summarize_strand_results(fit)
    R5[[a]] = res
  }
  
  R6 = NULL
  for (a in 1:length(R3)) {
    tmp = R5[[a]]$summary
    r = as.numeric(tmp[tmp$Variable %in% 'target effects coeffs (in-degree), Hairy',2])
    
    R6 = rbind(R6, data.frame('simulated effect' = effect[a], 'estimated effect' = r ))
  }
  
  p3 = ggplot(R6, aes(x = simulated.effect, y = estimated.effect))+ geom_point()+xlab('Simulated effect')+ylab('Estimated effect')
  
  p = ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
  return(list(R, R3, R5, p ))
}
estimations = estimation.evaluation() 



# 2. Comparing approaches estimations-----------------------------------------
## 2.1.Without exposure or censoring bias-------------------- 
R1 = parallel_sim(Reps = 5, 
                 N = seq(50, 90, by = 5), 
                 tie_effect = seq(-5, 5, by = 0.5),
                 detect_effect = 0,
                 censoring_effects =  0)

## 2.1.With exposure and no censoring bias-------------------- 
R2 = parallel_sim(Reps = 5, 
                 N = seq(50, 90, by = 5), 
                 tie_effect = seq(-5, 5, by = 0.5),
                 detect_effect = seq(-5, 5, by = 0.5),
                 censoring_effects =  0)

## 2.1.With exposure and censoring bias-------------------- 
R3 = parallel_sim(Reps = 5, 
                 N = seq(50, 90, by = 5), 
                 tie_effect = seq(-5, 5, by = 0.5),
                 detect_effect = 0,
                 censoring_effects =  seq(-5, 5, by = 0.5))


tmp =split(R1, R1$tie_effect)
tmp = do.call("rbind",lapply(tmp, function(x){
  r = data.frame(tapply(x$`scale(hair)`, x$approach, mean))
  r$ci5 = tapply(x$`5 %`, x$approach, mean)
  r$ci5 = tapply(x$`95 %`, x$approach, mean)
  r$p = tapply(x$`p-value`, x$approach, mean)
  r$tie_effect = tapply(x$tie_effect, x$approach, unique)
  r$approach = rownames(r) 
  rownames(r) = NULL
  colnames(r)[1] = "Estimate"
  r
}))

tmp2 = tmp[tmp$tie_effect %in% max(unique(tmp$tie_effect)),]
ggplot(tmp, aes(x = tie_effect, y = p, group = approach, color = approach, label = approach))+geom_line()+
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1)+
  geom_label_repel(data = tmp2, aes(x = tie_effect, y = p, colour = approach),
                   nudge_x = 0.5,
                   nudge_y = 0.02,
                   segment.curvature = -1e-20,
                   force = 10,
                   direction = "y",
                   hjust = 0)+
  theme(legend.position = 'None', text = element_text(size=20))+xlab("Error type")






## Generate appendices-------------
#render("1.Codes/Appendix 1.R",
#       output_format = pdf_document(),
#       output_file = "2.Results/Appendices/1/Appendix 1.pdf")
#
#
##render(input = "1.Codes/Appendix 2.R",
##       output_format = pdf_document(),
##       output_file = "2.Results/Appendices/2/Appendix 2.pdf")
#
#test = parallel_sim(Reps = 5)
### False positives rates ------------
#### No differences in sociality, no biases----------
#result1 = simulations(Reps = 100, ncores = 100, 
#                      N_id =  seq(30, 90, by = 10),
#                      sr_rho =0.5, dr_sigma = 0.5, dr_rho = 0.5, sr_sigma = c(1,1),
#                      hairy_tie_effect = seq(-0.25, 0.25, by = 0.01),
#                      hairy_detect_effect = seq(0, 0, by = 0.5),
#                      BISON = FALSE,
#                      STRAND = T, 
#                      simulate.interactions = F,
#                      simulate.censoring = F,
#                      cens_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
#                      cens_slope = c(Inf,Inf),# No effect of individuals attributes
#                      blockModel = TRUE) # No block model
#write.csv(result1, "2.Results/Simulations/No differences in sociality, no biases.csv", row.names = FALSE)
#
#### No differences in sociality, exposure bias----------
#result2 = simulations(Reps = 100, ncores = 100, 
#                      exposure_sigma = 1, 
#                      N_id =  seq(30, 90, by = 10), 
#                      sr_rho =0.5, dr_sigma = 0.5, dr_rho = 0.5, sr_sigma = c(1,1),
#                      hairy_tie_effect = seq(-0.25, 0.25, by = 0.01),
#                      hairy_detect_effect = seq(-0.40, 0.40, by = 0.1),
#                      BISON = FALSE,
#                      STRAND = T, 
#                      simulate.interactions = F,
#                      simulate.censoring = F,
#                      cens_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
#                      cens_slope = c(Inf,Inf),# No effect of individuals attributes
#                      blockModel = TRUE) # No block model
#
#write.csv(result2, "2.Results/Simulations/No differences in sociality, exposure bias.csv", row.names = FALSE)
#
#### No differences in sociality, censoring bias----------
#result3 = simulations(Reps = 100, ncores = 100, 
#                      exposure_sigma = 1, 
#                      N_id =  seq(30, 90, by = 10), 
#                      hairy_tie_effect = seq(-0.25, 0.25, by = 0.01),
#                      hairy_detect_effect = seq(0, 0, by = 0.5),
#                      BISON = FALSE,
#                      STRAND = T, 
#                      simulate.interactions = F,
#                      simulate.censoring = T,
#                      cens_intercept = c(4,4), #invert log of inf = 1 of prob to observe interaction for both focal and alter
#                      cens_slope = c(4, 4),# No effect of individuals attributes
#                      blockModel = TRUE) # No block model
#
#write.csv(result3, "2.Results/Simulations/No differences in sociality, censoring bias.csv", row.names = FALSE)
#
#
#
### False Negatives rates ------------
#### Differences in sociality, no biases----------
#result4 = simulations(Reps = 100, ncores = 100, 
#                      N_id =  seq(30, 90, by = 10), 
#                      sr_rho =0.5, dr_sigma = 0.5, dr_rho = 0.5, sr_sigma = c(1,1),
#                      hairy_tie_effect =  c(seq(-0.51, -0.26, by = 0.01), seq(0.26, 0.51, by = 0.01)),
#                      hairy_detect_effect = seq(0, 0, by = 0.5),
#                      BISON = FALSE,
#                      STRAND = T, 
#                      simulate.interactions = F, 
#                      simulate.censoring = F,
#                      cens_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
#                      cens_slope = c(-Inf,-Inf),# No effect of individuals attributes
#                      blockModel = TRUE) # No block model
#write.csv(result4, "2.Results/Simulations/Differences in sociality, no biases.csv", row.names = FALSE)
#
#
### Differences in sociality, exposure bias----------
#result5 = simulations(Reps = 100, ncores = 100, 
#                      N_id =  seq(30, 90, by = 10), 
#                      sr_rho =0.5, dr_sigma = 0.5, dr_rho = 0.5, sr_sigma = c(1,1),
#                      hairy_tie_effect = c(seq(-0.51, -0.26, by = 0.01), seq(0.26, 0.51, by = 0.01)),
#                      hairy_detect_effect = seq(-0.80, -0.20, by = 0.5),
#                      BISON = FALSE,
#                      STRAND = T, 
#                      simulate.interactions = F, 
#                      simulate.censoring = F,
#                      cens_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
#                      cens_slope = c(-Inf,-Inf),# No effect of individuals attributes
#                      blockModel = TRUE) # No block model
#write.csv(result5, "2.Results/Simulations/Differences in sociality, exposure bias.csv", row.names = FALSE)
#
#
### Differences in sociality, censoring bias----------
#result6 = simulations(Reps = 100, ncores = 100, 
#                      N_id =  seq(30, 90, by = 10), 
#                      sr_rho =0.5, dr_sigma = 0.5, dr_rho = 0.5, sr_sigma = c(1,1),
#                      hairy_tie_effect =  c(seq(-0.51, -0.26, by = 0.01), seq(0.26, 0.51, by = 0.01)),
#                      hairy_detect_effect = seq(0, 0, by = 0.5),
#                      BISON = FALSE,
#                      STRAND = T, 
#                      simulate.interactions = F, 
#                      simulate.censoring = T, 
#                      cens_intercept = c(4,4), #invert log of inf = 1 of prob to observe interaction for both focal and alter
#                      cens_slope = c(4,4),# No effect of individuals attributes
#                      blockModel = TRUE) # No block model
#write.csv(result6, "2.Results/Simulations/Differences in sociality, censoring bias.csv", row.names = FALSE)
#save.image("2.Results/Simulations/Simulations.RData")
#
### Get results-------------
#### Rates of false negatives and false positives------------------
#error.rates <- function(result, threshold){
#  
#  lower = result[result$tie_effect < threshold &  result$tie_effect  > -threshold,]$p.value
#  upper = result[result$tie_effect >= threshold | result$tie_effect <= -threshold,]$p.value
#  
#  if(length(lower) != 0){
#    t1 = tapply(result[result$tie_effect < threshold &  result$tie_effect  > -threshold,]$p.value, result[result$tie_effect < threshold & result$tie_effect > -threshold,]$approach, function(x){sum(x <= 0.05)/length(x)})
#  }else{t1 = NULL}
#  if(length(upper) != 0){
#    t2 = tapply(result[result$tie_effect >= threshold | result$tie_effect <= -threshold,]$p.value, result[result$tie_effect >= threshold | result$tie_effect <= -threshold,]$approach, function(x){sum(x >= 0.05)/length(x)})
#  }else{t2 = NULL}
#  
#  if(!is.null(t1) & !is.null(t2)){
#    summary = data.frame("Approaches" = c(names(t1),  names(t2)), 
#                         "Error type" = c(rep('False positives', length(t1)),
#                                          rep('False negatives', length(t2))) ,
#                         "Percent" = c(t1, t2))
#    
#    summary$Percent = summary$Percent * 100
#    return(summary)
#  }
#  if(!is.null(t1) != 0 & is.null(t2)){
#    summary = data.frame("Approaches" = names(t1), 
#                         "Error type" = rep('False positives', length(t1)) ,
#                         "Percent" = t1)
#    
#    summary$Percent = summary$Percent * 100
#    return(summary)
#  }
#  if(is.null(t1) & !is.null(t2)){
#    summary = data.frame("Approaches" = names(t2), 
#                         "Error type" = rep('False negatives', length(t2)) ,
#                         "Percent" = t2)
#    
#    summary$Percent = summary$Percent * 100
#    return(summary)
#  }
#}
#get.rates <- function(path, threshold = 0.2){
#  files = list.files(path)
#  results = errors = NULL
#  for(a in 1:length(files)){
#    if(grepl('.csv', files[a], fixed = TRUE)){
#      tmp = read.csv(paste(path,files[a], sep = '/'))
#      type = gsub('.csv', '', files[a])
#      tmp$type = type
#      results = rbind(results, tmp)
#      
#      tmp =  error.rates(tmp, threshold = threshold)
#      tmp$type = type
#      errors = rbind(errors, tmp)
#    }
#  }
#  rownames(errors) = NULL
#  return(list(errors, results))
#}
#error = get.rates(path = '~/ASNA_reliability/2.Results/Simulations', threshold = 0.26)
#
### Plots------------------
#### Rates------------------
#tmp = error[[1]]
#tmp$Approaches = gsub("2.","",tmp$Approaches)
#tmp$Approaches = gsub("3.","",tmp$Approaches)
#
#### False positives----------------
#tmp2 = tmp[tmp$type %in% c("No differences in sociality, censoring bias", "No differences in sociality, exposure bias", "No differences in sociality, no biases"),]
#tmp2$type = ifelse(tmp2$type == "No differences in sociality, censoring bias", 'Censoring biases', tmp2$type)
#tmp2$type = ifelse(tmp2$type == "No differences in sociality, exposure bias", 'Exposure biases', tmp2$type)
#tmp2$type = ifelse(tmp2$type == "No differences in sociality, no biases", 'No biases', tmp2$type)
#
#tmp3 = tmp2[tmp2$type =='No biases',]
#ggplot(tmp2, aes(x = type, y = Percent, group = Approaches, colour  = Approaches, label = Approaches))+geom_point(aes(size = 2))+geom_line()+facet_grid(~ Error.type) +
#  geom_label_repel(data = tmp3, aes(x = type, y = Percent, colour = Approaches),
#                   nudge_x = 1,
#                   nudge_y = 8,
#                   segment.curvature = -1e-20,
#                   force = 10,
#                   direction = "y",
#                   hjust= 0)+
#  theme(legend.position = 'None', text = element_text(size=20))+xlab("Error type")+geom_hline(yintercept = 5, linetype = "dashed", color = 'red')
#
#tmp = error[[2]]
#tmp$approach = gsub("2.","",tmp$approach)
#tmp$approach = gsub("3.","",tmp$approach)
#
#ggplot(tmp[tmp$type == "No biases", ], aes(x = tie_effect,  y = z, group = sim, label = z))+
#  geom_point(aes(color = sr_rho, size = detect_effect), show.legend = TRUE, alpha = 0.5) +
#  geom_hline(yintercept = 0, linetype = "dashed")+
#  facet_grid( . ~ approach, space="free") +
#  theme(legend.position = 'none')+
#  ylab("Estimated effect size (z-score)") +
#  xlab("True effect size") +
#  theme(axis.text = element_text(size = 12),axis.text.x = element_text(angle=45),strip.text = element_text(size = 12))
#
#
#
#ggplot(tmp[tmp$approach == '2.Rates weigthed' & tmp$type == 'Exposure bias',], aes(x = tie_effect, y = `p.value`, group = approach))+
#  geom_point(aes(size = N_id,  color = detect_effect), show.legend = TRUE, position=position_jitter(0.2))+
#  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1)+
#  xlab("True effect size")+
#  geom_vline(xintercept = 0.25, linetype = "dashed")+
#  geom_vline(xintercept = -0.25, linetype = "dashed")+
#  facet_grid( . ~ approach, space="free")
#
#### False negatives----------------
#tmp = error[[1]]
#tmp4 = tmp[tmp$type %in% c("Differences in sociality, censoring bias", "Differences in sociality, exposure bias", "Differences in sociality, no biases"),]
#tmp4$type = ifelse(tmp4$type == "Differences in sociality, censoring bias", 'Censoring biases', tmp4$type)
#tmp4$type = ifelse(tmp4$type == "Differences in sociality, exposure bias", 'Exposure biases', tmp4$type)
#tmp4$type = ifelse(tmp4$type == "Differences in sociality, no biases", 'No biases', tmp4$type)
#
#tmp5 = tmp4[tmp4$type =='No biases',]
#ggplot(tmp4, aes(x = type, y = Percent, group = Approaches, colour  = Approaches, label = Approaches))+geom_point(aes(size = 2))+geom_line()+facet_grid(~ Error.type) +
#  geom_label_repel(data = tmp5, aes(x = type, y = Percent, colour = Approaches),
#                   nudge_x = 1,
#                   nudge_y = 8,
#                   segment.curvature = -1e-20,
#                   force = 10,
#                   direction = "y",
#                   hjust= 0)+
#  theme(legend.position = 'None', text = element_text(size=20))+xlab("Error type")+geom_hline(yintercept = 5, linetype = "dashed", color = 'red')
#
#tmp = error[[2]]
#tmp$approach = gsub("2.","",tmp$approach)
#tmp$approach = gsub("3.","",tmp$approach)
#
#ggplot(tmp[tmp$type == "No biases", ], aes(x = tie_effect,  y = z, group = sim, label = z))+
#  geom_point(aes(color = sr_rho, size = detect_effect), show.legend = TRUE, alpha = 0.5) +
#  geom_hline(yintercept = 0, linetype = "dashed")+
#  facet_grid( . ~ approach, space="free") +
#  theme(legend.position = 'none')+
#  ylab("Estimated effect size (z-score)") +
#  xlab("True effect size") +
#  theme(axis.text = element_text(size = 12),axis.text.x = element_text(angle=45),strip.text = element_text(size = 12))
#
#
#
#ggplot(tmp[tmp$approach == '2.Rates weigthed' & tmp$type == 'Exposure bias',], aes(x = tie_effect, y = `p.value`, group = approach))+
#  geom_point(aes(size = N_id,  color = detect_effect), show.legend = TRUE, position=position_jitter(0.2))+
#  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1)+
#  xlab("True effect size")+
#  geom_vline(xintercept = 0.25, linetype = "dashed")+
#  geom_vline(xintercept = -0.25, linetype = "dashed")+
#  facet_grid( . ~ approach, space="free")
#
######
#ggplot
#
#