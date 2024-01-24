###############################
##  Simulations parameters ####
###############################
###############################
## Individuals parameters
#Reps = 100
#N_id =  seq(30, 50, by = 20)
#hairy_tie_effect = seq(-2, 2, by = 0.5)
#hairy_detect_effect = seq(-4, 4, by = 1)
## Block parameter
#blocks = TRUE
#B = NULL
#V = 1
#G = 3
#diagB = -6.5
#groups=NULL
## Sender-receiver parameter
#sr_sigma = c(1.4, 0.8)
#sr_rho = 0.5
## Dyadic parameter
#dr_sigma = 1.2
#dr_rho = 0.8
## Observation bias parameter
#exposure.bias = TRUE
#exposure_sigma = 2.9
#exposure_baseline = 40
## Interactions bias parameter
#simulate.interactions = TRUE
#int_intercept =c(Inf,Inf) #invert log of inf = 1 of prob to observe interaction for both focal and alte
#int_slope = c(-Inf,-Inf) # No effect of individuals attributes on missing interactio
## Simulation parameter
#BISON = TRUE
#STRAND = TRUE # use STRAN mode
#ncores = 1 # number of core
#blockModel = TRUE

simulations <- function(
    # Individuals parameters
    Reps = 100,
    N_id =  seq(30, 50, by = 20),
    hairy_tie_effect = seq(-2, 2, by = 0.5),
    hairy_detect_effect = seq(-4, 4, by = 1),
    
    # Block parameters
    blocks = TRUE,
    B = NULL,
    V = 1,
    G = 3,
    diagB = -6.5,
    groups=NULL,
    
    # Sender-receiver parameters
    sr_sigma = c(1.4, 0.8),
    sr_rho = 0.5,
    
    # Dyadic parameters
    dr_sigma = 1.2,
    dr_rho = 0.8,
    
    # Observation bias parameters
    exposure.bias = TRUE,
    exposure_sigma = 2.9,
    exposure_baseline = 40,
    
    # Interactions bias parameters
    simulate.interactions = TRUE,
    int_intercept =c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
    int_slope = c(-Inf,-Inf), # No effect of individuals attributes on missing interaction
    
    # Simulation parameters
    BISON = TRUE,
    STRAND = TRUE, # use STRAN model
    ncores = 1, # number of cores
    blockModel = TRUE
){
  #require(parallel)
  #require(doParallel)
  ########################################
  #######  Create grid of simulations ####
  ########################################
  grid = expand.grid(N_id, hairy_tie_effect, hairy_detect_effect)
  if(nrow(grid) <= Reps){
    grid_subsample = grid
    ncores = nrow(grid)
  }else{
    grid_subsample = grid[sample(1:nrow(grid), Reps, replace = FALSE),]
  }
  colnames(grid_subsample) = c("N_id", "tie_effect", "detect_effect")

  #############################
  #######  Parallelisation ####
  #############################
  require(foreach)
  library(doSNOW)
  cl <- ncores
  
  cl <- makeCluster(cl, outfile="")
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = nrow(grid_subsample), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  #cl <-makeCluster(cl, type="PSOCK")
  #registerDoParallel(cores=cl)
  
  on.exit(registerDoSEQ())
  #############################
  #######  Testing methods ####
  #############################
  r <- foreach(i = 1:nrow(grid_subsample),.options.snow = opts) %dopar% {
    print(i)
    library(ANTs)
    library(STRAND)
    source("./1.Codes/2.data_simulation.R")
    RESULTS = NULL

    # pre-network permutation-----------
    preNetPerm = function (df, focal = 'focal', alters = 'receiver', ctrl = 'focal', nperm = 10, progress = T, index = 'sri') 
    {
      col.alters <- ANTs:::df.col.findId(df, alters)
      col.focal <- ANTs:::df.col.findId(df, focal)
      col.ctrl <- ANTs:::df.col.findId(df, ctrl)
      
      ctrl <- c(col.ctrl, col.focal)
      df <- df.ctrlFactor(df, control = ctrl)
      df$control <- as.factor(df$control)
      Vecids <- unique(c(as.character(df[, col.focal]), as.character(df[, col.alters])))
      group_scan <- unique(df$control)
      
      GBI2 <- ANTs:::df_to_gbi(df, ncol(df), col.focal, Vecids, group_scan)
      GBI <- ANTs:::df_to_gbi(df, ncol(df), col.alters, Vecids, group_scan) 

      result <- ANTs:::perm_dataStream1_focal(GBI, GBI2, nperm = nperm, 
                                       progress = progress, method = index)
      result <- lapply(seq_along(result), function(x, Vecids, i) {
        colnames(x[[i]]) <- Vecids
        rownames(x[[i]]) <- Vecids
        attr(x[[i]], "permutation") <- i
        return(x[[i]])
      }, x = result, Vecids = Vecids)
      return(result)
    }
    
    # "Bayesian P value".---------------
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

    # Build table-------------
    get_res = function(x) {
      y = c(coef(x)[2],confint(x, level = 0.9)[2,])
      return(y)
    }
    get.result <- function(test1.1, strand = TRUE, name){
      if(strand){
        # "Bayesian P value"
        strand_est = P_se(test1.1$sample$srm_model_samples$focal_coeffs[,1])
        a = test1.1$sample$srm_model_samples$focal_coeffs[,1]
        strand_se = sqrt(sum((a-mean(a))^2/(length(a)-1)))/sqrt(length(a))

        result = as.data.frame(t(data.frame(as.numeric(unlist(c(unlist(test1.1$summary[2,2:4]), grid_subsample[i,], strand_est, strand_se))))))
        rownames(result) = NULL
        result$approach = 'strand'
        result$sim = i
        result$sr_sigma = paste(picked_sr_sigma[1], picked_sr_sigma[2], sep = "/")
        result$sr_rho  = picked_sr_rho
        result$dr_sigma = picked_dr_sigma
        result$dr_rho = picked_dr_rho
        result$exposure_sigma  = picked_exposure_sigma
        result$exposure_baseline = picked_exposure_baseline
        result$NG = ifelse(is.null(NG), NA, NG)
        result$mean.between.GR = ifelse(is.null(mean.between.GR), NA, mean.between.GR)
        result$mean.within.GR = ifelse(is.null(mean.within.GR), NA, mean.within.GR)

        colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value', 'se',
                             'approach', 'sim', 'sr_sigma', 'sr_rho', 'dr_sigma','dr_rho','exposure_sigma','exposure_baseline',
                             'NG', 'mean.between.GR', 'mean.within.GR')
        return(result)
      }else{
        ants_p = summary(test1.1)$coefficients[2,4]
        ants_se = summary(test1.1)$coefficients[2,2]
        result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,],ants_p, ants_se)))))
        rownames(result) = NULL
        result$approach = name
        result$sim = i
        result$sr_sigma = paste(picked_sr_sigma[1], picked_sr_sigma[2], sep = "/")
        result$sr_rho  = picked_sr_rho
        result$dr_sigma = picked_dr_sigma
        result$dr_rho = picked_dr_rho
        result$exposure_sigma  = picked_exposure_sigma
        result$exposure_baseline = picked_exposure_baseline
        result$NG = ifelse(is.null(NG), NA, NG)
        result$mean.between.GR = ifelse(is.null(mean.between.GR), NA, mean.between.GR)
        result$mean.within.GR = ifelse(is.null(mean.within.GR), NA, mean.within.GR)

        colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value', 'se',
                             'approach', 'sim', 'sr_sigma', 'sr_rho', 'dr_sigma','dr_rho','exposure_sigma','exposure_baseline',
                             'NG', 'mean.between.GR', 'mean.within.GR')
        return(result)
      }
    }

    # Make data--------------------------------
    N_id = grid_subsample$N_id[i]
    Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
    indiv =  data.frame(Hairy = Hairy)
    ## Variation in sr ad dr --------------
    if(is.null(sr_sigma)){
      range_sr_sigma = seq(from = 0.5, to = 3, by = 0.2)
      picked_sr_sigma = sample(range_sr_sigma, 2)
    }else{
      picked_sr_sigma = sr_sigma
    }
    
    if(is.null(sr_rho)){
      range_sr_rho = seq(from = -0.9, to = 0.9, by = 0.1)
      picked_sr_rho = sample(range_sr_rho, 1)
    }else{
      picked_sr_rho = sr_rho
    }
    
    if(is.null(dr_sigma)){
      range_dr_sigma = seq(from = 0.5, to = 3, by = 0.2)
      picked_dr_sigma = sample(range_dr_sigma, 1)
    }else{
      picked_dr_sigma = dr_sigma
    }
    
    if(is.null(dr_rho)){
      range_dr_rho = seq(from = -0.9, to = 0.9, by = 0.1)
      picked_dr_rho = sample(range_dr_rho, 1)
    }else{
      picked_dr_rho = dr_rho
    }
    
    if(is.null(exposure_sigma)){
      range_exposure_sigma  = seq(from = 0.5, to = 3, by = 0.2)
      picked_exposure_sigma = sample(range_dr_rho, 1)
    }else{
      picked_exposure_sigma = exposure_sigma
    }
    
    if(is.null(exposure_baseline)){
      range_exposure_baseline  = seq(from = 0.5, to = 3, by = 0.2)
      picked_exposure_baseline  = sample(range_exposure_baseline, 1)
    }else{
      picked_exposure_baseline = exposure_baseline
    }
    
    if(exposure.bias){
      sample.percent = sample(seq(from = 0.05, to = 0.5, by = 0.01),1)
      grid_subsample$detect_effect[i] =  grid_subsample$tie_effect[i]*sample.percent
      exposure_predictors = cbind(rep(1,N_id),Hairy)
    }else{
      grid_subsample$detect_effect[i] = 0
      exposure_predictors = NULL
    }

    ## Block data -----------------
    if(blockModel){
      NG = sample(c(1,3,7), 1) # Random number of groups
      clique = sample(1:NG, N_id, replace = TRUE)

      mean.within.GR = sample(c(seq(from = -9, to = 9, by = 1)), 1) # Probability of random ties within a group.
      B = matrix(rnorm(NG*NG, mean.within.GR, sd = 1), NG, NG)

      mean.between.GR = sample(c(seq(from = 0, to = 9, by = 1)), 1) # Increase randomly the probability of  ties within groups.
      diag(B) = diag(B) + rnorm(NG, mean.between.GR, sd = 1)

      block = data.frame(Clique=factor(clique))

    }else{
      NG = V = 1
      B = mean.within.GR = mean.between.GR = block = NULL
    }

    # Simulated data -----------------
    A = simulate_sbm_plus_srm_network_with_measurement_bias(
      N_id = N_id,
      B=list(B=B),
      V=V,
      groups=block,
      individual_predictors=Hairy,
      individual_effects=matrix(c(grid_subsample$tie_effect[i], grid_subsample$tie_effect[i]),ncol=1, nrow=2),

      sr_sigma = picked_sr_sigma,
      sr_rho = picked_sr_rho,
      dr_sigma = picked_dr_sigma,
      dr_rho = picked_dr_rho,
      exposure_predictors = cbind(rep(1,N_id),Hairy),
      exposure_effects = c(-1, grid_subsample$detect_effect[i]),
      exposure_sigma = picked_exposure_sigma,
      exposure_baseline = picked_exposure_baseline,
      simulate.interactions = simulate.interactions,
      int_intercept = int_intercept, # Prob of miss interactions
      int_slope = int_slope
    )
    A$network[is.na(A$network)] = 0
    # BISON------------------------
    if(simulate.interactions & BISON){
      library(bisonR)
      A$interactions$ego = as.factor(A$interactions$sender)
      A$interactions$alter = as.factor(A$interactions$receiver)
      priors <- get_default_priors("binary")
      fit_edge <- bison_model(
        (interaction | exposure) ~ dyad(ego, alter),
        data=A$interaction,
        model_type="binary",
        priors=priors
      )

      cv_samples <- extract_metric(fit_edge, "node_strength", num_draws = 1)
      df = data.frame('id' = names(fit_edge$node_to_idx), 'strength' = cv_samples[1,])
      df = df[order(as.numeric(df$id)),]
      df$hair = indiv$Hairy
      test1.1 = lm(strength ~ hair, data = df)
      result = get.result(test1.1, strand = FALSE, name = 'BISON')
      RESULTS = rbind(RESULTS, result)
    }

    # STRAND--------------------------------
    if(STRAND){
      
      nets = list(Grooming = A$network)
      
      if(exposure.bias){
        exposure_nets = list(Exposure = A$true_samps)
      }else{
        exposure_nets = NULL
      }


      model_dat = make_strand_data(outcome = nets,
                                   individual_covariates = indiv,
                                   block_covariates = block,
                                   outcome_mode = "binomial",
                                   exposure = exposure_nets
      )

      if(NG == 1){
        fit =  fit_social_relations_model(data=model_dat,
                                          focal_regression = ~ Hairy,
                                          target_regression = ~ Hairy,
                                          dyad_regression = ~  1,
                                          mode="mcmc",return_predicted_network = TRUE,
                                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                      iter_warmup = 1000, iter_sampling = 1000,
                                                                      max_treedepth = NULL, adapt_delta = .98))
      }else{
        fit =  fit_block_plus_social_relations_model(data=model_dat,
                                                     block_regression = ~ Clique,
                                                     focal_regression = ~ Hairy,
                                                     target_regression = ~ Hairy,
                                                     dyad_regression = ~  1,
                                                     mode="mcmc",return_predicted_network = TRUE,
                                                     stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                                 iter_warmup = 1000, iter_sampling = 1000,
                                                                                 max_treedepth = NULL, adapt_delta = .98))
      }
      
      res = summarize_strand_results(fit)
      result = get.result(res, strand = TRUE)
      RESULTS = rbind(RESULTS, result)
      
      # Compare estimated network
      est_net <- round(apply(res$samples$predicted_network_sample,2:3,mean ))
      m = matrix(0, ncol = N_id, nrow = N_id)
      for(a in 1:N_id){
        x = res$samples$predicted_network_sample[,,a]
        m[,a] = apply(x,2,mean)
      }
      
      plot(m,A$network/(1+A$true_samps))
    }

    # Rates of interactions unweighted--------------------------------
    tie_strength = A$network/(1+A$true_samps)
    colnames(tie_strength) = rownames(tie_strength) = 1:ncol(tie_strength)
    df = ANTs:::df.create(tie_strength)
    df$strength = met.strength(tie_strength)
    df$hair = indiv$Hairy
    test1.1 = lm(strength ~ hair, data = df)
    result = get.result(test1.1, strand = FALSE, name = '2.Rates')
    result
    #ggplot(df, aes(x = hair, y = strength))+geom_point()
    RESULTS = rbind(RESULTS, result)

    # Node label permutations
    perm = perm.net.nl(df, labels = 'hair', nperm = 10000, progress = FALSE)
    ptest = stat.lm(perm, formula = strength ~ hair, progress = FALSE)
    ptest = ANTs:::stat.p(c(ptest$Original.model$coefficients[2,1], ptest$permutations$hair))[3]
    result$approach = paste(result$approach, 'permuted', sep = " ")
    result$`p-value` = ptest
    RESULTS = rbind(RESULTS, result)

    # Rates of interactions weighted--------------------------------
    test1.1 = lm(strength  ~ hair, data = df, weights = A$true_exposure)
    result = get.result(test1.1, strand = FALSE, name = '2.Rates weigthed')
    RESULTS = rbind(RESULTS, result)

    # Node label permutations
    ptest = unlist(lapply(perm, function(x, exposure){
      m = lm(data = x, formula = strength ~ hair, weights = exposure)
      summary(m)$coefficients[2,1]
    }, exposure = A$true_exposure))
    ptest = ANTs:::stat.p(ptest)[3]
    result$approach = paste(result$approach, 'permuted', sep = " ")
    result$`p-value` = ptest
    RESULTS = rbind(RESULTS, result)

    # SRI of interactions unweighted--------------------------------
    m = A$network
    fa = fb = matrix(0, ncol = nrow(m), nrow = nrow(m))
    for (a in 1:nrow(fa)) {
      fa[a,] = A$true_exposure
      fb[,a] = A$true_exposure
    }
    ya = abs(fa - m)
    yb = abs(fb - m)

    sri <- ((m) /(m + ya + yb ))
    sri[is.na(sri)] = 0
    colnames(sri)= rownames(sri) = 1:ncol(sri)
    df = ANTs:::df.create(sri)
    df$strength = met.strength(sri)
    df$hair = indiv$Hairy
    test1.1 = lm(strength ~ hair, data = df, weights = A$true_exposure)
    result = get.result(test1.1, strand = FALSE, name = '3.SRI')
    RESULTS = rbind(RESULTS, result)

    # Node label permutations
    perm = perm.net.nl(df, labels = 'hair', nperm = 10000, progress = FALSE)
    ptest = stat.lm(perm, formula = strength ~ hair, progress = FALSE)
    ptest = ANTs:::stat.p(c(ptest$Original.model$coefficients[2,1], ptest$permutations$hair))[3]
    result$approach = paste(result$approach, 'permuted', sep = " ")
    result$`p-value` = ptest
    RESULTS = rbind(RESULTS, result)

    # SRI of interactions weighted--------------------------------
    test1.1 = lm(strength ~ hair, data = df, weights = A$true_exposure)
    result = get.result(test1.1, strand = FALSE, name = '3.SRI weigthed')
    RESULTS = rbind(RESULTS, result)

    # Node label permutations
    ptest = unlist(lapply(perm, function(x, exposure){
      m = lm(data = x, formula = strength ~ hair, weights = exposure)
      summary(m)$coefficients[2,1]
    }, exposure = A$true_exposure))
    ptest = ANTs:::stat.p(ptest)[3]
    result$approach = paste(result$approach, 'permuted', sep = " ")
    result$`p-value` = ptest
    RESULTS = rbind(RESULTS, result)
    
    ## SRI pre-network weighted--------------------------------
    #d = A$interactions
    #d = d[d$interaction == 1, ] #keeping only interactions
    #d$sender = as.character(d$sender)
    #d$receiver = as.character(d$receiver)
    ##d$ctrl = paste(d$follow, d$focal, sep = '_')
    ### Pre-net permutation
    #perms =preNetPerm(d,focal = 'follow', alters='receiver', ctrl = 'focal', nperm=1000, progress=F, index='sri')
    #tmp = perms[[1]]
    #r = met.strength(x,df,1)
    #r$id = as.numeric(r$id)
    #r = r[order(r$id),]
    #r$sampling = A$true_exposure
    #r$hair = indiv$Hairy
    #test = lm(strength ~ hair, data = r, weights = r$sampling)
    #result = get.result(test, strand = FALSE, name = 'SRI weigthed preNet')
    #
    #coefs = unlist(lapply(perms, function(x){
    #  r = met.strength(x,df,1)
    #  r$id = as.numeric(r$id)
    #  r = r[order(r$id),]
    #  r$sampling = A$true_exposure
    #  r$hair = indiv$Hairy
    #  test = coef(lm(strength ~ hair, data = r, weights = r$sampling))[2]
    #}))
    #ptest = ANTs:::stat.p(coefs)[3]
    #result$`p-value` = ptest
    #RESULTS = rbind(RESULTS, result)
    #
    ## Node label permutations
    #ptest = unlist(lapply(perm, function(x, exposure){
    #  m = lm(data = x, formula = strength ~ hair, weights = exposure)
    #  summary(m)$coefficients[2,1]
    #}, exposure = A$true_exposure))
    #ptest = ANTs:::stat.p(ptest)[3]
    #result$approach = paste(result$approach, 'permuted', sep = " ")
    #result$`p-value` = ptest
    #RESULTS = rbind(RESULTS, result)
  }

  #############################
  #######  Return #############
  #############################
  parallel::stopCluster(cl)
  gc()
  registerDoSEQ()
  result = do.call('rbind', r)
  result$z = result$`scale(hair)`/(result$`95 %` - result$`5 %`)
  result$MOE = result$z * result$se
  result$ci5 = result$z - abs(result$MOE)
  result$ci95 = result$z + abs(result$MOE)
  result$resid = result$tie_effect - result$z
  return(result)
}

# Plot ----------
plots <- function(result){
  require(ggplot2)
  # Removing permuted approaches, as the effect sizes are identical.
  p1 = ggplot(result[!result$approach %in% c("2.Rates permuted", "3.SRI permuted",'3.SRI weigthed permuted', "2.Rates weigthed permuted"),], aes(x = tie_effect,  y = z, group = sim, label = z))+
    #geom_linerange(aes(xmin=result[,2], xmax=result[,3])) +
    geom_point(aes(color = sr_rho, size = detect_effect), show.legend = TRUE, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    facet_grid( . ~ approach, space="free") +
    theme(legend.position = 'none')+
    ylab("Estimated effect size (z-score)") +
    xlab("True effect size") +
    theme(axis.text = element_text(size = 14))

 p4 = ggplot(result, aes(x = tie_effect, y = `p-value`, group = approach))+
   geom_point(aes(size = detect_effect,  color = sr_rho), alpha = 0.5, show.legend = TRUE, position=position_jitter(0.2))+
   geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1)+
   xlab("True effect size")+
   geom_vline(xintercept = 0.5, linetype = "dashed")+
   geom_vline(xintercept = -0.5, linetype = "dashed")+
   facet_grid( . ~ approach, space="free")
  return(list(p1, p4))
}




