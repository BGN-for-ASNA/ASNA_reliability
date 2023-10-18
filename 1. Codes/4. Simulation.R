###############################
##  Simulations parameters ####
###############################
###############################
simulations <- function(
    Reps = 100,
    N_id =  seq(30, 100, by = 20),
    hairy_tie_effect = seq(-4, 4, by = 1),
    hairy_detect_effect = seq(-4, 4, by = 1),
    blocks = TRUE,
    B = NULL,
    V = 1,
    G = 3,
    diagB = -6.5,
    groups=NULL,
    sr_sigma = c(1.4, 0.8),
    sr_rho = 0.5,
    dr_sigma = 1.2,
    dr_rho = 0.8,
    exposure_sigma = 2.9,
    exposure_baseline = 40,
    strand = T, # use STRAN model
    ncores = 1, # number of cores
    run_blocks_model = T
){
  require(parallel)
  require(foreach)
  require(doParallel)
  ########################################
  #######  Create grid of simulations ####
  ########################################
  grid = expand.grid(N_id, hairy_tie_effect, hairy_detect_effect)
  grid_subsample = grid[sample(1:nrow(grid), Reps, replace = FALSE),]
  colnames(grid_subsample) = c("N_id", "tie_effect", "detect_effect")

  #############################
  #######  Parallelisation ####
  #############################
  cl <- ncores
  cl <-makeCluster(cl, type="PSOCK")
  registerDoParallel(cores=cl)
  #############################
  #######  Testing methods ####
  #############################
  r <- foreach(i = 1:nrow(grid_subsample)) %dopar% {
    library(ANTs)
    library(STRAND)
    source("./1. Codes/3. Data simulation.R")
    RESULTS = NULL
    get_res = function(x) {
      y = c(coef(x)[2],confint(x, level = 0.9)[2,])
      return(y)
    }

    #"Bayesian P value".---------------
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

    # Variation --------------
    range_sr_sigma = seq(from = 0.5, to = 3, by = 0.2)
    range_sr_rho = seq(from = -0.9, to = 0.9, by = 0.1)

    range_dr_sigma = seq(from = 0.5, to = 3, by = 0.2)
    range_dr_rho = seq(from = -0.9, to = 0.9, by = 0.1)

    range_exposure_sigma  = seq(from = 0.5, to = 3, by = 0.2)
    range_exposure_baseline = seq(from = 10, to = 100, by = 20)

    picked_sr_sigma = sample(range_sr_sigma, 2)
    picked_sr_rho = sample(range_sr_rho, 1)

    picked_dr_sigma = sample(range_dr_sigma, 1)
    picked_dr_rho = sample(range_dr_rho, 1)

    picked_exposure_sigma  = sample(range_exposure_sigma, 1)
    picked_exposure_baseline  = sample(range_exposure_baseline, 1)


    # Make data--------------------------------
    N_id = grid_subsample$N_id[i]

    ## Block data -----------------
    NG = sample(c(1,3,7), 1) # Random number of groups
    clique = sample(1:G, N_id, replace = TRUE)
    mean.within.GR = sample(c(seq(from = -9, to = 9, by = 1)), 1) # Probability of random ties within a group.
    m = matrix(rnorm(length(m), mean.within.GR, sd = 1), NG, NG)

    mean.between.GR = sample(c(seq(from = -9, to = 0, by = 1)), 1) # Reduce randomly the probability of  ties between groups.
    diag(m) = diag(m) + rnorm(NG, mean.between.GR, sd = 1)


    block = data.frame(Clique=factor(clique))

    Clique = rep(1, N_id)
    Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)

    # Simulated data -----------------
    A = simulate_sbm_plus_srm_network_with_measurement_bias(
      N_id = N_id,
      B=list(B=B),
      V=V,
      groups=groups,
      individual_predictors=Hairy,
      individual_effects=matrix(c(grid_subsample$tie_effect[i], grid_subsample$tie_effect[i]),ncol=1, nrow=2),

      sr_sigma = picked_sr_sigma,
      sr_rho = picked_sr_rho,
      dr_sigma = picked_dr_sigma,
      dr_rho = picked_dr_rho,
      exposure_predictors = cbind(rep(1,N_id),Hairy),
      exposure_effects = c(-1, grid_subsample$hairy_detect_effect[i]),
      exposure_sigma = picked_exposure_sigma,
      exposure_baseline = picked_exposure_baseline,
      int_bias = T,
      return.network = FALSE,
    )

    # Zero-inflated Poisson model ----------------------
    y = A$interaction
    exposure = A$exposure
    library(rethinking)
    m <- map(
      alist(
        y ~ dzipois( p , lambda ),
        logit(p) <- ap, #Probability of true zero
        log(lambda) <- al + log(exposure+1) , #  logarithm of the exposure to address observation bias
        ap ~ dnorm(0,1),
        al ~ dnorm(0,10)
      ) ,data=list(y=y, exposure = exposure))

    r = precis(m)
    logistic(r$mean[1]) # probability false zero
    exp(r$mean[2]) # rate interaction

    # STRAND--------------------------------
    indiv =  data.frame(Hairy = Hairy)
    if(strand){
      nets = list(Grooming = A$network)
      exposure_nets = list(Exposure = A$true_samps)


      model_dat = make_strand_data(outcome = nets,
                                   individual_covariates = indiv,
                                   block_covariates = block,
                                   outcome_mode = "binomial",
                                   exposure = exposure_nets
      )

      # !!  I replaced fit_block_plus_social_relations_model by fit_social_relations_model
      if(run_blocks_model){
        fit =  fit_block_plus_social_relations_model(data=model_dat,
                                                     block_regression = ~ Clique,
                                                     focal_regression = ~ Hairy,
                                                     target_regression = ~ Hairy,
                                                     dyad_regression = ~  1,
                                                     mode="mcmc",
                                                     stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                                 iter_warmup = 1000, iter_sampling = 1000,
                                                                                 max_treedepth = NULL, adapt_delta = .98)
        )
      }else{
        fit =  fit_social_relations_model(data=model_dat,
                                          focal_regression = ~ Hairy,
                                          target_regression = ~ Hairy,
                                          dyad_regression = ~  1,
                                          mode="mcmc",
                                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                      iter_warmup = 1000, iter_sampling = 1000,
                                                                      max_treedepth = NULL, adapt_delta = .98)
        )

      }




      res = summarize_strand_results(fit)

      # "Bayesian P value"
      strand_est = P_se(res$sample$srm_model_samples$focal_coeffs[,1])

      result = as.data.frame(t(data.frame(as.numeric(unlist(c(unlist(res$summary[2,2:4]), grid_subsample[i,], strand_est))))))
      rownames(result) = NULL
      result$approach = 'strand'
      result$sim = i
      result$sr_sigma = paste(picked_sr_sigma[1], picked_sr_sigma[2], sep = "/")
      result$sr_rho  = picked_sr_rho
      result$dr_sigma = picked_dr_sigma
      result$dr_rho = picked_dr_rho
      result$exposure_sigma  = picked_exposure_sigma
      result$exposure_baseline = picked_exposure_baseline

      colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value',
                           'approach', 'sim', 'sr_sigma', 'sr_rho', 'dr_sigma','dr_rho','exposure_sigma','exposure_baseline')
      RESULTS = rbind(RESULTS, result)
    }

    #  Rates of interactions unweighted--------------------------------
    tie_strength = A$network/(1+A$true_samps)
    colnames(tie_strength) = rownames(tie_strength) = 1:ncol(tie_strength)
    df = ANTs:::df.create(tie_strength)
    df$strength = met.strength(tie_strength)
    df$hair = indiv$Hairy
    test1.1 = lm(strength ~ hair, data = df)

    ants_est = summary(test1.1)$coefficients[2,4]
    ants_se = summary(test1.1)$coefficients[2,2]
    result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,],ants_est, ants_se)))))
    rownames(result) = NULL
    result$approach = 'Rates unweighted'
    result$sim = i
    result$sr_sigma = paste(picked_sr_sigma[1], picked_sr_sigma[2], sep = "/")
    result$sr_rho  = picked_sr_rho
    result$dr_sigma = picked_dr_sigma
    result$dr_rho = picked_dr_rho
    result$exposure_sigma  = picked_exposure_sigma
    result$exposure_baseline = picked_exposure_baseline

    colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value', 'se',
                         'approach', 'sim', 'sr_sigma', 'sr_rho', 'dr_sigma','dr_rho','exposure_sigma','exposure_baseline')
    RESULTS = rbind(RESULTS, result)


    #  Rates of interactions weighted--------------------------------
    test1.1 = lm(strength  ~ hair, data = df, weights = A$true_exposure)
    ants_est = summary(test1.1)$coefficients[2,4]
    ants_se = summary(test1.1)$coefficients[2,2]

    result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,],ants_est, ants_se)))))
    rownames(result) = NULL
    result$approach = 'Rates weighted'
    result$sim = i
    result$sr_sigma = paste(picked_sr_sigma[1], picked_sr_sigma[2], sep = "/")
    result$sr_rho  = picked_sr_rho
    result$dr_sigma = picked_dr_sigma
    result$dr_rho = picked_dr_rho
    result$exposure_sigma  = picked_exposure_sigma
    result$exposure_baseline = picked_exposure_baseline

    colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value', 'se',
                         'approach', 'sim', 'sr_sigma', 'sr_rho', 'dr_sigma','dr_rho','exposure_sigma','exposure_baseline')
    RESULTS = rbind(RESULTS, result)

    #  SRI of interactions unweighted--------------------------------
    m = A$network
    fa = fb =matrix(0, ncol = nrow(m), nrow = nrow(m))
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
    ants_est = summary(test1.1)$coefficients[2,4]
    ants_se = summary(test1.1)$coefficients[2,2]

    result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,],ants_est, ants_se)))))
    rownames(result) = NULL
    result$approach = 'SRI unweigthed'
    result$sim = i
    result$sr_sigma = paste(picked_sr_sigma[1], picked_sr_sigma[2], sep = "/")
    result$sr_rho  = picked_sr_rho
    result$dr_sigma = picked_dr_sigma
    result$dr_rho = picked_dr_rho
    result$exposure_sigma  = picked_exposure_sigma
    result$exposure_baseline = picked_exposure_baseline

    colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value', 'se',
                         'approach', 'sim', 'sr_sigma', 'sr_rho', 'dr_sigma','dr_rho','exposure_sigma','exposure_baseline')
    RESULTS = rbind(RESULTS, result)

    #  SRI of interactions weighted--------------------------------
    test1.1 = lm(strength ~ hair, data = df, weights = A$true_exposure)
    ants_est = summary(test1.1)$coefficients[2,4]
    ants_se = summary(test1.1)$coefficients[2,2]

    result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,],ants_est, ants_se)))))
    rownames(result) = NULL
    result$approach = 'SRI weighted'
    result$sim = i
    result$sr_sigma = paste(picked_sr_sigma[1], picked_sr_sigma[2], sep = "/")
    result$sr_rho  = picked_sr_rho
    result$dr_sigma = picked_dr_sigma
    result$dr_rho = picked_dr_rho
    result$exposure_sigma  = picked_exposure_sigma
    result$exposure_baseline = picked_exposure_baseline

    colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value', 'se',
                         'approach', 'sim', 'sr_sigma', 'sr_rho', 'dr_sigma','dr_rho','exposure_sigma','exposure_baseline')
    RESULTS = rbind(RESULTS, result)
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

  result$approach = ifelse(result$approach == "strand", '1.Bayesian', result$approach)
  result$approach = ifelse(result$approach == "Rates unweighted", '2.Rates', result$approach)
  result$approach = ifelse(result$approach == "SRI unweigthed", '3.SRI', result$approach)

  return(result)
}


# Plot ----------
plots <- function(result){
  require(ggplot2)
  #p1 = ggplot(result, aes(x= z,  y= tie_effect, color = sim, group = sim, label = z))+
  #  #geom_linerange(aes(xmin=result[,2], xmax=result[,3])) +
  #  geom_point(aes(color = sim, size = 1), show.legend = FALSE, alpha = 0.5) +
  #  facet_grid( . ~ approach, space="free") +
  #  theme(legend.position = 'none')+
  #  xlim(min(result$z), max(result$z)) +
  #  ylim(min(result$z), max(result$z)) +
  #  geom_abline()+
  #  ylab("True efect size") +
  #  xlab("Estimated effect size") +
  #  theme(axis.text = element_text(size = 14))  +
  #  theme(axis.title = element_text(size = 14))

  p1 = ggplot(result, aes(x= z,  y= tie_effect))+
    geom_point(aes(color = sim, size = 1), show.legend = FALSE, alpha = 0.5) +
    geom_linerange(aes(xmin=ci5, xmax=ci95)) +
    facet_grid( . ~ approach, space="free") +
    theme(legend.position = 'none')+
    xlim(min(result$ci5), max(result$ci95)) +
    ylim(min(result$ci5), max(result$ci95)) +
    geom_abline()+
    ylab("True efect size") +
    xlab("Estimated effect size") +
    theme(axis.text = element_text(size = 14))  +
    theme(axis.title = element_text(size = 14))



  p2 = ggplot(result, aes(x = resid, y = tie_effect, color = approach))+
    geom_point(aes(size = N_id),  alpha = 0.5, show.legend = FALSE)+
    geom_vline(xintercept = 0, linetype = "dashed")+
    xlab("Difference between true effect size and estimated effect size")+
    ylab("True effect size")+
    theme(axis.text = element_text(size = 14))+
    theme(axis.title = element_text(size = 14))+
    theme(legend.position = 'none')

  p3 = ggplot(result, aes(x= approach, y = resid, color = approach, group = approach))+
    geom_violin()+
    stat_summary(fun= function(i){r = quantile(i)[2:4]},  shape=23, size=10, geom = 'point', color = 'black')+
    geom_jitter(aes(size = 1), alpha = 0.5, show.legend = FALSE, position=position_jitter(0.2)) +
    theme(axis.text = element_text(size = 14))+
    theme(axis.title = element_text(size = 14))+
    theme(legend.position = 'none')
  return(list(p1, p2, p3))
}


library(ggpubr)
result = simulations(Reps = 10, strand = F, ncores = 1)
tmp = result[result$approach %in% c('1.Bayesian', '2.Rates', '3.SRI'),]
p = plots(result)
ggarrange(p[[1]], p[[3]], ncol = 1, nrow = 2)

# Analysis ----------
library(lmerTest)
library(emmeans)
m = lmer(formula = resid ~ approach + (1|tie_effect), data = tmp)
summary(m)

emmeans(m, list(pairwise ~ approach), adjust = "tukey")
