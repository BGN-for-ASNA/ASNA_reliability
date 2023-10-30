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
    run_blocks_model = TRUE
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
        strand_se = P_se(test1.1$sample$srm_model_samples$focal_coeffs[,1])

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
        result$NG = NG
        result$mean.between.GR = mean.between.GR
        result$mean.within.GR = mean.within.GR

        colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value', 'se',
                             'approach', 'sim', 'sr_sigma', 'sr_rho', 'dr_sigma','dr_rho','exposure_sigma','exposure_baseline',
                             'NG', 'mean.between.GR', 'mean.within.GR')
        return(result)
      }else{
        ants_est = summary(test1.1)$coefficients[2,4]
        ants_se = summary(test1.1)$coefficients[2,2]
        result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,],ants_est, ants_se)))))
        rownames(result) = NULL
        result$approach = name
        result$sim = i
        result$sr_sigma = paste(picked_sr_sigma[1], picked_sr_sigma[2], sep = "/")
        result$sr_rho  = picked_sr_rho
        result$dr_sigma = picked_dr_sigma
        result$dr_rho = picked_dr_rho
        result$exposure_sigma  = picked_exposure_sigma
        result$exposure_baseline = picked_exposure_baseline
        result$NG = NG
        result$mean.between.GR = mean.between.GR
        result$mean.within.GR = mean.within.GR

        colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value', 'se',
                             'approach', 'sim', 'sr_sigma', 'sr_rho', 'dr_sigma','dr_rho','exposure_sigma','exposure_baseline',
                             'NG', 'mean.between.GR', 'mean.within.GR')
        return(result)
      }
    }

    # Make data--------------------------------
    N_id = grid_subsample$N_id[i]

    ## Variation in sr ad dr --------------
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

    ## Block data -----------------
    NG = sample(c(1,3,7), 1) # Random number of groups
    clique = sample(1:G, N_id, replace = TRUE)

    mean.within.GR = sample(c(seq(from = -9, to = 9, by = 1)), 1) # Probability of random ties within a group.
    B = matrix(rnorm(NG*NG, mean.within.GR, sd = 1), NG, NG)

    mean.between.GR = sample(c(seq(from = 0, to = 9, by = 1)), 1) # Increase randomly the probability of  ties within groups.
    diag(B) = diag(B) + rnorm(NG, mean.between.GR, sd = 1)

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
      return.network = TRUE,
    )
    indiv =  data.frame(Hairy = Hairy)
    # Zero-inflated Poisson model ----------------------
    #y = A$interaction
    #exposure = A$exposure
    #library(rethinking)
    #m <- map(
    #  alist(
    #    y ~ dzipois( p , lambda ),
    #    logit(p) <- ap, #Probability of true zero
    #    log(lambda) <- al + log(exposure+1) , #  logarithm of the exposure to address observation bias
    #    ap ~ dnorm(0,1),
    #    al ~ dnorm(0,10)
    #  ) ,data=list(y=y, exposure = exposure))
    #r = precis(m)
    #logistic(r$mean[1]) # probability false zero
    #exp(r$mean[2]) # rate interaction

    # STRAND--------------------------------
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
      result = get.result(res, strand = TRUE)
      RESULTS = rbind(RESULTS, result)
    }

    # Rates of interactions unweighted--------------------------------
    tie_strength = A$network/(1+A$true_samps)
    colnames(tie_strength) = rownames(tie_strength) = 1:ncol(tie_strength)
    df = ANTs:::df.create(tie_strength)
    df$strength = met.strength(tie_strength)
    df$hair = indiv$Hairy
    test1.1 = lm(strength ~ hair, data = df)
    result = get.result(test1.1, strand = FALSE, name = '2.Rates')
    RESULTS = rbind(RESULTS, result)

    # Rates of interactions weighted--------------------------------
    test1.1 = lm(strength  ~ hair, data = df, weights = A$true_exposure)
    result = get.result(test1.1, strand = FALSE, name = '2.Rates weigthed')
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

    # SRI of interactions weighted--------------------------------
    test1.1 = lm(strength ~ hair, data = df, weights = A$true_exposure)
    result = get.result(test1.1, strand = FALSE, name = '3.SRI weigthed')
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


  return(result)
}

# Plot ----------
plots <- function(result){
  p1 = ggplot(result, aes(x = tie_effect,  y = z, color = sim, group = sim, label = z))+
    #geom_linerange(aes(xmin=result[,2], xmax=result[,3])) +
    geom_point(aes(color = sim, size = sr_rho), show.legend = FALSE, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    facet_grid( . ~ approach, space="free") +
    theme(legend.position = 'none')+
    ylab("Estimated effect size (z-score)") +
    xlab("True efect size") +
    theme(axis.text = element_text(size = 14))

  #p1 = ggplot(result, aes(x= z,  y= tie_effect))+
  #  geom_point(aes(color = sim, size = 1), show.legend = FALSE, alpha = 0.5) +
  #  geom_linerange(aes(xmin=ci5, xmax=ci95)) +
  #  facet_grid( . ~ approach, space="free") +
  #  theme(legend.position = 'none')+
  #  xlim(min(result$ci5), max(result$ci95)) +
  #  ylim(min(result$ci5), max(result$ci95)) +
  #  geom_abline()+
  #  ylab("True efect size") +
  #  xlab("Estimated effect size") +
  #  theme(axis.text = element_text(size = 14))  +
  #  theme(axis.title = element_text(size = 14))


 p2 = ggplot(result, aes(x = tie_effect, y = z, color = approach))+
   geom_point(aes(size = sr_rho),  show.legend = FALSE, alpha = 0.5)+
   geom_hline(yintercept = 0, linetype = "dashed")+
   xlab("True efect size")+
   ylab("Estimated effect size (z-score)")+
   theme(axis.text = element_text(size = 14))
   #theme(legend.position = 'none')+
   #facet_grid( . ~ detect_effect, space="free")

 #p3 = ggplot(result, aes(x= approach, y = resid, color = approach, group = approach))+
 #  geom_violin()+
 #  stat_summary(fun= function(i){r = quantile(i)[2:4]},  shape=23, size=10, geom = 'point', color = 'black')+
 #  geom_jitter(aes(size = sr_rho), alpha = 0.5, show.legend = FALSE, position=position_jitter(0.2)) +
 #  theme(axis.text = element_text(size = 14))+
 #  theme(axis.title = element_text(size = 14))+
 #  theme(legend.position = 'none')

 p4 = ggplot(result, aes(x = tie_effect, y = `p-value`, group = approach))+
   geom_point(aes(size = sr_rho,  color = approach), alpha = 0.5, show.legend = TRUE, position=position_jitter(0.2))+
   geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1)+
   xlab("True efect size")+
   geom_vline(xintercept = 0.5, linetype = "dashed")+
   geom_vline(xintercept = -0.5, linetype = "dashed")
  return(list(p1, p2, p4))
}

library(ggpubr)
result = simulations(Reps = 1, strand = F, ncores = 4)
p = plots(result)
ggarrange(p[[2]], p[[3]], ncol = 2, nrow = 1, common.legend = TRUE)

# Rates of false negatives -------------
t1 = sum(result[result$tie_effect > 0.5 | result$tie_effect < -0.5,]$`p-value` >= 0.05)/nrow(result[result$tie_effect > 0.5 | result$tie_effect < -0.5,])
t2 = tapply(result[result$tie_effect > 0.5 | result$tie_effect < -0.5,]$`p-value`, result[result$tie_effect > 0.5 | result$tie_effect < -0.5,]$approach, function(x){sum(x >= 0.05)/length(x)})

# Rates of false positives -------------
t3 = sum(result[result$tie_effect <= 0.5 & result$tie_effect >= -0.5,]$`p-value` <= 0.05)/nrow(result[result$tie_effect <= 0.5 & result$tie_effect >= -0.5,])
t4 = tapply(result[result$tie_effect <= 0.5 &  result$tie_effect  >= -0.5,]$`p-value`, result[result$tie_effect<= 0.5 &  result$tie_effect  >= -0.5,]$approach, function(x){sum(x <= 0.05)/length(x)})


summary = data.frame("Approaches" = c('all', names(t2), "all", names(t4)), "Error type" = c(rep('False negatives', 6), rep('False positives', 6)) ,"Percent" = c(t1, t2, t3, t4))
summary$Percent = summary$Percent * 100
summary

write.csv(summary, "2. Results/ Rates of Type I and Type II errors.csv", row.names = FALSE)

