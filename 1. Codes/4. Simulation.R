###############################
##  Simulations parameters ####
###############################
###############################
simulations <- function(
    Reps = 100,
    N_id =  seq(30, 50, by = 20),
    hairy_tie_effect = seq(-2, 2, by = 0.5),
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
    simulate.interactions = TRUE,
    int_intercept =c(5,5), # Prob of miss interactions
    int_slope = c(10,5)
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
      exposure_effects = c(-1, grid_subsample$detect_effect[i]),
      exposure_sigma = picked_exposure_sigma,
      exposure_baseline = picked_exposure_baseline,
      simulate.interactions = simulate.interactions,
      int_intercept = int_intercept, # Prob of miss interactions
      int_slope =int_slope
    )

    indiv =  data.frame(Hairy = Hairy)
    ## Zero-inflated Poisson model ----------------------
    #y = A$interaction$interaction
    #exposure = A$interaction$exposure
    ## Due to some individuals not being observed we will
    #x = unique(A$interactions$ego)
    #id = as.numeric(as.character(
    #    factor(A$interactions$ego, levels = unique(A$interactions$ego), labels = seq_along(unique(A$interactions$ego)))
    #  ))
    #unique(id)

    #library(rethinking)
    #m <- ulam (
    #  alist(
    #    y ~ dzipois( p , lambda ),
    #    logit(p) <- ap + a_id[id], #Probability of true zero
    #    log(lambda) <- al + a_id[id] + log(exposure+1) , #  logarithm of the exposure to address observation bias

    #    ap ~ dnorm(0,1),
    #    al ~ dnorm(0,10),

    #    a_id[id] ~ dnorm(0, sigma_id ),
    #    sigma_id ~ dcauchy(0,1)
    #  ) , data=list(y=y, exposure = exposure, id = id))

    #r = precis(m, depth=2)
    #logistic(r$mean[1]) # probability false zero
    #exp(r$mean[2]) # rate interaction

    #link.m <- link( m , data=list(y=y, exposure = exposure, id = id) )
    #pred.p <- apply( link.m , 2 , mean )
    #pred.p.PI <- apply( link.m , 2 , PI )

    #vcov(m) # variance-covariance
    #post <- extract.samples( m , n = N_id )
    #head(post)
    ### Generating implied observations -------------------
    # BISON------------------------
    if(simulate.interactions){
      library(bisonR)
      A$interactions$ego = as.factor(A$interactions$ego)
      A$interactions$alter = as.factor(A$interactions$alter)
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
    if(strand){
      nets = list(Grooming = A$network)
      exposure_nets = list(Exposure = A$true_samps)


      model_dat = make_strand_data(outcome = nets,
                                   individual_covariates = indiv,
                                   block_covariates = block,
                                   outcome_mode = "binomial",
                                   exposure = exposure_nets
      )

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
  p1 = ggplot(result[!result$approach %in% c("2.Rates permuted", "3.SRI permuted"),], aes(x = tie_effect,  y = z, group = sim, label = z))+
    #geom_linerange(aes(xmin=result[,2], xmax=result[,3])) +
    geom_point(aes(color = sr_rho, size = detect_effect), show.legend = TRUE, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed")+
    facet_grid( . ~ approach, space="free") +
    theme(legend.position = 'none')+
    ylab("Estimated effect size (z-score)") +
    xlab("True efect size") +
    theme(axis.text = element_text(size = 14))

 p4 = ggplot(result, aes(x = tie_effect, y = `p-value`, group = approach))+
   geom_point(aes(size = detect_effect,  color = sr_rho), alpha = 0.5, show.legend = TRUE, position=position_jitter(0.2))+
   geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1)+
   xlab("True efect size")+
   geom_vline(xintercept = 0.5, linetype = "dashed")+
   geom_vline(xintercept = -0.5, linetype = "dashed")+
   facet_grid( . ~ approach, space="free")
  return(list(p1, p4))
}

# Simulation ----------
result = simulations(Reps = 1, ncores = 1, strand = T, simulate.interactions = T, int_p = c(Inf,Inf)) # Simulate interactions with the probability of observing an interaction of 1.
p = plots(result)
library(ggpubr)
ggarrange(p[[1]], p[[2]], ncol = 2, nrow = 1, common.legend = TRUE)

# Rates of false negatives -------------
t1 = tapply(result[result$tie_effect > 0.5 | result$tie_effect < -0.5,]$`p-value`, result[result$tie_effect > 0.5 | result$tie_effect < -0.5,]$approach, function(x){sum(x >= 0.05)/length(x)})

# Rates of false positives -------------
t2 = tapply(result[result$tie_effect <= 0.5 &  result$tie_effect  >= -0.5,]$`p-value`, result[result$tie_effect<= 0.5 &  result$tie_effect  >= -0.5,]$approach, function(x){sum(x <= 0.05)/length(x)})

summary = data.frame("Approaches" = c(names(t1),  names(t2)), "Error type" = c(rep('False negatives', 7), rep('False positives', 7)) ,"Percent" = c(t1, t2))
summary$Percent = summary$Percent * 100
summary

write.csv(result, "SIM -2 to 2.csv", row.names = FALSE)
write.csv(summary, "SIM -2 to 2 rates of type I and type II errors.csv", row.names = FALSE)
save.image("SIM -2 to 2.RData")

ggplot(result[result$tie_effect>= -0.5 & result$tie_effect <= 0.5,], aes(x = tie_effect,  y = z,group = sim, label = z))+
  geom_point(aes(shape = approach, size = N_id, color = exposure_sigma ),  show.legend = TRUE, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  facet_grid( . ~ detect_effect, space="free") +
  theme(legend.position = 'none')+
  ylab("Estimated effect size (z-score)") +
  xlab("True efect size") +
  theme(axis.text = element_text(size = 9))



