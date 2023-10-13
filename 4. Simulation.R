###############################
##  Simulations parameters ####
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
    ncores = 1 # number of cores
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
    source("3. Data simulation.R")

    get_res = function(x) {
      y = c(coef(x)[2],confint(x, level = 0.9)[2,])
      return(y)
    }

    #' "Bayesian P value". 1)The integral over the posterior from negative infinity to zero,
    #'                    2) the integral over the posterior from zero to infinity. T
    #' This is easy to approximate from the samples:

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

    RESULTS = NULL
    # Make data--------------------------------
    N_id = grid_subsample$N_id[i]

    ## Block data -----------------
    if(is.null(B)){
      clique = sample(1:G, N_id, replace=TRUE)
      B = matrix(-9, nrow=G, ncol=G)
      diag(B) = diagB
      B[1,3] = -8.9
      B[3,2] = -7.9
    }


    block = data.frame(Clique=factor(clique))

    Clique = rep(1, N_id)
    Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)

    ## Simulated data -----------------
    A = simulate_sbm_plus_srm_network_with_measurement_bias(
      N_id = N_id,
      B=list(B=B),
      V=V,
      groups=groups,
      individual_predictors=Hairy,
      individual_effects=matrix(c(grid_subsample$tie_effect[i], grid_subsample$tie_effect[i]),ncol=1, nrow=2),
      sr_sigma = sr_sigma,
      sr_rho = sr_rho,
      dr_sigma = dr_sigma,
      dr_rho = dr_rho,
      exposure_predictors = cbind(rep(1,N_id),Hairy),
      exposure_effects = c(-1, grid_subsample$hairy_detect_effect[i]),
      exposure_sigma = exposure_sigma,
      exposure_baseline =exposure_baseline,
      int_bias = T
    )

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
      if(blocks){
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

      colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value', 'approach', 'sim')
      RESULTS = rbind(RESULTS, result)
    }

    #  Rates of interactions unweighted--------------------------------
    tie_strength = A$network/(1+A$true_samps)
    colnames(tie_strength) = rownames(tie_strength) = 1:ncol(tie_strength)
    df = ANTs:::df.create(tie_strength)
    df = met.strength(tie_strength, df = df, dfid = 1)
    df = met.degree(tie_strength, df = df, dfid = 1)
    df$hair = indiv$Hairy
    test1.1 = lm(strength ~ hair, data = df)

    ants_est = summary(test1.1)$coefficients[2,4]
    result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,],ants_est)))))
    rownames(result) = NULL
    result$approach = 'Rates unweighted'
    result$sim = i
    colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value', 'approach', 'sim')
    RESULTS = rbind(RESULTS, result)


    #  Rates of interactions weighted--------------------------------
    test1.1 = lm(strength  ~ hair, data = df, weights = A$true_exposure)
    result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,],ants_est)))))
    rownames(result) = NULL
    result$approach = 'Rates weighted'
    result$sim = i
    result$sim = i
    colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value', 'approach', 'sim')
    RESULTS = rbind(RESULTS, result)

    #  SRI of interactions unweighted--------------------------------
    m = A$network
    fa = matrix(0, ncol = nrow(m), nrow = nrow(m))
    for (a in 1:nrow(fa)) {
      fa[a,] = A$true_exposure
    }

    fb = matrix(0, ncol = nrow(m), nrow = nrow(m))
    for (a in 1:nrow(m)) {
      fb[,a] = A$true_exposure
    }

    ya = abs(fa - m)

    yb = abs(fb - m)

    sri <- ((m) /(m + ya + yb ))
    sri[is.na(sri)] = 0
    colnames(sri)= rownames(sri) = 1:ncol(sri)
    df = ANTs:::df.create(sri)
    df = met.strength(sri, df = df, dfid = 1)
    df = met.degree(sri, df = df, dfid = 1)
    df$hair = indiv$Hairy

    test1.1 = lm(strength ~ hair, data = df, weights = A$true_exposure)
    ants_est = summary(test1.1)$coefficients[2,4]
    result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,],ants_est)))))
    rownames(result) = NULL
    result$approach = 'SRI unweigthed'
    result$sim = i
    colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value', 'approach', 'sim')
    RESULTS = rbind(RESULTS, result)

    #  SRI of interactions weighted--------------------------------
    test1.1 = lm(strength ~ hair, data = df, weights = A$true_exposure)
    ants_est = summary(test1.1)$coefficients[2,4]
    result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,],ants_est)))))
    rownames(result) = NULL
    result$approach = 'SRI weighted'
    result$sim = i
    colnames(result) = c('scale(hair)','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'p-value', 'approach', 'sim')
    RESULTS = rbind(RESULTS, result)
  }

  parallel::stopCluster(cl)
  gc()
  registerDoSEQ()
  return(do.call('rbind', r))
}

result = simulations(Reps = 50, strand = F, ncores = 4, blocks = TRUE)
result$resid = result$tie_effect - result$`scale(hair)`
result$z = result$`scale(hair)`/(result$`95 %` - result$`5 %`)

# Plot ----------
library(ggpubr)
require(ggplot2)
p1 = ggplot(result, aes(x=z,  y= tie_effect , xmin=result[,2], xmax=result[,3],  group = sim, label = z))+
  #geom_errorbar() +
  geom_point(aes(color = N_id, size = detect_effect), alpha = 0.5) +
  facet_grid( . ~ approach, space="free") +
  xlim(c(-5,5))+
  ylim(c(-5,5))+
  geom_abline()+
  ylab("True efect size") +
  xlab("Estimated effect size") +
  theme(axis.text = element_text(size = 14))  +
  theme(axis.title = element_text(size = 14))
p1

p2 = ggplot(result, aes(x = resid, y = detect_effect, color = approach))+
  geom_point(aes(size = detect_effect), show.legend = FALSE, alpha = 0.5)+
  geom_vline(xintercept = 0, linetype = "dashed")+
  xlab("Difference between true effect size and estimated effect size")+
  ylab("True effect size")+
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14))


p3 = ggplot(result, aes(x = tie_effect,  y = `p-value`,  group = approach, color = approach))+
  geom_point()+
  geom_hline(yintercept=0.05, linetype="dashed",
             color = "red", size=1)
p3
ggarrange(p1, p2, ncol  = 1, nrow = 2, common.legend = FALSE)

