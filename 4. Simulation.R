###############################
#######  Function ############
###############################
# Objects to store the results.
get_res = function(x) {
  y = c(coef(x)[2],confint(x, level = 0.9)[2,])
  return(y)
}


###############################
##  Simulations parameters ####
###############################
simulations <- function(
    Reps = 100,
    N_id =  seq(30, 100, by = 20),
    hairy_tie_effect = seq(-4, 4, by = 1),
    hairy_detect_effect = seq(-4, 4, by = 1),
    B = NULL,
    V = 1,
    groups=NULL,
    sr_sigma = c(0.3, 1.5),
    sr_rho = 0.0,
    dr_sigma = 1,
    dr_rho = 0.0,
    exposure_sigma = 2.9,
    exposure_baseline = 40,
    strand = T, # use STRAN model
    ncores = 1 # number of cores
){

  require(parallel)
  require(foreach)
  require(doParallel)
  source("3. Data simulation.R")

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
        # Make data--------------------------------
        Clique = rep(1, grid_subsample$N_id[i])
        Hairy = matrix(rnorm(grid_subsample$N_id[i], 0, 1), nrow=grid_subsample$N_id[i], ncol=1)
        A = simulate_sbm_plus_srm_network_with_measurement_bias(
          N_id = grid_subsample$N_id[i],
          B=list(B=B),
          V=V,
          groups=groups,
          individual_predictors=Hairy,
          individual_effects=matrix(c(grid_subsample$tie_effect[i], grid_subsample$tie_effect[i]),ncol=1, nrow=2),
          sr_sigma = sr_sigma,
          sr_rho = sr_rho,
          dr_sigma = dr_sigma,
          dr_rho = dr_rho,
          exposure_predictors = cbind(rep(1,grid_subsample$N_id[i]),Hairy),
          exposure_effects = c(-1, grid_subsample$hairy_detect_effect[i]),
          exposure_sigma = exposure_sigma,
          exposure_baseline =exposure_baseline
        )

        # STRAND--------------------------------
        if(strand){
          nets = list(Grooming = A$network)
          exposure_nets = list(Exposure = A$true_samps)
          block = data.frame(Clique=factor(Clique))
          indiv =  data.frame(Hairy = Hairy)
          model_dat = make_strand_data(outcome = nets,
                                       individual_covariates = indiv,
                                       block_covariates = block,
                                       outcome_mode = "binomial",
                                       exposure = exposure_nets
          )

          # !!  I replaced fit_block_plus_social_relations_model by fit_social_relations_model
          fit =  fit_social_relations_model(data=model_dat,
                                            focal_regression = ~ Hairy,
                                            target_regression = ~ Hairy,
                                            dyad_regression = ~  1,
                                            mode="mcmc",
                                            stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                        iter_warmup = 1000, iter_sampling = 1000,
                                                                        max_treedepth = NULL, adapt_delta = .98)
          )
          res = summarize_strand_results(fit)

          result = as.data.frame(t(data.frame(as.numeric(unlist(c(unlist(res$summary[2,2:4]), grid_subsample[i,]))))))
          rownames(result) = NULL
          result$approach = 'strand'
          result$sim = i

          colnames(result) = c('hair','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'approach', 'sim')
          RESULTS = rbind(RESULTS, result)
        }


        #  Rates of interactions unweighted--------------------------------
        tie_strength = A$network/(1+A$true_samps)
        colnames(tie_strength) = rownames(tie_strength) = 1:grid_subsample$N_id[i]
        df = ANTs:::df.create(tie_strength)
        df = met.strength(tie_strength, df = df, dfid = 1)
        df$hair = 1:grid_subsample$N_id[i]
        test1.1 = lm(strength ~ hair, data = df)

        result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,])))))
        rownames(result) = NULL
        result$approach = 'Rates unweighted'
        result$sim = i
        colnames(result) = c('hair','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'approach', 'sim')
        RESULTS = rbind(RESULTS, result)


        #  Rates of interactions weighted--------------------------------
        test1.1 = lm(strength  ~ hair, data = df, weights = A$true_exposure)
        result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,])))))
        rownames(result) = NULL
        result$approach = 'Rates weighted'
        result$sim = i
        RESULTS = rbind(RESULTS, result)

        #  SRI of interactions unweighted--------------------------------
        m = A$network
        fa = matrix(0, ncol = nrow(m), nrow = nrow(m))
        for (a in 1:nrow(fa)) {
          fa[a,] = A$true_exposure
        }

        fb = matrix(0, ncol = nrow(m), nrow = nrow(m))
        for (a in 1:nrow(m)) {
          fb[a,] = A$true_exposure
        }

        ya = abs(fa - m)

        yb = abs(fb - m)

        sri <- ((m) /(m + ya + yb ))
        sri[is.na(sri)] = 0
        colnames(sri)= rownames(sri) = 1:ncol(sri)
        df = ANTs:::df.create(sri)
        df = met.strength(sri, df = df, dfid = 1)
        df$hair = 1:grid_subsample$N_id[i]

        test1.1 = lm(strength ~ hair, data = df, weights = A$true_exposure)
        result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,])))))
        rownames(result) = NULL
        result$approach = 'SRI unweigthed'
        result$sim = i
        RESULTS = rbind(RESULTS, result)

        #  SRI of interactions weighted--------------------------------
        test1.1 = lm(strength ~ hair, data = df, weights = A$true_exposure)
        result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,])))))
        rownames(result) = NULL
        result$approach = 'SRI weighted'
        result$sim = i
        RESULTS = rbind(RESULTS, result)

      }
  return(do.call('rbind', r))
}

result = simulations(Reps = 10, strand = T)
result$resid = result$detect_effect - result$hair
# Plot ----------
library(ggpubr)
require(ggplot2)
p1 = ggplot(result, aes(x=hair,  y= detect_effect, xmin=result[,2], xmax=result[,3], color = approach, group = sim, label = hair))+
  geom_errorbar() +
  geom_point() +
  facet_grid( . ~ approach, space="free") +
  theme(legend.position = 'none')+
  xlim(c(-5,5))+
  ylim(c(-5,5))+
  geom_abline()+
  ylab("True efect size") +
  xlab("Estimated effect size") +
  theme(axis.text = element_text(size = 14))  +
  theme(axis.title = element_text(size = 14))


p2 = ggplot(result, aes(x = resid, y = detect_effect, color = approach))+
  geom_point()+
  geom_vline(xintercept = 0, linetype = "dashed")+
  xlab("Difference between True efect size and Estimated effect size")+
  ylab("True efect size")+
  theme(axis.text = element_text(size = 14))+
  theme(axis.title = element_text(size = 14))


ggarrange(p1, p2, ncol = 1, nrow = 2, common.legend = TRUE)

