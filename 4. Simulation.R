###############################
#######  Functions ############
###############################
# Objects to store the results.
get_res = function(x) {
  y = c(coef(x)[2],confint(x, level = 0.9)[2,])
  return(y)
}


###############################
##  Simulations parameters ####
###############################
# Network variables----------------------
#' @param N_id Number of individuals
#' @param B Tie probabilities
#' @param V Blocking variables (subgroups within the network)
#' @param groups Subgroup IDs

# Sender/receiver variables----------------------
#' @param sr_mu Average sender (cell 1) and receiver (cell 2) effect log odds
#' @param sr_sigma Sender (cell 1) and receiver (cell 2) effect variances
#' @param sr_rho Correlation of sender and receiver effects

# Dyadic variables---------------------------
#' @param dr_mu Average i to j dyad effect (cell 1) and j to i dyad effect (cell 2) log odds
#' @param dr_sigma Variance of dyad effects
#' @param dr_rho Correlation of i to j dyad effect and j to i dyad effect

# Individuals characteristics----------------------
#' @param individual_predictors A matrix of covariates
#' @param dyadic_predictors An array of covariates
#' @param individual_effects The effects of predictors on sender effects (row 1) and receiver effects (row 2)
#' @param dyadic_effects The effects of predictors on dyadic ties

# Biases variables----------------------
#' @param exposure_predictors A matrix of covariates
#' @param exposure_effects A vector of slopes
#' @param exposure_sigma Variance in exposure (observations) random effects
#' @param exposure_baseline Baseline exposure (observations) rate

simulations <- function(
    Reps = 100,
    N_id =  seq(30, 100, by = 20),
    hairy_tie_effect = seq(-4, 4, by = 0.5),
    hairy_detect_effect = seq(-4, 4, by = 0.5),
    B = NULL,
    V = 1,
    groups=NULL,
    sr_sigma = c(0.1, 0.1),
    sr_rho = 0.0,
    dr_sigma = 1.2,
    dr_rho = 0.0,
    exposure_sigma = 2.9,
    exposure_baseline = 40
){
  require(ANTs)
  require(ggplot2)
  source("3. Data simulation.R")
  RESULTS = NULL

  ########################################
  #######  Create grid of simulations ####
  ########################################
  grid = expand.grid(N_id, hairy_tie_effect, hairy_detect_effect)
  grid_subsample = grid[sample(1:nrow(grid), Reps, replace = FALSE),]
  colnames(grid_subsample) = c("N_id", "tie_effect", "detect_effect")

  #############################
  #######  Testing methods ####
  #############################
  for(i in 1:nrow(grid_subsample)){
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
    result$exposure =
    colnames(result) = c('hair','5 %','95 %', 'N_id', 'tie_effect', 'detect_effect', 'approach', 'sim')
    RESULTS = rbind(RESULTS, result)

    #  Rates of interactions unweighted--------------------------------
    tie_strength = A$network/(1+A$true_samps)
    colnames(tie_strength) = rownames(tie_strength) = 1:grid_subsample$N_id[i]
    df = ANTs:::df.create(tie_strength)
    df = met.strength(tie_strength, df = df, dfid = 1)
    df$hair = indiv$Hairy
    test1.1 = lm(strength ~ hair, data = df)

    result = as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid_subsample[i,])))))
    rownames(result) = NULL
    result$approach = 'Rates unweighted'
    result$sim = i
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
    df$hair = indiv$Hairy

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
  return(RESULTS)
}

result = simulations()

# Plot ----------

ggplot(result, aes(x=hair,  y= detect_effect, ymin=result[,2], ymax=result[,3], color = as.factor(sim), group = sim, label = hair))+
  geom_linerange() +
  geom_point(aes(size = N_id)) +
  facet_grid( . ~ approach, space="free") +
  #theme(legend.position = 'none')+
  xlim(c(range(result[,2])[1], range(result[,3])[2]))+
  geom_abline()+
  xlab("True efect size") +
  ylab("Estimated effect size") +
  theme(axis.text = element_text(size = 14))  +
  theme(axis.title = element_text(size = 14))

#ggsave("Model Comparison.pdf",p1,width=12,height=3)

