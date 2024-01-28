#########################################################################################
#### Testing simulation sr and dyadic effects
###########################################################
library(STRAND)
source("./1.Codes/2.data_simulation.R")
source("./1.Codes/2.1.STRAND_censoring.R")
N_id = 30
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
indiv =  data.frame(Hairy = Hairy)
individual_predictors = Hairy
individual_effects=matrix(c(1, 1),ncol=1, nrow=2)



NG = sample(c(1,3,7), 1) # Random number of groups
clique = sample(1:NG, N_id, replace = TRUE)
mean.within.GR = sample(c(seq(from = -9, to = 9, by = 1)), 1) # Probability of random ties within a group.
B = matrix(rnorm(NG*NG, mean.within.GR, sd = 1), NG, NG)
mean.between.GR = sample(c(seq(from = 0, to = 9, by = 1)), 1) # Increase randomly the probability of  ties within groups.
diag(B) = diag(B) + rnorm(NG, mean.between.GR, sd = 1)
block = data.frame(Clique=factor(clique))
sr_mu = c(0,0)
sr_sigma = c(1, 1)
sr_rho = 0.5
dr_mu = c(0,0)
dr_sigma = 1
dr_rho = 0.75
V = 1
dyadic_predictors = NULL
dyadic_effects = NULL
exposure_predictors = NULL
exposure_effects = NULL
exposure_sigma = 1
exposure_baseline = 50
int_intercept = c(Inf,Inf)
int_slope = c(Inf,Inf)
simulate.interactions = TRUE


A = simulate_sbm_plus_srm_network_with_measurement_bias(N_id = N_id,
                                                        B =list(B=B),
                                                        V = V,
                                                        groups=block,
                                                        sr_mu = sr_mu,
                                                        sr_sigma = sr_sigma,
                                                        sr_rho = sr_rho,
                                                        dr_mu = dr_mu,
                                                        dr_sigma = dr_sigma,
                                                        dr_rho = dr_rho,
                                                        individual_predictors = individual_predictors,
                                                        dyadic_predictors = dyadic_predictors,
                                                        individual_effects = individual_effects,
                                                        dyadic_effects = dyadic_effects,
                                                        exposure_predictors = exposure_predictors,
                                                        exposure_effects = exposure_effects,
                                                        exposure_sigma = exposure_sigma,
                                                        exposure_baseline = exposure_baseline,
                                                        int_intercept = int_intercept,
                                                        int_slope = int_slope,
                                                        simulate.interactions = simulate.interactions)

nets = list(Grooming = A$network)
exposure_nets = list(Exposure = A$true_samps)

data = make_strand_data_censoring(outcome = nets,
                                  individual_covariates = indiv,
                                  block_covariates = block,
                                  outcome_mode = "binomial",
                                  exposure = exposure_nets,
                                  self_report = NULL,
                                  ground_truth = NULL,
                                  dyadic_covariates = NULL,
                                  censoring = indiv)



fit_social_relations_model_bias <- function (data, focal_regression, target_regression, dyad_regression, 
          mode = "mcmc", return_predicted_network = FALSE, 
          stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL, 
          iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL), 
          priors = NULL){
  if (attributes(data)$class != "STRAND Data Object") {
    stop("fit_latent_network_model() requires a data object of class: STRAND Data Object. Please use make_strand_data() to build your data list.")
  }
  if (!("SRM" %in% attributes(data)$supported_models)) {
    stop("The supplied data are not appropriate for an SRM model. Please ensure that self_report data is single-layer.")
  }
  if (data$N_individual_predictors == 0 & focal_regression != 
      ~1) {
    stop("No individual covariate data has been provided. focal_regression must equal ~ 1 ")
  }
  if (data$N_individual_predictors == 0 & target_regression != 
      ~1) {
    stop("No individual covariate data has been provided. target_regression must equal ~ 1 ")
  }
  if (data$N_dyadic_predictors == 0 & dyad_regression != ~1) {
    stop("No individual covariate data has been provided. dyad_regression must equal ~ 1 ")
  }
  ind_names = colnames(data$individual_predictors)
  dyad_names = names(data$dyadic_predictors)
  if (data$N_dyadic_predictors > 0) {
    dyad_dims = c(data$N_id, data$N_id, length(dyad_names))
    dyad_dat = list()
    for (i in 1:dyad_dims[3]) {
      dyad_dat[[i]] = c(data$dyadic_predictors[[i]])
    }
    dyad_dat = as.data.frame(do.call(cbind, dyad_dat))
    colnames(dyad_dat) = dyad_names
    dyad_model_matrix = model.matrix(dyad_regression, dyad_dat)
    dyad_dat_out = array(NA, c(dyad_dims[1], dyad_dims[2], 
                               ncol(dyad_model_matrix)))
    for (i in 1:ncol(dyad_model_matrix)) {
      dyad_dat_out[, , i] = matrix(dyad_model_matrix[, 
                                                     i], nrow = dyad_dims[1], ncol = dyad_dims[2])
    }
    dimnames(dyad_dat_out)[[3]] = colnames(dyad_model_matrix)
    data$dyad_set = dyad_dat_out
  }else {
    data$dyad_set = array(1, c(data$N_id, data$N_id, 1))
  }
  if (data$N_individual_predictors > 0) {
    data$focal_set = model.matrix(focal_regression, data$individual_predictors)
    data$target_set = model.matrix(target_regression, data$individual_predictors)
  }else {
    data$focal_set = matrix(1, nrow = data$N_id, ncol = 1)
    data$target_set = matrix(1, nrow = data$N_id, ncol = 1)
  }
  
  # censoring----------------
  if (data$N_censoring_predictors > 0) {
    data$censoring_set = model.matrix(focal_regression, data$censoring_predictors)
  }else {
    data$censoring_set = matrix(1, nrow = data$N_id, ncol = 1)
  }
  
  
  data$N_params = c(ncol(data$focal_set), ncol(data$target_set), 
                    dim(data$dyad_set)[3], ncol(data$censoring_set))
  
  data$export_network = ifelse(return_predicted_network == 
                                 TRUE, 1, 0)
  if (is.null(priors)) {
    data$priors = make_priors()
  }else {
    data$priors = priors
  }


  model = cmdstanr::cmdstan_model("1.Codes/2.3.STRAND_censoring.stan")
  if (mode == "mcmc") {
    fit = model$sample(data = unclass(data), seed = stan_mcmc_parameters$seed, 
                       chains = stan_mcmc_parameters$chain, parallel_chains = stan_mcmc_parameters$parallel_chains, 
                       refresh = stan_mcmc_parameters$refresh, iter_warmup = stan_mcmc_parameters$iter_warmup, 
                       iter_sampling = stan_mcmc_parameters$iter_sampling, 
                       max_treedepth = stan_mcmc_parameters$max_treedepth, 
                       adapt_delta = stan_mcmc_parameters$adapt_delta)
  }
  if (mode == "vb") {
    print("Variational inference is fast, but not always dependable. We recommend using vb only for test runs.")
    fit = model$variational(data = unclass(data), seed = 123, 
                            output_samples = 2000)
  }
  if (mode == "optim") {
    print("Optimazation is fast, but not always dependable. We recommend using optim only for test runs.")
    fit = model$optimize(data = unclass(data), seed = 123)
  }
  if (!mode %in% c("mcmc", "vb", "optim")) {
    stop("Must supply a legal mode value: mcmc, vb, or optim.")
  }
  bob = list(data = data, fit = fit, return_predicted_network = return_predicted_network)
  attr(bob, "class") = "STRAND Model Object"
  attr(bob, "fit_type") = mode
  attr(bob, "model_type") = "SRM"
  return(bob)
}


focal_regression = ~ Hairy
target_regression = ~ Hairy
dyad_regression = ~  1
mode="mcmc"
return_predicted_network = TRUE
stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                            iter_warmup = 1000, iter_sampling = 1000,
                            max_treedepth = NULL, adapt_delta = .98)
priors = NULL




