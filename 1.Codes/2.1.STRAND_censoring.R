library(STRAND)
source("./1.Codes/2.data_simulation.R")

N_id = 30
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
indiv =  data.frame(Hairy = Hairy)
individual_effects=matrix(c(1, 1),ncol=1, nrow=2)
A = simulate_sbm_plus_srm_network_with_measurement_bias(N_id = N_id, individual_predictors= Hairy, individual_effects=individual_effects)

nets = list(Grooming = A$network)
exposure_nets = list(Exposure = A$true_samps)

NG = sample(c(1,3,7), 1) # Random number of groups
clique = sample(1:NG, N_id, replace = TRUE)
mean.within.GR = sample(c(seq(from = -9, to = 9, by = 1)), 1) # Probability of random ties within a group.
B = matrix(rnorm(NG*NG, mean.within.GR, sd = 1), NG, NG)
mean.between.GR = sample(c(seq(from = 0, to = 9, by = 1)), 1) # Increase randomly the probability of  ties within groups.
diag(B) = diag(B) + rnorm(NG, mean.between.GR, sd = 1)
block = data.frame(Clique=factor(clique))

make_strand_data_censoring = function(outcome = nets,
                                       individual_covariates = indiv,
                                       block_covariates = block,
                                       outcome_mode = "binomial",
                                       exposure = exposure_nets,
                                       self_report = NULL,
                                       ground_truth = NULL,
                                       dyadic_covariates = NULL,
                                       censoring = indiv){

  outcome_mode_numeric = NULL
  if (outcome_mode == "bernoulli") {
    outcome_mode_numeric = 1
  }
  if (outcome_mode == "binomial") {
    outcome_mode_numeric = 2
  }
  if (outcome_mode == "poisson") {
    outcome_mode_numeric = 3
  }
  if (is.null(outcome_mode_numeric)) 
    stop("outcome_mode not supported")
  if (is.null(self_report)) {
    self_report = outcome
  }
  if (is.null(exposure)) {
    if (is.null(self_report)) 
      stop("self_report must be a list of matrices.")
    if (!is.list(self_report)) 
      stop("self_report must be a list of matrices.")
    if (!length(self_report) %in% c(1, 2)) 
      stop("self_report must be a list length 1 or 2.")
    for (i in 1:length(self_report)) {
      if (outcome_mode == "bernoulli") {
        if (!all(self_report[[i]] %in% c(0, 1))) 
          stop("self_report must be binary 0 or 1")
      }
    }
  }
  if (!is.null(exposure)) {
    if (length(self_report) != length(exposure)) 
      stop("self_report and exposure must be lists of matrices equal in length.")
  }
  if (!is.null(ground_truth)) {
    if (!is.list(ground_truth)) 
      stop("ground_truth must be a list of matrices.")
  }
  if (!is.null(block_covariates)) {
    if (!is.data.frame(block_covariates)) 
      stop("block_covariates must be a data frame.")
    for (i in 1:dim(block_covariates)[2]) {
      if (!is.factor(block_covariates[, i])) 
        stop("block_covariates must be factor variables.")
    }
  }
  if (!is.null(individual_covariates)) {
    if (!is.data.frame(individual_covariates)) 
      stop("individual_covariates must be a data frame.")
  }
  if (!is.null(dyadic_covariates)) {
    if (!is.list(dyadic_covariates)) 
      stop("dyadic_covariates must be a list of matrices.")
  }
  N_id = dim(self_report[[1]])[1]
  N_responses = length(self_report)
  outcomes = array(NA, c(N_id, N_id, N_responses))
  if (is.null(exposure)) {
    exposure_on = 0
    exposure_risk = array(0, c(N_id, N_id, N_responses))
  }else {
    exposure_on = 1
    exposure_risk = array(NA, c(N_id, N_id, N_responses))
    for (i in 1:length(exposure)) {
      exposure_risk[, , i] = exposure[[i]]
    }
  }
  for (i in 1:length(self_report)) {
    outcomes[, , i] = self_report[[i]]
  }
  if (is.null(block_covariates)) {
    block_covariates = rep(1, N_id)
    N_groups_per_type = 1
    N_block_types = 0
    group_ids_character = rep("Any", N_id)
    group_ids = block_covariates
    group_ids_levels = "No Blocks"
  }else {
    N_block_types = length(block_covariates[1, ])
    N_groups_per_type = rep(NA, N_block_types)
    group_ids_character = array(NA, c(N_id, N_block_types))
    group_ids = array(NA, c(N_id, N_block_types))
    group_ids_levels = vector("list", N_block_types)
    for (i in 1:N_block_types) {
      N_groups_per_type[i] = max(as.numeric(block_covariates[, 
                                                             i]))
      group_ids_character[, i] = as.character(block_covariates[, 
                                                               i])
      group_ids[, i] = as.numeric(block_covariates[, i])
      group_ids_levels[[i]] = levels(block_covariates[, 
                                                      i])
    }
    group_ids = data.frame(group_ids)
    colnames(group_ids) = colnames(block_covariates)
  }
  if (is.null(ground_truth)) {
    N_networktypes = N_responses
    N_periods = 0
    flows = 0
  }else {
    N_networktypes = N_responses + 1
    N_periods = length(ground_truth)
    flows = array(NA, c(N_id, N_id, N_periods))
    for (i in 1:length(ground_truth)) flows[, , i] = ground_truth[[i]]
  }
  if (is.null(individual_covariates)) {
    N_individual_predictors = 0
    individual_predictors = 0
  }else {
    N_individual_predictors = dim(individual_covariates)[2]
    individual_predictors = individual_covariates
  }
  if (is.null(dyadic_covariates)) {
    N_dyadic_predictors = 0
    dyadic_predictors = 0
  }else {
    N_dyadic_predictors = length(dyadic_covariates)
    dyadic_predictors = dyadic_covariates
  }
  if (N_responses == 1) {
    if (max(N_groups_per_type) > 1) {
      supported_models = c("SRM", "SBM", "SRM+SBM")
    }
    else {
      supported_models = c("SRM")
    }
  }
  if (N_responses == 2 & N_networktypes == 2) {
    supported_models = c("LNM")
  }
  if (N_responses == 2 & N_networktypes == 3) {
    supported_models = c("LNM", "LNM+Flows")
  }
  
  # Adding censoring bias
  if(is.null(censoring)){
    N_censoring_predictors = 0
    censoring_predictors = 0
    model_dat = list(N_networktypes = N_networktypes, N_id = N_id, 
                     N_responses = N_responses, N_periods = N_periods, N_individual_predictors = N_individual_predictors, 
                     N_dyadic_predictors = N_dyadic_predictors, outcomes = outcomes, 
                     flows = flows, individual_predictors = individual_predictors, 
                     dyadic_predictors = dyadic_predictors, N_block_predictors = N_block_types, 
                     N_groups_per_block_type = N_groups_per_type, block_predictors = group_ids, 
                     outcome_mode = outcome_mode_numeric, exposure = exposure_risk)
  }else{
    N_censoring_predictors = dim(censoring)[2]
    censoring_predictors = censoring
    model_dat = list(N_networktypes = N_networktypes, N_id = N_id, 
                     N_responses = N_responses, N_periods = N_periods, N_individual_predictors = N_individual_predictors, 
                     N_dyadic_predictors = N_dyadic_predictors, outcomes = outcomes, 
                     flows = flows, individual_predictors = individual_predictors, 
                     dyadic_predictors = dyadic_predictors, N_block_predictors = N_block_types, 
                     N_groups_per_block_type = N_groups_per_type, block_predictors = group_ids, 
                     outcome_mode = outcome_mode_numeric, exposure = exposure_risk,
                     N_censoring_predictors = N_censoring_predictors, censoring_predictors = censoring_predictors
    )
    attr(model_dat, "censoring") = N_censoring_predictors
  }
  

  
  attr(model_dat, "class") = "STRAND Data Object"
  attr(model_dat, "supported_models") = supported_models
  attr(model_dat, "group_ids_character") = group_ids_character
  attr(model_dat, "group_ids_levels") = group_ids_levels
  colnames(attr(model_dat, "group_ids_character")) = colnames(model_dat$block_predictors)
  names(attr(model_dat, "group_ids_levels")) = colnames(attr(model_dat, 
                                                             "group_ids_character"))
  return(model_dat)
  
}


model_dat = make_strand_data_censoring()

