source('1.Codes/2.1.STRAND_censoring.R')
source('1.Codes/2.2.STRAND fit_social_relations_model_bias.R')

# "Bayesian P value".---------------
P_se <- function(x) {
  M_x <- mean(x)

  if (M_x < 0) {
    N_x <- length(x)
    P_x <- length(which(x > 0))
  } else {
    N_x <- length(x)
    P_x <- length(which(x < 0))
  }

  return(P_x / N_x)
}

# Pre-network permutation-----------
preNetPerm <- function(df, focal = "focal", alters = "receiver", ctrl = "focal", nperm = 10, progress = T, index = "sri") {
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

  result <- ANTs:::perm_dataStream1_focal(GBI, GBI2,
    nperm = nperm,
    progress = progress, method = index
  )
  result <- lapply(seq_along(result), function(x, Vecids, i) {
    colnames(x[[i]]) <- Vecids
    rownames(x[[i]]) <- Vecids
    attr(x[[i]], "permutation") <- i
    return(x[[i]])
  }, x = result, Vecids = Vecids)
  return(result)
}

# Build results-------------
get_res <- function(x) {
  y <- c(coef(x)[2], confint(x, level = 0.9)[2, ])
  return(y)
}

get_result <- function(test1.1, strand = TRUE, name, original.matrix, exposure, grid) {
  if (strand) {

    # Compare estimated network
    N_id = grid$N_id
    m <- matrix(0, ncol = N_id, nrow = N_id)
    for (a in 1:N_id) {
      x <- test1.1$samples$predicted_network_sample[, , a]
      m[, a] <- apply(x, 2, mean)
    }
    M = original.matrix[,,1] / exposure[,,1]
    diag(M) = 0
    cor <- cor(met.strength(m), met.strength(M))

    # "Bayesian P value"
    strand_est <- P_se(test1.1$sample$srm_model_samples$focal_coeffs[, 1])
    a <- test1.1$sample$srm_model_samples$focal_coeffs[, 1]
    strand_se <- sqrt(sum((a - mean(a))^2 / (length(a) - 1))) / sqrt(length(a))

    result <- as.data.frame(t(data.frame(as.numeric(unlist(c(unlist(test1.1$summary[2, 2:4]), grid, strand_est, strand_se))))))

    rownames(result) <- NULL
    result$approach <- "strand"
    
    result$cor <- cor

    colnames(result) <- c("scale(hair)", "5 %", "95 %", "N_id", "tie_effect", "detect_effect",  "sim", "p-value", "se",
                          "approach",  'prediction')
    return(result)
  } else {
    cor <- cor(test1.1$fitted.values, test1.1$model$strength)
    ants_p <- summary(test1.1)$coefficients[2, 4]
    ants_se <- summary(test1.1)$coefficients[2, 2]
    result <- as.data.frame(t(data.frame(unlist(c(unlist(get_res(test1.1)), grid, ants_p, ants_se)))))
    rownames(result) <- NULL
    result$approach <- name
    result$cor <- cor
    colnames(result) <- c("scale(hair)", "5 %", "95 %", "N_id", "tie_effect", "detect_effect", "sim", "p-value", "se",
                          "approach",  'prediction')
    return(result)
  }
}

# make_all_data-------------
make_all_data <- function(
    # Other parameters
    N_id = 50,
    # Sender-receiver parameters
    individual_effects_sr = c(1.7, 0.3),
    sr_sigma = c(1.4, 0.8),
    sr_rho = 0.5,
    # Dyadic parameters
    dr_sigma = 1.2,
    dr_rho = 0.8,
    # Observation bias parameters
    exposure_bias = TRUE,
    exposure_mu = 1.9,
    exposure_sigma = 2.9,
    exposure_max = 40,
    exposure_effects = -1.1,
    # Censoring bias parameters
    simulate_censoring = TRUE,
    censoring_mu = -1.9,
    censoring_sigma = 1.0,
    censoring_effects = -1.1
    ) {
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

  ########################## And finally organize for each software package
  data_out <- NULL

  ###################################### STRAND---------------------
  nets <- list(Grooming = A$network)
  indiv <- data.frame(Hairy = Individual_Factor)
  block <- data.frame(Clique = factor(Clique))
  exposure_nets <- list(Exposure = A$true_samps)

  model_dat <- make_strand_data_censoring(
    outcome = nets,
    individual_covariates = indiv,
    block_covariates = block,
    outcome_mode = "poisson",
    exposure = exposure_nets,
    detected = A$detected,
    trials = A$trials
  )

  data_out[[1]] <- model_dat

  ###################################### ANTs Rate---------------------
  tie_strength <- A$network / (1 + A$true_samps)
  colnames(tie_strength) <- rownames(tie_strength) <- 1:ncol(tie_strength)
  df <- ANTs:::df.create(tie_strength)
  df$strength <- met.strength(tie_strength)
  df$hair <- indiv$Hairy
  df$exposure <- A$true_exposure

  data_out[[2]] <- df

  ###################################### ANTs Rate Permutation---------------------
  data_out[[3]] <- perm <- perm.net.nl(df, labels = "hair", nperm = 10000, progress = FALSE)

  ###################################### ANTs SRI---------------------
  m <- A$network
  fa <- fb <- matrix(0, ncol = nrow(m), nrow = nrow(m))

  for (a in 1:nrow(fa)) {
    fa[a, ] <- A$true_exposure
    fb[, a] <- A$true_exposure
  }

  ya <- abs(fa - m)
  yb <- abs(fb - m)

  sri <- (m / (m + ya + yb))
  sri[is.na(sri)] <- 0
  colnames(sri) <- rownames(sri) <- 1:ncol(sri)
  #df <- ANTs:::df.create(sri)
  df = data.frame('id' = colnames(A$network))

  df$strength <- met.strength(sri)
  df$hair <- indiv$Hairy
  df$exposure <- A$true_exposure

  data_out[[4]] <- df

  ###################################### ANTs SRI NL Permutation---------------------
  data_out[[5]] <- perm.net.nl(df, labels = "hair", nperm = 10000, progress = FALSE)

  ###################################### ANTs SRI DS Permutation---------------------
  d = A$interactions
  d = d[d$outcome > 0,]
  d = d[rep(row.names(d), d$outcome), ]
  d$outcome = 1
  rownames(d) = NULL
  newD = NULL
  for( a in 1:nrow(d)){
    ID = c(d$focal[a], d$target[a])
    duration = d$duration[a]
    obs = paste(paste(ID, collapse = "_"), as.character(a), sep = '_')
    newD = rbind(newD, data.frame(ID, duration,obs))
  }

  data_out[[6]] <- list('individual_predictors' = A$individual_predictors,
                        'exposure' = A$true_exposure,
                        'oda' = newD,
                        'perms' = perm.ds.grp(df = newD, scan ='obs',  nperm = 10000, index = 'sri', progress = FALSE))

  ###################################### BiSON---------------------
  data_out[[7]] <- A

  ##################################### Return object---------------------
  names(data_out) <- c("STRAND", "ANTs_Rate", "ANTs_Rate_Perm", "ANTs_SRI", "ANTS_SRI_Perm", "ANTS_SRI_PreNet", "BiSON")
  return(data_out)
}

# Run models-------------
run_STRAND_model <- function(model_dat, grid) {
  fit <- fit_social_relations_model_bias(
    data = model_dat,
    block_regression = ~Clique,
    focal_regression = ~Hairy,
    target_regression = ~Hairy,
    censoring_regression = ~Hairy,
    dyad_regression = ~1,
    mode = "mcmc", return_predicted_network = TRUE,
    stan_mcmc_parameters = list(
      chains = 1, parallel_chains = 1, refresh = 1000,
      iter_warmup = 1000, iter_sampling = 1000,
      max_treedepth = NULL, adapt_delta = 0.98
    )
  )

  res <- summarize_strand_results(fit)

  result <- get_result(res,
    strand = TRUE,
    original.matrix = model_dat$outcomes, 
    exposure = model_dat$exposure,
    grid = grid
  )

  return(result)
}

# For each ANTs models we run the weigthed and unweigthed version
run_ANTs_Rate_model <- function(model_dat, grid) {
  # Unweigthed
  test1.1 <- lm(strength ~ hair, data = model_dat)
  result1 <- get_result(test1.1, 
                        strand = FALSE,
                        name = "ANTs Rates",
                        grid = grid)

  # Weigthed
  test1.2 <- lm(strength ~ hair, data = model_dat, weights = model_dat$exposure)
  result2 <- get_result(test1.2, 
                        strand = FALSE, 
                        name = "ANTs Rates Weighted",
                        grid = grid)

  return(rbind(result1, result2))
}

run_ANTs_Rate_Perm_model <- function(model_dat, grid) {
  # Unweighted
  test1.1 <- lm(strength ~ hair, data = model_dat[[1]])
  result1 <- get_result(test1.1, 
                        strand = FALSE, 
                        name = "ANTs Rates",
                        grid = grid)

  # perm = perm.net.nl(model_dat, labels = 'hair', nperm = 10000, progress = FALSE)
  ptest <- stat.lm(model_dat, formula = strength ~ hair, progress = FALSE)
  ptest <- ANTs:::stat.p(c(ptest$Original.model$coefficients[2, 1], ptest$permutations$hair))[3]

  result1$approach <- paste(result1$approach, "permuted", sep = " ")
  result1$`p-value` <- ptest

  # Weighted
  test1.2 <- lm(strength ~ hair, data = model_dat[[1]], weights = model_dat[[1]]$exposure)
  result2 <- get_result(test1.2, 
                        strand = FALSE, 
                        name = "ANTs Rates Weighted",
                        grid = grid)

  ptest <- unlist(lapply(model_dat, function(x, exposure) {
    m <- lm(data = x, formula = strength ~ hair, weights = exposure)
    summary(m)$coefficients[2, 1]
  }, exposure = model_dat$exposure))

  ptest <- ANTs:::stat.p(ptest)[3]
  result2$approach <- paste(result2$approach, "permuted", sep = " ")
  result2$`p-value` <- ptest

  return(rbind(result1, result2))
}

run_ANTs_SRI_model <- function(model_dat, grid) {
  # Unweigthed
  test1 = lm(strength ~ hair, data = model_dat)
  result1 = get_result(test1,
                       strand = FALSE, 
                       name = 'ANTs SRI',
                       grid = grid)
  
  # Weigthed
  test2 = lm(strength ~ hair, data = model_dat, weights = model_dat$true_exposure)
  result2 = get_result(test2, 
                       strand = FALSE, 
                       name = 'ANTs SRI weigthed',
                       grid = grid)
  
  return(rbind(result1, result2))
}

run_ANTS_SRI_Perm_model <- function(model_dat, grid) {
  # Unweigthed
  test1 = lm(strength ~ hair, data = model_dat[[1]]) # Run original model to get results
  result1 = get_result(test1, 
                       strand = FALSE, 
                       name = 'ANTs SRI NL',
                       grid = grid)
  
  ptest = stat.lm(model_dat, formula = strength ~ hair, progress = FALSE) # Run permuted models
  ptest = ANTs:::stat.p(c(ptest$Original.model$coefficients[2,1], ptest$permutations$hair))[3]
  result1$approach = paste(result1$approach, 'permuted', sep = " ")
  result1$`p-value` = ptest
  
  # Weigthed
  test2 = lm(strength ~ hair, data = model_dat[[1]],  weights = model_dat[[1]]$exposure) # Run original model to get results
  result2 = get_result(test2, 
                       strand = FALSE, 
                       name = 'ANTs SRI NL',
                       grid = grid)
  
  ptest = unlist(lapply(model_dat, function(x, exposure){# Run permuted models
    m = lm(data = x, formula = strength ~ hair, weights = exposure)
    summary(m)$coefficients[2,1]
  }, exposure = model_dat[[1]]$true_exposure))
  ptest = ANTs:::stat.p(ptest)[3]
  result2$approach = paste(result2$approach, 'permuted', sep = " ")
  result2$`p-value` = ptest
  
  return(rbind(result1, result2))
}

run_ANTS_SRI_PreNet_model <- function(model_dat, grid) {
  # Unweigthed
  ds = model_dat$perms
  # As we work with interactions some individual smay no be present
  d = data.frame('id' = colnames(model_dat$perms[[1]]))
  d2 = data.frame('id' = 1:nrow(model_dat$individual_predictors), "hair" =model_dat$individual_predictors[,1])
  d = merge.data.frame(d,d2, by = "id", all = T)

  perms = met.strength(ds, df = d, dfid = 'id')
  test1 = lm(strength ~ hair, data = perms[[1]]) # Run original model to get results

  result1 = get_result(test1, 
                       strand = FALSE, 
                       name = 'ANTs SRI prenet',
                       grid = grid)
  ptest = stat.lm(perms, formula = strength ~ hair, progress = FALSE, oda = model_dat$oda) # Run permuted models
  ptest = ANTs:::stat.p(c(ptest$Original.model$coefficients[2,1], ptest$permutations$hair))[3]
  result1$approach = paste(result1$approach, 'DS', sep = " ")
  result1$`p-value` = ptest
  
  # Weigthed
  test2 = lm(strength ~ hair, data = perms[[1]], weights = model_dat$exposure)
  result2 = get_result(test2, 
                       strand = FALSE, 
                       name = 'ANTs SRI prenet weigthed',
                       grid = grid)
  
  ptest = unlist(lapply(perms, function(x, exposure){# Run permuted models
    m = lm(data = x, formula = strength ~ hair, weights = exposure)
    summary(m)$coefficients[2,1]
  }, exposure = model_dat$exposure))
  ptest = ANTs:::stat.p(ptest)[3]
  result2$approach = paste(result2$approach, 'permuted', sep = " ")
  result2$`p-value` = ptest
  
  return(rbind(result1, result2))
}

run_BiSON_model <- function(model_dat, grid){
  model_dat$interactions$ego = as.factor(model_dat$interactions$focal)
  model_dat$interactions$alter = as.factor(model_dat$interactions$target)
  priors <- get_default_priors("count_conjugate")
  priors$edge <- "gamma(0.1, 0.1)"
  fit_edge <- bison_model(
    (outcome | duration) ~ dyad(focal, target),
    data=model_dat$interaction,
    model_type="count_conjugate",
    priors=priors
  )

  d = data.frame('id' = 1:nrow(model_dat$individual_predictors),
                 'char' = model_dat$individual_predictors[,1])
  

  fit_brm <- bison_brm(
    bison(node_strength(id)) ~ char,
    fit_edge,
    d
  )
  
  t = summary(fit_brm)
  t2 = p_map(fit_brm)

  result <- data.frame(t$fixed$Estimate[2], 
                      t$fixed$`l-95% CI`[2], 
                      t$fixed$`u-95% CI`[2],
                      grid,
                      t2$p_MAP[2],
                      t$fixed$Est.Error[2]
                      )

  rownames(result) <- NULL
  result$approach <- 'BISON'
  result$cor <- NA

  colnames(result) <- c("scale(hair)", "5 %", "95 %", "N_id", "tie_effect", "detect_effect", "sim", "p-value", "se",
                        "approach", 'prediction')
  return(result)
}


# generate_grid -------------
generate_grid <- function(Reps = 100,
                          N_id = seq(20, 100, by = 20),
                          tie_effect = seq(-2, 2, by = 0.5),
                          detect_effect = seq(-4, 4, by = 1)){
  # Create grid of paramters for fitting models
  grid <- expand.grid(N_id, tie_effect, detect_effect)
  if (nrow(grid) <= Reps) {
    grid_subsample <- grid
    ncores <- nrow(grid)
  } else {
    grid_subsample <- grid[sample(1:nrow(grid), Reps, replace = FALSE), ]
  }
  colnames(grid_subsample) <- c("N_id", "tie_effect", "detect_effect")
  return(grid_subsample)
}

# Run all models -------------
run_models <- function(model_names, model_dat, grid_params) {
  output <- NULL
  for (a in 1:length(model_names)) {
    print(paste("Working on model: ", model_names[a]))
    output[[a]] <- do.call(model_names[a], list(model_dat[[a]], grid_params))
  }

  output2 <- do.call(rbind, output)
  return(output2)
}


# Generate data and run all models on a single set of parameters -------------
simulation <- function(grid_params,...) {
  print('Generating data')
  opt <- list(...)
  opt$N_id = grid_params$N_id
  opt$individual_effects_sr = grid_params$tie_effect
  opt$exposure_effects = grid_params$detect_effect
  model_dat = do.call(make_all_data, opt)

  print('Runing models data')
  model_names <- c("run_STRAND_model", "run_ANTs_Rate_model", "run_ANTs_Rate_Perm_model", "run_ANTs_SRI_model", "run_ANTS_SRI_Perm_model", "run_ANTS_SRI_PreNet_model", "run_BiSON_model")
  output <- run_models(model_names, model_dat, grid_params)
  for(i in 1:length(opt)){
    output[names(opt)[i]] = paste(opt[[i]])
  }
  return(output)
}

# Test -----------------
# temp = generate_grid(1)
# tmp = simulation(
#     individual_effects=matrix(c(-1, temp$tie_effect[1]), ncol=1, nrow=2),
#     exposure_effects = c(-1, temp$detect_effect[1]),
#     dr_sigma  = 20,
#    )


# Parallelized simulations ----------------
parallel_sim <- function(Reps = 100,
                         N = seq(20, 100, by = 20),
                         tie_effect = seq(-2, 2, by = 0.5),
                         detect_effect = seq(-4, 4, by = 1),
                         censoring_effect = seq(-4, 4, by = 1),
                         ...) {
  # Setup grid ------------------
  grid_subsample <- generate_grid(Reps, N, tie_effect, detect_effect, censoring_effect)
  # Setup cores ------------------
  require(foreach)
  library(doSNOW)
  cl <- nrow(grid_subsample)
  cl <- makeCluster(cl, outfile="")
  on.exit(registerDoSEQ())
  r <- foreach(i = 1:nrow(grid_subsample)) %dopar% {
    source('1.Codes/1.load_packages.R')
    source('1.Codes/2.data_simulation.R')
    tmp = grid_subsample[i,]
    tmp$sim = i
    #opt <- list(...)
    simulation(
      N_id =  grid_subsample$N_id[i],
      individual_effects_sr = c(-1, grid_subsample$tie_effect[i]),
      exposure_effects = c(-1, grid_subsample$detect_effect[i]), 
      censoring_effects = grid_subsample$censoring_effect[i], 
      grid_params = tmp, ...
    )
  }
  print('multi proc done')
  print(r)
  # Close cores ------------------
  parallel::stopCluster(cl)
  gc()
  registerDoSEQ()

  # Stack results ------------------
  result <- do.call("rbind", r)

  result$z <- result$`scale(hair)` / (result$`95 %` - result$`5 %`)
  result$MOE <- result$z * result$se
  result$ci5 <- result$z - abs(result$MOE)
  result$ci95 <- result$z + abs(result$MOE)
  result$resid <- result$tie_effect - result$z
  return(result)
}


#R = parallel_sim(Reps = 5, 
#                 N = seq(20, 30, by = 5), 
#                 tie_effect = seq(-1, 1, by = 0.1),
#                 sr_sigma = c(1.4, 0.8),
#                 sr_rho = 0.5,
#                 # Dyadic parameters
#                 dr_sigma = 1.2,
#                 dr_rho = 0.8,
#                 # Observation bias parameters
#                 exposure_bias = TRUE,
#                 exposure_mu = 1.9,
#                 exposure_sigma = 2.9,
#                 exposure_max = 40,
#                 # Censoring bias parameters
#                 simulate_censoring = TRUE,
#                 censoring_mu = 0,
#                 censoring_sigma = 1.0,
#                 censoring_effects = 0)
#