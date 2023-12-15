#########################################################################################
#### Testing simulation sr and dyadic effects
###########################################################
source("1. Codes/3. Data simulation.R")
test.strand <- function(att = NULL,
                        N_id = 30,
                        B = NULL,
                        V = 1,
                        groups=NULL,
                        
                        sr_mu = c(0,0),
                        sr_sigma = c(1, 1),
                        sr_rho = 0.5,
                        
                        dr_mu = c(0,0),
                        dr_sigma = 1,
                        dr_rho = 0.75,
                        
                        individual_predictors = NULL,
                        dyadic_predictors = NULL,
                        individual_effects = NULL,
                        dyadic_effects = NULL,
                        exposure_predictors = NULL,
                        exposure_effects = NULL,
                        exposure_sigma = 1,
                        exposure_baseline = 50,
                        int_intercept = c(Inf,Inf),
                        int_slope = c(Inf,Inf),
                        simulate.interactions = FALSE){
  
  #if(is.null(att)){stop("Individuals attributes (att) can't be NULL")}
  data = simulate_sbm_plus_srm_network_with_measurement_bias(N_id = N_id,
                                                             B = B,
                                                             V = V,
                                                             groups=groups,
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
  if(is.null(att)){
    model_dat = make_strand_data(outcome = list(Grooming = data$network),
                                 individual_covariates = NULL,
                                 block_covariates = NULL,
                                 outcome_mode = "binomial",
                                 exposure = list(data$true_samps)
    )
    
    fit =  fit_social_relations_model(data=model_dat,
                                      focal_regression = ~ 1,
                                      target_regression = ~ 1,
                                      dyad_regression = ~  1,
                                      mode="mcmc",
                                      stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                  iter_warmup = 1000, iter_sampling = 1000,
                                                                  max_treedepth = NULL, adapt_delta = .98),
                                      return_predicted_network = TRUE)
  }else{
    model_dat = make_strand_data(outcome = list(Grooming = data$network),
                                 individual_covariates = data.frame(Hairy = Hairy),
                                 block_covariates = NULL,
                                 outcome_mode = "binomial",
                                 exposure = list(data$true_samps)
    )
    
    fit = fit_social_relations_model(data=model_dat,
                                     focal_regression = ~ Hairy,
                                     target_regression = ~ Hairy,
                                     dyad_regression = ~  1,
                                     mode="mcmc",
                                     stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                 iter_warmup = 1000, iter_sampling = 1000,
                                                                 max_treedepth = NULL, adapt_delta = .98),
                                     return_predicted_network = TRUE)
  }
  
  r = summarize_strand_results(fit)
  require(rstan)
  # Get rhat and n_eff
  tmp = read_stan_csv(fit$fit$output_files())
  tmp = monitor(extract(tmp, permuted = FALSE, inc_warmup = TRUE))
  
  # Get probability matrix
  
  est_net <- round(apply(r$samples$predicted_network_sample,2:3,mean ))
  m = matrix(0, ncol = N_id, nrow = N_id)
  for(a in 1:N_id){
    x = r$samples$predicted_network_sample[,,a]
    m[,a] = apply(x,2,mean)
  }
  m = m/data$true_samps
  
  return(list('arguments'=  as.list( sys.call()), 'data' = data, 'model_dat' = model_dat, 'fit' = fit,
              'summary' = r, 'diag' = tmp, 'matrix' = m))
}

#' @param var A list of variation of a single argument
#' @param list A list of arguments for test.strand
test.scenarios <- function(name = "",
                           var = list(dr_rho = seq(0.01, 0.8, length.out=10)),
                           list = list(att = NULL,
                                       N_id = 50,
                                       B = NULL,
                                       V = 1,
                                       groups=NULL,
                                       
                                       sr_mu = c(0,0),
                                       sr_sigma = c(1, 1),
                                       sr_rho = 0.5,
                                       
                                       dr_mu = c(0,0),
                                       dr_sigma = 1,
                                       
                                       individual_predictors = NULL,
                                       dyadic_predictors = NULL,
                                       individual_effects = NULL,
                                       dyadic_effects = NULL,
                                       exposure_predictors = NULL,
                                       exposure_effects = NULL,
                                       exposure_sigma = 1,
                                       exposure_baseline = 50,
                                       int_intercept = c(Inf,Inf),
                                       int_slope = c(Inf,Inf),
                                       simulate.interactions = FALSE)){
  if(length(var)>1){stop()}
  require(parallel)
  require(foreach)
  require(doParallel)
  cl <- length(var[[1]])
  registerDoParallel(cores=cl)
  on.exit(registerDoSEQ())
  result <- foreach(i = 1:length(var[[1]]), .export = ls(globalenv())) %dopar% {
    #result <- NULL
    #for(i in 1:length(var[[1]])){
    set.seed(1)
    if(names(var) %in% c('sr_mu', 'sr_sigma', 'dr_mu')){
      tmp = list(c(1.7*var[[1]][i], 1*var[[1]][i]))
    }else{ tmp = var[[1]][i]}
    
    if(names(var) %in% names(list)){
      if(is.list(tmp)){
        list[[which(names(list) %in% names(var))]] = tmp[[1]]
        L = list
      }else{
        list[[which(names(list) %in% names(var))]] = tmp
        L = list
      }
      
    }else{
      L = c(tmp, list)
      names(L) = c(names(var), names(list))
    }
    #result[[i]] = do.call('test.strand',L)
    r = do.call('test.strand',L)
    return(r)
  }
  
  result$tested = var
  save(result, file = paste('2. Results/Appendices/', names(var), name, '.Rdata',sep = "",  collapse = " "))
  return(result)
}

plot.function <- function(result){
  data =  NULL
  for(a in 1:(length(result)-1)){
    tmp = result[[a]]
    
    t1 = data.frame(Estimated = tmp$diag[1,1],
                    Simulated = 1,
                    name = "B")
    
    t2 = data.frame(Estimated = as.numeric(tmp$summary$summary$Median[1]),
                    Simulated = tmp$arguments$sr_sigma[1],
                    name = "sr_sigma1")
    
    t3 = data.frame(Estimated = as.numeric(tmp$summary$summary$Median[2]),
                    Simulated = tmp$arguments$sr_sigma[2],
                    name = "sr_sigma2")
    
    t4 = data.frame(Estimated = as.numeric(tmp$summary$summary$Median[4]),
                    Simulated = tmp$arguments$sr_rho,
                    name = "sr_rho")
    
    t5 = data.frame(Estimated = as.numeric(tmp$summary$summary$Median[3]),
                    Simulated = tmp$arguments$dr_sigma,
                    name = "dr_sigma")
    
    t6 = data.frame(Estimated = as.numeric(tmp$summary$summary$Median[5]),
                    Simulated =  tmp$arguments$dr_rho,
                    name = "dr_rho")
    
    t = rbind(t1, t2, t3, t4, t5, t6)
    t$sim = a
    
    
    data = rbind(data, t)
  }
  
  data$tested = names(result[[length(result)]])
  
  library(ggplot2)
  p = ggplot(data, aes(x = Simulated, y = Estimated, group = name))+
    geom_boxplot()+geom_point()+facet_grid(~name, scales="free" )+
    ggtitle(paste(names(result$tested)[1], ' variation from ', result$tested[[1]][1], ' to ', result$tested[[1]][length(result$tested[[1]])]))
  return(p)
}

test.scenarios(var = list(dr_sigma = seq(0.1, 5, length.out = 4)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           sr_rho = 0.5,
                           dr_mu = c(0,0),
                           dr_rho = 0.75)
)


test.scenarios(var = list(sr_rho = seq(0.01, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1,
                           sr_rho = 0)
)

test.scenarios(var = list(dr_rho = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1,
                           sr_rho = 0)
)


test.scenarios(var = list(sr_mu = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1,
                           sr_rho = 0)
)


test.scenarios(var = list(sr_sigma = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1,
                           sr_rho = 0)
)

test.scenarios(var = list(dr_mu = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1,
                           sr_rho = 0)
)


load("2. Results/Appendices/dr_mu.Rdata")
p.dr_mu = plot.function(result)
p.dr_mu
ggsave('2. Results/Appendices/p.dr_mu.png')

load("2. Results/Appendices/dr_rho.Rdata")
p.dr_rho = plot.function(result)
p.dr_rho
ggsave('2. Results/Appendices/p.dr_rho.png')

load("2. Results/Appendices/dr_sigma.Rdata")
p.dr_sigma = plot.function(result)
p.dr_sigma
ggsave('2. Results/Appendices/p.dr_sigma.png')

load("2. Results/Appendices/sr_mu.Rdata")
p.sr_mu = plot.function(result)
p.sr_mu
ggsave('2. Results/Appendices/p.sr_mu.png')

load("2. Results/Appendices/sr_rho.Rdata")
p.sr_rho = plot.function(result) 
p.sr_rho
ggsave('2. Results/Appendices/p.sr_rho.png')

load("2. Results/Appendices/sr_sigma.Rdata")
p.sr_sigma = plot.function(result)
p.sr_sigma
ggsave('2. Results/Appendices/p.sr_sigma.png')


#########################################################################################
#### Testing simulation sr and dyadic effects through interactions simulations
#########################################################################################
test.scenarios(var = list(dr_sigma = seq(0.1, 5, length.out=5)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           sr_rho = 0.5,
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           simulate.interactions = TRUE),
               name = " interactions"
)


test.scenarios(var = list(sr_rho = seq(0.01, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)

test.scenarios(var = list(dr_rho = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)


test.scenarios(var = list(sr_mu = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)


test.scenarios(var = list(sr_sigma = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)

test.scenarios(var = list(dr_mu = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)


load("2. Results/Appendices/dr_mu interactions.Rdata")
p.dr_mu_interactions = plot.function(result)
p.dr_mu_interactions
ggsave('2. Results/Appendices/p.dr_mu_interactions.png')

load("2. Results/Appendices/dr_rho interactions.Rdata")
p.dr_rho_interactions = plot.function(result)
p.dr_rho_interactions
ggsave('2. Results/Appendices/p.dr_rho_interactions.png')

load("2. Results/Appendices/dr_sigma interactions.Rdata")
p.dr_sigma_interactions = plot.function(result)
p.dr_sigma_interactions
ggsave('2. Results/Appendices/p.dr_sigma_interactions.png')

load("2. Results/Appendices/sr_mu interactions.Rdata")
p.sr_mu_interactions = plot.function(result)
p.sr_mu_interactions
ggsave('2. Results/Appendices/p.sr_mu_interactions.png')

load("2. Results/Appendices/sr_rho interactions.Rdata")
p.sr_rho_interactions = plot.function(result)
p.sr_rho_interactions
ggsave('2. Results/Appendices/p.sr_rho_interactions.png')

load("2. Results/Appendices/sr_sigma interactions.Rdata")
p.sr_sigma_interactions = plot.function(result)
p.sr_sigma_interactions
ggsave('2. Results/Appendices/p.sr_sigma_interactions.png')


#########################################################################################
#### Testing simulation sr and dyadic effects with bias
#########################################################################################
Hairy = matrix(rnorm(50, 0, 1), nrow=50, ncol=10)
test.scenarios(var = list(dr_sigma = seq(0.1, 5, length.out=2)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           sr_rho = 0.5,
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           exposure_predictors = cbind(rep(1,50),Hairy),
                           exposure_effects = c(-1, -4),# strong negative bais effect
                           exposure_sigma = 2.9,
                           exposure_baseline = 40),
               name = " with bias"
)


test.scenarios(var = list(sr_rho = seq(0.01, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           exposure_predictors = cbind(rep(1,50),Hairy),
                           exposure_effects = c(-1, -4),# strong negative bais effect
                           exposure_sigma = 2.9,
                           exposure_baseline = 40),
               name = " with bias"
)

test.scenarios(var = list(dr_rho = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           exposure_predictors = cbind(rep(1,50),Hairy),
                           exposure_effects = c(-1, -4),# strong negative bais effect
                           exposure_sigma = 2.9,
                           exposure_baseline = 40),
               name = " with bias"
)


test.scenarios(var = list(sr_mu = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           exposure_predictors = cbind(rep(1,50),Hairy),
                           exposure_effects = c(-1, -4),# strong negative bais effect
                           exposure_sigma = 2.9,
                           exposure_baseline = 40),
               name = " with bias"
)


test.scenarios(var = list(sr_sigma = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           exposure_predictors = cbind(rep(1,50),Hairy),
                           exposure_effects = c(-1, -4),# strong negative bais effect
                           exposure_sigma = 2.9,
                           exposure_baseline = 40),
               name = " with bias"
)

test.scenarios(var = list(dr_mu = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           exposure_predictors = cbind(rep(1,50),Hairy),
                           exposure_effects = c(-1, -4),# strong negative bais effect
                           exposure_sigma = 2.9,
                           exposure_baseline = 40),
               name = " with bias"
               
)


# A faire !!!!!!!!!!!!!!!!!!!!!!
load("2. Results/Appendices/dr_mu with bias.Rdata")
p.dr_mu_withBias = plot.function(result)
p.dr_mu_withBias
ggsave('2. Results/Appendices/p.dr_mu_withBias.png')

load("2. Results/Appendices/dr_rho with bias.Rdata")
p.dr_rho_withBias = plot.function(result)
p.dr_rho_withBias
ggsave('2. Results/Appendices/p.dr_rho_withBias.png')

load("2. Results/Appendices/dr_sigma with bias.Rdata")
p.dr_sigma_withBias = plot.function(result)
p.dr_sigma_withBias
ggsave('2. Results/Appendices/p.dr_sigma_withBias.png')

load("2. Results/Appendices/sr_mu with bias.Rdata")
p.sr_mu_withBias = plot.function(result)
p.sr_mu_withBias
ggsave('2. Results/Appendices/p.sr_mu_withBias.png')

load("2. Results/Appendices/sr_rho with bias.Rdata")
p.sr_rho_withBias = plot.function(result)
p.sr_rho_withBias
ggsave('2. Results/Appendices/p.sr_rho_withBias.png')

load("2. Results/Appendices/sr_sigma with bias.Rdata")
p.sr_sigma_withBias = plot.function(result)
p.sr_sigma_withBias
ggsave('2. Results/Appendices/p.sr_sigma_withBias.png')

#########################################################################################
#### Testing simulation sr and dyadic effects with bias through interactions simulations
#########################################################################################
Hairy = matrix(rnorm(50, 0, 1), nrow=50, ncol=1)
test.scenarios(var = list(dr_sigma = seq(0.1, 5, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           sr_rho = 0.5,
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           exposure_predictors = cbind(rep(1,50),Hairy),
                           exposure_effects = c(-1, -4),# strong negative bias effect
                           exposure_sigma = 2.9,
                           exposure_baseline = 40,
                           simulate.interactions = TRUE),
               name = " with bias interactions"
)


test.scenarios(var = list(sr_rho = seq(0.01, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           exposure_predictors = cbind(rep(1,50),Hairy),
                           exposure_effects = c(-1, -4),# strong negative bias effect
                           exposure_sigma = 2.9,
                           exposure_baseline = 40,
                           simulate.interactions = TRUE),
               name = " with bias interactions"
)

test.scenarios(var = list(dr_rho = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           exposure_predictors = cbind(rep(1,50),Hairy),
                           exposure_effects = c(-1, -4),# strong negative bias effect
                           exposure_sigma = 2.9,
                           exposure_baseline = 40,
                           simulate.interactions = TRUE),
               name = " with bias interactions"
)


test.scenarios(var = list(sr_mu = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           exposure_predictors = cbind(rep(1,50),Hairy),
                           exposure_effects = c(-1, -4),# strong negative bias effect
                           exposure_sigma = 2.9,
                           exposure_baseline = 40,
                           simulate.interactions = TRUE),
               name = " with bias interactions"
)


test.scenarios(var = list(sr_sigma = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           exposure_predictors = cbind(rep(1,50),Hairy),
                           exposure_effects = c(-1, -4),# strong negative bias effect
                           exposure_sigma = 2.9,
                           exposure_baseline = 40,
                           simulate.interactions = TRUE),
               name = " with bias interactions"
)

test.scenarios(var = list(dr_mu = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(2.1, 0.7),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1.4,
                           sr_rho = 0,
                           exposure_predictors = cbind(rep(1,50),Hairy),
                           exposure_effects = c(-1, -4),# strong negative bias effect
                           exposure_sigma = 2.9,
                           exposure_baseline = 40,
                           simulate.interactions = TRUE),
               name = " with bias interactions"
               
)

load("2. Results/Appendices/dr_muwith bias interactions.Rdata")
p.dr_mu_interactions_bias = plot.function(result)
p.dr_mu_interactions_bias
ggsave('2. Results/Appendices/p.dr_mu_interactions_bias.png')

load("2. Results/Appendices/dr_rhowith bias interactions.Rdata")
p.dr_rho_interactions_bias = plot.function(result)
p.dr_rho_interactions_bias
ggsave('2. Results/Appendices/p.dr_rho_interactions_bias.png')

load("2. Results/Appendices/dr_sigmawith bias interactions.Rdata")
p.dr_sigma_interactions_bias = plot.function(result)
p.dr_sigma_interactions_bias
ggsave('2. Results/Appendices/p.dr_sigma_interactions_bias.png')

load("2. Results/Appendices/sr_muwith bias interactions.Rdata")
p.sr_mu_interactions_bias = plot.function(result)
p.sr_mu_interactions_bias
ggsave('2. Results/Appendices/p.sr_mu_interactions_bias.png')

load("2. Results/Appendices/sr_rhowith bias interactions.Rdata")
p.sr_rho_interactions_bias = plot.function(result)
p.sr_rho_interactions_bias
ggsave('2. Results/Appendices/p.sr_rho_interactions_bias.png')

load("2. Results/Appendices/sr_sigmawith bias interactions.Rdata")
p.sr_sigma_interactions_bias = plot.function(result)
p.sr_sigma_interactions_bias
ggsave('2. Results/Appendices/p.sr_sigma_interactions_bias.png')



