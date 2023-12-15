# Functions -------
knitr::opts_chunk$set(warning = FALSE, message = FALSE)
#' This function uses the "3. Data simulation.R" function to generate simulated data. 
#' Once the data is generated, the function evaluates the relationship between individual characteristics and sociality, exposure, and censoring. 
#' We tested under different variations of sociality, exposure, and censoring to ensure the simulation generates accurate data, for example, when the individual_effects parameter is set to 0.
#' Additionally, we conducted tests with different values of the individual_effects parameter to identify when it generates non-significant relationships between individual characteristics and sociality. 
#' This helps in determining a threshold for false positive and false negative rates when testing different methods.
library(sjPlot)
library(ggplot2)
library(lmerTest)
library(ANTs)
set.seed(1)
source("1. Codes/3. Data simulation.R")
test.function <- function(att = NULL,
                          
                          N_id = 30,
                          B = NULL,
                          V = 1,
                          groups=NULL,
                          
                          sr_mu = c(0,0),
                          sr_sigma = c(1, 1),
                          sr_rho = 0,
                          
                          dr_mu = c(0,0),
                          dr_sigma = 1,
                          dr_rho = 0,
                          
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
                          simulate.interactions = FALSE,
                          test = TRUE,
                          print = TRUE,
                          legend = ''){
  require(ggplot2)
  require(ggpubr)
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

  if(test & !is.null(att)){
    if(simulate.interactions){
      df = cbind(Hairy, met.strength(data$network), 
                 met.strength(data$network/(1+data$true_samps)), data$true_exposure, tapply(data$interactions$s_i, data$interactions$sender  , mean))
      df = as.data.frame(df)
    }else{
      df = cbind(Hairy, met.strength(data$network),
                 met.strength(data$network/(1+data$true_samps)), data$true_exposure)
      df= as.data.frame(df)
      df$Censoring = 1
    }
    
    colnames(df) = c('att','Strength', 'Strength.corrected', 'Exposure', 'Censoring')

    p1 = ggplot(df, aes( x= att, y = Strength))+geom_point(aes(size = 1), alpha = 0.5, show.legend=c(size=FALSE, alpha = FALSE))+xlab("Individuals characteristics")+theme(text = element_text(size=13))+ labs(tag = "(a)")
    p2 = ggplot(df, aes( x= att, y = Strength.corrected))+geom_point(aes(size = 1), alpha = 0.5, show.legend =  c(size=FALSE, alpha = FALSE))+xlab("Individuals characteristics")+theme(text = element_text(size=13))+labs(tag = "(b)")
    p3 = ggplot(df, aes( x= att, y = Exposure))+geom_point(aes(size = 1), alpha = 0.5, show.legend =  c(size=FALSE, alpha = FALSE))+xlab("Individuals characteristics")+theme(text = element_text(size=13))+ labs(tag = "(c)")
    p4 = ggplot(df, aes( x= att, y = Censoring))+geom_point(aes(size = 1),alpha = 0.5,  show.legend =  c(size=FALSE, alpha = FALSE))+xlab("Individuals characteristics")+theme(text = element_text(size=13))+ labs(tag = "(d)")

    result = NULL
    test = lm(Strength~att, data = df)
    s = summary(test)
    estimate = s$coefficients[2,1]
    se = s$coefficients[2,2]
    sig = s$coefficients[2,4]
    result[[1]] = test
    if(print){
      cat("Relationship between individuals characteristics and strength none corrected---------------------------------", '\n')
      print(s)
    }
    
    test = lm(Strength.corrected~att, data = df)
    s = summary(test)
    estimate = s$coefficients[2,1]
    se = s$coefficients[2,2]
    sig = s$coefficients[2,4]
    result[[2]] = test
    if(print){
      cat("Relationship between individuals characteristics and strength corrected ---------------------------------", '\n')
      print(s)
    }
    
    test = lm(Strength.corrected~att, data = df, weights =Exposure)
    s = summary(test)
    estimate = s$coefficients[2,1]
    se = s$coefficients[2,2]
    sig = s$coefficients[2,4]
    result[[3]] = test
    if(print){
      cat("Relationship between individuals characteristics and strength corrected and lm with weigth---------------------------------", '\n')
      print(s)
    }

    test = lm(Exposure~att, data = df)
    s = summary(test)
    estimate = s$coefficients[2,1]
    se = s$coefficients[2,2]
    sig = s$coefficients[2,4]
    result[[4]] = test
    if(print){
      cat("Relationship between individuals characteristics and exposure ---------------------------------", '\n')
      print(s)
    }

    test = lm(Censoring~att, data = df)
    s = summary(test)
    estimate = s$coefficients[2,1]
    se = s$coefficients[2,2]
    sig = s$coefficients[2,4]
    result[[5]] = test
    if(print){
      cat("Relationship between individuals characteristics and censoring ---------------------------------", '\n')
      print(s)
    }
    
    names(result) = c('Strength', 'Strength.corrected', 'Strength.corrected.weigthed', 'Exposure', 'Censoring')

    return(list('data' = data, 'result' = result, 
                'plots' = annotate_figure(ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2), 
                                          bottom  = text_grob(legend, color = "black",
                                                             hjust = 1, x = 1, face = "italic", size = 10)
                                          )))
  }else{
    return(data)
  }
}

error.rates <- function(d, threshold = 0.15, legend = ''){
  # Rates of false negatives -------------
  tmp = d[d$effect >= threshold | d$effect <= -threshold,]
  if(nrow(tmp) != 0){
    t1 = unlist(lapply(split(tmp, tmp$approach), function(x){
      n = nrow(x)
      p = sum(x$p >= 0.05)
      return(p/n)
    }))
  }else{t1 = NULL}

  # Rates of false positives -------------
  tmp = d[d$effect <= threshold &  d$effect  >= -threshold,]
  if(nrow(tmp) != 0){
    t2 = unlist(lapply(split(tmp, tmp$approach), function(x){
      n = nrow(x)
      p = sum(x$p <= 0.05)
      return(p/n)
    }))
  }else{t2 = NULL}


  summary = data.frame(t1*100,t2*100, names(t1))
  colnames(summary) = c('false negatives', 'false positives', 'approaches')
  rownames(summary) = NULL
  
  p = ggplot(d[!d$approach %in% c('Censoring', 'Exposure'),], aes(x = effect, y = p))+geom_jitter()+
    geom_hline(yintercept = 0.05, linetype="dashed", color = "red")+
    geom_vline(xintercept = threshold, linetype="dashed", color = "blue")+
    facet_grid(~approach)+ylab('p-value')+xlab('individual_effects value') + labs(caption = legend) +
    theme(plot.caption = element_text(size = 12, hjust = 0, margin = margin(15,0,0,0)))
  print(p)
  return(list(summary, p))
}


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


# 1. Testing simulation sr and dyadic effects-----
test.scenarios(var = list(dr_sigma = seq(0.1, 5, length.out = 10)),
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


# 2. Testing simulation sr and dyadic effects through interactions simulations------
test.scenarios(var = list(dr_sigma = seq(0.1, 5, length.out=5)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           sr_rho = 0.5,
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           simulate.interactions = TRUE),
               name = " interactions"
)


test.scenarios(var = list(sr_rho = seq(0.01, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)

test.scenarios(var = list(dr_rho = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)


test.scenarios(var = list(sr_mu = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)


test.scenarios(var = list(sr_sigma = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)

test.scenarios(var = list(dr_mu = seq(0.1, 0.8, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1,
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

# 3. Testing individual characteristics on sociality, exposure, and censoring  -------
## 1.1. Individual characteristics do not impact sociality, exposure, or censoring
test1 = test.function(att = Hairy,
                      N_id = N_id,
                      individual_predictors=NULL, # individuals characteristics
                      individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                      sr_mu =  c(0.3, 1.5), sr_sigma =  c(1,1), # no sender-receiver effect
                      dr_mu = c(0,0), dr_sigma = 1, # no dyadic effect
                      exposure_predictors = NULL,
                      exposure_effects = c(0, 0), exposure_sigma = 1, # exposure effect
                      int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf),#no censoring effect
                      simulate.interactions = T,
                      legend = "Figure 1. No Relationship between individuals characteristics (a) sociality, (b) sociality corrected by exposure, (c) exposure, or (d) censoring.") 
test1$plots

#' The results of the regressions show, as expected, no significant effect in the relationship between individual characteristics, sociality, exposure, or censoring.
#'
## 1.2. There is a relationship between individual characteristics and sociality, but there is no relationship between individual characteristics, observation bias, and censoring
test2 = test.function(att = Hairy,
                      N_id = N_id,
                      individual_predictors=Hairy, # individuals characteristics
                      individual_effects=matrix(c(0.4,0.4),ncol=1, nrow=2), # individuals characteristics on interaction probability
                      sr_mu = c(0,0), sr_sigma =  c(1,1), # no sender-receiver effect
                      dr_mu = c(0,0), dr_sigma = 1, # no dyadic effect
                      exposure_predictors = NULL,
                      exposure_effects = NULL,exposure_sigma = 1, #no exposure effect
                      int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf), #no censoring effect
                      simulate.interactions = T, 
                      legend = "Figure 2. Relationship between individuals characteristics and (a) sociality, (b) sociality corrected by exposure, 
                      but no relationship betweenindividuals characteristics (c) exposure, (d) censoring, or (d) censoring.") 
test2$plots

#' The results of the regressions show, as expected, a significant effect in the relationship between individual characteristics and sociality, but no significant effect between  individuals characteristics exposure, and censoring.
#'
## 1.3. There is no relationship between individual characteristics,  sociality and censoring, but there is a relationship between individual characteristics and exposure
test3 = test.function(att = Hairy,
                      N_id = N_id,
                      individual_predictors=Hairy, # individuals characteristics
                      individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                      sr_mu =  c(0, 0), sr_sigma =  c(1,1), # no sender-receiver effect
                      dr_mu = c(0,0), dr_sigma = 1, # no dyadic effect
                      exposure_predictors = cbind(rep(1,N_id),Hairy),
                      exposure_effects = c(-1, 4), exposure_sigma = 1, # exposure effect
                      int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf),#no censoring effect
                        simulate.interactions = TRUE, 
                      legend = "Figure 3. No relationship between individuals characteristics and (a) sociality, (b) sociality corrected by exposure, (d) censoring, 
                      but precense of relationship betweenindividuals characteristics and (c) exposure.") 
test3$plots

#' The results of the regressions show, as expected, a significant effect in the relationship between individual characteristics and exposure 
#' which lead to a significant effect between individuals characteristics and (a) sociality and near significant effect between individuals characteristics and (b) correct ed sociality.


## 1.4. There is no relationship between individual characteristics, sociality and exposure but there is a relationship between individual characteristics and censoring
test4 = test.function(att = Hairy,
                      N_id = N_id,
                      individual_predictors=Hairy, # individuals characteristics
                      individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                      sr_mu =  c(0, 0), sr_sigma =  c(1,1), # no sender-receiver effect
                      dr_mu = c(0,0), dr_sigma = 1, # no dyadic effect
                      exposure_predictors = NULL,
                      exposure_effects = c(0, 0), exposure_sigma = 1, # exposure effect
                      int_intercept = c(0,0), int_slope = c(0.4,0.4),# censoring effect
                      simulate.interactions = T, 
                      legend = "Figure 5. No relationship between individuals characteristics and (a) sociality, (b) sociality corrected by exposure,
                      (c) exposure, but precense of relationship between individuals characteristics and (d) censoring.") 
test4$plots
#' 
#' The results of the regressions show, as expected, a significant effect in the relationship between individual characteristics and censoring 
#' which lead to a significant effect between individuals characteristics, (a) sociality and (b) correct ed sociality.
#' 



# 4. Testing when the coefficient of individual characteristics (individual_effects parameter) results in a significant effect on simulated data -------
N_id = 30
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
TEST = seq(from = 0, to = 0.6, by = 0.01)
length(TEST)
r = NULL
a = 1
for (a in a:length(TEST)) {
  for(b in 1:10){
    r[[length(r)+1]] = test.function(att = Hairy,
                           N_id = N_id,
                           individual_predictors=Hairy, # individuals characteristics
                           individual_effects=matrix(c(TEST[a],TEST[a]),ncol=1, nrow=2), # individuals characteristics on interaction probability
                           sr_mu =  c(0, 0), sr_sigma =  c(1,1), # no sender-receiver effect
                           dr_mu = c(0,0), dr_sigma = 1, # no dyadic effect
                           exposure_predictors = NULL,
                           exposure_effects = c(0, 0), exposure_sigma = 1, # exposure effect
                           int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf),#no censoring effect
                           simulate.interactions = TRUE, print = FALSE) 
  }
}
d = NULL
test = rep(TEST, each = 10)
for(a in 1: length(r)){
  for (b in 1:length(r[[a]]$result)) {
    s = summary(r[[a]]$result[[b]])
    p = s$coefficients[2,4]
    c = s$coefficients[2,1]
    d = rbind(d, data.frame('coef' = c, 'p' = p, 'effect' = test[a], 'approach' = names(r[[a]]$result)[b], 'sim' = a))
  }
}

error.rates(d, threshold = 0.20)

#'
#' From a visual perspective and error rates we can see that bellow a value of 0.20 for individual_effects parameters, we obtain no or or null effects.
#' We will use values of individual_effects ranging from 0 to 0.19 for simulations without sociality effect and values ranging 0.2 to 0.4 for simulations with sociality effect.
#'

## 2.1. An example of individual_effects being equal to 0.2 in simulated data
N_id = 50
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
test = test.function(att = Hairy,
                      N_id = N_id,
                      individual_predictors=Hairy, # individuals characteristics
                      individual_effects=matrix(c(0.2,0.2),ncol=1, nrow=2), # individuals characteristics on interaction probability
                      sr_mu =  c(0, 0), sr_sigma =  c(1,1), # no sender-receiver effect
                      dr_mu = c(0,0), dr_sigma = 1, # no dyadic effect
                      exposure_predictors = NULL,
                      exposure_effects = c(0, 0), exposure_sigma = 1, # exposure effect
                      int_intercept = c(Inf,Inf), int_slope = c(-Inf,-Inf),
                      simulate.interactions = T) #no censoring effect
test$plots

## 2.2. An example of individual_effects being equal to 0.4 in simulated data
N_id = 50
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
test = test.function(att = Hairy,
                     N_id = N_id,
                     individual_predictors=Hairy, # individuals characteristics
                     individual_effects=matrix(c(0.4,0.4),ncol=1, nrow=2), # individuals characteristics on interaction probability
                     sr_mu =  c(0, 0), sr_sigma =  c(1,1), # no sender-receiver effect
                     dr_mu = c(0,0), dr_sigma = 1, # no dyadic effect
                     exposure_predictors = NULL,
                     exposure_effects = c(0, 0), exposure_sigma = 1, # exposure effect
                     int_intercept = c(Inf,Inf), int_slope = c(Inf, Inf),
                     simulate.interactions = T) #no censoring effect
test$plots


# 5. Testing when the coefficient of exposure (exposure_effects parameter) lead to significant effect on simulated data -------
N_id = 30
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
TEST = seq(from = 0, to = 0.6, by = 0.01)
length(TEST)
r = NULL
a = 1
for (a in a:length(TEST)) {
  for(b in 1:10){
    r[[length(r)+1]] = test.function(att = Hairy,
                                     N_id = N_id,
                                     individual_predictors=Hairy, # individuals characteristics
                                     individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                                     sr_mu =  c(0, 0), sr_sigma =  c(1,1), # no sender-receiver effect
                                     dr_mu = c(0,0), dr_sigma = 1, # no dyadic effect
                                     exposure_predictors = cbind(rep(1,N_id),Hairy),
                                     exposure_effects = c(-1, TEST[a]), exposure_sigma = 1, # exposure effect
                                     int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf),#no censoring effect
                                     simulate.interactions = TRUE, print = FALSE) 
  }
}
d = NULL
test = rep(TEST, each = 10)
for(a in 1: length(r)){
  for (b in 1:length(r[[a]]$result)) {
    s = summary(r[[a]]$result[[b]])
    p = s$coefficients[2,4]
    c = s$coefficients[2,1]
    d = rbind(d, data.frame('coef' = c, 'p' = p, 'effect' = test[a], 'approach' = names(r[[a]]$result)[b], 'sim' = a))
  }
}

error.rates(d, threshold = 0.30)

#'
#' From a visual perspective and error rates we can see that above a value of 0.30 for individual_effects parameters, we start to observe increase of false positive.
#' We will use values of exposure_effects ranging from 0 to 0.20 for simulations without exposure bias and values ranging 0.4 to 0.6 for simulations with exposure bias.
#'

## 2.1. An example of exposure_effects being equal to 0.2 in simulated data
N_id = 50
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
test = test.function(att = Hairy,
                     N_id = N_id,
                     individual_predictors=Hairy, # individuals characteristics
                     individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                     sr_mu =  c(0, 0), sr_sigma =  c(1,1), # no sender-receiver effect
                     dr_mu = c(0,0), dr_sigma = 1, # no dyadic effect
                     exposure_predictors = cbind(rep(1,N_id),Hairy),
                     exposure_effects = c(-1, 0.2), exposure_sigma = 1, # exposure effect
                     int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf),
                     simulate.interactions = T) #no censoring effect
test$plots

## 2.2. An example of individual_effects being equal to 0.4 in simulated data
test = test.function(att = Hairy,
                     N_id = N_id,
                     individual_predictors=Hairy, # individuals characteristics
                     individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                     sr_mu =  c(0, 0), sr_sigma =  c(1,1), # no sender-receiver effect
                     dr_mu = c(0,0), dr_sigma = 1, # no dyadic effect
                     exposure_predictors = cbind(rep(1,N_id),Hairy),
                     exposure_effects = c(-1, 0.4), exposure_sigma = 1, # exposure effect
                     int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf),
                     simulate.interactions = T) #no censoring effect
test$plots

# 6. Testing when the coefficient of censoring (int_slope parameter) lead to significant effect on simulated data -------
N_id = 30
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
TEST = seq(from = 0, to = 0.2, by = 0.003)
length(TEST)
r = NULL
a = 1
for (a in a:length(TEST)) {
  for(b in 1:10){
    r[[length(r)+1]] = test.function(att = Hairy,
                                     N_id = N_id,
                                     individual_predictors=Hairy, # individuals characteristics
                                     individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                                     sr_mu =  c(0, 0), sr_sigma =  c(1,1), # no sender-receiver effect
                                     dr_mu = c(0,0), dr_sigma = 1, # no dyadic effect
                                     exposure_predictors = NULL,
                                     exposure_effects = c(0, 0), exposure_sigma = 1, # exposure effect
                                     int_intercept = c(TEST[a],TEST[a]), int_slope = c(TEST[a],TEST[a]),#no censoring effect
                                     simulate.interactions = TRUE, print = FALSE) 
  }
}
d = NULL
test = rep(TEST, each = 10)
for(a in 1: length(r)){
  for (b in 1:length(r[[a]]$result)) {
    s = summary(r[[a]]$result[[b]])
    p = s$coefficients[2,4]
    c = s$coefficients[2,1]
    d = rbind(d, data.frame('coef' = c, 'p' = p, 'effect' = test[a], 'approach' = names(r[[a]]$result)[b], 'sim' = a))
  }
}

error.rates(d, threshold = 0.10)

#'
#' From a visual perspective and error rates we can see that above a value of 0.30 for individual_effects parameters, we start to observe increase of false positive.
#' We will use values of exposure_effects ranging from 0 to 0.20 for simulations without exposure bias and values ranging 0.4 to 0.6 for simulations with exposure bias.
#'

## 2.1. An example of censoring intercept and slope are equal to 0.1 in simulated data
N_id = 50
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
test = test.function(att = Hairy,
                     N_id = N_id,
                     individual_predictors=Hairy, # individuals characteristics
                     individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                     sr_mu =  c(0, 0), sr_sigma =  c(1,1), # no sender-receiver effect
                     dr_mu = c(0,0), dr_sigma = 1, # no dyadic effect
                     exposure_predictors = NULL,
                     exposure_effects = c(0, 0), exposure_sigma = 1, # exposure effect
                     int_intercept = c(0.1,0.1), int_slope = c(0.1,0.1),
                     simulate.interactions = T) #no censoring effect
test$plots

## 2.2. An example of individual_effects being equal to 0.4 in simulated data
test = test.function(att = Hairy,
                     N_id = N_id,
                     individual_predictors=Hairy, # individuals characteristics
                     individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                     sr_mu =  c(0, 0), sr_sigma =  c(1,1), # no sender-receiver effect
                     dr_mu = c(0,0), dr_sigma = 1, # no dyadic effect
                     exposure_predictors = NULL,
                     exposure_effects = c(0, 0), exposure_sigma = 1, # exposure effect
                     int_intercept = c(0.2,0.2), int_slope = c(0.2,0.2),
                     simulate.interactions = T) #no censoring effect
test$plots


save.image(file='2. Results/Appendices/Appendix.RData')
