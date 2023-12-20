# Functions -------
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
source("1.Codes/2.data_simulation.R")
test.function <- function(att = NULL,
                          
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
                           var = list(dr_rho = seq(0.5, 1, length.out=10)),
                           list = list(att = NULL,
                                       N_id = 50,
                                       B = NULL,
                                       V = 1,
                                       groups=NULL,
                                       
                                       sr_mu = c(0,0),
                                       sr_sigma = c(1, 1),
                                       sr_rho = 0,
                                       
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
  #registerDoParallel(cores=cl)
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
  save(result, file = paste('2.Results/Appendices/2/', names(var), name, '.Rdata',sep = "",  collapse = " "))
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


# 1. Testing simulation sr and dyadic effects------
test.scenarios(var = list(dr_sigma = seq(0.5, 1, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           sr_rho = 0,
                           dr_mu = c(0,0),
                           dr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)


test.scenarios(var = list(sr_rho = seq(0.5, 1, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0,
                           dr_sigma = 1,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)

test.scenarios(var = list(dr_rho = seq(0.5, 1, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0,
                           dr_sigma = 1,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)


test.scenarios(var = list(sr_mu = seq(0.5, 1, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)


test.scenarios(var = list(sr_sigma = seq(0.5, 1, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)

test.scenarios(var = list(dr_mu = seq(0.5, 1, length.out=10)),
               list = list(sr_mu = c(0,0),
                           sr_sigma  = c(1, 1),
                           dr_mu = c(0,0),
                           dr_rho = 0.75,
                           dr_sigma = 1,
                           sr_rho = 0,
                           simulate.interactions = TRUE),
               name = " interactions"
)


load("2.Results/Appendices/2/dr_mu interactions.Rdata")
p.dr_mu_interactions = plot.function(result)
p.dr_mu_interactions
ggsave('2.Results/Appendices/2/p.dr_mu_interactions.png')

load("2.Results/Appendices/2/dr_rho interactions.Rdata")
p.dr_rho_interactions = plot.function(result)
p.dr_rho_interactions
ggsave('2.Results/Appendices/2/p.dr_rho_interactions.png')

load("2.Results/Appendices/2/dr_sigma interactions.Rdata")
p.dr_sigma_interactions = plot.function(result)
p.dr_sigma_interactions
ggsave('2.Results/Appendices/2/p.dr_sigma_interactions.png')

load("2.Results/Appendices/2/sr_mu interactions.Rdata")
p.sr_mu_interactions = plot.function(result)
p.sr_mu_interactions
ggsave('2.Results/Appendices/2/p.sr_mu_interactions.png')

load("2.Results/Appendices/2/sr_rho interactions.Rdata")
p.sr_rho_interactions = plot.function(result)
p.sr_rho_interactions
ggsave('2.Results/Appendices/2/p.sr_rho_interactions.png')

load("2.Results/Appendices/2/sr_sigma interactions.Rdata")
p.sr_sigma_interactions = plot.function(result)
p.sr_sigma_interactions
ggsave('2.Results/Appendices/2/p.sr_sigma_interactions.png')
