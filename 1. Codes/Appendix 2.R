# Testing simulation approach
test.function <- function(att = NULL,
                          N_id = 30,
                          B = NULL,
                          V = 1,
                          groups=NULL,
                          sr_mu = c(0,0),
                          sr_sigma = c(0, 0),
                          sr_rho = 0,
                          dr_mu = c(0,0),
                          dr_sigma = 0,
                          dr_rho = 0,
                          individual_predictors = NULL,
                          dyadic_predictors = NULL,
                          individual_effects = NULL,
                          dyadic_effects = NULL,
                          exposure_predictors = NULL,
                          exposure_effects = NULL,
                          exposure_sigma = 0,
                          exposure_baseline = 50,
                          int_intercept = c(0,0),
                          int_slope = c(0,0),
                          simulate.interactions = TRUE,
                          test = TRUE){
  require(ggplot2)
  require(ggpubr)
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



  if(test & !is.null(att)){
    df = cbind(Hairy, met.strength(data$network), data$true_exposure, tapply(data$interactions$s_i , data$interactions$ego, mean))
    colnames(df) = c('att', 'Strength', 'Exposure', 'Censoring')
    df = as.data.frame(df)

    p1 = ggplot(df, aes( x= att, y = Strength))+geom_point(aes(size = 1), alpha = 0.5, show.legend = F)+xlab("Individuals characteristics")+theme(text = element_text(size=13))
    p2 = ggplot(df, aes( x= att, y = Exposure))+geom_point(aes(size = 1), alpha = 0.5, show.legend = F)+xlab("Individuals characteristics")+theme(text = element_text(size=13))
    p3 = ggplot(df, aes( x= att, y = Censoring))+geom_point(aes(size = 1),alpha = 0.5,  show.legend = F)+xlab("Individuals characteristics")+theme(text = element_text(size=13))

    result = NULL
    cat("Relationship between individuals characteristics and strength ---------------------------------", '\n')
    test = lm(Strength~att, data = df)
    s = summary(test)
    estimate = s$coefficients[2,1]
    se = s$coefficients[2,2]
    sig = s$coefficients[2,4]
    result[[1]] = test
    print(s)

    cat("Relationship between individuals characteristics and exposure ---------------------------------", '\n')
    test = lm(Exposure~att, data = df)
    s = summary(test)
    estimate = s$coefficients[2,1]
    se = s$coefficients[2,2]
    sig = s$coefficients[2,4]
    result[[2]] = test

    print(s)

    cat("Relationship between individuals characteristics and censoring ---------------------------------", '\n')
    test = lm(Censoring~att, data = df)
    s = summary(test)
    estimate = s$coefficients[2,1]
    se = s$coefficients[2,2]
    sig = s$coefficients[2,4]
    result[[3]] = test

    print(s)

    return(list('data' = data, 'result' = result, 'plots' = ggarrange(p1, p2, p3, ncol = 3, nrow = 1)))
  }else{
    return(data)
  }
}

# No relationship between individuals characteristics and sociality, observation bias and interaction bias

test.strand <- function(att = NULL,
                        N_id = 30,
                        B = NULL,
                        V = 1,
                        groups=NULL,

                        sr_mu = c(0,0),
                        sr_sigma = c(2.1, 0.7),
                        sr_rho = 0.5,

                        dr_mu = c(0,0),
                        dr_sigma = 1.4,
                        dr_rho = 0.75,

                        individual_predictors = NULL,
                        dyadic_predictors = NULL,
                        individual_effects = NULL,
                        dyadic_effects = NULL,
                        exposure_predictors = NULL,
                        exposure_effects = NULL,
                        exposure_sigma = 0,
                        exposure_baseline = 50,
                        int_intercept = c(0,0),
                        int_slope = c(0,0),
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
  if(is.null(att)){Hairy = matrix(1, nrow=N_id, ncol=1)} # If no attributes are declared, generate identical values for everyone.
  model_dat = make_strand_data(outcome = list(Grooming = data$network),
                               individual_covariates = data.frame(Hairy = Hairy),
                               block_covariates = NULL,
                               outcome_mode = "binomial",
                               exposure = list(data$true_samps)
  )
  fit_latent_network_model
  fit =  fit_social_relations_model(data=model_dat,
                                    focal_regression = ~ Hairy,
                                    target_regression = ~ Hairy,
                                    dyad_regression = ~  1,
                                    mode="mcmc",
                                    stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                iter_warmup = 1000, iter_sampling = 1000,
                                                                max_treedepth = NULL, adapt_delta = .98),
                                    return_predicted_network = TRUE)
  r = summarize_strand_results(fit)

  # Get rhat and n_eff
  tmp = rstan::read_stan_csv(fit$fit$output_files())
  tmp
  rhat = tmp@.MISC$summary$rhat[,1]
  rhat
  n_eff = tmp@.MISC$summary$ess[,1]
  n_eff
  r$summary$rhat = mean(rhat, na.rm= TRUE)
  r$summary$n_eff =  mean(n_eff, na.rm= TRUE)

  # Get probability matrix
  m = matrix(0, ncol = N_id, nrow = N_id)
  for(a in 1:N_id){
    x = r$samples$predicted_network_sample[,,a]
    m[,a] = apply(x,2,mean)
  }
  m = m/data$true_samps

  return(list('strand' = r, 'summary' = r$summary, 'matrix' = m))
}

N_id = 10
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
range = seq(from = 0.1, to = 1, by = 0.2)
t.s = list()
for(a in range){
  t.s[[length(t.s)+1]] = test.strand(att = NULL, N_id = N_id,
                                     sr_rho = 0, sr_mu = c(a,0), sr_sigma = c(0,0),
                                     simulate.interactions = FALSE)
}
for(a in 1:length(t.s)){
  cat(range[a], '---------------------------------------','\n')
  print(t.s[[a]]$strand$summary)
}


t.s2 = list()
for(a in range){
  t.s2[[length(t.s2)+1]] = test.strand(att = NULL, N_id = N_id,
                                       sr_rho = 0, sr_mu = c(0,a), sr_sigma = c(0,0),
                                       simulate.interactions = FALSE)
}
for(a in 1:length(t.s2)){
  cat(range[a], '---------------------------------------','\n')
  print(t.s2[[a]]$strand$summary)
}


t.srho = list()
for(a in range){
  t.srho[[length(t.srho)+1]] = test.strand(att = NULL, N_id = N_id,
                                           sr_rho = a, sr_mu = c(0,0), sr_sigma = c(0,0),
                                           simulate.interactions = FALSE)
}
for(a in 1:length(t.s2)){
  cat(range[a], '---------------------------------------','\n')
  print(t.srho[[a]]$strand$summary)
}

t.s.sigma = list()
for(a in range){
  t.s.sigma[[length(t.s.sigma)+1]] = test.strand(att = NULL, N_id = N_id,
                                                 sr_rho = 0, sr_mu = c(0,0), sr_sigma = c(a,0),
                                                 simulate.interactions = FALSE)
}
for(a in 1:length(t.s2)){
  cat(range[a], '---------------------------------------','\n')
  print(t.s.sigma[[a]]$strand$summary)
}

t.s.sigma2 = list()
for(a in range){
  t.s.sigma2[[length(t.s.sigma2)+1]] = test.strand(att = NULL, N_id = N_id,
                                                   sr_rho = 0, sr_mu = c(0,0), sr_sigma = c(0,a),
                                                   simulate.interactions = FALSE)
}
for(a in 1:length(t.s2)){
  cat(range[a], '---------------------------------------','\n')
  print(t.s.sigma2[[a]]$strand$summary)
}

# -------------------------------------------------------------------------------------------------------
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
test3 = test.function(att = Hairy,
                      N_id = N_id,
                      individual_predictors=NULL, # individuals characteristics
                      individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                      sr_mu =  c(0.3, 1.5), sr_sigma =  c(0,0), # no sender-receiver effect
                      dr_mu = c(0,0), dr_sigma = 0, # no dyadic effect
                      exposure_predictors = NULL,
                      exposure_effects = c(0, 0), exposure_sigma = 0, # exposure effect
                      int_intercept = c(Inf,Inf), int_slope = c(-Inf,-Inf),
                      return_predicted_network = TRUE) #no censoring effect
test3$plots
test1 = test.function(att = Hairy,
                      N_id = N_id,
                      individual_predictors=Hairy, # individuals characteristics
                      individual_effects=matrix(c(0,0),ncol=1, nrow=2), #no individuals characteristics on interaction probability
                      sr_mu = c(0,0), sr_sigma =  c(0,0), # no sender-receiver effect
                      dr_mu = c(0,0), dr_sigma = 0, # no dyadic effect
                      exposure_predictors = cbind(rep(Inf,N_id),rep(0, nrow(Hairy))),
                      exposure_effects = NULL,exposure_sigma = 0, #no exposure effect
                      int_intercept = c(Inf,Inf), int_slope = c(-Inf,-Inf)) #no censoring effect
test1$plots
tab1 = tab_model(test1$result[[1]], test1$result[[2]], test1$result[[3]], show.se = T, show.ci = F, file = "tab1.xls")

# Relationship between individuals characteristics and sociality but no relationship between individuals characteristics observation bias and interaction bias
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
test2 = test.function(att = Hairy,
                      N_id = N_id,
                      individual_predictors=Hairy, # individuals characteristics
                      individual_effects=matrix(c(0.4,0.4),ncol=1, nrow=2), # individuals characteristics on interaction probability
                      sr_mu = c(0,0), sr_sigma =  c(0,0), # no sender-receiver effect
                      dr_mu = c(0,0), dr_sigma = 0, # no dyadic effect
                      exposure_predictors = cbind(rep(1,N_id),rep(0, nrow(Hairy))),
                      exposure_effects = NULL,exposure_sigma = 0, #no exposure effect
                      int_intercept = c(Inf,Inf), int_slope = c(-Inf,-Inf)) #no censoring effect
test2$plots

# Relationship between individuals characteristics and sociality and relationship between individuals characteristics and exposure
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
test3 = test.function(att = Hairy,
                      N_id = N_id,
                      individual_predictors=NULL, # individuals characteristics
                      individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                      sr_mu =  c(0.3, 1.5), sr_sigma =  c(0,0), # no sender-receiver effect
                      dr_mu = c(0,0), dr_sigma = 0, # no dyadic effect
                      exposure_predictors = NULL,
                      exposure_effects = c(0, 0), exposure_sigma = 0, # exposure effect
                      int_intercept = c(Inf,Inf), int_slope = c(-Inf,-Inf)) #no censoring effect
test3$plots


# ggsave(".\3. Manuscript\Figures\Appendix 2\3.Social and observation relationship.png")
test3$result
tab3 = tab_model(test3$result[[1]], test3$result[[2]], test3$result[[3]], show.se = T, show.ci = F, file = "tab3.html")
tab3
webshot("tab3.html", "3. Manuscript/Figures/Appendix2/tab3.png")



# Relationship between individuals characteristics and censoring
N_id = 10
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
test4 = test.function(att = Hairy,
                      N_id = N_id, V = 1,
                      individual_predictors=Hairy, # individuals characteristics
                      individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                      sr_mu = c(0,0), sr_sigma =  c(0,0), # no sender-receiver effect
                      dr_mu = c(0,0), dr_sigma = 0, # no dyadic effect
                      exposure_predictors = cbind(rep(1,N_id),rep(0, nrow(Hairy))),
                      exposure_effects = NULL,exposure_sigma = 0, #no exposure effect
                      int_intercept = c(0,0), int_slope = c(0,0)) # censoring effect on focal individual
test4$plots
# ggsave(".\3. Manuscript\Figures\Appendix 2\4.social and interaction relationship.png")
test4$result
tab4 = tab_model(test4$result[[1]], test4$result[[2]], test4$result[[3]], show.se = T, show.ci = F, file = "tab4.html")
tab4
webshot("tab3.html", "3. Manuscript/Figures/Appendix2/tab4.png")

ggarrange(test1$plots, test2$plots, test3$plots, test4$plots, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave('Figure3.png', path = '3. Manuscript/Figures/Appendix2/')
