#' # 1. Run parallel functions: to assess simulation parameter accuracy.
options(warn = -1)
source("1.Codes/Appendix 2 parallel functions.R")
load("2.Results/Appendices/2/dr_rho interactions.Rdata")
p.dr_rho_interactions = plot.function(result)

load("2.Results/Appendices/2/dr_sigma interactions.Rdata")
p.dr_sigma_interactions = plot.function(result)

load("2.Results/Appendices/2/sr_rho interactions.Rdata")
p.sr_rho_interactions = plot.function(result)

load("2.Results/Appendices/2/sr_sigma interactions.Rdata")
p.sr_sigma_interactions = plot.function(result)

ggpubr::ggarrange(p.dr_rho_interactions, p.dr_sigma_interactions, p.sr_rho_interactions, p.sr_sigma_interactions[[1]], p.sr_sigma_interactions[[2]],  ncol = 3, nrow = 2)

#' # 2. Testing individual characteristics on sociality, exposure, and censoring  -------
N_id = 50
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
## Run parallel functions: to assess simulation parameter accuracy..1. Individual characteristics do not impact sociality, exposure, or censoring
test1 = test.function(att = Hairy,
                      N_id = N_id,
                      individual_predictors=NULL, # individuals characteristics
                      individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                      exposure_predictors = NULL,
                      exposure_effects = c(0, 0), exposure_sigma = 1, # exposure effect
                      int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf),#no censoring effect
                      simulate.interactions = T,
                      legend = "Figure 1. No Relationship between individuals characteristics (a) sociality, (b) sociality corrected by exposure, (c) exposure, or (d) censoring.") 
test1$plots

#' The results of the regressions show, as expected, no significant effect in the relationship between individual characteristics, sociality, exposure, or censoring.
#'
##' ## 2.1. There is a relationship between individual characteristics and sociality, but there is no relationship between individual characteristics, observation bias, and censoring
test2 = test.function(att = Hairy,
                      N_id = N_id,
                      individual_predictors=Hairy, # individuals characteristics
                      individual_effects=matrix(c(0.4,0.4),ncol=1, nrow=2), # individuals characteristics on interaction probability
                      exposure_predictors = NULL,
                      exposure_effects = NULL,exposure_sigma = 1, #no exposure effect
                      int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf), #no censoring effect
                      simulate.interactions = T, 
                      legend = "Figure 2. Relationship between individuals characteristics and (a) sociality, (b) sociality corrected by exposure, 
                      but no relationship betweenindividuals characteristics (c) exposure, (d) censoring, or (d) censoring.") 
test2$plots

#' The results of the regressions show, as expected, a significant effect in the relationship between individual characteristics and sociality, but no significant effect between  individuals characteristics exposure, and censoring.
#'
##' ## 2.2. There is no relationship between individual characteristics,  sociality and censoring, but there is a relationship between individual characteristics and exposure
test3 = test.function(att = Hairy,
                      N_id = N_id,
                      individual_predictors=Hairy, # individuals characteristics
                      individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                      exposure_predictors = cbind(rep(1,N_id),Hairy),
                      exposure_effects = c(-1, 4), exposure_sigma = 1, # exposure effect
                      int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf),#no censoring effect
                      simulate.interactions = TRUE, 
                      legend = "Figure 3. No relationship between individuals characteristics and (a) sociality, (b) sociality corrected by exposure, (d) censoring, 
                      but precense of relationship betweenindividuals characteristics and (c) exposure.") 
test3$plots

#' The results of the regressions show, as expected, a significant effect in the relationship between individual characteristics and exposure 
#' which lead to a significant effect between individuals characteristics and (a) sociality and near significant effect between individuals characteristics and (b) correct ed sociality.


##' ## 2.3. There is no relationship between individual characteristics, sociality and exposure but there is a relationship between individual characteristics and censoring
test4 = test.function(att = Hairy,
                      N_id = N_id,
                      individual_predictors=Hairy, # individuals characteristics
                      individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
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



#' # 3. Testing when the coefficient of individual characteristics (individual_effects parameter) results in a significant effect on simulated data
N_id = 30
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
TEST = seq(from = 0, to = 0.5, by = 0.05)
length(TEST)
r = NULL
a = 1
for (a in a:length(TEST)) {
  for(b in 1:10){
    r[[length(r)+1]] = test.function(att = Hairy,
                                     N_id = N_id,
                                     individual_predictors=Hairy, # individuals characteristics
                                     individual_effects=matrix(c(TEST[a],TEST[a]),ncol=1, nrow=2), # individuals characteristics on interaction probability
                                     exposure_predictors = NULL,
                                     exposure_effects = c(0, 0), exposure_sigma = 0.5, # exposure effect
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

error.rates(d, threshold = 0.15)

#'
#' From a visual perspective and error rates we can see that bellow a value of 0.20 for individual_effects parameters, we obtain no or or null effects.
#' We will use values of individual_effects ranging from 0 to 0.19 for simulations without sociality effect and values ranging 0.2 to 0.4 for simulations with sociality effect.
#'

##' ## 3.1. An example of individual_effects being equal to 0.2 in simulated data
N_id = 50
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
test = test.function(att = Hairy,
                     N_id = N_id,
                     individual_predictors=Hairy, # individuals characteristics
                     individual_effects=matrix(c(0.19,0.19),ncol=1, nrow=2), # individuals characteristics on interaction probability
                     exposure_predictors = NULL,
                     exposure_effects = c(0, 0), exposure_sigma = 1, # exposure effect
                     int_intercept = c(Inf,Inf), int_slope = c(-Inf,-Inf),
                     simulate.interactions = T) #no censoring effect
test$plots

##' ## 3.2. An example of individual_effects being equal to 0.4 in simulated data
N_id = 50
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
test = test.function(att = Hairy,
                     N_id = N_id,
                     individual_predictors=Hairy, # individuals characteristics
                     individual_effects=matrix(c(0.4,0.4),ncol=1, nrow=2), # individuals characteristics on interaction probability
                     exposure_predictors = NULL,
                     exposure_effects = c(0, 0), exposure_sigma = 1, # exposure effect
                     int_intercept = c(Inf,Inf), int_slope = c(Inf, Inf),
                     simulate.interactions = T) #no censoring effect
test$plots


#' # 4. Testing when the coefficient of exposure (exposure_effects parameter) lead to significant effect on simulated data
N_id = 30
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
TEST = seq(from = 0, to = 5, by = 0.2)
length(TEST)
r = NULL
a = 1
for (a in a:length(TEST)) {
  for(b in 1:10){
    r[[length(r)+1]] = test.function(att = Hairy,
                                     N_id = N_id,
                                     individual_predictors=Hairy, # individuals characteristics
                                     individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
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

error.rates(d, threshold =  1)

#'
#' From a visual perspective and error rates we can see that above a value of 0.30 for individual_effects parameters, we start to observe increase of false positive.
#' We will use values of exposure_effects ranging from 0 to 0.20 for simulations without exposure bias and values ranging 0.4 to 0.6 for simulations with exposure bias.
#'

#' ## 4.1. An example of exposure_effects being equal to 0.2 in simulated data
N_id = 50
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
test = test.function(att = Hairy,
                     N_id = N_id,
                     individual_predictors=Hairy, # individuals characteristics
                     individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                     exposure_predictors = cbind(rep(1,N_id),Hairy),
                     exposure_effects = c(-1, 0.2), exposure_sigma = 1, # exposure effect
                     int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf),
                     simulate.interactions = T) #no censoring effect
test$plots

#' # 4.2. An example of individual_effects being equal to 0.4 in simulated data
test = test.function(att = Hairy,
                     N_id = N_id,
                     individual_predictors=Hairy, # individuals characteristics
                     individual_effects=matrix(c(0,0),ncol=1, nrow=2), # individuals characteristics on interaction probability
                     exposure_predictors = cbind(rep(1,N_id),Hairy),
                     exposure_effects = c(1,  1), exposure_sigma = 1, # exposure effect
                     int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf),
                     simulate.interactions = T) #no censoring effect
test$plots

#' Sociality patterns observed in plot (a) are only due to exposure bias (plot (c)).

#' # 5. Testing when the coefficient of censoring (int_slope parameter) lead to significant effect on simulated data
N_id = 30
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
TEST = seq(from = 0, to = 0.5, by = 0.05)
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

error.rates(d, threshold = 0.20)

#'
#' From a visual perspective and error rates we can see that above a value of 0.30 for individual_effects parameters, we start to observe increase of false positive.
#' We will use values of exposure_effects ranging from 0 to 0.20 for simulations without exposure bias and values ranging 0.4 to 0.6 for simulations with exposure bias.
#'

#' ## 5.1. An example of censoring intercept and slope are equal to 0.1 in simulated data
N_id = 100
Hairy = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
test = test.function(att = Hairy,
                     N_id = N_id,
                     individual_predictors=Hairy, # individuals characteristics
                     individual_effects=matrix(c(0.2,0.2),ncol=1, nrow=2), # individuals characteristics on interaction probability
                     exposure_predictors = NULL,
                     exposure_effects = c(0, 0), exposure_sigma = 1, # exposure effect
                     int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf),
                     simulate.interactions = T) #no censoring effect
test$plots

#' ## 5.2. An example of individual_effects being equal to 0.4 in simulated data
test = test.function(att = Hairy,
                     N_id = N_id,
                     individual_predictors=Hairy, # individuals characteristics
                     individual_effects=matrix(c(0.18,0.18),ncol=1, nrow=2), # individuals characteristics on interaction probability
                     exposure_predictors = NULL,
                     exposure_effects = c(0, 0), exposure_sigma = 1, # exposure effect
                     int_intercept = c(Inf,Inf), int_slope = c(Inf,Inf),
                     simulate.interactions = T) #no censoring effect
test$plots


save.image(file='2.Results/Appendices/Appendix.RData')
