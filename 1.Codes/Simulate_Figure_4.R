####################################################################################### Setup
# Simulator function
simulate_series = function(x,y,z){                                                               
set.seed(666)

V = 1            # One blocking variable
G = 3            # Three categories in this variable
N_id = 65        # Number of people

Group = sample(1:3, N_id, replace=TRUE)
B = matrix(-11, nrow=G, ncol=G)
diag(B) = -7.2 # Block matrix

B[1,3] = -8.1
B[3,2] = -8.9

alpha = 0.75
Coloration = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
Shyness = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)*alpha + Coloration*(1-alpha)
SizeDiff = array(rnorm(N_id*N_id, 0, 1), c(N_id, N_id, 1))
                                                               
A = simulate_sbm_plus_srm_network_with_measurement_bias(N_id = N_id, 
                                                   B=list(B=B),
                                                   V=V,
                                                   groups=data.frame(Group=factor(Group)),
                                                   individual_predictors=Coloration,
                                                   individual_effects=matrix(c(x, 0.55),ncol=1, nrow=2),
                                                   dyadic_predictors = SizeDiff,
                                                   dyadic_effects = c(-1.5),
                                                   sr_mu = c(0,0),
                                                   sr_sigma = c(1.4, 1.8),
                                                   sr_rho = 0.5,
                                                   dr_mu = c(0,0),
                                                   dr_sigma = 1.2,
                                                   dr_rho = 0.75,
                                                   exposure_mu = 5.0,
                                                   exposure_sigma = 0.25,
                                                   exposure_max = 20,
                                                   censoring_mu = -3.75,
                                                   censoring_sigma = 0.25,
                                                   exposure_predictors = Coloration,
                                                   censoring_predictors = Coloration,
                                                   exposure_effects = c(y),
                                                   censoring_effects = c(z)
                                                   )
# Prep data
grooming = list(Grooming = A$net)
exposure = list(Exposure = A$true_samps)
dyad = list(SizeDiff = SizeDiff)
block = data.frame(Group = as.factor(Group))
indiv =  data.frame(Coloration = Coloration, Shyness=Shyness)

model_dat = make_strand_data(outcome = grooming,
                             individual_covariates = indiv, 
                             block_covariates = block,
                             dyadic_covariates = dyad,
                             outcome_mode = "binomial",
                             exposure = exposure
                             )

model_dat$sampled = A$true_exposure
model_dat$sampled_exposure = rep(20, N_id)

model_dat$detected = A$detected
model_dat$detected_exposure = A$trials 

#################################### All methods functions go here
res1 = fit_lm_Rate_Unweighted(model_dat,x)
res3 = fit_ANTs_Rate_Unweighted(model_dat,x)
res4 = fit_ANTs_Rate_Perm_Unweighted(model_dat,x)
res5 = fit_asnipe(model_dat, x)
res6 = fit_AMEN(model_dat, x)
res7 = fit_bison(model_dat, x)
res8 = fit_STRAND_basic(model_dat, x)
res9 = fit_STRAND_ME(model_dat, x)

res = rbind(res1, res3, res4, res5, res6, res7, res8, res9)
print(x)
return(res)
}




#################################### Now run the simulations. No bias model
set.seed(3494)
 N_out = 13
 samps = 5

deploy_on_cores = function(z){
 slopes = rep(seq(-2.0, 2.0, length.out=N_out),samps)
 res = simulate_series(slopes[z], 0, 0) # no bias
 return(res)
}

fit = mclapply(1:(N_out*samps), function(z){
                       deploy_on_cores(z)
                       }, mc.cores = N_out*samps, mc.preschedule=FALSE)



# ##################################### OR run the simulations without supercomputer. Basic model
# N_out = 17
# results_basic = vector("list", N_out)
# slopes = seq(-2.0, 2.0, length.out=N_out)

# for(i in 1:N_out){
#  results_basic[[i]] = simulate_series(slopes[i], 0, 0) # no biases
#  }

#  results_basic_nb = results_basic

################################
results_basic_nb = fit
res_all_basic_nb = prepare_series_plots(results_basic_nb)
res_all_basic_nb = aggregate_p_values(res_all_basic_nb)

cols = plvs_vltra("mystic_mausoleum")

p3 = ggplot(res_all_basic_nb, aes(x=Setting, y=P, group=Model, color=Model, fill=Model)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid( . ~ Variable ) + xlab("True effect of Coloration on out-degree, \u03ba") + ylab("P") +
  geom_hline(yintercept=0.05, color="darkred", linetype="dashed") + 
  #geom_hline(yintercept=-0.05, color="darkred", linetype="dashed") +
  theme(legend.position="bottom") +  scale_fill_manual(
    values=c("STRAND Basic"= cols[5], 
             "STRAND ME"=cols[7],
             "ANTs Rates"=cols[12], 
             "ANTs Rates Permuted"=cols[11],
             "lm Rates"=cols[3], 
             "AMEN" = cols[10],
             "asnipe" = cols[9],
             "bison" = cols[1] 
             ) 
    ) +
  scale_color_manual(
    values=c("STRAND Basic"= cols[5], 
             "STRAND ME"=cols[7],
             "ANTs Rates"=cols[12], 
             "ANTs Rates Permuted"=cols[11],
             "lm Rates"=cols[3], 
             "AMEN" = cols[10],
             "asnipe" = cols[9],
             "bison" = cols[1] 
             )
    )
p3
 ggsave("Basic_model.pdf", p3, width=11, height=4, device = cairo_pdf)


##################################### Now run the simulations. Sampling model
set.seed(3494)
 N_out = 13
 samps = 5

deploy_on_cores = function(z){
 slopes = rep(seq(-2.0, 2.0, length.out=N_out),samps)
 res = simulate_series(slopes[z], -2.75, 0) # sampling bias
 return(res)
}

fit = mclapply(1:(N_out*samps), function(z){
                       deploy_on_cores(z)
                       }, mc.cores = N_out*samps, mc.preschedule=FALSE)


library(doSNOW)
library(foreach)
library(doParallel)
cl <- makeCluster(4, outfile="")
registerDoSNOW(cl)
registerDoParallel(cores=cl)
slopes = rep(seq(-2, 2, length.out=N_out),samps)
res = foreach(i = 1:(N_out*samps), .export = ls(globalenv())) %dopar% {
  library(igraph)
  library(ggplot2)
  library(cmdstanr)
  library(ANTs)
  library(mice)
  library(parallel)
  
  library(bisonR)
  library(STRAND)
  library(amen)
  library(asnipe)
  
  source("Support_Functions.R")
  res = simulate_series(slopes[i], -2.75, 0) # sampling bias
  res
}
##################################### OR run the simulations without supercomputer. Sampling model
# N_out = 13
# results_basic = vector("list", N_out)
# slopes = seq(-2.0, 2.0, length.out=N_out)

# for(i in 1:N_out){
#  results_basic[[i]] = simulate_series(slopes[i], -2.75, 0) # sampling bias
#  }

# results_basic_samp = results_basic

####################################
results_basic_samp = fit
res_all_basic_samp = prepare_series_plots(results_basic_samp)
res_all_basic_samp = aggregate_p_values(res_all_basic_samp)


p2 = ggplot(res_all_basic_samp, aes(x=Setting, y=P, group=Model, color=Model, fill=Model)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid( . ~ Variable ) + xlab("True effect of Coloration on out-degree, \u03ba") + ylab("P") +
  geom_hline(yintercept=0.05, color="darkred", linetype="dashed") + 
  #geom_hline(yintercept=-0.05, color="darkred", linetype="dashed") +
  theme(legend.position="bottom") +  scale_fill_manual(
    values=c("STRAND Basic"= cols[5], 
             "STRAND ME"=cols[7],
             "ANTs Rates"=cols[12], 
             "ANTs Rates Permuted"=cols[11],
             "lm Rates"=cols[3], 
             "AMEN" = cols[10],
             "asnipe" = cols[9],
             "bison" = cols[1] 
             ) 
    ) +
  scale_color_manual(
    values=c("STRAND Basic"= cols[5], 
             "STRAND ME"=cols[7],
             "ANTs Rates"=cols[12], 
             "ANTs Rates Permuted"=cols[11],
             "lm Rates"=cols[3], 
             "AMEN" = cols[10],
             "asnipe" = cols[9],
             "bison" = cols[1] 
             )
    )

 ggsave("Sampling_model.pdf", p2, width=11, height=4, device = cairo_pdf)





##################################### Now run the simulations. Censoring model
 set.seed(3494)
 N_out = 13
 samps = 5

deploy_on_cores = function(z){
 slopes = rep(seq(-2.0, 2.0, length.out=N_out),samps)
 res = simulate_series(slopes[z], 0, 2.75) # sampling bias
 return(res)
}

fit = mclapply(1:(N_out*samps), function(z){
                       deploy_on_cores(z)
                       }, mc.cores = N_out*samps, mc.preschedule=FALSE)

##################################### OR run the simulations without supercomputer. Censoring model
# N_out = 13
# results_basic = vector("list", N_out)
# slopes = seq(-2.0, 2.0, length.out=N_out)

# for(i in 1:N_out){
#  results_basic[[i]] = simulate_series(slopes[i], 0, 2.75) # censoring bias
#  }

#  res_basic_cens = results_basic

#########################
results_basic_cens = fit
res_all_basic_cens = prepare_series_plots(results_basic_cens)
res_all_basic_cens = aggregate_p_values(res_all_basic_cens)


p1 = ggplot(res_all_basic_cens, aes(x=Setting, y=P, group=Model, color=Model, fill=Model)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid( . ~ Variable ) + xlab("True effect of Coloration on out-degree, \u03ba") + ylab("P") +
  geom_hline(yintercept=0.05, color="darkred", linetype="dashed") + 
  #geom_hline(yintercept=-0.05, color="darkred", linetype="dashed") +
  theme(legend.position="bottom") +  scale_fill_manual(
    values=c("STRAND Basic"= cols[5], 
             "STRAND ME"=cols[7],
             "ANTs Rates"=cols[12], 
             "ANTs Rates Permuted"=cols[11],
             "lm Rates"=cols[3], 
             "AMEN" = cols[10],
             "asnipe" = cols[9],
             "bison" = cols[1] 
             ) 
    ) +
  scale_color_manual(
    values=c("STRAND Basic"= cols[5], 
             "STRAND ME"=cols[7],
             "ANTs Rates"=cols[12], 
             "ANTs Rates Permuted"=cols[11],
             "lm Rates"=cols[3], 
             "AMEN" = cols[10],
             "asnipe" = cols[9],
             "bison" = cols[1] 
             )
    )

 ggsave("Censoring_model.pdf", p1, width=11, height=4, device = cairo_pdf)

