####################################################################################### Censoring model                                   
simulate_series_cens = function(x){                                                               
set.seed(420)

V = 1            # One blocking variable
G = 3            # Three categories in this variable
N_id = 65        # Number of people

Group = sample(1:3, N_id, replace=TRUE)
B = matrix(-12, nrow=G, ncol=G)
diag(B) = -8.2 # Block matrix

B[1,3] = -9.1
B[3,2] = -9.9

Coloration = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
SizeDiff = array(rnorm(N_id*N_id, 0, 1), c(N_id, N_id, 1))
                                                               
A = simulate_sbm_plus_srm_network_with_measurement_bias(N_id = N_id, 
                                                   B=list(B=B),
                                                   V=V,
                                                   groups=data.frame(Group=factor(Group)),
                                                   individual_predictors=Coloration,
                                                   individual_effects=matrix(c(1.9, 0.55),ncol=1, nrow=2),
                                                   dyadic_predictors = SizeDiff,
                                                   dyadic_effects = c(-1.5),
                                                   sr_mu = c(0,0),
                                                   sr_sigma = c(1.4, 0.8),
                                                   sr_rho = 0.6,
                                                   dr_mu = c(0,0),
                                                   dr_sigma = 1.0,
                                                   dr_rho = 0.75,
                                                   exposure_mu = 5.0,
                                                   exposure_sigma = 0.001,
                                                   exposure_max = 20,
                                                   censoring_mu = -3.5,
                                                   censoring_sigma = 0.5,
                                                   exposure_predictors = Coloration,
                                                   censoring_predictors = Coloration,
                                                   exposure_effects = c(0.0),
                                                   censoring_effects = c(x)
                                                   )
# Prep data
grooming = list(Grooming = A$net)
exposure = list(Exposure = A$true_samps)
dyad = list(SizeDiff = SizeDiff)
block = data.frame(Group = as.factor(Group))
indiv =  data.frame(Coloration = Coloration)

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

# Basic model
fit =  fit_block_plus_social_relations_model(data=model_dat,
                              block_regression = ~ Group,
                              focal_regression = ~ Coloration,
                              target_regression = ~ Coloration,
                              dyad_regression = ~  SizeDiff,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 600, iter_sampling = 600,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

res1 = summarize_strand_results(fit)
res1$summary$Model = "Basic"

# Censoring model with correct specification                                                                
fit =  fit_block_plus_social_relations_model_with_measurement_bias(data=model_dat,
                              block_regression = ~ Group,
                              focal_regression = ~ Coloration,
                              target_regression = ~ Coloration,
                              sampling_regression = ~ Coloration,
                              censoring_regression = ~ Coloration,
                              dyad_regression = ~  SizeDiff,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 600, iter_sampling = 600,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

res2 = summarize_bsrm_results_with_measurement_bias(fit)
res2$summary$Model = "Basic+ME"

res = rbind(res1$summary,res2$summary)
print(x)
return(res)
}

####################################################################################### Run simulations and model fit  
Nout = 11
results = vector("list", Nout)
slopes = seq(-4.5, 4.5, length.out=Nout)

for(i in 1:Nout){
 results[[i]] = simulate_series_cens(slopes[i])
 }

results_backup = results


for(i in 1:Nout){
 results[[i]]$Setting = slopes[i]
 }


res_all_censoring = do.call(rbind, results)

###### Plot results
colnames(res_all_censoring) = c("Variable", "Median", "L", "H", "Mean", "SD", "Model", "Setting")
res_all_censoring$Median = as.numeric(res_all_censoring$Median)
res_all_censoring$L = as.numeric(res_all_censoring$L)
res_all_censoring$H = as.numeric(res_all_censoring$H)
in_set = c("focal effects coeffs (out-degree), Coloration","target effects coeffs (in-degree), Coloration","censoring effects coeffs, Coloration")

res_all_censoring2 = res_all_censoring[which(res_all_censoring$Variable %in% in_set),]

dummy1 = data.frame(Variable = in_set, Z = c(1.9, 0.5, 0), Slope=c(0,0,1))

p1 = ggplot(res_all_censoring2, aes(x=Setting, y=Median, group=Model, color=Model, fill=Model)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=L, ymax=H), color=NA, alpha=0.5) +
  facet_wrap(~ Variable, nrow=3, scale="free") + xlab("Censoring bias") + ylab("Posterior effect size")+
  geom_abline(data = dummy1, aes(intercept = Z, slope = Slope), linetype="dashed") + 
  theme(legend.position="bottom") + scale_fill_manual(values=c("Basic"= "#be9a3f", "Basic+ME"="#447c81")) +
  scale_color_manual(values=c("Basic"= "#be9a3f", "Basic+ME"="#447c81"))

# ggsave("Censoring_model.pdf", p1, height=11, width=4)



####################################################################################### Sampling bias model                                   
simulate_series_exposure = function(x){                                                               
set.seed(420)

V = 1            # One blocking variable
G = 3            # Three categories in this variable
N_id = 65        # Number of people

Group = sample(1:3, N_id, replace=TRUE)
B = matrix(-12, nrow=G, ncol=G)
diag(B) = -8.2 # Block matrix

B[1,3] = -9.1
B[3,2] = -9.9

Coloration = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
SizeDiff = array(rnorm(N_id*N_id, 0, 1), c(N_id, N_id, 1))
                                                               
A = simulate_sbm_plus_srm_network_with_measurement_bias(N_id = N_id, 
                                                   B=list(B=B),
                                                   V=V,
                                                   groups=data.frame(Group=factor(Group)),
                                                   individual_predictors=Coloration,
                                                   individual_effects=matrix(c(1.9, 0.55),ncol=1, nrow=2),
                                                   dyadic_predictors = SizeDiff,
                                                   dyadic_effects = c(-1.5),
                                                   sr_mu = c(0,0),
                                                   sr_sigma = c(1.4, 0.8),
                                                   sr_rho = 0.6,
                                                   dr_mu = c(0,0),
                                                   dr_sigma = 1.0,
                                                   dr_rho = 0.75,
                                                   exposure_mu = 3.5,
                                                   exposure_sigma = 0.5,
                                                   exposure_max = 20,
                                                   censoring_mu = -5,
                                                   censoring_sigma = 0.001,
                                                   exposure_predictors = Coloration,
                                                   censoring_predictors = Coloration,
                                                   exposure_effects = c(x),
                                                   censoring_effects = c(0)
                                                   )
# Prep data
grooming = list(Grooming = A$net)
exposure = list(Exposure = A$true_samps)
dyad = list(SizeDiff = SizeDiff)
block = data.frame(Group = as.factor(Group))
indiv =  data.frame(Coloration = Coloration)

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

# Basic model
fit =  fit_block_plus_social_relations_model(data=model_dat,
                              block_regression = ~ Group,
                              focal_regression = ~ Coloration,
                              target_regression = ~ Coloration,
                              dyad_regression = ~  SizeDiff,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 600, iter_sampling = 600,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

res1 = summarize_strand_results(fit)
res1$summary$Model = "Basic"

# Censoring model with correct specification                                                                
fit =  fit_block_plus_social_relations_model_with_measurement_bias(data=model_dat,
                              block_regression = ~ Group,
                              focal_regression = ~ Coloration,
                              target_regression = ~ Coloration,
                              sampling_regression = ~ Coloration,
                              censoring_regression = ~ Coloration,
                              dyad_regression = ~  SizeDiff,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 600, iter_sampling = 600,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

res2 = summarize_bsrm_results_with_measurement_bias(fit)
res2$summary$Model = "Basic+ME"

res = rbind(res1$summary,res2$summary)
print(x)
return(res)
}

####################################################################################### Run simulations and model fit 
Nout = 11
results = vector("list", Nout)
slopes = seq(-4.5, 4.5, length.out=Nout)

for(i in 1:Nout){
 results[[i]] = simulate_series_exposure(slopes[i])
 }

results_backup = results


for(i in 1:Nout){
 results[[i]]$Setting = slopes[i]
 }


res_all_sampling = do.call(rbind, results)


########### Plot results
colnames(res_all_sampling) = c("Variable", "Median", "L", "H", "Mean", "SD", "Model", "Setting")
res_all_sampling$Median = as.numeric(res_all_sampling$Median)
res_all_sampling$L = as.numeric(res_all_sampling$L)
res_all_sampling$H = as.numeric(res_all_sampling$H)
in_set = c("sampling effects coeffs, Coloration","focal effects coeffs (out-degree), Coloration","target effects coeffs (in-degree), Coloration")
res_all_sampling2 = res_all_sampling[which(res_all_sampling$Variable %in% in_set),]
res_all_sampling2$Variable = factor(res_all_sampling2$Variable)
res_all_sampling2$Variable = factor(res_all_sampling2$Variable,levels=c("sampling effects coeffs, Coloration","focal effects coeffs (out-degree), Coloration","target effects coeffs (in-degree), Coloration"))

dummy2 = data.frame(Variable = in_set, Z = c(0, 1.9, 0.5), Slope=c(1,0,0))
dummy2$Variable = factor(dummy2$Variable)
dummy2$Variable = factor(dummy2$Variable,levels=c("sampling effects coeffs, Coloration","focal effects coeffs (out-degree), Coloration","target effects coeffs (in-degree), Coloration"))

p2 = ggplot(res_all_sampling2, aes(x=Setting, y=Median, group=Model, color=Model, fill=Model)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=L, ymax=H), color=NA, alpha=0.5) +
  facet_wrap(~ Variable, nrow=3, scale="free") + xlab("Sampling bias") + ylab("Posterior effect size")+
  geom_abline(data = dummy2, aes(intercept = Z, slope = Slope), linetype="dashed") + 
  theme(legend.position="bottom") + scale_fill_manual(values=c("Basic"= "#be9a3f", "Basic+ME"="#447c81")) +
  scale_color_manual(values=c("Basic"= "#be9a3f", "Basic+ME"="#447c81")) + coord_cartesian(xlim =c(-4.5, 4.5), ylim = c(-4.5, 4.5))

# ggsave("sampling_model.pdf", p2, height=11, width=4)



####################################################################### Merged plot
res_all_sampling2$Pred="Sampling"
res_all_censoring2$Pred="Censoring"
res_all = do.call(rbind, list(res_all_censoring2,res_all_sampling2))

res_all$Variable = as.character(res_all$Variable)
res_all$Variable[which(res_all$Variable=="sampling effects coeffs, Coloration")] = "Effect of Coloration on sampling/censoring"
res_all$Variable[which(res_all$Variable=="censoring effects coeffs, Coloration")] = "Effect of Coloration on sampling/censoring"
res_all$Variable[which(res_all$Variable=="focal effects coeffs (out-degree), Coloration")] = "Effect of Coloration on out-degree"
res_all$Variable[which(res_all$Variable=="target effects coeffs (in-degree), Coloration")] = "Effect of Coloration on in-degree"
res_all$Variable = as.factor(res_all$Variable)

dummy2$Pred="Sampling"
dummy1$Pred="Censoring"
dummy0=rbind(dummy1, dummy2)

dummy0$Pred = factor(dummy0$Pred)
dummy0$Pred = factor(dummy0$Pred,levels=c("Sampling", "Censoring"))

res_all$Pred = factor(res_all$Pred)
res_all$Pred = factor(res_all$Pred,levels=c("Sampling", "Censoring"))

dummy0$Variable = c("Effect of Coloration on out-degree", "Effect of Coloration on in-degree", 
                    "Effect of Coloration on sampling/censoring","Effect of Coloration on sampling/censoring", 
                    "Effect of Coloration on out-degree", "Effect of Coloration on in-degree" ) 

res_all$Variable = factor(res_all$Variable, levels=c("Effect of Coloration on sampling/censoring", 
                    "Effect of Coloration on out-degree", "Effect of Coloration on in-degree"))
dummy0$Variable = factor(dummy0$Variable, levels=c("Effect of Coloration on sampling/censoring", 
                    "Effect of Coloration on out-degree", "Effect of Coloration on in-degree"))

cols = plvs_vltra("mystic_mausoleum")

res_all$Model[which(res_all$Model=="Basic")] = "STRAND Basic"
res_all$Model[which(res_all$Model=="Basic+ME")] = "STRAND ME"


p3 = ggplot(res_all, aes(x=Setting, y=Median, group=Model, color=Model, fill=Model)) + 
geom_ribbon(aes(ymin=L, ymax=H), color=NA, alpha=0.45) +
  geom_line() + geom_point() + theme_bw() + 
   geom_hline(yintercept=0) +
  facet_grid(Variable~Pred , scale="free") + xlab("Effect of Coloration on sampling bias/censoring") + ylab("Posterior effect size")+
  geom_abline(data = dummy0, aes(intercept = Z, slope = Slope), linetype="dashed", color="black") + 
  theme(legend.position="bottom") + scale_fill_manual(values=c("STRAND Basic"= cols[5], "STRAND ME"=cols[7])) +
  scale_color_manual(values=c("STRAND Basic"= cols[5], "STRAND ME"=cols[7])) 

ggsave("Parameter_Recovery.pdf", p3, height=11, width=6)























