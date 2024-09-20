####################################################################################### Setup
# Simulator function
simulate_series = function(y, z, Q, S, R){                                                               
set.seed(666)
x = z+y
V = 1            # One blocking variable
G = 3            # Three categories in this variable
N_id = 65        # Number of people

Group = sample(1:3, N_id, replace=TRUE)
B = matrix(-11 + Q[1], nrow=G, ncol=G)
diag(B) = -7.2 + Q[2]# Block matrix

B[1,3] = -8.1 + Q[3]
B[3,2] = -8.9 + Q[4]

alpha = R[1]
Coloration = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
Shyness = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)*alpha + Coloration*(1-alpha)
SizeDiff = array(rnorm(N_id*N_id, 0, 1), c(N_id, N_id, 1))
                                                               
A = simulate_sbm_plus_srm_network_with_measurement_bias(N_id = N_id, 
                                                   B=list(B=B),
                                                   V=V,
                                                   groups=data.frame(Group=factor(Group)),
                                                   individual_predictors=Coloration,
                                                   individual_effects=matrix(c(1.75, 0.75),ncol=1, nrow=2),
                                                   dyadic_predictors = SizeDiff,
                                                   dyadic_effects = c(Q[5]),
                                                   sr_mu = c(0,0),
                                                   sr_sigma = c(S[1], S[2]),
                                                   sr_rho = R[2],
                                                   dr_mu = c(0,0),
                                                   dr_sigma = S[3],
                                                   dr_rho = R[3],
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

##################################### Now run the simulations. Sampling model
N_out = 200
set.seed(3494)

deploy_on_cores = function(z){
 set.seed(z)
 Q = rnorm(5, 0, 2.5) 
 R = rbeta(3, 2, 2)
 S = abs(rnorm(3, 0, 2))
 res_fit_0 = simulate_series(0.00, 0, Q, S, R)
 res_fit_1 = simulate_series(-0.75, 0, Q, S, R)
 res_fit_2 = simulate_series(-1.75, 0, Q, S, R)
 res_fit_3 = simulate_series(-2.75, 0, Q, S, R)
 res_all = rbind(res_fit_0, res_fit_1, res_fit_2, res_fit_3)
 return(res_all)
}

fit = mclapply(1:N_out, function(z){
                       deploy_on_cores(z)
                       }, mc.cores = 100)

res_basic_samp = do.call(rbind, fit)

##################################### Prepare figures
res_all_basic = prepare_robustness_plots(res_basic_samp)

cols = plvs_vltra("mystic_mausoleum")

##### Figure type 1
p1 = ggplot(res_all_basic, aes(y=Model, x=P, group=Model, color=Model, fill=Model)) + 
   geom_jitter(alpha=0.5, height=0.25, width=0.025, shape=20) + geom_boxplot(outliers = FALSE, col="black",alpha=0.7)+
  facet_grid( Variable ~ Setting ) + xlab("P value") + ylab("Model") +
  geom_vline(xintercept=0.05, color=cols[9], linetype="dashed") + 
  #geom_hline(yintercept=-0.05, color="darkred", linetype="dashed") +
  theme(legend.position="none") +scale_fill_manual(
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
    )+ theme(panel.spacing = unit(0.7, "lines"))+ coord_cartesian(xlim = c(0, 1))

 p1
ggsave("Sampling_Robustness_data.pdf", p1, width=11, height=8, device = cairo_pdf)


##### Figure type 2
res_all_basic_agg = aggregate_error_rates(res_all_basic)

p3 = ggplot(res_all_basic_agg, aes(x = Model, y = ResultWrong, color = Model)) +
  geom_segment(
    aes(x = Model, xend = Model, y = 0, yend = ResultWrong)
    ) + coord_flip()+ ylab("Frequency of incorrrect inference")+xlab("")+ theme_bw()+
  geom_point(aes(color = Model), size = 3, shape=18) +
  facet_grid( Variable ~  Setting) +
     theme(legend.position="none")  +
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
    ) +ylim(0,1)+ theme(panel.spacing = unit(0.7, "lines"))
p3
   ggsave("Sampling_Robustness_Freq.pdf", p3, width=11, height=8, device = cairo_pdf)

res_all_basic_agg$Outcome = "Sampling"
df_to_merg_sampling = res_all_basic_agg
save(df_to_merg_sampling, file = "df_to_merg_sampling.RData")



##################################### Now run the simulations. Censoring model
N_out = 200
set.seed(3494)

deploy_on_cores = function(z){
 set.seed(z)
 Q = rnorm(5, 0, 2.5) 
 R = rbeta(3, 2, 2)
 S = abs(rnorm(3, 0, 2))
 res_fit_0 = simulate_series(0, 0.00, Q, S, R)
 res_fit_1 = simulate_series(0, 0.75, Q, S, R)
 res_fit_2 = simulate_series(0, 1.75, Q, S, R)
 res_fit_3 = simulate_series(0, 2.75, Q, S, R)
 res_all = rbind(res_fit_0, res_fit_1, res_fit_2, res_fit_3)
 return(res_all)
}

fit = mclapply(1:N_out, function(z){
                       deploy_on_cores(z)
                       }, mc.cores = 100)

res_basic_cens = do.call(rbind, fit)

##################################### Prepare figures
res_all_basic = prepare_robustness_plots(res_basic_cens)

##### Figure type 1
p1 = ggplot(res_all_basic, aes(y=Model, x=P, group=Model, color=Model, fill=Model)) + 
   geom_jitter(alpha=0.5, height=0.25, width=0.025, shape=20) + geom_boxplot(outliers = FALSE, col="black",alpha=0.7)+
  facet_grid( Variable ~ Setting ) + xlab("P value") + ylab("Model") +
  geom_vline(xintercept=0.05, color=cols[9], linetype="dashed") + 
  #geom_hline(yintercept=-0.05, color="darkred", linetype="dashed") +
  theme(legend.position="none") + scale_fill_manual(
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
    )+ theme(panel.spacing = unit(0.7, "lines"))+ coord_cartesian(xlim = c(0, 1))

 p1
ggsave("Censoring_Robustness_data.pdf", p1, width=11, height=8, device = cairo_pdf)

##### Figure type 2
res_all_basic_agg = aggregate_error_rates(res_all_basic)

p3 = ggplot(res_all_basic_agg, aes(x = Model, y = ResultWrong, color = Model)) +
  geom_segment(
    aes(x = Model, xend = Model, y = 0, yend = ResultWrong)
    ) + coord_flip()+ ylab("Frequency of incorrrect inference")+xlab("")+ theme_bw()+
  geom_point(aes(color = Model), size = 3, shape=18) +
  facet_grid( Variable ~  Setting) +
     theme(legend.position="none")  +
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
    )+ylim(0,1)+ theme(panel.spacing = unit(0.7, "lines"))
p3
   ggsave("Censoring_Robustness_Freq.pdf", p3, width=11, height=8, device = cairo_pdf)

res_all_basic_agg$Outcome = "Censoring"
df_to_merg_censoring = res_all_basic_agg
save(df_to_merg_censoring, file = "df_to_merg_censoring.RData")



##################################################################### Merged plot
load("df_to_merg_censoring.RData")
load("df_to_merg_sampling.RData")

df_all_performance = rbind(df_to_merg_censoring, df_to_merg_sampling)

df_all_performance$Setting = as.character(df_all_performance$Setting)
df_all_performance$Setting[which(df_all_performance$Setting=="κ = -0.75")] = "|κ| = 0.75"
df_all_performance$Setting[which(df_all_performance$Setting=="κ = -1.75")] = "|κ| = 1.75"
df_all_performance$Setting[which(df_all_performance$Setting=="κ = -2.75")] = "|κ| = 2.75"

df_all_performance$Setting[which(df_all_performance$Setting=="κ = 0")] = "|κ| = 0.00"

df_all_performance$Setting[which(df_all_performance$Setting=="κ = 0.75")] = "|κ| = 0.75"
df_all_performance$Setting[which(df_all_performance$Setting=="κ = 1.75")] = "|κ| = 1.75"
df_all_performance$Setting[which(df_all_performance$Setting=="κ = 2.75")] = "|κ| = 2.75"
df_all_performance$Setting = factor(df_all_performance$Setting)


p4 = ggplot(df_all_performance, aes(x = Model, y = ResultWrong, color = Model,  group=Outcome)) +
       geom_linerange(aes(x = Model, ymin = 0, ymax = ResultWrong, colour = Model, linetype = Outcome), 
                   position = position_dodge(width = 0.5))+
    geom_point(aes(x = Model, y = ResultWrong, colour = Model),
               position = position_dodge(width = 0.5))+
  coord_flip()+ ylab("Frequency of incorrrect inference")+xlab("")+ theme_bw()+
  facet_grid( Variable ~  Setting) +
     theme(legend.position="bottom")  +
  scale_color_manual(
    values=c("STRAND Basic"= cols[5], 
             "STRAND ME"=cols[7],
             "ANTs Rates"=cols[12], 
             "ANTs Rates Permuted"=cols[11],
             "lm Rates"=cols[3], 
             "AMEN" = cols[10],
             "asnipe" = cols[9],
             "bison" = cols[1] 
             ),
    guide = "none"
    ) +ylim(0,1)+ geom_hline(yintercept=0.05, color=cols[9], linetype="dashed")+
  theme(panel.spacing = unit(0.7, "lines")) + scale_linetype_manual(values=c("solid","twodash"))
p4
   ggsave("Robustness_Freq.pdf", p4, width=11, height=8, device = cairo_pdf)
