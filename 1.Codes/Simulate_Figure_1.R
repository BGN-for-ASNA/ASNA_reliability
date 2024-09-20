####################################################### First Simulate Data with Bias
 set.seed(420)

 V = 1            # One blocking variable
 G = 3            # Three categories in this variable
 N_id = 65        # Number of people

 Group = sample(1:3, N_id, replace=TRUE)
 B = matrix(-12, nrow=G, ncol=G)
 diag(B) = -8.2 # Block matrix

 B[1,3] = -9.1
 B[3,2] = -9.9

 Coloration = matrix(rbinom(N_id, size=3, prob=0.65), nrow=N_id, ncol=1)
 SizeDiff = array(rnorm(N_id*N_id, 0, 1), c(N_id, N_id, 1))
                                                               
 A = simulate_sbm_plus_srm_network_with_measurement_bias(N_id = N_id, 
                                                   B=list(B=B),
                                                   V=V,
                                                   groups=data.frame(Group=factor(Group)),
                                                   individual_predictors=Coloration,
                                                   individual_effects=matrix(c(3.7, 5.3),ncol=1, nrow=2),
                                                   dyadic_predictors = SizeDiff,
                                                   dyadic_effects = c(-0.1),
                                                   sr_mu = c(0,0),
                                                   sr_sigma = c(0.2, 0.2), 
                                                   sr_rho = 0.6,
                                                   dr_mu = c(0,0),
                                                   dr_sigma = 1.2,
                                                   dr_rho = 0.5,
                                                   exposure_mu = 3.5,
                                                   exposure_sigma = 0.5,
                                                   exposure_max = 10,
                                                   censoring_mu = -5,
                                                   censoring_sigma = 0.001,
                                                   exposure_predictors = Coloration,
                                                   censoring_predictors = Coloration,
                                                   exposure_effects = c(-4.5),
                                                   censoring_effects = c(0)
                                                   )

########################################################### Plot Figure 1
 diag(A$network)=NA
 diag(A$tie_strength)=NA
 diag(A$true_samps)=NA
 A$tie_strength_ii = A$network/A$true_samps

 cols = plvs_vltra("mystic_mausoleum")

 A2 = data.frame(tie_strength=c(A$tie_strength), tie_strength_ii=c(A$tie_strength_ii), N=c(A$true_samps), Obs=c(A$network))
 A2 = A2[complete.cases(A2),]

 p1 = ggplot(A2, aes(x=tie_strength, y=tie_strength_ii)) + theme_bw() + 
     geom_point(size=1) + 
     ylab("Inferred tie strength") + 
     xlab("True tie strength") + 
     geom_abline(intercept = 0, slope = 1, color=cols[9], linetype="dashed")

 p2 = ggplot(A2, aes(x=tie_strength, y=tie_strength_ii, alpha=N)) + theme_bw() + 
     geom_point(size=1) + scale_alpha_continuous() + 
     ylab("Inferred tie strength") + 
     xlab("True tie strength") + 
     geom_abline(intercept = 0, slope = 1, color=cols[9], linetype="dashed") + 
     theme(legend.position = "none")

 ggsave("ScatterFrameA.png",p1,width=3, height=3)
 ggsave("ScatterFrameB.png",p2,width=3, height=3)


