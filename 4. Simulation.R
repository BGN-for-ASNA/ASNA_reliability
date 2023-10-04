###############################
#######  Iterate over grid ####
###############################

# Varying parameters of the simulations.
N_id = 30:100
hairy_tie_effect = seq(-3.5, 3.5, by = 0.5)
hairy_detect_effect = seq(-3.5, 3.5, by = 0.5)
grid = expand.grid(N_id, hairy_tie_effect, hairy_detect_effect)
grid_subsample = grid[sample(1:nrow(grid), 100, replace = FALSE),]
colnames(grid_subsample) = c("N_id", "hairy_tie_effect", "hairy_detect_effect")

# SUbgroups variables ----------
B = NULL
V = 1
groups=NULL


# Objects to store the results.
res_array_lmw = res_array_lm = res_array = array(NA,c(Reps, Reps,3))
hte_array = hde_array = array(NA,c(Reps, Reps))

get_res = function(x) {
  y = c(coef(x)[2],confint(x, level = 0.9)[2,])
  return(y)
}

#######################
#######  Make data ####
#######################
for(i in 1:nrow(grid_subsample)){
  Clique = rep(1, grid_subsample$N_id[i])
  Hairy = matrix(rnorm(grid_subsample$N_id[i], 0, 1), nrow=grid_subsample$N_id[i], ncol=1)
  A = simulate_sbm_plus_srm_network_with_measurement_bias(
    N_id = grid_subsample$N_id[i],
    B=list(B=B),
    V=V,
    groups=groups,
    individual_predictors=Hairy,
    individual_effects=matrix(c(grid_subsample$hairy_tie_effect[i], grid_subsample$hairy_tie_effect[i]),ncol=1, nrow=2),
    sr_sigma = c(0.1, 0.1),
    sr_rho = 0.0,
    dr_sigma = 1.2,
    dr_rho = 0.0,
    exposure_predictors = cbind(rep(1,grid_subsample$N_id[i]),Hairy),
    exposure_effects = c(-1, grid_subsample$hairy_detect_effect[i]),
    exposure_sigma = 2.9,
    exposure_baseline = 40
  )

  ###############################
  #######  STRAND ###############
  ###############################
  nets = list(Grooming = A$network)
  exposure_nets = list(Exposure = A$true_samps)
  block = data.frame(Clique=factor(clique))
  indiv =  data.frame(Hairy = Hairy)
  model_dat = make_strand_data(outcome = nets,
                               individual_covariates = indiv,
                               block_covariates = block,
                               outcome_mode = "binomial",
                               exposure = exposure_nets
  )

  # !!  I replaced fit_block_plus_social_relations_model by fit_social_relations_model
  fit =  fit_social_relations_model(data=model_dat,
                                    focal_regression = ~ Hairy,
                                    target_regression = ~ Hairy,
                                    dyad_regression = ~  1,
                                    mode="mcmc",
                                    stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                iter_warmup = 1000, iter_sampling = 1000,
                                                                max_treedepth = NULL, adapt_delta = .98)
  )
  res = summarize_strand_results(fit)
  res_array[i,j,] = c(unlist(res$summary[2,2:4]))

  ###############################
  #######  ANTs unweighted ######
  ###############################
  tie_strength = A$network/(1+A$true_samps)
  colnames(tie_strength) = rownames(tie_strength) = 1:50
  df = ANTs:::df.create(tie_strength)
  df = met.strength(tie_strength, df = df, dfid = 1)
  df = met.degree(tie_strength, df = df, dfid = 1)
  df = met.eigen(tie_strength, df = df, dfid = 1)
  df$hair = indiv$Hairy
  test1.1 = lm(degree ~ hair, data = df)
  res_array_lm[i,j,] = c(unlist(get_res(test1.1)))

  ###############################
  #######  ANTs weighted ######
  ###############################
  tie_strength = A$network/(1+A$true_samps)
  colnames(tie_strength) = rownames(tie_strength) = 1:50
  df = ANTs:::df.create(tie_strength)
  df = met.strength(tie_strength, df = df, dfid = 1)
  df = met.degree(tie_strength, df = df, dfid = 1)
  df = met.eigen(tie_strength, df = df, dfid = 1)
  df$hair = indiv$Hairy
  test1.1 = lm(degree ~ hair, data = df, weights = A$true_exposure)
  res_array_lmw[i,j,] = c(unlist(get_res(test1.1)))
}


for(i in 1:Reps){
  for(j in 1:Reps){
    hte_array[i,j] = hairy_tie_effect[i]
    hde_array[i,j] = hairy_detect_effect[j]
  }}

df_all_true = data.frame(tie_effect = c(hte_array), detect_effect = c(hde_array), M = c(hte_array), L = NA, H = NA, model="True")
df_all_strand = data.frame(tie_effect = c(hte_array), detect_effect = c(hde_array), M = c(res_array[,,1]), L = c(res_array[,,2]), H = c(res_array[,,3]), model="STRAND")
df_all_ants = data.frame(tie_effect = c(hte_array), detect_effect = c(hde_array), M = c(res_array_lm[,,1]), L = c(res_array_lm[,,2]), H = c(res_array_lm[,,3]), model="ANTs")
df_all_antsw = data.frame(tie_effect = c(hte_array), detect_effect = c(hde_array), M = c(res_array_lmw[,,1]), L = c(res_array_lmw[,,2]), H = c(res_array_lmw[,,3]), model="ANTsW")

df_all = rbind(df_all_true, df_all_strand, df_all_ants, df_all_antsw)
df_all$M = as.numeric(df_all$M)
df_all$L = as.numeric(df_all$L)
df_all$H = as.numeric(df_all$H)
df_all$detect_effect = round(df_all$detect_effect,1)
df_all$tie_effect = round(df_all$tie_effect,0)

dodgewidth = 0.5
p1 = ggplot(df_all, aes(y=factor(tie_effect),  x=M, xmin=L, xmax=H , color=model)) +
  geom_linerange(position = position_dodge(width = dodgewidth)) +  geom_vline(xintercept = 0,linetype="dashed") +
  geom_point(position = position_dodge(width = dodgewidth)) +  facet_grid( . ~ detect_effect  , space="free") + theme(axis.text = element_text(size = 14))  +
  theme(axis.title = element_text(size = 14))  + ylab("True efect size") + xlab("Estimated effect size") + theme(strip.text.x = element_text(size = 14)) + coord_flip() +
  scale_colour_manual(values = c("True" = "black", "STRAND" = "#008080", "ANTs" = "darkred", "ANTsW" = "goldenrod4"))
p1

ggsave("ModelCompareANS.pdf",p1,width=12,height=3)

