################################### Need development version of STRAND to simulate the network data
 library(devtools)
 library(igraph)
 install_github('ctross/STRAND@measurement_error')

 library(STRAND)
 library(bisonR)

 ################################# Helper functions
# Bayesian P value (prob of sign error)
 P_se = function(x){
 M_x = mean(x)

 if(M_x<0){
   N_x = length(x)
   P_x = length(which(x>0))
  } else{
   N_x = length(x)
   P_x = length(which(x<0))
  }
 return(P_x/N_x)
    }

# Summarize samples from MCMC
process_samples = function(z, y1, y2, x){
 res = c(y1, median(z), HPDI(z, 0.95), mean(z), sd(z), y2, P_se(z), x)
 return(res)
}

################################## Create data
set.seed(123)

V = 1            # One blocking variable
G = 3            # Three categories in this variable
N_id = 65        # Number of people

Group = sample(1:3, N_id, replace=TRUE)
B = matrix(-13, nrow=G, ncol=G)
diag(B) = -8.2 # Block matrix

B[1,3] = -9.1
B[3,2] = -9.9

alpha = 0.75
Coloration = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
Shyness = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)*alpha + Coloration*(1-alpha)
SizeDiff = array(rnorm(N_id*N_id, 0, 1), c(N_id, N_id, 1))
                                                               
A = simulate_sbm_plus_srm_network_with_measurement_bias(N_id = N_id, 
                                  B=list(B=B),
                                  V=V,
                                  groups=data.frame(Group=factor(Group)),
                                  individual_predictors=Coloration,   # Coloration affects prob of individual "i" emiting ties (1.9) and "j" receiving ties (0.55). Shyness has no effect.
                                  individual_effects=matrix(c(1.9, 0.55),ncol=1, nrow=2),
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


########################################################### Quick Viz
Net = graph_from_adjacency_matrix(A$net, mode = c("directed"))
V(Net)$color = c("turquoise4","gray13", "goldenrod3")[Group]

plot(Net, edge.arrow.size =0.1, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)


############################## Now BisonR
# All relevant data is contained in model_dat, and is sent to the bisonR module
# There are 2 key variables of interest in the regression: Coloration and Shyness
# Coloration affects prob of individual "i" emiting ties (beta = 1.9) and "j" receiving ties (beta = 0.55). Shyness has no effects (beta=0 for both effects).

fit_bison = function(model_dat, x){
 # STRAND and other methods use matrices, rather than edge-lists, so we need to organize the data for bisonR

 # Create matrix form of N-vector data objects
 node_i_id = node_j_id = Shyness_j = Coloration_j = Shyness_i = Coloration_i = matrix(NA, nrow=model_dat$N_id, ncol=model_dat$N_id)
 
 # Loop over individuals, and fill matrices
 for(i in 1:model_dat$N_id){
  node_i_id[,i] = 1:model_dat$N_id
  node_j_id[i,] = 1:model_dat$N_id

  Coloration_i[,i] = model_dat$individual_predictors$Coloration
  Coloration_j[i,] = model_dat$individual_predictors$Coloration

  Shyness_i[,i] = model_dat$individual_predictors$Shyness
  Shyness_j[i,] = model_dat$individual_predictors$Shyness
  }

 # Turn all of the matrix data into long-form
 df_bison = data.frame(event=c(model_dat$outcomes), duration=c(model_dat$exposure), size_diff=c(model_dat$dyadic_predictors$SizeDiff),
                      node_i_id = c(node_i_id), node_j_id=c(node_j_id), coloration_i=c(Coloration_i), coloration_j=c(Coloration_j),
                      shyness_i=c(Shyness_i), shyness_j=c(Shyness_j)
                        )

 # Drop the rows that refer to ego-ego relationships
 df_bison = df_bison[which(df_bison$node_i_id != df_bison$node_j_id),]


#################################################### This comes from the tutorial
 # Fit Bison Model
  fit_edge = bison_model(
   (event | duration) ~ dyad(node_i_id, node_j_id),
   data=df_bison,
   model_type="binary", 
   directed = TRUE
  )

 fit_brm = bison_brm(
  bison(edge_weight(node_i_id, node_j_id)) ~ size_diff + coloration_i + coloration_j + shyness_i + shyness_j,
  fit_edge,
  df_bison,
  num_draws = 20
  )

  res_df = matrix(NA, nrow=4, ncol=9) # store results
  colnames(res_df) = c("Variable", "Median", "L", "H", "Mean", "SD", "Model", "P","Setting")

  res_df[1,] = process_samples(c(rstan::extract(fit_brm$fit,pars="b_coloration_i")$b_coloration_i),"Focal, Coloration","bison", x)
  res_df[2,] = process_samples(c(rstan::extract(fit_brm$fit,pars="b_coloration_j")$b_coloration_j),"Target, Coloration","bison", x)
  res_df[3,] = process_samples(c(rstan::extract(fit_brm$fit,pars="b_shyness_i")$b_shyness_i),"Focal, Shyness","bison", x)
  res_df[4,] = process_samples(c(rstan::extract(fit_brm$fit,pars="b_shyness_j")$b_shyness_j),"Target, Shyness","bison", x)

  res_df = data.frame(res_df)
  return(res_df)
}

res = fit_bison(model_dat, "test 1")

res # Looks about right now!
