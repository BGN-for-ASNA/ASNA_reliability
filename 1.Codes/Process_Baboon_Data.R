##################################################################################### Setup
# Load libraries 
 library(STRAND)
 library(psych)

# Load data
setwd("C:\\Users\\pound\\Dropbox\\Open Papers\\Sampling Bias in Animal Networks")
d1 = read.table("RFID_data.txt", sep="\t", header=TRUE) # RFID data
d2 = read.table("OBS_data.txt", sep="\t", header=TRUE)  # Observation data

# Prune data sets in order to match the time-points when data were collected. 
# This makes the two data sets be alternative measures of the same set of real-world interactions.
d_1 = d1[which(d1$DateTime %in% d2$DateTime),]
d_2 = d2[which(d2$DateTime %in% d_1$DateTime),]

# Pull the names of the 13 Baboons
Names = sort(unique(c(d_1$i, d_1$j)))

# Prune data for the 13 focal Baboons
d_2 = d_2[which(d_2$Actor %in% Names),]
d_2b = d_2[which(d_2$Recipient %in% Names),]

# Now, we can get an estimate of censoring-rate from comparing detections of the same individuals
# using RFID versus focal follow methods.
table(c(d_1$i,d_1$j))

table(c(d_2b$Recipient))
table(c(d_2b$Actor))

Detect_Rate_Absolute = table(c(d_2b$Recipient, d_2b$Actor))/table(c(d_1$i,d_1$j))         # Absolute size of these rates is arbitrary, and depends on relative number of datapoints
Detect_Rate_Relative_1 = 3.2*Detect_Rate_Absolute                                         # 3.2 normalizes rates, so that most well-observed baboon has nearly no censoring 
Detect_Rate_Relative_2 = 2.0*Detect_Rate_Absolute                                         # 2.0 normalizes rates, so that most well-observed baboon has about 40% censoring 

#################################################### Prep data for models
################################### Collar
A_1 = matrix(0, nrow=13, ncol=13)
colnames(A_1) = Names
rownames(A_1) = Names

# Loop over rows and build adjacency matrix 
for(t in 1:length(d_1$DateTime)){
  # Keeps a tally of directed ties from RFID collars
  A_1[d_1$i[t],d_1$j[t]] = 1 + A_1[d_1$i[t],d_1$j[t]]
 }

# Data are undirected, so add to transpose to get a full matrix
CollarTies = (A_1 + t(A_1))

# Exposure is the same for all individuals, number of minutes in sample, times 3 intervals per minute
Exposure_CollarTies = matrix(length(unique(d_1$DateTime))*3, nrow=13, ncol=13) # 3 20-s intervals per min

# Print tie prob during observation
round(CollarTies/Exposure_CollarTies,3)

################################### Observation data
A_2 = matrix(0, nrow=13, ncol=13)
Dur = matrix(0, nrow=13, ncol=1)
Invis = matrix(0, nrow=13, ncol=1)
Exposure_Observation = matrix(NA, nrow=13, ncol=13)

colnames(A_2) = Names
rownames(A_2) = Names
rownames(Dur) = Names
rownames(Invis) = Names

# Loop over rows and build adjacency matrix and exposure
 for(t in 1:length(d_2$DateTime)){
  if(d_2$Actor[t] %in% Names){
   
   # Keep a running tally of total duration for each actor
   Dur[d_2$Actor[t],1] = d_2$Duration[t] + Dur[d_2$Actor[t],1]         
   
   # Keep a running tally of invisible duration for each actor
   if(d_2$Behavior[t]=="Invisible"){
     Invis[d_2$Actor[t],1] = d_2$Duration[t] + Invis[d_2$Actor[t],1]   
    }
   
   # Keep a running tally of interaction duration for each dyad (directed)
   if(d_2$Recipient[t] %in% Names){
  		A_2[d_2$Actor[t],d_2$Recipient[t]] = d_2$Duration[t] + A_2[d_2$Actor[t],d_2$Recipient[t]]
  	}
  }

 }

for(i in 1:13){
 # Directed exposure is duration minus duration-invisible, so map to all columns 
 Exposure_Observation[,i] = Dur-Invis
}

# Need undirected to match the RFID collar data, so add transpose and divide
ObservationTies = round((A_2 + t(A_2))/2,0)           
Exposure_ObservationTies = round((Exposure_Observation+t(Exposure_Observation))/2,0)

# Print tie prob during observation
round(ObservationTies/Exposure_ObservationTies,3)


############################################ Now we fit some models
############################## RFID collar model
# Prep data
outcome_collar = list(Ties = CollarTies)
exposure_collar = list(Ties = Exposure_CollarTies)
Group = rep("Group 1", 13)
block = data.frame(Group = as.factor(Group), Name = as.factor(rownames(Dur)))

model_dat_collar = make_strand_data(outcome = outcome_collar,
                             individual_covariates = NULL, 
                             block_covariates = NULL,
                             dyadic_covariates = NULL,
                             outcome_mode = "binomial",
                             exposure = exposure_collar
                             )

# Collar model
fit_collar =  fit_block_plus_social_relations_model(data=model_dat_collar,
                              block_regression = ~ 1,
                              focal_regression = ~ 1,
                              target_regression = ~ 1,
                              dyad_regression = ~  1,
                              mode="mcmc",
                              return_predicted_network = TRUE,
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 600, iter_sampling = 600,
                                                          max_treedepth = 12, adapt_delta = .98)
                              )

res_collar = summarize_strand_results(fit_collar)

pred_collar = res_collar$samples$predicted_network_sample

############################## Observation model
# Prep data
outcome_obs = list(Ties = ObservationTies)
exposure_obs = list(Ties = Exposure_ObservationTies)

model_dat_obs = make_strand_data(outcome = outcome_obs,
                             individual_covariates = NULL, 
                             block_covariates = NULL,
                             dyadic_covariates = NULL,
                             outcome_mode = "binomial",
                             exposure = exposure_obs
                             )

# Observation model
fit_obs =  fit_block_plus_social_relations_model(data=model_dat_obs,
                              block_regression = ~ 1,
                              focal_regression = ~ 1,
                              target_regression = ~ 1,
                              dyad_regression = ~  1,
                              mode="mcmc",
                              return_predicted_network = TRUE,
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 600, iter_sampling = 600,
                                                          max_treedepth = 12, adapt_delta = .98)
                              )


res_obs = summarize_strand_results(fit_obs)

pred_obs = res_obs$samples$predicted_network_sample

############################## Observation model with Censoring
model_dat_obs$sampled = c(round(Dur/60,0))
model_dat_obs$sampled_exposure = c(rep(230, 13))

# Set sample size = 20
model_dat_obs$detected_exposure = c(rep(20, 13))
model_dat_obs$detected = as.vector(round(Detect_Rate_Relative_1*20,0))

# Censoring model with correct specification                                                                
fit_cens =  fit_block_plus_social_relations_model_with_measurement_bias(data=model_dat_obs,
                              block_regression = ~ 1,
                              focal_regression = ~ 1,
                              target_regression = ~ 1,
                              sampling_regression = ~ 1,
                              censoring_regression = ~ 1,
                              dyad_regression = ~  1,
                              return_predicted_network = TRUE,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 600, iter_sampling = 600,
                                                          max_treedepth = 12, adapt_delta = .98)
)

res_cens = summarize_bsrm_results_with_measurement_bias(fit_cens)

pred_cens = res_cens$samples$predicted_network_sample


############################## Observation model with Censoring, v2
model_dat_obs_2 = model_dat_obs
model_dat_obs_2$sampled = c(round(Dur/60,0))
model_dat_obs_2$sampled_exposure = c(rep(230, 13))

# Set sample size = 20
model_dat_obs_2$detected_exposure = c(rep(20, 13))
model_dat_obs_2$detected = as.vector(round(Detect_Rate_Relative_2*20,0))

# Censoring model with correct specification                                                                
fit_cens_2 =  fit_block_plus_social_relations_model_with_measurement_bias(data=model_dat_obs_2,
                              block_regression = ~ 1,
                              focal_regression = ~ 1,
                              target_regression = ~ 1,
                              sampling_regression = ~ 1,
                              censoring_regression = ~ 1,
                              dyad_regression = ~  1,
                              return_predicted_network = TRUE,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 600, iter_sampling = 600,
                                                          max_treedepth = 12, adapt_delta = .98)
)

res_cens_2 = summarize_bsrm_results_with_measurement_bias(fit_cens_2)

pred_cens_2 = res_cens_2$samples$predicted_network_sample

########################################################
par(mfrow=c(1,4))
image(apply(pred_collar,2:3,mean))
image(apply(pred_obs,2:3,mean))
image(apply(pred_cens,2:3,mean))
image(apply(pred_cens_2,2:3,mean))

Frobenius_norm = function(x,y){
 z = sqrt(tr(t(x-y) %*% (x-y)))
 return(z)
}

Div_Collar_Perm = Div_Collar_Obs = Div_Collar_Cens = Div_Collar_Cens_2 = rep(NA, 600)

pred_collar_rs = pred_collar
pred_obs_rs = pred_obs
pred_cens_rs = pred_cens
pred_cens_2_rs = pred_cens_2

pred_perm_rs = pred_obs

for(i in 1:600){
 pred_collar_rs[i,,] = pred_collar[i,,]/sum(pred_collar[i,,])
 pred_obs_rs[i,,] = pred_obs[i,,]/sum(pred_obs[i,,])
 pred_cens_rs[i,,] = pred_cens[i,,]/sum(pred_cens[i,,])
 pred_cens_2_rs[i,,] = pred_cens_2[i,,]/sum(pred_cens_2[i,,])
 
 swap = sample(1:13)
 pred_perm_rs[i,,] = pred_obs[i,swap,swap]/sum(pred_obs[i,swap,swap])

 Div_Collar_Perm[i] = Frobenius_norm(pred_collar_rs[i,,], pred_perm_rs[i,,])
 Div_Collar_Obs[i] = Frobenius_norm(pred_collar_rs[i,,], pred_obs_rs[i,,])
 Div_Collar_Cens[i] = Frobenius_norm(pred_collar_rs[i,,], pred_cens_rs[i,,])
 Div_Collar_Cens_2[i] = Frobenius_norm(pred_collar_rs[i,,], pred_cens_2_rs[i,,])
}

df = data.frame(
  Median = c(
median(Div_Collar_Perm),
median(Div_Collar_Obs),
median(Div_Collar_Cens),
median(Div_Collar_Cens_2)),

  L = c(
HPDI(Div_Collar_Perm)[1],
HPDI(Div_Collar_Obs)[1],
HPDI(Div_Collar_Cens)[1],
HPDI(Div_Collar_Cens_2)[1]),

  H = c(
HPDI(Div_Collar_Perm)[2],
HPDI(Div_Collar_Obs)[2],
HPDI(Div_Collar_Cens)[2],
HPDI(Div_Collar_Cens_2)[2]),

  Model = c("Random", "Obs.", "Cens., 0.99", "Cens., 0.60")
  )

pd = position_dodge(0.1)
cols = c(plvs_vltra("mystic_mausoleum"),"black")

p1 = ggplot(df, aes(x=Model, y=Median, colour=Model)) + 
    geom_errorbar(aes(ymin=L, ymax=H), width=0.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd) + theme_bw() + 
  scale_color_manual(values = cols[c(3, 2, 7, 13)]) + 
  theme(legend.position="bottom") + guides(color=guide_legend(nrow=2,byrow=TRUE)) + ylab("Frobenius norm") + xlab("Model contrast") +
  theme(text=element_text(size=18)) +  theme(legend.title=element_blank()) 

ggsave("Baboons_Frob.pdf", p1, height=6.5, width=6.5)

########################################## Nets
library(GGally)
library(PlvsVltra)
library(stringr)


Net_Collar = as.matrix(apply(pred_collar_rs,2:3,mean))
rownames(Net_Collar) = colnames(Net_Collar) = str_to_title(Names)

Net_Obs = as.matrix(apply(pred_obs_rs,2:3,mean))
rownames(Net_Obs) = colnames(Net_Obs) = str_to_title(Names)

Net_Cens_1 = as.matrix(apply(pred_cens_rs,2:3,mean))
rownames(Net_Cens_1) = colnames(Net_Cens_1) = str_to_title(Names)

Net_Cens_2 = as.matrix(apply(pred_cens_2_rs,2:3,mean))
rownames(Net_Cens_2) = colnames(Net_Cens_2) = str_to_title(Names)


palf = colorRampPalette(c("#fff6eb",plvs_vltra("robin_feathers", rev=TRUE, elements=NULL, show=FALSE)) )

#Net_Collar[lower.tri(Net_Collar,diag=TRUE)] = NA
g1 = pheatmap(Net_Collar, cluster_rows = FALSE,  cluster_cols = FALSE,   legend=FALSE, col = palf(100), fontsize = 18 )
ggsave("Baboons_Collar.pdf", g1, height=6, width=6)

#Net_Obs[lower.tri(Net_Obs,diag=TRUE)] = NA
g2 = pheatmap(Net_Obs, cluster_rows = FALSE,  cluster_cols = FALSE,   legend=FALSE, col = palf(100), fontsize = 18 )
ggsave("Baboons_Obs.pdf", g2, height=6, width=6)

#Net_Cens_1[lower.tri(Net_Cens_1,diag=TRUE)] = NA
g3 = pheatmap(Net_Cens_1, cluster_rows = FALSE,  cluster_cols = FALSE,   legend=FALSE, col = palf(100), fontsize = 18)
ggsave("Baboons_ME_1.pdf", g3, height=6, width=6)

#Net_Cens_2[lower.tri(Net_Cens_2, diag=TRUE)] = NA
g4 = pheatmap(Net_Cens_2, cluster_rows = FALSE,  cluster_cols = FALSE,   legend=FALSE, col = palf(100), fontsize = 18 )
ggsave("Baboons_ME_2.pdf", g4, height=6, width=6)



############################################################################################### Eigen
g = graph_from_adjacency_matrix(as.matrix(apply(pred_collar_rs,2:3,mean)), weighted=TRUE)
ord = order(eigen_centrality(g)$vector)
ev1 = cumsum(eigen_centrality(g)$vector[ord])

g = graph_from_adjacency_matrix(as.matrix(apply(pred_obs_rs,2:3,mean)), weighted=TRUE)
ev2 = cumsum(eigen_centrality(g)$vector[ord])

g = graph_from_adjacency_matrix(as.matrix(apply(pred_cens_rs,2:3,mean)), weighted=TRUE)
ev3 = cumsum(eigen_centrality(g)$vector[ord])

g = graph_from_adjacency_matrix(as.matrix(apply(pred_cens_2_rs,2:3,mean)), weighted=TRUE)
ev4 = cumsum(eigen_centrality(g)$vector[ord])

ev1 = ev1/max(ev1)
ev2 = ev2/max(ev2)
ev3 = ev3/max(ev3)
ev4 = ev4/max(ev4)

cols = c(plvs_vltra("mystic_mausoleum"),"black")

df2 = data.frame(Centrality = c(ev1, ev2, ev3, ev4), Model = rep(c("RFID", "Obs.", "Cens., 0.99", "Cens., 0.60"), each = 13), Individual = rep(1:13, 4))

p2 = ggplot(df2, aes(x=Individual, y=Centrality, group=Model)) +
  geom_line(aes(color=Model))+
  geom_point(aes(color=Model)) + theme_bw() + 
  scale_color_manual(values = cols[c(3, 2, 7, 13)]) + 
  theme(legend.position="bottom") + guides(color=guide_legend(nrow=2,byrow=TRUE)) + ylab("Cumulative eigenvector centrality") + xlab("Individuals sorted by RFID eigenvector centrality") +
  theme(text=element_text(size=18)) +  theme(legend.title=element_blank())

ggsave("Baboons_Centrality.pdf", p2, height=6.5, width=6.5)
