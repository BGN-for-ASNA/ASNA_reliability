#' Create an interaction matrix according to the ties probability matrix and observation bias.
#' Network variables----------------------
#' @param N_id Number of individuals
#' @param individual_predictors A matrix of covariates
#' @param dyadic_predictors A array of covariates
#' @param exposure_predictors A matrix of covariates
#' @param censoring_predictors An array of covariates

# Block variables----------------------
#' @param B Tie probabilities
#' @param V Blocking variables (subgroups within the network)
#' @param groups Subgroup IDs

# Sender/receiver variables----------------------
#' @param sr_mu Average sender (cell 1) and receiver (cell 2) effect log odds
#' @param sr_sigma Sender (cell 1) and receiver (cell 2) effect variances
#' @param sr_rho Correlation of sender and receiver effects

# Dyadic variables---------------------------
#' @param dr_mu Average i to j dyad effect (cell 1) and j to i dyad effect (cell 2) log odds
#' @param dr_sigma Variance of dyad effects
#' @param dr_rho Correlation of i to j dyad effect and j to i dyad effect

# Exposure variables---------------------------
#' @param exposure_mu Average log odds scale
#' @param exposure_sigma Variance of random effects
#' @param exposure_max Max count of observations per focal

# Censoring variables---------------------------
#' @param censoring_mu Average log odds scale
#' @param censoring_sigma Variance of random effects

# Effects----------------------
#' @param individual_effects The effects of predictors on sender effects (row 1) and receiver effects (row 2)
#' @param dyadic_effects The effects of predictors on dyadic ties
#' @param exposure_effects The effects of predictors on observation count
#' @param censoring_effects The effects of predictors on censoring

# Exposure variables----------------------
#' @param exposure_predictors A matrix of covariates
#' @param exposure_effects A vector of slopes
#' @param exposure_sigma Variance in exposure (observations) random effects
#' @param exposure_baseline Baseline exposure (observations) rate

# Censoring variables----------------------
#' @param simulate_censoring A boolean to generate censoring bias or not

# Censoring variables----------------------
#' @param simulate_interactions A boolean to generate interaction, edge list or matrix of interactions (Require censoring to be siumulated)


#' @return

simulate_sbm_plus_srm_network_with_measurement_bias = function(N_id = 30,
                                                               individual_predictors = NULL,
                                                               dyadic_predictors = NULL,                                                               
                                                               exposure_predictors = NULL,
                                                               censoring_predictors = NULL,
                                                               
                                                               B = NULL,
                                                               V = 3,
                                                               groups=NULL,

                                                               sr_mu = c(0,0),
                                                               sr_sigma = c(0.3, 1.5),
                                                               sr_rho = 0.6,

                                                               dr_mu = c(0,0),
                                                               dr_sigma = 1,
                                                               dr_rho = 0.7,

                                                               exposure_mu = 1.9,
                                                               exposure_sigma = 1.0,
                                                               exposure_max = 50,
                                                               
                                                               censoring_mu = 1.9,
                                                               censoring_sigma = 1.0,
                                                               N_trials = 20,

                                                               individual_effects = NULL,
                                                               dyadic_effects = NULL,
                                                               exposure_effects = NULL,
                                                               censoring_effects = NULL, 

                                                               simulate_censoring = TRUE,
                                                               simulate_interactions = TRUE){  
  require(STRAND)
  ###############################
  ####### Run some checks #######
  ###############################
   if(!is.null(individual_predictors)){
     if(is.null(individual_effects)){
      stop("If individual_predictors is supplied, a matching matrix of individual_effects must be supplied.")
     }

    if(length(individual_effects[1,]) != length(individual_predictors[1,])){
      stop("The number of columns of individual_effects must match that of individual_predictors.")
    }

    if(nrow(individual_effects) != 2 ){
      stop("The number of rows of individual_effects must be 2: one for sender effects, the other for receiver. Un-needed slopes can be set to 0.")
    }

    if(nrow(individual_predictors) != N_id ){
     stop("The number of rows of individual_predictors must be N_id.")
    }
   }

  ###############################
  ####### Model true network ####
  ###############################
  # True networks make reference to the network without bias.
  # Create correlation matrices (aka matrixes). It defines the reciprocity of interactions.
  Rho_sr = Rho_dr = Rho_int = diag(c(1,1))
  Rho_sr[1,2] = Rho_sr[2,1] = sr_rho
  Rho_dr[1,2] = Rho_dr[2,1] = dr_rho

  # Varying effects on individuals
  # ## Determine for each dyad its baseline interaction frequencies (rmvnorm2) +
  # ## individuals' characteristics effect on interactions sum(individual_effects[1,] * individual_predictors[i,])
  sr = matrix(NA, nrow=N_id, ncol=2)
  for(i in 1:N_id){
    sr[i,] = rmvnorm2(1 , Mu=sr_mu, sigma=sr_sigma, Rho= Rho_sr)

    if(!is.null(individual_predictors)){
      sr[i,1] = sr[i,1] + sum(individual_effects[1,]*individual_predictors[i,])
      sr[i,2] = sr[i,2] + sum(individual_effects[2,]*individual_predictors[i,])
    }
  }

  # Build network
  dr = p = y = matrix(0, N_id, N_id)
  # Loop over upper triangle and create ties from i to j, and j to i
  # Determine for each interaction its baseline interaction frequencies (rmvnorm2) +
  # dyads' characteristics effect on interactions sum(dyadic_effects * dyadic_predictors[i, j,])
  for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){
      # Dyadic effects
       dr_scrap = rmvnorm2(1, Mu=c(dr_mu), sigma=rep(dr_sigma,2), Rho=Rho_dr)

       if(!is.null(dyadic_predictors)){
         dr_scrap[1] = dr_scrap[1] + sum(dyadic_effects*dyadic_predictors[i,j,])
         dr_scrap[2] = dr_scrap[2] + sum(dyadic_effects*dyadic_predictors[j,i,])
        }

       # ## If subgroups are declared, determine within and between group link frequencies.
       B_i_j = B_j_i = c()
       for(v in 1:V){
          B_i_j[v] =  B[[v]][groups[i,v] , groups[j,v] ]
          B_j_i[v] =  B[[v]][groups[j,v] , groups[i,v] ]
       }
       dr[i,j] = dr_scrap[1] + sum(B_i_j)
       dr[j,i] = dr_scrap[2] + sum(B_j_i)
    }
  }

  # Sum dyad and interaction weights and create tie probability matrix.
  for (i in 1:(N_id-1) ){
    for (j in (i+1):N_id){
      p[i,j] = inv_logit(sr[i,1] + sr[j,2] + dr[i,j])
      p[j,i] = inv_logit(sr[j,1] + sr[i,2] + dr[j,i])
    }
  }

  ##################################
  #######  Model exposure bias   ###
  ##################################
  ideal_samps = matrix(NA, nrow=N_id, ncol=N_id) # Sample without bias.
  true_samps = matrix(NA, nrow=N_id, ncol=N_id)  # Sample with bias.
  exposure_prob = rep(NA, N_id)                  # The probability of sample based on individual characteristics.
  true_exposure = rep(NA, N_id)                  # Observation time with bias.
  exposure_offset = rep(NA, N_id)                # Observation
  exposure_factors = rep(0, N_id)               
  diag(true_samps) = 0
  diag(ideal_samps) = 0
  
  # For each individual, determine:
  for(i in 1:N_id){
        
    if(!is.null(exposure_predictors)){
    exposure_factors[i] = exposure_factors[i] + sum(exposure_effects*exposure_predictors[i,])
    }

    exposure_offset[i] = rnorm(1, 0, exposure_sigma)                                             # Observation bias.
    exposure_prob[i] = inv_logit(exposure_mu + exposure_factors[i] + exposure_offset[i])         # Its observation probabilities based on individual characteristics.
    true_exposure[i] = rbinom(1, size=exposure_max, prob=exposure_prob[i])                       # Its observations with bias.
  }

  # For each dyad, determine its sampling with and without bias according to previous steps. 
  for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){
      ideal_samps[i,j] = sum(exposure_max + exposure_max)
      ideal_samps[j,i] = ideal_samps[i,j]

      true_samps[i,j] = sum(true_exposure[i] + true_exposure[j])
      true_samps[j,i] = true_samps[i,j]
    }
  }
  
  ###################################
  #######  Model censoring bias   ###
  ###################################
   censoring_offset = rep(NA, N_id)                
   censoring_prob = rep(0, N_id) 
   censoring_factors = rep(0, N_id)               

   if(simulate_censoring == TRUE){
   # For each individual, determine:
   for(i in 1:N_id){

    if(!is.null(censoring_predictors)){
    censoring_factors[i] = censoring_factors[i] + sum(censoring_effects*censoring_predictors[i,])
    }

    censoring_offset[i] = rnorm(1, 0, censoring_sigma)                                                                 
    censoring_prob[i] = inv_logit(censoring_mu + censoring_factors[i] + censoring_offset[i])                                      
    }
     }

  eta = 1 - censoring_prob             # Flip to represent NOT censoring

  ###############################
  #######  Model outcomes #######
  ###############################
  # Create an interaction matrix according to the ties probability matrix and observation bias.
  for ( i in 1:(N_id-1) ){
    for ( j in (i+1):N_id){
      y[i,j] = rbinom( 1 , size=true_samps[i,j], prob = p[i,j]*eta[i]*eta[j] ) 
      y[j,i] = rbinom( 1 , size=true_samps[j,i], prob = p[j,i]*eta[j]*eta[i] )
    }
  }
  
  
  trials = rep(N_trials, N_id)
  detected = NULL
  for(i in 1: N_id){
    detected[i]  = rbinom(1, size = trials[i], prob = eta[i])
  }

  diag(y) = 0
  diag(p) = 0
  diag(dr) = 0
  
  colnames(y) = rownames(y) = 1:ncol(y)
  interactions = NULL

 if(simulate_interactions==TRUE){
    for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){
          temp_dat_ij = data.frame(focal=i, target=j, outcome=y[i,j], duration=true_samps[i,j])
          temp_dat_ji = data.frame(focal=j, target=i, outcome=y[j,i], duration=true_samps[j,i])

         interactions = rbind(interactions, temp_dat_ij, temp_dat_ji)
      }}
    }
  
  return(list(interactions = interactions,
              network= y,
              tie_strength=p,
              group_ids=groups,
              individual_predictors=individual_predictors,
              dyadic_predictors=dyadic_predictors,
              exposure_predictors=exposure_predictors,
              sr=sr,
              dr=dr,
              true_samps=true_samps,
              ideal_samps=ideal_samps,
              true_exposure=true_exposure,
              censoring = simulate_censoring,
              detected = detected,
              trials = trials))

}


                                                               
                                                               
                                                             
                                                                                       
                                                 

                                                                

# Test ----------------------------------------------------------------------
#V = 1           # One blocking variable
#G = 3           # Three categories in this variable
#N_id = 30       # Number of people

#clique = sample(1:3, N_id, replace=TRUE)
#B = matrix(-8, nrow=G, ncol=G)
#diag(B) = -4.5 # Block matrix

#B[1,3] = -5.9
#B[3,2] = -6.9


#                                                                
#A = simulate_sbm_plus_srm_network_with_measurement_bias(N_id = N_id, 
#                                                    B=list(B=B),
#                                                    V=V,
#                                                    groups=data.frame(clique=factor(clique)),
#                                                    individual_predictors=matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1),
#                                                    individual_effects=matrix(c(1.7, 0.3),ncol=1, nrow=2),
#                                                    dyadic_predictors = array(rnorm(N_id*N_id, 0, 1), c(N_id, N_id, 1)),
#                                                    dyadic_effects = c(0.03),
#                                                    sr_mu = c(0,0),
#                                                    sr_sigma = c(1.4, 0.8),
#                                                    sr_rho = 0.5,
#                                                    dr_mu = c(0,0),
#                                                    dr_sigma = 1.2,
#                                                    dr_rho = 0.8,
#                                                    exposure_mu = 1.9,
#                                                    exposure_sigma = 1.0,
#                                                    exposure_max = 50,
#                                                    censoring_mu = -1.9,
#                                                    censoring_sigma = 1.0,
#                                                    exposure_predictors = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1),
#                                                    censoring_predictors = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1),
#                                                    simulate_censoring = TRUE,
#                                                    simulate_interactions = TRUE,
#                                                    exposure_effects = c(-1.1),
#                                                    censoring_effects = c(-1.1)
#                                                    )
#                               

#Net = graph_from_adjacency_matrix(A$network, mode = c("directed"))
#V(Net)$color = c("turquoise4","gray13", "goldenrod3")[A$group_ids$clique]
#
#plot(Net, edge.arrow.size =0.1, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)

