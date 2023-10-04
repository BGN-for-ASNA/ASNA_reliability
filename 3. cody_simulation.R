# Network variables----------------------
#' @param N_id Number of individuals
#' @param B Tie probabilities
#' @param V  Blocking variables (subgroups within the network)
#' @param groups Subgroup IDs

# Individuals variables----------------------
#' @param sr_mu  Average sender (cell 1) and reciever (cell 2) effect log odds
#' @param sr_sigma Sender (cell 1) and reciever (cell 2) effect variances
#' @param sr_rho Correlation of sender and reciever effects

# Dyadic variables---------------------------
#' @param dr_mu  Average i to j dyad effect (cell 1) and j to i dyad effect (cell 2) log odds
#' @param dr_sigma  Variance of dyad effects
#' @param dr_rho Correlation of i to j dyad effect and j to i dyad effect

#' @param individual_predictors A matrix of covariates
#' @param dyadic_predictors An array of covariates
#' @param individual_effects  The effects of predictors on sender effects (row 1) and receiver effects (row 2)
#' @param dyadic_effects The effects of predictors on dyadic ties

# Biases variables----------------------
#' @param exposure_predictors  A matrix of covariates
#' @param exposure_effects A vector of slopes
#' @param exposure_sigma Variance in exposure random effects
#' @param exposure_baseline Baseline exposure rate

#' @return
#' A list with:
#'  - network,
#'  - tie_strength
#'  - group_ids
#'  - individual_predictors
#'  - dyadic_predictors
#'  - exposure_predictors
#'  - sr
#'  - dr
#'  - true_samps
#'  - ideal_samps
#'  - ideal_exposure
#'  - true_exposure

#' @example
#' simulate_sbm_plus_srm_network_with_measurement_bias()

simulate_sbm_plus_srm_network_with_measurement_bias = function(N_id = 99,
                                                               B = NULL,
                                                               V = 3,
                                                               groups=NULL,
                                                               sr_mu = c(0,0),
                                                               sr_sigma = c(0.3, 1.5),
                                                               sr_rho = 0.6,
                                                               dr_mu = c(0,0),
                                                               dr_sigma = 1,
                                                               dr_rho = 0.7,
                                                               individual_predictors = NULL,
                                                               dyadic_predictors = NULL,
                                                               individual_effects = NULL,
                                                               dyadic_effects = NULL,
                                                               exposure_predictors = NULL,
                                                               exposure_effects = NULL,
                                                               exposure_sigma = 1.9,
                                                               exposure_baseline = 50
                                                                ){
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
# !! True networks make reference to the network without bias.
# Create correlation matrices (aka matrixes). !! It defines the reciprocity of interactions.
Rho_sr = Rho_dr = diag(c(1,1))
Rho_sr[1,2] = Rho_sr[2,1] = sr_rho
Rho_dr[1,2] = Rho_dr[2,1] = dr_rho

# Varying effects on individuals
# !! ## Determine for each dyad its baseline interaction frequencies (rmvnorm2) +
# !! ## individuals' characteristics effect on interactions sum(individual_effects[1,] * individual_predictors[i,])
sr = matrix(NA, nrow=N_id, ncol=2)
for( i in 1:N_id){
 sr[i,] = rmvnorm2(1 , Mu=sr_mu, sigma=sr_sigma, Rho= Rho_sr)

 if(!is.null(individual_predictors)){
  sr[i,1] = sr[i,1] + sum(individual_effects[1,]*individual_predictors[i,])
  sr[i,2] = sr[i,2] + sum(individual_effects[2,]*individual_predictors[i,])
  }
 }

# Build true network
dr = p = y_true = matrix(NA, N_id, N_id)
# Loop over upper triangle and create ties from i to j, and j to i
# !! Determine for each interaction its baseline interaction frequencies (rmvnorm2) +
# !! dyads' characteristics effect on interactions sum(dyadic_effects * dyadic_predictors[i, j,])
for ( i in 1:(N_id-1) ){
  for ( j in (i+1):N_id){
    # Dyadic effects
     dr_scrap = rmvnorm2(1, Mu=c(dr_mu), sigma=rep(dr_sigma,2), Rho=Rho_dr)

     if(!is.null(dyadic_predictors)){
      dr_scrap[1] = dr_scrap[1] + sum(dyadic_effects*dyadic_predictors[i,j,])
      dr_scrap[2] = dr_scrap[2] + sum(dyadic_effects*dyadic_predictors[j,i,])
      }

     # !! ## If subgroups are declared, determine within and between group link frequencies.
     B_i_j = B_j_i = c()
     for(v in 1:V){
        B_i_j[v] =  B[[v]][groups[i,v] , groups[j,v] ]
        B_j_i[v] =  B[[v]][groups[j,v] , groups[i,v] ]
      }

     dr[i,j] = dr_scrap[1] + sum(B_i_j)
     dr[j,i] = dr_scrap[2] + sum(B_j_i)
  }
}

# !! Sum dyad and interaction probabilities and create ties probabilities matrix.
for ( i in 1:(N_id-1) ){
  for ( j in (i+1):N_id){
    p[i,j] = inv_logit( sr[i,1] + sr[j,2] + dr[i,j])
    p[j,i] = inv_logit( sr[j,1] + sr[i,2] + dr[j,i])
  }
}

###############################
#######  Model measurements ###
###############################
 ideal_samps = matrix(NA, nrow=N_id, ncol=N_id) #!! Sample without bias.
 true_samps = matrix(NA, nrow=N_id, ncol=N_id) #!! Sample with bias.
 exposure_prob = rep(NA, N_id) #!! ??
 ideal_exposure = rep(NA, N_id) #!! Observation time without bias.
 true_exposure = rep(NA, N_id) #!! Observation time with bias.
 exposure_offset = rep(NA, N_id) #!! ??

 diag(true_samps) = 0
 diag(ideal_samps) = 0


# !!  For each individual, determine:
 # !! 1) Observation baseline.
 # !! 2) Observation bias.
 # !! 3) Its observation probabilities based on individual characteristics.
 # !! 4) Its observations with bias.
 for( i in 1:N_id){
  ideal_exposure[i] = rpois(1, lambda=exposure_baseline)
  exposure_offset[i] = rnorm(1,0,exposure_sigma)
  exposure_prob[i] = inv_logit(sum(exposure_effects*exposure_predictors[i,]) + exposure_offset[i])
  true_exposure[i] = rbinom(1, size=ideal_exposure[i], prob=exposure_prob[i])
 }

# !! For each dyad, determine its sampling with and without bias according to previous steps.
  for ( i in 1:(N_id-1) ){
    for ( j in (i+1):N_id){
      ideal_samps[i,j] = sum(ideal_exposure[i] + ideal_exposure[j])
      ideal_samps[j,i] = ideal_samps[i,j]
      true_samps[i,j] = sum(true_exposure[i] + true_exposure[j])
      true_samps[j,i] = true_samps[i,j]
    }
  }

###############################
#######  Model outcomes #######
###############################
# !! Create an interaction matrix according to the ties probability matrix and observation bias.
  for ( i in 1:(N_id-1) ){
    for ( j in (i+1):N_id){
      y_true[i,j] = rbinom( 1 , size=true_samps[i,j], prob = p[i,j] )
      y_true[j,i] = rbinom( 1 , size=true_samps[j,i], prob = p[j,i] )
    }
  }

  for ( i in 1:N_id ){
    y_true[i,i] = 0
    p[i,i] = 0
    dr[i,i] = 0
  }

return(list(network=y_true,
            tie_strength=p,
            group_ids=groups,
            individual_predictors=individual_predictors,
            dyadic_predictors=dyadic_predictors,
            exposure_predictors=exposure_predictors,
            sr=sr,
            dr=dr,
            true_samps=true_samps,
            ideal_samps=ideal_samps,
            ideal_exposure=ideal_exposure,
            true_exposure=true_exposure))
}


