###############################
#######  Model outcomes #######
###############################
# Create an interaction matrix according to the ties probability matrix and observation bias.
# Network variables----------------------
#' @param N_id Number of individuals
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

# Individuals characteristics----------------------
#' @param individual_predictors A matrix of covariates
#' @param dyadic_predictors An array of covariates
#' @param individual_effects The effects of predictors on sender effects (row 1) and receiver effects (row 2)
#' @param dyadic_effects The effects of predictors on dyadic ties

# Exposure variables----------------------
#' @param exposure_predictors A matrix of covariates
#' @param exposure_effects A vector of slopes
#' @param exposure_sigma Variance in exposure (observations) random effects
#' @param exposure_baseline Baseline exposure (observations) rate

# Censoring variables----------------------
#' @param simulate.censoring A boolean to generate censoring bias or not
#' @param cens_intercept Intercept value
#' @param cens_slope Slope value

# Censoring variables----------------------
#' @param simulate.interactions A boolean to generate interaction, edge list or matrix of interactions (Require censoring to be siumulated)


#' @return

simulate_sbm_plus_srm_network_with_measurement_bias <- function(N_id = 30,

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
                                                               exposure_baseline = 50,
                                                                
                                                               simulate.censoring = TRUE,
                                                               cens_intercept = c(Inf,Inf),
                                                               cens_slope = c(Inf,Inf), #Change this to affect characteristics effect on detectability
                                                               simulate.interactions = TRUE){
  
    #N_id = 30
    #B = NULL
    #V = 3
    #groups=NULL
    #sr_mu = c(0,0)
    #sr_sigma = c(0.3, 1.5)
    #sr_rho = 0.6
    #dr_mu = c(0,0)
    #dr_sigma = 1
    #dr_rho = 0.7
    #individual_predictors = NULL
    #dyadic_predictors = NULL
    #individual_effects = NULL
    #dyadic_effects = NULL
    #exposure_predictors = NULL
    #exposure_effects = NULL
    #exposure_sigma = 1.9
    #exposure_baseline = 50
    #simulate.censoring = TRUE
    #cens_intercept = c(Inf,Inf)
    #cens_slope = c(Inf,Inf)
    #simulate.interactions = TRUE
    
  
  require(STRAND)
  require(ANTs)
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
  # Determine for each interaction its baseline interaction frequencies (rmvnorm2) +
  # dyads' characteristics effect on interactions sum(dyadic_effects * dyadic_predictors[i, j,])
  for ( i in 1:(N_id-1) ){
    for ( j in (i+1):N_id){
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

  # Sum dyad and interaction probabilities and create ties probabilities matrix.
  for ( i in 1:(N_id-1) ){
    for ( j in (i+1):N_id){
      p[i,j] = inv_logit( sr[i,1] + sr[j,2] + dr[i,j])
      p[j,i] = inv_logit( sr[j,1] + sr[i,2] + dr[j,i])
    }
  }

  ##################################
  #######  Model exposure bias###
  ##################################
  ideal_samps = matrix(NA, nrow=N_id, ncol=N_id) #Sample without bias.
  true_samps = matrix(NA, nrow=N_id, ncol=N_id) #Sample with bias.
  exposure_prob = rep(NA, N_id) #The probability of sample based on individual characteristics.
  ideal_exposure = rep(NA, N_id) #Observation time without bias.
  true_exposure = rep(NA, N_id) #Observation time with bias.
  exposure_offset = rep(NA, N_id) #Observation
  diag(true_samps) = 0
  diag(ideal_samps) = 0
  mu = rep(NA, N_id)
  
  ##  For each individual, determine:
  #for( i in 1:N_id){
  #  ideal_exposure[i] = rpois(1, lambda=exposure_baseline)
  #  exposure_offset[i] = rnorm(1, 0.001, exposure_sigma)/10# Reduce variance as scale is exponential
  #  mu[i] = exp(sum(exposure_effects*exposure_predictors[i,]) + exposure_offset[i])
  #  true_exposure[i] = ceiling(rgamma(n = 1, shape = mu[i]*ideal_exposure[i]/10, scale = ideal_exposure[i]/10))
  #}
  #hist(true_exposure)
  ##  For each individual, determine:
  #for( i in 1:N_id){
  #  ideal_exposure[i] = rpois(1, lambda=exposure_baseline)
  #  exposure_offset[i] = rnorm(1, 0, exposure_sigma)
  #  mu[i] = sum(exposure_effects*exposure_predictors[i,]) + 1
  #  true_exposure[i] = ceiling(rnorm(n = 1, mean = mu[i]*ideal_exposure[i], sd = abs(exposure_offset[i])))
  #}
  #hist(true_exposure)
  #
  #  For each individual, determine:
  for( i in 1:N_id){
    ideal_exposure[i] = rpois(1, lambda=exposure_baseline) # Observation baseline.
    exposure_offset[i] = rnorm(1,0,exposure_sigma) # Observation bias.
    exposure_prob[i] = inv_logit(sum(exposure_effects*exposure_predictors[i,]) + exposure_offset[i]) # Its observation probabilities based on individual characteristics.
    true_exposure[i] = rbinom(1, size=ideal_exposure[i], prob=exposure_prob[i]) # Its observations with bias.
  }

  # For each dyad, determine its sampling with and without bias according to previous steps.
  for ( i in 1:(N_id-1) ){
    for ( j in (i+1):N_id){
      ideal_samps[i,j] = sum(ideal_exposure[i] + ideal_exposure[j])
      ideal_samps[j,i] = ideal_samps[i,j]
      true_samps[i,j] = sum(true_exposure[i] + true_exposure[j])
      true_samps[j,i] = true_samps[i,j]
    }
  }
  
  ###################################
  #######  Model censoring bias###
  ###################################
  if(simulate.censoring){
    censoring = matrix(0, N_id, N_id)
    cens = rep(0, N_id)
    interactions = NULL
    for(a in 1:N_id){
      for(b in 1:N_id){
        if(a == b){next}
        if(!is.null(individual_predictors) &
           !is.infinite(cens_intercept[1]) &
           !is.infinite(cens_intercept[2])){
          p.ego =  cens_intercept[1] + cens_slope[1]*individual_predictors[a]
          p.alter = cens_intercept[2] + cens_slope[2]*individual_predictors[b]
        }else{
          p.ego =  cens_intercept[1]
          p.alter = cens_intercept[2]
        }

        censoring[a,b] = inv_logit(p.ego)*inv_logit(p.alter)
        #censoring[b,a] = inv_logit(p.ego)*inv_logit(p.alter)
        
        #if(true_exposure[a] == 0){next}
        #
        ##Probability of censor data
        #observedWithoutCensoring = rbinom(true_exposure[a] , size = 1, prob = p[a,b])
        #Nobservations = sum(observedWithoutCensoring == 1)
        #
        #observed = observedWithoutCensoring
        #observed[which(observed == 1)] = rbinom(Nobservations, size = 1, prob = censoring[a,b]) 
        #
        #cens[a] = cens[a] + sum(observedWithoutCensoring - observed)
        #cens[b] = cens[b] + sum(observedWithoutCensoring - observed)
        #interactions = rbind(interactions, data.frame('follow' = a, 'focal' = 1:true_exposure[a], 'sender' = a,  'receiver' = b, 'interaction' = observed,
        #                                             'exposure' = true_exposure[a], 'sr_p' = p[a,b], 's_i' = inv_logit(p.ego), 'r_i' = inv_logit(p.alter)))
      }
    }

    diag(censoring) = 0
  }else{
    censoring = matrix(1, N_id, N_id)
    diag(censoring) = 0
    interactions = NULL
  }


  ###############################
  #######  Model outcomes #######
  ###############################
  # Create an interaction matrix according to the ties probability matrix and observation bias.
  for ( i in 1:(N_id-1) ){
    for ( j in (i+1):N_id){
      y_true[i,j] = rbinom( 1 , size=true_samps[i,j], prob = p[i,j]*censoring[i,j] ) # This concider undirected behaviours
      y_true[j,i] = rbinom( 1 , size=true_samps[j,i], prob = p[j,i]*censoring[j,i])
    }
  }

  diag(y_true) = 0
  diag(p) = 0
  diag(dr) = 0
  y_true[is.na(y_true)] = 0
  
  colnames(y_true) = rownames(y_true) = 1:ncol(y_true)
  interactions = NULL

  for(a in 1:N_id){
    for(b in a:N_id){
      if(a != b){
        samp = true_samps[a,b]
        if(y_true[a,b] != 0){
          interactions = rbind(interactions, data.frame('i' = a, 'j' = b,
                                                        'int' = rep(1, y_true[a,b])))
        }
        if(y_true[b,a] != 0){
          interactions = rbind(interactions, data.frame('i' = b, 'j' = a,
                                                        'int' = rep(1,y_true[b,a])))
        }else{receive = NULL}

      }
    }
  }
  #tmp = df.to.mat(interactions, 'i', 'j', weighted = 'int', sym = F)
  #tmp = tmp[order(as.numeric(colnames(tmp))),order(as.numeric(colnames(tmp)))]
  #all(tmp == y_true)
  
  
  return(list(interactions = interactions,
              network= y_true,
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
              true_exposure=true_exposure,
              censoring = censoring))

}

