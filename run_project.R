library(rmarkdown)
library(ggpubr)
# Generate appendices-------------
render("1.Codes/Appendix 1.R",
       output_format = pdf_document(),
       output_file = "2.Results/Appendices/1/Appendix 1.pdf")


render(input = "1.Codes/Appendix 2.R",
        output_format = pdf_document(),
        output_file = "2.Results/Appendices/2/Appendix 2.pdf")

# Generate analysis-------------
## No differences in sociality, no biases----------
source('1.Codes/2.data_simulation.R')
source('1.Codes/3.simulation.R')
result1 = simulations(Reps = 100, ncores = 100, 
                      sr_rho = 0, sr_sigma =  c(1,1), # no sender-receiver effect
                      dr_rho = 0, dr_sigma = 1, # no dyadic effect
                      exposure_sigma = 1, 
                      N_id =  seq(30, 90, by = 10), 
                      hairy_tie_effect = seq(-0.19, 0.19, by = 0.01),
                      hairy_detect_effect = seq(0, 0, by = 0.5),
                      BISON = FALSE,
                      STRAND = T, 
                      simulate.interactions = F, 
                      int_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                      int_slope = c(-Inf,-Inf),# No effect of individuals attributes
                      blockModel = TRUE) # No block model

write.csv(result1, "2.Results/Simulations/No differences in sociality, no biases.csv", row.names = FALSE)
## Differences in sociality, no biases----------
result2 = simulations(Reps = 100, ncores = 100, 
                      sr_rho = 0, sr_sigma =  c(1,1), # no sender-receiver effect
                      dr_rho = 0, dr_sigma = 1, # no dyadic effect
                      exposure_sigma = 1,
                      N_id =  seq(30, 90, by = 10), 
                      hairy_tie_effect =  c(seq(-0.40, -0.20, by = 0.01), seq(0.20, 0.40, by = 0.01)),
                      hairy_detect_effect = seq(0, 0, by = 0.5),
                      BISON = FALSE,
                      STRAND = T, 
                      simulate.interactions = F, 
                      int_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                      int_slope = c(-Inf,-Inf),# No effect of individuals attributes
                      blockModel = TRUE) # No block model
write.csv(result2, "2.Results/Simulations/Differences in sociality, no biases.csv", row.names = FALSE)


## Differences in sociality, exposure bias----------
result3 = simulations(Reps = 100, ncores = 100, 
                      sr_rho = 0, sr_sigma =  c(1,1), # no sender-receiver effect
                      dr_rho = 0, dr_sigma = 1, # no dyadic effect
                      exposure_sigma = 1,
                      N_id =  seq(30, 90, by = 10), 
                      hairy_tie_effect =  c(seq(-0.40, -0.20, by = 0.01), seq(0.20, 0.40, by = 0.01)),
                      hairy_detect_effect = seq(0, 0, by = 0.5),
                      BISON = FALSE,
                      STRAND = T, 
                      simulate.interactions = F, 
                      int_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                      int_slope = c(-Inf,-Inf),# No effect of individuals attributes
                      blockModel = TRUE) # No block model
write.csv(result3, "2.Results/Simulations/Differences in sociality, exposure bias.csv", row.names = FALSE)


## Differences in sociality, censorign bias----------
result4 = simulations(Reps = 100, ncores = 100, 
                      sr_rho = 0, sr_sigma =  c(1,1), # no sender-receiver effect
                      dr_rho = 0, dr_sigma = 1, # no dyadic effect
                      exposure_sigma = 1,
                     N_id =  seq(30, 90, by = 10), 
                     hairy_tie_effect =  c(seq(-0.40, -0.20, by = 0.01), seq(0.20, 0.40, by = 0.01)),
                     hairy_detect_effect = seq(0, 0, by = 0.5),
                     BISON = FALSE,
                     STRAND = T, 
                     simulate.interactions = F, 
                     int_intercept = c(4,4), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                     int_slope = c(-1,-1),# No effect of individuals attributes
                     blockModel = TRUE) # No block model
write.csv(result4, "2.Results/Simulations/Differences in sociality, censorign bias.csv", row.names = FALSE)
save.image("2.Results/Simulations/Simulations.RData")

# Get results-------------
## Rates of false negatives and false positives------------------
error.rates <- function(result, threshold){
  t1 = tapply(result[result$tie_effect > threshold | result$tie_effect < -threshold,]$`p-value`, result[result$tie_effect > threshold | result$tie_effect < -threshold,]$approach, function(x){sum(x >= 0.05)/length(x)})
  t2 = tapply(result[result$tie_effect <= threshold &  result$tie_effect  >= -threshold,]$`p-value`, result[result$tie_effect<= threshold &  result$tie_effect  >= -threshold,]$approach, function(x){sum(x <= 0.05)/length(x)})
  
  if(length(t1) != 0 & length(t2) != 0){
    summary = data.frame("Approaches" = c(names(t1),  names(t2)), 
                         "Error type" = c(rep('False negatives', length(t1)),
                                          rep('False positives', length(t2))) ,
                         "Percent" = c(t1, t2))
    
    summary$Percent = summary$Percent * 100
    return(summary)
  }
  if(length(t1) != 0 & length(t2) == 0){
    summary = data.frame("Approaches" = names(t1), 
                         "Error type" = rep('False negatives', length(t1)) ,
                         "Percent" = t1)
    
    summary$Percent = summary$Percent * 100
    return(summary)
  }
  if(length(t1) == 0 & length(t2) != 0){
    summary = data.frame("Approaches" = names(t2), 
                         "Error type" = rep('False positives', length(t2)) ,
                         "Percent" = t2)
    
    summary$Percent = summary$Percent * 100
    return(summary)
  }

  
}
rates1 = error.rates(result1, threshold = 0.20)
rates2 = error.rates(result2, threshold = 0.19)
rates3 = error.rates(result3, threshold = 0.19)
rates4 = error.rates(result4, threshold = 0.19)

rates1$type = "Sim1"
rates2$type = "Sim2"
rates3$type = "Sim3"
rates4$type = "Sim4"

rates = rbind(rates1, rates2, rates3, rates4)

## Plots------------------
### Rates------------------
ggplot(rates, aes(x = type, y= Percent, group = Approaches, color =Approaches))+geom_point()+geom_line()

### Coefficients------------------
plots(result1)[[1]]
plots(result2)[[1]]
plots(result3)[[1]]
plots(result4)[[1]]
