library(rmarkdown)
# Generate appendices-------------

render("D:/OneDrive/Travail/Max Planck/Projects/ASNA_reliability/1.Codes/Appendix 1.R",
       output_format = pdf_document(),
       output_file = "D:/OneDrive/Travail/Max Planck/Projects/ASNA_reliability/2.Results/Appendices/1/Appendix 1.pdf")


render(input = "D:/OneDrive/Travail/Max Planck/Projects/ASNA_reliability/1.Codes/Appendix 2.R",
        output_format = pdf_document(),
        output_file = "D:/OneDrive/Travail/Max Planck/Projects/ASNA_reliability/2.Results/Appendices/2/Appendix 2.pdf")

# Generate analysis-------------
## No differences in sociality, no biases----------
result = simulations(Reps = 1, ncores = 1, 
                     sr_rho = 0, sr_sigma =  c(0,0), # no sender-receiver effect
                     dr_rho = 0, dr_sigma = 0, # no dyadic effect
                     N_id =  seq(30, 90, by = 10), 
                     hairy_tie_effect = seq(-0.1, 0.1, by = 0.01),
                     hairy_detect_effect = seq(0, 0, by = 0.5),
                     BISON = FALSE,
                     STRAND = T, 
                     simulate.interactions = F, 
                     int_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                     int_slope = c(-Inf,-Inf),# No effect of individuals attributes
                     blockModel = TRUE) # No block model
p = plots(result)
library(ggpubr)
ggarrange(p[[1]], p[[2]], ncol = 2, nrow = 1, common.legend = TRUE)

# Rates of false negatives and false positives
t1 = tapply(result[result$tie_effect > 0.5 | result$tie_effect < -0.5,]$`p-value`, result[result$tie_effect > 0.5 | result$tie_effect < -0.5,]$approach, function(x){sum(x >= 0.05)/length(x)})
t2 = tapply(result[result$tie_effect <= 0.5 &  result$tie_effect  >= -0.5,]$`p-value`, result[result$tie_effect<= 0.5 &  result$tie_effect  >= -0.5,]$approach, function(x){sum(x <= 0.05)/length(x)})
summary = data.frame("Approaches" = c(names(t1),  names(t2)), "Error type" = c(rep('False negatives', 7), rep('False positives', 7)) ,"Percent" = c(t1, t2))
summary$Percent = summary$Percent * 100
summary

write.csv(result, "No differences in sociality, no biases.csv", row.names = FALSE)
write.csv(summary, " No differences in sociality, no biases rates of type I and type II errors.csv", row.names = FALSE)
save.image("SIM -2 to 2.RData")

## No biases, no differences in sociality

## Differences in sociality, no biases----------
result = simulations(Reps = 1, ncores = 1, 
                     sr_rho = 0, sr_sigma =  c(0,0), # no sender-receiver effect
                     dr_rho = 0, dr_sigma = 0, # no dyadic effect
                     N_id =  seq(30, 90, by = 10), 
                     hairy_tie_effect = seq(-0.1, 0.1, by = 0.01),
                     hairy_detect_effect = seq(0, 0, by = 0.5),
                     BISON = FALSE,
                     STRAND = T, 
                     simulate.interactions = F, 
                     int_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                     int_slope = c(-Inf,-Inf),# No effect of individuals attributes
                     blockModel = TRUE) # No block model
p = plots(result)
library(ggpubr)
ggarrange(p[[1]], p[[2]], ncol = 2, nrow = 1, common.legend = TRUE)

### Rates of false negatives & false positives
t1 = tapply(result[result$tie_effect > 0.5 | result$tie_effect < -0.5,]$`p-value`, result[result$tie_effect > 0.5 | result$tie_effect < -0.5,]$approach, function(x){sum(x >= 0.05)/length(x)})
t2 = tapply(result[result$tie_effect <= 0.5 &  result$tie_effect  >= -0.5,]$`p-value`, result[result$tie_effect<= 0.5 &  result$tie_effect  >= -0.5,]$approach, function(x){sum(x <= 0.05)/length(x)})

summary = data.frame("Approaches" = c(names(t1),  names(t2)), "Error type" = c(rep('False negatives', 7), rep('False positives', 7)) ,"Percent" = c(t1, t2))
summary$Percent = summary$Percent * 100
summary

write.csv(result, "Differences in sociality, no biases.csv", row.names = FALSE)
write.csv(summary, "Differences in sociality, no biases rates of type I and type II errors.csv", row.names = FALSE)
save.image("SIM -2 to 2.RData")

## Differences in sociality, exposure bias----------
result = simulations(Reps = 1, ncores = 1, 
                     sr_rho = 0, sr_sigma =  c(0,0), # no sender-receiver effect
                     dr_rho = 0, dr_sigma = 0, # no dyadic effect
                     N_id =  seq(30, 90, by = 10), 
                     hairy_tie_effect = seq(-0.1, 0.1, by = 0.01),
                     hairy_detect_effect = seq(0, 0, by = 0.5),
                     BISON = FALSE,
                     STRAND = T, 
                     simulate.interactions = F, 
                     int_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                     int_slope = c(-Inf,-Inf),# No effect of individuals attributes
                     blockModel = TRUE) # No block model
p = plots(result)
library(ggpubr)
ggarrange(p[[1]], p[[2]], ncol = 2, nrow = 1, common.legend = TRUE)

### Rates of false negatives and false positives
t1 = tapply(result[result$tie_effect > 0.5 | result$tie_effect < -0.5,]$`p-value`, result[result$tie_effect > 0.5 | result$tie_effect < -0.5,]$approach, function(x){sum(x >= 0.05)/length(x)})
t2 = tapply(result[result$tie_effect <= 0.5 &  result$tie_effect  >= -0.5,]$`p-value`, result[result$tie_effect<= 0.5 &  result$tie_effect  >= -0.5,]$approach, function(x){sum(x <= 0.05)/length(x)})

summary = data.frame("Approaches" = c(names(t1),  names(t2)), "Error type" = c(rep('False negatives', 7), rep('False positives', 7)) ,"Percent" = c(t1, t2))
summary$Percent = summary$Percent * 100
summary

write.csv(result, "Differences in sociality, exposure bias.csv", row.names = FALSE)
write.csv(summary, "Differences in sociality, exposure bias rates of type I and type II errors.csv", row.names = FALSE)
save.image("SIM -2 to 2.RData")
## Differences in sociality, censorign bias----------
result = simulations(Reps = 1, ncores = 1, 
                     sr_rho = 0, sr_sigma =  c(0,0), # no sender-receiver effect
                     dr_rho = 0, dr_sigma = 0, # no dyadic effect
                     N_id =  seq(30, 90, by = 10), 
                     hairy_tie_effect = seq(-0.1, 0.1, by = 0.01),
                     hairy_detect_effect = seq(0, 0, by = 0.5),
                     BISON = FALSE,
                     STRAND = T, 
                     simulate.interactions = F, 
                     int_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                     int_slope = c(-Inf,-Inf),# No effect of individuals attributes
                     blockModel = TRUE) # No block model
p = plots(result)
library(ggpubr)
ggarrange(p[[1]], p[[2]], ncol = 2, nrow = 1, common.legend = TRUE)

### Rates of false negatives and false positives
t1 = tapply(result[result$tie_effect > 0.5 | result$tie_effect < -0.5,]$`p-value`, result[result$tie_effect > 0.5 | result$tie_effect < -0.5,]$approach, function(x){sum(x >= 0.05)/length(x)})
t2 = tapply(result[result$tie_effect <= 0.5 &  result$tie_effect  >= -0.5,]$`p-value`, result[result$tie_effect<= 0.5 &  result$tie_effect  >= -0.5,]$approach, function(x){sum(x <= 0.05)/length(x)})

summary = data.frame("Approaches" = c(names(t1),  names(t2)), "Error type" = c(rep('False negatives', 7), rep('False positives', 7)) ,"Percent" = c(t1, t2))
summary$Percent = summary$Percent * 100
summary

write.csv(result, "Differences in sociality, censorign bias.csv", row.names = FALSE)
write.csv(summary, "Differences in sociality, censorign bias rates of type I and type II errors.csv", row.names = FALSE)
save.image("SIM -2 to 2.RData")

