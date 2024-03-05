list.of.packages <- c("rmarkdown", "ggpubr", 'ggplot2','ggrepel', 'ANTs', 'STRAND', 'lmerTest', 'bisonR', 'rstan')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(rmarkdown)
library(ggpubr)
library(ggplot2)
library(ggrepel)
# Generate appendices-------------
render("1.Codes/Appendix 1.R",
       output_format = pdf_document(),
       output_file = "2.Results/Appendices/1/Appendix 1.pdf")


render(input = "1.Codes/Appendix 2.R",
       output_format = pdf_document(),
       output_file = "2.Results/Appendices/2/Appendix 2.pdf")

# Generate analysis-------------
source('1.Codes/2.data_simulation.R')
source('1.Codes/3.simulation.R')
## False positives rates ------------
### No differences in sociality, no biases----------
result1 = simulations(Reps = 10, ncores = 10, 
                      N_id =  seq(30, 90, by = 10), 
                      hairy_tie_effect = seq(-0.20, 0.20, by = 0.01),
                      hairy_detect_effect = seq(0, 0, by = 0.5),
                      BISON = FALSE,
                      STRAND = T, 
                      simulate.interactions = F,
                      simulate.censoring = F,
                      cens_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                      cens_slope = c(Inf,Inf),# No effect of individuals attributes
                      blockModel = TRUE) # No block model
write.csv(result1, "2.Results/Simulations/No differences in sociality, no biases.csv", row.names = FALSE)

### No differences in sociality, exposure bias----------
result2 = simulations(Reps = 10, ncores = 10, 
                      exposure_sigma = 1, 
                      N_id =  seq(30, 90, by = 10), 
                      hairy_tie_effect = seq(-0.20, 0.20, by = 0.01),
                      hairy_detect_effect = seq(-0.40, 0.40, by = 0.1),
                      BISON = FALSE,
                      STRAND = T, 
                      simulate.interactions = F,
                      simulate.censoring = F,
                      cens_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                      cens_slope = c(Inf,Inf),# No effect of individuals attributes
                      blockModel = TRUE) # No block model

write.csv(result2, "2.Results/Simulations/No differences in sociality, exposure bias.csv", row.names = FALSE)

### No differences in sociality, censoring bias----------
result3 = simulations(Reps = 10, ncores = 10, 
                      exposure_sigma = 1, 
                      N_id =  seq(30, 90, by = 10), 
                      hairy_tie_effect = seq(-0.25, 0.25, by = 0.01),
                      hairy_detect_effect = seq(0, 0, by = 0.5),
                      BISON = FALSE,
                      STRAND = T, 
                      simulate.interactions = F,
                      simulate.censoring = F,
                      cens_intercept = c(4,4), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                      cens_slope = c(-1,-1),# No effect of individuals attributes
                      blockModel = TRUE) # No block model

write.csv(result3, "2.Results/Simulations/No differences in sociality, censoring bias.csv", row.names = FALSE)



## False Negatives rates ------------
### Differences in sociality, no biases----------
result4 = simulations(Reps = 100, ncores = 100, 
                      sr_rho = 0.5, sr_sigma =  c(1.7,0.8),
                      dr_rho = 0.8, dr_sigma = 1.2,
                      exposure_sigma = 1,
                      N_id =  seq(30, 90, by = 10), 
                      hairy_tie_effect =  c(seq(-1, -0.40, by = 0.01), seq(0.40, 1, by = 0.01)),
                      hairy_detect_effect = seq(0, 0, by = 0.5),
                      BISON = FALSE,
                      STRAND = T, 
                      simulate.interactions = F, 
                      simulate.censoring = F,
                      cens_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                      cens_slope = c(-Inf,-Inf),# No effect of individuals attributes
                      blockModel = TRUE) # No block model
write.csv(result4, "2.Results/Simulations/Differences in sociality, no biases.csv", row.names = FALSE)


## Differences in sociality, exposure bias----------
result5 = simulations(Reps = 100, ncores = 100, 
                      sr_rho = 0.5, sr_sigma =  c(1.7,0.8),
                      dr_rho = 0.8, dr_sigma = 1.2,
                      exposure_sigma = 1,
                      N_id =  seq(30, 90, by = 10), 
                      hairy_tie_effect =  c(seq(-1, -0.40, by = 0.01), seq(0.40, 1, by = 0.01)),
                      hairy_detect_effect = seq(-0.80, -0.20, by = 0.5),
                      BISON = FALSE,
                      STRAND = T, 
                      simulate.interactions = F, 
                      simulate.censoring = F,
                      cens_intercept = c(Inf,Inf), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                      cens_slope = c(-Inf,-Inf),# No effect of individuals attributes
                      blockModel = TRUE) # No block model
write.csv(result5, "2.Results/Simulations/Differences in sociality, exposure bias.csv", row.names = FALSE)


## Differences in sociality, censoring bias----------
result6 = simulations(Reps = 100, ncores = 100, 
                      sr_rho = 0.5, sr_sigma =  c(1.7,0.8),
                      dr_rho = 0.8, dr_sigma = 1.2,
                      exposure_sigma = 1,
                      N_id =  seq(30, 90, by = 10), 
                      hairy_tie_effect =  c(seq(-1, -0.40, by = 0.01), seq(0.40, 1, by = 0.01)),
                      hairy_detect_effect = seq(0, 0, by = 0.5),
                      BISON = FALSE,
                      STRAND = T, 
                      simulate.interactions = F, 
                      simulate.censoring = F, 
                      cens_intercept = c(4,4), #invert log of inf = 1 of prob to observe interaction for both focal and alter
                      cens_slope = c(-1,-1),# No effect of individuals attributes
                      blockModel = TRUE) # No block model
write.csv(result6, "2.Results/Simulations/Differences in sociality, censoring bias.csv", row.names = FALSE)
save.image("2.Results/Simulations/Simulations.RData")

## Get results-------------
### Rates of false negatives and false positives------------------
error.rates <- function(result, threshold){
  
  lower = result[result$tie_effect < threshold &  result$tie_effect  > -threshold,]$p.value
  upper = result[result$tie_effect >= threshold | result$tie_effect <= -threshold,]$p.value
  
  if(length(lower) != 0){
    t1 = tapply(result[result$tie_effect < threshold &  result$tie_effect  > -threshold,]$p.value, result[result$tie_effect < threshold & result$tie_effect > -threshold,]$approach, function(x){sum(x <= 0.05)/length(x)})
  }else{t1 = NULL}
  if(length(upper) != 0){
    t2 = tapply(result[result$tie_effect >= threshold | result$tie_effect <= -threshold,]$p.value, result[result$tie_effect >= threshold | result$tie_effect <= -threshold,]$approach, function(x){sum(x >= 0.05)/length(x)})
  }else{t2 = NULL}
  
  if(!is.null(t1) & !is.null(t2)){
    summary = data.frame("Approaches" = c(names(t1),  names(t2)), 
                         "Error type" = c(rep('False positives', length(t1)),
                                          rep('False negatives', length(t2))) ,
                         "Percent" = c(t1, t2))
    
    summary$Percent = summary$Percent * 100
    return(summary)
  }
  if(!is.null(t1) != 0 & is.null(t2)){
    summary = data.frame("Approaches" = names(t1), 
                         "Error type" = rep('False positives', length(t1)) ,
                         "Percent" = t1)
    
    summary$Percent = summary$Percent * 100
    return(summary)
  }
  if(is.null(t1) & !is.null(t2)){
    summary = data.frame("Approaches" = names(t2), 
                         "Error type" = rep('False negatives', length(t2)) ,
                         "Percent" = t2)
    
    summary$Percent = summary$Percent * 100
    return(summary)
  }
}
get.rates <- function(path, threshold = 0.2){
  files = list.files(path)
  results = errors = NULL
  for(a in 1:length(files)){
    if(grepl('.csv', files[a], fixed = TRUE)){
      tmp = read.csv(paste(path,files[a], sep = '/'))
      type = gsub('.csv', '', files[a])
      tmp$type = type
      results = rbind(results, tmp)
      
      tmp =  error.rates(tmp, threshold = threshold)
      tmp$type = type
      errors = rbind(errors, tmp)
    }
  }
  rownames(errors) = NULL
  return(list(errors, results))
}
error = get.rates(path = '~/ASNA_reliability/2.Results/Simulations', threshold = 0.20)

## Plots------------------
### Rates------------------
tmp = error[[1]]
tmp$Approaches = gsub("2.","",tmp$Approaches)
tmp$Approaches = gsub("3.","",tmp$Approaches)

### False positives----------------
tmp2 = tmp[tmp$type %in% c("No differences in sociality, censoring bias", "No differences in sociality, exposure bias", "No differences in sociality, no biases"),]
tmp2$type = ifelse(tmp2$type == "No differences in sociality, censoring bias", 'Censoring biases', tmp2$type)
tmp2$type = ifelse(tmp2$type == "No differences in sociality, exposure bias", 'Exposure biases', tmp2$type)
tmp2$type = ifelse(tmp2$type == "No differences in sociality, no biases", 'No biases', tmp2$type)

tmp3 = tmp2[tmp2$type =='No biases',]
ggplot(tmp2, aes(x = type, y = Percent, group = Approaches, colour  = Approaches, label = Approaches))+geom_point(aes(size = 2))+geom_line()+facet_grid(~ Error.type) +
  geom_label_repel(data = tmp3, aes(x = type, y = Percent, colour = Approaches),
                   nudge_x = 1,
                   nudge_y = 8,
                   segment.curvature = -1e-20,
                   force = 10,
                   direction = "y",
                   hjust= 0)+
  theme(legend.position = 'None', text = element_text(size=20))+xlab("Error type")+geom_hline(yintercept = 5, linetype = "dashed", color = 'red')

tmp = error[[2]]
tmp$approach = gsub("2.","",tmp$approach)
tmp$approach = gsub("3.","",tmp$approach)

ggplot(tmp[tmp$type == "No biases", ], aes(x = tie_effect,  y = z, group = sim, label = z))+
  geom_point(aes(color = sr_rho, size = detect_effect), show.legend = TRUE, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  facet_grid( . ~ approach, space="free") +
  theme(legend.position = 'none')+
  ylab("Estimated effect size (z-score)") +
  xlab("True effect size") +
  theme(axis.text = element_text(size = 12),axis.text.x = element_text(angle=45),strip.text = element_text(size = 12))



ggplot(tmp[tmp$approach == '2.Rates weigthed' & tmp$type == 'Exposure bias',], aes(x = tie_effect, y = `p.value`, group = approach))+
  geom_point(aes(size = N_id,  color = detect_effect), show.legend = TRUE, position=position_jitter(0.2))+
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1)+
  xlab("True effect size")+
  geom_vline(xintercept = 0.25, linetype = "dashed")+
  geom_vline(xintercept = -0.25, linetype = "dashed")+
  facet_grid( . ~ approach, space="free")

### False negatives----------------
tmp = error[[1]]
tmp4 = tmp[tmp$type %in% c("Differences in sociality, censorign bias", "Differences in sociality, exposure bias", "Differences in sociality, no biases"),]
tmp4$type = ifelse(tmp4$type == "Differences in sociality, censorign bias", 'Censoring biases', tmp4$type)
tmp4$type = ifelse(tmp4$type == "Differences in sociality, exposure bias", 'Exposure biases', tmp4$type)
tmp4$type = ifelse(tmp4$type == "Differences in sociality, no biases", 'No biases', tmp4$type)

tmp5 = tmp4[tmp4$type =='No biases',]
ggplot(tmp4, aes(x = type, y = Percent, group = Approaches, colour  = Approaches, label = Approaches))+geom_point(aes(size = 2))+geom_line()+facet_grid(~ Error.type) +
  geom_label_repel(data = tmp5, aes(x = type, y = Percent, colour = Approaches),
                   nudge_x = 1,
                   nudge_y = 8,
                   segment.curvature = -1e-20,
                   force = 10,
                   direction = "y",
                   hjust= 0)+
  theme(legend.position = 'None', text = element_text(size=20))+xlab("Error type")+geom_hline(yintercept = 5, linetype = "dashed", color = 'red')

tmp = error[[2]]
tmp$approach = gsub("2.","",tmp$approach)
tmp$approach = gsub("3.","",tmp$approach)

ggplot(tmp[tmp$type == "No biases", ], aes(x = tie_effect,  y = z, group = sim, label = z))+
  geom_point(aes(color = sr_rho, size = detect_effect), show.legend = TRUE, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  facet_grid( . ~ approach, space="free") +
  theme(legend.position = 'none')+
  ylab("Estimated effect size (z-score)") +
  xlab("True effect size") +
  theme(axis.text = element_text(size = 12),axis.text.x = element_text(angle=45),strip.text = element_text(size = 12))



ggplot(tmp[tmp$approach == '2.Rates weigthed' & tmp$type == 'Exposure bias',], aes(x = tie_effect, y = `p.value`, group = approach))+
  geom_point(aes(size = N_id,  color = detect_effect), show.legend = TRUE, position=position_jitter(0.2))+
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1)+
  xlab("True effect size")+
  geom_vline(xintercept = 0.25, linetype = "dashed")+
  geom_vline(xintercept = -0.25, linetype = "dashed")+
  facet_grid( . ~ approach, space="free")

