 ######################### All libraries
 library(igraph)
 library(ggplot2)
 library(cmdstanr) # remotes::install_github("stan-dev/cmdstanr")
 library(ANTs)  # devtools :: install_github('https://github.com/SebastianSosa/ANTs')
 library(mice)
 library(parallel)
 library(PlvsVltra) # devtools :: install_github('https://github.com/ctross/PlvsVltra.git')

 library(bisonR) # devtools :: install_github('https://github.com/JHart96/bisonR/tree/dev') ! Dev version is mandatory!
 library(STRAND) #  devtools :: install_github('ctross/STRAND@measurement_error')
 library(amen)
 library(asnipe)



setwd(normalizePath(paste0(getwd(), "/1.Codes")))

######################### Appendix 1
rmarkdown::render("Appendix 1.R", "pdf_document")

######################### Helper functions
source("Support_Functions.R")

########################## Replicate analyses
# Effects of dividing by sample size
source("Simulate_Figure_1.R")

# Visualize a network
source("Simulate_Figure_2.R")

# Parameter recovery for STRAND models
source("Simulate_Figure_3.R")

# Main model comparisons
source("Simulate_Figure_4.R")

# Robustness sweeps. Need a supercomputer with 100 cores.
source("Simulate_Figure_5.R")
