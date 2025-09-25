source("~/MATH532/Network model.R")

A <- makeNetwork()#("smw",NULL,1000)
R0 <- 3
model <- runSimulation(A,R0)

#Practicals:
#Test the R0 threshold
#Loop through a network parameter
#Test k vs. k_nbr calibration
#Node removal - RMG: target a single module or inter-module links only;
#PA: rank by degree