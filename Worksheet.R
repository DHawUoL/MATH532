source("~/MATH532/Network model.R")

A <- makeNetwork()#("smw",NULL,1000)
R0 <- 3
model <- runSimulation(A,R0,modelType="SIS")

simulateEpidemic <- function(network = "rand", R0 = 3, N = 1000, frac = 0, seed = NULL,
                             modelType = "SIR", plotOutput = TRUE, savePlot = FALSE, filename = NULL) {
  A <- makeNetwork(network, No = N)
  df <- runSimulation(A, R0, frac = frac, plotOutput = plotOutput, modelType = modelType,
                      seed = seed, savePlot = savePlot, filename = filename)
  return(df)
}

simulateEpidemic(
  network = "smw",
  R0 = 2.5,
  seed = 123,
  modelType = "SEIR",
  savePlot = TRUE
)

#Practicals:
#Test the R0 threshold
#Loop through a network parameter
#Test k vs. k_nbr calibration
#Compare outbreaks on RG and RL for same R0, and explain why
#- hint: both are basically degree uncorrelated
#Node removal - RMG: target a single module or inter-module links only;
#PA: rank by degree