# Y_D = The state variable of the DPM
# Y_G = The state variable of  Grain
# Y_W = The state variable of Wheat
# Y_P = The state variable of Pasture
# Y_IOM = The state variable of the IOM
# Y_B = The state variable of Microbial biomass


# obs_R = The observation value

# Path for calling dataset (should be changed)
#setwd("/home/n10159495")
set.seed(12433)
library(invgamma)
library(truncnorm)
library(mvtnorm)

# Initial values of the model parameters
params <- c(0.065, 0.45, 0.41, 1.1, 0.1, 2.1, 0.54, 2.62, 1.35, 0.5, 0.42, 0.006, 0.2, 0.01, 0.042, 0.025, 0.023, 0.133, 0.067, 0.01, 12, 60, 60, 60, 0.096, 0.154, 0.1, 0.743, 0.076, 0.047, 0.86, 0.9, 0.1, 0.424, 0.003, 0.84, 0.5, 0.5, 0.54, 0.11, 0.5, 0.5, 0.01, 0.5, 0.5, 0.2, 2, 0.023, 0.133, 0.05, 0.05)
names(params) <- c("KD", "c", "r_W", "r_P", "p", "h", "mu_G", "mu_W", "mu_P", "rho_G", "rho_P", "s_2nu", "s_2G", "s_2W", "s_2P", "s_2TOC", "s_2eG", "s_2eW", "s_2eP", "s_2eIOM", "Y_IOM", "Y_D1", "Y_D2", "Y_D3", "s_2B", "KB", "pi_DB", "pi_BB", "KR", "KH", "PD", "s_2POC", "s_2eH", "s_2R", "s_2H", "pi_DH", "pi_RH", "pi_HH", "pi_BH", "pi_RB", "pi_HB", "r_S", "s_2S", "mu_GS", "rho_GS", "s_2GS", "h_S", "s_2eGS", "s_2eS", "pi_BC", "pi_CB") 

# Number of particles required for particle filter
numParticles = 150

# Number of the MCMC samples and burn-in
numSamples = 200000
burnIn = 100000

# Names of observation error parameters
observationErrorParamNames <- c("s_2TOC", "s_2eG", "s_2eW", "s_2eGS", "s_2eS", "s_2eIOM", "s_2eH", "s_2POC")

# Names of state variables
stateNames <- c("Y_D", "Y_G", "Y_W", "Y_GS", "Y_S", "Y_IOM", "Y_B", "Y_H", "Y_R")

# Call the dataset
data <- read.csv(file = "BrigalowDataSet.csv", header = TRUE, stringsAsFactors = FALSE)
data[(data[, "Crop"] == "w"), "Crop"] <- "Wheat"
data[(data[, "Crop"] == "f"), "Crop"] <- "Fallow"
data[(data[, "Crop"] == "s"), "Crop"] <- "Sorghum"
data[(data[, "Crop"] == "cleared"), "Crop"] <- "cleared"
data[(data[, "Crop"] == "s1"),"TDM"] <- data[(data[, "Crop"] == "s1"),"TDM"] + data[(data[, "Crop"] == "s2"),"TDM"]
data[(data[, "Crop"] == "s1"),"Grain"] <- data[(data[, "Crop"] == "s1"),"Grain"] + data[(data[, "Crop"] == "s2"),"Grain"]
data[(data[, "Crop"] == "s1" & data[, "year"] %in% 1995), "Crop"] <- "Sorghum"

data <- data[(data[, "Crop"] != "s2"),]

# Number of fields
numFields = 3
observationMatrix_list = list()
forcingMatrix_list = list()

for(i in 1:numFields)
{
  thisFieldInds <- which(data[, "field"] == i)
  thisFieldData <- data[thisFieldInds, c("TOC", "Grain", "TDM", "Grain", "TDM", "IOM", "HUM", "POC")]
  colnames(thisFieldData) <- c("TOC", "Grain_W", "Wheat", "Grain_S", "Sorghum", "IOM", "HUM", "POC")
  thisFieldForcing <- data[thisFieldInds, c("year", "Crop")]
  
  
  #Set the years that were not Sorghum to have missing TDM for Pasture
  noPastureInds <- (thisFieldForcing[, "Crop"] == "Wheat" | thisFieldForcing[, "Crop"] == "Fallow")
  if(length(noPastureInds) > 0)
  {
    thisFieldData[noPastureInds, "Sorghum"] <- NA
  }
  
  #Set the years that were not wheat to have missing TDM for Wheat
  noPastureInds <- (thisFieldForcing[, "Crop"] == "Sorghum" | thisFieldForcing[, "Crop"] == "Fallow")
  if(length(noPastureInds) > 0)
  {
    thisFieldData[noPastureInds, "Wheat"] <- NA
  }
  
  #Set the years that were not wheat to have missing Grain for Wheat
  noPastureInds <- (thisFieldForcing[, "Crop"] == "Grain_S"| thisFieldForcing[, "Crop"] == "Sorghum" | thisFieldForcing[, "Crop"] == "Fallow")
  if(length(noPastureInds) > 0)
  {
    thisFieldData[noPastureInds, "Grain_W"] <- NA
  }
  
  #Set the years that were not Sorghum to have missing Grain for Pasture
  noPastureInds <- (thisFieldForcing[, "Crop"] == "Grain_W" | thisFieldForcing[, "Crop"] == "Wheat" | thisFieldForcing[, "Crop"] == "Fallow")
  if(length(noPastureInds) > 0)
  {
    thisFieldData[noPastureInds, "Grain_S"] <- NA
  }
  
  observationMatrix_list[[i]] = thisFieldData[,]
  forcingMatrix_list[[i]] = thisFieldForcing[,]
  
}

#A stable function for computing log(sum(exp(x)))
logsumexp <- function(x){
  myMax <- max(x)
  x <- x - myMax
  return (log(sum(exp(x))) + myMax)
}


# A function to compute the carbon inputs
computeCarbonInputs <- function(params, cropType, wheatTDM, grainWTDM, SorghumTDM, grainSTDM)
{
  inputs <- rep(0, length(wheatTDM))
  ################# "Wheat" ###############################
  if(cropType == "Wheat")
  {
    inputs <- params["c"]*(wheatTDM - grainWTDM) + params["c"]*params["r_W"]*wheatTDM
  }
  else if(cropType == "Wheat For Hay")
  {
    inputs <- params["c"]*params["p"]*wheatTDM + params["c"]*params["r_W"]*wheatTDM
  }
  ################# "Sorghum" ###############################
  if(cropType == "Sorghum")
  {
    inputs <- params["c"]*(SorghumTDM - grainSTDM) + params["c"]*params["r_S"]*SorghumTDM
  }
  else if(cropType == "Sorghum For Hay")
  {
    inputs <- params["c"]*params["p"]*SorghumTDM + params["c"]*params["r_S"]*SorghumTDM
  }
  
  return(inputs)
}

# Auxiliary random numbers used in the CPM algorithm
tMax <- nrow(observationMatrix_list[[1]])
samp <- tMax*numParticles + tMax
RN <- rnorm(samp, mean=0, sd=1)
R1 <- RN


# Transition process model
transitionProcessModel <- function(params, stateVectors, forcingVariables, RV)
{
  #stateVector is a matrix with dim(state) columns and numParticles rows
  
  newStateVectors <- stateVectors
  cropType <- forcingVariables[, "Crop"]
  
  newStateVectors["Y_IOM", ] <- newStateVectors["Y_IOM", ]
  CarbonInput <- computeCarbonInputs(params, cropType, newStateVectors["Y_W", ], newStateVectors["Y_G", ], newStateVectors["Y_S", ], newStateVectors["Y_GS", ])
  
  # In case we have negative amount of carbon which it should not be, hence, remove inputs
  # as they have negative values
  mm1 = which(newStateVectors["Y_D", ]*exp(-1*params["KD"]) + CarbonInput + stateVectors["Y_B", ]*(1-exp(-1*params["KB"])) * params["pi_BC"] < 0)
  if (length(mm1) > 0)
  {
    CarbonInput[mm1] <- 0
  }
  
  carbonMean <- newStateVectors["Y_D", ]*exp(-1*params["KD"]) + CarbonInput + stateVectors["Y_B", ]*(1-exp(-1*params["KB"])) * params["pi_BC"]    
  
  newStateVectors["Y_D", ] <- exp(log(carbonMean) + sqrt(params["s_2nu"])*RV)
  
  
  
  mm2 = which(stateVectors["Y_B", ]*exp(-1*params["KB"]) + stateVectors["Y_D", ]*(1 - exp(-1*params["KD"])) * params["pi_CB"] + stateVectors["Y_B", ]*(1 - exp(-1*params["KB"])) * params["pi_BB"]< 0)
  if (length(mm2) > 0)
  {
    
    stateVectors["Y_B", mm2] <- stateVectors["Y_B", ]*exp(-1*params["KB"])
  }
  
  
  BMean <- stateVectors["Y_B", ]*exp(-1*params["KB"]) + stateVectors["Y_D", ]*(1 - exp(-1*params["KD"])) * params["pi_CB"]
  + stateVectors["Y_B", ]*(1 - exp(-1*params["KB"])) * params["pi_BB"]
  
  
  
  newStateVectors["Y_B", ] <- exp(log(BMean) + sqrt(params["s_2B"])*RV)
  
  
  return(newStateVectors)
}


# Initial values of the process model
initialiseProcessModel <- function(params, Mean, numParticles)
{
  initialStateVector <- matrix(NA, nrow = 3, ncol =  numParticles)
  rownames(initialStateVector) <- c("Y_D", "Y_IOM", "Y_B")
  initialStateVector["Y_D",] <- rtruncnorm(numParticles, a = 0, b = Inf, mean = Mean, sd = 15)                
  initialStateVector["Y_IOM",] <- rep(params["Y_IOM"],numParticles)
  initialStateVector["Y_B",] <- runif(1,0,0.05)*(initialStateVector["Y_D",] + initialStateVector["Y_IOM",])
  
  return(initialStateVector)
}


# The log-likelihood function
logLikelihoodFunction <- function(params, obs_t, y_est, observationErrorParamNames)
{
  numParticles <- nrow(y_est)
  # There is no observation for "Y_B" therefore to match "i" with the dim observations, numStates should be ncol(y_est) - 1
  numStates <- ncol(y_est) - 1
  logLike <- rep(0, numParticles)
  
  for(i in 1:numStates)
  {
    if(!is.na(obs_t[i]))
    {
      if(i == 1)
      {
        LogDens <- dlnorm(as.numeric(obs_t[1]), log(y_est[,1] + y_est[,2] + y_est[,3]), 
                          sqrt(params[observationErrorParamNames[1]]),
                          log = TRUE)
      }
      else if(i == 2)
      {
        LogDens <- dlnorm(as.numeric(obs_t[2]), log(y_est[,2]), 
                          sqrt(params[observationErrorParamNames[6]]),   
                          log = TRUE)
      }
      logLike <- logLike + LogDens
    }
  } 
  
  
  return(logLike)
}

# Proposal values for inital state variable  
rProposalY <- function(currentY)
{
  theta_propY = currentY
  for (i in 1:length(currentY))
  {
    theta_propY[i] = rtruncnorm(1, a = 0, b = Inf, mean = currentY[i], sd = 2)
  }
  
  return(theta_propY)
}

# Proposal of the model parameters
rProposal <- function(currentTheta)
{
  
  theta_prop = currentTheta
  theta_prop[1] = rnorm(1, currentTheta[1], 0.001*10)
  theta_prop[2] = rtruncnorm(1, a = 0, b = 1, mean = currentTheta[2], sd = 0.005)
  theta_prop[3] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[3], sd = 0.05)
  
  theta_prop[5] = rtruncnorm(1, a = 0, b = 1, mean = currentTheta[5], sd = 0.005)
  theta_prop[6] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[6], sd = 0.05*2)
  theta_prop[7] = rnorm(1, currentTheta[7], 0.05)
  theta_prop[8] = rnorm(1, currentTheta[8], 0.2)
  
  theta_prop[10] = rtruncnorm(1, a = -1, b = 1, mean = currentTheta[10], sd = 0.05*5)
  
  theta_prop[12] = rtruncnorm(1, a  =0, b = Inf, mean = currentTheta[12], sd = 0.001)
  theta_prop[13] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[13], sd = sqrt(currentTheta[13]/20))
  theta_prop[14] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[14], sd = 0.001*10)
  
  
  theta_prop[21] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[21], sd = 0.05)
  
  theta_prop[22] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[22], sd = 2)
  theta_prop[23] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[23], sd = 2)
  theta_prop[24] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[24], sd = 2)
  
  theta_prop[25] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[25], sd = 0.01*10)
  theta_prop[26] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[26], sd = 0.09)
  theta_prop[27] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[27], sd = 0.05*2)
  theta_prop[28] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[28], sd = 0.1)
  theta_prop[29] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[29], sd = 0.05*2)
  
  theta_prop[42] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[42], sd = 0.05)
  theta_prop[43] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[43], sd = 0.001*10)
  theta_prop[44] = rnorm(1, currentTheta[44], 0.05*2)
  theta_prop[45] = rtruncnorm(1, a = -1, b = 1, mean = currentTheta[45], sd = 0.05*5)
  theta_prop[46] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[46], sd = sqrt(currentTheta[46]/20))
  theta_prop[47] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[47], sd = 0.05*5)
  theta_prop[50] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[50], sd = 0.05*2)
  theta_prop[51] = rtruncnorm(1, a = 0, b = Inf, mean = currentTheta[51], sd = 0.05)
  return(theta_prop)
}

# Log density propsal of the intials of the state variables used in MH ratio
logDensityProposalY <- function(fromY, toY)
{
  d_Y_prop = rep(0, length(fromY))
  d_Y_prop = log(dtruncnorm(toY, a = 0, b = Inf, mean = fromY, sd = 2))
  
  return(sum(d_Y_prop))
}


# Log density propsal of the model parameters used in MH ratio
logDensityProposal <- function(fromTheta, toTheta)
{
  d_theta_prop = rep(0, length(fromTheta))
  d_theta_prop[1] = dnorm(toTheta[1], fromTheta[1], 0.001*10, log = TRUE)
  d_theta_prop[2] = log(dtruncnorm(toTheta[2], a = 0, b = 1, mean = fromTheta[2] , sd = 0.005))
  d_theta_prop[3] = log(dtruncnorm(toTheta[3], a = 0, b = Inf, mean = fromTheta[3], sd = 0.05))
  
  d_theta_prop[5] = log(dtruncnorm(toTheta[5], a = 0, b = 1, mean = fromTheta[5], sd = 0.005))
  d_theta_prop[6] = log(dtruncnorm(toTheta[6], a = 0, b = Inf, mean = fromTheta[6], sd = 0.05*2))
  d_theta_prop[7] = dnorm(toTheta[7], fromTheta[7], 0.05, log = TRUE)
  d_theta_prop[8] = dnorm(toTheta[8], fromTheta[8], 0.2, log = TRUE)
  
  d_theta_prop[10] = log(dtruncnorm(toTheta[10], a = -1, b = 1, mean = fromTheta[10], sd = 0.05*5))
  
  d_theta_prop[12] = log(dtruncnorm(toTheta[12], a = 0, b = Inf, mean = fromTheta[12], sd = 0.001))
  d_theta_prop[13] = log(dtruncnorm(toTheta[13], a = 0, b = Inf, mean = fromTheta[13], sd = sqrt(fromTheta[13]/20)))
  d_theta_prop[14] = log(dtruncnorm(toTheta[14], a = 0, b = Inf, mean = fromTheta[14], sd = 0.001*10))
  
  
  d_theta_prop[21] = log(dtruncnorm(toTheta[21], a = 0, b = Inf, mean = fromTheta[21], sd = 0.05))
  d_theta_prop[22] = log(dtruncnorm(toTheta[22], a = 0, b = Inf, mean = fromTheta[22], sd = 2))
  d_theta_prop[23] = log(dtruncnorm(toTheta[23], a = 0, b = Inf, mean = fromTheta[23], sd = 2))
  d_theta_prop[24] = log(dtruncnorm(toTheta[24], a = 0, b = Inf, mean = fromTheta[24], sd = 2))
  d_theta_prop[25] = log(dtruncnorm(toTheta[25], a = 0, b = Inf, mean = fromTheta[25], sd = 0.01*10))
  d_theta_prop[26] = log(dtruncnorm(toTheta[26], a = 0, b = Inf, mean = fromTheta[26], sd = 0.09))
  d_theta_prop[27] = log(dtruncnorm(toTheta[27], a = 0, b = Inf, mean = fromTheta[27], sd = 0.05*2))
  d_theta_prop[28] = log(dtruncnorm(toTheta[28], a = 0, b = Inf, mean = fromTheta[28], sd = 0.1))
  d_theta_prop[29] = log(dtruncnorm(toTheta[29], a = 0, b = Inf, mean = fromTheta[29], sd = 0.05*2))
  
  
  d_theta_prop[42] = log(dtruncnorm(toTheta[42], a = 0, b = Inf, mean = fromTheta[42], sd = 0.05))
  d_theta_prop[43] = log(dtruncnorm(toTheta[43], a = 0, b = Inf, mean = fromTheta[43], sd = 0.001*10))
  d_theta_prop[44] = dnorm(toTheta[44], fromTheta[44], 0.05*2, log = TRUE)
  d_theta_prop[45] = log(dtruncnorm(toTheta[45], a = -1, b = 1, mean = fromTheta[45], sd = 0.05*5))
  d_theta_prop[46] = log(dtruncnorm(toTheta[46], a = 0, b = Inf, mean = fromTheta[46], sd = sqrt(fromTheta[46]/20)))
  d_theta_prop[47] = log(dtruncnorm(toTheta[47], a = 0, b = Inf, mean = fromTheta[47], sd = 0.05*5))
  d_theta_prop[50] = log(dtruncnorm(toTheta[50], a = 0, b = Inf, mean = fromTheta[50], sd = 0.05*2))
  d_theta_prop[51] = log(dtruncnorm(toTheta[51], a = 0, b = Inf, mean = fromTheta[51], sd = 0.05))
  
  
  return(sum(d_theta_prop))
}

# Log density prior of the intials of the state variables used in MH ratio
logDensityPriorY <- function(inintialValueY)
{
  d_Y_prior = rep(0, length(inintialValueY))
  d_Y_prior = log(dtruncnorm(inintialValueY, a = 0, b = Inf, mean = 60, sd = 15))
  
  return(sum(d_Y_prior))
}

# Log density prior of the intials of the model parameters used in MH ratio
logDensityPrior <- function(theta)
{
  d_theta_prior = rep(0, length(theta))
  d_theta_prior[1] = dlnorm(theta[1], -2.71 , 0.127, log=TRUE)  
  d_theta_prior[2] = dnorm(theta[2], 0.45 , 0.01, log=TRUE) 
  d_theta_prior[3] = dnorm(theta[3], 0.5 , 0.067, log=TRUE) 
   
  d_theta_prior[5] = dbeta(theta[5], 89.9, 809.1, log=TRUE)
  d_theta_prior[6] = dlnorm(theta[6], 0.825, 0.36, log=TRUE)  
  d_theta_prior[7] = dnorm(theta[7], 0.46 , 1.18, log=TRUE) 
  d_theta_prior[8] = dnorm(theta[8], 1.24 , 1.12, log=TRUE) 
  
  d_theta_prior[10] = dunif(theta[10], -1, 1, log=TRUE)
  
  d_theta_prior[12] = dinvgamma(theta[12], 0.001, 0.001, log=TRUE)
  d_theta_prior[13] = dinvgamma(theta[13], 0.001, 0.001, log=TRUE)
  d_theta_prior[14] = dinvgamma(theta[14], 0.001, 0.001, log=TRUE)
  
  
  d_theta_prior[21] = log(dtruncnorm(theta[21], a = 0, b = Inf, mean = 12, sd = 2))
  d_theta_prior[22] = log(dtruncnorm(theta[22], a = 0, b = Inf, mean = 60, sd = 5))
  d_theta_prior[23] = log(dtruncnorm(theta[23], a = 0, b = Inf, mean = 60, sd = 5))
  d_theta_prior[24] = log(dtruncnorm(theta[24], a = 0, b = Inf, mean = 60, sd = 5))
  d_theta_prior[25] = dinvgamma(theta[25], 0.001, 0.001, log=TRUE)
  d_theta_prior[26] = dunif(theta[26], 0, 1, log=TRUE)
  d_theta_prior[27] = dunif(theta[27], 0, 1, log=TRUE)
  d_theta_prior[28] = dunif(theta[28], 0, 1, log=TRUE)
  d_theta_prior[29] = dunif(theta[29], 0, 1, log=TRUE)
  
  d_theta_prior[42] = dnorm(theta[42], 0.5 , 0.067, log=TRUE)
  d_theta_prior[43] = dinvgamma(theta[43], 0.01, 0.01, log=TRUE)
  d_theta_prior[44] = dnorm(theta[44], 0.42 , 1.18, log=TRUE)
  d_theta_prior[45] = dunif(theta[45], -1, 1, log=TRUE)
  d_theta_prior[46] = dinvgamma(theta[46], 0.001, 0.001, log=TRUE)
  d_theta_prior[47] = dlnorm(theta[47], 0.46, 0.6, log=TRUE) #log(dtruncnorm(theta[47], a = 0, b = Inf, mean = 0.46, sd = 1.5))
  d_theta_prior[50] = dunif(theta[50], 0, 1, log=TRUE)
  d_theta_prior[51] = dunif(theta[51], 0, 1, log=TRUE)
  
  return(sum(d_theta_prior))
}


# Kalman filter
KFtrans <- function(params,observationMatrix,numParticles, RV)
{
  # Observations related to the linear models
  obs_t <- observationMatrix[,2:5]
  
  # Number of linear states
  stateDimension <- ncol(obs_t)
  tMax <- nrow(obs_t)
  
  y <-  matrix(0,tMax, stateDimension)
  y_hat <- matrix(0,tMax, stateDimension)
  
  Z_density <- rep(0, tMax)  
  y_est3 <- array(0,dim = c(tMax, stateDimension, numParticles))
  
  
  A <- matrix(c(params["rho_G"], 0, 0, 0, params["rho_G"], 0, 0, 0, 0, 0, params["rho_GS"], 0, 0, 0, params["rho_GS"], 0), nrow = stateDimension, ncol = stateDimension, byrow = TRUE)
  B <- matrix(c(params["mu_G"]*(1- params["rho_G"]), params["mu_G"]*(1- params["rho_G"]) + log(params["h"]), params["mu_GS"]*(1- params["rho_GS"]), params["mu_GS"]*(1- params["rho_GS"]) + log(params["h_S"])), nrow = stateDimension, ncol = 1, byrow = TRUE)
  
  # Covariance matrix of the state variables
  covMatrix_y <- matrix(c(params["s_2G"], params["s_2G"], 0, 0, params["s_2G"], params["s_2W"] + params["s_2G"], 0, 0, 0, 0, params["s_2GS"], params["s_2GS"], 0, 0, params["s_2GS"], params["s_2S"] + params["s_2GS"]), nrow = stateDimension, ncol = stateDimension, byrow = TRUE)            
  
  # Covariance matrix of the measurement variables
  covMatrix_z <- matrix(c(params["s_2eG"],0,0,0, 0,params["s_2eW"],0,0, 0,0,params["s_2eGS"],0, 0,0,0,params["s_2eS"]), nrow = stateDimension, ncol = stateDimension, byrow = TRUE)
  
  
  
  p <- array(0, dim=c(nrow(covMatrix_y), ncol(covMatrix_y), tMax))
  p_hat <- array(0, dim=c(nrow(covMatrix_y), ncol(covMatrix_y), tMax))
  
  phat_zero <- matrix(c(4*params["s_2G"], 0,0,0, 0, params["s_2W"],0,0 ,0,0,4*params["s_2GS"],0, 0,0,0,params["s_2S"]), nrow = stateDimension, ncol = stateDimension, byrow = TRUE)
 
  
  yhat_zero <- c(params["mu_G"] - (params["rho_G"]*params["mu_G"]), log(params["h"]) + params["mu_G"], params["mu_GS"] - (params["rho_GS"]*params["mu_GS"]), log(params["h_S"]) + params["mu_GS"])
  
  
  
  for (i in 1:tMax)
  {
    # Initialise
    if(i == 1)
    {
      y[i,] <- A%*% yhat_zero + B
      p[,,i] <- A %*% phat_zero %*% t(A) + covMatrix_y
    }
    else
    {
      y[i,] <- A %*% y_hat[i-1,] + B           #Predict next state 
      p[,,i] <- A %*% p_hat[,,i-1] %*% t(A) + covMatrix_y         #Predict next covariance
    }
    
    
    if (any(!is.na(obs_t[i,])))
    {
      count <- which(as.numeric(obs_t[i,]) %in% NA)
      
      C = diag(stateDimension) 
      C <- C[-count,]
      
      
      
      Zt <- obs_t[i, -count]
      Zt <- as.matrix(log(as.numeric(Zt)))
      
      Zt_edit <- t(Zt)%*%C
      
      
      Rt_edit <- C%*%covMatrix_z%*%t(C)
      
      Pt <- p[,,i]
      Yt <- y[i,]
      
      Kt <-  Pt %*% t(C)%*%(solve(Rt_edit + C%*%Pt%*%t(C)))
      
      
      
      # Update the state estimate
      Yhat_t <- Yt + Kt%*%t((Zt_edit[,-count] - Yt%*%t(C)))
      y_hat[i,] <- Yhat_t
      
      
      # Update covariance estimation
      Phat_t <- (diag(stateDimension) - Kt%*%C)%*%Pt 
      p_hat[,,i] <- Phat_t
      
      
      
      # Index related to the random numbers (RV)
      ad = 1 + (i-1)*numParticles
      bd = i*numParticles
      
      # Draw a sample 
      for(j in 1:stateDimension)
      {
        y_est3[i,j,] <-  exp(y[i,j] + sqrt(p[j,j,i])*RV[ad:bd])
      }
      

      
      y_hat_edit <- y_hat[i, -count]
      p_hat_edit <- covMatrix_z[-count, -count]       
      
      # Compute the density of Z
      if (length(Zt) > 1)
      {
        eps <- Zt - C%*%Yt
        
        Z_density[i] <- dmvnorm(t(eps), matrix(0,1,2), p_hat_edit + C%*%Pt%*%t(C), log = TRUE)
        Z_density[i] <- Z_density[i] - sum(Zt)
      }
      else
      {
        eps <- Zt - C%*%Yt
        Z_density[i] <- dnorm(eps, 0, sqrt(p_hat_edit + Pt[-count,-count]), log = TRUE)
        Z_density[i] <- Z_density[i] - Zt 
      }
      
    }
    else
    {
      y_hat[i,] = y[i,]
      p_hat[,,i] <- p[,,i]
      
      # Index related to the random numbers (RV)
      ad = 1 + (i-1)*numParticles
      bd = i*numParticles
      
      
      for(j in 1:stateDimension)
      {
        y_est3[i,j,] <-  exp(y[i,j] + sqrt(p[j,j,i])*RV[ad:bd])
      }
      Z_density[i] <- 0
    }
    
  }  
  
  
  
  
  output <- list()
  output[["y_est3"]] = y_est3
  output[["y_hat"]] = y_hat
  output[["p_hat"]] = p_hat
  output[["Z_density"]] = Z_density
  return(output)
}


# Bootstrap particle filter
BF <- function(params, KFtrans, numParticles, stateNames, observationMatrix, forcingMatrix, observationErrorParamNames, transitionFunction, initialiseFunction, likelihoodFunction, initialY, RV)
{
  stateDimension <- ncol(observationMatrix) - 1
  tMax <- nrow(observationMatrix)
  weight <- matrix(nrow = tMax, ncol = numParticles)
  sampledParticleIndices <- matrix(nrow = tMax, ncol = numParticles)
  y_est_new <- array(NA, dim=c(stateDimension, 1, tMax))
  
  weights <- rep(0, tMax)
  y_est<- array(NA, dim=c(stateDimension, numParticles, tMax))
  theseStateVectors <- array(NA, dim=c(stateDimension, numParticles))
  Dens3Ys <- rep(0, tMax)
  
  # Sampling initial value of non-linear sub models
  y_est[c(1,6,7), ,1] <- initialY
  
  # Sort the samples
  for(i in c(1,6,7))
  {
    sorting = order(y_est[i, , 1], decreasing = FALSE)
    y_est[i, , 1] <- y_est[i, sorting, 1]
  }
  
  # Random Numbers
  R <- RV
  
  # Call the Kalman function
  KFestimate <- KFtrans(params,observationMatrix,numParticles,R)
  
  # The Density of linear submodels
  Dens3Ys <- KFestimate$Z_density
  
  # Sampling from Kalam filter
  y_est[2:5, , 1] <- KFestimate$y_est3[1,,]
  
  
  if(any(!is.na(observationMatrix[1, ])))
  {
    
    #Compute particle weights for this timestep
    weight[1, ] <- likelihoodFunction(params, observationMatrix[1, c(1,6)], t(y_est[c(1,6,7), , 1]), observationErrorParamNames)
  }
  
  weights[1] = -log(numParticles) + logsumexp(weight[1, ]) + Dens3Ys[1]
  
  #Normalise weights
  weight[1, ] = weight[1, ] - max(weight[1, ])
  weight[1, ] = exp(weight[1, ])
  weight[1, ] = weight[1, ]/sum(weight[1, ])
  
  
  #Sample Particles using weights
  
  particleIndices <- 1:numParticles
  sampledParticleIndices[1, ] <- sample(particleIndices, size = numParticles, replace = TRUE, prob = weight[1, ])
  y_est_copy <- y_est
  j = 1
  SumW = weight[1,j]
  tim1 <- tMax*numParticles + 1
  U = pnorm(R[tim1], 0, 1)/numParticles
  for (k in 1:numParticles)
  {
    
    while(SumW < U && j < numParticles)
    {
      j = j+1
      SumW = SumW + weight[1,j]
    }
    y_est[c(1,6,7), k, 1] = y_est_copy[c(1,6,7), j, 1]
    U = U + 1/numParticles
    
  }
  
  
  for (t in 2:tMax) 
  {
    #Sampling 4 Y's from KF function for this timestep
    y_est[2:5, , t] <- KFestimate$y_est3[t,,]
    
    # Index of random numbers
    ad = 1 + (t-1)*numParticles
    bd = t*numParticles
    
    #Transition part
    theseStateVectors[c(1,6,7),] <- y_est[c(1,6,7), , t-1]
    theseStateVectors[2:5,] <- y_est[2:5, , t]
    rownames(theseStateVectors) <- stateNames[1:stateDimension]
    y_est[, , t] <- transitionFunction(params, theseStateVectors, forcingMatrix[t, ], R[ad:bd])
    
    # Sorting the samples
    for(i in c(1,6,7))
    {
      sorting = order(y_est[i, , t], decreasing = FALSE)
      y_est[i, , t] <- y_est[i, sorting, t]
    }
    
    
    
    if(any(!is.na(observationMatrix[t, ])))
    {
      
      #Compute particle weights for this timestep
      weight[t, ] <- likelihoodFunction(params, observationMatrix[t, c(1,6)], t(y_est[c(1,6,7), , t]), observationErrorParamNames)
      weights[t] <- -log(numParticles) + logsumexp(weight[t, ]) + Dens3Ys[t]
      
    }
    else 
    {
      weight[t, ] <- log(1/numParticles) 
      y_est[c(1,6,7), , t] <- y_est[c(1,6,7), , t]
      weights[t] <- 0
    }
    
    
    
    
    #Normalise weights
    weight[t, ] = weight[t, ] - max(weight[t, ])
    weight[t, ] = exp(weight[t, ])
    weight[t, ] = weight[t, ]/sum(weight[t, ])
    
    
    #Resampling step for this timestep
    particleIndices <- 1:numParticles
    sampledParticleIndices[t, ] <- sample(particleIndices, size = numParticles, replace = TRUE, prob = weight[t, ])
    y_est_copy <- y_est
    j = 1
    SumW = weight[t,j]
    
    timt <- tMax*numParticles + t
    U = pnorm(R[timt], 0, 1)/numParticles
    for (k in 1:numParticles)
    {
      while (SumW < U && j < numParticles)
      {
        j = j+1
        SumW = SumW + weight[t,j]
      }
      
      
      
      y_est[c(1,6,7), k, t] = y_est_copy[c(1,6,7), j, t]
      U = U + 1/numParticles
      
    }
    
  }
  
  # Generate indeces through multinomial resampling to get the trajectory of the state variables
  Ind <- rep(0,tMax)
  Ind[tMax] <- sample(particleIndices, size = 1, replace = TRUE, prob = weight[tMax, ])   
  y_est_new[, ,tMax] <- y_est[, Ind[tMax], tMax]
  for (i in tMax:2)
  {
    Ind[i-1] <- sampledParticleIndices[i, Ind[i]]
    y_est_new[, ,i-1] <- y_est[, Ind[i-1], i-1]
  }
  
  
  
  # The estimated log-likelihood
  loglike <- sum(weights)
  
  
  output <- list()
  output$Y_est <- y_est
  output$loglike <- loglike
  output$weight <- weight
  output$weights <- weights
  output$Y_hat_new <- y_est_new
  return(output)
}

# The CPM algorithm
MH <- function(numSamples, burnIn, numParticles, KFtrans, initialParams, stateNames, observationMatrix_list, forcingMatrix_list, observationErrorParamNames, transitionFunction, initialiseFunction, logLikelihoodFunction, rProposal, logDensityProposal, logDensityPrior, rProposalY, logDensityProposalY, logDensityPriorY, R1, R_initial)
{
  stateDimension <- 7
  numberOfBlocks <- length(observationMatrix_list)
  theta = matrix(NA, nrow = numSamples, ncol = length(initialParams))
  colnames(theta) <- names(initialParams)
  theta[1, ] <- initialParams
  storedLogLikes <- rep(NA, numSamples)
  storedLogLikesHitRun <- rep(NA, numSamples)
  storedLogLikesY <- rep(NA, numSamples)
  tMax <- nrow(observationMatrix_list[[1]])
  storedrproposeY <- array(NA, dim = c(1,numParticles, numberOfBlocks))
  
  Y_hat_new <- array(NA, dim = c(numberOfBlocks, numSamples, stateDimension, tMax))
  Y_hat_temporary <- array(NA, dim = c(numberOfBlocks, stateDimension, tMax))
  
  
  
  num1 <- tMax*numParticles + tMax
  currentR = array(NA, dim = c(num1, numSamples))
  
  currentR[,1] <- R1
  # initial values that should be estimated, Y_D 
  numInitials <- 1
  
  # For storing all initial values
  y_samples <- array(NA, dim=c(numInitials, numSamples, numParticles, numberOfBlocks))
  MeanInitial = params[22]
  Ini <- initialiseFunction(initialParams, MeanInitial, numParticles)
  rownames(Ini) <- c("Y_D", "Y_IOM", "Y_B")  
  y_samples[,1,, ] <- Ini["Y_D",]
  
  # A is a counter for counting acceptance rate
  A <- rep(0, numSamples-1)
  AA <- rep(0, numSamples-1)
  ll <- 0
  for(b in 1:numberOfBlocks)
  {
    MeanInitial = params[21+b]
    initialY1 <- initialiseFunction(initialParams, MeanInitial, numParticles)
    
    thisObservationMatrix <- observationMatrix_list[[b]]
    thisForcingMatrix <- forcingMatrix_list[[b]]
    
    
    
    bfOut <- BF(initialParams, KFtrans, numParticles, stateNames, thisObservationMatrix, thisForcingMatrix, observationErrorParamNames, transitionFunction, initialiseFunction, logLikelihoodFunction, initialY1, R1)
    # The trajectory of estimated state variables
    Y_hat_new[b, 1,,] <- bfOut$Y_hat_new
    
    ll <- ll + bfOut$loglike
  }
  storedLogLikes[1] <- ll
  storedLogLikesHitRun[1] <- ll
  storedLogLikesY[1] <- ll
  
  
  for(i in 2:numSamples)
  {
    #Propose new parameter values and calculate the important quantities
    currentTheta <- theta[i-1, ]
    names(currentTheta) <- names(initialParams)
    
    
    proposedTheta <- rProposal(currentTheta)
    names(proposedTheta) <- names(initialParams)
    proposalLogNumerator <- logDensityProposal(currentTheta, proposedTheta)
    proposalLogDenominator <- logDensityProposal(proposedTheta, currentTheta)
    logPriorNumerator <- logDensityPrior(proposedTheta)
    logPriorDenominator <- logDensityPrior(currentTheta)
    
    
    
    
    ll <- 0
    for(b in 1:numberOfBlocks)
    {
      #Propose new initial "Y_D"  and calculate the important quantities
      currentY <- y_samples[,i-1,, b]
      names(currentY) <- c("Y_D")
      
      #Generating random numbers
      R <- 0.8 * currentR[,i-1] + sqrt(1 - (0.8)^2)*rnorm(num1,mean=0,sd=1)
      
      proposedY <- rProposalY(currentY) 
      names(proposedY) <- names(currentY)
      proposalLogNumeratorY <- logDensityProposalY(currentY, proposedY)
      proposalLogDenominatorY <- logDensityProposalY(proposedY, currentY)
      logPriorNumeratorY <- logDensityPriorY(proposedY)
      logPriorDenominatorY <- logDensityPriorY(currentY)
      
      pas <- Ini
      rownames(pas) <- rownames(Ini)
      pas["Y_IOM",] <- rep(proposedTheta["Y_IOM"],numParticles)
      pas["Y_D",] <-  proposedY
      
      thisObservationMatrix <- observationMatrix_list[[b]]
      thisForcingMatrix <- forcingMatrix_list[[b]]
      
      bfOut <- BF(proposedTheta, KFtrans, numParticles, stateNames, thisObservationMatrix, thisForcingMatrix, observationErrorParamNames, transitionFunction, initialiseFunction, logLikelihoodFunction, pas, R)
      
      ll <- ll + bfOut$loglike
      # The trajectory of estimated state variables
      Y_hat_temporary[b,,] <- bfOut$Y_hat_new
      # This is to check if the mass of BIO is less than 0.05 of the total mass of SOC
      for(j in 1:nrow(Y_hat_temporary[b,,]))
      {
        
        Y_b <- Y_hat_temporary[b,7,]
        for (t in 1:tMax)
        {
          SumStateVar <- sum(Y_hat_temporary[b,,t])
          tem <- (Y_b[t]/SumStateVar)
          if(tem > 0.05)
          {
            ll <- ll - (tem - 0.05)
          }
        }
      }
      storedrproposeY[,,b] <- proposedY
    }
    
    
    
    num_r_Y = ll + logPriorNumeratorY + proposalLogNumeratorY
    denom_r_Y = storedLogLikesY[i-1] + logPriorDenominatorY + proposalLogDenominatorY
    
    r_Y <- exp(num_r_Y - denom_r_Y)
    #test the condition for acceptance of the proposed "Y_D"
    storedLogLikesY[i] <- storedLogLikesY[i-1]
    if(!is.na(r_Y))
    {
      if(!is.infinite(r_Y))
      {
        if (runif(1) < r_Y)
        {
          currentY <- storedrproposeY
          storedLogLikesY[i] <- ll
        }
      }
    }
    
    y_samples[,i,, ] <- currentY
    
    
    num_r = ll + logPriorNumerator + proposalLogNumerator
    denom_r = storedLogLikes[i-1] + logPriorDenominator + proposalLogDenominator
    
    r <- exp(num_r - denom_r)
    
    #test the condition for acceptance of the proposed theta
    storedLogLikes[i] <- storedLogLikes[i-1]
    Y_hat_new[, i,,] <-  Y_hat_new[, i-1,,]
    
    
    
    
    
    if(!is.na(r))
    {
      if(!is.infinite(r))
      {
        if (runif(1) < r)
        {
          currentTheta <- proposedTheta #[-HitRun1 ]
          R1 <- R
          storedLogLikes[i] <- ll
          
          Y_hat_new[, i,,] <- Y_hat_temporary
          A[i] <- 1
        }
      }
    }
    AcceptRate <-sum(A)/numSamples
    
    theta[i, ] = currentTheta 
    currentR[,i] = R1
    
    
    if(i %in% seq(0, numSamples, by = 100)) print(i)
  }
  
  output=list()
  output[["loglike"]] = storedLogLikes
  output[["theta"]] = theta[(burnIn + 1):numSamples, ]
  output[["AcceptRate"]] = AcceptRate
  output[["initial"]] = y_samples
  output[["Y_hat_new"]] = Y_hat_new
  
  return(output)
}


theta_post_final = MH(numSamples = numSamples, burnIn = burnIn, numParticles = numParticles, KFtrans = KFtrans, initialParams = params, stateNames = stateNames, observationMatrix_list = observationMatrix_list, forcingMatrix_list = forcingMatrix_list, observationErrorParamNames = observationErrorParamNames, transitionFunction = transitionProcessModel, initialiseFunction = initialiseProcessModel, logLikelihoodFunction = logLikelihoodFunction, rProposal = rProposal, logDensityProposal = logDensityProposal, logDensityPrior = logDensityPrior, rProposalY = rProposalY, logDensityProposalY = logDensityProposalY, logDensityPriorY = logDensityPriorY, R1 = R1, R_initial)



