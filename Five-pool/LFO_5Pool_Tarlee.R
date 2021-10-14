logLikelihoodFunctionLFO <- function(params, obs_t, y_est, observationErrorParamNames)
{
  numParticles <- nrow(y_est)
  numStates <- ncol(y_est) - 1
  numobs <- ncol(obs_t)
  logLike <- rep(0, numParticles)
  for(i in 1:numStates)
  {
    if(!is.na(obs_t[i]))
    {
      if(i == 1)
      {
        
        logLike <- logLike + dlnorm(as.numeric(obs_t[i]), log(as.numeric(y_est["Y_D"]) + as.numeric(y_est["Y_R"]) + as.numeric(y_est["Y_H"]) + as.numeric(y_est["Y_B"]) + as.numeric(y_est["Y_IOM"])), 
                                    sqrt(as.numeric(params[observationErrorParamNames[1]])),   
                                    log = TRUE)
      }
      else if(i == 7)
      {
        logLike <- logLike + dlnorm(as.numeric(obs_t[i]), log(as.numeric(y_est["Y_D"]) + as.numeric(y_est["Y_R"]) + as.numeric(y_est["Y_B"])), 
                                    sqrt(as.numeric(params[observationErrorParamNames[7]])),   
                                    log = TRUE)
        
      }
      else if(i == 6)
      {
        logLike <- logLike + dlnorm(as.numeric(obs_t[i]), log(as.numeric(y_est["Y_H"])), 
                                    sqrt(as.numeric(params[observationErrorParamNames[6]])),   
                                    log = TRUE) 
      }
      else if(2<= i & i<= 5)
      {
        logLike <- logLike + dlnorm(as.numeric(obs_t[i]), log(as.numeric(y_est[i])), 
                                    sqrt(as.numeric(params[observationErrorParamNames[i]])),
                                    log = TRUE)
      }
    }
    
  }
  
  return(logLike)
}


LPD = array(0, dim = c(3, numSamples - burnIn))
ELPD = rep(0,7)

obs = list()
forcing = list()
for(Ts in 13:19){
  
  for(i in 1:3){
    obs[[i]] = observationMatrix_list[[i]][1:Ts-1,]
    forcing[[i]] = forcingMatrix_list[[i]][1:Ts-1,]
    
  }
  
  tMax <- nrow(obs[[1]])
  samp <- tMax*numParticles + tMax
  RN <- rnorm(samp, mean=0, sd=1)
  R1 <- RN
  
  theta_post_final = MH(numSamples = numSamples, burnIn = burnIn, numParticles = numParticles, KFtrans = KFtrans, initialParams = params, stateNames = stateNames, observationMatrix_list = obs, forcingMatrix_list = forcing, observationErrorParamNames = observationErrorParamNames, transitionFunction = transitionProcessModel, initialiseFunction = initialiseProcessModel, logLikelihoodFunction = logLikelihoodFunction, rProposal = rProposal, logDensityProposal = logDensityProposal, logDensityPrior = logDensityPrior, rProposalY = rProposalY, logDensityProposalY = logDensityProposalY, logDensityPriorY = logDensityPriorY, R1 = R1, R_initial)
  theta_post <- theta_post_final$theta
  Y_hat_new <- theta_post_final$Y_hat_new
  for(b in 1:3){
    
    obser = observationMatrix_list[[i]][1:Ts,]
    thisForcingMatrix <- forcingMatrix_list[[i]][1:Ts,]
    y_est <- colMeans(Y_hat_new[b,,,Ts])
    for(i in 1:(numSamples - burnIn))
    {
      LPD[b, i] <- logLikelihoodFunctionLFO(theta_post[i,], obser[Ts,], t(y_est), observationErrorParamNames)
    }
  }
  
  AllFs = rep(0, numSamples - burnIn)
  
  for(i in 1:(numSamples - burnIn))
  {
    AllFs[i] = LPD[1, i] + LPD[2, i] + LPD[3, i]
  }
  
  expAllFs = rep(0, (numSamples - burnIn)/30)
  expAllFs = exp(AllFs[seq(1,numSamples - burnIn,30)] - max(AllFs))
  
  ELPD[Ts-12] = log(mean(expAllFs))
  
  
}

LFO = sum(ELPD)