###### Function of Change of GV Among All Generations ######
######' @param gvLst [list] list of GV and allele frequency of all simulation of list of all breeding generations, length of nGen, nInd x nTrait matirx
######' @param gvResQuantile [integer] number of quantile (xx% quantile)
######' @param nTrait [integer] number of traits
######' @param nSimu [integer] number of simulations, = length of gvLst
######' @param selTraitNames [vector] names of traits selected for selection index
######' @param gvInitMat [matrix] nInd x nTrait matrix of GV of the initial population


CalcGvChange <- function(gvLst, gvResQuantile, nTrait, 
                         nSimu, selTraitNames, gvInitMat) {
  
  ### parameters
  nGen <- length(gvLst[[1]])
  nInd <- nrow(gvLst[[1]][[1]])
  nIndAboveQuant <- (1 - gvResQuantile) * nInd ### top n% individuals
  
  resNameVec <- c("mean", 
                  "max", 
                  paste0(gvResQuantile * 100, "%quantile"), 
                  paste0("Top", round(nIndAboveQuant), "Mean"), 
                  "variance")
  
  ### index about GV for the initial population
  ### one trait
  gvInitMeanVec <- apply(X = gvInitMat[, selTraitNames], MARGIN = 2, FUN = mean)
  gvInitMaxVec <- apply(X = gvInitMat[, selTraitNames], MARGIN = 2, FUN = max)
  gvInitQuantVec <- apply(X = gvInitMat[, selTraitNames], MARGIN = 2, 
                          FUN = quantile, probs = gvResQuantile)
  gvInitTopMeanVec <- apply(X = gvInitMat[, selTraitNames], MARGIN = 2, 
                            FUN = function(x) {
                              mean(sort(x, decreasing = TRUE)[1:nIndAboveQuant])
                            })
  gvInitVarVec <- apply(X = gvInitMat[, selTraitNames], MARGIN = 2, FUN = var)
  gvInitSdVec <- apply(X = gvInitMat[, selTraitNames], MARGIN = 2, FUN = sd)
  gvInitIndexLst <- list(array(data = rbind(matrix(data = rep(0, 8), nrow = 4), 
                                            gvInitSdVec), 
                               dim = c(length(resNameVec), nTrait, 1), 
                               dimnames = list(resNameVec, selTraitNames)))
  
  ### selection index
  selIndexInitVec <- apply(X = scale(gvInitMat)[, selTraitNames], MARGIN = 1, FUN = sum)
  
  selIndexInitMean <- mean(selIndexInitVec)
  selIndexInitMax <- max(selIndexInitVec)
  selIndexInitQuant <- quantile(selIndexInitVec, probs = gvResQuantile)
  selIndexInitTopMean <- mean(sort(selIndexInitVec, decreasing = TRUE)[1:nIndAboveQuant])
  selIndexInitVar <- var(selIndexInitVec)
  selIndexInitSd <- sd(selIndexInitVec)
  selIndexInitIndexLst <- list(c(0, 0, 0, 0, selIndexInitSd))
  
  
  ### calculate index about GV for each repetition
  resAllSimuLst <- lapply(X = gvLst, FUN = function(gvEachSimuLst) {
    
    resOneSimuArray <- array(data = NA, dim = c(length(resNameVec), nTrait, (nGen + 1), 1), 
                             dimnames = list(resNameVec, selTraitNames, paste0("G", 0:nGen)))
    resOneSimuArray[, , , 1] <- do.call(what = abind::abind, 
                                        args = c(gvInitIndexLst, lapply(X = gvEachSimuLst, FUN = function(gvEachGenMat) {
                                          
                                          resOneGenArray <- array(data = NA, dim = c(length(resNameVec), nTrait, 1), 
                                                                  dimnames = list(resNameVec, selTraitNames))
                                          ### mean
                                          eachTraitMeanVec <- apply(X = gvEachGenMat, MARGIN = 2, FUN = mean)
                                          eachTraitMeanScaledVec <- 
                                            (eachTraitMeanVec - gvInitMeanVec) / gvInitSdVec
                                          ### max
                                          eachTraitMaxVec <- apply(X = gvEachGenMat, MARGIN = 2, FUN = max)
                                          eachTraitMaxScaledVec <- 
                                            (eachTraitMaxVec - gvInitMaxVec) / gvInitSdVec
                                          ### n% quantile
                                          eachTraitQuantileVec <- apply(X = gvEachGenMat, MARGIN = 2, 
                                                                        FUN = quantile, probs = gvResQuantile)
                                          eachTraitQuantileScaledVec <- 
                                            (eachTraitQuantileVec - gvInitQuantVec) / gvInitSdVec
                                          ### mean of top n%
                                          eachTraitTopMeanVec <- apply(X = gvEachGenMat, MARGIN = 2, 
                                                                       FUN = function(x) {
                                                                         mean(sort(x, decreasing = TRUE)[1:nIndAboveQuant])
                                                                       })
                                          eachTraitTopMeanScaledVec <- (eachTraitTopMeanVec - gvInitTopMeanVec) / gvInitSdVec[selTraitNames]
                                          ### variance
                                          eachTraitVarVec <- apply(X = gvEachGenMat, MARGIN = 2, FUN = var)
                                          eachTraitVarScaledVec <- eachTraitVarVec / gvInitSdVec
                                          
                                          resOneGenArray[, , 1] <- rbind(eachTraitMeanScaledVec, 
                                                                         eachTraitMaxScaledVec, 
                                                                         eachTraitQuantileScaledVec, 
                                                                         eachTraitTopMeanScaledVec, 
                                                                         eachTraitVarScaledVec)
                                          
                                          return(resOneGenArray)
                                        })))
    
    return(resOneSimuArray)
  })
  resAllSimuArray <- do.call(what = abind::abind, args = resAllSimuLst)
  dimnames(resAllSimuArray)[[4]] <- paste0("simuNo_", 1:nSimu)
  
  ### selection index (1:1)
  selIndexAllSimuLst <- lapply(X = gvLst, FUN = function(gvEachSimuLst) {
    
    selIndexOneSimuArray <- array(data = NA, dim = c(length(resNameVec), (nGen + 1), 1), 
                                  dimnames = list(resNameVec, paste0("G", 0:nGen)))
    selIndexOneSimuArray[, , 1] <- 
      do.call(what = cbind, 
              args = c(selIndexInitIndexLst, lapply(X = gvEachSimuLst, FUN = function(gvEachGenMat) {
                
                gvSelIndex <- 
                  do.call(what = "+", 
                          args = lapply(X = 1:length(selTraitNames), FUN = function(selTNo) {
                            
                            selTraitNameNow <- selTraitNames[selTNo]
                            scaledGvOneTrait <- 
                              (gvEachGenMat[, selTraitNameNow] - gvInitMeanVec[selTraitNameNow]) / 
                              gvInitSdVec[selTraitNameNow]
                            
                            return(scaledGvOneTrait)
                          }))
                ### mean
                selIndexMean <- (mean(gvSelIndex) - selIndexInitMean) / selIndexInitSd
                ### max
                selIndexMax <- (max(gvSelIndex) - selIndexInitMax) / selIndexInitSd
                ### n% quantile
                selIndexQuant <- (quantile(x = gvSelIndex, probs = gvResQuantile) - selIndexInitQuant) / selIndexInitSd
                ### mean of top n%
                selIndexTopMean <- (mean(sort(x = gvSelIndex, decreasing = TRUE)[1:nIndAboveQuant]) - selIndexInitTopMean) / selIndexInitSd
                ### variance
                selIndexVar <- var(gvSelIndex) / selIndexInitSd
                
                resVec <- c(selIndexMean, selIndexMax, selIndexQuant, selIndexTopMean, selIndexVar)
                names(resVec) <- resNameVec
                
                return(resVec)
              })))
    
    return(selIndexOneSimuArray)
  })
  selIndexAllSimuArray <- do.call(what = abind::abind, args = selIndexAllSimuLst)
  dimnames(selIndexAllSimuArray)[[3]] <- paste0("simuNo_", 1:nSimu)
  
  ### mean and sd of all repetitions (each trait)
  resMeanOfAllSimuArray <- apply(X = resAllSimuArray, MARGIN = c(1, 2, 3), 
                                 FUN = mean) ### for plot mean
  resSdOfAllSimuArray <- apply(X = resAllSimuArray, MARGIN = c(1, 2, 3), 
                               FUN = sd) ### for error bar
  
  ### mean and sd of all repetitions (selection index)
  selIndexMeanOfAllSimuArray <- apply(X = selIndexAllSimuArray, MARGIN = c(1, 2), 
                                      FUN = mean) ### for plot mean
  selIndexSdOfAllSimuArray <- apply(X = selIndexAllSimuArray, MARGIN = c(1, 2), 
                                    FUN = sd) ### for error bar
  
  resLst <- list(oneTrait = list(All = resAllSimuArray, 
                                 mean = resMeanOfAllSimuArray, 
                                 sd = resSdOfAllSimuArray), 
                 selIndex = list(All = selIndexAllSimuArray, 
                                 mean = selIndexMeanOfAllSimuArray, 
                                 sd = selIndexSdOfAllSimuArray))
  
  return(resLst)
}