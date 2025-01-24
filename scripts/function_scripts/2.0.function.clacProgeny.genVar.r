### calculating the recombination frequency
###' @param linkInfoDf [data.frame] each column contains: CHROM: chromosome number, MapDis: map distance, RR: recombinationrate
###' @param k [integer] F'k' generation, number of generation of selfing
###'


CalcD <- function(linkInfoDf, k) {
  ### calculate the distance of each marker set
  myDist <- sapply(X = 1:nrow(linkInfoDf),
                   FUN = function(mrkNo) abs(linkInfoDf$MapDis[mrkNo] - linkInfoDf$MapDis))
  
  ### check the marker on the same chr or not
  myCHR1 <- do.call(what = rbind,
                    args = lapply(1:nrow(linkInfoDf),
                                  function(mrkNo) rep(linkInfoDf$CHROM[mrkNo],
                                                      nrow(linkInfoDf))))
  myCHR2 <- t(myCHR1)
  
  ### convert the distance to the recombination frequency (Haldane)
  c1 <- 0.5 * (1 - exp(-2 * (myDist / 100)))
  ### r of the markers on the diff chr is 0.5
  c1[myCHR1 != myCHR2] <- 0.5
  
  ### r for the future F"k" generations (for selfing)
  l <- k - 1
  ck <- 2 * c1 * (1 - ((0.5)^l) * (1 - 2 * c1)^l) / (1 + 2 * c1)
  
  
  ### LD parameter which is common to all cross pairs
  D_12 <- 0.25 * (1 - ck) * (1 - 2 * c1)
  # diag(D_12) <- 0.25 / 4 ### 1/16 same as UCPC
  diag(D_12) <- 0.25
  D_34 <- D_12
  
  if (k == 1) {
    
    t <- 0
    D_13 <- matrix(data = 0, nrow = nrow(myDist), ncol = nrow(myDist))
    D_14 <- D_23 <- D_24 <- D_13
    
  } else {
    
    t <- k - 2
  D_13 <- 0.25 * (1 - 2 * ck - (0.5 * (1 - 2 * c1))^(t))
  # diag(D_13) <- 0.25 / 4 ### 1/16
  diag(D_13) <- 0.25
  D_14 <- D_23 <- D_24 <- D_13
  
  }
  
  return(list(D_12 = D_12,
              D_34 = D_34,
              D_13 = D_13,
              D_14 = D_14,
              D_23 = D_23,
              D_24 = D_24))
}




### Calculating the variance-covariance between selected 2 genotypes
GenCovProgeny <- function(ObjectGenoP1, ObjectGenoP2, DLst, MarkEffect1, MarkEffect2)
{
  
  beta1_1 <- ObjectGenoP1[1, ] * MarkEffect1
  beta1_2 <- ObjectGenoP1[2, ] * MarkEffect1
  beta2_1 <- ObjectGenoP2[1, ] * MarkEffect2
  beta2_2 <- ObjectGenoP2[2, ] * MarkEffect2
  
  sigma_12 <- t(beta1_1 - beta1_2) %*% (DLst$D_12) %*% (beta1_1 - beta1_2)
  sigma_34 <- t(beta2_1 - beta2_2) %*% (DLst$D_34) %*% (beta2_1 - beta2_2)
  sigma_13 <- t(beta1_1 - beta2_1) %*% (DLst$D_13) %*% (beta1_1 - beta2_1)
  sigma_14 <- t(beta1_1 - beta2_2) %*% (DLst$D_14) %*% (beta1_1 - beta2_2)
  sigma_23 <- t(beta1_2 - beta2_1) %*% (DLst$D_23) %*% (beta1_2 - beta2_1)
  sigma_24 <- t(beta1_2 - beta2_2) %*% (DLst$D_24) %*% (beta1_2 - beta2_2)
  sigma <- sigma_12 + sigma_34 + sigma_13 + sigma_14 + sigma_23 + sigma_24
  
  return(sigma)
}

### For different populations with different marker effects
###' @param ObjectHaploP1 [array] nMrk x nPop x nHaplo, haplotype of parent1
###' @param mrkEffLst [list] list length of nTrait, each contains nMrk x nPop matrix, marker effects of each population binded
###' 
GenCovProgenyMultiPop <- function(ObjectHaploP1, ObjectHaploP2, DLst, mrkEffLst) {
  
  sigmaAlltraitVec <- sapply(X = 1:length(mrkEffLst), FUN = function(traitNo) {
    
    mrkEffMat <- mrkEffLst[[traitNo]]
    
    beta1_1 <- apply(X = (ObjectHaploP1[, , 1] * mrkEffMat), MARGIN = 1, FUN = sum)
    beta1_2 <- apply(X = (ObjectHaploP1[, , 2] * mrkEffMat), MARGIN = 1, FUN = sum)
    beta2_1 <- apply(X = (ObjectHaploP2[, , 1] * mrkEffMat), MARGIN = 1, FUN = sum)
    beta2_2 <- apply(X = (ObjectHaploP2[, , 2] * mrkEffMat), MARGIN = 1, FUN = sum)
    
    sigma_12 <- t(beta1_1 - beta1_2) %*% (DLst$D_12) %*% (beta1_1 - beta1_2)
    sigma_34 <- t(beta2_1 - beta2_2) %*% (DLst$D_34) %*% (beta2_1 - beta2_2)
    sigma_13 <- t(beta1_1 - beta2_1) %*% (DLst$D_13) %*% (beta1_1 - beta2_1)
    sigma_14 <- t(beta1_1 - beta2_2) %*% (DLst$D_14) %*% (beta1_1 - beta2_2)
    sigma_23 <- t(beta1_2 - beta2_1) %*% (DLst$D_23) %*% (beta1_2 - beta2_1)
    sigma_24 <- t(beta1_2 - beta2_2) %*% (DLst$D_24) %*% (beta1_2 - beta2_2)
    sigma <- sigma_12 + sigma_34 + sigma_13 + sigma_14 + sigma_23 + sigma_24
    
    return(sigma)
  })
  names(sigmaAlltraitVec) <- names(mrkEffLst)
  
  return(sigmaAlltraitVec)
}
