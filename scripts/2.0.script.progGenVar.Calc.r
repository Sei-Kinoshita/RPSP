###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load packages #####
require(stringr)
require(pbmcapply)


##### 0.3. Define parameters #####
scriptID <- "2.0."

nGenVec <- 1:15 ### 4: cross hetero lines --> selfing 3 generation

nCoresComb <- 15

popNames <- c("S827", "S840", "S844")
nPop <- length(popNames)

selMethod <- "kMedoids" ### selection method of parents


##### 0.4. Define directories #####
haploDir <- "data/haplotype/"
linkInfoDir <- "data/linkageInfo/" ### rr, map distance
mrkEffDir <- "data/mrkEff/"
crsCombDir <- "data/crossComb/kMedoids_nClust_20/"
crsCombAllF4Dir <- "data/crossComb/F4895_allCross/"

saveDir <- paste0("result/", scriptID, "GeneticVariance_Calculation/")
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}


##### 0.5. Load Function #####
source("scripts/function_scripts/2.0.function.clacProgeny.genVar.r")




###### 1. Load Data ######
##### 1.1. Haplotype Data #####
fileNameF4Haplo <- paste0(haploDir, "0.1.F4895_Haplotype.rds")
f4HaploLst <- readRDS(file = fileNameF4Haplo)


##### 1.2. Linkage Information data #####
fileNameRR <- paste0(linkInfoDir, "0.1.mapDistance_rr.csv")
linkInfoDf <- read.csv(file = fileNameRR, header = TRUE, row.names = 1)


##### 1.3. Marker effects data #####
### single-trait GP
fileNameMrkEff <- paste0(mrkEffDir, "0.2.F4_allPops_bindedMrkEffLst.rds")
mrkEffLst <- readRDS(file = fileNameMrkEff)

nTrait <- dim(mrkEffLst[[1]])[1]
traitNames <- dimnames(mrkEffLst[[1]])[[1]]


##### 1.4. Selected parents data #####
### selected by k-medoids
fileNameSelIndName <- paste0(crsCombDir, "0.3.2Pops_selectedIndName_kMedoids_nClust_20.rds")
selectedIndName <- readRDS(file = fileNameSelIndName)
fileNameCrsComb <- paste0(crsCombDir, "0.3.crossCombTable.csv")
selectedCrsCombNameDf <- read.csv(file = fileNameCrsComb, header = TRUE, row.names = 1)

nSelectInd <- length(selectedIndName)
nComb <- nrow(selectedCrsCombNameDf)




###### 2. Modify data ######
##### 2.1. Haplotype data of selected individuals #####
f4SelectedHaplo1 <- f4HaploLst[[1]][, selectedIndName]
f4SelectedHaplo2 <- f4HaploLst[[2]][, selectedIndName]


##### 2.2. Marker effect of selected individuals #####
### for selected by k-medoids
isS827 <- str_detect(string = selectedIndName, pattern = "S827")
isS840 <- str_detect(string = selectedIndName, pattern = "S840")

mrkEffNoVec <- rep(NA, nSelectInd)
names(mrkEffNoVec) <- selectedIndName

mrkEffNoVec[isS827] <- 1
mrkEffNoVec[isS840] <- 2




###### 3. Genotypic variance of progeny ######
##### 3.1. For Selected Individuals #####
genVarAllCombAllGenArray <- array(data = NA, 
                                  dim = c(nTrait, nTrait, nComb, length(nGenVec)), 
                                  dimnames = list(traitNames, 
                                                  traitNames, 
                                                  apply(X = selectedCrsCombNameDf[, 1:2], 
                                                        MARGIN = 1, 
                                                        FUN = paste0, collapse = " x "), 
                                                  paste0("G", nGenVec)))

for (genNo in nGenVec) {
  
  cat("------ Generation No.", genNo, " will be performed ------\n")
  if (genNo == 1) {
    k <- 1
  } else {
    k <- genNo + 1
  }
  ### LD parameter common to all the crosses
  DLst <- CalcD(linkInfoDf = linkInfoDf, k = k)
  
  ### calculate genetic variance
  genVarLst <- pbmcapply::pbmclapply(X = 1:nComb, FUN = function(combNo) {
    
    ### Name of Parents Now
    sireName <- selectedCrsCombNameDf$sireID[combNo]
    damName <- selectedCrsCombNameDf$damID[combNo]
    
    ### Haplotype of selected parents
    genoMatP1 <- rbind(Haplo1 = f4SelectedHaplo1[, sireName],
                       Haplo2 = f4SelectedHaplo2[, sireName])
    genoMatP2 <- rbind(Haplo1 = f4SelectedHaplo1[, damName],
                       Haplo2 = f4SelectedHaplo2[, damName])
    
    ### Marker effects of selected parents
    mrkEffNoP1 <- mrkEffNoVec[sireName]
    mrkEffNoP2 <- mrkEffNoVec[damName]
    mrkEffP1 <- t(mrkEffLst[[mrkEffNoP1]])
    mrkEffP2 <- t(mrkEffLst[[mrkEffNoP2]]) ### single-trait GP
    
    ### genotypic variance of one cross combination
    genVarOneCombArray <- array(data = NA, dim = c(nTrait, nTrait, 1), 
                                dimnames = list(traitNames, traitNames))
    
    genVarOneCombArray[, , 1] <- GenCovProgeny(ObjectGenoP1 = genoMatP1,
                                               ObjectGenoP2 = genoMatP2,
                                               DLst = DLst,
                                               MarkEffect1 = mrkEffP1,
                                               MarkEffect2 = mrkEffP2)
    
    return(genVarOneCombArray)
    
  }, mc.cores = nCoresComb)
  
  genVarCovMatAllCombArray <- do.call(what = abind::abind, args = genVarLst)
  
  genNameNow <- paste0("G", genNo)
  genVarAllCombAllGenArray[, , , genNameNow] <- genVarCovMatAllCombArray
  
}

### save genotypic variance
### single-trait GP
fileNameGenVar <- paste0(saveDir, scriptID,
                         "progGenVarAllGenerations_nCrs_",
                         nComb, ".rds")
saveRDS(object = genVarAllCombAllGenArray, file = fileNameGenVar)