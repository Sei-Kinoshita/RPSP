###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)


##### 0.2. Load Packages and Functions #####
source("scripts/function_scripts/1.0.function.simuProgeny.r")


##### 0.3. Define Parameters #####
scriptID <- "1.0."

nProgMate <- 1000
nSimu <- 30
nGenTot <- 7 ### total number of generation in simulation

nCores <- parallel::detectCores() - 3


##### 0.4. Define Directories #####
crsCombDir <- "data/crossComb/kMedoids_nClust_20/"
haploDir <- "data/haplotype/"
mrkEffDir <- "data/mrkEff/"
linkageDir <- "data/linkageInfo/"

saveDir <- paste0("result/", scriptID, "GeneticVariance_Simulation/")
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}




###### 1. Load Data ######
##### 1.1. Selected Cross Combination #####
fileNameSelIndName <- paste0(crsCombDir, "0.3.2Pops_selectedIndName_kMedoids_nClust_20.rds")
selIndName <- readRDS(file = fileNameSelIndName)
fileNameSelCrs <- paste0(crsCombDir, "0.3.crossCombTable.csv")
crsCombDf <- read.csv(file = fileNameSelCrs, header = TRUE, row.names = 1)

nComb <- nrow(crsCombDf)


##### 1.2. Haplotype Data #####
fileNameHaplo <- paste0(haploDir, "0.2.F4895_Haplotype_4ParentsArray.rds")
f4AllHaploArray <- readRDS(file = fileNameHaplo)
f4SelHaploArray <- f4AllHaploArray[, , selIndName, ]


##### 1.3. Marker Effects Data #####
fileNameMrkEff <- paste0(mrkEffDir, "0.2.F4_allTraits_bindedMrkEffLst.rds")
mrkEffAllTraitLst <- readRDS(file = fileNameMrkEff)

nTrait <- length(mrkEffAllTraitLst)
traitNames <- names(mrkEffAllTraitLst)


##### 1.4. Linkage Information #####
### marker name in each chromosome
fileNameMrkNameEachChr <- paste0(linkageDir, "0.2.F4_mrkName_eachChromLst.rds")
mrkNameEachChrLst <- readRDS(file = fileNameMrkNameEachChr)
### recombination rate in each chromosome
fileNameRREachChr <- paste0(linkageDir, "0.2.F4_rr_eachChromLst.rds")
rrEachChrLst <- readRDS(file = fileNameRREachChr)




###### 2. Simulate Genetic Variance of Progeny ######
genVarAllCrsLst <- list()

for (crsCombNo in 1:nComb) {
  
  # tictoc::tic()
  
  cat("------ Cross Combination No.", crsCombNo, "/", nComb, " will be performed ------\n")
  
  sireNameNow <- crsCombDf$sireID[crsCombNo]
  damNameNow <- crsCombDf$damID[crsCombNo]
  
  genVarAllSimuArray <- 
    do.call(what = abind::abind, 
            args = pbmcapply::pbmclapply(X = 1:nSimu, FUN = function(simuNo) {
              
              genVarOneSimuArray <- array(data = NA, dim = c(nTrait, nGenTot, 1), 
                                          dimnames = list(traitNames, 
                                                          paste0("G", 1:nGenTot)))
              genVarOneSimuArray[, , 1] <- 
                matingChain4Parents(haploArrayP1 = f4SelHaploArray[, , sireNameNow, ],
                                    haploArrayP2 = f4SelHaploArray[, , damNameNow, ],
                                    mrkEffLst = mrkEffAllTraitLst, 
                                    nProgMate = nProgMate, 
                                    nGenTot = nGenTot)
              
              return(genVarOneSimuArray)
            }, mc.cores = nCores))
  
  genVarAllSimuMeanMat <- apply(X = genVarAllSimuArray, MARGIN = c(1, 2), 
                                FUN = mean)
  
  genVarAllCrsLst <- c(genVarAllCrsLst, list(genVarAllSimuMeanMat))
  
  # tictoc::toc()
}

names(genVarAllCrsLst) <- apply(X = crsCombDf[, 1:2], MARGIN = 1, 
                                FUN = paste0, collapse = " x ")

fileNameProgGenVarAllCrs <- paste0(saveDir, scriptID, 
                                   "progGenVarAllGenerations_nProgMate_", nProgMate, 
                                   "_nCrs_", nComb, 
                                   "_nSimu_", nSimu, ".rds")
saveRDS(object = genVarAllCrsLst, file = fileNameProgGenVarAllCrs)