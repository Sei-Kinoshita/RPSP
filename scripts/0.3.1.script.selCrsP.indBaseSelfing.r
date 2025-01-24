###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)


##### 0.2. Load Packages #####
require(stringr)
require(pbmcapply)


##### 0.3. Define Parameters #####
scriptID <- "0.3.1."

nPops <- 2

nSelInd <- 9 ### inidvidual selection based on BV, only selfing
nSelInd27 <- 5 ### select parents from both populations
nSelInd40 <- 4 ### select parents from both populations

nSelIndForCrs <- 9 ### cross pair selection based on BV, all-combination
nSelCrsForProgVar <- 10 ### cross pair selection based on progeny
selTraitNames <- c("RosmarinicAcid", "Perillaldehyde")

genVec <- 1:15 ### generation for calculated progeny variance

k <- qnorm(p = 0.975) ### top 5%

nCores <- parallel::detectCores() - 3


##### 0.4. Define Directories #####
dataDir <- "data/"
resDir <- "result/"

haploDir <- paste0(dataDir, "haplotype/")
mrkEffDir <- paste0(dataDir, "mrkEff/")
genVarDir <- paste0(resDir, "2.0.GeneticVariance_Calculation/")
crsCombDir <- paste0(dataDir, "crossComb/")

gvDir <- paste0(dataDir, "genotypicValue/")
if (!dir.exists(gvDir)) {
  dir.create(gvDir, recursive = TRUE)
}

indBaseCrsCombDir <- paste0(crsCombDir, "F4_indBase/")
if (!dir.exists(indBaseCrsCombDir)) {
  dir.create(indBaseCrsCombDir, recursive = TRUE)
}
crsPairBaseCrsCombDir <- paste0(crsCombDir, "F4_crsPairBase/")
if (!dir.exists(crsPairBaseCrsCombDir)) {
  dir.create(crsPairBaseCrsCombDir, recursive = TRUE)
}




###### 1. Load Data ######
##### 1.1. Haplotype Data #####
fileNameF4AllHaplo <- paste0(haploDir, "0.1.F4895_Haplotype.rds")
fileName27Haplo <- paste0(haploDir, "0.1.F4_S827_Haplotype.rds")
fileName40Haplo <- paste0(haploDir, "0.1.F4_S840_Haplotype.rds")
f4AllHaploLst <- readRDS(file = fileNameF4AllHaplo)
s827HaploLst <- readRDS(file = fileName27Haplo)
s840HaploLst <- readRDS(file = fileName40Haplo)

f4AllGenoMat <- f4AllHaploLst[[1]] + f4AllHaploLst[[2]]
s827GenoMat <- s827HaploLst[[1]] + s827HaploLst[[2]]
s840GenoMat <- s840HaploLst[[1]] + s840HaploLst[[2]]

f4AllGenoMatLst <- list(S827 = s827GenoMat, S840 = s840GenoMat)


##### 1.2. Marker Effects Data #####
fileNameMrkEff <- paste0(mrkEffDir, "0.2.F4_allPops_bindedMrkEffLst.rds")
mrkEffLst <- readRDS(file = fileNameMrkEff)

nTrait <- dim(mrkEffLst[[1]])[1]
traitNames <- dimnames(mrkEffLst[[1]])[[1]]


##### 1.3. Genetic Variance of Progeny #####
fileNameProgGenVar <- paste0(genVarDir, "2.0.progGenVarAllGenerations_F4AllCross_G1-G15.rds")
genVarArray <- readRDS(file = fileNameProgGenVar)




###### 2. Select Individual with high BV ######
### calculate BV
gvMat <- do.call(what = rbind, 
                 args = lapply(X = 1:nPops, FUN = function(popNo) {
                   t(f4AllGenoMatLst[[popNo]]) %*% t(mrkEffLst[[popNo]])
                 }))

gvScaledMat <- apply(X = gvMat, MARGIN = 2, FUN = scale) ### scaling
rownames(gvScaledMat) <- rownames(gvMat)

gvSelTraitMat <- gvScaledMat[, selTraitNames] ### GV for selected traits

### add up BV of selected traits
gvSumVec <- apply(X = gvSelTraitMat, MARGIN = 1, FUN = sum)

### order sum of GV of selected traits
gvSumOrderVec <- sort(gvSumVec, decreasing = TRUE)
gvSumOrderVec27 <- gvSumOrderVec[str_detect(string = names(gvSumOrderVec), pattern = "S827")]
gvSumOrderVec40 <- gvSumOrderVec[str_detect(string = names(gvSumOrderVec), pattern = "S840")]

### select individuals with top BV
selIndNameHighBv <- names(gvSumOrderVec[1:nSelInd])
selIndNameHighBv2Pop <- c(names(gvSumOrderVec27[1:nSelInd27]), 
                          names(gvSumOrderVec40[1:nSelInd40]))



###### 3. Cross Table of each method ######
##### 3.1. individual-base, all combinations #####
### selection based on 2 population integrated
crsTable <- expand.grid(rep(list(1:nSelIndForCrs), 2))
crsTable <- crsTable[(crsTable[, 1] < crsTable[, 2]), ]

### selection from both populations
selIndHighBv2PopCrsTable <- do.call(what = rbind, 
                                    args = apply(X = crsTable, MARGIN = 1, FUN = function(x) {
                                      
                                      sireName <- selIndNameHighBv2Pop[as.numeric(x[1])]
                                      damName <- selIndNameHighBv2Pop[as.numeric(x[2])]
                                      
                                      sirePop <- str_sub(string = sireName, start = 1, end = 4)
                                      damPop <- str_sub(string = damName, start = 1, end = 4)
                                      if (sireName == damName) {
                                        crsPat <- "Selfing"
                                      } else if (sirePop == damPop) {
                                        crsPat <- "Within"
                                      } else {
                                        crsPat <- "Across"
                                      }
                                      c(sireName, damName, crsPat)
                                    }, simplify = FALSE))
selIndHighBv2PopCrsTable <- as.data.frame(selIndHighBv2PopCrsTable)
colnames(selIndHighBv2PopCrsTable) <- c("sireID", "damID", "crossPattern")


##### 3.2. cross pairs-base #####
### mean & variance of GV of parent population
gvMeanF4AllTraitVec <- apply(X = gvMat, MARGIN = 2, FUN = mean)
gvVarF4AllTraitVec <- apply(X = gvMat, MARGIN = 2, FUN = var)
### name of each cross pair
allCrsNameMat <- data.frame(do.call(what = rbind, 
                                    args = str_split(string = dimnames(genVarArray)[[3]], 
                                                     pattern = " x ")))
nCrsAll <- dim(genVarArray)[3]

### calculate GV of progeny
progGvArray <- do.call(what = abind::abind, 
                       args = pbmcapply::pbmclapply(X = 1:nCrsAll, 
                                                    FUN = function(crsNo) {
                                                      
                                                      progGvOneCrsArray <- array(data = NA, 
                                                                                 dim = c(nTrait, length(genVec), 1), 
                                                                                 dimnames = dimnames(genVarArray)[c(2, 4)])
                                                      
                                                      ### name of parent
                                                      sireNameNow <- allCrsNameMat[crsNo, 1]
                                                      damNameNow <- allCrsNameMat[crsNo, 2]
                                                      ### GV of parent
                                                      sireGvVec <- gvMat[sireNameNow, ]
                                                      damGvVec <- gvMat[damNameNow, ]
                                                      ### GV mean of progeny (corrected by mean of parent population)
                                                      # progGvMeanVec <- (sireGvVec + damGvVec) / 2 - gvMeanF4AllTraitVec
                                                      progGvMeanVec <- 
                                                        ((sireGvVec + damGvVec) / 2 - gvMeanF4AllTraitVec) / sqrt(gvVarF4AllTraitVec)
                                                      
                                                      ### GV sd of progeny (corrected by sd of parent population)
                                                      progGvSdMat <- matrix(data = NA, nrow = length(genVec), ncol = nTrait, 
                                                                            dimnames = c(dimnames(genVarArray)[4], list(traitNames)))
                                                      for (genNo in 1:length(genVec)) {
                                                        progGvSdNow <- 
                                                          sqrt(diag(genVarArray[, , crsNo, genNo])) / sqrt(gvVarF4AllTraitVec)
                                                        progGvSdMat[genNo, ] <- progGvSdNow
                                                      }
                                                      
                                                      ### (mean + kappa * sigma) of progeny
                                                      progGvOneCrsArray[, , 1] <- do.call(what = cbind, 
                                                                                          args = apply(X = progGvSdMat, MARGIN = 1, FUN = function(x) {
                                                                                            progGvMeanVec + k * x
                                                                                          }, simplify = FALSE))
                                                      
                                                      return(progGvOneCrsArray)
                                                    }, mc.cores = nCores))
progGvArray <- aperm(a = progGvArray, perm = c(3, 1, 2))
dimnames(progGvArray)[[1]] <- dimnames(genVarArray)[[3]]

### sum up selected traits (ros & peril)
progGvSumAllGenMat <- 
  progGvArray[, selTraitNames[1], ] + progGvArray[, selTraitNames[2], ]
progGvSumAllGenOrderedLst <- apply(X = progGvSumAllGenMat, 
                                   MARGIN = 2, 
                                   FUN = sort, decreasing = TRUE, 
                                   simplify = FALSE) ### order by each generation

### name of selected cross pair
selCrsPairNameForProgVarLst <- lapply(X = progGvSumAllGenOrderedLst, 
                                      FUN = function(x) {
                                        names(x)[1:nSelCrsForProgVar]
                                      })

selCrsTableForProgVarLst <- 
  lapply(X = selCrsPairNameForProgVarLst, 
         FUN = function(eachGenVec) {
           
           crsPairMat <- do.call(what = rbind, 
                                 args = str_split(string = eachGenVec, pattern = " x "))
           crsPatVec <- apply(X = crsPairMat, MARGIN = 1, FUN = function(x) {
             sirePop <- str_sub(string = x[1], start = 1, end = 4)
             damPop <- str_sub(string = x[2], start = 1, end = 4)
             if (x[1] == x[2]) {
               crsPat <- "Selfing"
             } else if (sirePop == damPop) {
               crsPat <- "Within"
             } else {
               crsPat <- "Across"
             }
             return(crsPat)
           })
           crsTableRes <- data.frame(crsPairMat, crsPatVec)
           colnames(crsTableRes) <- c("sireID", "damID", "crossPattern")
           
           return(crsTableRes)
         })




###### 3. Save ######
### GV of all trait
fileNameGv <- paste0(gvDir, scriptID, "gv_allTrait.rds")
saveRDS(object = gvMat, file = fileNameGv)

### Mean and sd of GV of the initial population
initPopInfoLst <- list(Mean = gvMeanF4AllTraitVec, Sd = sqrt(gvVarF4AllTraitVec))
fileNameInitPopInfo <- paste0(gvDir, scriptID, "gv_initialPop_meanSd.rds")
saveRDS(object = initPopInfoLst, file = fileNameInitPopInfo)

### cross table for individual-base, all combination
### selection from both populations
fileNameCrsTableHighBvCrs2Pop <- paste0(indBaseCrsCombDir, scriptID, 
                                    "crsCombDf_highBv_allCombination_2Pop_nCrs_", 
                                    nrow(selIndHighBv2PopCrsTable), ".rds")
saveRDS(object = selIndHighBv2PopCrsTable, file = fileNameCrsTableHighBvCrs2Pop)

### cross table for cross pair-base
fileNameCrsTableForProgVar <- paste0(crsPairBaseCrsCombDir, scriptID, 
                                     "crsCombDf_progVarBase_nCrs_", nSelCrsForProgVar, ".rds")
saveRDS(object = selCrsTableForProgVarLst, file = fileNameCrsTableForProgVar)

