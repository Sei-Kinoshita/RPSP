###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load packages #####
require(stringr)


##### 0.3. Define parameters #####
scriptID <- "3.0."

nGenVec <- 1:7 ### 4: cross hetero lines --> selfing 3 generation

selMethod <- "kMedoids" ### selection method of parents


##### 0.4. Define directories #####
progVarSimuDir <- "result/1.0.GeneticVariance_Simulation/"
progVarCalcDir <- "result/2.0.GeneticVariance_Calculation/"
crsCombDir <- "data/crossComb/kMedoids_nClust_20/"

saveDir <- paste0("result/", scriptID, "GeneticVarinace_SimuVsCalc/")
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}




###### 1. Load Data ######
##### 1.1. Simulation Data #####
fileNameSimuVar <- paste0(progVarSimuDir, "1.0.progGenVarAllGenerations_nProgMate_1000_nCrs_820_nSimu_30.rds")
progVarSimuLst <- readRDS(file = fileNameSimuVar)

##### 1.2. Calculation Data #####
fileNameCalcVar <- paste0(progVarCalcDir, "2.0.progGenVarAllGenerations_nCrs_820.rds")
progVarCalcArray <- readRDS(file = fileNameCalcVar)

traitNames <- dimnames(progVarCalcArray)[[1]]
nTrait <- length(traitNames)


##### 1.3. Cross Combination Data #####
fileNameCrsComb <- paste0(crsCombDir, "0.3.crossCombTable.csv")
selectedCrsCombNameDf <- read.csv(file = fileNameCrsComb, header = TRUE, row.names = 1)




###### 2. Check simulation and calculation ######
fileNameSimuCalcMatchPlot <- paste0(saveDir, scriptID, "progGenVar_simuVsCalc_matchPlot.pdf")
pdf(file = fileNameSimuCalcMatchPlot, width = 7, height = 7)
for (traitNo in 1:nTrait) {
  
  traitNamesNow <- traitNames[traitNo]
  
  for (genNo in nGenVec) {
    
    genNameNow <- paste0("G", genNo)
    
    ### genetic variance by simulation
    genVarSimuNow <- unlist(lapply(X = progVarSimuLst, 
                                   FUN = function(x) {
                                     x[traitNo, genNameNow]
                                   }))
    
    ### genetic variance by calculation
    genVarCalcNow <- apply(X = progVarCalcArray[, , , genNameNow], 
                           MARGIN = 3, 
                           FUN = function(x) {
                             x[traitNo, traitNo]
                           })
    
    ### plot simulation v.s. calculation
    mainTitle <- paste0("Simulation v.s. Calculation, ", traitNamesNow, ", ", genNameNow)
    plot(x = genVarSimuNow, y = genVarCalcNow, 
         col = as.factor(selectedCrsCombNameDf$crossPattern), 
         xlab = "Simulation", ylab = "Calculation", 
         main = mainTitle) ## Across:1, Selfing:2, Within:3
    # legend("topleft", legend = c("Across", "Selfing", "within"), col = 1:3, 
    #        pch = 1)
    abline(a = 0, b = 1, col = "red")
    
    corNow <- format(cor(x = genVarSimuNow, y = genVarCalcNow), digits = 4)
    legend("bottomright", legend = paste0("R2 = ", corNow))
    
  }
}
dev.off()