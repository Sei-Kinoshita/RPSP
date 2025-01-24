###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)


##### 0.2. Load Packages and Functions #####
require(stringr)
require(ggplot2)

source("scripts/function_scripts/6.1.3.function.calculateChangeOfGv.progeny.r")


##### 0.3. Define Parameters #####
dirID <- "6.1.4."
scriptID <- "6.1.4.1."

selTraitNames <- c("RosmarinicAcid", "Perillaldehyde")
nTrait <- length(selTraitNames)

nGenTot <- 15

nMethod <- 2 ### number of breeding methods

nSimu <- 50

gvResQuantile <- 0.99 ### extract 95% quantile of each generation

##### 0.4. Define Directories #####

dataDir <- "data/"
resDir <-"result/"

haploDir <- paste0(dataDir, "haplotype/")
mrkEffDir <- paste0(dataDir, "mrkEff/")
gvDir <- paste0(dataDir, "genotypicValue/")

saveDir <- paste0(resDir, dirID, "BreedProgram_Compare/GV/")
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Load Data ######
##### 1.1. Genotypic Value of Initial F4 Population #####
fileNameGvInit <- paste0(gvDir, "0.3.1.gv_allTrait.rds")
gvInitMat <- readRDS(file = fileNameGvInit)


##### 1.2. GV of each generation #####
### parameter set for selection generation
fileNameSelGenParSet <- paste0(saveDir, "6.0.3.parameterSet.csv")
selGenParSetDf <- read.csv(file = fileNameSelGenParSet, header = TRUE, row.names = 1)
nParSet <- nrow(selGenParSetDf)
nSelGenNo <- length(unique(selGenParSetDf$selGenNo))


##### 1.4. Marker effect #####
fileNameMrkEff <- paste0(mrkEffDir, "0.2.F4_allTraits_bindedMrkEffLst.rds")
mrkEffLst <- readRDS(fileNameMrkEff)[selTraitNames]




###### 3. Change of GV through All Generations (index) ######
##### 3.0. Calculate the change of GV through all generations #####
#### 3.0.1. method 2 (highBV, Cross) ####
calcGvResM2AllParLst <- list()
for (parSetNo in 1:nSelGenNo) {
  
  selGenNoFstCrs <- selGenParSetDf$selGenNoFstCrs[parSetNo]
  selGenNo <- selGenParSetDf$selGenNo[parSetNo]
  
  fileNameGvAllGenM2 <- paste0(saveDir, scriptID, 
                               "gvAllSimuAllGenLst_highBv_selGenNoFstCrs_", selGenNoFstCrs, 
                               "_selGenNo_", 
                               selGenNo, 
                               "_nSimu_", nSimu, "_nGen_15.rds")
  m2GvAllGenLst <- readRDS(file = fileNameGvAllGenM2)
  
  ### calculate mean through all simulations
  calcGvResM2Lst <- CalcGvChange(gvLst = m2GvAllGenLst,  
                                 gvResQuantile = gvResQuantile, 
                                 nTrait = nTrait, 
                                 nSimu = nSimu, 
                                 selTraitNames = selTraitNames,
                                 gvInitMat = gvInitMat)
  calcGvResM2AllParLst <- c(calcGvResM2AllParLst, list(calcGvResM2Lst))
}
names(calcGvResM2AllParLst) <- paste0("selGenNo_", 1:nSelGenNo)




#### 3.0.2. method 3 (progeny base, cross) ####
### genoArray to gv each generation
calcGvResM3AllParLst <- list()
for (parSetNo in 1:nSelGenNo) {
  
  selGenNoFstCrs <- selGenParSetDf$selGenNoFstCrs[parSetNo]
  selGenNo <- selGenParSetDf$selGenNo[parSetNo]
  
  fileNameGvAllSimuAllGenLst <- paste0(saveDir, 
                                       "6.1.3.BreedProgram_Compare/GV/6.1.3.gvAllSimuAllGenLst_progVarBase_", 
                                       "selGenNoFstCrs_", selGenNoFstCrs, 
                                       "_selGenNo_", selGenNo, 
                                       "_nSimu_", nSimu, 
                                       "_nGen_", nGenTot, ".rds")
  
  
  gvAllSimuBindParNowLst <- readRDS(file = fileNameGvAllSimuAllGenLst)
  
  
  ### calculate mean through all simulations
  calcGvResM3Lst <- CalcGvChange(gvLst = gvAllSimuBindParNowLst,  
                                 gvResQuantile = gvResQuantile, 
                                 nTrait = nTrait, 
                                 nSimu = nSimu, 
                                 selTraitNames = selTraitNames,
                                 gvInitMat = gvInitMat)
  
  calcGvResM3AllParLst <- c(calcGvResM3AllParLst, list(calcGvResM3Lst))
}
names(calcGvResM3AllParLst) <- paste0("ParSetNo_", 1:4)



##### 3.1. Plot Change of GV #####
indexName <- dimnames(calcGvResM2AllParLst$selGenNo_1$oneTrait$All)[[1]]
nIndex <- length(indexName)
nGenChange <- dim(calcGvResM2AllParLst$selGenNo_1$oneTrait$All)[3]

plotColMethod <- c("#808000", "#6a1e48")
plotLtyMethod <- c("73", "solid")
plotLtyPar <- c(1, 5, 2, 3) ### line type for different parameter

#### 3.1.1. For One Trait ####
fileNameGvChangeCompare <- paste0(saveDir, scriptID, 
                                  "ChangeOfGV_CompareEachMethod_nSimu_", 
                                  nSimu, "_thesis", ".pdf")
pdf(file = fileNameGvChangeCompare, width = 15, height = 6)
for (resNo in 1:nIndex) {
  for (traitNo in 1:length(selTraitNames)) {
    
    traitNamesNow <- selTraitNames[traitNo]
    
    ### high GV, individual base, all combination
    plotIndexM2Lst <- 
      lapply(X = 1:nSelGenNo, FUN = function(parNo) {
        
        meanNow <- calcGvResM2AllParLst[[parNo]]$oneTrait$mean[resNo, traitNamesNow, ]
        sdNow <- calcGvResM2AllParLst[[parNo]]$oneTrait$sd[resNo, traitNamesNow, ]
        erUpperNow <- meanNow + sdNow
        erLowerNow <- meanNow - sdNow
        
        return(rbind(meanNow, sdNow, erUpperNow, erLowerNow))
      })
    
    ### progeny high GV, cross pair base
    plotIndexM3Lst <- 
      lapply(X = 1:nSelGenNo, FUN = function(parNo) {
        
        meanNow <- calcGvResM3AllParLst[[parNo]]$oneTrait$mean[resNo, traitNamesNow, ]
        sdNow <- calcGvResM3AllParLst[[parNo]]$oneTrait$sd[resNo, traitNamesNow, ]
        erUpperNow <- meanNow + sdNow
        erLowerNow <- meanNow - sdNow
        
        return(rbind(meanNow, sdNow, erUpperNow, erLowerNow))
      })
    
    ### plot
    plotRange <- range(c(unlist(lapply(X = plotIndexM2Lst, FUN = function(x) {x[3:4, ]})), 
                         c(unlist(lapply(X = plotIndexM3Lst, FUN = function(x) {x[3:4, ]})))))
    mainTitle <- paste0("Change of GV ", indexName[resNo], " ", traitNamesNow)
    
    ### nSelGenCrs in different plot
    par(mfrow = c(1, 4), mar = c(5, 2.3, 4, 0.5))
    for (parNo in 1:nSelGenNo) {
      
      plot(NULL, type = "l", ylim = plotRange, cex.axis = 1.7, 
           xlim = c(0, (nGenChange - 1)), xlab = "Generation", ylab = " ", 
           main = mainTitle)
      
      plotIndexM2Now <- plotIndexM2Lst[[parNo]]
      plotIndexM3Now <- plotIndexM3Lst[[parNo]]
      
      ### dark grey lines for other nSelGenCrs
      if (parNo != 1) {
        for (x in 1:(parNo - 1)) {
          ### high GV, individual base, all combination
          lines(x = 0:(nGenChange - 1), y = plotIndexM2Lst[[x]][1, ],
                type = "l", lwd = 2, lty = plotLtyMethod[1], 
                col = "darkgrey")
          points(x = 0:(nGenChange - 1), y = plotIndexM2Lst[[x]][1, ], 
                 col = "darkgrey", pch = 19)
          arrows(x0 = 0:(nGenChange - 1), y0 = plotIndexM2Lst[[x]][3, ],
                 x1 = 0:(nGenChange - 1), y1 = plotIndexM2Lst[[x]][4, ],
                 code = 3, angle = 90, length = 0.01, col = "darkgrey")
          ### progeny high GV, cross pair base
          lines(x = 0:(nGenChange - 1), y = plotIndexM3Lst[[x]][1, ],
                type = "l", lwd = 2, lty = plotLtyMethod[2], 
                col = "darkgrey")
          points(x = 0:(nGenChange - 1), y = plotIndexM3Lst[[x]][1, ], 
                 col = "darkgrey", pch = 19)
          arrows(x0 = 0:(nGenChange - 1), y0 = plotIndexM3Lst[[x]][3, ],
                 x1 = 0:(nGenChange - 1), y1 = plotIndexM3Lst[[x]][4, ],
                 code = 3, angle = 90, length = 0.01, col = "darkgrey")
        }
      }
      
      ### color lines for nSelGenCrsNow
      ### high GV, individual base, all combination
      lines(x = 0:(nGenChange - 1), y = plotIndexM2Now[1, ],
            type = "l", lwd = 3, lty = plotLtyMethod[1], 
            col = plotColMethod[1])
      points(x = 0:(nGenChange - 1), y = plotIndexM2Now[1, ], 
             col = plotColMethod[1], pch = 19)
      arrows(x0 = 0:(nGenChange - 1), y0 = plotIndexM2Now[3, ],
             x1 = 0:(nGenChange - 1), y1 = plotIndexM2Now[4, ],
             code = 3, angle = 90, length = 0.01)
      ### progeny high GV, cross pair base
      lines(x = 0:(nGenChange - 1), y = plotIndexM3Now[1, ],
            type = "l", lwd = 3, lty = plotLtyMethod[2], 
            col = plotColMethod[2])
      points(x = 0:(nGenChange - 1), y = plotIndexM3Now[1, ], 
             col = plotColMethod[2], pch = 19)
      arrows(x0 = 0:(nGenChange - 1), y0 = plotIndexM3Now[3, ],
             x1 = 0:(nGenChange - 1), y1 = plotIndexM3Now[4, ],
             code = 3, angle = 90, length = 0.01)
    }
    par(mfrow = c(1, 1))
  }
}
dev.off()


#### 3.1.2. For Selection Index ####
fileNameGvChangeCompSelIndex <- paste0(saveDir, scriptID, 
                                       "ChangeOfGV_SelIndex_CompareEachMethod_nSimu_", 
                                       nSimu, "_thesis", ".pdf")
pdf(file = fileNameGvChangeCompSelIndex, width = 15, height = 6)
for (resNo in 1:nIndex) {
  
  ### high GV, individual base, all combination
  plotIndexM2Lst <- 
    lapply(X = 1:nSelGenNo, FUN = function(parNo) {
      
      meanNow <- calcGvResM2AllParLst[[parNo]]$selIndex$mean[resNo, ]
      sdNow <- calcGvResM2AllParLst[[parNo]]$selIndex$sd[resNo, ]
      erUpperNow <- meanNow + sdNow
      erLowerNow <- meanNow - sdNow
      
      return(rbind(meanNow, sdNow, erUpperNow, erLowerNow))
    })
  
  ### progeny high GV, cross pair base
  plotIndexM3Lst <- 
    lapply(X = 1:nSelGenNo, FUN = function(parNo) {
      
      meanNow <- calcGvResM3AllParLst[[parNo]]$selIndex$mean[resNo, ]
      sdNow <- calcGvResM3AllParLst[[parNo]]$selIndex$sd[resNo, ]
      erUpperNow <- meanNow + sdNow
      erLowerNow <- meanNow - sdNow
      
      return(rbind(meanNow, sdNow, erUpperNow, erLowerNow))
    })
  
  ### plot
  plotRange <- range(c(unlist(lapply(X = plotIndexM2Lst, FUN = function(x) {x[3:4, ]})), 
                       c(unlist(lapply(X = plotIndexM3Lst, FUN = function(x) {x[3:4, ]})))))
  mainTitle <- paste0("Change of GV ", indexName[resNo], " Selection Index")
  
  ### nSelGenCrs in different plot
  par(mfrow = c(1, 4), mar = c(5, 2.3, 4, 0.5))
  for (parNo in 1:nSelGenNo) {
    
    plot(NULL, type = "l", ylim = plotRange, cex.axis = 1.7, 
         xlim = c(0, (nGenChange - 1)), xlab = "Generation", ylab = " ", 
         main = mainTitle)
    
    plotIndexM2Now <- plotIndexM2Lst[[parNo]]
    plotIndexM3Now <- plotIndexM3Lst[[parNo]]
    
    ### dark grey lines for other nSelGenCrs
    if (parNo != 1) {
      for (x in 1:(parNo - 1)) {
        ### high GV, individual base, all combination
        lines(x = 0:(nGenChange - 1), y = plotIndexM2Lst[[x]][1, ],
              type = "l", lwd = 2, lty = plotLtyMethod[1], 
              col = "darkgrey")
        points(x = 0:(nGenChange - 1), y = plotIndexM2Lst[[x]][1, ], 
               col = "darkgrey", pch = 19)
        arrows(x0 = 0:(nGenChange - 1), y0 = plotIndexM2Lst[[x]][3, ],
               x1 = 0:(nGenChange - 1), y1 = plotIndexM2Lst[[x]][4, ],
               code = 3, angle = 90, length = 0.01, col = "darkgrey")
        ### progeny high GV, cross pair base
        lines(x = 0:(nGenChange - 1), y = plotIndexM3Lst[[x]][1, ],
              type = "l", lwd = 2, lty = plotLtyMethod[2], 
              col = "darkgrey")
        points(x = 0:(nGenChange - 1), y = plotIndexM3Lst[[x]][1, ], 
               col = "darkgrey", pch = 19)
        arrows(x0 = 0:(nGenChange - 1), y0 = plotIndexM3Lst[[x]][3, ],
               x1 = 0:(nGenChange - 1), y1 = plotIndexM3Lst[[x]][4, ],
               code = 3, angle = 90, length = 0.01, col = "darkgrey")
      }
    }
    
    ### color lines for nSelGenCrsNow
    ### high GV, individual base, all combination
    lines(x = 0:(nGenChange - 1), y = plotIndexM2Now[1, ],
          type = "l", lwd = 3, lty = plotLtyMethod[1], 
          col = plotColMethod[1])
    points(x = 0:(nGenChange - 1), y = plotIndexM2Now[1, ], 
           col = plotColMethod[1], pch = 19)
    arrows(x0 = 0:(nGenChange - 1), y0 = plotIndexM2Now[3, ],
           x1 = 0:(nGenChange - 1), y1 = plotIndexM2Now[4, ],
           code = 3, angle = 90, length = 0.01)
    ### progeny high GV, cross pair base
    lines(x = 0:(nGenChange - 1), y = plotIndexM3Now[1, ],
          type = "l", lwd = 3, lty = plotLtyMethod[2], 
          col = plotColMethod[2])
    points(x = 0:(nGenChange - 1), y = plotIndexM3Now[1, ], 
           col = plotColMethod[2], pch = 19)
    arrows(x0 = 0:(nGenChange - 1), y0 = plotIndexM3Now[3, ],
           x1 = 0:(nGenChange - 1), y1 = plotIndexM3Now[4, ],
           code = 3, angle = 90, length = 0.01)
  }
  par(mfrow = c(1, 1))
}
dev.off()
