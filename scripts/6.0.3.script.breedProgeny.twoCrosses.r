###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)


##### 0.2. Load Packages and Functions #####
source("scripts/function_scripts/6.0.3.function.breedProgramTwoCrosses.r")
source("scripts/function_scripts/2.0.function.clacProgeny.genVar.r")


##### 0.3. Define Parameters #####
scriptID <- "6.0.3."

### common to all three methods
nGenTot <- 15 ### total number of generations
selGenNoVec <- c(2, 3, 4, 5) ### at which generation selection will be done

### individual selection based on BV, all-combination
nSelInd <- 9 ### how many individuals to select at the second cross generation
nProgG0M2 <- 25 
nProgGSelM2 <- 25

### cross pair selection based on progeny
nSelParent <- 450 ### how many parents to narrow down at selGenNo
nSelCrs <- 10
nProgG0M3 <- 90
nProgGSelM3 <- 90

selTraitNames <- c("RosmarinicAcid", "Perillaldehyde")

nSimu <- 50
nCores <- 1


##### 0.4. Define Directories #####
dataDir <- "data/"
resDir <- "result/"

haploDir <- paste0(dataDir, "haplotype/")
mrkEffDir <- paste0(dataDir, "mrkEff/")
linkageDir <- paste0(dataDir, "linkageInfo/")

gvDir <- paste0(dataDir, "genotypicValue/")

crsCombDir <- paste0(dataDir, "crossComb/")
crsPairBaseCrsCombDir <- paste0(crsCombDir, "F4_crsPairBase/")

### saving directory
scrDir <- paste0(resDir, scriptID, "Selection_result_Progeny_TwoCross/")
if (!dir.exists(scrDir)) {
  dir.create(scrDir, recursive = TRUE)
}
progVarBaseDir <- paste0(scrDir, "progVarBase/")
if (!dir.exists(progVarBaseDir)) {
  dir.create(progVarBaseDir, recursive = TRUE)
}




###### 1. Load Data ######
##### 1.1. Haplotype Data #####
fileNameHaplo <- paste0(haploDir, "0.2.F4895_Haplotype_4ParentsArray.rds")
f4AllHaploArray <- readRDS(file = fileNameHaplo)


##### 1.2. Marker Effects Data #####
fileNameMrkEff <- paste0(mrkEffDir, "0.2.F4_allTraits_bindedMrkEffLst.rds")
mrkEffAllTraitLst <- readRDS(file = fileNameMrkEff)
mrkEffLst <- mrkEffAllTraitLst[selTraitNames]


##### 1.3. GV mean and Sd of the Initial Population #####
fileNameInitPopInfo <- paste0(gvDir, "0.3.1.gv_initialPop_meanSd.rds")
initPopInfoLst <- readRDS(file = fileNameInitPopInfo)


##### 1.4. Selected Individual Data ######
### progeny high GV, cross pair base
fileNameProgBaseCrs <- paste0(crsPairBaseCrsCombDir, 
                              "0.3.1.crsCombDf_progVarBase_nCrs_10.rds")
progBaseCrsCombDf <- readRDS(file = fileNameProgBaseCrs)


##### 1.5. Linkage Information #####
### marker name in each chromosome
fileNameMrkNameEachChr <- paste0(linkageDir, "0.2.F4_mrkName_eachChromLst.rds")
mrkNameEachChrLst <- readRDS(file = fileNameMrkNameEachChr)
### recombination rate in each chromosome
fileNameRREachChr <- paste0(linkageDir, "0.2.F4_rr_eachChromLst.rds")
rrEachChrLst <- readRDS(file = fileNameRREachChr)
### map distance for calculating LD parameter
fileNameRR <- paste0(linkageDir, "0.1.mapDistance_rr.csv")
linkInfoDf <- read.csv(file = fileNameRR, header = TRUE, row.names = 1)




###### 2. Preparation for Generate Progeny ######
##### 2.1. create parameter set for breeding simulation #####
selGenParDf <- data.frame(selGenNo = rep(selGenNoVec, 2), 
                          selGenNoFstCrs = c(selGenNoVec, 
                                             rep(nGenTot, length(selGenNoVec))))
nParSet <- nrow(selGenParDf)

fileNameSelGenPar <- paste0(scrDir, scriptID, "parameterSet.csv")
if (!file.exists(fileNameSelGenPar)) {
  write.csv(x = selGenParDf, file = fileNameSelGenPar)
}




###### 3. Generate Progeny & Calculate GV ######
##### 3.1. Progeny High GV, Cross Pair Base #####
for (parSetNo in 1:nParSet) {

  cat("------ Parameter Set No.", parSetNo, "/", nParSet, " will be performed ------\n")
  
  ### generation of 2nd cross
  selGenNo <- selGenParDf$selGenNo[parSetNo]
  ### first cross pairs selected based on which generation
  selGenNoFstCrs <- selGenParDf$selGenNoFstCrs[parSetNo]
  
  ### LD parameter (DLst) for method 3
  DLst <- CalcD(linkInfoDf = linkInfoDf, k = (nGenTot - selGenNo))
  ### haplotype data
  progBaseGNowCrsCombDf <- progBaseCrsCombDf[[paste0("G", selGenNoFstCrs)]]
  progBaseCrsIndName <- unique(c(progBaseGNowCrsCombDf[, 1], progBaseGNowCrsCombDf[, 2]))
  progBaseCrsSelHaploArray <- f4AllHaploArray[, , progBaseCrsIndName, ]
  
  for (simuNo in 1:nSimu) {
    
    ### simulation
    dirNameAllResProgBaseGv <- paste0(progVarBaseDir,  
                                      "nGenTot_", nGenTot, 
                                      "_selGenNoFstCrs_", selGenNoFstCrs, 
                                      "_selGenNo_", selGenNo, 
                                      "_nSelParent_", nSelParent, 
                                      "_nSelCrs_", nSelCrs, 
                                      "_simu_", simuNo, "/")
    
    if (!dir.exists(dirNameAllResProgBaseGv)) {
      
      dir.create(dirNameAllResProgBaseGv, recursive = TRUE)
      
      if ((simuNo %% 10) == 1) {
        cat("------ Simulation No.", simuNo, "/", nSimu, " will be performed ------\n")
      }
      
      gvOneSimuLst <- BreedSimuUcTwoCrs(crsTbl = progBaseGNowCrsCombDf[, 1:2],
                                        haploArray = progBaseCrsSelHaploArray,
                                        mrkEffLst = mrkEffLst,
                                        DLst = DLst,
                                        initPopInfoLst = initPopInfoLst,
                                        nGenTot = nGenTot,
                                        selGenNo = selGenNo,
                                        nSelParent = nSelParent,
                                        nSelCrs = nSelCrs,
                                        nProgG0 = nProgG0M3,
                                        nProgGSel = nProgGSelM3, 
                                        nCores = nCores, 
                                        saveDir = dirNameAllResProgBaseGv)
    }
    
    gc(reset = TRUE);gc(reset = TRUE)
  }
}
