###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)


##### 0.2. Load Packages #####
require(cluster)


##### 0.3. Define Directories #####
haploDir <- "data/haplotype/"

scrDir <- "data/crossComb/"
if (!dir.exists(scrDir)) {
  dir.create(scrDir, recursive = TRUE)
}


##### 0.4. Define Parameters #####
scriptID <- "0.3."

nClust <- 20
selMethod <- "kMedoids"




###### 1. Load Data ######
##### 1.1. Haplotype Data #####
fileName27Haplo <- paste0(haploDir, "0.1.F4_S827_Haplotype.rds")
fileName40Haplo <- paste0(haploDir, "0.1.F4_S840_Haplotype.rds")
s827HaploLst <- readRDS(file = fileName27Haplo)
s840HaploLst <- readRDS(file = fileName40Haplo)

s827GenoMat <- s827HaploLst[[1]] + s827HaploLst[[2]]
s840GenoMat <- s840HaploLst[[1]] + s840HaploLst[[2]]




###### 2. Compute k-medoids ######
### distance matrix
genoDist27 <- dist(x = t(s827GenoMat), method = "euclidean")
genoDist40 <- dist(x = t(s840GenoMat), method = "euclidean")

### k-medoids
kMedois27 <- pam(x = genoDist27, k = nClust)
kMedois40 <- pam(x = genoDist40, k = nClust)

### select representative lines
selIndName27 <- kMedois27$medoids
selIndName40 <- kMedois40$medoids
selIndName <- c(selIndName27, selIndName40)




###### 3. Create Cross Combination Table ######
nSelInd <- length(selIndName)

selectIndGrid <- as.matrix(expand.grid(rep(list(1:nSelInd), 2)))
crossCombGrid <- selectIndGrid[selectIndGrid[, 1] <= selectIndGrid[, 2], ]

### name of each cross combination
crossCombNameMat <- t(apply(X = crossCombGrid, MARGIN = 1, FUN = function(combNo) {
  
  sireName <- selIndName[combNo[1]]
  damName <- selIndName[combNo[2]]
  
  crossCombName <- c(sireName, damName)
  
  return(crossCombName)
}))

crossPatternVec <- sapply(X = 1:nrow(crossCombNameMat), FUN = function(combNo) {
  if (crossCombNameMat[combNo, 1] == crossCombNameMat[combNo, 2]) {
    crsPattern <- "Selfing"
  } else {
    whichPop <- str_sub(string = crossCombNameMat[combNo, ], start = 1, end = 4)
    if (whichPop[1] == whichPop[2]) {
      crsPattern <- "Within"
    } else {
      crsPattern <- "Across"
    }
  }
  return(crsPattern)
})
crossCombDf <- data.frame(crossCombNameMat, crossPatternVec)
colnames(crossCombDf) <- c("sireID", "damID", "crossPattern")




###### 4. Save ######
saveDir <- paste0(scrDir, selMethod, "_nClust_", nClust, "/")
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}
### name of selected individuals
fileNameSelNameAll <- paste0(saveDir, scriptID, "2Pops_selectedIndName_", 
                             selMethod, "_nClust_", nClust, ".rds")
saveRDS(object = selIndName, file = fileNameSelNameAll)

### cross combinations
fileNameCrsComb <- paste0(saveDir, scriptID, "crossCombTable.csv")
write.csv(x = crossCombDf, file = fileNameCrsComb)
