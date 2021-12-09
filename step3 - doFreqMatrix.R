### SCRIPT3.R: New script, continue here when we have "genotrans.txt" and 

#Merge CM code and animaID into genotrans matrix
setwd("C:/Users/Arends/Downloads/R.chisquare_01.11.2021_Junling/R.chisquare_01.11.2021_Junling/")
genotrans<- read.table("genotrans.txt", sep = "\t",header = TRUE)
cmdata<- read.table("CM dataset.txt", sep = "\t",header = TRUE)

genotrans <- cbind(genotrans, CM = NA)
for(animal in rownames(genotrans)){
  CMcode <- unique(cmdata[which(cmdata[, "AnimalID"] ==  animal), "CM"])
  if(length(CMcode) == 1){
    genotrans[animal, "CM"] <- CMcode
  }
}

# Create a matrix with NAs to hold the computed frequencies 10 rows (markers 1 to 10, named by the column names of cmdata), 3 columns called (AA,AB,BB)
outputFreqMatrixSick <- matrix(NA, 10, 3, dimnames = list(colnames(genotrans)[1:10], c("AA", "AB", "BB")))
outputFreqMatrixHealthy <- matrix(NA, 10, 3, dimnames = list(colnames(genotrans)[1:10], c("AA", "AB", "BB")))

for(marker in rownames(outputFreqMatrixSick)) {
  # Now compute the frequency using the 28 lines of code you used before
  AA <- genotrans[which(genotrans[, marker]== "AA"), ]
  AB <- genotrans[which(genotrans[, marker]== "AB"), ]
  BB <- genotrans[which(genotrans[, marker]== "BB"), ]
  sickAA <- AA[which(AA[, "CM"]==1), ]
  sickAB <- AB[which(AB[, "CM"]==1), ]
  sickBB <- BB[which(BB[, "CM"]==1), ]
  NALL <- nrow(sickAA) + nrow(sickAB) + nrow(sickBB)
  ##Sick Animal frequency 
  outputFreqMatrixSick[marker, "AA"] <- nrow(sickAA)/NALL
  outputFreqMatrixSick[marker, "AB"] <- nrow(sickAB)/NALL
  outputFreqMatrixSick[marker, "BB"] <- nrow(sickBB)/NALL
  
  healthyAA <- AA[which(AA[, "CM"]==0), ]
  healthyAB <- AB[which(AB[, "CM"]==0), ]
  healthyBB <- BB[which(BB[, "CM"]==0), ]
  NALL <- nrow(healthyAA) + nrow(healthyAB) + nrow(healthyBB)

  outputFreqMatrixHealthy[marker, "AA"] <- nrow(healthyAA)/NALL
  outputFreqMatrixHealthy[marker, "AB"] <- nrow(healthyAB)/NALL
  outputFreqMatrixHealthy[marker, "BB"] <- nrow(healthyBB)/NALL
}

rownames(outputFreqMatrixSick)<- c("SaS_M2:Chr06_rs41588957", "SaS_M3:Chr06_rs110707460", "HAS_M6:Chr13_rs109934030", "HAS_M3:Chr13_rs41634110", "HAS_M4:Chr13_rs109441194", "HAS_M8:Chr19_rs41257403", "HAS_M9:Chr19_rs41636878", "HAS_M1:Chr05_rs41257360", "HAS_M7:Chr18_rs29020544", "HAS_M10:ChrX_rs41629005")
write.table(outputFreqMatrixSick, file = "genotypeFrequenciesSickCows.txt", sep = "\t")

rownames(outputFreqMatrixHealthy)<- c("SaS_M2:Chr06_rs41588957", "SaS_M3:Chr06_rs110707460", "HAS_M6:Chr13_rs109934030", "HAS_M3:Chr13_rs41634110", "HAS_M4:Chr13_rs109441194", "HAS_M8:Chr19_rs41257403", "HAS_M9:Chr19_rs41636878", "HAS_M1:Chr05_rs41257360", "HAS_M7:Chr18_rs29020544", "HAS_M10:ChrX_rs41629005")
write.table(outputFreqMatrixHealthy, file = "genotypeFrequenciesHealthyCows.txt", sep = "\t")
