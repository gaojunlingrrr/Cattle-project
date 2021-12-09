### SCRIPT2.R: Continue script from here, once we have the "CM dataset.txt" file

#load CM data
setwd("C:/Users/Arends/Downloads/R.chisquare_01.11.2021_Junling/R.chisquare_01.11.2021_Junling/")
cmdata<- read.table("CM dataset.txt", sep = "\t",header = TRUE)

####Get season
calvingmonth<- unlist(lapply(strsplit(as.character(cmdata[, "Calvingyear"]), "-"), "[", 2))
calvingmonth<- as.numeric(calvingmonth)
getSeason <- function(calvingmonth) {
  ret <- rep(NA, length(calvingmonth))
  ret[calvingmonth >= 3 & calvingmonth <= 5] <- "Spring"
  ret[calvingmonth >= 6 & calvingmonth <= 8] <- "Summer"
  ret[calvingmonth >= 9 & calvingmonth <= 11] <- "Fall"
  ret[calvingmonth == 12 | calvingmonth == 1 | calvingmonth == 2] <- "Winter"
  return(ret)
}
Calvingseason<- getSeason(calvingmonth)

cmdata<-cbind(cmdata, Calvingseason)

##delete NA value
cmdata<- na.omit(cmdata)

####select animal with birth year >2008
birthyear <- as.numeric(unlist(lapply(strsplit(as.character(cmdata[,"Birthdate"]), "-"), "[", 1)))
cmdata<- cbind(cmdata, birthyear)
cmdata<- cmdata[which(cmdata[, "birthyear"] > 2008), ]

### select lactation 1-3
cmdata<- cmdata[which(cmdata[, "Lactation"] < 4 & cmdata[, "Lactation"] > 0), ]

###Delete elements with less than 3 reccords
frequency<- as.data.frame(table(cmdata[, "AnimalID"]))
frequency2<- frequency[which(frequency[, "Freq"]>2),]
qfanimal<- frequency2[, "Var1"]
cmdata<- cmdata[which(cmdata[, "AnimalID"] %in% qfanimal), ]

##get calving year
Calvingyear_N<- unlist(lapply(strsplit(as.character(cmdata[,"Calvingyear"]), "-"), "[", 1))
cmdata<- cbind(cmdata, Calvingyear_N)

# Make sure that we only use an animal once, so that we compute the real frequency
uniqueAnimals <- unique(cmdata[,"AnimalID"])
genotypes <- NULL
for(animal in uniqueAnimals){
  i <- which(cmdata[, "AnimalID"] == animal)[1]
  genotypes <- rbind(genotypes, cmdata[i,9:18])
}

rownames(genotypes) <- uniqueAnimals

#transform homozygous major alleles to "AA", homozygous minor alleleS to "BB", heterozygous alleles to "AB"
# markers in column, ind in rows
genotrans<- NULL
for(x in 1:ncol(genotypes)){
  marker <- genotypes[,x]
  splitted <- unlist(strsplit(as.character(marker), ""))
  mtab <- table(splitted)
  minorAllele <- names(mtab)[which.min(mtab)]
  majorAllele<- names(mtab)[which.max(mtab)]
  getgeno<- function(marker){
    v1<- rep(NA, length(marker))
    v1[marker== paste0(majorAllele,majorAllele)] <- "AA"
    v1[marker== paste0(majorAllele,minorAllele) | marker== paste0(minorAllele,majorAllele)] <- "AB"
    v1[marker== paste0(minorAllele,minorAllele)] <- "BB"
    return(v1)
  }
   genotrans<- cbind(genotrans, getgeno(marker))
   
}
colnames(genotrans)<- colnames(genotypes)
rownames(genotrans)<- rownames(genotypes)

write.table(genotrans, "genotrans.txt", sep = "\t")
