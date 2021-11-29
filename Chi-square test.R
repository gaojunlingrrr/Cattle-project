#load genotyping data
mysample<- read.table("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep = "\t", na.strings = c("","-" ,"NA","no"), fill = TRUE, header = TRUE)
#get all genotped animals
mysample <- mysample[, -c(13,14)]
myanimals<- mysample[,"ID"]

#Load raw data of mastitis
mydata <- read.table("krankData-Dedelow.txt", sep = "\t",fill = TRUE, header = TRUE)

#sick cows
ismastitis = which(grepl("1.13", mydata[, "STAUFEN"], fixed = TRUE))
Mastitisdata = mydata[ismastitis, ]
Mastitisdata= Mastitisdata[which(Mastitisdata[, "OHR"] %in% myanimals), ]
Mastitisdata<- Mastitisdata[, -(3:33)]
CM<- rep(1, length = nrow(Mastitisdata))
Mastitisdata<- cbind(Mastitisdata, CM)
sickanimal<- Mastitisdata[, "OHR"]
##healthy cows
healthydata = mydata[-ismastitis, ]
healthydata= healthydata[which(healthydata[, "OHR"] %in% myanimals), ]
healthydata<-healthydata[, -(3:33)]
healthydata<-healthydata[which(!healthydata[, "OHR"] %in% sickanimal), ]
CM<- rep(0, length = nrow(healthydata))
healthydata<- cbind(healthydata, CM)
CMdata<- rbind(Mastitisdata, healthydata)


###load herd data
library(foreign)
bestand<- read.dbf("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/BESTAND.DBF")

#load genotyping data
gts<- read.table("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep="\t",  na.strings = c("","-" ,"NA","no"), fill = TRUE,header = TRUE,colClasses = "character")

#merge data
mmatrix <- NULL
for(r in 1:nrow(CMdata)){
  ohr <- CMdata[r, "OHR"]
  inBestand <- which(as.character(bestand[, "OHR"]) == as.character(ohr))
  birthdate <- NA
  firstcalf <- NA
  if(length(inBestand) == 1){ 
    birthdate <- as.character(bestand[inBestand, "GEBURT"])
    firstcalf <- as.character(bestand[inBestand, "KALBUNG1"])
    firstcalfindays <- as.character(as.Date(firstcalf, format="%Y-%m-%d") - as.Date(birthdate, format="%Y-%m-%d"))
    Calvingy<- as.character(bestand[inBestand, "KALBUNGN"])
    fatherid<- as.character(bestand[inBestand, "VATER"])
  }
  
  mlpdata <- c(as.character(CMdata[r, "OHR"]), birthdate, firstcalfindays,firstcalf,fatherid,
               as.character(CMdata[r, "LAKTATION"]),
               as.character(CMdata[r, "CM"]))
  idx = which(as.character(gts[, "ID"]) == as.character(CMdata[r, "OHR"]))
  mlpdata <- c(mlpdata, as.character(gts[idx[1], -1]))
  
  mmatrix <- rbind(mmatrix, mlpdata)
  
}


colnames(mmatrix) <- c("AnimalID", "Birthdate", "FirstCalfIndays","Calvingyear","Father", "Lactation","CM" , colnames(gts)[-1])
rownames(mmatrix) <- 1:nrow(mmatrix)


#setwd("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/Rcode update-05-05-2021/")
#write.table(mmatrix, "CM dataset.txt", sep = "\t", quote= FALSE) ###quote 是否为字符型变量添加双引号


#load CM data
cmdata<- read.table("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/Rcode update-05-05-2021/CM dataset.txt", sep = "\t",header = TRUE)

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




####Load SNP data
genotypes <- cmdata[,9:18]

#transform homozygous major alleles to "AA", homozygous minor alleleS to "BB", heterozygous alleles to "AB"

#  markers in column, ind in rows
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
    v1[marker== paste0(majorAllele,minorAllele) |marker== paste0(minorAllele,majorAllele)] <- "AB"
    v1[marker== paste0(minorAllele,minorAllele)] <- "BB"
    return(v1)
  }
   genotrans<- cbind(genotrans, getgeno(marker))
   
}
  colnames(genotrans)<- colnames(genotypes)
  rownames(genotrans)<- 1:nrow(genotrans)
  

write.table(genotrans, "genotrans.txt", sep = "\t")

#Merge CM code and animaID into genotrans matrix
genotrans<- read.table("genotrans.txt", sep = "\t",header = TRUE)
CMcode<- cmdata[, "CM"]
Animal<- cmdata[, "AnimalID"]
chrdata<-cbind(CMcode, Animal, genotrans)

write.table(chrdata, "chrda.txt", sep = "\t")
# Create a matrix with NAs to hold the computed frequencies 10 rows (markers 1 to 10, named by the column names of cmdata), 3 columns called (AA,AB,BB)
outputFreqMatrix <- matrix(NA, 10, 3, dimnames = list(colnames(genotrans), c("AA", "AB", "BB")))

for(marker in colnames(chrdata[, 3:12])) {
  
  # Now compute the frequency using the 28 lines of code you used before
  AA<- chrdata[which(chrdata[, marker]== "AA"), ]
  AB<- chrdata[which(chrdata[, marker]== "AB"), ]
  BB<- chrdata[which(chrdata[, marker]== "BB"), ]
  sickAA<- AA[which(AA[, "CMcode"]==1), ]
  sickAB<- AB[which(AB[, "CMcode"]==1), ]
  sickBB<- BB[which(BB[, "CMcode"]==1), ]
  NAA<-length(unique(AA[, "Animal"]))
  NAB<-length(unique(AB[, "Animal"]))
  NBB<-length(unique(BB[, "Animal"]))
  NALL<-NAA+NAB+NBB
  ##Sick Animal frequency 
  freqAA<-NAA/NALL
  freqAB<-NAB/NALL
  freqBB<-NBB/NALL
  
  
  
  # Add the code to the output matrix
  outputFreqMatrix[marker, "AA"] <- freqAA
  outputFreqMatrix[marker, "AB"] <- freqAB
  outputFreqMatrix[marker, "BB"] <- freqBB
}

rownames(outputFreqMatrix)<- c("SaS_M2:Chr06_rs41588957", "SaS_M3:Chr06_rs110707460", "HAS_M6:Chr13_rs109934030", "HAS_M3:Chr13_rs41634110", "HAS_M4:Chr13_rs109441194", "HAS_M8:Chr19_rs41257403", "HAS_M9:Chr19_rs41636878", "HAS_M1:Chr05_rs41257360", "HAS_M7:Chr18_rs29020544", "HAS_M10:ChrX_rs41629005")
write.table(outputFreqMatrix, file = "genotypeFrequenciesSickCows.txt", sep = "\t")







###Pathogen cases with given genotypes

#load raw data of bacteria
bacteria<- read.table("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Pathogen data of Dedelow.txt", sep = "\t",header = TRUE)

#load herd data
library(foreign)
bestand<- read.dbf("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/BESTAND.DBF")
#order by date
ordering<- order(bestand[, "GEBURT"])
stalldata<- bestand[ordering, ]
stalldata<- cbind(stalldata, "year" = unlist(lapply(strsplit(as.character(stalldata[, "GEBURT"]),"-" ),"[", 1)))
#select birthyear > 2008
birthrange<- which(as.numeric(as.character(stalldata[, "year"])) > 2008)
stalldata<- stalldata[birthrange,]
#stallnumber
stallnumer<-bacteria[, "Stallnummer"]
stalldata<- stalldata[which(stalldata[, "STALL"] %in% stallnumer), ]
Animals<- stalldata[, "STALL"]
bacteria<- bacteria[which(bacteria[, "Stallnummer"] %in% Animals), ]
####get animalID and bacteria data
pathogendata<- NULL
for (r in 1:nrow(bacteria)){
  stall<- bacteria[r, "Stallnummer"]
  install<-which(as.character(stalldata[, "STALL"]) == as.character(stall))
  AnimalID<-NA
  if(length(install) ==1){
    AnimalID<- as.character(stalldata[install, "OHR"])
  }
  bacteriadata<- c(AnimalID,
                   as.character(bacteria[r, "Stallnummer"]),
                   as.character(bacteria[r, "Betrieb"]),
                   as.character(bacteria[r, "bv"]),
                   as.character(bacteria[r, "KNS"]),
                   as.character(bacteria[r, "Sc_pl"]),
                   as.character(bacteria[r, "Sc_min"]),
                   as.character(bacteria[r, "STA"]),
                   as.character(bacteria[r, "Ecoli"]),
                   as.character(bacteria[r, "Galt"]))
  pathogendata<-rbind(pathogendata, bacteriadata)
  
  
}
colnames(pathogendata)<- c("AnimalID", "STALL", "Farm", "Bacteria", "CNS", "SC+", "SC-", "STA", "E.coli", "GALT")
rownames(pathogendata)<- 1:nrow(pathogendata)


####Remove NA 
pathogendata<-pathogendata[-which(is.na(pathogendata[, "AnimalID"])), ]

###Merge Gene matrix into pathogen data
#load genotrans data
genotrans<- read.table("genotrans.txt", sep = "\t",header = TRUE)

#Merge data
pathogendata<- pathogendata[which(pathogendata[, "AnimalID"] %in% gts[, "ID"]), ]
genematrix<- NULL
for (r in 1:nrow(gts)) {
  ohr<- gts[r, "ID"]
  inPathogen<- which(as.character(pathogendata[, "AnimalID"]) == as.character(ohr))
  if(length(inPathogen)== 1){
    
    CNS<- as.character(pathogendata[inPathogen, "CNS"])
    Scplus<- as.character(pathogendata[inPathogen, "SC+"])
    Scmin<- as.character(pathogendata[inPathogen, "SC-"])
    STA<- as.character(pathogendata[inPathogen, "STA"])
    Ecoli<- as.character(pathogendata[inPathogen, "E.coli"])
    Galt<- as.character(pathogendata[inPathogen, "GALT"])
    
  }
  dataset<- c(as.character(gts[r, "ID"]), CNS, Scplus, Scmin, STA, Ecoli, Galt,
              as.character(gts[r, "SaS_M2"]),
              as.character(gts[r, "SaS_M3"]),
              as.character(gts[r, "HAS_M6"]),
              as.character(gts[r, "HAS_M3"]),
              as.character(gts[r, "HAS_M4"]),
              as.character(gts[r, "HAS_M8"]),
              as.character(gts[r, "HAS_M9"]),
              as.character(gts[r, "HAS_M1"]),
              as.character(gts[r, "HAS_M7"]),
              as.character(gts[r, "HAS_M10"]))
  
  genematrix<- rbind(genematrix, dataset)
}
colnames(genematrix)<- c("AnimalID", "CNS", "SC+", "SC-", "STA", "E.coli", "GALT", "CHR6-1", "CHR6-2","CHR13-1", "CHR13-2", "CHR13-3", "CHR19-1", "CHR19-2", "CHR5","CHR18", "CHRX")
rownames(genematrix)<- 1:nrow(genematrix)

#load CM data
cmdata<- read.table("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/Rcode update-05-05-2021/CM dataset.txt", sep = "\t",header = TRUE)

#only select mastitis cows
cmcowdata<- cmdata[which(cmdata[, "CM"]==1), ]
cmcowsid<- cmcowdata[, "AnimalID"]
genematrix<- genematrix[which(genematrix[, "AnimalID"] %in% cmcowsid), ]

write.table(genematrix, "genematrix.txt", sep = "\t", col.names= TRUE, quote= FALSE)



#loading dataset

genematrix<- read.table("genematrix.txt", sep = "\t",check.names= F, header = TRUE)

#transform homozygous major alleles to "AA", homozygous minor alleleS to "BB", heterozygous alleles to "AB"

#  markers in column, ind in rows
pathotrans<- NULL
for(id in colnames(genematrix[, 8:17])){
  marker <- genematrix[,id]
  splitted <- unlist(strsplit(as.character(marker), ""))
  mtab <- table(splitted)
  minorAllele <- names(mtab)[which.min(mtab)]
  majorAllele<- names(mtab)[which.max(mtab)]
  getgeno<- function(marker){
    v1<- rep(NA, length(marker))
    v1[marker== paste0(majorAllele,majorAllele)] <- "AA"
    v1[marker== paste0(majorAllele,minorAllele) |marker== paste0(minorAllele,majorAllele)] <- "AB"
    v1[marker== paste0(minorAllele,minorAllele)] <- "BB"
    return(v1)
  }
  pathotrans<- cbind(pathotrans, getgeno(marker))
  
}
colnames(pathotrans)<- colnames(genematrix[, 8:17])

#Merging pathogen data from genematrix
pathontrans1<- cbind(genematrix[, 1:7], pathotrans)


write.table(pathontrans1, "pathontans1.txt", sep = "\t")

### calculating pathogen cases of different SNPs

#creating empty matrix
bacfreq<- c()
for (pathogen in colnames(pathontrans1[, 2:7])) {
  
  idx<- which(pathontrans1[, pathogen]==1)
  pathoSubset<- pathontrans1[idx, ]
  
  for (marker in colnames(pathontrans1[, 8:17])) {
    snpdata<- pathoSubset[, marker]
  
  nAA<- length(which(snpdata == "AA"))
  nAB<- length(which(snpdata == "AB"))
  nBB<- length(which(snpdata == "BB"))
  bacfreq<-rbind(bacfreq, c(pathogen, marker, nAA, nAB, nBB ))
  
  }
}
colnames(bacfreq)<- c("Bacteria", "Chr", "nAA", "nAB", "nBB")
write.table(bacfreq, "bacfreq.txt", sep = "\t")






##chi-square test
#load number of pathogen cases


bacfreq<- read.table("bacfreq.txt", sep = "\t", header = TRUE)

#load sick cow frequency
sickfreq<- read.table("genotypeFrequenciesSickCows.txt", sep = "\t", header = TRUE)

#only pick CNS
cnsdata<- bacfreq[which(bacfreq[, "Bacteria"] == "CNS"), ]

pvals1<- c()
for (x in 1:10) {
    
  pval<- chisq.test(c(cnsdata[x,3], cnsdata[x,4],cnsdata[x,5]), p=c(sickfreq[x, 1], sickfreq[x, 2], sickfreq[x, 3]))$p.value
  
    pvals1<- round(p.adjust(c(pvals1, pval), method = "BH"), 4)
}

pvals1
[1] 0.7188 0.7188 0.7188 0.7188 0.7188 0.7188 0.7188 0.7188 0.7188 0.7188


# only pick SC+ (Streptococcus uberis)
SCplus<- bacfreq[which(bacfreq[, "Bacteria"] == "SC+"), ]

pvals2<- c()
for (x in 1:10) {
  
  pval<- chisq.test(c(SCplus[x,3], SCplus[x,4],SCplus[x,5]), p=c(sickfreq[x, 1], sickfreq[x, 2], sickfreq[x, 3]))$p.value
  
  pvals2<- round(p.adjust(c(pvals2, pval), method= "BH"),4)
}

pvals2
[1] 0.9544 0.9544 0.9544 0.0000 0.9544 0.9544 0.9544 0.9544 0.9544 0.9544


#only pick SC-(S. dysgalactiae)

SCminus<- bacfreq[which(bacfreq[, "Bacteria"] == "SC-"), ]

pvals3<- c()
for (x in 1:10) {
  
  pval<- chisq.test(c(SCminus[x,3], SCminus[x,4],SCminus[x,5]), p=c(sickfreq[x, 1], sickfreq[x, 2], sickfreq[x, 3]))$p.value
  
  pvals3<- round(p.adjust(c(pvals3, pval), method= "BH"),4)
}

pvals3
[1] 0.9707 0.9707 0.9707 0.9707 0.9707 0.9707 0.9707 0.9707 0.9707 0.9707



#only pick STA (S. aureus)

STA<- bacfreq[which(bacfreq[, "Bacteria"] == "STA"), ]

pvals4<- c()
for (x in 1:10) {
  
  pval<- chisq.test(c(STA[x,3], STA[x,4],STA[x,5]), p=c(sickfreq[x, 1], sickfreq[x, 2], sickfreq[x, 3]))$p.value
  
  pvals4<- round(p.adjust(c(pvals4, pval), method= "BH"),4)
}

pvals4
[1] 0.7339 0.7339 0.7339 0.7339 0.7339 0.7339 0.7339 0.7339 0.7339 0.7339


#only pick E.COLI 

ECO<- bacfreq[which(bacfreq[, "Bacteria"] == "E.coli"), ]

pvals5<- c()
for (x in 1:10) {
  
  pval<- chisq.test(c(ECO[x,3], ECO[x,4],ECO[x,5]), p=c(sickfreq[x, 1], sickfreq[x, 2], sickfreq[x, 3]))$p.value
  
  pvals5<- round(p.adjust(c(pvals5, pval), method= "BH"),4)
}


pvals5
[1] 0.8828 0.8828 0.8828 0.8828 0.8828 0.8828 0.8828 0.8828 0.8828 0.0422

#only pick GALT (Streptococcus agalactiae)

Galt<- bacfreq[which(bacfreq[, "Bacteria"] == "GALT"), ]

pvals6<- c()
for (x in 1:10) {
  
  pval<- chisq.test(c(Galt[x,3], Galt[x,4],Galt[x,5]), p=c(sickfreq[x, 1], sickfreq[x, 2], sickfreq[x, 3]))$p.value
  
  pvals6<- round(p.adjust(c(pvals6, pval), method= "BH"),4)
}

pvals6
[1] 0.9214 0.9214 0.9214 0.9214 0.9214 0.9214 0.9214 0.9214 0.9214 0.9214


psum<-cbind(pvals1, pvals2, pvals3, pvals4, pvals5, pvals6)
rownames(psum)<- rownames(sickfreq)
colnames(psum)<- c("CNS", "SC+", "SC-","STA","E.COLI","GALT")

#results
psum

                            CNS    SC+    SC-    STA  E.COLI  GALT
SaS_M2:Chr06_rs41588957  0.7188 0.9544 0.9707 0.7339 0.8828 0.9214
SaS_M3:Chr06_rs110707460 0.7188 0.9544 0.9707 0.7339 0.8828 0.9214
HAS_M6:Chr13_rs109934030 0.7188 0.9544 0.9707 0.7339 0.8828 0.9214
HAS_M3:Chr13_rs41634110  0.7188 0.0000 0.9707 0.7339 0.8828 0.9214
HAS_M4:Chr13_rs109441194 0.7188 0.9544 0.9707 0.7339 0.8828 0.9214
HAS_M8:Chr19_rs41257403  0.7188 0.9544 0.9707 0.7339 0.8828 0.9214
HAS_M9:Chr19_rs41636878  0.7188 0.9544 0.9707 0.7339 0.8828 0.9214
HAS_M1:Chr05_rs41257360  0.7188 0.9544 0.9707 0.7339 0.8828 0.9214
HAS_M7:Chr18_rs29020544  0.7188 0.9544 0.9707 0.7339 0.8828 0.9214
HAS_M10:ChrX_rs41629005  0.7188 0.9544 0.9707 0.7339 0.0422 0.9214




