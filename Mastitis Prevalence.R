mmatrix<-read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/phenotypes2020.Dedelow.txt", sep = "\t",header = TRUE)
Mergingdata<-read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Mergingdata.txt",sep = "\t",header = TRUE)

phenotypes<- NULL
for (i in 1:nrow(mmatrix)){
  ohr1<- mmatrix[i,"AnimalID"]
  inMerging<-which(as.character(Mergingdata[, "ID"]) == as.character(ohr1))
  Code<- NA
  if(length(inMerging) == 1){
    Code<- as.numeric(Mergingdata[inMerging,"Animalcode"])
  }
  mastitiscode<- c(as.character(mmatrix[i, "AnimalID"]), Code )
  phenotypes<- rbind(phenotypes, mastitiscode)
}
colnames(phenotypes)<-c("AnimalIDD", "Code")
row.names(phenotypes)<- 1:nrow(phenotypes)
Animalcode<-as.numeric(phenotypes[,"Code"])
mmatrix<-mmatrix[,-13]
mmatrix1<- cbind(mmatrix,Animalcode)
####Get season
samplemonth<- unlist(lapply(strsplit(as.character(mmatrix1[, "Sampledate"]), "-"), "[", 2))
samplemonth<- as.numeric(samplemonth)
getSeason <- function(samplemonths) {
  ret <- rep(NA, length(samplemonths))
  ret[samplemonths >= 3 & samplemonths <= 5] <- "Spring"
  ret[samplemonths >= 6 & samplemonths <= 8] <- "Summer"
  ret[samplemonths >= 9 & samplemonths <= 11] <- "Fall"
  ret[samplemonths == 12 | samplemonths == 1 | samplemonths == 2] <- "Winter"
  return(ret)
}
Samplingseason<- getSeason(samplemonth)

calvingmonth<- unlist(lapply(strsplit(as.character(mmatrix1[, "calvingyear"]), "-"), "[", 2))
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

mmatrix1<-cbind(mmatrix1, Samplingseason, Calvingseason)
setwd("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/")
write.table(mmatrix1, "phenotypes2023.Dedelow.txt", sep = "\t", quote= FALSE) ###
####Delete animals leave before 2017
mmatrix1<-read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/phenotypes2023.Dedelow.txt",sep = "\t",header = TRUE)
mmatrix2<- mmatrix1[-which(unlist(lapply(strsplit(as.character(mmatrix1[,"Leaveherdtime"]),"-"), "[", 1)) <2017), ]### nrow()=0
# Get the phenotype we are analyzing
mmatrix1<- mmatrix1[which(!is.na(mmatrix1[, "Animalcode"])), ]##去除NA行

### select lactation 1-3
mmatrix1<- mmatrix1[which(mmatrix1[, "Lactation"] < 4 & mmatrix1[, "Lactation"] > 0), ]
####select animal with birth year >2007
birthyearofcows <- as.numeric(unlist(lapply(strsplit(as.character(mmatrix1[,"Birthdate"]), "-"), "[", 1)))
mmatrix1<- cbind(mmatrix1, birthyearofcows)
mmatrix1<- mmatrix1[which(mmatrix1[, "birthyearofcows"] > 2007), ]

###Delete elements with less than 3 reccords
frequency<- as.data.frame(table(mmatrix1[, "AnimalID"]))
frequency2<- frequency[which(frequency[, "Freq"]>2),]
qfanimal<- frequency2[, "Var1"]
mmatrix1<- mmatrix1[which(mmatrix1[, "AnimalID"] %in% qfanimal), ]

SCS <- log2(as.numeric(as.character(mmatrix1[, "SCC"])) / 100) + 3
SCS[which(!is.finite(SCS))] <- 0
# Get the covariates we want to investigate in model building
animal <- mmatrix1[, "AnimalID"]
lactation <- as.factor(mmatrix1[, "Lactation"])
year <- unlist(lapply(strsplit(as.character(mmatrix1[,"Sampledate"]), "-"), "[", 1))
birthyear <- unlist(lapply(strsplit(as.character(mmatrix1[,"Birthdate"]), "-"), "[", 1))
seasonsampling <- mmatrix1[, "Samplingseason"]
seasoncalving<- mmatrix1[, "Calvingseason"]
yearcalving <-  unlist(lapply(strsplit(as.character(mmatrix1[,"calvingyear"]), "-"), "[", 1))
firstCalf <- mmatrix1[, "FirstCalfIndays"]
father<- mmatrix1[, "Father"]
animalcode<- mmatrix1[, "Animalcode"] 

####select animal with birth year >2007

# Combine all data into a matrix
mdata <- data.frame(Y = as.numeric(SCS), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Samplingyear = as.factor(year),
                    Birthyear = as.factor(birthyear),
                    Samplingseason = as.factor(seasonsampling),
                    Calvingseason = as.factor(seasoncalving),
                    Calvingyear = as.factor(yearcalving),
                    Firstcalf = as.numeric(firstCalf),
                    Father= as.factor(father)
)
# Remove rows with missing data
mdata<- na.omit(mdata)

#Improve the association by removing fathers which do not have 10 offspring
offspring <- c()
for(f in unique(mdata[, "Father"])){
  nOffspring = length(unique(mdata[which(mdata[, "Father"] == f), "Animal"]))
  offspring <- c(offspring, nOffspring)
}
names(offspring) <- unique(mdata[, "Father"])

enoughoffspring <- names(offspring)[which(offspring >= 10)]
dim(mdata)
mdata <- mdata[which(mdata[, "Father"] %in% enoughoffspring),]
dim(mdata)

###pick animals used in association study

Cowinas<- mdata[, "Animal"]


####Calculate Mastitis Prevalence
mysample<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep = "\t", fill = TRUE, header = TRUE)
mysample <- mysample[, -c(13,14)]
myanimals<- mysample[,"ID"]
setwd("/Users/Pomelo/Documents/R Practising/Dedelow/")
mydata <- read.table("krankData-Dedelow.txt", sep = "\t",fill = TRUE, header = TRUE)
ismastitis = which(grepl("1.13", mydata[, "STAUFEN"], fixed = TRUE))
Mastitisdata = mydata[ismastitis, ]
Mastitisdata= Mastitisdata[which(Mastitisdata[, "OHR"] %in% myanimals), ]
# Covert dates to R format, and store as new column
Rdaterows = as.Date(as.character(Mastitisdata[, "DATUM"]), format="%y/%m/%d")
Mastitisdata = cbind(Mastitisdata, "Rdate"= Rdaterows)

#order by date
ordering = order(Mastitisdata[,"Rdate"])
Mastitisdata = Mastitisdata[ordering, ]
Mastitisdata = cbind(Mastitisdata, "year" = unlist(lapply(strsplit(as.character(Mastitisdata[,"Rdate"]), "-"), "[", 1)))


##
Mastitisdata<- Mastitisdata[which(Mastitisdata[, "OHR"] %in% Cowinas), ]

# Select all the entries  between 2009 and 2017  (the average productive life of a Holstein is approximately four years).
Timerange = which(as.numeric(as.character(Mastitisdata[, "year"])) > 2008 & as.numeric(as.character(Mastitisdata[, "year"])) < 2018 )
Mastitisdata = Mastitisdata[Timerange, ]

in2009 <- which(grepl("2009", Mastitisdata[,"Rdate"]))
in2010 <- which(grepl("2010", Mastitisdata[,"Rdate"]))
in2011 <- which(grepl("2011", Mastitisdata[,"Rdate"]))
in2012 <- which(grepl("2012", Mastitisdata[,"Rdate"]))
in2013 <- which(grepl("2013", Mastitisdata[,"Rdate"]))
in2014 <- which(grepl("2014", Mastitisdata[,"Rdate"]))
in2015 <- which(grepl("2015", Mastitisdata[,"Rdate"]))
in2016 <- which(grepl("2016", Mastitisdata[,"Rdate"]))
in2017 <- which(grepl("2017", Mastitisdata[,"Rdate"]))

generateInfo = function(Mastitisdata){
  Mastitistable = table(as.character(Mastitisdata[, "OHR"]))
  MastitisInfo <- NULL
  Mastitisentry <- vector("list", length(names(Mastitistable)))
  names(Mastitisentry) = names(Mastitistable)
  for (cow in names(Mastitistable)){
    cowID = which(Mastitisdata[, "OHR"] == cow)
    Mastitisentry[[cow]] = Mastitisdata[cowID, "Rdate"]
    nDiagnosis = length(Mastitisentry[[cow]])
    Diagnosisdate = Mastitisentry[[cow]][1]
    if (nDiagnosis > 1){
      Diagnosis = c("N", rep(NA, nDiagnosis-1))
      for(x in 2:nDiagnosis){
        if ((Mastitisentry[[cow]][x] - Diagnosisdate) < 14){
          Diagnosis[x] <- "F"
        }
        else{
          Diagnosis[x] <- "N"
        }
        Diagnosisdate = Mastitisentry[[cow]][x]
      }
      Mastitisentry[[cow]] = rbind(as.character(Mastitisentry[[cow]]), Diagnosis)
      nMastitis = length(which(Diagnosis =="N"))
      MastitisInfo = rbind(MastitisInfo, c(cow,nDiagnosis,nMastitis))
    }else{
      MastitisInfo = rbind(MastitisInfo, c(cow,1,1))
    }
  }
  
  colnames(MastitisInfo) <- c("OHR", "nDiagnosis", "nMastitis")
  return(MastitisInfo)
}
MastitisInfo = generateInfo(Mastitisdata)
sum(as.numeric(MastitisInfo[,"nDiagnosis"])) / nrow(MastitisInfo)
sum(as.numeric(MastitisInfo[,"nMastitis"])) / nrow(MastitisInfo)

####
#Per Lactation
inL1 <- which(Mastitisdata[,"LAKTATION"] == 1)
inL2 <- which(Mastitisdata[,"LAKTATION"] == 2)
inL3 <- which(Mastitisdata[,"LAKTATION"] == 3)



lactation = 1
for(individuals in list(inL1,inL2,inL3)){
  MastitisInfo = generateInfo(Mastitisdata[individuals, ])
  cat("Lactation:", lactation)
  cat(", Cows", nrow(MastitisInfo))
  cat(", nDiagnois", sum(as.numeric(MastitisInfo[, "nDiagnosis"])))
  cat(", nMastitis", sum(as.numeric(MastitisInfo[, "nMastitis"])))
  cat(", nDiagnosis:", sum(as.numeric(MastitisInfo[, "nDiagnosis"]))/nrow(MastitisInfo))
  cat(", nMastitis:", sum(as.numeric(MastitisInfo[, "nMastitis"]))/nrow(MastitisInfo),"\n")
  lactation = lactation + 1
  
}



###number of sick and healthy animal
mmatrix1<- mmatrix1[which(mmatrix1[, "AnimalID"] %in% Cowinas), ]
healthygroup<- mmatrix1[which(mmatrix1[, "Animalcode"] == 0), ]
sickgroup<- mmatrix1[which(mmatrix1[, "Animalcode"] == 1), ]
length(unique(healthygroup[, "AnimalID"]))
length(unique(sickgroup[, "AnimalID"]))






















