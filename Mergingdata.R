mysample<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep = "\t", na.strings = c("", "NA","no"), fill = TRUE, header = TRUE)
mysample <- mysample[, -c(13,14)]



###mastitis data
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

# Select all the entries  between 2009 and 2017
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

#Per Lactation
inL1 <- which(Mastitisdata[,"LAKTATION"] == 1)
inL2 <- which(Mastitisdata[,"LAKTATION"] == 2)
inL3 <- which(Mastitisdata[,"LAKTATION"] == 3)
inL4 <- which(Mastitisdata[,"LAKTATION"] > 3)

ML1data<- Mastitisdata[inL1,]
ML2data<- Mastitisdata[inL2,]
ML3data<- Mastitisdata[inL3,]
ML4data<- Mastitisdata[inL4,]

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





#### Create mergingdata dataset
Mergingdata<-MastitisInfo[, c("OHR", "nMastitis")]
Mergingdata[1:5,]

getCode<- function(Mastitiscases){
  nmastitis<- as.numeric(Mastitiscases)
  ret<- rep(NA, length(nmastitis))
  ret[nmastitis == 1] <-  "1"
  ret[nmastitis > 1]  <-  "1"
  return(ret)
  
}

Mastitiscases<- Mergingdata[, "nMastitis"]
animalcode<- getCode(Mastitiscases)
Mergingdata<- cbind(Mergingdata, animalcode)
colnames(Mergingdata)<- c("AnimalID", "nMastitis", "Code")
Mergingdata<- Mergingdata[, c("AnimalID", "Code")]


####My healthy animal
#Delete sick cows from 305 file 

myanimals<- mysample[,"ID"]

lactData <- read.table("/Users/Pomelo/Documents/R Practising/Dedelow/305data-Dedelow.txt", sep = "\t", header = TRUE)
lactData<- lactData[which(lactData[, "OHR"] %in% myanimals),]
# Covert dates to R format, and store as new column
KdateObjects = as.Date(as.character(lactData[, "KALBUNG"]), format="%y/%m/%d")
lactData = cbind(lactData, "Kdate" = KdateObjects) ### lapply function return a list 
lactData = cbind(lactData, "year" = unlist(lapply(strsplit(as.character(lactData[,"Kdate"]), "-"), "[", 1)))
lactData = cbind(lactData, "month" = unlist(lapply(strsplit(as.character(lactData[,"Kdate"]), "-", fixed = TRUE), "[", 2)))

## dates to season

getSeason<- function(lactData){
  dates<- as.numeric(lactData[, "month"]) ###factor->numeric
  ret<-rep(NA, length(dates))
  ret[dates>=4 & dates<=5] <- "Spring"
  ret[dates>=6 & dates<=8] <- "Summer"
  ret[dates>=9 & dates<=10] <- "Fall"
  ret[dates==11 | dates==12 | dates==1 | dates==2 | dates==3 ] <- "Winter"
  
  return(ret)
}

Seasonobjects <- getSeason(lactData)

lactData <- cbind(lactData, "Season" =Seasonobjects )

healthyCows = lactData[which(!(lactData[,"OHR"] %in% Mastitisdata[, "OHR"])), ]

healthycowid<- as.matrix(unique(healthyCows[,"OHR"])) ###必须先定义为矩阵后才能使用cbind

healthycode<- as.matrix(rep(0, length=length(healthycowid) ))
healthydata<- cbind(healthycowid, healthycode)
colnames(healthydata)<- c("AnimalID", "Code")

Mergingdata<- rbind(Mergingdata, healthydata)
Animalcode<- Mergingdata[, "Code"]

Animals<- Mergingdata[, "AnimalID"]
mysample<-mysample[which(unique(myanimals) %in% Animals), ]
Mergingdata<- cbind(mysample, Animalcode)




###Bestand file
library(foreign)
Bestanddata<- read.dbf("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/BESTAND.DBF")
ordering <- order(Bestanddata[, "GEBURT"])
Bestanddata<- Bestanddata[ordering, ]
Bestanddata<- cbind(Bestanddata, "year"= unlist(lapply(strsplit(as.character(Bestanddata[,"GEBURT"]), "-"), "[", 1)))
Bestanddata = cbind(Bestanddata, "month" = unlist(lapply(strsplit(as.character(Bestanddata[,"GEBURT"]), "-", fixed = TRUE), "[", 2)))




## dates to season

getSeason<- function(Bestanddata){
  dates<- as.numeric(Bestanddata[, "month"]) ###factor->numeric
  ret<-rep(NA, length(dates))
  ret[dates>=4 & dates<=5] <- "Spring"
  ret[dates>=6 & dates<=8] <- "Summer"
  ret[dates>=9 & dates<=10] <- "Fall"
  ret[dates==11 | dates==12 | dates==1 | dates==2 | dates==3 ] <- "Winter"
  
  return(ret)
}

Seasonobjects <- getSeason(Bestanddata)

Bestanddata <- cbind(Bestanddata, "Season" =Seasonobjects )


Bestanddata<- Bestanddata[which(Bestanddata[, "OHR"] %in% Animals), ]
Birthyear<- Bestanddata[, "year"]
Birthseason<- Bestanddata[, "Season"]

Mergingdata<- cbind(Mergingdata, Birthyear, Birthseason)

setwd("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/")
write.table(Mergingdata, file= "Mergingdata.txt", sep = "\t")

#### Age at first of calving
Bestanddata<- cbind(Bestanddata, "Afc"= unlist(lapply(strsplit(as.character(Bestanddata[,"KALBUNG1"]), "-"), "[", 1)))
Bestanddata = cbind(Bestanddata, "Afcmonth" = unlist(lapply(strsplit(as.character(Bestanddata[,"KALBUNG1"]), "-", fixed = TRUE), "[", 2)))

getSeason<- function(Bestanddata){
  dates<- as.numeric(Bestanddata[, "Afcmonth"]) ###factor->numeric
  ret<-rep(NA, length(dates))
  ret[dates>=4 & dates<=5] <- "Spring"
  ret[dates>=6 & dates<=8] <- "Summer"
  ret[dates>=9 & dates<=10] <- "Fall"
  ret[dates==11 | dates==12 | dates==1 | dates==2 | dates==3 ] <- "Winter"
  
  return(ret)
}

Seasonobjects <- getSeason(Bestanddata)

Bestanddata <- cbind(Bestanddata, "AfcSeason" =Seasonobjects )

Afc<-Bestanddata[, "Afc"]
Afcs<-Bestanddata[, "AfcSeason"]
Mergingdata<-cbind(Mergingdata,Afc,Afcs)
###father
Father<-Bestanddata[,"VATER"]
Mergingdata<- cbind(Mergingdata, Father)




setwd("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/")
write.table(Mergingdata, file= "Mergingdata.txt", sep = "\t")



