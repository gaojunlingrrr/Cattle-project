library(foreign)

#load genotyping data
setwd("C:/Users/Arends/Downloads/R.chisquare_01.11.2021_Junling/R.chisquare_01.11.2021_Junling/")
mysample<- read.table("Genotyping data- Dedelow.txt", sep = "\t", na.strings = c("","-" ,"NA","no"), fill = TRUE, header = TRUE)
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
bestand<- read.dbf("BESTAND.DBF")

#load genotyping data
gts<- read.table("Genotyping data- Dedelow.txt", sep="\t",  na.strings = c("","-" ,"NA","no"), fill = TRUE,header = TRUE,colClasses = "character")

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
  
    mlpdata <- c(as.character(CMdata[r, "OHR"]), birthdate, firstcalfindays, firstcalf, fatherid, as.character(CMdata[r, "LAKTATION"]), as.character(CMdata[r, "CM"]))
    idx = which(as.character(gts[, "ID"]) == as.character(CMdata[r, "OHR"]))
    mlpdata <- c(mlpdata, as.character(gts[idx[1], -1]))

    mmatrix <- rbind(mmatrix, mlpdata)
  }
}


colnames(mmatrix) <- c("AnimalID", "Birthdate", "FirstCalfIndays","Calvingyear","Father", "Lactation","CM" , colnames(gts)[-1])
rownames(mmatrix) <- 1:nrow(mmatrix)

write.table(mmatrix, "CM dataset.txt", sep = "\t", quote= FALSE) ###quote 是否为字符型变量添加双引号