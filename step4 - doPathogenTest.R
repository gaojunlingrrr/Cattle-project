#### Script4.R: Pathogen cases with given genotypes
setwd("C:/Users/Arends/Downloads/R.chisquare_01.11.2021_Junling/R.chisquare_01.11.2021_Junling/")
#load raw data of bacteria
bacteria<- read.table("Pathogen data of Dedelow.txt", sep = "\t",header = TRUE)

#load genotrans data & CM data to add the CM code to the genotrans matrix
genotrans<- read.table("genotrans.txt", sep = "\t",header = TRUE)
cmdata<- read.table("CM dataset.txt", sep = "\t",header = TRUE)
genotrans <- cbind(genotrans, CM = NA)
for(animal in rownames(genotrans)){
  CMcode <- unique(cmdata[which(cmdata[, "AnimalID"] ==  animal), "CM"])
  if(length(CMcode) == 1){
    genotrans[animal, "CM"] <- CMcode
  }
}

#load herd data
library(foreign)
bestand<- read.dbf("BESTAND.DBF")

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
}
colnames(pathogendata)<- c("AnimalID", "STALL", "Farm", "Bacteria", "CNS", "SC+", "SC-", "STA", "E.coli", "GALT")
rownames(pathogendata)<- 1:nrow(pathogendata)


pdata <- matrix(NA, length(unique(pathogendata[, "AnimalID"])), 18, dimnames = list(unique(pathogendata[, "AnimalID"]), c("Bacteria", "CNS", "SC+", "SC-", "STA", "E.coli", "GALT", colnames(genotrans))))

# Fill in the genotype data (if any)
for(animal in rownames(pdata)){
  pdata[animal, colnames(genotrans)] <- as.character(genotrans[animal, ])
}
# Remove animals not genotyped
pdata <- pdata[-which(pdata[,"CM"] == "NA"),]

# Fill in the bacterial data
for(animal in rownames(pdata)){
  iix <- which(pathogendata[, "AnimalID"] == animal)
  for(i in iix){
    for(colu in c("Bacteria", "CNS", "SC+", "SC-", "STA", "E.coli", "GALT")){
      if(!is.na(pathogendata[i, colu])) pdata[animal, colu] <- pathogendata[i, colu]
    }
  }
}

pdata <-  pdata[which(pdata[, "CM"] == "1"),]

markers <- colnames(pdata)[8:17]

mmatrix <- array(NA, dim=c(7, 10, 3), dimnames=list(c("Bacteria", "CNS", "SC+", "SC-", "STA", "E.coli", "GALT"), markers, c("AA", "AB", "BB")))

### Fill in the matrix
for(x in c("Bacteria", "CNS", "SC+", "SC-", "STA", "E.coli", "GALT")){
  for(y in markers){
    N <- length(which(pdata[,x] == "1"))
    mmatrix[x,y,"AA"] <- length(which(pdata[,y] == "AA" & pdata[,x] == "1"))
    mmatrix[x,y,"AB"] <- length(which(pdata[,y] == "AB" & pdata[,x] == "1"))
    mmatrix[x,y,"BB"] <- length(which(pdata[,y] == "BB" & pdata[,x] == "1"))
  }
}

sickfreq<- read.table("genotypeFrequenciesSickCows.txt", sep = "\t", header = TRUE)
rownames(sickfreq) <- lapply(strsplit(rownames(sickfreq), ":"), "[", 1)

mres <- c()
for(x in c("Bacteria", "CNS", "SC+", "SC-", "STA", "E.coli", "GALT")){
  for(y in markers){
    p <- chisq.test(c(mmatrix[x, y, "AA"], mmatrix[x, y, "AB"],mmatrix[x, y, "BB"]), p=as.numeric(sickfreq[y,]))$p.value
    cat(x, " ", y, " ", p, "\n")
    mres <- rbind(mres, c(x, y, p, NA, c(mmatrix[x, y, "AA"], mmatrix[x, y, "AB"],mmatrix[x, y, "BB"]), as.numeric(sickfreq[y,])))
  }
}

colnames(mres) <- c("Pathogen", "Marker", "pvalue", "pAdjusted", "nObs_AA", "nObs_AB", "nObs_BB", "freq_AA_sick", "freq_AB_sick", "freq_BB_sick")

# Padjust per bacteria (use BH on 10 pvalues)
for(x in c("Bacteria", "CNS", "SC+", "SC-", "STA", "E.coli", "GALT")){
  iix <- which(mres[, "Pathogen"] == x)
  mres[iix, "pAdjusted"] <- p.adjust(as.numeric(mres[iix, "pvalue"]), method = "BH")
}

write.table(mres, "pvalues_pathogens.txt", sep = "\t", quote = FALSE, row.names = FALSE)
