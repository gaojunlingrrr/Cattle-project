library(foreign)
setwd("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/")
gts<- read.table("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep="\t",  na.strings = c("","-" ,"NA","no"), fill = TRUE,header = TRUE,colClasses = "character")
bestand<- read.dbf("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/BESTAND.DBF")
mp<- read.table("/Users/gaojunling/Documents/R Practising/Dedelow/305data-Dedelow.txt", sep = "\t",fill = TRUE, header = TRUE)

# Figure out MKG == NA or MKG == 0, and remove those from the data 
noData <- which(is.na(mp[, "MKG"]) | mp[, "MKG"] == 0 | is.na(mp[, "MKG100"]) | mp[, "MKG100"] ==0 | is.na(mp[, "MKG200_N"]) | mp[, "MKG200_N"]==0)
mp<- mp[-noData,]
mp <- mp[which(mp[, "OHR"] %in% gts[, "ID"]),]
# How many genotyped cows do not have milk production data ?
length(which(!(gts[, "ID"] %in% mp[, "OHR"])))
# Remove the cows without milk production data
gts <- gts[-which(!(gts[, "ID"] %in% mp[, "OHR"])),]
mmatrix <- NULL
for(r in 1:nrow(mp)){
  ohr <- mp[r, "OHR"]
  inBestand <- which(as.character(bestand[, "OHR"]) == as.character(ohr))
  birthdate <- NA
  firstcalf <- NA
  if(length(inBestand) == 1){ 
    birthdate <- as.character(bestand[inBestand, "GEBURT"])
    firstcalf <- as.character(bestand[inBestand, "KALBUNG1"])
    firstcalfindays <- as.character(as.Date(firstcalf, format="%Y-%m-%d") - as.Date(birthdate, format="%Y-%m-%d"))
    fatherid<- as.character(bestand[inBestand, "VATER"])
    
  }
  
  mpdata <- c(as.character(mp[r, "OHR"]), birthdate, firstcalfindays,fatherid,
              as.character(mp[r, "KALBUNG"]),
              as.character(mp[r, "LAKTATION"]), 
              as.character(mp[r, "MKG100"]),as.character(mp[r, "MKG200_N"]),as.character(mp[r, "MKG"]), as.character(mp[r, "FETTKG100"]), as.character(mp[r, "FETTKG200"]),as.character(mp[r, "FETTKG"]),as.character(mp[r, "EIWEISSKG1"]),as.character(mp[r, "EIWEISSKG2"]),as.character(mp[r, "EIWEISSKGG"]))
  idx = which(as.character(gts[, "ID"]) == as.character(mp[r, "OHR"]))
  mpdata <- c(mpdata, as.character(gts[idx[1], -1]))
  
  mmatrix <- rbind(mmatrix, mpdata)
}

colnames(mmatrix) <- c("AnimalID", "Birthdate", "FirstCalfIndays","Father","Calvingdate", "Lactation", "MKG100", "MKG200_N", "MKG", "FATKG100","FATKG200","FATKG","PROTEINKG100", "PROTEINKG200" , "PROTEINKG"  ,colnames(gts)[-1])
rownames(mmatrix) <- 1:nrow(mmatrix)

setwd("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/Rcode update-05-05-2021/")
write.table(mmatrix, "mpdataset.txt", sep = "\t", quote= FALSE) ###quote 是否为字符型变量添加双引号

#####Qulaity control
mmatrix<-read.table("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/Rcode update-05-05-2021/mpdataset.txt", sep = "\t",header = TRUE)
###select animals 1-3
mmatrix<- mmatrix[which(mmatrix[, "Lactation"]<4 & mmatrix[, "Lactation"] > 0), ]
####select animal with birth year >2008
birthyearofcows <- as.numeric(unlist(lapply(strsplit(as.character(mmatrix[,"Birthdate"]), "-"), "[", 1)))
mmatrix<- cbind(mmatrix, birthyearofcows)
mmatrix<- mmatrix[which(mmatrix[, "birthyearofcows"] > 2008), ]
###Delete elements with less than 3 reccords
frequency<- as.data.frame(table(mmatrix[, "AnimalID"]))
frequency2<- frequency[which(frequency[, "Freq"]>2),]
qfanimal<- frequency2[, "Var1"]
mmatrix<- mmatrix[which(mmatrix[, "AnimalID"] %in% qfanimal), ]

####Get season
calvingdates<- as.Date(as.character(mmatrix[, "Calvingdate"]), format="%y/%m/%d")
calvingmonth<- unlist(lapply(strsplit(as.character(calvingdates), "-"), "[", 2))
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
Calvingyear<- unlist(lapply(strsplit(as.character(calvingdates), "-"), "[", 1))
mmatrix1<- cbind(mmatrix, Calvingseason, Calvingyear)

##
setwd("/Users/gaojunling/Desktop/PhD thesis/Manuscript-21.05-Junling/")
productiondata<- write.table(mmatrix1, file= "productiondata.txt", sep = "\t")

# Get the covariates we want to investigate in model building
animal <- mmatrix1[, "AnimalID"]
lactation <- mmatrix1[, "Lactation"]
birthyear <- unlist(lapply(strsplit(as.character(mmatrix1[,"Birthdate"]), "-"), "[", 1))
calvlingyear <- mmatrix1[, "Calvingyear"]
calvingseason<- mmatrix1[, "Calvingseason"]
afc <- mmatrix1[, "FirstCalfIndays"]
father<- mmatrix1[, "Father"]
mkg100<- mmatrix1[, "MKG100"]
mkg200<- mmatrix1[, "MKG200_N"]
mkg305<- mmatrix1[, "MKG"]
fat100<- 100*(mmatrix1[, "FATKG100"]/ mmatrix1[, "MKG100"])
fat200<- 100*(mmatrix1[, "FATKG200"]/ mmatrix1[, "MKG200_N"])
fat305<-100*(mmatrix1[, "FATKG"]/ mmatrix1[, "MKG"])
protein100<- 100*(mmatrix1[, "PROTEINKG100"]/ mmatrix1[, "MKG100"])
protein200<-100* (mmatrix1[, "PROTEINKG200"]/ mmatrix1[, "MKG200_N"])
protein305<- 100*(mmatrix1[, "PROTEINKG"]/ mmatrix1[, "MKG"])

###Set  Reference allele as default
genotypes <- mmatrix1[,17:26]
snp1<-genotypes[, "SaS_M2"]
snp1<-factor(snp1, levels = c("TT", "CT", "CC"))
snp2<- genotypes[, "SaS_M3"]
snp2<- factor(snp2, levels = c("AA", "AG", "GG"))
###
snp3<- genotypes[, "HAS_M6"]
snp3<- factor(snp3, levels = c("CC", "CT", "TT"))
###
snp4<- genotypes[, "HAS_M3"]
snp4<- factor(snp4, levels = c("AA", "AG", "GG"))
###
snp5<- genotypes[, "HAS_M4"]
snp5<- factor(snp5, levels = c("CC", "CT", "TT"))
###
snp6<- genotypes[, "HAS_M8"]
snp6<- factor(snp6, levels = c("GG", "AG", "AA"))
###
snp7<- genotypes[, "HAS_M9"]
snp7<- factor(snp7, levels = c("CC", "CT", "TT"))
###
snp8<- genotypes[, "HAS_M1"]
snp8<- factor(snp8, levels = c("AA", "AG", "GG"))
###
snp9<- genotypes[, "HAS_M7"]
snp9<- factor(snp9, levels = c("TT", "GT", "GG"))
###
snp10<- genotypes[, "HAS_M10"]
snp10<- factor(snp10, levels = c("TT", "CT", "CC"))


# Combine all data into a matrix
mdata <- data.frame(Y = as.numeric(mkg100), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Calvingyear = as.factor(calvlingyear),
                    Calvingseason = as.factor(calvingseason),
                    Afc = as.numeric(afc),
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


##animal selected
mpanimal<- unique(mdata[, "Animal"])
mpsire<- unique(mdata[, "Father"])
#per lactation
m1<- mdata[which(mdata[, "Lact"] == 1), ]
m2<- mdata[which(mdata[, "Lact"] == 2), ]
m3<- mdata[which(mdata[, "Lact"] == 3), ]
mean(m1[, "Y"])
min(m1[, "Y"])
max(m1[, "Y"])
sd(m1[, "Y"])
sd(m1[, "Y"])/sqrt(nrow(m1))

##
mean(m2[, "Y"])
min(m2[, "Y"])
max(m2[, "Y"])
sd(m2[, "Y"])
sd(m2[, "Y"])/sqrt(nrow(m2))

##
mean(m3[, "Y"])
min(m3[, "Y"])
max(m3[, "Y"])
sd(m3[, "Y"])
sd(m3[, "Y"])/sqrt(nrow(m3))

##
mean(mdata[, "Y"])
min(mdata[, "Y"])
max(mdata[, "Y"])
sd(mdata[, "Y"])
sd(mdata[, "Y"])/sqrt(nrow(mdata))



# Which covariates influence the phenotype?
#### 
library(lme4)
null.model <- lmer(Y ~ (1|Animal), data = mdata,REML = FALSE)

model1 <- lmer(Y ~ Lact + (1|Animal), data = mdata,REML = FALSE)
anova(model1, null.model) # < 2.2e-16 ***

model2 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata,REML = FALSE)
anova(model2, null.model) # < 2.2e-16 ***

model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata,REML = FALSE)
anova(model3, null.model) # < 0.03003 *


model4 <- lmer(Y ~ Afc + (1|Animal), data = mdata,REML = FALSE)
anova(model4, null.model) # <  0.3139 NS

model6 <- lmer(Y ~ (1|Father) + (1|Animal), data = mdata, REML = FALSE)
anova(model6, null.model) # <  0.9978


genotypes <- mmatrix1[,17:26]
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata1 <- data.frame(Y = as.numeric(mkg100), 
                       Marker = as.factor(snp5),
                       Lact = as.factor(lactation),
                       Animal = as.factor(animal),
                       
                       Calvingyear = as.factor(calvlingyear),
                       Calvingseason = as.factor(calvingseason),
                       
                       Father= as.factor(father)
  )
}
# Remove rows with missing data
hasMissing <- which(apply(apply(mdata1,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata1 <- mdata1[-hasMissing,] }

#Remove fathers with a small amount of children
mdata1 <- mdata1[which(mdata1[, "Father"] %in% enoughoffspring),]

# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- lmer(Y ~  Lact + Calvingseason +Calvingyear  + (1|Animal) , data = mdata1, REML = FALSE)
markermodel <- lmer(Y ~ Lact  + Calvingseason + Calvingyear  + Marker + (1|Animal), data = mdata1, REML = FALSE)
anova(null.model, markermodel)
markermodel


###Genetic variance
library("MuMIn")

variance.explained.raw <- function(m){
  model.anova = anova(m)
  FE.names = rownames(model.anova)
  RE.var = summary(m)$sigma^2
  FE.var = anova(m)[,"Sum Sq"] / nobs(m)
  names(FE.var) = FE.names
  VarExp.Model = r.squaredGLMM(m)[1, "R2c"]
  VarExp.FE = (VarExp.Model * FE.var) / (RE.var +sum(FE.var))
  VarExp.left = 1 - sum(VarExp.FE)
  VarExp.All = c(VarExp.FE, "Unexplained" = VarExp.left)
  return(VarExp.All)
}

variance.explained.scaled <- function(m){
  variance.explained <- variance.explained.raw(m)
  scaled <- variance.explained / sum(variance.explained)
  return(scaled)
}

100 * (variance.explained.scaled(null.model)["Unexplained"] - variance.explained.scaled(markermodel)["Unexplained"])



##fat100
mdata <- data.frame(Y = as.numeric(fat100), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Calvingyear = as.factor(calvlingyear),
                    Calvingseason = as.factor(calvingseason),
                    Afc = as.numeric(afc),
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


# Which covariates influence the phenotype?
#### 
library(lme4)
null.model <- lmer(Y ~ (1|Animal), data = mdata,REML = FALSE)

model1 <- lmer(Y ~ Lact + (1|Animal), data = mdata,REML = FALSE)
anova(model1, null.model) # < 3.875e-11 ***

model2 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata,REML = FALSE)
anova(model2, null.model) # < 5.067e-08 ***

model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata,REML = FALSE)
anova(model3, null.model) # < 1.226e-08 ***


model4 <- lmer(Y ~ Afc + (1|Animal), data = mdata,REML = FALSE)
anova(model4, null.model) # <  0.9054

model6 <- lmer(Y ~ (1|Father) + (1|Animal), data = mdata, REML = FALSE)
anova(model6, null.model) # < 0.04376 *



genotypes <- mmatrix1[,17:26]
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata1 <- data.frame(Y = as.numeric(fat100), 
                       Marker = as.factor(snp10),
                       Lact = as.factor(lactation),
                       Animal = as.factor(animal),
                       
                       Calvingyear = as.factor(calvlingyear),
                       Calvingseason = as.factor(calvingseason),
                       
                       Father= as.factor(father)
  )
}
# Remove rows with missing data
hasMissing <- which(apply(apply(mdata1,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata1 <- mdata1[-hasMissing,] }

#Remove fathers with a small amount of children
mdata1 <- mdata1[which(mdata1[, "Father"] %in% enoughoffspring),]

# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- lmer(Y ~  Lact + Calvingseason + Calvingyear + Father + (1|Animal) , data = mdata1, REML = FALSE)
markermodel <- lmer(Y ~ Lact  + Calvingseason + Calvingyear + Father + Marker + (1|Animal), data = mdata1, REML = FALSE)
anova(null.model, markermodel)
markermodel


###Genetic variance
library("MuMIn")

variance.explained.raw <- function(m){
  model.anova = anova(m)
  FE.names = rownames(model.anova)
  RE.var = summary(m)$sigma^2
  FE.var = anova(m)[,"Sum Sq"] / nobs(m)
  names(FE.var) = FE.names
  VarExp.Model = r.squaredGLMM(m)[1, "R2c"]
  VarExp.FE = (VarExp.Model * FE.var) / (RE.var +sum(FE.var))
  VarExp.left = 1 - sum(VarExp.FE)
  VarExp.All = c(VarExp.FE, "Unexplained" = VarExp.left)
  return(VarExp.All)
}

variance.explained.scaled <- function(m){
  variance.explained <- variance.explained.raw(m)
  scaled <- variance.explained / sum(variance.explained)
  return(scaled)
}

100 * (variance.explained.scaled(null.model)["Unexplained"] - variance.explained.scaled(markermodel)["Unexplained"])


##protein100
mdata <- data.frame(Y = as.numeric(protein100), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Calvingyear = as.factor(calvlingyear),
                    Calvingseason = as.factor(calvingseason),
                    Afc = as.numeric(afc),
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


# Which covariates influence the phenotype?
#### 
library(lme4)
null.model <- lmer(Y ~ (1|Animal), data = mdata,REML = FALSE)

model1 <- lmer(Y ~ Lact + (1|Animal), data = mdata,REML = FALSE)
anova(model1, null.model) # < 0.07389 

model2 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata,REML = FALSE)
anova(model2, null.model) # < 0.003506 **

model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata,REML = FALSE)
anova(model3, null.model) # < 2.2e-16 ***


model4 <- lmer(Y ~ Afc + (1|Animal), data = mdata,REML = FALSE)
anova(model4, null.model) # <  0.2977

model6 <- lmer(Y ~ (1|Father) + (1|Animal), data = mdata, REML = FALSE)
anova(model6, null.model) # < 1.939e-06 ***



genotypes <- mmatrix1[,17:26]
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata1 <- data.frame(Y = as.numeric(protein100), 
                       Marker = as.factor(snp10),
                       Lact = as.factor(lactation),
                       Animal = as.factor(animal),
                       
                       Calvingyear = as.factor(calvlingyear),
                       Calvingseason = as.factor(calvingseason),
                       
                       Father= as.factor(father)
  )
}
# Remove rows with missing data
hasMissing <- which(apply(apply(mdata1,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata1 <- mdata1[-hasMissing,] }

#Remove fathers with a small amount of children
mdata1 <- mdata1[which(mdata1[, "Father"] %in% enoughoffspring),]

# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- lmer(Y ~  Lact + Calvingseason + Calvingyear + Father + (1|Animal) , data = mdata1, REML = FALSE)
markermodel <- lmer(Y ~ Lact  + Calvingseason + Calvingyear + Father + Marker + (1|Animal), data = mdata1, REML = FALSE)
anova(null.model, markermodel)
markermodel


###Genetic variance
library("MuMIn")

variance.explained.raw <- function(m){
  model.anova = anova(m)
  FE.names = rownames(model.anova)
  RE.var = summary(m)$sigma^2
  FE.var = anova(m)[,"Sum Sq"] / nobs(m)
  names(FE.var) = FE.names
  VarExp.Model = r.squaredGLMM(m)[1, "R2c"]
  VarExp.FE = (VarExp.Model * FE.var) / (RE.var +sum(FE.var))
  VarExp.left = 1 - sum(VarExp.FE)
  VarExp.All = c(VarExp.FE, "Unexplained" = VarExp.left)
  return(VarExp.All)
}

variance.explained.scaled <- function(m){
  variance.explained <- variance.explained.raw(m)
  scaled <- variance.explained / sum(variance.explained)
  return(scaled)
}

100 * (variance.explained.scaled(null.model)["Unexplained"] - variance.explained.scaled(markermodel)["Unexplained"])



##mkg200
mdata <- data.frame(Y = as.numeric(mkg200), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Calvingyear = as.factor(calvlingyear),
                    Calvingseason = as.factor(calvingseason),
                    Afc = as.numeric(afc),
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


# Which covariates influence the phenotype?
#### 
library(lme4)
null.model <- lmer(Y ~ (1|Animal), data = mdata,REML = FALSE)

model1 <- lmer(Y ~ Lact + (1|Animal), data = mdata,REML = FALSE)
anova(model1, null.model) # < 2.2e-16 ***

model2 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata,REML = FALSE)
anova(model2, null.model) # < 2.2e-16 ***

model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata,REML = FALSE)
anova(model3, null.model) # < 0.005113 **


model4 <- lmer(Y ~ Afc + (1|Animal), data = mdata,REML = FALSE)
anova(model4, null.model) # <  0.3283

model6 <- lmer(Y ~ (1|Father) + (1|Animal), data = mdata, REML = FALSE)
anova(model6, null.model) # < 0.3646




genotypes <- mmatrix1[,17:26]
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata1 <- data.frame(Y = as.numeric(mkg200), 
                       Marker = as.factor(snp5),
                       Lact = as.factor(lactation),
                       Animal = as.factor(animal),
                       
                       Calvingyear = as.factor(calvlingyear),
                       Calvingseason = as.factor(calvingseason),
                       
                       Father= as.factor(father)
  )
}
# Remove rows with missing data
hasMissing <- which(apply(apply(mdata1,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata1 <- mdata1[-hasMissing,] }

#Remove fathers with a small amount of children
mdata1 <- mdata1[which(mdata1[, "Father"] %in% enoughoffspring),]

# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- lmer(Y ~  Lact + Calvingseason + Calvingyear + (1|Animal) , data = mdata1, REML = FALSE)
markermodel <- lmer(Y ~ Lact  + Calvingseason + Calvingyear  + Marker + (1|Animal), data = mdata1, REML = FALSE)
anova(null.model, markermodel)
markermodel


###Genetic variance
library("MuMIn")

variance.explained.raw <- function(m){
  model.anova = anova(m)
  FE.names = rownames(model.anova)
  RE.var = summary(m)$sigma^2
  FE.var = anova(m)[,"Sum Sq"] / nobs(m)
  names(FE.var) = FE.names
  VarExp.Model = r.squaredGLMM(m)[1, "R2c"]
  VarExp.FE = (VarExp.Model * FE.var) / (RE.var +sum(FE.var))
  VarExp.left = 1 - sum(VarExp.FE)
  VarExp.All = c(VarExp.FE, "Unexplained" = VarExp.left)
  return(VarExp.All)
}

variance.explained.scaled <- function(m){
  variance.explained <- variance.explained.raw(m)
  scaled <- variance.explained / sum(variance.explained)
  return(scaled)
}

100 * (variance.explained.scaled(null.model)["Unexplained"] - variance.explained.scaled(markermodel)["Unexplained"])


##fat200
mdata <- data.frame(Y = as.numeric(fat200), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Calvingyear = as.factor(calvlingyear),
                    Calvingseason = as.factor(calvingseason),
                    Afc = as.numeric(afc),
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


# Which covariates influence the phenotype?
#### 
library(lme4)
null.model <- lmer(Y ~ (1|Animal), data = mdata,REML = FALSE)

model1 <- lmer(Y ~ Lact + (1|Animal), data = mdata,REML = FALSE)
anova(model1, null.model) # < 2.2e-16 ***

model2 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata,REML = FALSE)
anova(model2, null.model) # < 2.2e-16 ***

model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata,REML = FALSE)
anova(model3, null.model) # < 2.86e-07 ***


model4 <- lmer(Y ~ Afc + (1|Animal), data = mdata,REML = FALSE)
anova(model4, null.model) # <  0.9709

model6 <- lmer(Y ~ (1|Father) + (1|Animal), data = mdata, REML = FALSE)
anova(model6, null.model) # < 0.006544 **




genotypes <- mmatrix1[,17:26]
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata1 <- data.frame(Y = as.numeric(fat200), 
                       Marker = as.factor(snp10),
                       Lact = as.factor(lactation),
                       Animal = as.factor(animal),
                       
                       Calvingyear = as.factor(calvlingyear),
                       Calvingseason = as.factor(calvingseason),
                       
                       Father= as.factor(father)
  )
}
# Remove rows with missing data
hasMissing <- which(apply(apply(mdata1,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata1 <- mdata1[-hasMissing,] }

#Remove fathers with a small amount of children
mdata1 <- mdata1[which(mdata1[, "Father"] %in% enoughoffspring),]

# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- lmer(Y ~  Lact + Calvingseason + Calvingyear + Father + (1|Animal) , data = mdata1, REML = FALSE)
markermodel <- lmer(Y ~ Lact  + Calvingseason + Calvingyear + Father + Marker + (1|Animal), data = mdata1, REML = FALSE)
anova(null.model, markermodel)
markermodel


###Genetic variance
library("MuMIn")

variance.explained.raw <- function(m){
  model.anova = anova(m)
  FE.names = rownames(model.anova)
  RE.var = summary(m)$sigma^2
  FE.var = anova(m)[,"Sum Sq"] / nobs(m)
  names(FE.var) = FE.names
  VarExp.Model = r.squaredGLMM(m)[1, "R2c"]
  VarExp.FE = (VarExp.Model * FE.var) / (RE.var +sum(FE.var))
  VarExp.left = 1 - sum(VarExp.FE)
  VarExp.All = c(VarExp.FE, "Unexplained" = VarExp.left)
  return(VarExp.All)
}

variance.explained.scaled <- function(m){
  variance.explained <- variance.explained.raw(m)
  scaled <- variance.explained / sum(variance.explained)
  return(scaled)
}

100 * (variance.explained.scaled(null.model)["Unexplained"] - variance.explained.scaled(markermodel)["Unexplained"])



##protein200
mdata <- data.frame(Y = as.numeric(protein200), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Calvingyear = as.factor(calvlingyear),
                    Calvingseason = as.factor(calvingseason),
                    Afc = as.numeric(afc),
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


# Which covariates influence the phenotype?
#### 
library(lme4)
null.model <- lmer(Y ~ (1|Animal), data = mdata,REML = FALSE)

model1 <- lmer(Y ~ Lact + (1|Animal), data = mdata,REML = FALSE)
anova(model1, null.model) # < 3.438e-05 ***

model2 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata,REML = FALSE)
anova(model2, null.model) # < 3.933e-07 ***

model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata,REML = FALSE)
anova(model3, null.model) # < 1.898e-11 ***


model4 <- lmer(Y ~ Afc + (1|Animal), data = mdata,REML = FALSE)
anova(model4, null.model) # <  0.3729

model6 <- lmer(Y ~ (1|Father) + (1|Animal), data = mdata, REML = FALSE)
anova(model6, null.model) # < 3.756e-08 ***




genotypes <- mmatrix1[,17:26]
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata1 <- data.frame(Y = as.numeric(protein200), 
                       Marker = as.factor(snp10),
                       Lact = as.factor(lactation),
                       Animal = as.factor(animal),
                       
                       Calvingyear = as.factor(calvlingyear),
                       Calvingseason = as.factor(calvingseason),
                       
                       Father= as.factor(father)
  )
}
# Remove rows with missing data
hasMissing <- which(apply(apply(mdata1,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata1 <- mdata1[-hasMissing,] }

#Remove fathers with a small amount of children
mdata1 <- mdata1[which(mdata1[, "Father"] %in% enoughoffspring),]

# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- lmer(Y ~  Lact + Calvingseason + Calvingyear + Father + (1|Animal) , data = mdata1, REML = FALSE)
markermodel <- lmer(Y ~ Lact  + Calvingseason + Calvingyear + Father + Marker + (1|Animal), data = mdata1, REML = FALSE)
anova(null.model, markermodel)
markermodel


###Genetic variance
library("MuMIn")

variance.explained.raw <- function(m){
  model.anova = anova(m)
  FE.names = rownames(model.anova)
  RE.var = summary(m)$sigma^2
  FE.var = anova(m)[,"Sum Sq"] / nobs(m)
  names(FE.var) = FE.names
  VarExp.Model = r.squaredGLMM(m)[1, "R2c"]
  VarExp.FE = (VarExp.Model * FE.var) / (RE.var +sum(FE.var))
  VarExp.left = 1 - sum(VarExp.FE)
  VarExp.All = c(VarExp.FE, "Unexplained" = VarExp.left)
  return(VarExp.All)
}

variance.explained.scaled <- function(m){
  variance.explained <- variance.explained.raw(m)
  scaled <- variance.explained / sum(variance.explained)
  return(scaled)
}

100 * (variance.explained.scaled(null.model)["Unexplained"] - variance.explained.scaled(markermodel)["Unexplained"])



##mkg305
mdata <- data.frame(Y = as.numeric(mkg305), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Calvingyear = as.factor(calvlingyear),
                    Calvingseason = as.factor(calvingseason),
                    Afc = as.numeric(afc),
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


# Which covariates influence the phenotype?
#### 
library(lme4)
null.model <- lmer(Y ~ (1|Animal), data = mdata,REML = FALSE)

model1 <- lmer(Y ~ Lact + (1|Animal), data = mdata,REML = FALSE)
anova(model1, null.model) # < 2.2e-16 ***

model2 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata,REML = FALSE)
anova(model2, null.model) # < 2.2e-16 ***

model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata,REML = FALSE)
anova(model3, null.model) # < 0.1374


model4 <- lmer(Y ~ Afc + (1|Animal), data = mdata,REML = FALSE)
anova(model4, null.model) # <  0.3102

model6 <- lmer(Y ~ (1|Father) + (1|Animal), data = mdata, REML = FALSE)
anova(model6, null.model) # < 0.1999




genotypes <- mmatrix1[,17:26]
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata1 <- data.frame(Y = as.numeric(mkg305), 
                       Marker = as.factor(snp5),
                       Lact = as.factor(lactation),
                       Animal = as.factor(animal),
                       
                       Calvingyear = as.factor(calvlingyear),
                       Calvingseason = as.factor(calvingseason),
                       
                       Father= as.factor(father)
  )
}
# Remove rows with missing data
hasMissing <- which(apply(apply(mdata1,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata1 <- mdata1[-hasMissing,] }

#Remove fathers with a small amount of children
mdata1 <- mdata1[which(mdata1[, "Father"] %in% enoughoffspring),]

# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- lmer(Y ~  Lact  + Calvingyear  + (1|Animal) , data = mdata1, REML = FALSE)
markermodel <- lmer(Y ~ Lact  + Calvingyear  + Marker + (1|Animal), data = mdata1, REML = FALSE)
anova(null.model, markermodel)
markermodel


###Genetic variance
library("MuMIn")

variance.explained.raw <- function(m){
  model.anova = anova(m)
  FE.names = rownames(model.anova)
  RE.var = summary(m)$sigma^2
  FE.var = anova(m)[,"Sum Sq"] / nobs(m)
  names(FE.var) = FE.names
  VarExp.Model = r.squaredGLMM(m)[1, "R2c"]
  VarExp.FE = (VarExp.Model * FE.var) / (RE.var +sum(FE.var))
  VarExp.left = 1 - sum(VarExp.FE)
  VarExp.All = c(VarExp.FE, "Unexplained" = VarExp.left)
  return(VarExp.All)
}

variance.explained.scaled <- function(m){
  variance.explained <- variance.explained.raw(m)
  scaled <- variance.explained / sum(variance.explained)
  return(scaled)
}

100 * (variance.explained.scaled(null.model)["Unexplained"] - variance.explained.scaled(markermodel)["Unexplained"])



##fat305
mdata <- data.frame(Y = as.numeric(fat305), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Calvingyear = as.factor(calvlingyear),
                    Calvingseason = as.factor(calvingseason),
                    Afc = as.numeric(afc),
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


# Which covariates influence the phenotype?
#### 
library(lme4)
null.model <- lmer(Y ~ (1|Animal), data = mdata,REML = FALSE)

model1 <- lmer(Y ~ Lact + (1|Animal), data = mdata,REML = FALSE)
anova(model1, null.model) # < 2.2e-16 ***

model2 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata,REML = FALSE)
anova(model2, null.model) # < 2.2e-16 ***

model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata,REML = FALSE)
anova(model3, null.model) # < 1.364e-05 ***


model4 <- lmer(Y ~ Afc + (1|Animal), data = mdata,REML = FALSE)
anova(model4, null.model) # <  0.8464

model6 <- lmer(Y ~ (1|Father) + (1|Animal), data = mdata, REML = FALSE)
anova(model6, null.model) # < 0.002602 **




genotypes <- mmatrix1[,17:26]
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata1 <- data.frame(Y = as.numeric(fat305), 
                       Marker = as.factor(snp10),
                       Lact = as.factor(lactation),
                       Animal = as.factor(animal),
                       
                       Calvingyear = as.factor(calvlingyear),
                       Calvingseason = as.factor(calvingseason),
                       
                       Father= as.factor(father)
  )
}
# Remove rows with missing data
hasMissing <- which(apply(apply(mdata1,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata1 <- mdata1[-hasMissing,] }

#Remove fathers with a small amount of children
mdata1 <- mdata1[which(mdata1[, "Father"] %in% enoughoffspring),]

# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- lmer(Y ~  Lact  + Calvingyear + Calvingseason + Father + (1|Animal) , data = mdata1, REML = FALSE)
markermodel <- lmer(Y ~ Lact  + Calvingyear + Calvingseason + Father + Marker + (1|Animal), data = mdata1, REML = FALSE)
anova(null.model, markermodel)
markermodel


###Genetic variance
library("MuMIn")

variance.explained.raw <- function(m){
  model.anova = anova(m)
  FE.names = rownames(model.anova)
  RE.var = summary(m)$sigma^2
  FE.var = anova(m)[,"Sum Sq"] / nobs(m)
  names(FE.var) = FE.names
  VarExp.Model = r.squaredGLMM(m)[1, "R2c"]
  VarExp.FE = (VarExp.Model * FE.var) / (RE.var +sum(FE.var))
  VarExp.left = 1 - sum(VarExp.FE)
  VarExp.All = c(VarExp.FE, "Unexplained" = VarExp.left)
  return(VarExp.All)
}

variance.explained.scaled <- function(m){
  variance.explained <- variance.explained.raw(m)
  scaled <- variance.explained / sum(variance.explained)
  return(scaled)
}

100 * (variance.explained.scaled(null.model)["Unexplained"] - variance.explained.scaled(markermodel)["Unexplained"])



##protein305
mdata <- data.frame(Y = as.numeric(protein305), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Calvingyear = as.factor(calvlingyear),
                    Calvingseason = as.factor(calvingseason),
                    Afc = as.numeric(afc),
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


# Which covariates influence the phenotype?
#### 
library(lme4)
null.model <- lmer(Y ~ (1|Animal), data = mdata,REML = FALSE)

model1 <- lmer(Y ~ Lact + (1|Animal), data = mdata,REML = FALSE)
anova(model1, null.model) # < 0.03396 *

model2 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata,REML = FALSE)
anova(model2, null.model) # <  0.49

model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata,REML = FALSE)
anova(model3, null.model) # <  0.02398 *


model4 <- lmer(Y ~ Afc + (1|Animal), data = mdata,REML = FALSE)
anova(model4, null.model) # <   0.3258

model6 <- lmer(Y ~ (1|Father) + (1|Animal), data = mdata, REML = FALSE)
anova(model6, null.model) # < 4.861e-07 ***




genotypes <- mmatrix1[,17:26]
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata1 <- data.frame(Y = as.numeric(protein305), 
                       Marker = as.factor(snp8),
                       Lact = as.factor(lactation),
                       Animal = as.factor(animal),
                       
                       Calvingyear = as.factor(calvlingyear),
                       Calvingseason = as.factor(calvingseason),
                       
                       Father= as.factor(father)
  )
}
# Remove rows with missing data
hasMissing <- which(apply(apply(mdata1,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata1 <- mdata1[-hasMissing,] }

#Remove fathers with a small amount of children
mdata1 <- mdata1[which(mdata1[, "Father"] %in% enoughoffspring),]

# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- lmer(Y ~  Calvingseason  + Lact + Father + (1|Animal) , data = mdata1, REML = FALSE)
markermodel <- lmer(Y ~ Calvingseason  + Lact + Father + Marker + (1|Animal), data = mdata1, REML = FALSE)
anova(null.model, markermodel)
markermodel


###Genetic variance
library("MuMIn")

variance.explained.raw <- function(m){
  model.anova = anova(m)
  FE.names = rownames(model.anova)
  RE.var = summary(m)$sigma^2
  FE.var = anova(m)[,"Sum Sq"] / nobs(m)
  names(FE.var) = FE.names
  VarExp.Model = r.squaredGLMM(m)[1, "R2c"]
  VarExp.FE = (VarExp.Model * FE.var) / (RE.var +sum(FE.var))
  VarExp.left = 1 - sum(VarExp.FE)
  VarExp.All = c(VarExp.FE, "Unexplained" = VarExp.left)
  return(VarExp.All)
}

variance.explained.scaled <- function(m){
  variance.explained <- variance.explained.raw(m)
  scaled <- variance.explained / sum(variance.explained)
  return(scaled)
}

100 * (variance.explained.scaled(null.model)["Unexplained"] - variance.explained.scaled(markermodel)["Unexplained"])



#### find closest gene

library(biomaRt)


bio.mart <- useMart(biomart="ensembl", dataset="btaurus_gene_ensembl")


res.biomart <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "chromosome_name", "start_position", "end_position", "strand"), 
                     filters = c("chromosomal_region", "biotype"), 
                     values = list(c("6:83553915:84053915"),"protein_coding"), mart = bio.mart)

res.biomart



####P-VALUE correction
p.adjust(c(0.078,
           0.7574,
           0.0644,
           0.3887,
           0.4793,
           0.6212,
           0.9424,
           0.5763,
           0.9314,
           0.7458), method = "BH")

