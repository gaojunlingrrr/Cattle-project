mysample<- read.table("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep = "\t", na.strings = c("","-" ,"NA","no"), fill = TRUE, header = TRUE)
mysample <- mysample[, -c(13,14)]
myanimals<- mysample[,"ID"]
setwd("/Users/gaojunling/Documents/R Practising/Dedelow/")
mydata <- read.table("krankData-Dedelow.txt", sep = "\t",fill = TRUE, header = TRUE)
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


###
library(foreign)
bestand<- read.dbf("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/BESTAND.DBF")
gts<- read.table("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep="\t",  na.strings = c("","-" ,"NA","no"), fill = TRUE,header = TRUE,colClasses = "character")
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


setwd("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/Rcode update-05-05-2021/")
write.table(mmatrix, "CM dataset.txt", sep = "\t", quote= FALSE) ###quote 是否为字符型变量添加双引号

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
### select lactation 1-3
cmdata<- cmdata[which(cmdata[, "Lactation"] < 4 & cmdata[, "Lactation"] > 0), ]


####select animal with birth year >2008
birthyear <- as.numeric(unlist(lapply(strsplit(as.character(cmdata[,"Birthdate"]), "-"), "[", 1)))
cmdata<- cbind(cmdata, birthyear)
cmdata<- cmdata[which(cmdata[, "birthyear"] > 2008), ]

###Delete elements with less than 3 reccords
frequency<- as.data.frame(table(cmdata[, "AnimalID"]))
frequency2<- frequency[which(frequency[, "Freq"]>2),]
qfanimal<- frequency2[, "Var1"]
cmdata<- cmdata[which(cmdata[, "AnimalID"] %in% qfanimal), ]


##get calving year
Calvingyear_N<- unlist(lapply(strsplit(as.character(cmdata[,"Calvingyear"]), "-"), "[", 1))

cmdata<- cbind(cmdata, Calvingyear_N)

##
setwd("/Users/gaojunling/Desktop/PhD thesis/Manuscript-21.05-Junling/")
CMDATA<- write.table(cmdata, file= "CMDATA.txt", sep = "\t")

# Get the covariates we want to investigate in model building
animal <- cmdata[, "AnimalID"]
lactation <- cmdata[, "Lactation"]
seasoncalving<- cmdata[, "Calvingseason"]
yearcalving <-  cmdata[, "Calvingyear_N"]
firstCalf <- cmdata[, "FirstCalfIndays"]
father<- cmdata[, "Father"]
animalcode<- cmdata[, "CM"] 


##combine all data into a matrix
###CM

mdata1 <- data.frame(Y = as.numeric(animalcode), 
                     Lact = as.factor(lactation),
                     Animal = as.factor(animal),
                     Calvingseason = as.factor(seasoncalving),
                     Calvingyear = as.factor(yearcalving),
                     Firstcalf = as.numeric(firstCalf),
                     Father= as.factor(father)
)



# Remove rows with missing data
mdata1<- na.omit(mdata1)

#Improve the association by removing fathers which do not have 10 offspring
offspring <- c()
for(f in unique(mdata1[, "Father"])){
  nOffspring = length(unique(mdata1[which(mdata1[, "Father"] == f), "Animal"]))
  offspring <- c(offspring, nOffspring)
}
names(offspring) <- unique(mdata1[, "Father"])

enoughoffspring <- names(offspring)[which(offspring >= 10)]
dim(mdata1)
mdata1 <- mdata1[which(mdata1[, "Father"] %in% enoughoffspring),]
dim(mdata1)

###number of cows and observations per lactation
Lact1<- mdata1[which(mdata1[, "Lact"] == 1), ]
Lact2<- mdata1[which(mdata1[, "Lact"] == 2), ]
Lact3<- mdata1[which(mdata1[, "Lact"] == 3), ]
length(unique(Lact1[, "Animal"]))#1197
length(unique(Lact2[, "Animal"]))#840
length(unique(Lact3[, "Animal"]))#448
nrow(Lact1)
nrow(Lact2)
nrow(Lact3)

###Maastitis incidence per lactation
###number of sick and healthy animal
healthygroup<- mdata1[which(mdata1[, "Y"] == 0), ]
sickgroup<- mdata1[which(mdata1[, "Y"] == 1), ]
length(unique(healthygroup[, "Animal"])) ##644
length(unique(sickgroup[, "Animal"])) ##941
healthylact1 = healthygroup[which(healthygroup[, "Lact"]==1), ]
healthylact2 = healthygroup[which(healthygroup[, "Lact"]==2), ]
healthylact3 = healthygroup[which(healthygroup[, "Lact"]==3), ]
sicklact1 = sickgroup[which(sickgroup[, "Lact"]==1), ]
sicklact2 = sickgroup[which(sickgroup[, "Lact"]==2), ]
sicklact3 = sickgroup[which(sickgroup[, "Lact"]==3), ]
length(unique(healthylact1[, "Animal"]))#644
length(unique(healthylact2[, "Animal"]))#288
length(unique(healthylact3[, "Animal"]))
length(unique(sicklact1[, "Animal"]))
length(unique(sicklact2[, "Animal"]))
length(unique(sicklact3[, "Animal"]))


###mastitis incidence per lactation
length(unique(sicklact1[, "Animal"]))/length(unique(Lact1[, "Animal"]))#34.9%
length(unique(sicklact2[, "Animal"]))/length(unique(Lact2[, "Animal"]))#34.8%
length(unique(sicklact3[, "Animal"]))/length(unique(Lact3[, "Animal"]))#23.1%
length(unique(healthygroup[, "Animal"]))/length(unique(mdata1[, "Animal"]))#40.6%
###Animal selected
cmanimal<- unique(mdata1[, "Animal"])
cmsire<- unique(mdata1[, "Father"])

# Which covariates influence the phenotype?
library(lme4)
null.model <- glmer(Y ~ (1|Animal), data = mdata1, family = binomial,nAGQ = 0)

model1 <- glmer(Y ~ Lact + (1|Animal), data = mdata1, family = binomial, nAGQ = 0)
anova(model1, null.model, test="Chisq") # < 2.2e-16 ***

model2 <- glmer(Y ~ Calvingseason + (1|Animal), data = mdata1,family = binomial, nAGQ = 0)　　　
anova(model2, null.model, test="Chisq") # < 4.833e-05 ***

model3 <- glmer(Y ~ Calvingyear + (1|Animal), data = mdata1, family = binomial, nAGQ = 0)
anova(model3, null.model, test="Chisq") # < 2.2e-16 ***

model4 <- glmer(Y ~ Firstcalf + (1|Animal), data = mdata1, family = binomial, nAGQ = 0)
anova(model4, null.model, test="Chisq") # < 3.102e-12 *** 

model5 <- glmer(Y ~ Father + (1|Animal), data = mdata1, family = binomial, nAGQ = 0)
anova(model5, null.model, test="Chisq") # < 2.2e-16 ***

###genotypes
genotypes <- cmdata[,9:18]
snp1<-genotypes[, "SaS_M2"]
snp1<-factor(snp1, levels = c("TT", "CT", "CC"))
####
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






# Association analysis
mdata2 <- data.frame(Y = as.numeric(animalcode), 
                     Marker = as.factor(snp6),
                     Animal = as.factor(animal),
                     Lact = as.factor(lactation),
                     Calvingseason = as.factor(seasoncalving),                 
                     Calvingyear = as.factor(yearcalving),
                     Firstcalf = as.numeric(firstCalf),
                     Father= as.factor(father))

# Remove rows with missing data
hasMissing <- which(apply(apply(mdata2,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata2 <- mdata2[-hasMissing,] }

#Remove fathers with a small amount of children
mdata2 <- mdata2[which(mdata2[, "Father"] %in% enoughoffspring),]

# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- glmer(Y ~ Lact + Calvingseason +Calvingyear + Firstcalf + Father + (1|Animal), data = mdata2, family = binomial, nAGQ=0)
markermodel <- glmer(Y ~ Lact + Calvingseason +Calvingyear + Firstcalf + Father+ Marker+ (1|Animal) , data = mdata2, family = binomial,nAGQ=0)  
anova(null.model, markermodel)
markermodel



##Genetic Variance
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


###Another way to calculate variance explained 
(r.squaredGLMM(markermodel)[1, "R2c"]-r.squaredGLMM(null.model)[1, "R2c"])*100




####Mastitis incidence of Ref/Alt per SNP
###chr5
chrdata <- data.frame(Y = as.numeric(animalcode), 
                     Marker = as.factor(snp8),
                     Animal = as.factor(animal),
                     Father= as.factor(father))

# Remove rows with missing data
hasMissing <- which(apply(apply(chrdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { chrdata <- chrdata[-hasMissing,] }

#Remove fathers with a small amount of children
chrdata <- chrdata[which(chrdata[, "Father"] %in% enoughoffspring),]
AA<- chrdata[which(chrdata[, "Marker"]== "AA"), ]
AB<- chrdata[which(chrdata[, "Marker"]== "AG"), ]
BB<- chrdata[which(chrdata[, "Marker"]== "GG"), ]
sickAA<-AA[which(AA[, "Y"]==1),]
sickAB<-AA[which(AB[, "Y"]==1),]
sickBB<-AA[which(BB[, "Y"]==1),]

length(unique(AA[, "Animal"]))#800
length(unique(AB[, "Animal"]))#650
length(unique(BB[, "Animal"]))#135

length(unique(sickAA[, "Animal"]))/length(unique(chrdata[, "Animal"]))#29.6%
length(unique(sickAB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#25.6%
length(unique(sickBB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#9.2%


###chr6-1
chrdata <- data.frame(Y = as.numeric(animalcode), 
                      Marker = as.factor(snp1),
                      Animal = as.factor(animal),
                      Father= as.factor(father))

# Remove rows with missing data
hasMissing <- which(apply(apply(chrdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { chrdata <- chrdata[-hasMissing,] }

#Remove fathers with a small amount of children
chrdata <- chrdata[which(chrdata[, "Father"] %in% enoughoffspring),]
AA<- chrdata[which(chrdata[, "Marker"]== "TT"), ]
AB<- chrdata[which(chrdata[, "Marker"]== "CT"), ]
BB<- chrdata[which(chrdata[, "Marker"]== "CC"), ]
sickAA<-AA[which(AA[, "Y"]==1),]
sickAB<-AA[which(AB[, "Y"]==1),]
sickBB<-AA[which(BB[, "Y"]==1),]

length(unique(AA[, "Animal"]))#606
length(unique(AB[, "Animal"]))#730
length(unique(BB[, "Animal"]))#249

length(unique(sickAA[, "Animal"]))/length(unique(chrdata[, "Animal"]))#22.0%
length(unique(sickAB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#34.3%
length(unique(sickBB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#12.4%


###chr6-2
chrdata <- data.frame(Y = as.numeric(animalcode), 
                      Marker = as.factor(snp2),
                      Animal = as.factor(animal),
                      Father= as.factor(father))

# Remove rows with missing data
hasMissing <- which(apply(apply(chrdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { chrdata <- chrdata[-hasMissing,] }

#Remove fathers with a small amount of children
chrdata <- chrdata[which(chrdata[, "Father"] %in% enoughoffspring),]
AA<- chrdata[which(chrdata[, "Marker"]== "AA"), ]
AB<- chrdata[which(chrdata[, "Marker"]== "AG"), ]
BB<- chrdata[which(chrdata[, "Marker"]== "GG"), ]
sickAA<-AA[which(AA[, "Y"]==1),]
sickAB<-AA[which(AB[, "Y"]==1),]
sickBB<-AA[which(BB[, "Y"]==1),]

length(unique(AA[, "Animal"]))#896
length(unique(AB[, "Animal"]))#730
length(unique(BB[, "Animal"]))#88

length(unique(sickAA[, "Animal"]))/length(unique(chrdata[, "Animal"]))#32.1%
length(unique(sickAB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#27.8%
length(unique(sickBB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#7.6%


###chr13-1
chrdata <- data.frame(Y = as.numeric(animalcode), 
                      Marker = as.factor(snp3),
                      Animal = as.factor(animal),
                      Father= as.factor(father))

# Remove rows with missing data
hasMissing <- which(apply(apply(chrdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { chrdata <- chrdata[-hasMissing,] }

#Remove fathers with a small amount of children
chrdata <- chrdata[which(chrdata[, "Father"] %in% enoughoffspring),]
AA<- chrdata[which(chrdata[, "Marker"]== "CC"), ]
AB<- chrdata[which(chrdata[, "Marker"]== "CT"), ]
BB<- chrdata[which(chrdata[, "Marker"]== "TT"), ]
sickAA<-AA[which(AA[, "Y"]==1),]
sickAB<-AA[which(AB[, "Y"]==1),]
sickBB<-AA[which(BB[, "Y"]==1),]

length(unique(AA[, "Animal"]))#1243
length(unique(AB[, "Animal"]))#320
length(unique(BB[, "Animal"]))#22

length(unique(sickAA[, "Animal"]))/length(unique(chrdata[, "Animal"]))#46.4%
length(unique(sickAB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#18.5%
length(unique(sickBB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#2.1%

###chr13-2
chrdata <- data.frame(Y = as.numeric(animalcode), 
                      Marker = as.factor(snp4),
                      Animal = as.factor(animal),
                      Father= as.factor(father))

# Remove rows with missing data
hasMissing <- which(apply(apply(chrdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { chrdata <- chrdata[-hasMissing,] }

#Remove fathers with a small amount of children
chrdata <- chrdata[which(chrdata[, "Father"] %in% enoughoffspring),]
AA<- chrdata[which(chrdata[, "Marker"]== "AA"), ]
AB<- chrdata[which(chrdata[, "Marker"]== "AG"), ]
BB<- chrdata[which(chrdata[, "Marker"]== "GG"), ]
sickAA<-AA[which(AA[, "Y"]==1),]
sickAB<-AA[which(AB[, "Y"]==1),]
sickBB<-AA[which(BB[, "Y"]==1),]

length(unique(AA[, "Animal"]))#610
length(unique(AB[, "Animal"]))#752
length(unique(BB[, "Animal"]))#218

length(unique(sickAA[, "Animal"]))/length(unique(chrdata[, "Animal"]))#23.6%
length(unique(sickAB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#28.4%
length(unique(sickBB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#10.2%


###chr13-3
chrdata <- data.frame(Y = as.numeric(animalcode), 
                      Marker = as.factor(snp5),
                      Animal = as.factor(animal),
                      Father= as.factor(father))

# Remove rows with missing data
hasMissing <- which(apply(apply(chrdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { chrdata <- chrdata[-hasMissing,] }

#Remove fathers with a small amount of children
chrdata <- chrdata[which(chrdata[, "Father"] %in% enoughoffspring),]
AA<- chrdata[which(chrdata[, "Marker"]== "CC"), ]
AB<- chrdata[which(chrdata[, "Marker"]== "CT"), ]
BB<- chrdata[which(chrdata[, "Marker"]== "TT"), ]
sickAA<-AA[which(AA[, "Y"]==1),]
sickAB<-AA[which(AB[, "Y"]==1),]
sickBB<-AA[which(BB[, "Y"]==1),]

length(unique(AA[, "Animal"]))#278
length(unique(AB[, "Animal"]))#776
length(unique(BB[, "Animal"]))#531

length(unique(sickAA[, "Animal"]))/length(unique(chrdata[, "Animal"]))#11.0%
length(unique(sickAB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#17.6%
length(unique(sickBB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#17.6%


###chr18
chrdata <- data.frame(Y = as.numeric(animalcode), 
                      Marker = as.factor(snp9),
                      Animal = as.factor(animal),
                      Father= as.factor(father))

# Remove rows with missing data
hasMissing <- which(apply(apply(chrdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { chrdata <- chrdata[-hasMissing,] }

#Remove fathers with a small amount of children
chrdata <- chrdata[which(chrdata[, "Father"] %in% enoughoffspring),]
AA<- chrdata[which(chrdata[, "Marker"]== "TT"), ]
AB<- chrdata[which(chrdata[, "Marker"]== "GT"), ]
BB<- chrdata[which(chrdata[, "Marker"]== "GG"), ]
sickAA<-AA[which(AA[, "Y"]==1),]
sickAB<-AA[which(AB[, "Y"]==1),]
sickBB<-AA[which(BB[, "Y"]==1),]

length(unique(AA[, "Animal"]))#460
length(unique(AB[, "Animal"]))#825
length(unique(BB[, "Animal"]))#300

length(unique(sickAA[, "Animal"]))/length(unique(chrdata[, "Animal"]))#17.9%
length(unique(sickAB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#29.1%
length(unique(sickBB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#10.6%


###chr19-1
chrdata <- data.frame(Y = as.numeric(animalcode), 
                      Marker = as.factor(snp6),
                      Animal = as.factor(animal),
                      Father= as.factor(father))

# Remove rows with missing data
hasMissing <- which(apply(apply(chrdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { chrdata <- chrdata[-hasMissing,] }

#Remove fathers with a small amount of children
chrdata <- chrdata[which(chrdata[, "Father"] %in% enoughoffspring),]
AA<- chrdata[which(chrdata[, "Marker"]== "GG"), ]
AB<- chrdata[which(chrdata[, "Marker"]== "AG"), ]
BB<- chrdata[which(chrdata[, "Marker"]== "AA"), ]
sickAA<-AA[which(AA[, "Y"]==1),]
sickAB<-AA[which(AB[, "Y"]==1),]
sickBB<-AA[which(BB[, "Y"]==1),]

length(unique(AA[, "Animal"]))#1203
length(unique(AB[, "Animal"]))#353
length(unique(BB[, "Animal"]))#29

length(unique(sickAA[, "Animal"]))/length(unique(chrdata[, "Animal"]))#46.1%
length(unique(sickAB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#17.9%
length(unique(sickBB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#2.6%

###chr19-2
chrdata <- data.frame(Y = as.numeric(animalcode), 
                      Marker = as.factor(snp7),
                      Animal = as.factor(animal),
                      Father= as.factor(father))

# Remove rows with missing data
hasMissing <- which(apply(apply(chrdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { chrdata <- chrdata[-hasMissing,] }

#Remove fathers with a small amount of children
chrdata <- chrdata[which(chrdata[, "Father"] %in% enoughoffspring),]
AA<- chrdata[which(chrdata[, "Marker"]== "CC"), ]
AB<- chrdata[which(chrdata[, "Marker"]== "CT"), ]
BB<- chrdata[which(chrdata[, "Marker"]== "TT"), ]
sickAA<-AA[which(AA[, "Y"]==1),]
sickAB<-AA[which(AB[, "Y"]==1),]
sickBB<-AA[which(BB[, "Y"]==1),]

length(unique(AA[, "Animal"]))#940
length(unique(AB[, "Animal"]))#554
length(unique(BB[, "Animal"]))#91

length(unique(sickAA[, "Animal"]))/length(unique(chrdata[, "Animal"]))#36.2%
length(unique(sickAB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#23.3%
length(unique(sickBB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#5.8%


###chrX
chrdata <- data.frame(Y = as.numeric(animalcode), 
                      Marker = as.factor(snp10),
                      Animal = as.factor(animal),
                      Father= as.factor(father))

# Remove rows with missing data
hasMissing <- which(apply(apply(chrdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { chrdata <- chrdata[-hasMissing,] }

#Remove fathers with a small amount of children
chrdata <- chrdata[which(chrdata[, "Father"] %in% enoughoffspring),]
AA<- chrdata[which(chrdata[, "Marker"]== "TT"), ]
AB<- chrdata[which(chrdata[, "Marker"]== "CT"), ]
BB<- chrdata[which(chrdata[, "Marker"]== "CC"), ]
sickAA<-AA[which(AA[, "Y"]==1),]
sickAB<-AA[which(AB[, "Y"]==1),]
sickBB<-AA[which(BB[, "Y"]==1),]

length(unique(AA[, "Animal"]))#1143
length(unique(AB[, "Animal"]))#408
length(unique(BB[, "Animal"]))#34

length(unique(sickAA[, "Animal"]))/length(unique(chrdata[, "Animal"]))#41.5%
length(unique(sickAB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#20.7%
length(unique(sickBB[, "Animal"]))/length(unique(chrdata[, "Animal"]))#4.2%



###MAF
Animaluse<- mdata1[, "Animal"]
gts<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep="\t",  na.strings = c("","-" ,"NA","no"), fill = TRUE,header = TRUE,colClasses = "character")
genodata<- gts[which(gts[, "ID"] %in% Animaluse), ]
##SNP1-rs41588957-chr6

TT<- subset(genodata, SaS_M2== "TT")
CT<- subset(genodata, SaS_M2== "CT")
CC<- subset(genodata, SaS_M2== "CC")
N11<-nrow(TT)
N12<-nrow(CT)
N22<-nrow(CC)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q


###SNP2-rs110707460-chr6
AA<- subset(genodata, SaS_M3== "AA")
AG<- subset(genodata, SaS_M3== "AG")
GG<- subset(genodata, SaS_M3== "GG")
N11<-nrow(AA)
N12<-nrow(AG)
N22<-nrow(GG)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q

###SNP3-rs109934030-chr13
TT<- subset(genodata, HAS_M6== "TT")
CT<- subset(genodata, HAS_M6== "CT")
CC<- subset(genodata, HAS_M6== "CC")
N11<-nrow(TT)
N12<-nrow(CT)
N22<-nrow(CC)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q
###SNP4-rs41634110-chr13
AA<- subset(genodata, HAS_M3== "AA")
AG<- subset(genodata, HAS_M3== "AG")
GG<- subset(genodata, HAS_M3== "GG")
N11<-nrow(AA)
N12<-nrow(AG)
N22<-nrow(GG)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q
###SNP5-rs109441194-chr13
TT<- subset(genodata, HAS_M4== "TT")
CT<- subset(genodata, HAS_M4== "CT")
CC<- subset(genodata, HAS_M4== "CC")
N11<-nrow(TT)
N12<-nrow(CT)
N22<-nrow(CC)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q
###SNP6-rs41257403-chr19

AA<- subset(genodata, HAS_M8== "AA")
AG<- subset(genodata, HAS_M8== "AG")
GG<- subset(genodata, HAS_M8== "GG")
N11<-nrow(AA)
N12<-nrow(AG)
N22<-nrow(GG)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q

###SNP7-rs41636878-chr19
TT<- subset(genodata, HAS_M9== "TT")
CT<- subset(genodata, HAS_M9== "CT")
CC<- subset(genodata, HAS_M9== "CC")
N11<-nrow(TT)
N12<-nrow(CT)
N22<-nrow(CC)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q
###SNP8-rs41257360-chr5

AA<- subset(genodata, HAS_M1== "AA")
AG<- subset(genodata, HAS_M1== "AG")
GG<- subset(genodata, HAS_M1== "GG")
N11<-nrow(AA)
N12<-nrow(AG)
N22<-nrow(GG)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q
###SNP9-rs29020544-chr18

TT<- subset(genodata, HAS_M7== "TT")
GT<- subset(genodata, HAS_M7== "GT")
GG<- subset(genodata, HAS_M7== "GG")
N11<-nrow(TT)
N12<-nrow(GT)
N22<-nrow(GG)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22 
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q
###SNP10-rs41629005-chr10

TT<- subset(genodata, HAS_M10== "TT")
CT<- subset(genodata, HAS_M10== "CT")
CC<- subset(genodata, HAS_M10== "CC")
N11<-nrow(TT)
N12<-nrow(CT)
N22<-nrow(CC)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22  
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q

#### find closest gene

library(biomaRt)


bio.mart <- useMart(biomart="ensembl", dataset="btaurus_gene_ensembl")


res.biomart <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "chromosome_name", "start_position", "end_position", "strand"), 
                     filters = c("chromosomal_region", "biotype"), 
                     values = list(c("13:78865467:79865467"),"protein_coding"), mart = bio.mart)

res.biomart



####P-VALUE correction
p.adjust(c(0.0024,
           0.1364,
           0.0239,
           0.7439,
           0.0593,
           0.3163,
           0.1255,
           0.0049,
           0.0543,
           0.8421), method = "BH")





