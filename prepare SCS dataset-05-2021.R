##read SCS dataset
library(foreign)
setwd("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/")
gts<- read.table("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep="\t",  na.strings = c("","-" ,"NA","no"), fill = TRUE,header = TRUE,colClasses = "character")
mlp<- read.dbf("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/MLP.DBF")
bestand<- read.dbf("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/BESTAND.DBF")

# Figure out MKG == NA or MKG == 0, and remove those from the data (no milk, no reliable SCS)
noData <- which(is.na(mlp[, "MKG"]) | mlp[, "MKG"] == 0)
mlp<- mlp[-noData,]
# Remove the NA cell counts
mlp <- mlp[-which(is.na(mlp[,"ZELLZAHL"])),]
mlp <- mlp[which(mlp[, "OHR"] %in% gts[, "ID"]),]
# How many genotyped cows do not have SCS data ?
length(which(!(gts[, "ID"] %in% mlp[, "OHR"])))
# Remove the cows without SCS
gts <- gts[-which(!(gts[, "ID"] %in% mlp[, "OHR"])),]
mmatrix <- NULL
for(r in 1:nrow(mlp)){
  ohr <- mlp[r, "OHR"]
  inBestand <- which(as.character(bestand[, "OHR"]) == as.character(ohr))
  birthdate <- NA
  firstcalf <- NA
  if(length(inBestand) == 1){ 
    birthdate <- as.character(bestand[inBestand, "GEBURT"])
    firstcalf <- as.character(bestand[inBestand, "KALBUNG1"])
    firstcalfindays <- as.character(as.Date(firstcalf, format="%Y-%m-%d") - as.Date(birthdate, format="%Y-%m-%d"))
    fatherid<- as.character(bestand[inBestand, "VATER"])
    
  }
  
  mlpdata <- c(as.character(mlp[r, "OHR"]), birthdate, firstcalfindays,firstcalf,fatherid,
               as.character(mlp[r, "DATUM"]), 
               as.character(mlp[r, "LAKTATION"]), 
                as.character(mlp[r, "ZELLZAHL"]))
  idx = which(as.character(gts[, "ID"]) == as.character(mlp[r, "OHR"]))
  mlpdata <- c(mlpdata, as.character(gts[idx[1], -1]))
  
  mmatrix <- rbind(mmatrix, mlpdata)
  
}


colnames(mmatrix) <- c("AnimalID", "Birthdate", "FirstCalfIndays","calvingyear","Father","Sampledate", "Lactation",  "SCC", colnames(gts)[-1])
rownames(mmatrix) <- 1:nrow(mmatrix)


setwd("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/")
write.table(mmatrix, "phenotypes2020.Dedelow.txt", sep = "\t", quote= FALSE) ###quote 是否为字符型变量添加双引号






mmatrix1<-read.table("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/phenotypes2020.Dedelow.txt", sep = "\t",header = TRUE)

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

### select lactation 1-3
mmatrix1<- mmatrix1[which(mmatrix1[, "Lactation"] < 4 & mmatrix1[, "Lactation"] > 0), ]
####select animal with birth year >2008
birthyear <- as.numeric(unlist(lapply(strsplit(as.character(mmatrix1[,"Birthdate"]), "-"), "[", 1)))
mmatrix1<- cbind(mmatrix1, birthyear)
mmatrix1<- mmatrix1[which(mmatrix1[, "birthyear"] > 2008), ]

###Delete elements with less than 3 reccords
frequency<- as.data.frame(table(mmatrix1[, "AnimalID"]))
frequency2<- frequency[which(frequency[, "Freq"]>2),]
qfanimal<- frequency2[, "Var1"]
mmatrix1<- mmatrix1[which(mmatrix1[, "AnimalID"] %in% qfanimal), ]

SCS <- log2(as.numeric(as.character(mmatrix1[, "SCC"])) / 100) + 3
SCS[which(!is.finite(SCS))] <- 0

mmatrix1<- cbind(mmatrix1, SCS)

##
Calvingyear<- unlist(lapply(strsplit(as.character(mmatrix1[,"calvingyear"]), "-"), "[", 1))

mmatrix1<- cbind(mmatrix1, Calvingyear)
# Get the covariates we want to investigate in model building
animal <- mmatrix1[, "AnimalID"]
lactation <- mmatrix1[, "Lactation"]
birthyear <- mmatrix1[, "birthyear"]
seasonsampling <- mmatrix1[, "Samplingseason"]
seasoncalving<- mmatrix1[, "Calvingseason"]
yearcalving<- mmatrix1[, "Calvingyear"]
firstCalf <- mmatrix1[, "FirstCalfIndays"]
father<- mmatrix1[, "Father"]
SCS<- mmatrix1[, "SCS"]

# Combine all data into a matrix
mdata <- data.frame(Y = as.numeric(SCS), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Samplingseason = as.factor(seasonsampling),
                    Calvingseason = as.factor(seasoncalving),
                    Calvingyear = as.factor(yearcalving),
                    Afc = as.numeric(firstCalf),
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
length(unique(mdata[, "Animal"])) ##1844

##
setwd("/Users/gaojunling/Desktop/PhD thesis/Manuscript-21.05-Junling/")
table1data<- write.table(mdata, file= "table1data.txt", sep = "\t")
setwd("/Users/gaojunling/Desktop/PhD thesis/Manuscript-21.05-Junling/")
scsdata<-write.table(mmatrix1, file= "scsdata.txt", sep = "\t")
###number of animal per lactation, table1
lact1<- mdata[which(mdata[, "Lact"] ==1), ]
lact2<- mdata[which(mdata[, "Lact"] ==2), ]
lact3<- mdata[which(mdata[, "Lact"] ==3), ]
length(unique(lact1[, "Animal"]))#1844
length(unique(lact2[, "Animal"]))#1333
length(unique(lact3[, "Animal"]))#712
nrow(lact1)#18856
nrow(lact2)#11520
nrow(lact3)#5856
mean(lact1[, "Y"])
mean(lact2[, "Y"])
mean(lact3[, "Y"])
mean(mdata[, "Y"])
min(lact1[, "Y"])
min(lact2[, "Y"])
min(lact3[, "Y"])
min(mdata[, "Y"])
max(lact1[, "Y"])
max(lact2[, "Y"])
max(lact3[, "Y"])
max(mdata[, "Y"])
sd(lact1[, "Y"])
sd(lact2[, "Y"])
sd(lact3[, "Y"])
sd(mdata[, "Y"])
sd(lact1[, "Y"])/sqrt(nrow(lact1))
sd(lact2[, "Y"])/sqrt(nrow(lact2))
sd(lact3[, "Y"])/sqrt(nrow(lact3))
sd(mdata[, "Y"])/sqrt(nrow(mdata))

##Animal selected
scsanimal<- unique(mdata[, "Animal"])
Totalanimal<- c(mpanimal, cmanimal, scsanimal)
length(unique(Totalanimal)) ###2174
scssire<- unique(mdata[, "Father"])
Totalsire<- c(mpsire, cmsire, scssire)
length(unique(Totalsire)) ##158

# Which covariates influence the phenotype?
library(lme4)
null.model <- lmer(Y ~ (1|Animal), data = mdata, REML = FALSE)
model1 <- lmer(Y ~ Lact + (1|Animal), data = mdata, REML = FALSE)
anova(model1, null.model) # < 2.2e-16 ***

model2 <- lmer(Y ~ Samplingseason + (1|Animal), data = mdata, REML = FALSE)
anova(model2, null.model) # < 5.805e-06 ***

model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata, REML = FALSE)
anova(model3, null.model) # < 0.01067 *

model4 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata, REML = FALSE)
anova(model4, null.model) # < 0.0001964 ***

model5 <- lmer(Y ~ Afc + (1|Animal), data = mdata, REML = FALSE)
anova(model5, null.model) # <  NS

model6 <- lmer(Y ~ Father + (1|Animal), data = mdata, REML = FALSE)
anova(model6, null.model) # < 3.233e-12 ***

genotypes <- mmatrix1[,10:19]
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



###Association analysis
genotypes <- mmatrix1[,10:19]
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata <- data.frame(Y = as.numeric(SCS), 
                      Marker = as.factor(snp5),
                      Samplingseason = as.factor(seasonsampling), 
                      Lact = as.factor(lactation),
                      Calvingseason = as.factor(seasoncalving), 
                      Calvingyear = as.factor(yearcalving),
                      Animal = as.factor(animal),
                      Father= as.factor(father))}

# Remove rows with missing data
hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata <- mdata[-hasMissing,] }

#Remove fathers with a small amount of children
mdata <- mdata[which(mdata[, "Father"] %in% enoughoffspring),]

# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- lmer(Y ~ Samplingseason + Lact + Calvingseason + Father + Calvingyear  + (1|Animal), data = mdata, REML = FALSE)
markermodel <- lmer(Y ~ Samplingseason + Lact + Calvingseason + Father+ Calvingyear + Marker + (1|Animal), data = mdata, REML = FALSE)  
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

# To get variance explained we need to use lmer model, since lme doesn't allow us to get sum of squares
100 * (variance.explained.scaled(null.model)["Unexplained"] - variance.explained.scaled(markermodel)["Unexplained"])




####P-VALUE correction
p.adjust(c(0.0076,
           0.8260,
           0.3422,
           0.6745,
           0.2512,
           0.5270,
           0.9211,
           0.7262,
           0.9249,
           0.8642), method = "BH")



