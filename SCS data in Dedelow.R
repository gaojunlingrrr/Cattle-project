##Dedelow analysis

library(foreign)
setwd("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/")
gts<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep="\t",  na.strings = c("","-" ,"NA","no"), fill = TRUE,header = TRUE,colClasses = "character")
mlp<- read.dbf("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/MLP.DBF")
bestand<- read.dbf("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/BESTAND.DBF")

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
    leavedate<- as.character(bestand[inBestand, "ABGANG"])
  }
  
  mlpdata <- c(as.character(mlp[r, "OHR"]), birthdate, firstcalfindays,firstcalf,leavedate,fatherid,
               as.character(mlp[r, "DATUM"]), 
               as.character(mlp[r, "LAKTATION"]), 
               as.character(mlp[r, "MKG"]),as.character(mlp[r, "FETT"]),as.character(mlp[r, "EIWEISS"]), as.character(mlp[r, "ZELLZAHL"]))
  idx = which(as.character(gts[, "ID"]) == as.character(mlp[r, "OHR"]))
  mlpdata <- c(mlpdata, as.character(gts[idx[1], -1]))
  
  mmatrix <- rbind(mmatrix, mlpdata)
  
}


colnames(mmatrix) <- c("AnimalID", "Birthdate", "FirstCalfIndays","calvingyear","Leaveherdtime","Father","Sampledate", "Lactation", "MKG", "Fat", "Protein", "SCC", colnames(gts)[-1])
rownames(mmatrix) <- 1:nrow(mmatrix)


setwd("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/")
write.table(mmatrix, "phenotypes2020.Dedelow.txt", sep = "\t", quote= FALSE) ###quote 是否为字符型变量添加双引号
###
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

# Which covariates influence the phenotype?
library(lme4)
null.model <- lmer(Y ~ (1|Animal), data = mdata, REML = FALSE)
model1 <- lmer(Y ~ Lact + (1|Animal), data = mdata, REML = FALSE)
anova(model1, null.model) # < 2.2e-16 ***

model2 <- lmer(Y ~ Samplingseason + (1|Animal), data = mdata, REML = FALSE)
anova(model2, null.model) # < 1.336e-09 ***

model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata, REML = FALSE)
anova(model3, null.model) # < 0.001351 **

model4 <- lmer(Y ~ Samplingyear + (1|Animal), data = mdata, REML = FALSE)
anova(model4, null.model) # < 2.2e-16 ***

model5 <- lmer(Y ~ Birthyear + (1|Animal), data = mdata, REML = FALSE)
anova(model5, null.model) # < 7.018e-05 ***

model6 <- lmer(Y ~ Firstcalf + (1|Animal), data = mdata, REML = FALSE)
anova(model6, null.model) # < 0.1267 NS

model7 <- lmer(Y ~ Father + (1|Animal), data = mdata, REML = FALSE)
anova(model7, null.model) # < 4.2e-12 ***

model8 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata, REML = FALSE)
anova(model8, null.model) # 0.0002043 ***

genotypes <- mmatrix1[,13:22]
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
snp5<- factor(snp5, levels = c("TT", "CT", "CC"))
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

###
genotypes <- mmatrix1[,13:22]
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata <- data.frame(Y = as.numeric(SCS), 
                 Marker = as.factor(snp10),
                 
                 Samplingseason = as.factor(seasonsampling), 
                 Lact = as.factor(lactation),
                 Season = as.factor(seasoncalving), 
                 Calvingyear = as.factor(yearcalving),
                 Animal = as.factor(animal),
                 Father= as.factor(father))}
  
  # Remove rows with missing data
  hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
  if (length(hasMissing) > 0) { mdata <- mdata[-hasMissing,] }

  #Remove fathers with a small amount of children
  mdata <- mdata[which(mdata[, "Father"] %in% enoughoffspring),]
  
  # Run a null model (all significant and suggestive covariates) and another model including the marker
  null.model <- lmer(Y ~ Samplingseason + Lact + Season + Father + Calvingyear  + (1|Animal), data = mdata, REML = FALSE)
  markermodel <- lmer(Y ~ Samplingseason + Lact + Season + Father+ Calvingyear + Marker + (1|Animal), data = mdata, REML = FALSE)  
  anova(null.model, markermodel)
  markermodel
  
  
  
  
  ###
  r.squaredGLMM(null.model)
  r.squaredGLMM(markermodel)
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
  
  
  
  # Pvalue for the model
  pval <- as.numeric(na.omit(anova(markermodel, null.model)[, "Pr(>Chisq)"]))
  pvals <- c(pvals, pval)
}
names(pvals) <- colnames(genotypes)
pvals

#SaS_M2     SaS_M3     HAS_M6     HAS_M3     HAS_M4     HAS_M8     HAS_M9     HAS_M1     HAS_M7 
#0.62827522 0.49011224 0.01046332 0.03788797 0.02306631 0.36745351 0.20981009 0.01106936 0.90059805 
#HAS_M10 
#0.49516825
###CM

mdata1 <- data.frame(Y = as.numeric(animalcode), 
                    Lact = as.factor(lactation),
                    Animal = as.factor(animal),
                    Birthyear = as.factor(birthyear),
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



# Which covariates influence the phenotype?
library(lme4)
null.model <- glmer(Y ~ (1|Animal), data = mdata1, family = binomial(link = "logit"),nAGQ=0)
#null.model <- lm(Y ~ 1, data = mdata1)

model1 <- glmer(Y ~ Lact + (1|Animal), data = mdata1, family = binomial(link = "logit"))
anova(model1, null.model, test="Chisq") # < NS 
　　　　
model2 <- glmer(Y ~ Calvingseason + (1|Animal), data = mdata1,family = binomial(link = "logit"))　　　
anova(model2, null.model, test="Chisq") # < 2.2e-16 ***

model3 <- glmer(Y ~ Birthyear + (1|Animal), data = mdata1, family = binomial(link = "logit"))
anova(model3, null.model, test="Chisq") # < NS

model4 <- glmer(Y ~ Firstcalf + (1|Animal), data = mdata1, family = binomial(link = "logit"))
anova(model4, null.model, test="Chisq") # < NS

model5 <- glmer(Y ~ Father + (1|Animal), data = mdata1, family = binomial(link = "logit"),nAGQ=0)
anova(model5, null.model, test="Chisq") # < 2.2e-16 ***

 
model6 <- glmer(Y ~ Calvingyear + (1|Animal), data = mdata1, family = binomial(link = "logit"),nAGQ=0)
anova(model6, null.model, test="Chisq") # < 2.2e-16 ***



####

  # Combine all data into a matrix
  mdata1 <- data.frame(Y = as.numeric(animalcode), 
                 Marker = as.factor(snp8),
                 Animal = as.factor(animal),
                 Calvingseason = as.factor(seasoncalving),                 
                 Calvingyear = as.factor(yearcalving),
              
                 Father= as.factor(father))
  
  # Remove rows with missing data
  hasMissing <- which(apply(apply(mdata1,1,is.na),2,sum) != 0)
  if (length(hasMissing) > 0) { mdata1 <- mdata1[-hasMissing,] }
  
  #Remove fathers with a small amount of children
  mdata1 <- mdata1[which(mdata1[, "Father"] %in% enoughoffspring),]
  
  # Run a null model (all significant and suggestive covariates) and another model including the marker
  null.model <- glmer(Y ~ Calvingseason +Calvingyear + Father + (1|Animal), data = mdata1, family = binomial(link = "logit"), nAGQ=0)
  markermodel <- glmer(Y ~ Calvingseason +Calvingyear + Father+ Marker+ (1|Animal) , data = mdata1, family = binomial(link = "logit"),nAGQ=0)  
  anova(null.model, markermodel)
  markermodel
  
  
  ###
  r.squaredGLMM(null.model)
  r.squaredGLMM(markermodel)
  
  
  
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
  
  
  
  
  
  
  
  
  
  
  
  
  # Pvalue for the model
  pval <- as.numeric(na.omit(anova(markermodel, null.model, test="Chisq")[, "Pr(>Chi)"]))
  pvals <- c(pvals, pval)
}
names(pvals) <- colnames(genotypes)
pvals

####

yearofbirth<- unlist(lapply(strsplit(as.character(mmatrix1[,"Birthdate"]), "-"), "[", 1))


range(yearofbirth)

#### find closest gene
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")


library("biomaRt")


bio.mart <- useMart(biomart="ensembl", dataset="btaurus_gene_ensembl")


res.biomart <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "chromosome_name", "start_position", "end_position", "strand"), 
                     filter = c("chromosomal_region", "biotype"), 
                     values = list(c("6:83553915:84053915"),"protein_coding"), mart = bio.mart)


