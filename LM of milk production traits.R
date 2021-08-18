library(foreign)
  setwd("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/")
  gts<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep="\t",  na.strings = c("","-" ,"NA","no"), fill = TRUE,header = TRUE,colClasses = "character")
  bestand<- read.dbf("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/BESTAND.DBF")
  mp<- read.table("/Users/Pomelo/Documents/R Practising/Dedelow/305data-Dedelow.txt", sep = "\t",fill = TRUE, header = TRUE)
  
  # Figure out MKG == NA or MKG == 0, and remove those from the data 
  noData <- which(is.na(mp[, "MKG"]) | mp[, "MKG"] == 0 | is.na(mp[, "MKG100"]) | mp[, "MKG100"] ==0 | is.na(mp[, "MKG200_N"]) | mp[, "MKG200_N"]==0)
  mp<- mp[-noData,]
  mp <- mp[which(mp[, "OHR"] %in% gts[, "ID"]),]
# How many genotyped cows do not have milk production data ?
length(which(!(gts[, "ID"] %in% mp[, "OHR"])))
# Remove the cows without SCS
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
    leavedate<- as.character(bestand[inBestand, "ABGANG"])
  }
  
  mpdata <- c(as.character(mp[r, "OHR"]), birthdate, firstcalfindays,firstcalf,leavedate,fatherid,
               as.character(mp[r, "KALBUNG"]), 
               as.character(mp[r, "LAKTATION"]), 
               as.character(mp[r, "MKG100"]),as.character(mp[r, "MKG200_N"]),as.character(mp[r, "MKG"]), as.character(mp[r, "FETTKG100"]), as.character(mp[r, "FETTKG200"]),as.character(mp[r, "FETTKG"]),as.character(mp[r, "EIWEISSKG1"]),as.character(mp[r, "EIWEISSKG2"]),as.character(mp[r, "EIWEISSKGG"]))
  idx = which(as.character(gts[, "ID"]) == as.character(mp[r, "OHR"]))
  mpdata <- c(mpdata, as.character(gts[idx[1], -1]))
  
  mmatrix <- rbind(mmatrix, mpdata)
}

colnames(mmatrix) <- c("AnimalID", "Birthdate", "FirstCalfIndays","Afc","Leaveherdtime","Father","Calvingdate", "Lactation", "MKG100", "MKG200_N", "MKG", "FATKG100","FATKG200","FATKG","PROTEINKG100", "PROTEINKG200" , "PROTEINKG"  ,colnames(gts)[-1])
rownames(mmatrix) <- 1:nrow(mmatrix)

setwd("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/")
write.table(mmatrix, "phenotypes-milkproduction.Dedelow.txt", sep = "\t", quote= FALSE) ###quote 是否为字符型变量添加双引号

#####
mmatrix<-read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/phenotypes-milkproduction.Dedelow.txt", sep = "\t",header = TRUE)
Mergingdata<-read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Mergingdata.txt",sep = "\t",header = TRUE)
mmatrix<- mmatrix[which(mmatrix[, "Lactation"]<4 & mmatrix[, "Lactation"] > 0), ]
####select animal with birth year >2007
birthyearofcows <- as.numeric(unlist(lapply(strsplit(as.character(mmatrix[,"Birthdate"]), "-"), "[", 1)))
mmatrix<- cbind(mmatrix, birthyearofcows)
mmatrix<- mmatrix[which(mmatrix[, "birthyearofcows"] > 2007), ]
###Delete elements with less than 3 reccords
frequency<- as.data.frame(table(mmatrix[, "AnimalID"]))
frequency2<- frequency[which(frequency[, "Freq"]>2),]
qfanimal<- frequency2[, "Var1"]
mmatrix<- mmatrix[which(mmatrix[, "AnimalID"] %in% qfanimal), ]
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

mmatrix<-mmatrix[,-18]
mmatrix1<- cbind(mmatrix,Animalcode)

####Get season
calvingdates<- as.Date(as.character(mmatrix1[, "Calvingdate"]), format="%y/%m/%d")
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
mmatrix1<- cbind(mmatrix1, Calvingseason, Calvingyear)


# Get the covariates we want to investigate in model building
animal <- mmatrix1[, "AnimalID"]
lactation <- as.factor(mmatrix1[, "Lactation"])
birthyear <- unlist(lapply(strsplit(as.character(mmatrix1[,"Birthdate"]), "-"), "[", 1))
calvlingyear <- mmatrix1[, "Calvingyear"]
calvingseason<- mmatrix1[, "Calvingseason"]
afc <- mmatrix1[, "FirstCalfIndays"]
father<- mmatrix1[, "Father"]
mkg100<- mmatrix1[, "MKG100"]
mkg200<- mmatrix1[, "MKG200_N"]
mkg<- mmatrix1[, "MKG"]
fatkg100<- 100*(mmatrix1[, "FATKG100"]/ mmatrix1[, "MKG100"])
fatkg200<- 100*(mmatrix1[, "FATKG200"]/ mmatrix1[, "MKG200_N"])
fatkg<-100*(mmatrix1[, "FATKG"]/ mmatrix1[, "MKG"])
proteinkg100<- 100*(mmatrix1[, "PROTEINKG100"]/ mmatrix1[, "MKG100"])
proteinkg200<-100* (mmatrix1[, "PROTEINKG200"]/ mmatrix1[, "MKG200_N"])
proteinkg<- 100*(mmatrix1[, "PROTEINKG"]/ mmatrix1[, "MKG"])

###Set  Reference allele as default
genotypes <- mmatrix1[,18:27]
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


# Combine all data into a matrix
mdata <- data.frame(Y = as.numeric(proteinkg), 
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


# ###
# weights <- c()
# for(x in 1:nrow(mdata)){
#   nInMdata <- length(which(mdata[, "Animal"] == mdata[x, "Animal"]))
#   weights <- c(weights, 1.0 / nInMdata)
# }
# mdata <- cbind(mdata, weight = weights)
# library(lme4)
# # Not deal with the repeated measurement, will over estimate the effect
# anova(lm(Y ~ Lact, data = mdata))
# 
# # Deal with repeat measurements by weighing them less in the regression, 
# # will correctly estimate the effect, but overestimate the significance
# anova(lm(Y ~ Lact, weights = mdata[, "weight"], data = mdata))
# 
# # Deal with repeat measurements by fitting multiple regression lines per individual 
# # will correctly estimate the effect, will correctly estimate the significance
# m1 <- lmer(Y ~ Lact + (1|Animal), data = mdata)
# mNull <- lmer(Y ~ (1|Animal), data = mdata)
# anova(m1,mNull)[[8]]



# Which covariates influence the phenotype?
#### 
library(lme4)
null.model <- lmer(Y ~ (1|Animal), data = mdata,REML = FALSE)

model1 <- lmer(Y ~ Lact + (1|Animal), data = mdata,REML = FALSE)
anova(model1, null.model) # < 2.2e-16 ***

model2 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata,REML = FALSE)
anova(model2, null.model) # < 3.805e-09 ***

model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata,REML = FALSE)
anova(model3, null.model) # < 0.0006105 ***


model4 <- lmer(Y ~ Afc + (1|Animal), data = mdata,REML = FALSE)
anova(model4, null.model) # <  0.6333 NS

model6 <- lmer(Y ~ Father + (1|Animal), data = mdata, REML = FALSE)
anova(model6, null.model) # < 2.2e-16 ***



###
genotypes <- mmatrix1[,18:27]
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata1 <- data.frame(Y = as.numeric(proteinkg), 
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
  null.model <- lmer(Y ~ Lact + Calvingseason + Calvingyear + Father + (1|Animal) , data = mdata1, REML = FALSE)
  markermodel <- lmer(Y ~ Lact + Calvingseason + Calvingyear + Father + Marker + (1|Animal), data = mdata1, REML = FALSE)
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
  
  
  
  
  
  
  
  
  
  
  # Pvalue for the model
  pval <- as.numeric(na.omit(anova(markermodel, null.model)[, "Pr(>Chisq)"]))
  pvals <- c(pvals, pval)
}
names(pvals) <- colnames(genotypes)
pvals










