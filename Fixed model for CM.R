############CM
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

##### Divided into lactation 1, 2, 3
inL1<- which(cmdata[, "Lactation"] ==1)
inL2<- which(cmdata[, "Lactation"] ==2)
inL3<- which(cmdata[, "Lactation"] ==3)


lactation=1

for(individuals in list(inL1,inL2,inL3)){
  
  mmatrix<- cmdata[individuals, ]

# Get the covariates we want to investigate in model building
animal <- mmatrix[, "AnimalID"]
seasoncalving<- mmatrix[, "Calvingseason"]
yearcalving <-  mmatrix[, "Calvingyear_N"]
firstCalf <- mmatrix[, "FirstCalfIndays"]
AFC<- round(firstCalf/30) #21-37
father<- mmatrix[, "Father"]
animalcode<- mmatrix[, "CM"] 


##combine all data into a matrix


mdata1 <- data.frame(Y = as.numeric(animalcode), 
                     
                     Animal = as.factor(animal),
                     Calvingseason = as.factor(seasoncalving),
                     Calvingyear = as.factor(yearcalving),
                     Firstcalf = as.numeric(AFC),#21-37
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

mdata1 <- mdata1[which(mdata1[, "Father"] %in% enoughoffspring),]


# Which covariates influence the phenotype?
library(lme4)
null.model <- glmer(Y ~ (1|Animal), data = mdata1, family = binomial,nAGQ = 0)


model1 <- glmer(Y ~ Calvingseason + (1|Animal), data = mdata1,family = binomial, nAGQ = 0)　　　
anova(model1, null.model, test="Chisq") # < 4.833e-05 ***

model2 <- glmer(Y ~ Calvingyear + (1|Animal), data = mdata1, family = binomial, nAGQ = 0)
anova(model2, null.model, test="Chisq") # < 2.2e-16 ***

model3 <- glmer(Y ~ Firstcalf + (1|Animal), data = mdata1, family = binomial, nAGQ = 0)
anova(model3, null.model, test="Chisq") # < 1.839e-12 ***

model4 <- glmer(Y ~ Father + (1|Animal), data = mdata1, family = binomial, nAGQ = 0)
anova(model4, null.model, test="Chisq") # < 2.2e-16 ***


cat("Lactation:", lactation)
cat(", Calvingseason:", as.numeric(na.omit(anova(model1, null.model)[, "Pr(>Chisq)"])))
cat(", Calvingyear:", as.numeric(na.omit(anova(model2, null.model)[, "Pr(>Chisq)"])))
cat(", Afc:", as.numeric(na.omit(anova(model3, null.model)[, "Pr(>Chisq)"])))
cat(", Father:", as.numeric(na.omit(anova(model4, null.model)[, "Pr(>Chisq)"])), "\n")

lactation= lactation + 1

}


#Pvalues of covarites for CM1, CM2, and CM3
Lactation: 1, Calvingseason: 0.0002309639, Calvingyear: 8.53591e-36, Afc: 7.642581e-05, Father: 3.546783e-23 
Lactation: 2, Calvingseason: 0.1658807, Calvingyear: 3.07656e-33, Afc: 0.07593098, Father: 1.702671e-17 
Lactation: 3, Calvingseason: 0.6485439, Calvingyear: 2.875661e-11, Afc: 0.4954228, Father: 3.951342e-07



###Association analysis for CM1
#loding dataset
mmatrix3<- cmdata[inL1, ]

# Get the covariates we want to investigate in model building
animal <- mmatrix3[, "AnimalID"]
seasoncalving<- mmatrix3[, "Calvingseason"]
yearcalving <-  mmatrix3[, "Calvingyear_N"]
firstCalf <- mmatrix3[, "FirstCalfIndays"]
AFC<- round(firstCalf/30) #21-37
father<- mmatrix3[, "Father"]
animalcode<- mmatrix3[, "CM"] 

#Loading SNP data
genotypes <- mmatrix3[,9:18]
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



pvals <- NULL
for (x in 1:ncol(genotypes)) {
mdata <- data.frame(Y = as.numeric(animalcode), 
                     Marker = as.factor(genotypes[, x]), 
                     Animal = as.factor(animal),
                     Calvingseason = as.factor(seasoncalving),                 
                     Calvingyear = as.factor(yearcalving),
                     Afc = as.numeric(AFC),
                     Father= as.factor(father))

# Remove rows with missing data
hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata <- mdata[-hasMissing,] }

#Remove fathers with a small amount of children
offspring <- c()
for(f in unique(mdata[, "Father"])){
  nOffspring = length(unique(mdata[which(mdata[, "Father"] == f), "Animal"]))
  offspring <- c(offspring, nOffspring)
}
names(offspring) <- unique(mdata[, "Father"])

enoughoffspring <- names(offspring)[which(offspring >= 10)]

mdata <- mdata[which(mdata[, "Father"] %in% enoughoffspring),]

# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- glmer(Y ~ Calvingseason + Calvingyear + Afc + Father + (1|Animal), data = data.frame(mdata), family = binomial, nAGQ=0)
markermodel <- glmer(Y ~ Calvingseason + Calvingyear + Afc + Father + Marker+ (1|Animal) , data = data.frame(mdata), family = binomial,nAGQ=0)  
anova(null.model, markermodel)
markermodel

# Pvalue for the model
pval <- as.numeric(na.omit(anova(markermodel, null.model)[, "Pr(>Chisq)"]))
pvals <- c(pvals, pval)

}
names(pvals) <- colnames(genotypes)

pvals

##P-values

SaS_M2       SaS_M3       HAS_M6       HAS_M3       HAS_M4       HAS_M8       HAS_M9       HAS_M1       HAS_M7      HAS_M10 
3.342524e-01 1.535136e-01 6.064183e-01 5.352446e-02 2.504998e-01 3.794348e-03 1.071277e-02 1.901639e-05 8.561648e-02 8.440390e-01 


####P-VALUE correction
round(p.adjust(c(3.342524e-01, 1.535136e-01, 6.064183e-01, 5.352446e-02, 2.504998e-01, 3.794348e-03, 1.071277e-02, 1.901639e-05, 8.561648e-02, 8.440390e-01  ), method = "BH"), 4)

[1] 0.4178 0.2559 0.6738 0.1338 0.3579 0.0190 0.0357 0.0002 0.1712 0.8440




###Association analysis for CM2
#loding dataset
mmatrix3<- cmdata[inL2, ]

# Get the covariates we want to investigate in model building
animal <- mmatrix3[, "AnimalID"]
seasoncalving<- mmatrix3[, "Calvingseason"]
yearcalving <-  mmatrix3[, "Calvingyear_N"]
firstCalf <- mmatrix3[, "FirstCalfIndays"]
AFC<- round(firstCalf/30) #21-37
father<- mmatrix3[, "Father"]
animalcode<- mmatrix3[, "CM"] 

#Loading SNP data
genotypes <- mmatrix3[,9:18]
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

pvals <- NULL
for (x in 1:ncol(genotypes)) {
  mdata <- data.frame(Y = as.numeric(animalcode), 
                      Marker = as.factor(genotypes[, x]), 
                      Animal = as.factor(animal),
                      Calvingseason = as.factor(seasoncalving),                 
                      Calvingyear = as.factor(yearcalving),
                      Afc = as.numeric(AFC),
                      Father= as.factor(father))
  
  # Remove rows with missing data
  hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
  if (length(hasMissing) > 0) { mdata <- mdata[-hasMissing,] }
  
  #Remove fathers with a small amount of children
  offspring <- c()
  for(f in unique(mdata[, "Father"])){
    nOffspring = length(unique(mdata[which(mdata[, "Father"] == f), "Animal"]))
    offspring <- c(offspring, nOffspring)
  }
  names(offspring) <- unique(mdata[, "Father"])
  
  enoughoffspring <- names(offspring)[which(offspring >= 10)]
  
  mdata <- mdata[which(mdata[, "Father"] %in% enoughoffspring),]
  
  # Run a null model (all significant and suggestive covariates) and another model including the marker
  null.model <- glmer(Y ~ Calvingyear + Afc + Father + (1|Animal), data = data.frame(mdata), family = binomial, nAGQ=0)
  markermodel <- glmer(Y ~ Calvingyear + Afc + Father + Marker+ (1|Animal) , data = data.frame(mdata), family = binomial,nAGQ=0)  
  anova(null.model, markermodel)
  markermodel
  
  # Pvalue for the model
  pval <- as.numeric(na.omit(anova(markermodel, null.model)[, "Pr(>Chisq)"]))
  pvals <- c(pvals, pval)
  
}
names(pvals) <- colnames(genotypes)

pvals

##P-values

SaS_M2     SaS_M3     HAS_M6     HAS_M3     HAS_M4     HAS_M8     HAS_M9     HAS_M1     HAS_M7    HAS_M10 
0.74017562 0.03530723 0.62503202 0.11389714 0.68719151 0.08235440 0.06257772 0.89761951 0.16091558 0.44679163 



####P-VALUE correction
round(p.adjust(c(0.74017562, 0.03530723, 0.62503202, 0.11389714, 0.68719151, 0.08235440, 0.06257772, 0.89761951, 0.16091558, 0.44679163  ), method = "BH"), 4)

[1] 0.8224 0.2745 0.8224 0.2847 0.8224 0.2745 0.2745 0.8976 0.3218 0.7447




###Association analysis for CM3
#loding dataset
mmatrix3<- cmdata[inL3, ]

# Get the covariates we want to investigate in model building
animal <- mmatrix3[, "AnimalID"]
seasoncalving<- mmatrix3[, "Calvingseason"]
yearcalving <-  mmatrix3[, "Calvingyear_N"]
firstCalf <- mmatrix3[, "FirstCalfIndays"]
AFC<- round(firstCalf/30) #21-37
father<- mmatrix3[, "Father"]
animalcode<- mmatrix3[, "CM"] 

#Loading SNP data
genotypes <- mmatrix3[,9:18]
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

pvals <- NULL
for (x in 1:ncol(genotypes)) {
  mdata <- data.frame(Y = as.numeric(animalcode), 
                      Marker = as.factor(genotypes[, x]), 
                      Animal = as.factor(animal),
                      Calvingseason = as.factor(seasoncalving),                 
                      Calvingyear = as.factor(yearcalving),
                      Afc = as.numeric(AFC),
                      Father= as.factor(father))
  
  # Remove rows with missing data
  hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
  if (length(hasMissing) > 0) { mdata <- mdata[-hasMissing,] }
  
  #Remove fathers with a small amount of children
  offspring <- c()
  for(f in unique(mdata[, "Father"])){
    nOffspring = length(unique(mdata[which(mdata[, "Father"] == f), "Animal"]))
    offspring <- c(offspring, nOffspring)
  }
  names(offspring) <- unique(mdata[, "Father"])
  
  enoughoffspring <- names(offspring)[which(offspring >= 10)]
  
  mdata <- mdata[which(mdata[, "Father"] %in% enoughoffspring),]
  
  # Run a null model (all significant and suggestive covariates) and another model including the marker
  null.model <- glmer(Y ~ Calvingyear + Father + (1|Animal), data = data.frame(mdata), family = binomial, nAGQ=0)
  markermodel <- glmer(Y ~ Calvingyear + Father + Marker+ (1|Animal) , data = data.frame(mdata), family = binomial,nAGQ=0)  
  anova(null.model, markermodel)
  markermodel
  
  # Pvalue for the model
  pval <- as.numeric(na.omit(anova(markermodel, null.model)[, "Pr(>Chisq)"]))
  pvals <- c(pvals, pval)
  
}
names(pvals) <- colnames(genotypes)

pvals

##P-values

SaS_M2     SaS_M3     HAS_M6     HAS_M3     HAS_M4     HAS_M8     HAS_M9     HAS_M1     HAS_M7    HAS_M10 
0.01421638 0.71426130 1.00000000 0.83255152 0.73612422 0.21421950 0.50755032 0.65298718 0.98483044 0.17439675 



####P-VALUE correction
round(p.adjust(c(0.01421638, 0.71426130, 1.00000000, 0.83255152, 0.73612422, 0.21421950, 0.50755032, 0.65298718, 0.98483044, 0.17439675 ), method = "BH"), 4)

[1] 0.1422 1.0000 1.0000 1.0000 1.0000 0.7141 1.0000 1.0000 1.0000 0.7141







