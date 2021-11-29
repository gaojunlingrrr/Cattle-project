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
  
 
  
  # Which covariates influence the phenotype?
  library(lme4)
  null.model <- glmer(Y ~ (1|Animal), data = mdata1, family = binomial,nAGQ = 0)
  
  
  model1 <- glmer(Y ~ Calvingseason + (1|Animal), data = mdata1,family = binomial, nAGQ = 0)　　　
  anova(model1, null.model, test="Chisq") 
  
  model2 <- glmer(Y ~ Calvingyear + (1|Animal), data = mdata1, family = binomial, nAGQ = 0)
  anova(model2, null.model, test="Chisq") 
  
  model3 <- glmer(Y ~ Firstcalf + (1|Animal), data = mdata1, family = binomial, nAGQ = 0)
  anova(model3, null.model, test="Chisq") 
  
  model4 <- glmer(Y ~ (1|Father) + (1|Animal), data = mdata1, family = binomial, nAGQ = 0)
  anova(model4, null.model, test="Chisq") 
  
  
  cat("Lactation:", lactation)
  cat(", Calvingseason:", as.numeric(na.omit(anova(model1, null.model)[, "Pr(>Chisq)"])))
  cat(", Calvingyear:", as.numeric(na.omit(anova(model2, null.model)[, "Pr(>Chisq)"])))
  cat(", Afc:", as.numeric(na.omit(anova(model3, null.model)[, "Pr(>Chisq)"])))
  cat(", Father:", as.numeric(na.omit(anova(model4, null.model)[, "Pr(>Chisq)"])), "\n")
  
  lactation= lactation + 1
  
}


#Pvalues of covarites for CM1, CM2, and CM3
Lactation: 1, Calvingseason: 3.708009e-07, Calvingyear: 4.869603e-81, Afc: 1.819017e-08, Father: 2.601951e-30 
Lactation: 2, Calvingseason: 0.07235081, Calvingyear: 4.343314e-63, Afc: 0.07033527, Father: 7.914664e-25 
Lactation: 3, Calvingseason: 0.06688018, Calvingyear: 6.099554e-20, Afc: 0.4424046, Father: 0.0003487974 



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
  
  
  
  # Run a null model (all significant and suggestive covariates) and another model including the marker
  null.model <- glmer(Y ~ Calvingseason + Calvingyear + Afc + (1|Father) + (1|Animal), data = data.frame(mdata), family = binomial, nAGQ=0)
  markermodel <- glmer(Y ~ Calvingseason + Calvingyear + Afc + (1|Father) + Marker+ (1|Animal) , data = data.frame(mdata), family = binomial,nAGQ=0)  
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
0.69672002 0.39738992 0.58032179 0.41403504 0.04866959 0.24353723 0.28186180 0.05326573 0.25114750 0.77236158 


####P-VALUE correction
round(p.adjust(c(0.69672002, 0.39738992, 0.58032179, 0.41403504, 0.04866959, 0.24353723, 0.28186180, 0.05326573, 0.25114750, 0.77236158  ), method = "BH"), 4)

[1] 0.7724 0.5915 0.7254 0.5915 0.2663 0.5637 0.5637 0.2663 0.5637 0.7724



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
  
  
  
  # Run a null model (all significant and suggestive covariates) and another model including the marker
  null.model <- glmer(Y ~ Calvingseason + Calvingyear + Afc + (1|Father) + (1|Animal), data = data.frame(mdata), family = binomial, nAGQ=0)
  markermodel <- glmer(Y ~ Calvingseason + Calvingyear + Afc + (1|Father) + Marker+ (1|Animal) , data = data.frame(mdata), family = binomial,nAGQ=0)  
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
0.9052578837 0.0006430409 0.9185649270 0.9419640343 0.5949729782 0.0307849975 0.2725001345 0.4741004689 0.0551572275 0.1327023946 


####P-VALUE correction
round(p.adjust(c(0.9052578837, 0.0006430409, 0.9185649270, 0.9419640343, 0.5949729782, 0.0307849975, 0.2725001345, 0.4741004689, 0.0551572275, 0.1327023946   ), method = "BH"), 4)

[1] 0.9420 0.0064 0.9420 0.9420 0.8500 0.1539 0.5450 0.7902 0.1839 0.3318




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
  
  
  
  # Run a null model (all significant and suggestive covariates) and another model including the marker
  null.model <- glmer(Y ~ Calvingseason + Calvingyear + (1|Father) + (1|Animal), data = data.frame(mdata), family = binomial, nAGQ=0)
  markermodel <- glmer(Y ~ Calvingseason + Calvingyear + (1|Father) + Marker+ (1|Animal) , data = data.frame(mdata), family = binomial,nAGQ=0)  
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
0.40043029 0.34168917 0.18797921 0.13295155 0.07012412 0.28191273 0.66911961 0.10905844 0.99305641 0.05239808 


####P-VALUE correction
round(p.adjust(c(0.40043029, 0.34168917, 0.18797921, 0.13295155, 0.07012412, 0.28191273, 0.66911961, 0.10905844, 0.99305641, 0.05239808 ), method = "BH"), 4)

[1] 0.5005 0.4881 0.3760 0.3324 0.3324 0.4699 0.7435 0.3324 0.9931 0.3324




