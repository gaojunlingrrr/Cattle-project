####SCS
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

####select animal with birth year >2008
birthyear <- as.numeric(unlist(lapply(strsplit(as.character(mmatrix1[,"Birthdate"]), "-"), "[", 1)))
mmatrix1<- cbind(mmatrix1, birthyear)
mmatrix1<- mmatrix1[which(mmatrix1[, "birthyear"] > 2008), ]


### select lactation 1-3
mmatrix1<- mmatrix1[which(mmatrix1[, "Lactation"] < 4 & mmatrix1[, "Lactation"] > 0), ]

###Delete elements with less than 3 reccords
frequency<- as.data.frame(table(mmatrix1[, "AnimalID"]))
frequency2<- frequency[which(frequency[, "Freq"]>2),]
qfanimal<- frequency2[, "Var1"]
mmatrix1<- mmatrix1[which(mmatrix1[, "AnimalID"] %in% qfanimal), ]

# SCC TRANSFORMED TO SCS
SCS <- log2(as.numeric(as.character(mmatrix1[, "SCC"])) / 100) + 3
SCS[which(!is.finite(SCS))] <- 0

mmatrix1<- cbind(mmatrix1, SCS)

## ADD THE COULUMN OF CALVING YEAR
Calvingyear<- unlist(lapply(strsplit(as.character(mmatrix1[,"calvingyear"]), "-"), "[", 1))

mmatrix1<- cbind(mmatrix1, Calvingyear)



## Divided into lactation 1, 2, 3
inL1<- which(mmatrix1[, "Lactation"] == 1)
inL2<- which(mmatrix1[, "Lactation"] == 2)
inL3<- which(mmatrix1[, "Lactation"] == 3)

lactation=1

for(individuals in list(inL1,inL2,inL3)){
  
  mmatrix<- mmatrix1[individuals, ]
  
  # Get the covariates we want to investigate in model building
  animal <- mmatrix[, "AnimalID"]
  seasonsampling <- mmatrix[, "Samplingseason"]
  seasoncalving<- mmatrix[, "Calvingseason"]
  yearcalving<- mmatrix[, "Calvingyear"]
  firstCalf <- mmatrix[, "FirstCalfIndays"]
  AFC<- round(firstCalf/30) ## 21-37
  DIM<- mmatrix[, "DIM"]
  father<- mmatrix[, "Father"]
  SCS<- mmatrix[, "SCS"]
  
  
  # Combine all data into a matrix
  mdata <- data.frame(Y = as.numeric(SCS), 
                      
                      Animal = as.factor(animal),
                      Samplingseason = as.factor(seasonsampling),
                      Calvingseason = as.factor(seasoncalving),
                      Calvingyear = as.factor(yearcalving),
                      Afc = as.numeric(AFC), #21-37
                      DIM = as.numeric(DIM),
                      Father= as.factor(father)
  )
  # Remove rows with missing data
  mdata<- na.omit(mdata)
  
  
  
  # Which covariates influence the phenotype?
  library(lme4)
  null.model <- lmer(Y ~ (1|Animal), data = mdata, REML = FALSE)
  model1 <- lmer(Y ~ DIM + (1|Animal), data = mdata, REML = FALSE)
  anova(model1, null.model) 
  
  model2 <- lmer(Y ~ Samplingseason + (1|Animal), data = mdata, REML = FALSE)
  anova(model2, null.model)  
  
  model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata, REML = FALSE)
  anova(model3, null.model) 
  
  model4 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata, REML = FALSE)
  anova(model4, null.model) 
  
  model5 <- lmer(Y ~ Afc + (1|Animal), data = mdata, REML = FALSE)
  anova(model5, null.model) 
  
  model6 <- lmer(Y ~ (1|Father) + (1|Animal), data = mdata, REML = FALSE)
  anova(model6, null.model) 
  
  cat("Lactation:", lactation)
  cat(", DIM:", as.numeric(na.omit(anova(model1, null.model)[, "Pr(>Chisq)"])))
  cat(", Samplingseason:", as.numeric(na.omit(anova(model2, null.model)[, "Pr(>Chisq)"])))
  cat(", Calvingseason:", as.numeric(na.omit(anova(model3, null.model)[, "Pr(>Chisq)"])))
  cat(", Calvingyear:", as.numeric(na.omit(anova(model4, null.model)[, "Pr(>Chisq)"])))
  cat(", Afc:", as.numeric(na.omit(anova(model5, null.model)[, "Pr(>Chisq)"])))
  cat(", Father:", as.numeric(na.omit(anova(model6, null.model)[, "Pr(>Chisq)"])), "\n")
  
  lactation= lactation + 1
  
}


#Pvalues of covarites for SCS1, SCS2, and SCS3

Lactation: 1, DIM: 1.067434e-23, Samplingseason: 0.0001623048, Calvingseason: 3.140637e-06, Calvingyear: 1.916381e-11, Afc: 0.06462796, Father: 6.476886e-13 
Lactation: 2, DIM: 4.075e-263, Samplingseason: 0.0004825206, Calvingseason: 0.8920693, Calvingyear: 3.787526e-11, Afc: 0.8761192, Father: 0.0003214427 
Lactation: 3, DIM: 1.892904e-205, Samplingseason: 4.493681e-07, Calvingseason: 0.1290677, Calvingyear: 5.691753e-14, Afc: 0.5313662, Father: 0.0001011047 


###Association analysis for SCS1

#loding dataset
mmatrix2<- mmatrix1[inL1, ]
# Get the significant covariates  in model building
animal <- mmatrix2[, "AnimalID"]
seasonsampling <- mmatrix2[, "Samplingseason"]
seasoncalving<- mmatrix2[, "Calvingseason"]
yearcalving<- mmatrix2[, "Calvingyear"]
DIM<- mmatrix2[, "DIM"]
firstCalf <- mmatrix2[, "FirstCalfIndays"]
AFC<- round(firstCalf/30) ## 21-37
father<- mmatrix2[, "Father"]
SCS<- mmatrix2[, "SCS"]
#Loading SNP data
genotypes <- mmatrix2[,11:20]
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
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata <- data.frame(Y = as.numeric(SCS), 
                      Marker = as.factor(genotypes[,x]), 
                      Samplingseason = as.factor(seasonsampling), 
                      DIM = as.numeric(DIM),
                      Calvingseason = as.factor(seasoncalving), 
                      Calvingyear = as.factor(yearcalving),
                      Animal = as.factor(animal),
                      Afc = as.factor(AFC),
                      Father= as.factor(father))
  
  # Remove rows with missing data
  hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
  if (length(hasMissing) > 0) { mdata <- mdata[-hasMissing,] }
  
 
  
  # Run a null model (all significant and suggestive covariates) and another model including the marker
  null.model <- lmer(Y ~ Samplingseason   + Calvingseason + Calvingyear + DIM + Afc + (1|Father) + (1|Animal), data = data.frame(mdata), REML = FALSE)
  markermodel <- lmer(Y ~ Samplingseason  + Calvingseason + Calvingyear + DIM + Afc + (1|Father) + Marker + (1|Animal), data = data.frame(mdata), REML = FALSE)  
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
0.49919948 0.61397816 0.42329174 0.45101026 0.54060753 0.48083171 0.46523686 0.07290838 0.56105718 0.35452211 

####P-VALUE correction
round(p.adjust(c(0.49919948, 0.61397816, 0.42329174, 0.45101026, 0.54060753, 0.48083171, 0.46523686, 0.07290838, 0.56105718, 0.35452211 ), method = "BH"), 4)

[1] 0.614 0.614 0.614 0.614 0.614 0.614 0.614 0.614 0.614 0.614



###Association analysis for SCS2

#loding dataset
mmatrix2<- mmatrix1[inL2, ]
# Get the significant covariates  in model building
animal <- mmatrix2[, "AnimalID"]
seasonsampling <- mmatrix2[, "Samplingseason"]
seasoncalving<- mmatrix2[, "Calvingseason"]
yearcalving<- mmatrix2[, "Calvingyear"]
DIM<- mmatrix2[, "DIM"]
firstCalf <- mmatrix2[, "FirstCalfIndays"]
AFC<- round(firstCalf/30) ## 21-37
father<- mmatrix2[, "Father"]
SCS<- mmatrix2[, "SCS"]
#Loading SNP data
genotypes <- mmatrix2[,11:20]
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
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata <- data.frame(Y = as.numeric(SCS), 
                      Marker = as.factor(genotypes[,x]), 
                      Samplingseason = as.factor(seasonsampling), 
                      DIM = as.numeric(DIM),
                      Calvingseason = as.factor(seasoncalving), 
                      Calvingyear = as.factor(yearcalving),
                      Animal = as.factor(animal),
                      Afc = as.factor(AFC),
                      Father= as.factor(father))
  
  # Remove rows with missing data
  hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
  if (length(hasMissing) > 0) { mdata <- mdata[-hasMissing,] }
  
  
  
  # Run a null model (all significant and suggestive covariates) and another model including the marker
  null.model <- lmer(Y ~ Samplingseason  + Calvingyear + DIM  + (1|Father) + (1|Animal), data = data.frame(mdata), REML = FALSE)
  markermodel <- lmer(Y ~ Samplingseason + Calvingyear + DIM  + (1|Father) + Marker + (1|Animal), data = data.frame(mdata), REML = FALSE)  
  anova(null.model, markermodel)
  markermodel
  
  # Pvalue for the model
  pval <- as.numeric(na.omit(anova(markermodel, null.model)[, "Pr(>Chisq)"]))
  pvals <- c(pvals, pval)
  
  
}

names(pvals) <- colnames(genotypes)

pvals

##P-values

SaS_M2      SaS_M3      HAS_M6      HAS_M3      HAS_M4      HAS_M8      HAS_M9      HAS_M1      HAS_M7     HAS_M10 
0.443984531 0.008623234 0.672925113 0.493416100 0.025288883 0.043105392 0.162732612 0.249753564 0.132479829 0.368760002 

####P-VALUE correction
round(p.adjust(c(0.443984531, 0.008623234, 0.672925113, 0.493416100, 0.025288883, 0.043105392, 0.162732612, 0.249753564, 0.132479829, 0.368760002), method = "BH"), 4)

[1] 0.5482 0.0862 0.6729 0.5482 0.1264 0.1437 0.3255 0.4163 0.3255 0.5268



###Association analysis for SCS3

#loding dataset
mmatrix2<- mmatrix1[inL3, ]
# Get the significant covariates  in model building
animal <- mmatrix2[, "AnimalID"]
seasonsampling <- mmatrix2[, "Samplingseason"]
seasoncalving<- mmatrix2[, "Calvingseason"]
yearcalving<- mmatrix2[, "Calvingyear"]
DIM<- mmatrix2[, "DIM"]
firstCalf <- mmatrix2[, "FirstCalfIndays"]
AFC<- round(firstCalf/30) ## 21-37
father<- mmatrix2[, "Father"]
SCS<- mmatrix2[, "SCS"]
#Loading SNP data
genotypes <- mmatrix2[,11:20]
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
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata <- data.frame(Y = as.numeric(SCS), 
                      Marker = as.factor(genotypes[,x]), 
                      Samplingseason = as.factor(seasonsampling), 
                      DIM = as.numeric(DIM),
                      Calvingseason = as.factor(seasoncalving), 
                      Calvingyear = as.factor(yearcalving),
                      Animal = as.factor(animal),
                      Afc = as.factor(AFC),
                      Father= as.factor(father))
  
  # Remove rows with missing data
  hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
  if (length(hasMissing) > 0) { mdata <- mdata[-hasMissing,] }
  
  
  
  # Run a null model (all significant and suggestive covariates) and another model including the marker
  null.model <- lmer(Y ~ Samplingseason  + Calvingyear + DIM  + (1|Father) + (1|Animal), data = data.frame(mdata), REML = FALSE)
  markermodel <- lmer(Y ~ Samplingseason + Calvingyear + DIM  + (1|Father) + Marker + (1|Animal), data = data.frame(mdata), REML = FALSE)  
  anova(null.model, markermodel)
  markermodel
  
  # Pvalue for the model
  pval <- as.numeric(na.omit(anova(markermodel, null.model)[, "Pr(>Chisq)"]))
  pvals <- c(pvals, pval)
  
  
}

names(pvals) <- colnames(genotypes)

pvals

##P-values

SaS_M2    SaS_M3    HAS_M6    HAS_M3    HAS_M4    HAS_M8    HAS_M9    HAS_M1    HAS_M7   HAS_M10 
0.8081961 0.2252467 0.9656387 0.6202122 0.1437456 0.6419212 0.6251156 0.6519740 0.5859615 0.8742098 

####P-VALUE correction
round(p.adjust(c(0.8081961, 0.2252467, 0.9656387, 0.6202122, 0.1437456, 0.6419212, 0.6251156, 0.6519740, 0.5859615, 0.8742098 ), method = "BH"), 4)

[1] 0.9656 0.9314 0.9656 0.9314 0.9314 0.9314 0.9314 0.9314 0.9314 0.9656






