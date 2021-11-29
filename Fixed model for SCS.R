
##Build SCS dataset
library(foreign)
setwd("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/")
#load genotyping data
gts<- read.table("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep="\t",  na.strings = c("","-" ,"NA","no"), fill = TRUE,header = TRUE,colClasses = "character")
#load SCS data
mlp<- read.dbf("/Users/gaojunling/Documents/R Practising/Reanalying for Dedelow/MLP.DBF")
#load herd data
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
               as.character(mlp[r, "ZELLZAHL"]),
               as.character(mlp[r, "MT"]) )
  idx = which(as.character(gts[, "ID"]) == as.character(mlp[r, "OHR"]))
  mlpdata <- c(mlpdata, as.character(gts[idx[1], -1]))
  
  mmatrix <- rbind(mmatrix, mlpdata)
  
}


colnames(mmatrix) <- c("AnimalID", "Birthdate", "FirstCalfIndays","calvingyear","Father","Sampledate", "Lactation",  "SCC", "DIM",colnames(gts)[-1])
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
  
  #Improve the association by removing fathers which do not have 10 offspring
  offspring <- c()
  for(f in unique(mdata[, "Father"])){
    nOffspring = length(unique(mdata[which(mdata[, "Father"] == f), "Animal"]))
    offspring <- c(offspring, nOffspring)
  }
  names(offspring) <- unique(mdata[, "Father"])
  
  enoughoffspring <- names(offspring)[which(offspring >= 10)]
  
  mdata <- mdata[which(mdata[, "Father"] %in% enoughoffspring),]
  
  
  
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
  
  model6 <- lmer(Y ~ Father + (1|Animal), data = mdata, REML = FALSE)
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

Lactation: 1, DIM: 4.643452e-18, Samplingseason: 0.006731092, Calvingseason: 9.991956e-05, Calvingyear: 6.545157e-06, Afc: 0.1380045, Father: 4.210681e-13 
Lactation: 2, DIM: 2.067163e-171, Samplingseason: 3.975553e-05, Calvingseason: 0.8602455, Calvingyear: 3.856069e-07, Afc: 0.4520431, Father: 0.0002501813 
Lactation: 3, DIM: 3.398343e-117, Samplingseason: 1.516974e-09, Calvingseason: 0.09883874, Calvingyear: 3.258244e-08, Afc: 0.7202842, Father: 4.555995e-05 
  
  
  
###Association analysis for SCS1

#loding dataset
mmatrix2<- mmatrix1[inL1, ]
# Get the significant covariates  in model building
animal <- mmatrix2[, "AnimalID"]
seasonsampling <- mmatrix2[, "Samplingseason"]
seasoncalving<- mmatrix2[, "Calvingseason"]
yearcalving<- mmatrix2[, "Calvingyear"]
DIM<- mmatrix2[, "DIM"]
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
  null.model <- lmer(Y ~ Samplingseason   + Calvingseason + Calvingyear + DIM + Father + (1|Animal), data = data.frame(mdata), REML = FALSE)
  markermodel <- lmer(Y ~ Samplingseason  + Calvingseason + Calvingyear + DIM + Father + Marker + (1|Animal), data = data.frame(mdata), REML = FALSE)  
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
0.417229751 0.802025336 0.008989096 0.021668756 0.450427123 0.462721201 0.276498043 0.111182343 0.931168978 0.832494718 

####P-VALUE correction
round(p.adjust(c(0.417229751, 0.802025336, 0.008989096, 0.021668756, 0.450427123, 0.462721201, 0.276498043, 0.111182343, 0.931168978, 0.832494718 ), method = "BH"), 4)

[1] 0.6610 0.9250 0.0899 0.1083 0.6610 0.6610 0.6610 0.3706 0.9312 0.9250



###Association analysis for SCS2

#loding dataset
mmatrix2<- mmatrix1[inL2, ]
# Get the significant covariates  in model building
animal <- mmatrix2[, "AnimalID"]
seasonsampling <- mmatrix2[, "Samplingseason"]

yearcalving<- mmatrix2[, "Calvingyear"]
DIM<- mmatrix2[, "DIM"]
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
                   
                      Calvingyear = as.factor(yearcalving),
                      Animal = as.factor(animal),
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
  null.model <- lmer(Y ~ Samplingseason   + Calvingyear + DIM + Father + (1|Animal), data = data.frame(mdata), REML = FALSE)
  markermodel <- lmer(Y ~ Samplingseason  + Calvingyear + DIM + Father + Marker + (1|Animal), data = data.frame(mdata), REML = FALSE)  
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
0.443439583 0.009934231 0.687446697 0.590051913 0.014680796 0.052267701 0.294685042 0.060805958 0.143716957 0.409548618 

####P-VALUE correction
round(p.adjust(c(0.443439583, 0.009934231, 0.687446697, 0.590051913, 0.014680796, 0.052267701, 0.294685042, 0.060805958, 0.143716957, 0.409548618  ), method = "BH"), 4)

[1] 0.5543 0.0734 0.6874 0.6556 0.0734 0.1520 0.4911 0.1520 0.2874 0.5543




###Association analysis for SCS3

#loding dataset
mmatrix2<- mmatrix1[inL3, ]
# Get the significant covariates  in model building
animal <- mmatrix2[, "AnimalID"]
seasonsampling <- mmatrix2[, "Samplingseason"]
seasoncalving<- mmatrix2[, "Calvingseason"]
yearcalving<- mmatrix2[, "Calvingyear"]
DIM<- mmatrix2[, "DIM"]
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
  null.model <- lmer(Y ~ Samplingseason   + Calvingseason + Calvingyear + DIM + Father + (1|Animal), data = data.frame(mdata), REML = FALSE)
  markermodel <- lmer(Y ~ Samplingseason  + Calvingseason + Calvingyear + DIM + Father + Marker + (1|Animal), data = data.frame(mdata), REML = FALSE)  
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
0.87058830 0.32827470 0.97253815 0.80646058 0.09571419 0.88286953 0.27173677 0.80835907 0.46334173 0.61968680 

####P-VALUE correction
round(p.adjust(c(0.87058830, 0.32827470, 0.97253815, 0.80646058, 0.09571419, 0.88286953, 0.27173677, 0.80835907, 0.46334173, 0.61968680 ), method = "BH"), 4)

[1] 0.9725 0.9725 0.9725 0.9725 0.9571 0.9725 0.9725 0.9725 0.9725 0.9725












