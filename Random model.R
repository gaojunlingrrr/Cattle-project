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

SCS <- log2(as.numeric(as.character(mmatrix1[, "SCC"])) / 100) + 3
SCS[which(!is.finite(SCS))] <- 0

mmatrix1<- cbind(mmatrix1, SCS)

##
Calvingyear<- unlist(lapply(strsplit(as.character(mmatrix1[,"calvingyear"]), "-"), "[", 1))

mmatrix1<- cbind(mmatrix1, Calvingyear)



## Divided into lactation 1, 2, 3
mmatrix1<- mmatrix1[which(mmatrix1[, "Lactation"] == 1), ]
mmatrix1<- mmatrix1[which(mmatrix1[, "Lactation"] == 2), ]
mmatrix1<- mmatrix1[which(mmatrix1[, "Lactation"] == 3), ]

# Get the covariates we want to investigate in model building
animal <- mmatrix1[, "AnimalID"]
lactation <- mmatrix1[, "Lactation"]
birthyear <- mmatrix1[, "birthyear"]
seasonsampling <- mmatrix1[, "Samplingseason"]
seasoncalving<- mmatrix1[, "Calvingseason"]
yearcalving<- mmatrix1[, "Calvingyear"]
firstCalf <- mmatrix1[, "FirstCalfIndays"]
AFC<- round(firstCalf/30) ## 21-37
DIM<- mmatrix1[, "DIM"]
father<- mmatrix1[, "Father"]
SCS<- mmatrix1[, "SCS"]
mmatrix1<-cbind(mmatrix1, AFC)

# Combine all data into a matrix
mdata <- data.frame(Y = as.numeric(SCS), 
                    Lact = as.factor(lactation),
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
anova(model1, null.model) # < 2.2e-16 ***

model2 <- lmer(Y ~ Samplingseason + (1|Animal), data = mdata, REML = FALSE)
anova(model2, null.model) # < 5.805e-06 ***

model3 <- lmer(Y ~ Calvingseason + (1|Animal), data = mdata, REML = FALSE)
anova(model3, null.model) # < 0.01067 *

model4 <- lmer(Y ~ Calvingyear + (1|Animal), data = mdata, REML = FALSE)
anova(model4, null.model) # < 0.0001964 ***

model5 <- lmer(Y ~ Afc + (1|Animal), data = mdata, REML = FALSE)
anova(model5, null.model) # <  NS

model6 <- lmer(Y ~ (1|Father) + (1|Animal), data = mdata, REML = FALSE)
anova(model6, null.model) # < 3.233e-12 ***


genotypes <- mmatrix1[,11:20]
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
genotypes <- mmatrix1[,11:20]
pvals <- NULL
for (x in 1:ncol(genotypes)) {
  # Combine all data into a matrix
  mdata <- data.frame(Y = as.numeric(SCS), 
                      Marker = as.factor(snp8), #order(snp8,snp1,snp2,snp3,snp4,snp5,snp9,snp6,snp7,snp10)
                      Samplingseason = as.factor(seasonsampling), 
                      DIM = as.numeric(DIM),
                      Calvingseason = as.factor(seasoncalving), 
                      Calvingyear = as.factor(yearcalving),
                      Animal = as.factor(animal),
                      Afc = as.numeric(AFC),
                      Father= as.factor(father))}

# Remove rows with missing data
hasMissing <- which(apply(apply(mdata,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata <- mdata[-hasMissing,] }


# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- lmer(Y ~ Samplingseason   + Calvingseason + Calvingyear + DIM + (1|Father) + (1|Animal), data = mdata, REML = FALSE)
markermodel <- lmer(Y ~ Samplingseason  + Calvingseason + Calvingyear + DIM + (1|Father) + Marker + (1|Animal), data = mdata, REML = FALSE)  
anova(null.model, markermodel)
markermodel

####P-VALUE correction
p.adjust(c(0.1091,
           0.4004,
           0.3417,
           0.188,
           0.0964,
           0.0701,
           0.9931,
           0.2819,
           0.6691,
           0.0524), method = "BH")









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

###CM1, CM2, CM3
cmdata<- cmdata[which(cmdata[, "Lactation"] ==1), ]
cmdata<- cmdata[which(cmdata[, "Lactation"] ==2), ]
cmdata<- cmdata[which(cmdata[, "Lactation"] ==3), ]

# Get the covariates we want to investigate in model building
animal <- cmdata[, "AnimalID"]
lactation <- cmdata[, "Lactation"]
seasoncalving<- cmdata[, "Calvingseason"]
yearcalving <-  cmdata[, "Calvingyear_N"]
firstCalf <- cmdata[, "FirstCalfIndays"]
AFC<- round(firstCalf/30) #21-37
father<- cmdata[, "Father"]
animalcode<- cmdata[, "CM"] 
cmdata<- cbind(cmdata, AFC)

##combine all data into a matrix


mdata1 <- data.frame(Y = as.numeric(animalcode), 
                     Lact = as.factor(lactation),
                     Animal = as.factor(animal),
                     Calvingseason = as.factor(seasoncalving),
                     Calvingyear = as.factor(yearcalving),
                     Firstcalf = as.numeric(AFC),#21-37
                     Father= as.factor(father)
)



# Remove rows with missing data
mdata1<- na.omit(mdata1)





library(lme4)
null.model <- glmer(Y ~ (1|Animal), data = mdata1, family = binomial,nAGQ = 0)


model2 <- glmer(Y ~ Calvingseason + (1|Animal), data = mdata1,family = binomial, nAGQ = 0)　　　
anova(model2, null.model, test="Chisq") # < 4.833e-05 ***

model3 <- glmer(Y ~ Calvingyear + (1|Animal), data = mdata1, family = binomial, nAGQ = 0)
anova(model3, null.model, test="Chisq") # < 2.2e-16 ***

model4 <- glmer(Y ~ Firstcalf + (1|Animal), data = mdata1, family = binomial, nAGQ = 0)
anova(model4, null.model, test="Chisq") # < 1.839e-12 ***

model5 <- glmer(Y ~ (1|Father) + (1|Animal), data = mdata1, family = binomial, nAGQ = 0)
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
                     Marker = as.factor(snp10), #order(snp8,snp1,snp2,snp3,snp4,snp5,snp9,snp6,snp7,snp10)
                     Animal = as.factor(animal),
                     Lact = as.factor(lactation),
                     Calvingseason = as.factor(seasoncalving),                 
                     Calvingyear = as.factor(yearcalving),
                     Afc = as.numeric(AFC),
                     Father= as.factor(father))

# Remove rows with missing data
hasMissing <- which(apply(apply(mdata2,1,is.na),2,sum) != 0)
if (length(hasMissing) > 0) { mdata2 <- mdata2[-hasMissing,] }



# Run a null model (all significant and suggestive covariates) and another model including the marker
null.model <- glmer(Y ~ Calvingseason + Calvingyear   + (1|Father) + (1|Animal), data = mdata2, family = binomial, nAGQ=0)
markermodel <- glmer(Y ~ Calvingseason + Calvingyear  + (1|Father) + Marker+ (1|Animal) , data = mdata2, family = binomial,nAGQ=0)  
anova(null.model, markermodel)
markermodel


####P-VALUE correction
p.adjust(c(0.1091,
           0.4004,
           0.3417,
           0.188,
           0.0964,
           0.0701,
           0.9931,
           0.2819,
           0.6691,
           0.0524), method = "BH")



