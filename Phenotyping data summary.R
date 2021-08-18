mysample<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep = "\t", fill = TRUE, header = TRUE)
mysample <- mysample[, -c(13,14)]
myanimals<- mysample[,"ID"]
setwd("/Users/Pomelo/Documents/R Practising/Dedelow/")
mydata <- read.table("krankData-Dedelow.txt", sep = "\t",fill = TRUE, header = TRUE)
ismastitis = which(grepl("1.13", mydata[, "STAUFEN"], fixed = TRUE))
Mastitisdata = mydata[ismastitis, ]
Mastitisdata= Mastitisdata[which(Mastitisdata[, "OHR"] %in% myanimals), ]
# Covert dates to R format, and store as new column
Rdaterows = as.Date(as.character(Mastitisdata[, "DATUM"]), format="%y/%m/%d")
Mastitisdata = cbind(Mastitisdata, "Rdate"= Rdaterows)

#order by date
ordering = order(Mastitisdata[,"Rdate"])
Mastitisdata = Mastitisdata[ordering, ]
Mastitisdata = cbind(Mastitisdata, "year" = unlist(lapply(strsplit(as.character(Mastitisdata[,"Rdate"]), "-"), "[", 1)))

# Select all the entries  between 2009 and 2017  (the average productive life of a Holstein is approximately four years).
Timerange = which(as.numeric(as.character(Mastitisdata[, "year"])) > 2009 & as.numeric(as.character(Mastitisdata[, "year"])) < 2018 )
Mastitisdata = Mastitisdata[Timerange, ]

in2010 <- which(grepl("2010", Mastitisdata[,"Rdate"]))
in2011 <- which(grepl("2011", Mastitisdata[,"Rdate"]))
in2012 <- which(grepl("2012", Mastitisdata[,"Rdate"]))
in2013 <- which(grepl("2013", Mastitisdata[,"Rdate"]))
in2014 <- which(grepl("2014", Mastitisdata[,"Rdate"]))
in2015 <- which(grepl("2015", Mastitisdata[,"Rdate"]))
in2016 <- which(grepl("2016", Mastitisdata[,"Rdate"]))
in2017 <- which(grepl("2017", Mastitisdata[,"Rdate"]))

generateInfo = function(Mastitisdata){
  Mastitistable = table(as.character(Mastitisdata[, "OHR"]))
  MastitisInfo <- NULL
  Mastitisentry <- vector("list", length(names(Mastitistable)))
  names(Mastitisentry) = names(Mastitistable)
  for (cow in names(Mastitistable)){
    cowID = which(Mastitisdata[, "OHR"] == cow)
    Mastitisentry[[cow]] = Mastitisdata[cowID, "Rdate"]
    nDiagnosis = length(Mastitisentry[[cow]])
    Diagnosisdate = Mastitisentry[[cow]][1]
    if (nDiagnosis > 1){
      Diagnosis = c("N", rep(NA, nDiagnosis-1))
      for(x in 2:nDiagnosis){
        if ((Mastitisentry[[cow]][x] - Diagnosisdate) < 14){
          Diagnosis[x] <- "F"
        }
        else{
          Diagnosis[x] <- "N"
        }
        Diagnosisdate = Mastitisentry[[cow]][x]
      }
      Mastitisentry[[cow]] = rbind(as.character(Mastitisentry[[cow]]), Diagnosis)
      nMastitis = length(which(Diagnosis =="N"))
      MastitisInfo = rbind(MastitisInfo, c(cow,nDiagnosis,nMastitis))
    }else{
      MastitisInfo = rbind(MastitisInfo, c(cow,1,1))
    }
  }
  
  colnames(MastitisInfo) <- c("OHR", "nDiagnosis", "nMastitis")
  return(MastitisInfo)
}
MastitisInfo = generateInfo(Mastitisdata)
sum(as.numeric(MastitisInfo[,"nDiagnosis"])) / nrow(MastitisInfo)
sum(as.numeric(MastitisInfo[,"nMastitis"])) / nrow(MastitisInfo)

year = 2010
for(individuals in list(in2010,in2011,in2012,in2013,in2014,in2015,in2016,in2017)){
  MastitisInfo = generateInfo(Mastitisdata[individuals, ])
  cat("Year:", year)
  cat(", Cows:", nrow(MastitisInfo))
  cat(", nDiagnosis:", sum(as.numeric(MastitisInfo[, "nDiagnosis"])))
  cat(", nMastitis:", sum(as.numeric(MastitisInfo[, "nMastitis"])))
  cat(", nDiagnosis:", sum(as.numeric(MastitisInfo[, "nDiagnosis"]))/nrow(MastitisInfo))
  cat(", nMastitis:", sum(as.numeric(MastitisInfo[, "nMastitis"]))/nrow(MastitisInfo),"\n")
  year = year +1
}
###plot sick animal distribution  and mastitis frequency among years
year<- c(2010,2011,2012,2013,2014,2015,2016,2017)
Ncows<- c(24,50,103,206,311,622,996,1153)
Frequency<- c(1.38,1.22,1.26,1.44,1.43,1.58,1.83,1.62)
plot(year,Ncows,yaxt= "n", ylim= c(0, 1200),xlab = "Year", ylab = "Number of sick cows", main = "Sick animal distribution among years of Farm Dedelow ")
axis(2, seq(0,1200,100), seq(0,1200,100), las=1)
plot(year, Frequency, xlab = "Year",yaxt= "n", ylim= c(1.2,1.9),ylab = "Mastitis incidence", main = "Mastitis incidence among years of Farm Dedelow ")
axis(2, seq(1.2,1.9,0.1), seq(1.2,1.9,0.1), las=1)

#Per Lactation
inL1 <- which(Mastitisdata[,"LAKTATION"] == 1)
inL2 <- which(Mastitisdata[,"LAKTATION"] == 2)
inL3 <- which(Mastitisdata[,"LAKTATION"] == 3)
inL4 <- which(Mastitisdata[,"LAKTATION"] > 3)


lactation = 1
for(individuals in list(inL1,inL2,inL3,inL4)){
  MastitisInfo = generateInfo(Mastitisdata[individuals, ])
  cat("Lactation:", lactation)
  cat(", Cows", nrow(MastitisInfo))
  cat(", nDiagnois", sum(as.numeric(MastitisInfo[, "nDiagnosis"])))
  cat(", nMastitis", sum(as.numeric(MastitisInfo[, "nMastitis"])))
  cat(", nDiagnosis:", sum(as.numeric(MastitisInfo[, "nDiagnosis"]))/nrow(MastitisInfo))
  cat(", nMastitis:", sum(as.numeric(MastitisInfo[, "nMastitis"]))/nrow(MastitisInfo),"\n")
  lactation = lactation + 1
  
}

###plot animal distribution among lactations
x<- c(907, 920, 661, 401)
label<- c("Lact1", "Lact2", "Lact3", "Lact3+")
piepercentage<- round(100*x/sum(x),1)
piepercentage<- paste(piepercentage, "%", sep= "" )###paste函数用来连接字符串

pie(x, labels = piepercentage, main = "Sick cows' distribution over lactations in Dedelow", col = rainbow(length(x)))
legend("topright", label, cex=0.8, fill= rainbow(length(x)))
###
Mf<-c(1.39,1.78,1.92,3.33)
nlactation<-c(1,2,3,4)
plot(nlactation,Mf, xaxt= 'n',yaxt= 'n', ylim= c(1.3,3.5),xlab = "Lactation" , ylab = "Mastitis incidence", col="red", main = "Mastitis incidence among lactations of Farm Dedelow")
axis(1, c(1,2,3,4), labels = c("1", "2", "3","3+"))
axis(2, seq(1.3,3.5,0.1),seq(1.3,3.5,0.1), las=1 )

####Mastitis incidence per lactation per year
Mastitisdata2010<- Mastitisdata[in2017,]
MastitisdataL1<-Mastitisdata2010[inL4, ]
generateInfo = function(MastitisdataL1){
  Mastitistable = table(as.character(MastitisdataL1[, "OHR"]))
  MastitisInfo <- NULL
  Mastitisentry <- vector("list", length(names(Mastitistable)))
  names(Mastitisentry) = names(Mastitistable)
  for (cow in names(Mastitistable)){
    cowID = which(MastitisdataL1[, "OHR"] == cow)
    Mastitisentry[[cow]] = MastitisdataL1[cowID, "Rdate"]
    nDiagnosis = length(Mastitisentry[[cow]])
    Diagnosisdate = Mastitisentry[[cow]][1]
    if (nDiagnosis > 1){
      Diagnosis = c("N", rep(NA, nDiagnosis-1))  ###无效值NA, NaN主要用于应付某操作没完成、结果未知的情况
      for(x in 2:nDiagnosis){
        if ((Mastitisentry[[cow]][x] - Diagnosisdate) < 14){
          Diagnosis[x] <- "F"
        }
        else{
          Diagnosis[x] <- "N"
        }
        Diagnosisdate = Mastitisentry[[cow]][x]
      }
      Mastitisentry[[cow]] = rbind(as.character(Mastitisentry[[cow]]), Diagnosis)
      nMastitis = length(which(Diagnosis =="N"))
      MastitisInfo = rbind(MastitisInfo, c(cow,nDiagnosis,nMastitis))
    }else{
      MastitisInfo = rbind(MastitisInfo, c(cow,1,1))
    }
  }
  
  colnames(MastitisInfo) <- c("OHR", "nDiagnosis", "nMastitis")
  return(MastitisInfo)
}
MastitisInfo = generateInfo(MastitisdataL1)
sum(as.numeric(MastitisInfo[,"nDiagnosis"])) / nrow(MastitisInfo)
sum(as.numeric(MastitisInfo[,"nMastitis"])) / nrow(MastitisInfo)
sd(as.numeric(MastitisInfo[,"nMastitis"]))/sqrt(nrow(MastitisInfo))##标准误差
###   plot mastitis frequency
###Year
Mfreq <- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/MF per lactatation per year.txt", sep = "\t", header = TRUE)
ggplot(Mfreq, aes(x = Year, y= MF, fill = Lactation) ) +  ###Error: Continuous value supplied to discrete scale---as.factor(Lactation)
  geom_bar(stat= "identity", position = "dodge" , color= "black", width = 0.8) +
  scale_x_continuous(breaks = seq(2010,2017,1)) +
  scale_fill_brewer(palette = 'Accent') +
  theme(legend.key=element_blank()) +
  labs(x= "Year", y= "Mastitis incidence", title = "Mastitis incidence of Farm Dedelow") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10, color = "black"), plot.title = element_text(hjust = 0.5, size = 20))


##绘制误差棒及显著性标记
###Among years
ggplot(data = Mfreq, aes(x = Year, y = MF, fill = Lactation)) +
  
  geom_bar(stat = "identity", position = "dodge",color= "black", width = 0.8) +
  scale_x_continuous(breaks = seq(2010,2017,1)) +
  
  geom_errorbar(aes(ymax = MF + SE, ymin = MF - SE),
                
                position = position_dodge(0.9), width = 0.3) +
  
  scale_fill_brewer(palette = "Set1") +
  labs(x= "Year", y= "Mastitis incidence", title = "Mastitis incidence of Farm Dedelow") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10, color = "black"), plot.title = element_text(hjust = 0.5, size = 20))

## Among lactations
Mfreq[, "Year"] <- factor(Mfreq[, "Year"])
ggplot(data = Mfreq, aes(x = Lactation, y = MF, fill = Year)) +
  
  geom_bar(stat = "identity", position = "dodge",color= "black", width = 0.8) +
  
  geom_errorbar(aes(ymax = MF + SE, ymin = MF - SE),
                
                position = position_dodge(width = 0.8), width = 0.3) +
  
  scale_fill_brewer(palette = "Set1") +
  labs(x= "Lactation", y= "Mastitis incidence", title = "Mastitis incidence of Farm Dedelow") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10, color = "black"), plot.title = element_text(hjust = 0.5, size = 20))

### Once and multiple mastitis cases
lactation = 1
for (individuals in list(inL1, inL2, inL3, inL4)){
  MastitisInfo = generateInfo(Mastitisdata[individuals, ])
  Multicases = which(as.numeric(MastitisInfo[,"nMastitis"]) ==1)
  MastitisInfo = MastitisInfo[Multicases, ]
  cat(", Cows:",nrow(MastitisInfo))
  lactation = lactation + 1
}

a<- c(907, 920, 661, 401)
b<- c(248, 418, 351, 286)
c<- c(659, 502, 310, 115)
b/a
c/a

year = 2010
for (individuals in list(in2010, in2011, in2012, in2013, in2014, in2015, in2016, in2017)){
  MastitisInfo = generateInfo(Mastitisdata[individuals, ])
  Multicases = which(as.numeric(MastitisInfo[,"nMastitis"]) >1)
  MastitisInfo = MastitisInfo[Multicases, ]
  cat(", Cows:",nrow(MastitisInfo))
  year = year + 1
}
a<- c(24,50,103,206,311,622,996,1153)
b<- c(17,40,81,143,218,389,521,681)
c<- c(7,10,22,63,93,233,475,472)
b/a
c/a




# Make a plot for multiple mastitis cases among lactations
MutipleData <- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/once and multiple mastitis cases.txt", sep = "\t", header = TRUE)
ggplot(MutipleData, aes(x = Lactation, y = Percentage, fill = Group)) +
  
  geom_bar(position = "dodge", stat = "identity", colour= "black")+
  scale_fill_brewer(palette="Pastel1")+
  geom_text(mapping = aes(label = Percentage), vjust = 0, colour = "black", position = position_dodge(.9), size = 5)+
  labs(x= "Lactation", y= "Percentage of cows (%)", title = "Percentage of sick cows with once or multiple cases among lactations ") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10, color = "black"), plot.title = element_text(hjust = 0.5, size = 12))
###among years
MutipleData <- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/sick cow percentage among years.txt", sep = "\t", header = TRUE)
ggplot(MutipleData, aes(x = Year, y = Percentage, fill = Group)) +
  
  geom_bar(position = "dodge", stat = "identity", colour= "black")+
  scale_fill_brewer(palette="Pastel1")+
  scale_x_continuous(breaks = seq(2010,2017,1)) +
  geom_text(mapping = aes(label = Percentage), vjust = 0, colour = "black", position = position_dodge(.9), size = 5)+
  labs(x= "Year", y= "Percentage of cows (%)", title = "Percentage of sick cows with once or multiple cases among years ") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10, color = "black"), plot.title = element_text(hjust = 0.5, size = 12))







####My healthy animal
#Delete sick cows from 305 file 
mysample<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep = "\t", fill = TRUE, header = TRUE)
myanimals<- mysample[,"ID"]

lactData <- read.table("/Users/Pomelo/Documents/R Practising/Dedelow/305data-Dedelow.txt", sep = "\t", header = TRUE)
lactData<- lactData[which(lactData[, "OHR"] %in% myanimals),]
# Covert dates to R format, and store as new column
KdateObjects = as.Date(as.character(lactData[, "KALBUNG"]), format="%y/%m/%d")
lactData = cbind(lactData, "Kdate" = KdateObjects) ### lapply function return a list 
lactData = cbind(lactData, "year" = unlist(lapply(strsplit(as.character(lactData[,"Kdate"]), "-"), "[", 1)))
lactData = cbind(lactData, "month" = unlist(lapply(strsplit(as.character(lactData[,"Kdate"]), "-", fixed = TRUE), "[", 2)))

## dates to season

getSeason<- function(lactData){
  dates<- as.numeric(lactData[, "month"]) ###factor->numeric
  ret<-rep(NA, length(dates))
  ret[dates>=4 & dates<=5] <- "Spring"
  ret[dates>=6 & dates<=8] <- "Summer"
  ret[dates>=9 & dates<=10] <- "Fall"
  ret[dates==11 | dates==12 | dates==1 | dates==2 | dates==3 ] <- "Winter"
  
  return(ret)
}

Seasonobjects <- getSeason(lactData)

lactData <- cbind(lactData, "Season" =Seasonobjects )
###select cows giving birth in 2010
lact2010 <- lactData[which(lactData[ , "year"] == 2017), ]
####select sick cows in 2010 
Mastitis2010 = Mastitisdata[which(Mastitisdata[, "year"]== 2017), ]
###select healthy cows in 2010
healthy2010 = lact2010[which(!lact2010[, "OHR"] %in% Mastitis2010[, "OHR"] ),  ]
#### number of sick and healthy cows
length(unique(Mastitis2010[, "OHR"]))
length(unique(healthy2010[, "OHR"]))
length(unique(lact2010[, "OHR"]))


###animal distribution among years
a<- c(56,121,229,466,860,1367,2015,2028) ###all cows
b<- c(17,40,81,143,218,389,521,681) ###once
c<- c(7,10,22,63,93,233,475,472) ###multiple
d<- c(32,71,126,260,549,745,1019,875)###healthy
b/a
c/a
d/a

b+c+d

## animal distribution among lactations
healthyCows = lactData[which(!(lactData[,"OHR"] %in% Mastitisdata[, "OHR"])), ]
healthylact1 = healthyCows[which(healthyCows[, "LAKTATION"]==1), ]
healthylact2 = healthyCows[which(healthyCows[, "LAKTATION"]==2), ]
healthylact3 = healthyCows[which(healthyCows[, "LAKTATION"]==3), ]

length(unique(healthylact1[,"OHR"]))
length(unique(healthylact2[,"OHR"]))
length(unique(healthylact3[,"OHR"]))

a<- c(1756,1336,781,441) ### all cows
b<- c(248, 418, 351, 286) ##once
c<- c(659, 502, 310, 115) ##multiple
d<- c(849,416,120,40)##healthy
b+c+d
b/a
c/a
d/a
###plot - year

MutipleData <- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/healthy and sick cow among years.txt", sep = "\t", header = TRUE)
ggplot(MutipleData, aes(x = Year, y = Percentage, fill = Group)) +
  
  geom_bar(position = "dodge", stat = "identity", colour= "black")+
  scale_fill_brewer(palette="Pastel1")+
  scale_x_continuous(breaks = seq(2010,2017,1)) +
  geom_text(mapping = aes(label = Percentage), vjust = 0, colour = "black", position = position_dodge(.9), size = 4)+
  labs(x= "Year", y= "Percentage of cows (%)", title = "Percentage of healthy and sick cows among years ") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10, color = "black"), plot.title = element_text(hjust = 0.5, size = 15))

### plot- lactation
MutipleData <- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/healthy and sick cows among lactations.txt", sep = "\t", header = TRUE)
ggplot(MutipleData, aes(x = Lactation, y = Percentage, fill = Group)) +
  
  geom_bar(position = "dodge", stat = "identity", colour= "black")+
  scale_fill_brewer(palette="Pastel1")+
  geom_text(mapping = aes(label = Percentage), vjust = 0, colour = "black", position = position_dodge(.9), size = 4)+
  labs(x= "Lactation", y= "Percentage of cows (%)", title = "Percentage of healthy and sick cows among lactations ") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10, color = "black"), plot.title = element_text(hjust = 0.5, size = 15))

###### Least square means
lactDataformean <- lactData[which(as.numeric(as.character(lactData[, "year"])) >2009 & as.numeric(as.character(lactData[, "year"])) < 2018),  ]
lactDataformean <- lactDataformean[which(lactDataformean[, "EIWEISSKG"]!=0), ]
EIWEISSKG<- lactDataformean[, "EIWEISSKG"]
y<- EIWEISSKG

x<- lactDataformean[ , "Season"]


mmodel <- lm(y ~  x)
summary(mmodel)
anova(mmodel)
mean(y, na.rm = TRUE)  ## na.rm- delete NA value
sd<-sd(y, na.rm = TRUE)
se<- sd/sqrt(length(y)) ## stand error
se 
confint(mmodel) ### confidence interval


plot(c(0.5, 4.5), y = c(250, 400), t = 'n', xaxt = 'n', xlab = "Season", ylab = " 305 days Protein yield",main= "Protein yield per lacatation among seasons" , las=2)
mymeans = c()
mystderrs = c()
at = 1
for(season in c("Winter", "Spring", "Summer", "Fall")){
  vals <- y[which(x == season)]
  mymean = mean(vals, na.rm = TRUE)
  mystderr = sd(vals, na.rm = TRUE) #/ sqrt(length(which(!is.na(vals))))
  mymeans = c(mymeans, mymean)
  mystderrs = c(mystderrs, mystderr)
  points(c(at, at), c(mymean + mystderr, mymean - mystderr), t = 'b', pch = "-")
  text(at, 250, paste0("n=", length(vals)))
  at = at + 1
}
points(1:4, mymeans, t = "b", pch = 19)###type= "b "  在图形中数据显示为点和连接线
axis(1, 1:4, c("Winter", "Spring", "Summer", "Fall"))

### somatic cell score
library(foreign)
somaticdata<- read.dbf("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/MLP.DBF")
mysample<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep = "\t", fill = TRUE, header = TRUE)
myanimals<- mysample[,"ID"]
somaticdata= somaticdata[which(somaticdata[, "OHR"] %in% myanimals), ]
ordering <- order(somaticdata[, "DATUM"])
somaticdata<- somaticdata[ordering, ]
somaticdata<- cbind(somaticdata, "year"= unlist(lapply(strsplit(as.character(somaticdata[,"DATUM"]), "-"), "[", 1)))
somaticdata <- somaticdata[which(as.numeric(as.character(somaticdata[, "year"])) >2009 & as.numeric(as.character(somaticdata[, "year"])) < 2018),  ]
somaticdata <- somaticdata[which(somaticdata[, "ZELLZAHL"]!=0), ]
SCS<- log2(somaticdata[,"ZELLZAHL"]/100)+3
somaticdata<- cbind(somaticdata, "SCS" = SCS)
plot(sort(somaticdata[, "ZELLZAHL"], decreasing = TRUE), col= "red", type = "h", xlab = "Individuals", ylab = "SCC (*1000)", main = "Somatic cell count of all samples")
summary(somaticdata[, "ZELLZAHL"])

#### Check subclinical mastitis cases in "healthy group"
healthyCows = lactData[which(!(lactData[,"OHR"] %in% Mastitisdata[, "OHR"])), ]
healthylact1 = healthyCows[which(healthyCows[, "LAKTATION"]==1), ]
healthylact2 = healthyCows[which(healthyCows[, "LAKTATION"]==2), ]
healthylact3 = healthyCows[which(healthyCows[, "LAKTATION"]==3), ]
healthylact4 = healthyCows[which(healthyCows[, "LAKTATION"]>3), ]
Animalobjects = healthyCows[, "OHR"]

somaticdata= somaticdata[which(somaticdata[, "OHR"] %in% Animalobjects), ]
ordering <- order(somaticdata[, "DATUM"])
somaticdata<- somaticdata[ordering, ]
somaticdata<- cbind(somaticdata, "year"= unlist(lapply(strsplit(as.character(somaticdata[,"DATUM"]), "-"), "[", 1)))
somaticdata <- somaticdata[which(as.numeric(as.character(somaticdata[, "year"])) >2009 & as.numeric(as.character(somaticdata[, "year"])) < 2018),  ]
somaticdata <- somaticdata[which(somaticdata[, "ZELLZAHL"]!=0), ]
SCS<- log2(somaticdata[,"ZELLZAHL"]/100)+3
somaticdata<- cbind(somaticdata, "SCS" = SCS)
plot(sort(somaticdata[, "ZELLZAHL"], decreasing = TRUE), col= "red", type = "h", xlab = "Individuals", ylab = "SCC (*1000)", main = "Somatic cell count of healthy samples")
summary(somaticdata[, "ZELLZAHL"])
