bacteria<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Pathogen data of Dedelow.txt", sep = "\t",header = TRUE)
library(foreign)
bestand<- read.dbf("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/BESTAND.DBF")
#order by date
ordering<- order(bestand[, "GEBURT"])
stalldata<- bestand[ordering, ]
stalldata<- cbind(stalldata, "year" = unlist(lapply(strsplit(as.character(stalldata[, "GEBURT"]),"-" ),"[", 1)))
#select birthyear > 2017
birthrange<- which(as.numeric(as.character(stalldata[, "year"])) > 2007)
stalldata<- stalldata[birthrange,]
#stallnumber
stallnumer<-bacteria[, "Stallnummer"]
stalldata<- stalldata[which(stalldata[, "STALL"] %in% stallnumer), ]
Animals<- stalldata[, "STALL"]
bacteria<- bacteria[which(bacteria[, "Stallnummer"] %in% Animals), ]
####
pathogendata<- NULL
for (r in 1:nrow(bacteria)){
  stall<- bacteria[r, "Stallnummer"]
  install<-which(as.character(stalldata[, "STALL"]) == as.character(stall))
  AnimalID<-NA
  if(length(install) ==1){
    AnimalID<- as.character(stalldata[install, "OHR"])
  }
  bacteriadata<- c(AnimalID,
                   as.character(bacteria[r, "Stallnummer"]),
                   as.character(bacteria[r, "Betrieb"]),
                   as.character(bacteria[r, "bv"]),
                   as.character(bacteria[r, "KNS"]),
                   as.character(bacteria[r, "Sc_pl"]),
                   as.character(bacteria[r, "Sc_min"]),
                   as.character(bacteria[r, "STA"]),
                   as.character(bacteria[r, "Ecoli"]),
                   as.character(bacteria[r, "Galt"]))
                   pathogendata<-rbind(pathogendata, bacteriadata)
  
  
}
colnames(pathogendata)<- c("AnimalID", "STALL", "Farm", "Bacteria", "CNS", "SC+", "SC-", "STA", "E.coli", "GALT")
rownames(pathogendata)<- 1:nrow(pathogendata)


####Remove NA Animal ID
pathogendata<-pathogendata[-which(is.na(pathogendata[, "AnimalID"])), ]
###Gene matrix
gts<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt", sep="\t",  na.strings = c("","-" ,"NA","no"), fill = TRUE,header = TRUE,colClasses = "character")
gts<- gts[, -1]
gts<- gts[which(gts[,"ID"] %in% pathogendata[, "AnimalID"]),]
pathogendata<- pathogendata[which(pathogendata[, "AnimalID"] %in% gts[, "ID"]), ]
genematrix<- NULL
for (r in 1:nrow(gts)) {
  ohr<- gts[r, "ID"]
  inPathogen<- which(as.character(pathogendata[, "AnimalID"]) == as.character(ohr))
  if(length(inPathogen)== 1){
    Bacteria<- as.character(pathogendata[inPathogen, "Bacteria"])
    CNS<- as.character(pathogendata[inPathogen, "CNS"])
    Scplus<- as.character(pathogendata[inPathogen, "SC+"])
    Scmin<- as.character(pathogendata[inPathogen, "SC-"])
    STA<- as.character(pathogendata[inPathogen, "STA"])
    Ecoli<- as.character(pathogendata[inPathogen, "E.coli"])
    Galt<- as.character(pathogendata[inPathogen, "GALT"])
    
  }
  dataset<- c(as.character(gts[r, "ID"]), Bacteria,CNS, Scplus, Scmin, STA, Ecoli, Galt,
              as.character(gts[r, "SaS_M2"]),
              as.character(gts[r, "SaS_M3"]),
              as.character(gts[r, "HAS_M6"]),
              as.character(gts[r, "HAS_M3"]),
              as.character(gts[r, "HAS_M4"]),
              as.character(gts[r, "HAS_M8"]),
              as.character(gts[r, "HAS_M9"]),
              as.character(gts[r, "HAS_M1"]),
              as.character(gts[r, "HAS_M7"]),
              as.character(gts[r, "HAS_M10"]))
  
  genematrix<- rbind(genematrix, dataset)
}
colnames(genematrix)<- c("AnimalID", "Bacteria","CNS", "SC+", "SC-", "STA", "E.coli", "GALT", "CHR6-1", "CHR6-2","CHR13-1", "CHR13-2", "CHR13-3", "CHR19-1", "CHR19-2", "CHR5","CHR18", "CHRX")
rownames(genematrix)<- 1:nrow(genematrix)


### pathogen cases
##CNS
##CHR5
chr5<- genematrix[, c(8,16)]
AA<- chr5[which(chr5[, "CHR5"] == "GG"), ]
AB<- chr5[which(chr5[, "CHR5"] == "AG"), ]
BB<- chr5[which(chr5[, "CHR5"] == "AA"), ]
table(AA[, "GALT"])
table(AB[, "GALT"])
table(BB[, "GALT"])
##CHR6-1
chr61<- genematrix[, c(8,9)]
AA<- chr61[which(chr61[, "CHR6-1"] == "CC"), ]
AB<- chr61[which(chr61[, "CHR6-1"] == "CT"), ]
BB<- chr61[which(chr61[, "CHR6-1"] == "TT"), ]
table(AA[, "GALT"])
table(AB[, "GALT"])
table(BB[, "GALT"])

##CHR6-2
  chr62<- genematrix[, c(8,10)]
  AA<- chr62[which(chr62[, "CHR6-2"] == "GG"), ]
  AB<- chr62[which(chr62[, "CHR6-2"] == "AG"), ]
  BB<- chr62[which(chr62[, "CHR6-2"] == "AA"), ]
  table(AA[, "GALT"])
  table(AB[, "GALT"])
  table(BB[, "GALT"])
##CHR13-1
chr131<- genematrix[, c(8,11)]
AA<- chr131[which(chr131[, "CHR13-1"] == "TT"), ]
AB<- chr131[which(chr131[, "CHR13-1"] == "CT"), ]
BB<- chr131[which(chr131[, "CHR13-1"] == "CC"), ]
table(AA[, "GALT"])
table(AB[, "GALT"])
table(BB[, "GALT"])
##CHR13-2
chr132<- genematrix[, c(8,12)]
AA<- chr132[which(chr132[, "CHR13-2"] == "GG"), ]
AB<- chr132[which(chr132[, "CHR13-2"] == "AG"), ]
BB<- chr132[which(chr132[, "CHR13-2"] == "AA"), ]
table(AA[, "GALT"])
table(AB[, "GALT"])
table(BB[, "GALT"])
##CHR13-3
  chr133<- genematrix[, c(8,13)]
  AA<- chr133[which(chr133[, "CHR13-3"] == "CC"), ]
  AB<- chr133[which(chr133[, "CHR13-3"] == "CT"), ]
  BB<- chr133[which(chr133[, "CHR13-3"] == "TT"), ]
  table(AA[, "GALT"])
  table(AB[, "GALT"])
  table(BB[, "GALT"])
##CHR18
chr18<- genematrix[, c(8,17)]
AA<- chr18[which(chr18[, "CHR18"] == "GG"), ]
AB<- chr18[which(chr18[, "CHR18"] == "GT"), ]
BB<- chr18[which(chr18[, "CHR18"] == "TT"), ]
table(AA[, "GALT"])
table(AB[, "GALT"])
table(BB[, "GALT"])
##CHR19-1
chr191<- genematrix[, c(8,14)]
AA<- chr191[which(chr191[, "CHR19-1"] == "AA"), ]
AB<- chr191[which(chr191[, "CHR19-1"] == "AG"), ]
BB<- chr191[which(chr191[, "CHR19-1"] == "GG"), ]
table(AA[, "GALT"])
table(AB[, "GALT"])
table(BB[, "GALT"])
##CHR19-2
chr192<- genematrix[, c(8,15)]
AA<- chr192[which(chr192[, "CHR19-2"] == "TT"), ]
AB<- chr192[which(chr192[, "CHR19-2"] == "CT"), ]
BB<- chr192[which(chr192[, "CHR19-2"] == "CC"), ]
table(AA[, "GALT"])
table(AB[, "GALT"])
table(BB[, "GALT"])
##CHRx
  chrx<- genematrix[, c(8,18)]
  AA<- chrx[which(chrx[, "CHRX"] == "CC"), ]
  AB<- chrx[which(chrx[, "CHRX"] == "CT"), ]
  BB<- chrx[which(chrx[, "CHRX"] == "TT"), ]
  table(AA[, "GALT"])
  table(AB[, "GALT"])
  table(BB[, "GALT"])

####plot pathogen distribution across genotype class
##Bacteria
mydata<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Pathogen across genotype classes/Bacteria.txt", sep = "\t",header = TRUE)
mydata<- data.frame(SNP= as.character(mydata[, "SNP"]), Cases= as.numeric(mydata[, "Cases"]), Geno= as.character(mydata[, "Geno"]))

ggplot(mydata, aes(x = SNP, y= Cases, fill= Geno))+
  geom_bar(stat = "identity", colour = "black")+
  scale_fill_brewer(palette = "Set3")+
  ylim(0, 15)+
  ggtitle("Bacteria cases across genotype class of significant SNP for CM") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))


###CNS
mydata1<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Pathogen across genotype classes/CNS.txt", sep = "\t",header = TRUE)
mydata1<- data.frame(SNP= as.character(mydata1[, "SNP"]), Cases= as.numeric(mydata1[, "Cases"]), Geno= as.character(mydata1[, "Geno"]))

ggplot(mydata1, aes(x = SNP, y= Cases, fill= Geno))+
  geom_bar(stat = "identity", colour = "black")+
  scale_fill_brewer(palette = "Set3")+
  ylim(0, 20)+
  ggtitle("CNS cases across genotype class of significant SNP for CM") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
###SC+
mydata1<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Pathogen across genotype classes/SC+.txt", sep = "\t",header = TRUE)
mydata1<- data.frame(SNP= as.character(mydata1[, "SNP"]), Cases= as.numeric(mydata1[, "Cases"]), Geno= as.character(mydata1[, "Geno"]))

ggplot(mydata1, aes(x = SNP, y= Cases, fill= Geno))+
  geom_bar(stat = "identity", colour = "black")+
  scale_fill_brewer(palette = "Set3")+
  ylim(0, 100)+
  ggtitle("SC+ cases across genotype class of significant SNP for CM") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

###SC-
mydata1<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Pathogen across genotype classes/SC-.txt", sep = "\t",header = TRUE)
mydata1<- data.frame(SNP= as.character(mydata1[, "SNP"]), Cases= as.numeric(mydata1[, "Cases"]), Geno= as.character(mydata1[, "Geno"]))

ggplot(mydata1, aes(x = SNP, y= Cases, fill= Geno))+
  geom_bar(stat = "identity", colour = "black")+
  scale_fill_brewer(palette = "Set3")+
  ylim(0, 6)+
  ggtitle("SC- cases across genotype class of significant SNP for CM") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

###STA
mydata1<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Pathogen across genotype classes/STA.txt", sep = "\t",header = TRUE)
mydata1<- data.frame(SNP= as.character(mydata1[, "SNP"]), Cases= as.numeric(mydata1[, "Cases"]), Geno= as.character(mydata1[, "Geno"]))

ggplot(mydata1, aes(x = SNP, y= Cases, fill= Geno))+
  geom_bar(stat = "identity", colour = "black")+
  scale_fill_brewer(palette = "Set3")+
  ylim(0, 15)+
  ggtitle("STA cases across genotype class of significant SNP for CM") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

###E.coli
mydata1<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Pathogen across genotype classes/E.coli.txt", sep = "\t",header = TRUE)
mydata1<- data.frame(SNP= as.character(mydata1[, "SNP"]), Cases= as.numeric(mydata1[, "Cases"]), Geno= as.character(mydata1[, "Geno"]))

ggplot(mydata1, aes(x = SNP, y= Cases, fill= Geno))+
  geom_bar(stat = "identity", colour = "black")+
  scale_fill_brewer(palette = "Set3")+
  ylim(0, 40)+
  ggtitle("E.coli cases across genotype class of significant SNP for CM") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

####GALT
mydata1<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Pathogen across genotype classes/GALT.txt", sep = "\t",header = TRUE)
mydata1<- data.frame(SNP= as.character(mydata1[, "SNP"]), Cases= as.numeric(mydata1[, "Cases"]), Geno= as.character(mydata1[, "Geno"]))

ggplot(mydata1, aes(x = SNP, y= Cases, fill= Geno))+
  geom_bar(stat = "identity", colour = "black")+
  scale_fill_brewer(palette = "Set3")+
  ylim(0, 15)+
  ggtitle("GALT cases across genotype class of significant SNP for CM") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))






