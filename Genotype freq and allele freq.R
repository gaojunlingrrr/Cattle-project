genodata<- read.table("/Users/Pomelo/Documents/R Practising/Reanalying for Dedelow/Myanimal in Dedelow/newdata/Genotyping data- Dedelow.txt",sep = "\t",  na.strings = c(" ","-" ,"NA","no"),fill = TRUE, header = TRUE)
##SNP1
TT<- subset(genodata, SaS_M2== "TT")
CT<- subset(genodata, SaS_M2== "CT")
CC<- subset(genodata, SaS_M2== "CC")
N11<-nrow(TT)
N12<-nrow(CT)
N22<-nrow(CC)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q


###SNP2
AA<- subset(genodata, SaS_M3== "AA")
AG<- subset(genodata, SaS_M3== "AG")
GG<- subset(genodata, SaS_M3== "GG")
N11<-nrow(AA)
N12<-nrow(AG)
N22<-nrow(GG)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q

###SNP3
TT<- subset(genodata, HAS_M6== "TT")
CT<- subset(genodata, HAS_M6== "CT")
CC<- subset(genodata, HAS_M6== "CC")
N11<-nrow(TT)
N12<-nrow(CT)
N22<-nrow(CC)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q
###SNP4
AA<- subset(genodata, HAS_M3== "AA")
AG<- subset(genodata, HAS_M3== "AG")
GG<- subset(genodata, HAS_M3== "GG")
N11<-nrow(AA)
N12<-nrow(AG)
N22<-nrow(GG)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q
###SNP5
TT<- subset(genodata, HAS_M4== "TT")
CT<- subset(genodata, HAS_M4== "CT")
CC<- subset(genodata, HAS_M4== "CC")
N11<-nrow(TT)
N12<-nrow(CT)
N22<-nrow(CC)
N<- N11+N12+N22
N11/N
N12/N
N22/N
N11
N12
N22
p<- (N11 + N12/2)/N
q<- (N22 + N12/2)/N
p
q
###SNP6

  AA<- subset(genodata, HAS_M8== "AA")
  AG<- subset(genodata, HAS_M8== "AG")
  GG<- subset(genodata, HAS_M8== "GG")
  N11<-nrow(AA)
  N12<-nrow(AG)
  N22<-nrow(GG)
  N<- N11+N12+N22
  N11/N
  N12/N
  N22/N
  N11
  N12
  N22
  p<- (N11 + N12/2)/N
  q<- (N22 + N12/2)/N
  p
  q
  
  ###SNP7
  TT<- subset(genodata, HAS_M9== "TT")
  CT<- subset(genodata, HAS_M9== "CT")
  CC<- subset(genodata, HAS_M9== "CC")
  N11<-nrow(TT)
  N12<-nrow(CT)
  N22<-nrow(CC)
  N<- N11+N12+N22
  N11/N
  N12/N
  N22/N
  N11
  N12
  N22
  p<- (N11 + N12/2)/N
  q<- (N22 + N12/2)/N
  p
  q
###SNP8
  
  AA<- subset(genodata, HAS_M1== "AA")
  AG<- subset(genodata, HAS_M1== "AG")
  GG<- subset(genodata, HAS_M1== "GG")
  N11<-nrow(AA)
  N12<-nrow(AG)
  N22<-nrow(GG)
  N<- N11+N12+N22
  N11/N
  N12/N
  N22/N
  N11
  N12
  N22
  p<- (N11 + N12/2)/N
  q<- (N22 + N12/2)/N
  p
  q
###SNP9

  TT<- subset(genodata, HAS_M7== "TT")
  GT<- subset(genodata, HAS_M7== "GT")
  GG<- subset(genodata, HAS_M7== "GG")
  N11<-nrow(TT)
  N12<-nrow(GT)
  N22<-nrow(GG)
  N<- N11+N12+N22
  N11/N
  N12/N
  N22/N
  N11
  N12
  N22 
  p<- (N11 + N12/2)/N
  q<- (N22 + N12/2)/N
  p
  q
###SNP10
  
  TT<- subset(genodata, HAS_M10== "TT")
  CT<- subset(genodata, HAS_M10== "CT")
  CC<- subset(genodata, HAS_M10== "CC")
  N11<-nrow(TT)
  N12<-nrow(CT)
  N22<-nrow(CC)
  N<- N11+N12+N22
  N11/N
  N12/N
  N22/N
  N11
  N12
  N22  
  p<- (N11 + N12/2)/N
  q<- (N22 + N12/2)/N
  p
  q







N22/N