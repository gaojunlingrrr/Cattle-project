##install gprofiler2 package
install.packages("gprofiler2")
library(gprofiler2)

# installing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
BiocManager::install(c("DESeq2", "airway"))
library(DESeq2) 
library(airway) 
library(gprofiler2)
data(airway)