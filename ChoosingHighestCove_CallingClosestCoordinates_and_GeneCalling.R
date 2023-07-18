rm(list = ls())
library(dplyr)
library(stringr)

#Enter the file names
Treatment <- read.csv("Chl_1Q_t60_exclusiveComplementPeaks.csv")
NDCexc <- read.csv("Chl_ndc_t60_exclusiveComplementPeaks.csv")
Shared <- read.csv("Chl_1Q_t60_vs_Chl_ndc_t60_SharedComplementPeaks.csv")

#trim the files and combine shared and NDC exclusive files
NDCexc <- select(NDCexc, c(3,4))
Shared <- select(Shared, c(1,2))
colnames(Shared) <- colnames(NDCexc)
NDCexc <- rbind(NDCexc, Shared)

Treatment <- select(Treatment, c(3,4))
colnames(Treatment) <- c("TreatmentHighestPeak", "TreatmentCoverage")

#Work on the NA sites of the shorter file
dim(NDCexc)
dim(Treatment)
rowdif = dim(NDCexc) - dim(Treatment)
fill <- data.frame((matrix(ncol = 2, nrow = rowdif)))
fill[, c(1,2)] <- "0"
colnames(fill) <- colnames(Treatment)
Treatment <- rbind(Treatment, fill)

Treatment[, c(3,4)] <- NDCexc[,c(1,2)]
colnames(Treatment)[3] <- "Shared_NDC_HighestPeak"
colnames(Treatment)[4] <- "Shared_NDC_Coverage"

Treatment$TreatmentHighestPeak <- as.numeric(Treatment$TreatmentHighestPeak)
Treatment$Shared_NDC_HighestPeak <- as.numeric(Treatment$Shared_NDC_HighestPeak)

#Calling the annotation file adjusting the gene UTRs
annot <- read.delim(file = "NC_003028.v3.17.ncrna.genes", header = FALSE)
annot <- annot[, c(2,3,5,6)]
colnames(annot) <- c("Start", "End", "GeneFunction", "GeneName")

for(w in 1:nrow(annot)){
  annot$Start[w] <- annot$Start[w] - 150
  annot$End[w] <- annot$End[w] + 150
}

annot$GeneFunction <- str_replace_all(annot$GeneFunction, ",", "_")

#Calling gene names wrt to their coordinates
Treatment$TreatmentGene <- ""
Treatment$NDCGene <- ""
Treatment$TreatmentGeneFuncion<- ""
Treatment$NDCGeneFunction <- ""


for(i in 1:nrow(Treatment)){
  for(j in 1:nrow(annot)){
    if((Treatment$Shared_NDC_HighestPeak[i] > annot$Start[j]) & (Treatment$Shared_NDC_HighestPeak[i] < annot$End[j])){
      Treatment$NDCGene[i] <- annot$GeneName[j]
      Treatment$NDCGeneFunction[i] <- annot$GeneFunction[j]
    }
  }
}

#choosing the highest cov and > 10% of it.
sample <- Treatment

listOfCats <- unique(sample$NDCGene)

tempFrame <- c()

outputFrame <- as.data.frame(matrix(nrow = 0, ncol = ncol(sample)))
colnames(outputFrame) <- colnames(sample)


for(x in 1:length(listOfCats)){
  
  tempFrame <- sample[sample$NDCGene %in% listOfCats[x],]
  
  tempMaximum <- max(tempFrame$Shared_NDC_Coverage)
  
  tempFrame <- tempFrame[tempFrame$Shared_NDC_Coverage >= (.1*tempMaximum),]
  
  outputFrame <- rbind(tempFrame, outputFrame)
  
}

rownames(outputFrame) <- 1:nrow(outputFrame)


#Calling the closest coordinates from the treatment file wrt to NDCShared coordinates
outputFrame$TreatmentPosCalled <- ""
inds = sapply(outputFrame$Shared_NDC_HighestPeak, function(x) which.min(abs(x - outputFrame$TreatmentHighestPeak)))
result <- transform(outputFrame, TreatmentPosCalled = outputFrame$TreatmentHighestPeak[inds])
result <- result[, -c(1,2)]
result$TreatmentCov <- ""
Treatment <- na.omit(Treatment)
for(x in 1:nrow(result)){
  for(y in 1:nrow(Treatment)){
    if(result$TreatmentPosCalled[x] == Treatment$TreatmentHighestPeak[y]){
      result$TreatmentCov[x] <- Treatment$TreatmentCoverage[y]
    }
  }
}

#Adjusting the output results and calling the gene annotations again now with the trimmed coordinates
output <- result[,c(1,2,7,8,4,6)]
output$TreatmentGene <- ""
output$TreatmentGeneFunction <- ""

for(i in 1:nrow(output)){
  for(j in 1:nrow(annot)){
    if((output$TreatmentPosCalled[i] > annot$Start[j]) & (output$TreatmentPosCalled[i] < annot$End[j])){
      output$TreatmentGene[i] <- annot$GeneName[j]
      output$TreatmentGeneFuncion[i] <- annot$GeneFunction[j]
    }
  }
}

output <- output[, -c(8)]
output$Distance <- ""
output$Distance <- output$TreatmentPosCalled - output$Shared_NDC_HighestPeak
output <- output[, c(1,2,3,4,9,5,6,7,8)]

output$AllEqual <- ""
output$AllEqual <- output$TreatmentGene==output$NDCGene

output <- arrange(output, NDCGene)

#Change the output file name
write.csv(output, file = "Chl_1Q_vs_NDCsharedSupplemented_T60_Compement_Newannot_TrimmedEarly.csv", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)



