rm(list = ls())
library(dplyr)
library(stringr)
Treatment <- read.csv("Ksg_3Q_t60_exclusiveTopPeaks.csv")
NDCexc <- read.csv("Ksg_ndc_t60_exclusiveTopPeaks.csv")
Shared <- read.csv("Ksg_3Q_t60_vs_Ksg_ndc_t60_SharedTopPeaks.csv")

NDCexc <- select(NDCexc, c(3,4))
Shared <- select(Shared, c(1,2))
colnames(Shared) <- colnames(NDCexc)
NDCexc <- rbind(NDCexc, Shared)

Treatment <- select(Treatment, c(3,4))
colnames(Treatment) <- c("TreatmentHighestPeak", "TreatmentCoverage")

dim(NDCexc)
dim(Treatment)
rowdif = dim(NDCexc) - dim(Treatment)
fill <- data.frame((matrix(ncol = 2, nrow = rowdif)))
fill[, c(1,2)] <- "NA"
colnames(fill) <- colnames(Treatment)
Treatment <- rbind(Treatment, fill)

Treatment[, c(3,4)] <- NDCexc[,c(1,2)]
colnames(Treatment)[3] <- "Shared_NDC_HighestPeak"
colnames(Treatment)[4] <- "Shared_NDC_Coverage"

Treatment$TreatmentHighestPeak <- as.numeric(Treatment$TreatmentHighestPeak)
Treatment$Shared_NDC_HighestPeak <- as.numeric(Treatment$Shared_NDC_HighestPeak)

inds = sapply(Treatment$Shared_NDC_HighestPeak, function(x) which.min(abs(x - Treatment$TreatmentHighestPeak)))
result <- transform(Treatment, TreatmentPosCalled = Treatment$TreatmentHighestPeak[inds])
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
output <- transform(result, Distance_Treatment_NDCShared = result$TreatmentPosCalled - result$Shared_NDC_HighestPeak)

annot <- read.delim(file = "NC_003028.v3.17.ncrna.genes", header = FALSE)
annot <- annot[, c(2,3,5,6)]
colnames(annot) <- c("Start", "End", "GeneFunction", "GeneName")

for(w in 1:nrow(annot)){
  annot$Start[w] <- annot$Start[w] - 150
  annot$End[w] <- annot$End[w] + 150
}

annot$GeneFunction <- str_replace_all(annot$GeneFunction, ",", "_")
output$TreatmentGene <- ""
output$NDCGene <- ""
output$TreatmentGeneFuncion<- ""
output$NDCGeneFunction <- ""

for(i in 1:nrow(output)){
  for(j in 1:nrow(annot)){
    if((output$TreatmentPosCalled[i] > annot$Start[j]) & (output$TreatmentPosCalled< annot$End[j])){
      output$TreatmentGene[i] <- annot$GeneName[j]
      output$TreatmentGeneFuncion[i] <- annot$GeneFunction[j]
    }
  }
}

for(i in 1:nrow(output)){
  for(j in 1:nrow(annot)){
    if((output$Shared_NDC_HighestPeak[i] > annot$Start[j]) & (output$Shared_NDC_HighestPeak[i] < annot$End[j])){
      output$NDCGene[i] <- annot$GeneName[j]
      output$NDCGeneFunction[i] <- annot$GeneFunction[j]
    }
  }
}

output$AllEqual <- ""
output$AllEqual <- output$TreatmentGene==output$NDCGene

output <- arrange(output, Distance_Treatment_NDCShared)
output <- arrange(output, TreatmentGeneFuncion)

#Change the output file name
write.csv(output, file = "Ksg_3Q_vs_NDCsharedSupplemented_T60_Top_Newannot.csv", quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

#choosing the highest cov
rm(list = ls())
sample <- read.csv("Ksg_3Q_vs_NDCsharedSupplemented_T60_Top_Newannot.csv")
sample2 <- filter(sample, AllEqual == FALSE)
sample2 <- arrange(sample2, Distance_Treatment_NDCShared)
sample <- filter(sample, AllEqual == TRUE)
sample <- arrange(sample, TreatmentGeneFuncion)


listOfCats <- unique(sample$TreatmentGeneFuncion)

tempFrame <- c()

outputFrame <- as.data.frame(matrix(nrow = 0, ncol = ncol(sample)))
colnames(outputFrame) <- colnames(sample)


for(x in 1:length(listOfCats)){
  
  tempFrame <- sample[sample$TreatmentGeneFuncion %in% listOfCats[x],]
  
  tempMaximum <- max(tempFrame$Shared_NDC_Coverage)
  
  tempFrame <- tempFrame[tempFrame$Shared_NDC_Coverage >= (.1*tempMaximum),]
  
  outputFrame <- rbind(tempFrame, outputFrame)
  
}

rownames(outputFrame) <- 1:nrow(outputFrame)
outputFrame <- arrange(outputFrame, Distance_Treatment_NDCShared)
outputFrame <- rbind(outputFrame, sample2)

write.csv(outputFrame, file = "Ksg_3Q_vs_NDCsharedSupplemented_T60_Top_Newannot_Trimmed.csv", quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)


