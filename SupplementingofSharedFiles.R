
share <- read.csv("Chl_1Q_t60_vs_Chl_ndc_t60_SharedComplementPeaks.csv")
fileone<- select(share, c(1,2,5,6))
fileone$allequal <- ""
fileone$allequal <- fileone$FileOnePeak==fileone$FileTwoPeak
fileone <- filter(fileone, allequal == FALSE)
rowfileone <- dim(fileone)
filetwo <- data.frame(matrix(nrow = rowfileone, ncol=10))
filetwo[,c(1,2)] <- select(fileone, c(3,4))
colnames(filetwo) <- colnames(share)
share<- rbind(share, filetwo)
