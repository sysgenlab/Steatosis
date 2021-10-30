library(DiffBind)

#Load data
NAFLD <- dba(sampleSheet="DiffBind/sheet/samplesheet.DiffBind.FRiP20.csv")
NAFLD$config$AnalysisMethod <- factor("edgeR")
NAFLD$minOverlap <- 1

#Counting reads
NAFLD <- dba.count(NAFLD)

#Differential analysis
NAFLD <- dba.contrast(NAFLD, categories=DBA_CONDITION)
NAFLD$contrasts[[1]] <- list(group1=NAFLD$masks$fNASH, group2=NAFLD$masks$Steatosis, name1="fNASH", name2="Steatosis")
NAFLD$contrasts[[2]] <- list(group1=NAFLD$masks$Steatosis, group2=NAFLD$masks$Normal, name1="Steatosis", name2="Normal")
NAFLD$contrasts[[3]] <- list(group1=NAFLD$masks$fNASH, group2=NAFLD$masks$Normal, name1="fNASH", name2="Normal")

NAFLD <- dba.analyze(NAFLD)

#Save results
contr <- c("FvsS", "SvsN", "FvsN")

for(i in seq(contr)){
  myData <- dba.report(NAFLD, contrast=i)

  up <- myData[myData$Fold > 0,]
  write.table(up, paste0("DiffBind/Result.DiffBind.", contr[i], ".upregulated.txt"), quote=FALSE, sep="\t", row.names=FALSE)
  df_up <- data.frame(chr=seqnames(up), start=start(up)-1, end=end(up))  
  df_up$name <- paste0(df_up$chr, ":", df_up$start, "-", df_up$end)
  write.table(df_up, paste0("DiffBind/Result.DiffBind.", contr[i], ".upregulated.bed"), quote=FALSE, sep="\t", row.names=FALSE, col.names=F)

  down <- myData[myData$Fold < 0,]
  write.table(down, paste0("DiffBind/Result.DiffBind.", contr[i], ".downregulated.txt"), quote=FALSE, sep="\t", row.names=FALSE)
  df_down <- data.frame(chr=seqnames(down), start=start(down)-1, end=end(down))  
  df_down$name <- paste0(df_down$chr, ":", df_down$start, "-", df_down$end)
  write.table(df_down, paste0("DiffBind/Result.DiffBind.", contr[i], ".downregulated.bed"), quote=FALSE, sep="\t", row.names=FALSE, col.names=F)
}
