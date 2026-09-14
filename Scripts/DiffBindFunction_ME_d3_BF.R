#cd Z:\NAS_Server_Shared\08-Cxt_RVFV\bam\Individual Bams\with duplicates
setwd("Z:/NAS_Server_Shared/08-Cxt_RVFV/bam/Round1_Round2_Bams/02.4_Cxt_BamMerged_with_replicate_count")
library(BiocParallel)
register(SerialParam())
library(DiffBind)
library(tidyverse)
library(edgeR)
library(GenomicRanges)

stage_bams_locally <- function(sampleFile, local_dir = "C:/Cxt_temp/bam") {
  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
  samplesOnly <- read.csv(sampleFile)
  bam_files <- unique(samplesOnly$bamReads)
  bai_files <- paste0(bam_files, ".bai")
  for (f in c(bam_files, bai_files)) {
    dest <- file.path(local_dir, basename(f))
    if (!file.exists(dest)) {
      message("Copying ", basename(f), " ...")
      file.copy(f, dest, overwrite = FALSE)
    } else {
      message(basename(f), " already staged, skipping.")
    }
  }
  samplesOnly$bamReads <- file.path(local_dir, basename(samplesOnly$bamReads))
  new_csv <- sub("\\.csv$", "_local.csv", sampleFile)
  write.csv(samplesOnly, new_csv, row.names = FALSE)
  message("Wrote updated sample sheet: ", new_csv)
  return(new_csv)
}

runDiffbind = function(sampleFile,samples1,samples2,plotting=FALSE,saving=TRUE,summits=FALSE,sampleName=FALSE,
                       namesOne="RVFV_d1Me",namesTwo="BF",fromFile=FALSE,CSV=TRUE,edger=TRUE){
  if(fromFile){
    contrastOnly=readRDS(sampleFile)
  } else{
    samplesOnly <- read.csv(sampleFile)
    dbOnly <- dba(sampleSheet=samplesOnly)
    print("Past dba")
    resultOnly <- dba.count(dbOnly,summits=summits)
    if(saving){
      saveRDS(resultOnly,file=paste(namesOne,"_vs_",namesTwo,"PeakCounts.RDS",sep = ""))
    }
    testPeaks=dba.peakset(resultOnly)
    normOnly <- dba.normalize(resultOnly)
    contrastOnly <- dba.contrast(normOnly,group1=samples1, group2=samples2,
                                 name1=namesOne, name2=namesTwo,minMembers = 2)
  }
  print("Here")
  analizOnly <- dba.analyze(contrastOnly,method=DBA_EDGER,bBlacklist = FALSE,bGreylist = FALSE)
  
  if(edger){
    res_deseq <- dba.report(analizOnly, method=DBA_EDGER, contrast = 1, th=1)
  } else{
    res_deseq <- dba.report(analizOnly, method=DBA_DESEQ2, contrast = 1, th=1)
  }
  
  resDF <- as.data.frame(res_deseq) %>% rename(chr = seqnames)
  
  if(saving){
    if(CSV){
      outname = paste(namesOne,"_vs_",namesTwo,".csv",sep="")
      print(paste("saving csv file:", outname))
      write.csv(resDF,outname,row.names=FALSE)
    }
    saveRDS(contrastOnly,paste(namesOne,"_vs_",namesTwo,".RDS",sep=""))
  }
  
  if(plotting){
    par(mar=c(4,4,4,4))
    print("Generating Global PCA...")
    dba.plotPCA(analizOnly, attributes=DBA_FACTOR, label=DBA_ID, mask=analizOnly$masks$All)
    print("Generating Contrast PCA...")
    dba.plotPCA(analizOnly, attributes=DBA_FACTOR, label=DBA_ID, contrast=1, th=1, method=DBA_EDGER)
    plot(analizOnly)
    dba.plotVolcano(analizOnly, method=DBA_EDGER)
    dba.plotMA(analizOnly, bXY=TRUE, method=DBA_EDGER)
  }
  return(resDF)
}

setwd("Z:/NAS_Server_Shared/08-Cxt_RVFV/diffbind_d3")

# confirm the exact ME csv filename before running
list.files(pattern = "\\.csv$")


local_csv_me <- stage_bams_locally("MPvBF_d3_ME123-withdupes_in.csv")  # use actual filename

MPBFd3_ME = runDiffbind(local_csv_me,
                        samples1 = c(1:3),
                        samples2 = c(4:6),
                        namesOne = "MP_ME_d3",
                        namesTwo = "BF-OUT",
                        plotting = TRUE)
