path <- "D:/idps/script/output/UniID/"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x) {
  paste(path, x, sep = '/')
})
data <- lapply(filePath, function(x) {
  read.table(x, header = T,sep = "\t",stringsAsFactors = FALSE,quote=NULL)
})
lapply(data,function(x){x[dim(x)[1],]})
#观察是否有NA，共13个
#VTDB_HUMAN、PLMN_HUMAN、MSHR_HUMAN、LCAT_HUMAN、HRG_HUMAN、DRA_HUMAN、CFAH_HUMAN、CBPB2_HUMAN、APOH_HUMAN、APOC3_HUMAN、AFAM_HUMAN、1B07_HUMAN、1A03_HUMAN