path <- "D:/idps/script/output/mobidb"
fileName <- dir(path)
filePath <- sapply(fileName,function(x){
  paste(path,x,sep = "/")
})
data <- lapply(filePath,function(x){
  read.table(x,header = T,sep = "\t",quote = NULL,stringsAsFactors = F)
})
data.disease <- lapply(data, function(x){
  subset(x,x$RelatedDisease!=0)
})

data.disorder <- lapply(data, function(x){
  subset(x,x$mobidb.disorder.predictors.mobidb.lite.score>0.501)
})
############疾病相关的PTM位点在所有PTM位点中的占比
a <- lapply(data.disorder, function(x){
  dim(x)[1]
})
sum <- 0
for (i in 1:143) {
  sum <- sum+a[[i]]
}
sum
#23775
##########所有PTM位点中无序的占比
data.PTM <- data.diease <- lapply(data, function(x){
  subset(x,x$Modification!=0)
})

data.PTM.disorder <- lapply(data.PTM, function(x){
  subset(x,x$mobidb.disorder.predictors.mobidb.lite.score>0.501)
})
b <- lapply(data.PTM.disorder, function(x){
  dim(x)[1]
})
sum <- 0
for (i in 1:143) {
  sum <- sum+b[[i]]
}
sum
###############疾病PTM位点中无序的占比
data.disease <- lapply(data, function(x){
  subset(x,x$RelatedDisease!=0)
})
data.disease.disorder <- lapply(data.disease, function(x){
  subset(x,x$mobidb.disorder.predictors.mobidb.lite.score>0.501)
})
c <- lapply(data.disease.disorder, function(x){
  dim(x)[1]
})
sum <- 0
for (i in 1:143) {
  sum <- sum+c[[i]]
}
sum
