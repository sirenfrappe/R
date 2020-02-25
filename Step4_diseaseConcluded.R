# 疾病相关数据总结

#设置工作目录
setwd("D:/idps/script")

#新建dataframe,存储疾病名称，蛋白ID和该疾病位点对应PTM类型
data <-
  data.frame(
    "diease" = character(),
    "proID" = character(),
    "PTMtype" = character(),
    stringsAsFactors = FALSE
  )

path <- "D:/idps/script/output/UniID"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x) {
  paste(path, x, sep = '/')
})

#mulu <- "D:/idps/script/output/UniID/1B07_HUMAN.txt"


for (m in 1:length(filePath)) {
 #读取文件
  a <- read.table(filePath[m],
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  header = T,
                  quote = NULL)
  #筛选出有疾病相关的行
  a <- subset(a, a$RelatedDisease != 0)
  #循环处理
  for (i in 1:dim(a)[1]) {
#无逗号，strsplit结果的长度为1直接添加即可
    if (length(strsplit(a[i, ]$RelatedDisease,",")[[1]]) == 1) {
      data[dim(data)[1]+ 1, 1] <- a[i, ]$RelatedDisease
      data[dim(data)[1], 2] <-
        strsplit(strsplit(filePath[m], "/")[[1]][6], ".txt")[[1]]
      data[dim(data)[1], 3] <- a[i, 3]
    } else{
      num <- length(strsplit(a[i, ]$RelatedDisease,",")[[1]])
      for (j in 1:num) {
        data[dim(data)[1] + 1, 1] <- strsplit(a[i, ]$RelatedDisease, ",")[[1]][j]
        data[dim(data)[1], 2] <-
          strsplit(strsplit(filePath[m], "/")[[1]][6], ".txt")[[1]]
        data[dim(data)[1], 3] <- a[i, 3]
      }
    }
  }
}
data$diease <- as.factor(data$diease)
data$proID <- as.factor(data$proID)
data$PTMtype<- as.factor(data$PTMtype)
