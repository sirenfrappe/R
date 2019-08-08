#DNCON2预测结果rr文件转化为邻接矩阵
library(tidyverse)
path <- dir("D:\\idps\\idps\\output")
for (i in path) {
  filename <- str_c("D:\\idps\\idps\\output\\",i,"\\",i,".rr.raw")
  data <- read.table(filename)[,c(1,2,5)]
  len <- data[dim(data)[1],2]
  outmatrix <- matrix(1,nrow = len,ncol = len)
  for (j in 1:dim(data)[1]) {
    num1 <- data[j,1]
    num2 <- data[j,2]
    p <- data[j,3]
    outmatrix[num1,num2] <- p
    outmatrix[num2,num1] <- p
  }
  write.table(outmatrix,file = str_c("D:\\idps\\script\\output\\rr2matrix\\",i,".txt"),sep="\t",quote = FALSE,
              row.names = FALSE,col.names = F)
}
