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

################################
#6.9
################################

#蛋白长度与PTM个数关系图
data.lengthANDPTM <- data.frame(matrix(NA,143,2))
colnames(data.lengthANDPTM) <- c("length","num")
for (i in 1:143) {
  data.lengthANDPTM[i,1] <- dim(data[[i]])[1]
  data.lengthANDPTM[i,2] <- dim(subset(data[[i]],data[[i]]$Modification!=0))[[1]]
}
library(ggplot2)
ggplot(data.lengthANDPTM,aes(x=data.lengthANDPTM$length,y=data.lengthANDPTM$num))+
  geom_point()+
  geom_text(label=data.lengthANDPTM$num,vjust=-0.5)+
  xlab("长度")+
  ylab("个数")+
  geom_smooth(method = 'auto')

#去除蛋白长度对位点个数的影响
data.rate.normaldisease <- data.frame(matrix(NA,143,2))
colnames(data.rate.normaldisease) <- c("normal","disease")
#分母：PTM个数*蛋白长度
for (i in 1:143) {
  data.rate.normaldisease[i,2] <- dim(data.disease[[i]])[1]/(dim(data[[i]])[1]*
                                                              dim(
                                                                subset(
                                                                  data[[i]],
                                                                  data[[i]]$Modification!=0)
                                                              )[1]
  )
  data.rate.normaldisease[i,1] <- (dim(subset(data[[i]],data[[i]]$Modification!=0))[1]-dim(data.disease[[i]])[1])/(dim(data[[i]])[1]*
                                                                                                                    dim(
                                                                                                                      subset(
                                                                                                                        data[[i]],
                                                                                                                        data[[i]]$Modification!=0)
                                                                                                                    )[1]
  )
}

#差异性检验
shapiro.test(data.rate.normaldisease[,1])
#p-value = 9.919e-12
shapiro.test(data.rate.normaldisease[,2])
#p-value < 2.2e-16
#方差齐性检验
bartlett.test(c(data.rate.normaldisease[,1],data.rate.normaldisease[,2])~factor(c(rep(1,143),rep(2,143))))
#p=1.817e-09

#远小于0.05，采用非参数检验
wilcox.test(data.rate.normaldisease[,1],data.rate.normaldisease[,2],paired=TRUE)
#p-value < 2.2e-16

data.PTM <- lapply(data, function(x){
  subset(x,x$Modification!=0)
})
#疾病、无序/结构，打分>0.501判断
data.disease.stru <- lapply(data.disease, function(x){
  subset(x,x$mobidb.disorder.predictors.mobidb.lite.score<0.501)
})
data.disease.diso <- lapply(data.disease, function(x){
  subset(x,x$mobidb.disorder.predictors.mobidb.lite.score>0.501)
})
#疾病无关、无序/结构，打分>0.501判断
data.normal <- lapply(data, function(x){
  subset(x,x$RelatedDisease==0 & x$Modification!=0 )
})
data.normal.stru <- lapply(data.normal, function(x){
  subset(x,x$mobidb.disorder.predictors.mobidb.lite.score<0.501)
})
data.normal.diso <- lapply(data.normal, function(x){
  subset(x,x$mobidb.disorder.predictors.mobidb.lite.score>0.501)
})
#疾病不相关PTM位点无序与否差异性检测
data.rate.NPTM <- data.frame(matrix(NA,143,2))
colnames(data.rate.NPTM) <- c("stru","diso")
#疾病无关
#分母：PTM个数*蛋白长度
for (i in 1:143) {
  data.rate.NPTM[i,1] <- dim(data.normal.stru[[i]])[1]/(dim(data[[i]])[1]*
                                                          dim(
                                                            subset(
                                                              data[[i]],
                                                              data[[i]]$Modification!=0)
                                                          )[1])
  data.rate.NPTM[i,2] <- dim(data.normal.diso[[i]])[1]/(dim(data[[i]])[1]*
                                                          dim(
                                                            subset(
                                                              data[[i]],
                                                              data[[i]]$Modification!=0)
                                                          )[1])
}
shapiro.test(data.rate.NPTM[,1])
shapiro.test(data.rate.NPTM[,2])
#p<0.05
bartlett.test(c(data.rate.NPTM[,1],data.rate.NPTM[,2])~factor(c(rep(1,143),rep(2,143))))
#p-value = 1.043e-08
wilcox.test(data.rate.NPTM[,1],data.rate.NPTM[,2],paired=TRUE)
#p-value = 6.57e-10，差异性显著

#疾病相关
data.rate.DPTM <- data.frame(matrix(NA,143,2))
colnames(data.rate.DPTM) <- c("stru","diso")
for (i in 1:143) {
  data.rate.DPTM[i,1] <- dim(data.disease.stru[[i]])[1]/(dim(data[[i]])[1]*
                                                          dim(
                                                            subset(
                                                              data[[i]],
                                                              data[[i]]$Modification!=0)
                                                          )[1])
  data.rate.DPTM[i,2] <- dim(data.disease.diso[[i]])[1]/(dim(data[[i]])[1]*
                                                          dim(
                                                            subset(
                                                              data[[i]],
                                                              data[[i]]$Modification!=0)
                                                          )[1])
}
shapiro.test(data.rate.DPTM[,1])
shapiro.test(data.rate.DPTM[,2])
#p<0.05
bartlett.test(c(data.rate.DPTM[,1],data.rate.DPTM[,2])~factor(c(rep(1,143),rep(2,143))))
# p-value < 2.2e-16
wilcox.test(data.rate.DPTM[,1],data.rate.DPTM[,2],paired=TRUE)
#p-value = 5.421e-07，差异性显著

#疾病无序，非疾病无序比较
shapiro.test(data.rate.DPTM[,2])
shapiro.test(data.rate.NPTM[,2])
#p<0.05
bartlett.test(c(data.rate.DPTM[,2],data.rate.NPTM[,2])~factor(c(rep(1,143),rep(2,143))))
#p-value = 0.01032
wilcox.test(data.rate.DPTM[,2],data.rate.NPTM[,2],paired=TRUE)
