path <- "D:/idps/script/output/mobidb"
fileName <- dir(path)
filePath <- sapply(fileName,function(x){
  paste(path,x,sep = "/")
})
data <- lapply(filePath,function(x){
  read.table(x,header = T,sep = "\t",quote = NULL,stringsAsFactors = F)
})
data.diease <- lapply(data, function(x){
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
  data.rate.normaldisease[i,2] <- dim(data.diease[[i]])[1]/(dim(data[[i]])[1]*
                                                              dim(
                                                                subset(
                                                                  data[[i]],
                                                                  data[[i]]$Modification!=0)
                                                              )[1]
  )
  data.rate.normaldisease[i,1] <- (dim(subset(data[[i]],data[[i]]$Modification!=0))[1]-dim(data.diease[[i]])[1])/(dim(data[[i]])[1]*
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