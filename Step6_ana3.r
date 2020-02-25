#dncon2 预测结果rr文件处理

path <- "D:/idps/script/output/mobidb"
fileName <- dir(path)
filePath <- sapply(fileName,function(x){
  paste(path,x,sep = "/")
})
data <- lapply(filePath,function(x){
  read.table(x,header = T,sep = "\t",quote = NULL,stringsAsFactors = F)
})


rr <- read.table("D:\\idps\\idps\\output\\P30443\\P30443.rr.raw")
colnames(rr) <- c("i","j","d1","d2","p")
head(rr)

leng <- 365

sequ<-  data[[1]][,2]
length(sequ)

dfP <- matrix(1,leng+1,leng+1)
dfP[1,2:(leng+1)] <- sequ
dfP[2:(leng+1),1] <- sequ

df01 <- matrix(1,leng,leng)


for (m in 1:dim(rr)[1]) {
  dfP[rr[m,1]+1,rr[m,2]+1] <- rr[m,5]
  dfP[rr[m,2]+1,rr[m,1]+1] <- rr[m,5]
  if(rr[m,5]>0.5){
    df01[rr[m,1],rr[m,2]] <- 1
    df01[rr[m,2],rr[m,1]] <- 1
  }else{
    df01[rr[m,1],rr[m,2]] <- 0
    df01[rr[m,2],rr[m,1]] <- 0
  }
}


library(pheatmap)

pheatmap(df01,cluster_row = FALSE, cluster_col = FALSE,legend = FALSE,show_rownames=F,show_colnames=F)
pheatmap(dfP,cluster_row = FALSE, cluster_col = FALSE,show_rownames=F,show_colnames=F)
write.table(dfP,"D:\\idps\\script\\output\\rrfile\\P30443_dfP.txt",quote = FALSE,row.names = FALSE,col.names = F)
write.table(df01,"D:\\idps\\script\\output\\rrfile\\P30443_df01.txt",quote = FALSE,row.names = FALSE,col.names = F)

#######################################
###########相邻残基注释为接触##########
#######################################

df01_2 <- df01
diag(df01_2[-leng,-1]) <- 1
diag(df01_2[-1,-leng]) <- 1
write.table(df01_2,"D:\\idps\\script\\output\\rrfile\\P30443_df01_2.txt",quote = FALSE,row.names = FALSE,col.names = F)
