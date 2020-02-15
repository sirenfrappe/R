#igraph参数比较

#读文件,将igraph输出文件与mobidb中矩阵PTM列、疾病列与无序打分列合并
library(tidyverse)
data1 <- tibble()
for (i in dir("D:/idps/igraph_shortest/igraph_out1")) {
  df <- read_tsv(str_c("D:/idps/igraph_shortest/igraph_out1/",i),col_names = c("position","residue","transitivity","degree","closeness","eigen_centrality","betweenness","page_rank"),skip = 1)
  data1 <- bind_rows(data1,df)
}
#data1将所有igraph输出结果按照文件名顺序依次排列放入tibble中
#这里读mobidb下的文件tibble会报错，改用dataframe
data2 <- data.frame()
for (i in dir("D:/idps/script/output/mobidb")) {
  df <- read.table(str_c("D:/idps/script/output/mobidb/",i),header = T,sep = "\t",quote = NULL,stringsAsFactors = F)
  data2 <- rbind(data2,df)
}
data2 <- data2[,c(3,6,21)]
data2 <- as_tibble(data2)
#data2只含有3列PTM列、疾病列与无序打分列
#合并
data <- bind_cols(data1,data2)
#筛选PTM
data_PTM <- filter(data,Modification != 0)
#共2965行，将分为四类：D0、D1、N0、N1
#D表示疾病相关，N表示疾病无关
#0表示无序，1表示有序
D <- filter(data_PTM,RelatedDisease !=0)
D0 <- filter(D,mobidb.disorder.predictors.mobidb.lite.score > 0.501)
D1 <- filter(D,mobidb.disorder.predictors.mobidb.lite.score < 0.501)
N <- filter(data_PTM,RelatedDisease == 0)
N0 <- filter(data_PTM,RelatedDisease == 0,mobidb.disorder.predictors.mobidb.lite.score > 0.501)
N1 <- filter(data_PTM,RelatedDisease == 0,mobidb.disorder.predictors.mobidb.lite.score < 0.501)

#采用Wilcox秩和检验比较参数之间是否有差别
#degree
wilcox.test(D$degree,N$degree)
wilcox.test(D0$degree,D1$degree)
wilcox.test(D0$degree,N0$degree)
wilcox.test(D0$degree,N1$degree)
wilcox.test(D1$degree,N0$degree)
wilcox.test(D1$degree,N1$degree)
wilcox.test(N0$degree,N1$degree)

#closeness
wilcox.test(D$closeness,N$closeness)
wilcox.test(D0$closeness,D1$closeness)
wilcox.test(D0$closeness,N0$closeness)
wilcox.test(D0$closeness,N1$closeness)
wilcox.test(D1$closeness,N0$closeness)
wilcox.test(D1$closeness,N1$closeness)
wilcox.test(N0$closeness,N1$closeness)

#eigen_centrality
wilcox.test(D$eigen_centrality,N$eigen_centrality)
wilcox.test(D0$eigen_centrality,D1$eigen_centrality)
wilcox.test(D0$eigen_centrality,N0$eigen_centrality)  #p=0.6217>0.05
wilcox.test(D0$eigen_centrality,N1$eigen_centrality)
wilcox.test(D1$eigen_centrality,N0$eigen_centrality)
wilcox.test(D1$eigen_centrality,N1$eigen_centrality)
wilcox.test(N0$eigen_centrality,N1$eigen_centrality)

#betweenness
wilcox.test(D$betweenness,N$betweenness)
wilcox.test(D0$betweenness,D1$betweenness)
wilcox.test(D0$betweenness,N0$betweenness)
wilcox.test(D0$betweenness,N1$betweenness)
wilcox.test(D1$betweenness,N0$betweenness)
wilcox.test(D1$betweenness,N1$betweenness)
wilcox.test(N0$betweenness,N1$betweenness)

#page_rank
wilcox.test(D$page_rank,N$page_rank)
wilcox.test(D0$page_rank,D1$page_rank)
wilcox.test(D0$page_rank,N0$page_rank)
wilcox.test(D0$page_rank,N1$page_rank)  #p=0.7347>0.05
wilcox.test(D1$page_rank,N0$page_rank)
wilcox.test(D1$page_rank,N1$page_rank)
wilcox.test(N0$page_rank,N1$page_rank)

