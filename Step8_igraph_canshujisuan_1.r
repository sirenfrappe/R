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
data_PTM  <-   filter(data,Modification != 0)
data_PTM <-  filter(data_PTM,is.na(mobidb.disorder.predictors.mobidb.lite.score)==F)
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

#表格整理，含网络参数，在每类总的均值、中位数
#将data_PTM最后两列整理，方便summarize()使用
#倒数第二列0表示疾病无关，1表示疾病有关
#最后一列0表示无序，1表示有序
for (i in 1:2933) {
  if(data_PTM[i,10]!=0){
    data_PTM[i,10] <- 1
  }
  if(data_PTM[i,11]>0.501){
    data_PTM[i,11] <- 0
  }else{
    data_PTM[i,11] <- 1
  }
}
#分组，均值
a <- group_by(data_PTM,RelatedDisease,mobidb.disorder.predictors.mobidb.lite.score) %>%
  summarise(mean(degree),mean(closeness),mean(eigen_centrality),mean(betweenness),mean(page_rank))
#表格第一列为参数
parameter <- tibble(parameter=c("degree","closeness","eigen_centrality","betweenness","page_rank"))
################################################################
###############################################################
D0_D1 <- mutate(parameter,
                mean_D0 = t(a[3,3:7])[,1],
                mean_D1 = t(a[4,3:7])[,1],
                p = c(wilcox.test(D0$degree,D1$degree)$p.value,
                      wilcox.test(D0$closeness,D1$closeness)$p.value,
                      wilcox.test(D0$eigen_centrality,D1$eigen_centrality)$p.value,
                      wilcox.test(D0$betweenness,D1$betweenness)$p.value,
                      wilcox.test(D0$page_rank,D1$page_rank)$p.value)
                )
write_csv(D0_D1,"./output/igraph_canshujisuan_1/D0_D1.csv")
###############################################################
D0_N0 <- mutate(parameter,
                mean_D0 = t(a[3,3:7])[,1],
                mean_N0 = t(a[1,3:7])[,1],
                p = c(wilcox.test(D0$degree,N0$degree)$p.value,
                      wilcox.test(D0$closeness,N0$closeness)$p.value,
                      wilcox.test(D0$eigen_centrality,N0$eigen_centrality)$p.value,
                      wilcox.test(D0$betweenness,N0$betweenness)$p.value,
                      wilcox.test(D0$page_rank,N0$page_rank)$p.value)
)
write_csv(D0_N0,"./output/igraph_canshujisuan_1/D0_N0.csv")
###############################################################
D0_N1 <- mutate(parameter,
                mean_D0 = t(a[3,3:7])[,1],
                mean_N1 = t(a[2,3:7])[,1],
                p = c(wilcox.test(D0$degree,N1$degree)$p.value,
                      wilcox.test(D0$closeness,N1$closeness)$p.value,
                      wilcox.test(D0$eigen_centrality,N1$eigen_centrality)$p.value,
                      wilcox.test(D0$betweenness,N1$betweenness)$p.value,
                      wilcox.test(D0$page_rank,N1$page_rank)$p.value)
)
write_csv(D0_N1,"./output/igraph_canshujisuan_1/D0_N1.csv")
###############################################################
D1_N0 <- mutate(parameter,
                mean_D1 = t(a[4,3:7])[,1],
                mean_N0 = t(a[1,3:7])[,1],
                p = c(wilcox.test(D1$degree,N0$degree)$p.value,
                      wilcox.test(D1$closeness,N0$closeness)$p.value,
                      wilcox.test(D1$eigen_centrality,N0$eigen_centrality)$p.value,
                      wilcox.test(D1$betweenness,N0$betweenness)$p.value,
                      wilcox.test(D1$page_rank,N0$page_rank)$p.value)
)
write_csv(D1_N0,"./output/igraph_canshujisuan_1/D1_N0.csv")
###############################################################
D1_N1 <- mutate(parameter,
                mean_D1 = t(a[4,3:7])[,1],
                mean_N1 = t(a[2,3:7])[,1],
                p = c(wilcox.test(D1$degree,N1$degree)$p.value,
                      wilcox.test(D1$closeness,N1$closeness)$p.value,
                      wilcox.test(D1$eigen_centrality,N1$eigen_centrality)$p.value,
                      wilcox.test(D1$betweenness,N1$betweenness)$p.value,
                      wilcox.test(D1$page_rank,N1$page_rank)$p.value)
)
write_csv(D1_N1,"./output/igraph_canshujisuan_1/D1_N1.csv")
###############################################################
N0_N1 <- mutate(parameter,
                mean_N0 = t(a[1,3:7])[,1],
                mean_N1 = t(a[2,3:7])[,1],
                p = c(wilcox.test(N0$degree,N1$degree)$p.value,
                      wilcox.test(N0$closeness,N1$closeness)$p.value,
                      wilcox.test(N0$eigen_centrality,N1$eigen_centrality)$p.value,
                      wilcox.test(N0$betweenness,N1$betweenness)$p.value,
                      wilcox.test(N0$page_rank,N1$page_rank)$p.value)
)
write_csv(N0_N1,"./output/igraph_canshujisuan_1/N0_N1.csv")
#################################################################
data_PTM_draw <- select(data_PTM,-c(position, residue,transitivity,Modification))
#加入新列class用于画图
data_PTM_draw <- mutate(data_PTM_draw,classification=NA)
for (i in 1:2933) {
  if(data_PTM_draw[i,6]== 0 && data_PTM_draw[i,7] == 0){
    data_PTM_draw[i,8] <- "Normal&Disorder"
  }
  if(data_PTM_draw[i,6]== 0 && data_PTM_draw[i,7] == 1){
    data_PTM_draw[i,8] <- "Normal&Structure"
  }
  if(data_PTM_draw[i,6]== 1 && data_PTM_draw[i,7] == 0){
    data_PTM_draw[i,8] <- "Disease&Disorder"
  }
  if(data_PTM_draw[i,6]== 1 && data_PTM_draw[i,7] == 1){
    data_PTM_draw[i,8] <- "Disease&Structure"
  }
}
###########################################################
#boxplot
library(ggsignif)
###########################################################
##drgree
ggplot(data_PTM_draw,aes(x=classification,y=degree))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("Normal&Disorder","Normal&Structure"),
                                 c("Normal&Disorder","Disease&Disorder"),
                                 c("Normal&Disorder","Disease&Structure"),
                                 c("Normal&Structure","Disease&Disorder"),
                                 c("Normal&Structure","Disease&Structure"),
                                 c("Disease&Disorder","Disease&Structure")
                                 ),
              test = wilcox.test,
              step_increase = 0.1
              )+
  ggsave("./output/igraph_canshujisuan_1/degree.tiff")
##closeness
ggplot(data_PTM_draw,aes(x=classification,y=closeness))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("Normal&Disorder","Normal&Structure"),
                                 c("Normal&Disorder","Disease&Disorder"),
                                 c("Normal&Disorder","Disease&Structure"),
                                 c("Normal&Structure","Disease&Disorder"),
                                 c("Normal&Structure","Disease&Structure"),
                                 c("Disease&Disorder","Disease&Structure")
  ),
  test = wilcox.test,
  step_increase = 0.1
  )+
  ggsave("./output/igraph_canshujisuan_1/closeness.tiff")
##eigen_centrality
ggplot(data_PTM_draw,aes(x=classification,y=eigen_centrality))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("Normal&Disorder","Normal&Structure"),
                                 c("Normal&Disorder","Disease&Disorder"),
                                 c("Normal&Disorder","Disease&Structure"),
                                 c("Normal&Structure","Disease&Disorder"),
                                 c("Normal&Structure","Disease&Structure"),
                                 c("Disease&Disorder","Disease&Structure")
  ),
  test = wilcox.test,
  step_increase = 0.1
  )+
  ggsave("./output/igraph_canshujisuan_1/eigen_centrality.tiff")
##betweenness
ggplot(data_PTM_draw,aes(x=classification,y=betweenness))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("Normal&Disorder","Normal&Structure"),
                                 c("Normal&Disorder","Disease&Disorder"),
                                 c("Normal&Disorder","Disease&Structure"),
                                 c("Normal&Structure","Disease&Disorder"),
                                 c("Normal&Structure","Disease&Structure"),
                                 c("Disease&Disorder","Disease&Structure")
  ),
  test = wilcox.test,
  step_increase = 0.1
  )+
  ggsave("./output/igraph_canshujisuan_1/betweenness.tiff")
##page_rank
ggplot(data_PTM_draw,aes(x=classification,y=page_rank))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("Normal&Disorder","Normal&Structure"),
                                 c("Normal&Disorder","Disease&Disorder"),
                                 c("Normal&Disorder","Disease&Structure"),
                                 c("Normal&Structure","Disease&Disorder"),
                                 c("Normal&Structure","Disease&Structure"),
                                 c("Disease&Disorder","Disease&Structure")
  ),
  test = wilcox.test,
  step_increase = 0.1
  )+
  ggsave("./output/igraph_canshujisuan_1/page_rank.tiff")
