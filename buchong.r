library(tidyverse)
library(ggpmisc)
#读取文件
data <- tibble()
for (i in dir("D:/idps/igraph_shortest/shortest_path2")) {
  df <- read_tsv(str_c("D:/idps/igraph_shortest/shortest_path2/",i)) %>%
    select(-(1:2))
  len <- dim(df)[1]
  mean_path <- sum(df)/(len*(len-1))
  data <- bind_rows(data,tibble(meanPath=mean_path,len=len))
}
ggplot(data,aes(x=len,y=meanPath))+
  geom_point()+
  geom_smooth(se=FALSE)+
  xlab("网络节点数")+
  ylab("平均最短路径")+
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')), formula = y ~ x,parse = T)
