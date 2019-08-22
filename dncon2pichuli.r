#查询未完成任务ID，假设已完成预测结果全部存放在path变量所指的文件夹下
#根据未完成任务，设置脚本
library(tidyverse)
path <- "D:\\idps\\idps\\output"
filename <- dir(path)

id <-
  read.table(
    "D:/idps/data/dbPTM/uniprot-yourlist_M20190417216DA2B77BFBD2E6699CA9B6D1C41EB202F439B.tab",
    header = T,
    sep = "\t",
    stringsAsFactors = FALSE
  )[,c(1,8)]
colnames(id) <- c("ID", "length")
id <- as_tibble(id)
doneName <- tibble(ID = filename)



undone <- id %>%
  anti_join(doneName) %>% 
  arrange(length)

for (i in 1:dim(undone)[1]) {
  str1 <- str_c("~/DNCON2/DNCON2/dncon2-v1.0.sh ~/DNCON2/protein/",undone[i,1],".fasta ~/DNCON2/output/",undone[i,1],"/ > ~/DNCON2/log/",undone[i,1],".log")
  cat(str1,"wait",append = T,sep = "\n",file = "D:\\idps\\script\\dncon2RUN.sh")
}
