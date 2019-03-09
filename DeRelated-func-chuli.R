


#包加载
library(rvest)
library(xml2)
library(magrittr)

ID <- NULL

for (i in 1:15) {
  URL <- paste("D:/idps/data/dbPTM/muluWeb/dbPTM", i, ".html", sep = "")
  page <- read_html(URL)
  table <-
    data.frame(html_table(
      page %>% html_nodes(".table-responsive") %>% html_nodes("table"),
      header = TRUE
    ))
  ID <- c(ID, unique(table[, 1]))
}
ID <- unique(ID)

#ID为含有疾病相关PTM的蛋白的uniprotID
#其中部分为dbPTM 中没有信息的ID
#下面进行无信息ID的筛除

temp <- NULL
for (i in 88:length(ID)) {
  cat(i, "\t", "开始", "\t")
  URL <-
    paste("http://dbptm.mbc.nctu.edu.tw/info.php?id=", ID[i], sep = "")
  
  page <- read_html(URL)
  cat("url wancheng\t")
  Nnodes <- page %>% html_nodes(".panel-body")
  
  cat("节点完毕\t")
  if (length(Nnodes) == 1) {
    temp <- c(temp, i)
    cat("空数据", "\t")
  }
  cat("jieshu \t")
  cat(ID[i], "\n")
  Sys.sleep(5)
}

ID <- ID[-temp]
#ID为存在的与疾病相关的蛋白uniprotID,共143个
cat(ID, file = "D:/idps/script/output/UniIDofDieaseassociated.txt")

library(stringr)
temp <- NULL
for (i in 1:length(ID)) {
  cat(i, "\t", ID[i], "\t")
  URL <-
    paste("http://dbptm.mbc.nctu.edu.tw/info.php?id=", ID[i], sep = "")
  page <- read_html(URL)
  cat("网页加载完成\t")
  DaPS <-
    as.data.frame(html_table(
      page %>% html_nodes(".table-responsive") %>% extract2(3) %>% html_nodes("table")
    ))[, c(1, 2, 3, 4, 6)]
  cat("数据框载入完成\t")
  temp <- str_count(DaPS[, 3], ";")
  cat("temp:", temp, "\t")
  for (j in 1:length(temp)) {
    if (temp[j] != 0) {
      cat(ID[i], "\t", DaPS[j, 3], "\t", DaPS[j, 4], "\t")
    }
  }
  cat("此条结束", "\n")
  Sys.sleep(2)
}
#ABI3_HUMAN,IQGA2_HUMAN,FA5_HUMAN,MSHR_HUMAN,TLR10_HUMAN,VTDB_HUMAN
#
#由于文件已经保存，可以直接从此处运行
ID <-
  readLines("D:/idps/script/output/UniIDofDieaseassociated.txt")
ID <- strsplit(ID, " ")[[1]]
#疾病相关蛋白全部处理
IDD <- ID
source("D:/idps/script/zhushi.R")
for (i in 1:length(IDD)) {
  lujing <-
    paste("D:/idps/script/output/UniID/",
          IDD[i],
          ".txt",
          sep = "")
  cat(i, "\t")
  func(IDD[i], lujing)
  cat("完成\n")
}
