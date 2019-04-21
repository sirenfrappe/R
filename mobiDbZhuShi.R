#mobiDb注释
#使用1A01_HUMAN为例，对应AC为P30443
#设置工作目录
setwd("D:/idps/srcipt")
#加载包，用于获取请求网页的得到的JSON文件
library("httr")
data <-
  content(GET("http://mobidb.bio.unipd.it/ws/Q13569/consensus"))
#sequ为mobidb数据库中对应的序列
sequ <-
  content(GET("http://mobidb.bio.unipd.it/ws/P01912/uniprot"))$sequence
#读取数据
table <-
  read.table(
    "D:/idps/script/output/UniID/2B13_HUMAN.txt",
    sep = "\t",
    header = T,
    stringsAsFactors = F
  )
sequDbptm <- paste(as.character(table$residue), collapse = "")
#判断mobidb数据库与dbPTM数据库中序列是否相同
if (sequ != sequDbptm) {
  cat("两序列不相同！\n")
  stop()
}
#来自手动添加的mobidb数据
#这里错误，db下的内容应该是data$mobidb_consensus$disorder$db，有的蛋白数据中不包含此项目
disorder.db <- data$mobidb_consensus$disorder$full[[1]][[1]]
tem <- data.frame(matrix(NA, dim(table)[1], 1))
for (i in 1:length(disorder.db)) {
  tem <- regions(disorder.db, i, tem)
}
colnames(tem) <- "disorder.db"
table <- cbind(table, tem)
#来自mobidb中derived数据
disorder.derived <- data$mobidb_consensus$disorder$derived
disorder.derived.length <-  length(disorder.derived)
for (i in 1:disorder.derived.length) {
  if (disorder.derived[[i]]$method == "full") {
    tem <- data.frame(matrix(NA, dim(table)[1], 1))
    for (j in 1:length(disorder.derived[[i]][[1]])) {
      tem <- regions(disorder.db = disorder.derived[[i]]$regions, j, tem)
    }
    colnames(tem) <- "disorder.derived.full"
    table <- cbind(table, tem)
  } else{
    tem  <- data.frame(matrix(NA, dim(table)[1], 2))
    method <- disorder.derived[[i]]$method
    colnames(tem) <-
      c(
        paste("disorder.derived", "regions", sep = "."),
        paste("disorder.derived", method, "scores", sep = ".")
      )
    #打分数据并入
    tem[, 2] <- as.character(disorder.derived[[1]]$scores_n)
    #打分结果并入
    for (j in 1:length(disorder.derived[[i]][[1]])) {
      tem <- regions(disorder.db = disorder.derived[[i]]$regions, j, tem)
    }
    table <- cbind(table,tem)
  }
}


#无scores的数据单元添加进上一步处理得到表格时用到的函数
#method = full\simple
regions <- function(disorder.db, i, tem) {
  qidian <- disorder.db[[i]][[1]]
  zhongdian  <- disorder.db[[i]][[2]]
  stru <- disorder.db[[i]][[3]]
  tem[qidian:zhongdian, 1] <- stru
  return(tem)
}
