#添加mobidb注释

mobidbZhushi <- function(UniID, UniAC, outputfile) {
  #设置工作目录
  setwd("D:/idps/srcipt")
  #加载httr包，用于获取请求网页的得到的JSON文件
  library("httr")
  data <-
    content(GET(
      paste("http://mobidb.bio.unipd.it/ws/", UniAC, "/consensus", sep = "")
    ))
  #sequ为mobidb数据库中对应的蛋白质序列
  sequ <-
    content(GET(
      paste("http://mobidb.bio.unipd.it/ws/", UniAC, "/uniprot", sep = "")
    ))$sequence
  #读取处理DBPTM数据库数据得到的表格，命名为table
  table <-
    read.table(
      paste("./output/UniID/", UniID, ".txt", sep = ""),
      sep = "\t",
      header = T,
      stringsAsFactors = F
    )
  #sequDbptm是dbptm中相应蛋白的序列
  sequDbptm <- paste(as.character(table$residue), collapse = "")
  #判断mobidb数据库与dbPTM数据库中序列是否相同，
  if (sequ != sequDbptm) {
    cat("两序列不相同！\n")
    stop()
  }
  #只做了一个例子，暂时还没发现有序列不同的情况
  #来自手动添加的mobidb数据
  #disorder下只有一个元素
  disorder.db <- data$mobidb_consensus$disorder$full[[1]][[1]]
  #新建tem数据框用于存储regions数据
  tem <- data.frame(matrix(NA, dim(table)[1], 1))
  #循环赋值
  for (i in 1:length(disorder.db)) {
    tem <- regions(disorder.db, i, tem)
  }
  #更改列名，合并。
  colnames(tem) <- "disorder.db"
  table <- cbind(table, tem)
  #用于处理derived和predictors数据的函数，两种情况下的结构相似，但区别仍存在，不新建方程处理
  #disorder.derived为来自mobidb中derived数据
  disorder.derived <- data$mobidb_consensus$disorder$derived
  disorder.derived.length <- length(disorder.derived)
  #方法不止一个，建立循环
  for (i in 1:disorder.derived.length) {
    #判断是否使用full方法
    if (disorder.derived[[i]]$method == "full") {
      #新建临时数据框
      tem <- data.frame(matrix(NA, dim(table)[1], 1))
      #根据regions的段数进行循环添加数据进入table
      for (j in 1:length(disorder.derived[[i]][[1]])) {
        tem <- regions(disorder.db = disorder.derived[[i]]$regions, j, tem)
      }
      #更改行名
      colnames(tem) <- "disorder.derived.full"
      table <- cbind(table, tem)
    } else {
      #方法为missing_residues, bfactor, mobile的数据
      tem <- data.frame(matrix(NA, dim(table)[1], 2))
      #获取方法名，后面命名列名时用到
      method <- disorder.derived[[i]]$method
      #tem的列名，一列regions、一列scores
      colnames(tem) <-
        c(
          paste("disorder.derived",method, "regions", sep = "."),
          paste("disorder.derived",method, "scores", sep = ".")
        )
      #打分数据并入
      tem[, 2] <- as.character(disorder.derived[[i]]$scores_n)
      #打分结果并入
      for (j in 1:length(disorder.derived[[i]][[1]])) {
        tem <- regions(disorder.db = disorder.derived[[i]]$regions, j, tem)
      }
      table <- cbind(table, tem)
    }
  }
  
  #来自predictior数据
  disorder.predictiors <- data$mobidb_consensus$disorder$predictors
  disorder.predictiors.length <- length(disorder.predictiors)
  for (i in 1:disorder.predictiors.length) {
    if (disorder.predictiors[[i]]$method == "simple") {
      tem <- data.frame(matrix(NA, dim(table)[1], 1))
      for (j in 1:length(disorder.predictiors[[i]][[1]])) {
        tem <-
          regions(disorder.db = disorder.predictiors[[i]]$regions, j, tem)
      }
      colnames(tem) <- "disorder.predictiors.simple"
      table <- cbind(table, tem)
    } else {
      tem <- data.frame(matrix(NA, dim(table)[1], 2))
      method <- disorder.predictiors[[i]]$method
      colnames(tem) <- c(
        paste("disorder.predictors", method,"regions", sep = "."),
        paste("disorder.predictors", method, "secores", sep = ".")
      )
      tem[, 2] <- as.character(disorder.predictiors[[i]]$scores)
      for (j in 1:length(disorder.predictiors[[i]][[1]])) {
        tem <-
          regions(disorder.db = disorder.predictiors[[i]]$regions, j, tem)
        
      }
      table <- cbind(table, tem)
      
    }
  }
  #此函数用来处理无评分的方法的数据（只含有regions）
  regions <- function(disorder.db, i, tem) {
    qidian <- disorder.db[[i]][[1]]
    zhongian <- disorder.db[[i]][[2]]
    stru <- disorder.db[[i]][[3]]
    tem[qidian:zhongian, 1] <- stru
    return(tem)
  }
  
  write.table(
    table,
    outputfile,
    quote = F,
    sep = "\t",
    row.names = FALSE
  )
  
}
mobidbZhushi("1A01_HUMAN","P30443","1A01_HUMAN(Mobidb).txt")
