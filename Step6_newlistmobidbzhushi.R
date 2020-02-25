#R 3.6.0 encoding=utf-8
#给定ID序列添加注释

#添加mobdbi注释
mobidbzhushi <- function(UniAC, outputfile) {
  #region处理函数,添加regions数据
  regions <- function(region, i, tem) {
    #确定起点终点以及是否无序
    qidian <- region[[i]][[1]]
    zhongian <- region[[i]][[2]]
    stru <- region[[i]][[3]]
    tem[qidian:zhongian, 1] <- stru
    return(tem)
  }
  #无score的region数据,disorder.data指传入的json格式数据，len为序列长度，用于拼接的数据框，
  region_noscore <- function(disorder.data, len) {
    #新建tem数据框用于存储regions数据
    tem <- data.frame(matrix(NA, len, 1))
    #循环赋值
    for (i in 1:length(disorder.data$regions)) {
      tem <- regions(disorder.data$regions, i, tem)
    }
    
    #返回tem
    return(tem)
  }
  
  #有score的region数据
  region_score <- function(disorder.data, len) {
    #新建tem数据框用于存储regions数据，由于有评分，tem为两列
    tem <- data.frame(matrix(NA, len, 2))
    #将打分数据添加进tem
    tem[, 2] <- as.character(disorder.data$scores)
    #有的情况下regions全空
    if (length(disorder.data$regions) == 0) {
      
    } else{
      #将region数据添加进表
      for (i in 1:length(disorder.data$regions)) {
        tem <- regions(disorder.data$regions, i, tem)
      }
    }
    
    #返回tem
    return(tem)
  }
  
  #加载httr包，用于获取请求网页的得到的JSON文件
  library("httr")
  #发送get请求，获取consensus内容，里面有disorder的3类数据
  data <-
    content(GET(
      paste("http://mobidb.bio.unipd.it/ws/", UniAC, "/consensus", sep = "")
    ))
  #判断请求是否成功，不成功返回两个节点的json结果
  if (!is.null(data$code)) {
    cat("consensus网页请求失败！\t")
  }
  #sequ为mobidb数据库中对应的蛋白质序列
  sequ <-
    content(GET(
      paste("http://mobidb.bio.unipd.it/ws/", UniAC, "/uniprot", sep = "")
    ))$sequence
  #同样进行判断请求是否成功
  if (is.null(sequ)) {
    cat("uniprot网页请求失败！\t")
  }
  
  
  #每个蛋白在mobidb数据不一，新建统一格式空数据框并入原表，根据request请求返回内容填入内容
  table <- data.frame(matrix(NA, nchar(sequ), 14))
  colnames(table) <- c(
    "position",
    "residue",
    "mobidb.disorder.full.regions",
    "mobidb.disorder.db.regions",
    "mobidb.disorder.derived.bfactor.regions",
    "mobidb.disorder.derived.bfactor.score",
    "mobidb.disorder.derived.mobile.regions",
    "mobidb.disorder.derived.mobile.score",
    "mobidb.disorder.derived.missing_residues.regions",
    "mobidb.disorder.derived.missing_residues.score",
    "mobidb.disorder.derived.full.regions",
    "mobidb.disorder.predictors.simple.regions",
    "mobidb.disorder.predictors.mobidb-lite.regions",
    "mobidb.disorder.predictors.mobidb-lite.score"
  )
  table[, 1] <- 1:nchar(sequ)
  table[, 2] <- strsplit(sequ, "")[[1]]
  
  #disorder.full,这部分内容不一定存在，所以需要提前判断。
  if (!is.null(data$mobidb_consensus$disorder$full)) {
    disorder.full <- data$mobidb_consensus$disorder$full[[1]]
    table[, 3] <-  region_noscore(disorder.full, dim(table)[1])
  }
  
  
  #disorder.db,这部分数据不一定存在，所以需要提前判断。
  
  if (!is.null(data$mobidb_consensus$disorder$db)) {
    disorder.db <- data$mobidb_consensus$disorder$db[[1]]
    table[, 4] <- region_noscore(disorder.db, dim(table)[1])
  }
  
  
  #disorder.derived,这部分不一定存在，需要提前判断。
  if (!is.null(data$mobidb_consensus$disorder$derived)) {
    disorder.derived <- data$mobidb_consensus$disorder$derived
    for (i in 1:length(disorder.derived)) {
      switch(
        disorder.derived[[i]]$method,
        bfactor = {
          table[, 5:6] <- region_score(disorder.derived[[i]], dim(table)[1])
        },
        mobile = {
          table[, 7:8] <- region_score(disorder.derived[[i]], dim(table)[1])
        },
        missing_residues = {
          table[, 9:10] <- region_score(disorder.derived[[i]], dim(table)[1])
        },
        full = {
          table[, 11] <- region_noscore(disorder.derived[[i]], dim(table)[1])
        }
      )
    }
  }
  
  
  #disorder.predictior,这部分不一定存在，需要判断
  if (!is.null(data$mobidb_consensus$disorder$predictors)) {
    disorder.predictors <- data$mobidb_consensus$disorder$predictors
    for (i in 1:length(disorder.predictors)) {
      switch(disorder.predictors[[i]]$method,
             "mobidb-lite" = {
               table[, 13:14] <-
                 region_score(disorder.predictors[[i]], dim(table)[1])
             },
             simple = {
               table[, 12] <-
                 region_noscore(disorder.predictors[[i]], dim(table)[1])
             })
    }
    
  }
  
  
  write.table(
    table,
    outputfile,
    quote = F,
    sep = "\t",
    row.names = FALSE
  )

}

#mobidbzhushi("P35609","D:\\idps\\data\\output\\P35609.txt")
