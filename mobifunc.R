#添加mobdbi注释
mobidbzhushi <- function(UniID,UniAC,outputfile){
  
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
  region_noscore <- function(disorder.data,len,type){
    #新建tem数据框用于存储regions数据
    tem <- data.frame(matrix(NA, len, 1))
    #循环赋值
    for (i in 1:length(disorder.data$regions)) {
      tem <- regions(disorder.data$regions, i, tem)
    }
    #更改列名，合并。
    colnames(tem) <- paste("disorder.",type,".",disorder.data$method,sep = "")
    table <- cbind(table, tem)
    return(table)
  }
  
  #有score的region数据
  region_score <- function(disorder.data,len,type){
    #新建tem数据框用于存储regions数据，由于有评分，tem为两列
    tem <- data.frame(matrix(NA,len,2))
    #tem的列名，一列regions、一列scores
    colnames(tem) <-
      c(
        paste("disorder",type, disorder.data$method, "regions", sep = "."),
        paste("disorder",type, disorder.data$method, "scores", sep = ".")
      )
    #将打分数据添加进tem
    tem[,2] <- as.character(disorder.data$scores)
    #有的情况下regions全空
    if(length(disorder.data$regions)==0){
      
    }else{
      #将region数据添加进表
      for (i in 1:length(disorder.data$regions)) {
        tem <- regions(disorder.data$regions,i,tem)
      }
    }
   
    table <- cbind(table,tem)
    return(table)
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
  
  #读取处理DBPTM数据库数据得到的表格，命名为table
  table <-
    read.table(
      paste("D:/idps/script/output/UniID/", UniID, ".txt", sep = ""),
      sep = "\t",
      header = T,
      stringsAsFactors = F
    )
  #sequDbptm是dbptm中相应蛋白的序列
  sequDbptm <- paste(as.character(table$residue), collapse = "")
  #判断mobidb数据库与dbPTM数据库中序列是否相同，
  if (sequ != sequDbptm) {
    cat("seq not same\n")
    stop()
  }
  
  #disorder.full,这部分内容不一定存在，所以需要提前判断。
  if(!is.null(data$mobidb_consensus$disorder$full)){
    disorder.full <- data$mobidb_consensus$disorder$full[[1]]
    table <-  region_noscore(disorder.full,dim(table)[1],"full")
  }
  
  
  #disorder.db,这部分数据不一定存在，所以需要提前判断。
  
  if(!is.null(data$mobidb_consensus$disorder$db)){
    disorder.db <- data$mobidb_consensus$disorder$db[[1]]
    table <- region_noscore(disorder.db,dim(table)[1],"db")
  }
  
  
  #disorder.derived,这部分不一定存在，需要提前判断。
  if(!is.null(data$mobidb_consensus$disorder$derived)){
    disorder.derived <- data$mobidb_consensus$disorder$derived
    table.list <-  lapply(disorder.derived,function(x){
      if(length(x)==3){
        region_score(x,dim(table)[1],"derived")
      }else{
        if(length(x)==2){
          region_noscore(x,dim(table)[1],"derived")
        }
      }
      
    })
    #得到的是list，需要处理一下,table.list长度为1~4不等
    if(length(table.list)==1){
      table <- table.list[[1]]
    }
    if(length(table.list)==2){
      table <- merge(table.list[[1]],table.list[[2]])
    }
    if(length(table.list)==3){
      table <- merge(merge(table.list[[1]],table.list[[2]]),table.list[[3]])
    }
    if(length(table.list)==4){
      table <- merge(merge(merge(table.list[[1]],table.list[[2]]),table.list[[3]]),table.list[[4]])
    }
  }
  #上面处理完的table会乱序
  table <- table[order(table$location),]
  #disorder.predictior,这部分不一定存在，需要判断
  if(!is.null(data$mobidb_consensus$disorder$predictors)){
    disorder.predictors <- data$mobidb_consensus$disorder$predictors
    table.list <-  lapply(disorder.predictors,function(x){
      if(length(x)==4){
        region_score(x,dim(table)[1],"predictors")
      }else{
        if(length(x)==2){
          region_noscore(x,dim(table)[1],"predictors")
        }
      }
      
    })
    #得到的是list，需要处理一下
    if(length(table.list)==1){
      table <- table.list[[1]]
    }
    if(length(table.list)==2){
      table <- merge(table.list[[1]],table.list[[2]])
    }
  }
  #上面处理完的table会乱序
  table <- table[order(table$location),]
  
  
  write.table(
    table,
    outputfile,
    quote = F,
    sep = "\t",
    row.names = FALSE
  )
  
}

#mobidbzhushi("TDG_HUMAN","Q13569","D:/idps/script/output/mobidb/TDG_HUMAN(Mobidb).txt")