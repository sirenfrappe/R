#ID为UniProt ID ，mulu为输出文件目录，均为字符串格式
func <- function(ID, mulu) {
  #包加载
  library(rvest)
  library(xml2)
  library(magrittr)
  library(stringr)
  
  #网页读取
  HTMLurl <-
    paste("http://dbptm.mbc.nctu.edu.tw/info.php?id=", ID, sep = "")
  page <- read_html(HTMLurl, encoding = "UTF-8")
  cat("网页读取完成\t")
  #序列获取，并分隔
  sequ <-
    page %>% html_nodes(".table-responsive") %>% extract2(1) %>% html_nodes("td") %>% extract2(10) %>% html_text()
  sequ <- as.data.frame(strsplit(sequ, ""))
  colnames(sequ) <- "residue"
  sequ <- data.frame(location = seq(1, dim(sequ)[1]), sequ)
  
  #位点的获取
  eps <-
    as.data.frame(html_table(
      page %>% html_nodes(".table-responsive") %>% extract2(2) %>% html_nodes("table")
    ))[, c(1, 2, 3)]
  eps <-
    data.frame(
      eps[, 1:2],
      SubstratePeptides = substr(eps[, 3], 1, 15),
      SecondaryStructure = substr(eps[, 3], 16, 30),
      stringsAsFactors = FALSE
    )
  
  #疾病相关数据
  DaPS <-
    as.data.frame(html_table(
      page %>% html_nodes(".table-responsive") %>% extract2(3) %>% html_nodes("table")
    ))[, c(1, 2, 3, 4, 6)]
  #有的突变点相关疾病不止一个，直接读取数据时未做出分别
  #处理多疾病情况
  for (i in 1:dim(DaPS)[1]) {
    numOfDisease <-
      length(
        page %>% html_nodes(".table-responsive") %>% extract2(3) %>% html_nodes("tr") %>% extract2(i +
                                                                                                     1) %>% html_nodes("li")
      )
    if (numOfDisease == 1) {
      next
    } else{
      DaPS[i, 5] <-
        paste(
          page %>% html_nodes(".table-responsive") %>% extract2(3) %>% html_nodes("tr") %>% extract2(i +
                                                                                                       1) %>% html_nodes("li") %>% extract(1:numOfDisease) %>% html_text(),
          collapse = ";"
        )
    }
  }
  #只保留突变结果
  DaPS$Residue.Change <-
    data.frame(sapply(strsplit(DaPS$Residue.Change, " "), function(x) {
      x[3]
    }), stringsAsFactors = FALSE)
  colnames(DaPS)[c(3, 4)] <- c("VariantPosition", "ResidueChange")
  
  #不能用数字开头命名，使用ac 来命名数据框
  UniID <- matrix(rep(0, dim(sequ)[1] * 9), ncol = 9)
  UniID <- as.data.frame(UniID)
  colnames(UniID) <-
    c(
      "location",
      "residue",
      "Modification",
      "seq",
      "secStruc",
      "RelatedDisease",
      "VariantPosition",
      "ResidueChange",
      "DaPSbySAP"
    )
  UniID[, c(1, 2)] <- sequ
  #eps整合
  for (i in 1:dim(eps)[1]) {
    UniID[eps[i, 1], 3:5] <- eps[i, 2:4]
  }
  #DaPS整合
  for (i in 1:dim(DaPS)[1]) {
    UniID[DaPS[i, 1], 6:7] <- DaPS[i, c(5, 3)]
    #有的附近变异点不止一个：ABI3_HUMAN,IQGA2_HUMAN,FA5_HUMAN,MSHR_HUMAN,TLR10_HUMAN,VTDB_HUMAN，共六个
    change <- as.character(DaPS[i, 4])
    pos <- strsplit(DaPS[i, 3], " ")[[1]][1]
    #突变
    if (!grepl(";", pos)) {
      if (UniID[pos, 8] == 0) {
        UniID[pos, 8] <- change
      }
      #突变点附近疾病相关PTM位点位置
      if (UniID[pos, 9] == 0) {
        UniID[pos, 9] <- DaPS[i, 1]
      } else{
        UniID[pos, 9] <- paste(UniID[pos, 9], DaPS[i, 1], sep = ",")
      }
    } else{
      #附近有两个突变点，突变情况相同，共三个蛋白FA5_HUMAN，TLR10_HUMAN，ABI3_HUMAN
      if (!grepl(";", change)) {
        pos1 <- strsplit(pos, ";")[[1]][1]
        pos2 <- strsplit(pos, ";")[[1]][2]
        if (UniID[pos1, 8] == 0) {
          UniID[pos1, 8] <- change
        }
        if (UniID[pos2, 8] == 0) {
          UniID[pos2, 8] <- change
        }
        if (UniID[pos1, 9] == 0) {
          UniID[pos1, 9] <- DaPS[i, 1]
        } else{
          UniID[pos1, 9] <- paste(UniID[pos1, 9], DaPS[i, 1], sep = ",")
        }
        if (UniID[pos2, 9] == 0) {
          UniID[pos2, 9] <- DaPS[i, 1]
        } else{
          UniID[pos2, 9] <- paste(UniID[pos2, 9], DaPS[i, 1], sep = ",")
        }
        
      } else{
        if (ID == "IQGA2_HUMAN") {
          pos1 <- strsplit(pos, ";")[[1]][1]
          pos2 <- strsplit(pos, ";")[[1]][2]
          change1 <- strsplit(change, ";")[[1]][1]
          change2 <- strsplit(change, ";")[[1]][2]
          UniID[pos1, 8] <- change1
          UniID[pos2, 8] <- change2
          if (UniID[pos1, 9] == 0) {
            UniID[pos1, 9] <- DaPS[i, 1]
          } else{
            UniID[pos1, 9] <- paste(UniID[pos1, 9], DaPS[i, 1], sep = ",")
          }
          if (UniID[pos2, 9] == 0) {
            UniID[pos2, 9] <- DaPS[i, 1]
          } else{
            UniID[pos2, 9] <- paste(UniID[pos2, 9], DaPS[i, 1], sep = ",")
          }
        }
        #此蛋白含有位点附近两个突变位点，其中一个突变位点有两种突变情况
        if (ID == "VTDB_HUMAN") {
          pos1 <- strsplit(pos, ";")[[1]][1]
          pos2 <- strsplit(pos, ";")[[1]][2]
          change1 <- strsplit(change, ";")[[1]][1]
          change2 <- substr(change, 3, 5)
          UniID[pos1, 8] <- change1
          UniID[pos2, 8] <- change2
          if (UniID[pos1, 9] == 0) {
            UniID[pos1, 9] <- DaPS[i, 1]
          } else{
            UniID[pos1, 9] <- paste(UniID[pos1, 9], DaPS[i, 1], sep = ",")
          }
          if (UniID[pos2, 9] == 0) {
            UniID[pos2, 9] <- DaPS[i, 1]
          } else{
            UniID[pos2, 9] <- paste(UniID[pos2, 9], DaPS[i, 1], sep = ",")
          }
        }
        #附近有三个突变位点，其中有两个位点残基相同，突变成为的残基不同
        if (ID == "MSHR_HUMAN") {
          pos1 <- strsplit(pos, ";")[[1]][1]
          pos2 <- strsplit(pos, ";")[[1]][2]
          pos3 <- strsplit(pos, ";")[[1]][3]
          change1 <- strsplit(change, ";")[[1]][1]
          change2 <- strsplit(change, ";")[[1]][2]
          change3 <- strsplit(change, ";")[[1]][3]
          UniID[pos1, 8] <- change1
          UniID[pos2, 8] <- change2
          UniID[pos3, 8] <- change3
          if (UniID[pos1, 9] == 0) {
            UniID[pos1, 9] <- DaPS[i, 1]
          } else{
            UniID[pos1, 9] <- paste(UniID[pos1, 9], DaPS[i, 1], sep = ",")
          }
          if (UniID[pos2, 9] == 0) {
            UniID[pos2, 9] <- DaPS[i, 1]
          } else{
            UniID[pos2, 9] <- paste(UniID[pos2, 9], DaPS[i, 1], sep = ",")
          }
          if (UniID[pos3, 9] == 0) {
            UniID[pos3, 9] <- DaPS[i, 1]
          } else{
            UniID[pos3, 9] <- paste(UniID[pos3, 9], DaPS[i, 1], sep = ",")
          }
        }
      }
    }
  }
  
  write.table(UniID,
              file = mulu,
              quote = FALSE,
              row.names = FALSE)
}


#func("TDG_HUMAN","C:/Users/Admin/Desktop/idps/script/output/TDG_HUMAN.txt")
