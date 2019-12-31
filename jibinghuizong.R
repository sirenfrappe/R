#用于统计每个疾病对应多少蛋白

#加载包
library(rvest)
library(xml2)
library(magrittr)
library(stringr)


#读取保存的目录本地网页，将信息进行汇总
for (i in 1:15) {
  html_url <-
    paste("D:/idps/data/dbPTM/muluWeb/dbPTM", i, ".html", sep = "")
  page <- read_html(html_url, encoding = "UTF-8")
  table <-
    data.frame(html_table(
      page %>% html_nodes(".table-responsive") %>% html_nodes("table"),
      header = TRUE
    ))
  for (j in (dim(table)[1] + 1):2) {
    num <-
      length(
        page %>% html_nodes(".table-responsive")  %>% html_nodes("tr") %>% extract2(j) %>% html_nodes("li")
      )
    if (num != 1) {
      for (m in 1:num) {
        disease <-
          page %>% html_nodes(".table-responsive")  %>% html_nodes("tr") %>% extract2(j) %>% html_nodes("li") %>% html_text()
        line <- table[j - 1, ]
        line$Related.Disease <- disease[m]
        table <- rbind(table, line)
        if (m == num) {
          table <- table[-j + 1, ]
        }
      }
    }
  }
  if (i == 1) {
    data <- table
  } else{
    data <- rbind(data, table)
  }
}

#去除ID不存在的蛋白信息
#读取之前处理过的ID信息
UniID <-
  readLines("D:/idps/script/output/pro122.txt")
UniID <- as.character(strsplit(UniID, " ")[[1]])
data <-  subset(data, is.na(match(data$ID, UniID)) == F)
#去除重复行
data <- unique(data)
data <- data[-c(4, 5, 7)]
data$Related.Disease <- as.factor(data$Related.Disease)
data$ID <- as.factor(data$ID)
data$PTM.Type <- as.factor(data$PTM.Type)
#按照相关位点个数，进行统计
library(plyr)
data1 <- count(data$Related.Disease)
data1 <- data1[order(data1[, 2], decreasing = T), ]
write.table(
  data1,
  "D:/idps/script/output/geshuweidian.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)
