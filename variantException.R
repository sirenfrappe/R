#用于统计由于不同isoform导致的残基突变情况不合的示例

library(rvest)
library(xml2)
library(magrittr)
library(stringr)

UniID <-
  readLines("D:/idps/script/output/UniIDofDieaseassociated.txt")
UniID <- strsplit(UniID, " ")[[1]]
UniID



#Outlier为异常值坐标的记录
Outlier <-
  data.frame("ID" = character(),
             Location = integer(),
             stringsAsFactors = FALSE)

for (i in 75:length(UniID)) {
  cat(i, "\t")
  ID <- UniID[i]
  cat(ID, "\n")
  URL <-
    paste("http://dbptm.mbc.nctu.edu.tw/info.php?id=", ID, sep = "")
  page <- read_html(URL, encoding = "UTF-8")
  sequence <-
    page %>% html_nodes(".table-responsive") %>% extract2(1) %>% html_nodes("td") %>% extract2(10) %>% html_text()
  SitesBySAP <-
    as.data.frame(html_table(
      page %>% html_nodes(".table-responsive") %>% extract2(3) %>% html_nodes("table")
    ))[, c(1, 3, 4)]
  
  for (j in 1:dim(SitesBySAP)[1]) {
    PosVariant <- strsplit(SitesBySAP[j, 2], " ")[[1]][1]
    cat(j, "\t", PosVariant, "\t")
    #位点附近有多个突变点，则跳过，此种情况手工验证֤
    if (grepl(";", PosVariant)) {
      cat(",附近有多个突变点，跳过\n\t")
      next
    }
    originaResidue <- strsplit(SitesBySAP[j, 3], " ")[[1]][1]
    if (originaResidue != substr(sequence, PosVariant, PosVariant)) {
      cat(ID, "\t", SitesBySAP[j, 1], "\t")
      Outlier[dim(Outlier)[1] + 1, 1] <- ID
      cat(dim(Outlier)[1], "aaaaaaaaaaaaaaaaaaaa\t")
      Outlier[dim(Outlier)[1], 2] <- PosVariant
      cat("添加新异常点\n")
    } else{
      cat("过\n")
    }
  }
  cat("一条结束\n")
}
#去重，并添加单ptm位点附近有多个突变点情况
Outlier <- unique(Outlier)
shoudongOutlier <-
  data.frame(
    ID = c("ABI3_HUMAN", "TLR10_HUMAN", "VTDB_HUMAN", "VTDB_HUMAN"),
    Location = c(209, 312, 451, 455)
  )
Outlier <- rbind(Outlier, shoudongOutlier)
write.table(Outlier,
            file = "D:/idps/script/output/yichangzhi.txt",
            quote = FALSE,
            row.names = FALSE)
