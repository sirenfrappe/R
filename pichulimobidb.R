id <- read.table("D:/idps/data/dbPTM/uniprot-yourlist_M20190417216DA2B77BFBD2E6699CA9B6D1C41EB202F439B.tab",header = T,sep = "\t",stringsAsFactors = FALSE)[,1:2]
source("D:/idps/script/mobifunc.R")
for (i in 142:143) {
  mulu <- paste("D:/idps/script/output/mobidb/",id[i,1],".txt",sep = "")
  mobidbzhushi(id[i,1],id[i,2],mulu)
  cat(i,"\t")
}
#P16860 	P00747 	Q13569 	Q9Y4G6 	P02788 	Q14191 
#6ge 

#6 ABI3_HUMAN	Q9P2A4  
#16 APOE_HUMAN	P02649
#19 ATM_HUMAN	Q13315
#46 FAS_HUMAN	P49327
#91 NUCB2_HUMAN	P80303
#99 PPR3B_HUMAN	Q86XI6
#100 PRC2C_HUMAN	Q9Y520
#141 WHRN_HUMAN	Q9P202

#8个序列不相同

pros <- c(6,16,19,46,91,99,100,141)
for (i in pros) {
  dbptm <- read.table(paste("D:/idps/script/output/UniID/", id[i,1], ".txt", sep = ""),
                          sep = "\t",
                          header = T,
                          stringsAsFactors = F)
  #dbptm中序列
  dbptm.seq <- paste(as.character(dbptm[,2]),collapse = "")
  #
  cat(">",id[i,1],"dbptm\n",sep = "",dbptm.seq,file = paste("D:/idps/script/output/diffSeqInDbptm&mobidb/",id[i,1],".fasta"))
  library(httr)
  sequ.mobidb <-
    content(GET(
      paste("http://mobidb.bio.unipd.it/ws/", id[i,2], "/uniprot", sep = "")
    ))$sequence
  cat("\n",">",id[i,1],"mobidb\n",sep = "",sequ.mobidb,file = paste("D:/idps/script/output/diffSeqInDbptm&mobidb/",id[i,1],".fasta"),append=TRUE)
}

#使用clustal进行比对,结果在diffSeqInDbptm&mobidb文件夹下

