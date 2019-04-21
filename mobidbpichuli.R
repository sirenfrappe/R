#读取uniprot上下载的id对照表,此处只需要前两列ID和AC
#id对照表中的检索用id是之前得到的output文件夹下UniIDofDieaseassociated.txt
id <- read.table("D:/idps/data/dbPTM/uniprot-yourlist_M20190417216DA2B77BFBD2E6699CA9B6D1C41EB202F439B.tab",header = T,sep = "\t",stringsAsFactors = FALSE)[,1:2]
for (i in 1:143) {
  source("D:/idps/script/func-mobidbZhushi.R")
  mulu <- paste("D:/idps/script/output/mobidb/",id[i,1],".txt",sep = "")
  mobidbZhushi(id[i,1],id[i,2],mulu)
  cat(i,"\t")
}
