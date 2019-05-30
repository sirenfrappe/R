library(seqinr)
seq <- read.fasta("D:\\idps\\data\\Uniprot_Sequence\\A0A0H2URD6.fasta", seqtype = "AA")[[1]]
seq
library("Biostrings")
#sequ <- Biostrings::readAAStringSet("D:\\idps\\data\\Uniprot_Sequence\\P45984.fasta")

path <- "D:/idps/data/Uniprot_Sequence/"
fileName <- dir(path)
filePath <- sapply(fileName, function(x) {
  paste(path, x, sep = "")
})
source("D:/idps/script/mobifunc.R")
for (i in 1:length(filePath)) {
  cat(i,"\t")
  sequ <- as.character(Biostrings::readAAStringSet(filePath[i])[[1]])
  id <- strsplit(names(filePath[i]), "\\.")[[1]][1]
  temp <- data.frame(matrix(NA, nchar(sequ), 9))
  temp[, 1] <- 1:nchar(sequ)
  colnames(temp)[1] <- "position"
  temp[, 2] <- unlist(strsplit(sequ, ""))
  colnames(temp)[2] <- "residue"
  write.table(temp,
        paste("D:/idps/script/output/UniID/", id, ".txt", sep = ""),
        quote = FALSE,
        row.names = FALSE,
        sep = "\t")

  cat("1st\t")

  mobidbzhushi(id, id, paste("D:/idps/data/output/", id, ".txt", sep = ""))
  cat("2nd\t")
  temp <- read.table(paste("D:/idps/data/output/", id, ".txt", sep = ""), 
      sep = "\t",
      header = T,
      stringsAsFactors = F,
      quote = "")
    temp <- temp[,-3:-9]
   write.table(temp,
        paste("D:/idps/data/output/", id, ".txt", sep = ""),
        quote = FALSE,
        row.names = FALSE,
        sep = "\t")
   cat(i)
}
