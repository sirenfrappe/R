#!/home/q/miniconda2/bin/Rscript
#将psiblast输出的结果转换为aln文件和fasta文件
#PSIBLAST max_target_num 分别设为100、500进行两次blast 
#此脚本将在所在目录根据蛋白ID建立结果文件夹
#protein_id为蛋白ID，protein_file蛋白序列文件
#blast地址可能需要做相应修改

blast_output2aln_Fasta <- function(protein_file,protein_id) {
    cat(protein_id," start.\n")
    if (is.na(match("tidyverse",(.packages())))) {
       library(tidyverse)
    }
    system(str_c("mkdir ","./",protein_id),intern=TRUE)
    filePrefix <- str_c("./",protein_id,"/",protein_id)
    psiblast_100 <- str_c("/mnt/d/DNCON2/ncbi-blast-2.2.25+/bin/psiblast -query ",protein_file," -evalue .001 -inclusion_ethresh .002 -db /mnt/d/DNCON2/databases/nr90-2012/nr90 -num_iterations 3 -num_threads 8 -seg yes -outfmt '7 std qseq sseq stitle' -max_target_seqs 100 -out " , filePrefix , "_100.output")
    psiblast_500 <- str_c("/mnt/d/DNCON2/ncbi-blast-2.2.25+/bin/psiblast -query ",protein_file," -evalue .001 -inclusion_ethresh .002 -db /mnt/d/DNCON2/databases/nr90-2012/nr90 -num_iterations 3 -num_threads 8 -seg yes -outfmt '7 std qseq sseq stitle' -max_target_seqs 500 -out ",filePrefix,"_500.output")
    cat("blast start\n")
    cat("running: ",psiblast_100,"\n")
    system(psiblast_100)
    cat("running: ",psiblast_500,"\n")
    system(psiblast_500)
    cat("blast done.\n")
    #blast运行完毕
    
    #用于blast output转换
    #max_target ：100、500
    tran2 <- function(max_target){
        text <- readLines(str_c(filePrefix,"_",max_target,".output"))
        firstStr <- grep("#",str_sub(text,1,1))
        start <- tail(firstStr,2)[1]+1
        end <- tail(firstStr,2)[2]-1
        cat(text[c(1:6,start:(end+1))],sep = "\n",file = str_c(filePrefix,"_",max_target,"_last.output"))
        system(str_c("/mnt/d/linux/mview/bin/mview -in blast ",filePrefix,"_",max_target,"_last.output -out aln > ",filePrefix,"_",max_target,".aln"))
        #提取id，用于输出序列
        id <- text[start:end] %>%
            str_split("\t",simplify = T) %>%
            .[,2] %>% 
            str_split("\\|",simplify = T) %>%
            .[,2]
        if (is.na(match("reutils",(.packages())))) {
            library("reutils")
        }

        reutils::epost(id,'protein') %>%
            efetch(retmode = "text" , rettype = "fasta",outfile = str_c(filePrefix,"_",max_target,".fasta"))
    }
    tran2(100)
    tran2(500)
    cat(protein_id," finished.\n")
}

if (is.na(match("tidyverse",(.packages())))) {
       library(tidyverse)
    }

path <- "/mnt/d/DNCON2/protein" 
fileNames <- dir(path)
filePath <- sapply(fileNames,function(x){
    str_c(path,x,sep='/')
})
for (i in 124:143) {
   blast_output2aln_Fasta(protein_file = filePath[[i]] , protein_id = str_split(names(filePath[i]),"\\.")[[1]][1] )
}
#blast_output2aln_Fasta(protein_file="/mnt/d/DNCON2/protein/1A03_HUMAN.fasta",protein_id="1A03_HUMAN")