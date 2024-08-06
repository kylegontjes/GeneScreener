# Libraries
library(tidyverse)

# File name
file_name <-  snakemake@input[[1]]
rowname <- gsub(".fna|.fasta","",file_name) %>% gsub("_blast.out","",.)
database <- readLines(snakemake@input[[2]]) 

get_contig_length <- function(file){
  entries <- 1:length(file) %>% .[lapply(.,"%%",2)==1]
  locus_tags <- file %>% subset(grepl(">",file)) %>% gsub(">","",.) %>% str_split(.," ",simplify=T) %>% .[,1]
  contigs <- 1:length(file) %>% .[lapply(.,"%%",2)==0] 
  get_length <- function(entry,file){
    file[[entry]] %>% nchar
  }
  contig_length <- sapply(contigs,get_length,file)
  results <- data.frame(locus_tag = locus_tags,length = contig_length)
  return(results)
}

database_length <- get_contig_length(database)  

blast_column_names <- c("qseqid" ,"sseqid" ,"slen" ,"length" ,"pident" ,"nident", "mismatch" ,"gapopen", "gaps", 'qstart', "qend" ,"sstart" ,'send' ,"evalue") 
blast_file =  read.delim(file_name,header = F,col.names = blast_column_names) 

get_entry_statistics <- function(loci,database_length,blast_file,reporting_vars){
  locus_length <- subset(database_length,locus_tag==loci) %>% .[,"length"]
  if(loci %in% blast_file$sseqid){
    blast_entry <- subset(blast_file,sseqid == loci)
    blast_entry$coverage <- (blast_entry$nident / blast_entry$slen) * 100
    blast_entry$blast_hit <- ifelse(blast_entry$coverage >40 & blast_entry$pident>80,1,0)
    
    blast_df <- blast_entry %>% select(any_of(reporting_vars)) %>% `colnames<-`(paste0(loci,"_",colnames(.)))
  } else {
    blast_df <-  data.frame(matrix(nrow=1,ncol=length(reporting_vars))) %>% `colnames<-`(paste0(loci,"_",reporting_vars)) 
  }
  return(blast_df)
}
reporting_vars <- c("seqid","blast_hit", "pident","coverage","qstart","qend")

blast_entry_statistics <- lapply(database_length$locus_tag,FUN=get_entry_statistics,database_length=database_length,blast_file=blast_file,reporting_vars=reporting_vars) %>% do.call(cbind,.)  %>% `rownames<-`(rowname)

write_delim(blast_entry_statistics,file = snakemake@output[[1]]) 

get_blast_matrix_entry <- function(entry_statistics){
  matrix <- entry_statistics %>% select_if(grepl("_blast_hit",colnames(entry_statistics))) %>% `colnames<-`(gsub("_blast_hut","",colnames(.)))
  matrix[is.na(matrix)] <- 0
  return(matrix)
}

blast_matrix <- get_blast_matrix_entry(blast_entry_statistics)
write_delim(blast_matrix,file = snakemake@output[[2]])