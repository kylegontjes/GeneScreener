# Libraries 
setwd(snakemake@params[[1]])
library(tidyverse)

# File name
file_name <-  snakemake@input[[1]]
isolate <- gsub(".fna|.fasta","",file_name) %>% gsub("_blast.out","",.) %>% gsub("results/","",.)
database <- readLines(snakemake@input[[2]],warn=F) 
locus_tags <- database %>% subset(grepl(">",database)) %>% gsub(">","",.) %>% str_split(.," ",simplify=T) %>% .[,1]

blast_column_names <- snakemake@params[[2]] %>% gsub("\\\"","",.) %>% str_split(.," ",simplify=T) %>% .[,2:length(.)]
blast_file =  read.delim(file_name,header = F,col.names = blast_column_names) 

curate_blast_entries <- function(loci,blast_file,name){ 
  if(loci %in% blast_file$sseqid){
    blast_entry <- subset(blast_file,sseqid == loci)
    blast_entry$coverage <- (blast_entry$nident / blast_entry$slen) * 100
    blast_entry$coverage <- round(blast_entry$coverage,3)
    blast_entry$isolate_no <- name
    return(blast_entry)
  } else {
    blast_entry <- NULL
  } 
}
  
blast_entry_statistics <- lapply(locus_tags,FUN=curate_blast_entries,blast_file=blast_file,name=isolate) %>% do.call(rbind,.) %>% as.data.frame   
if(nrow(blast_entry_statistics)>0){
  write_delim(blast_entry_statistics,file = snakemake@output[[1]]) 
} else {
   file.create(snakemake@output[[1]])
}   