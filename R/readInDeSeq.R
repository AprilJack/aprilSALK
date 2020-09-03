readInDeSeq <- function(pth,design){
  library(DESeq2)
  setwd(pth)
  raw <- read.delim("raw.txt",header=T,stringsAsFactors = F)
  row.names(raw) <- raw[,1]
  raw <- raw[,-c(1:8)]
  
  colnames(raw) <- paste(vapply(strsplit(colnames(raw),"[_.-]+"), `[`, 2, FUN.VALUE=character(1)),
                       vapply(strsplit(colnames(raw),"[_.-]+"), `[`, 3, FUN.VALUE=character(1)),
                       vapply(strsplit(colnames(raw),"[_.-]+"), `[`, 4, FUN.VALUE=character(1)),sep="_")
  
  if("PE.sequencing" %in% colnames(design))
  {
    design <- design[design$PE.sequencing=="R1",]
  }
  row.names(design) <- paste(vapply(strsplit(as.character(design$fastq.file),"[_-]+"), `[`, 1, FUN.VALUE=character(1)),
                             vapply(strsplit(as.character(design$fastq.file),"[_-]+"), `[`, 2, FUN.VALUE=character(1)),
                             vapply(strsplit(as.character(design$fastq.file),"[_-]+"), `[`, 3, FUN.VALUE=character(1)),sep="_")
  
  design$samp.uniq <- paste(vapply(strsplit(as.character(design$fastq.file),"[_-]+"), `[`, 1, FUN.VALUE=character(1)),
                             vapply(strsplit(as.character(design$fastq.file),"[_-]+"), `[`, 2, FUN.VALUE=character(1)),
                             vapply(strsplit(as.character(design$fastq.file),"[_-]+"), `[`, 3, FUN.VALUE=character(1)),sep="_")
  
  
  design$samp <- paste(vapply(strsplit(as.character(design$fastq.file),"[_-]+"), `[`, 1, FUN.VALUE=character(1)),
                       vapply(strsplit(as.character(design$fastq.file),"[_-]+"), `[`, 2, FUN.VALUE=character(1)),sep="_")
  
  design <- design[order(row.names(design)),]
  
  if(any(!(row.names(design) %in% colnames(raw))) == TRUE){
    raw[[row.names(design)[!(row.names(design) %in% colnames(raw))]]] <- 0 
    raw <- raw[,order(colnames(raw))]
  }
  #create deseq object 
  dds <- DESeqDataSetFromMatrix(countData = raw[,row.names(design)] * 2, colData = design, design = as.formula("~ Group"))
  #collapse resequencings
  dds <- collapseReplicates(dds, dds$samp, dds$samp.uniq, renameCols = TRUE)

  return(dds)
}