makeDesignAmandine <- function(design){
  
  if("PE.sequencing" %in% colnames(design))
  {
    design <- design[design$PE.sequencing=="R1",]
  }
  design <- design[design$Sequencing.round=="1",]
  
  row.names(design) <- paste(vapply(strsplit(as.character(design$fastq.file),"[_-]+"), `[`, 1, FUN.VALUE=character(1)),
                                    vapply(strsplit(as.character(design$fastq.file),"[_-]+"), `[`, 2, FUN.VALUE=character(1)),
                                    sep="_")
  design <- design[order(row.names(design)),]
  return(design)
}