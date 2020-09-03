#' A Ortholog Fetch Function
#'
#' This function gets human ortholog of refseq_mrna ID in rownames.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' getHumanOrtho()

getHumanOrtho <- function(df="data frame with refseq mrna id's to convert",
                          df_refseq_mrna="column in dataframe that contains the refseq mrna ids",
                          refseq_mrna="name of column in dataframe that contains the refseq mrna ids as a string"){
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package \"biomaRt\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  library("biomaRt")
  ensemblMouse = useMart("ensembl",dataset="mmusculus_gene_ensembl")

  test1 <- getBM(attributes=c('ensembl_gene_id','hsapiens_homolog_ensembl_gene','hsapiens_homolog_associated_gene_name'),
                 filters = 'refseq_mrna',
                 values = df_refseq_mrna,
                 mart = ensemblMouse)


  test <- getBM(attributes=c('ensembl_gene_id','refseq_mrna','external_gene_name'),
                filters = 'refseq_mrna',
                values = df_refseq_mrna,
                mart = ensemblMouse)
  test2<-merge(df,test,by.y="refseq_mrna",by.x=refseq_mrna,all.x=T)
  test3 <-merge(test2,test1,by="ensembl_gene_id",all.x=T)
  return(test3)

}
