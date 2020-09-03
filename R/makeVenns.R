#' A function to make Venn Diagrams
#'
#' This function allows you to make 2, 3, or 4 circle Venn diagrams. Outputs images and annotated lists.
#' @param no defaults.
#' @keywords venn
#' @export
#' @examples
#' makeVenns()

makeVenns <- function(overlap="list of elements to overlap",
                      nam="desired file base name",
                      cats="name of each element of your list",
                      geneAnno="gene annotation to add to overlap",
                      byy="what to merge the overlaps and gene annotation by"){
  if(length(cats) == 2){
    venn.diagram(x =overlap,category.names = cats,
                 filename = sprintf("%s_venn.png",nam),
                 main = sprintf("Overlap of Genes %s",nam),imagetype = 'png')

    names(overlap) <- cats
    VennDiagram_obj <- get.venn.partitions(x = overlap)

    make_venn_csv(VennDiagram_obj,nam,geneAnno,byy)

    tx <- shQuote(paste('python /home/april/venn.py ',sprintf("%s_venn_proportional.pdf",nam), ' [','\"',
                        VennDiagram_obj$..count..[3],'\"', ',','\"',
                        VennDiagram_obj$..count..[2],'\"',',','\"',
                        VennDiagram_obj$..count..[1],'\"',"]",' ["', cats[1], '", "', cats[2], '"]',sep=""),type="sh")
    cat(tx)

    system(command=paste('/gpfs/tools/anaconda/bin/python3.5 /home/april/venn.py ',
                         sprintf("%s_venn_proportional.pdf",nam), ' [','"',
                         VennDiagram_obj$..count..[3],'"', ',','"',
                         VennDiagram_obj$..count..[2],'"',',','"',
                         VennDiagram_obj$..count..[1],'"',"]",' ["', cats[1], '","', cats[2], '"]',sep=""))
  }else if(length(cats) == 3){
    venn.diagram(x = overlap,
                 category.names = cats,
                 filename = sprintf("%s_venn.png",nam),
                 main = sprintf("Overlap of Genes %s",nam),imagetype = 'png')

    names(overlap) <- cats
    VennDiagram_obj <- get.venn.partitions(x = overlap)

    make_venn_csv(VennDiagram_obj,nam,geneAnno,byy)

    tx <- shQuote(paste('python /home/april/venn.py ',sprintf("%s_venn_proportional.pdf",nam), ' [','\"',
                        VennDiagram_obj$..count..[7],'\"', ',','\"',
                        VennDiagram_obj$..count..[6],'\"', ',','\"',
                        VennDiagram_obj$..count..[5],'\"', ',','\"',
                        VennDiagram_obj$..count..[4],'\"', ',','\"',
                        VennDiagram_obj$..count..[3],'\"', ',','\"',
                        VennDiagram_obj$..count..[2],'\"',',','\"',
                        VennDiagram_obj$..count..[1],'\"',"]",' ["', cats[1], '", "', cats[2],'", "', cats[3], '"]',sep=""),type="sh")
    cat(tx)

    system(command=paste('/gpfs/tools/anaconda/bin/python3.5 /home/april/venn.py ', sprintf("%s_venn_proportional.pdf",nam),' [','"',
                         VennDiagram_obj$..count..[7],'"', ',','"',
                         VennDiagram_obj$..count..[6],'"', ',','"',
                         VennDiagram_obj$..count..[5],'"', ',','"',
                         VennDiagram_obj$..count..[4],'"', ',','"',
                         VennDiagram_obj$..count..[3],'"', ',','"',
                         VennDiagram_obj$..count..[2],'"', ',','"',
                         VennDiagram_obj$..count..[1],'"',"]",' ["', cats[1], '","', cats[2],'","', cats[3], '"]',sep=""))
  }else if(length(cats) == 4){
    venn.diagram(x = overlap,
                 category.names = cats,
                 filename = sprintf("%s_venn.png",nam),
                 main = sprintf("Overlap of Genes %s",nam),imagetype = 'png')

    names(overlap) <- cats
    VennDiagram_obj <- get.venn.partitions(x = overlap)

    make_venn_csv(VennDiagram_obj,nam,geneAnno,byy)
  }else{
    print("Please enter an appropriate number of categories for a 2, 3 or 4 circle Venn.")
  }
}
