#' A Venn Diagram Function
#'
#' This function takes in a Venn diagram object from VennDiagram package, and makes a CSV file of the overlaps with annotation.
#' @param no defaults
#' @keywords venn
#' @export
#' @examples
#' make_venn_csv()


make_venn_csv<- function(VennDiagram_obj,nam,ann,byy){
  if(length(VennDiagram_obj$..count..) == 7) {
    if(length(as.data.frame(VennDiagram_obj$..values..[1])[,1]) > 0){
      set1 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[1])[,1],
                         set = rep(VennDiagram_obj$..set..[1],length(VennDiagram_obj$..values..[1])))} else{
                           set1 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[1],length(VennDiagram_obj$..values..[1])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[2])[,1]) > 0){
      set2 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[2])[,1],
                         set = rep(VennDiagram_obj$..set..[2],length(VennDiagram_obj$..values..[2])))} else{
                           set2 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[2],length(VennDiagram_obj$..values..[2])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[3])[,1]) > 0){
      set3 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[3])[,1],
                         set = rep(VennDiagram_obj$..set..[3],length(VennDiagram_obj$..values..[3])))} else{
                           set3 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[3],length(VennDiagram_obj$..values..[3])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[4])[,1]) > 0){
      set4 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[4])[,1],
                         set = rep(VennDiagram_obj$..set..[4],length(VennDiagram_obj$..values..[4])))} else{
                           set4 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[4],length(VennDiagram_obj$..values..[4])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[5])[,1]) > 0){
      set5 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[5])[,1],
                         set = rep(VennDiagram_obj$..set..[5],length(VennDiagram_obj$..values..[5])))} else{
                           set5 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[5],length(VennDiagram_obj$..values..[5])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[6])[,1]) > 0){
      set6 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[6])[,1],
                         set = rep(VennDiagram_obj$..set..[6],length(VennDiagram_obj$..values..[6])))} else{
                           set6 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[6],length(VennDiagram_obj$..values..[6])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[7])[,1]) > 0){
      set7 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[7])[,1],
                         set = rep(VennDiagram_obj$..set..[7],length(VennDiagram_obj$..values..[7])))} else{
                           set7 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[7],length(VennDiagram_obj$..values..[7])))}

    venn_overlap <- rbind(set1,set2,set3,set4,set5,set6,set7)
    res <- merge(venn_overlap,ann,by.x=1,by.y=byy, all.x=T)

    write.csv(x = res, file = sprintf("%s_venn.csv",nam), quote = F, row.names = F)
    openxlsx::write.xlsx(x = res,file = sprintf("%s_venn.xlsx",nam),firstRow=T,keepNA=T)
  }
  else if(length(VennDiagram_obj$..count..) == 3){
    if(length(as.data.frame(VennDiagram_obj$..values..[1])[,1]) > 0){
      set1 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[1])[,1],
                         set = rep(VennDiagram_obj$..set..[1],length(VennDiagram_obj$..values..[1])))} else{
                           set1 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[1],length(VennDiagram_obj$..values..[1])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[2])[,1]) > 0){
      set2 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[2])[,1],
                         set = rep(VennDiagram_obj$..set..[2],length(VennDiagram_obj$..values..[2])))} else{
                           set2 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[2],length(VennDiagram_obj$..values..[2])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[3])[,1]) > 0){
      set3 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[3])[,1],
                         set = rep(VennDiagram_obj$..set..[3],length(VennDiagram_obj$..values..[3])))} else{
                           set3 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[3],length(VennDiagram_obj$..values..[3])))}
    venn_overlap <- rbind(set1,set2,set3)
    res <- merge(venn_overlap,ann,by.x=1,by.y=byy, all.x=T)

    write.csv(x = res, file = sprintf("%s_venn.csv",nam), quote = F, row.names = F)
    openxlsx::write.xlsx(x = res,file = sprintf("%s_venn.xlsx",nam),firstRow=T,keepNA=T)
  } else if(length(VennDiagram_obj$..count..) == 15){
    if(length(as.data.frame(VennDiagram_obj$..values..[1])[,1]) > 0){
      set1 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[1])[,1],
                         set = rep(VennDiagram_obj$..set..[1],length(VennDiagram_obj$..values..[1])))} else{
                           set1 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[1],length(VennDiagram_obj$..values..[1])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[2])[,1]) > 0){
      set2 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[2])[,1],
                         set = rep(VennDiagram_obj$..set..[2],length(VennDiagram_obj$..values..[2])))} else{
                           set2 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[2],length(VennDiagram_obj$..values..[2])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[3])[,1]) > 0){
      set3 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[3])[,1],
                         set = rep(VennDiagram_obj$..set..[3],length(VennDiagram_obj$..values..[3])))} else{
                           set3 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[3],length(VennDiagram_obj$..values..[3])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[4])[,1]) > 0){
      set4 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[4])[,1],
                         set = rep(VennDiagram_obj$..set..[4],length(VennDiagram_obj$..values..[4])))} else{
                           set4 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[4],length(VennDiagram_obj$..values..[4])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[5])[,1]) > 0){
      set5 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[5])[,1],
                         set = rep(VennDiagram_obj$..set..[5],length(VennDiagram_obj$..values..[5])))} else{
                           set5 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[5],length(VennDiagram_obj$..values..[5])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[6])[,1]) > 0){
      set6 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[6])[,1],
                         set = rep(VennDiagram_obj$..set..[6],length(VennDiagram_obj$..values..[6])))} else{
                           set6 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[6],length(VennDiagram_obj$..values..[6])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[7])[,1]) > 0){
      set7 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[7])[,1],
                         set = rep(VennDiagram_obj$..set..[7],length(VennDiagram_obj$..values..[7])))} else{
                           set7 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[7],length(VennDiagram_obj$..values..[7])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[8])[,1]) > 0){
      set8 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[8])[,1],
                         set = rep(VennDiagram_obj$..set..[8],length(VennDiagram_obj$..values..[8])))} else{
                           set8 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[8],length(VennDiagram_obj$..values..[8])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[9])[,1]) > 0){
      set9 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[9])[,1],
                         set = rep(VennDiagram_obj$..set..[9],length(VennDiagram_obj$..values..[9])))} else{
                           set9 <- data.frame(gene = "NA",
                                              set = rep(VennDiagram_obj$..set..[9],length(VennDiagram_obj$..values..[9])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[10])[,1]) > 0){
      set10 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[10])[,1],
                          set = rep(VennDiagram_obj$..set..[10],length(VennDiagram_obj$..values..[10])))} else{
                            set10 <- data.frame(gene = "NA",
                                                set = rep(VennDiagram_obj$..set..[10],length(VennDiagram_obj$..values..[10])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[11])[,1]) > 0){
      set11 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[11])[,1],
                          set = rep(VennDiagram_obj$..set..[11],length(VennDiagram_obj$..values..[11])))} else{
                            set11 <- data.frame(gene = "NA",
                                                set = rep(VennDiagram_obj$..set..[11],length(VennDiagram_obj$..values..[11])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[12])[,1]) > 0){
      set12 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[12])[,1],
                          set = rep(VennDiagram_obj$..set..[12],length(VennDiagram_obj$..values..[12])))} else{
                            set12 <- data.frame(gene = "NA",
                                                set = rep(VennDiagram_obj$..set..[12],length(VennDiagram_obj$..values..[12])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[13])[,1]) > 0){
      set13 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[13])[,1],
                          set = rep(VennDiagram_obj$..set..[13],length(VennDiagram_obj$..values..[13])))} else{
                            set13 <- data.frame(gene = "NA",
                                                set = rep(VennDiagram_obj$..set..[13],length(VennDiagram_obj$..values..[13])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[14])[,1]) > 0){
      set14 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[14])[,1],
                          set = rep(VennDiagram_obj$..set..[14],length(VennDiagram_obj$..values..[14])))} else{
                            set14 <- data.frame(gene = "NA",
                                                set = rep(VennDiagram_obj$..set..[14],length(VennDiagram_obj$..values..[14])))}

    if(length(as.data.frame(VennDiagram_obj$..values..[15])[,1]) > 0){
      set15 <- data.frame(gene = as.data.frame(VennDiagram_obj$..values..[15])[,1],
                          set = rep(VennDiagram_obj$..set..[15],length(VennDiagram_obj$..values..[15])))} else{
                            set15 <- data.frame(gene = "NA",
                                                set = rep(VennDiagram_obj$..set..[15],length(VennDiagram_obj$..values..[15])))}

    venn_overlap <- rbind(set1,set2,set3,set4,set5,set6,set7,set8,set9,set10,set11,set12,set13,set14,set15)
    res <- merge(venn_overlap,ann,by.x=1,by.y=byy, all.x=T)

    write.csv(x = res, file = sprintf("%s_venn.csv",nam), quote = F, row.names = F)
    openxlsx::write.xlsx(x = res,file = sprintf("%s_venn.xlsx",nam),firstRow=T,keepNA=T)

  }else{
    print("Invalid Venn diagram object.")
  }

}
