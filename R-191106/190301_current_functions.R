myExplore <- function(counts,design,colData1,colData2,formulaDE2,ref,name,anno) {
  suppressPackageStartupMessages(library("DESeq2"))
  suppressPackageStartupMessages(library("stringr"))
  suppressPackageStartupMessages(library("gplots"))
  suppressPackageStartupMessages(library("ggplot2"))
  suppressPackageStartupMessages(library("pheatmap"))
  suppressPackageStartupMessages(library("cluster"))
  suppressPackageStartupMessages(library("VennDiagram"))
  pdf(sprintf("%s_plots.pdf",name),width = 7,height = 7)

  print("#########################################################")
  print("Starting DESeq2 General QC")
  print("#########################################################")

  dds <- DESeqDataSetFromMatrix(countData = counts*2,
                                colData = design,
                                design = as.formula(formulaDE2))
  dds[[colData1]] <- relevel(dds[[colData1]], ref = ref)

  summary(rowSums(counts(dds)))
  hist(log10(rowSums(counts(dds))))
  abline(v=log10(10),col="red")

  keep <- rowSums(counts(dds)) >= 10
  print("Genes passing expression filter:")
  print(table(keep))

  cRef <- counts[,row.names(design[design$type == ref,])]
  cOth <- counts[,row.names(design[design$type != ref,])]

  makeVenns(overlap=list(row.names(cRef[rowSums(cRef) >= 10,]),row.names(cOth[rowSums(cOth) >= 10,])),
            nam=name,
            cats = c(levels(droplevels(design$type))[1], levels(droplevels(design$type))[2]),
            geneAnno = anno,
            byy = 0)

  dds <- dds[keep,]
  dds <- suppressMessages(DESeq(dds))

  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

  pca <- prcomp(t(assay(vsd)))
  summaryPCA <- summary(pca)
  print("PCA summary:")
  print(summaryPCA)

  loadings = as.data.frame(pca$x[])

  print("#########################################################")
  print("Starting QC graphing of VSD transformed raw count data.")
  print("#########################################################")

  print("...graphing PCA...")
  pca12 <- ggplot(data = loadings,
                  aes(x=PC1,
                      y=PC2,
                      color=as.factor(design[[colData1]]),
                      shape=as.factor(design[[colData2]])))
  print(pca12 +
          geom_point(size=17) +
          scale_shape_manual(values=seq(21, 21+length(levels(as.factor(design[[colData2]])))-1, 1)) +
          theme_bw() +
          labs(x=sprintf("PC1 (%s of variance)",paste(round(summaryPCA$importance[2,1]*100, 2), "%", sep="")),
               y=sprintf("PC2 (%s of variance)",paste(round(summaryPCA$importance[2,2]*100, 2), "%", sep=""))) +
          geom_text(aes(label=row.names(design))) +
          ggtitle(sprintf("PCA of normalized raw count data - %s",name)) + labs(color = colData1,shape=colData2))

  print("...graphing Euclidean distance matrix...")
  sampleDists <- dist(t(assay(vsd)))
  df <- as.data.frame(colData(vsd)[,c(colData1,colData2)])
  sampleDistMatrix <- as.matrix(sampleDists)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           show_rownames=FALSE,
           main="Euclidean distance of normalized counts",
           annotation_col=df,annotation_row=df
  )

  print("...graphing heatmap of top 100 expressed genes...")
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:100]

  df <- as.data.frame(colData(vsd)[,c(colData1,colData2)])
  pheatmap(assay(vsd)[select,],
           cluster_rows=TRUE,
           show_rownames=TRUE,
           cluster_cols=TRUE,
           annotation_col=df,
           main=sprintf("%s - Heatmap of Top 100 Expressed Genes",name),
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "correlation",
           clustering_method = "ward.D2",
           treeheight_row=25,
           fontsize = 6,
           color = colorRampPalette(c("green", "black", "red"))(50))
  dev.off()
  }




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

  make_venn_csv(VennDiagram_obj,sprintf("%s_venn.csv",nam),geneAnno,byy)

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

      make_venn_csv(VennDiagram_obj,sprintf("%s_venn.csv",nam),geneAnno,byy)

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

    make_venn_csv(VennDiagram_obj,sprintf("%s_venn.csv",nam),geneAnno,byy)
  }else{
    print("Please enter an appropriate number of categories for a 2, 3 or 4 circle Venn.")
  }
}



rmLaneNum <- function(test){
  #works to remove last 5 characters from string
  n <- nchar(test)
  n <- n-5
  rt <- substr(test,1,n)
  return(rt)
}

rmLastElementStrSplit <- function(test){
  #includes all elments of a stringsplit except the last
  return(paste(strsplit(test,"[._]+")[[1]][1:length(strsplit(test,"[._]+")[[1]])-1],collapse="_"))
}

#gets human ortholog of refseq_mrna ID in rownames
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

make_venn_csv<- function(VennDiagram_obj,filename,ann,byy){
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

    write.csv(x = res, file = filename, quote = F, row.names = F)
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

    write.csv(x = res, file = filename, quote = F, row.names = F)
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

    write.csv(x = res, file = filename, quote = F, row.names = F)

  }else{
    print("Invalid Venn diagram object.")
  }

}


meth25bp <- function(methylseqObj = "methy seq obj",nam = "name of exp",gene.obj = gene.obj,org="organism name for WebGestaltR",db="GO db's to query"){
  library("methylKit")
  library("genomation")
  library("ggplot2")
  library("gplots")
  library("pheatmap")
  library("cluster")
  pdf(sprintf("%s_plots.pdf",nam),width = 7,height = 7)
  print("---------Filtering methyl object")
  filtered.methylseqObj = filterByCoverage(methylseqObj,lo.count=10,lo.perc=NULL,
                                           hi.count=NULL,hi.perc=99.9)

  print("---------Normalizing methyl object")
  normed.methylseqObj <- normalizeCoverage(filtered.methylseqObj)

  print("---------Tiling with 25bp windows with 15bp steps")
  tiles <- tileMethylCounts(normed.methylseqObj,win.size=25,step.size=15)

  print("---------Uniting strands of methyl object")
  meth <- unite(tiles, destrand=FALSE)

  print("---------Calculating samples correlations")
  getCorrelation(meth,plot=T)
  print("---------Clustering samples by correlation")
  clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
  print("---------PCA analysis")
  PCASamples(meth, screeplot=TRUE)
  PCASamples(meth)
  print("---------Calculating methylation percentages")
  percmeth <- percMethylation(meth)

  print("---------Calculating differential methylation")
  myDiff <- calculateDiffMeth(meth,mc.cores=4)
  diffMethPerChr(myDiff,plot=T,qvalue.cutoff=0.01, meth.cutoff=25,exclude=c("chrM","chrX","chrY"))

  print("---------Annotating hyper and hypo methylated using Genomation")
  myDiff.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
  diffAnnHyper <- annotateWithGeneParts(as(myDiff.hyper,"GRanges"),gene.obj)
  genomation::plotTargetAnnotation(diffAnnHyper,precedence=TRUE,main=sprintf("%s\ndifferential methylation\nannotation hyper",nam))

  myDiff.hypo <- getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
  diffAnnHypo <- annotateWithGeneParts(as(myDiff.hypo,"GRanges"),gene.obj)
  genomation::plotTargetAnnotation(diffAnnHypo,precedence=TRUE,main=sprintf("%s\ndifferential methylation\nannotation hypo",nam))

  print("---------Writing out Genomation annotation results")
  diffAnnAll <- annotateWithGeneParts(as(myDiff,"GRanges"),gene.obj)
  annotatedmeth <- data.frame(myDiff,percmeth,diffAnnAll@members)
  write.table(x = annotatedmeth,file = sprintf("%s_annotatedDiffMeth.txt", nam),sep = "\t",quote = F,row.names = F)

  print("---------Annotating hyper and hypo methylated regions using HOMER")
  #annotate methyl regions using HOMER
  input <- sprintf("%s_annotatedDiffMeth.txt", nam)
  output <- sprintf("%s_annotatedDiffMeth_homerAnn.txt", nam)
  system(sprintf("source /gpfs/tools/.bashTools;
     annotatePeaks.pl %s mm10 > %s", input, output))
  annotatedmeth$PeakID <- paste("*-",row.names(annotatedmeth),sep="")
  annotatedHomermeth <- read.delim(output, header=T, stringsAsFactors=F, sep="\t")
  Final_annotatedHomermeth <- merge(annotatedmeth,annotatedHomermeth[,-c(2:4,6:7)],by.x="PeakID",by.y=1)
  write.table(x = Final_annotatedHomermeth,file = sprintf("%s_Final_annotatedDiffMeth.txt",nam),sep = "\t",quote = F,row.names = F)

  diffSigAnn <- subset(Final_annotatedHomermeth, qvalue < 0.01 & abs(meth.diff) > 25)
  getGenesNearCoordinateMouse

  sapply(db,webgX <- function(x) {
    webG(Final_annotatedHomermeth,diffSigAnn,org,x,nam)})

  graphics.off()

  Final_annotatedHomermeth50kb <- getGenesNearCoordinateMouse(Final_annotatedHomermeth,regionSize = 50000)
  write.table(x = Final_annotatedHomermeth50kb,file = sprintf("%s_Final_annotated_50kb.txt",nam),quote = FALSE,append = FALSE,sep = "\t",row.names = F)

  return(Final_annotatedHomermeth50kb)
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste0(toupper(substring(s, 1,1)), tolower(substring(s, 2)))
}

getGenesNearCoordinateMouse <- function(annMethResult, regionSize = 50000){
  if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)) {
    stop("Package \"TxDb.Mmusculus.UCSC.mm10.knownGene\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    stop("Package \"org.Mm.eg.db\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(dplyr)
  library(GenomicRanges)
  library(org.Mm.eg.db)

  #create a data frame of gene ids and symbols
  eid <- as.data.frame(org.Mm.egSYMBOL)
  #creat a coordinate identifier
  annMethResult$coord <- do.call(paste, c(annMethResult[,2:4],sep=":"))
  #create data frame with 50 kb intervals we want to search
  df <- data.frame(chrom=annMethResult$chr,start=as.numeric(annMethResult$start)-regionSize,end=as.numeric(annMethResult$end)+regionSize)
  #make a GRanges object from the data frame
  dfG <- makeGRangesFromDataFrame(df)
  rm(df)
  #merge the GRanges with mm10 UCSC known genes
  dfGgenes <- mergeByOverlaps(genes(TxDb.Mmusculus.UCSC.mm10.knownGene), dfG)
  rm(dfG)
  #create a data frame with the results
  dfGgenesDF <- data.frame(gene_id=dfGgenes$gene_id, methInterval = dfGgenes$dfG)
  rm(dfGgenes)
  #merge the GRanges gene results with symbol
  dfGgenesDFsym <- merge(dfGgenesDF,eid, by="gene_id",all.x=T)
  #add a column of unique coord
  dfGgenesDFsym$coord <- do.call(paste, c(list(dfGgenesDFsym$methInterval.seqnames,dfGgenesDFsym$methInterval.start+regionSize,dfGgenesDFsym$methInterval.end-regionSize),sep=":"))
  #collapse data frame by coord and list all genes found within
  result <- aggregate(symbol ~ coord, data = dfGgenesDFsym, paste, collapse = ",")

  annMethResult50kb <- merge(annMethResult,result,by="coord", all.x=T)
  return(annMethResult50kb)

}


