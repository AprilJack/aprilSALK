webG <- function(counts,resSig,org,db,name,gsea.df=NULL){
  #Use WebG
  #sapply(db,webgX <- function(x) {
  #  webG(counts,resSig,org,x,name)})
  #C: the number of reference genes in the category
  #O: the number of genes in the user gene list and also in the category
  #E: The expected number in the category
  #R: ratio of enrichment
  library("WebGestaltR")#,lib.loc = "/gpfs/analyses/april/April/lib/R/library/")
  print("#########################################################")
  print(sprintf("Starting GO Analysis for - %s using %s", name, db))
  print("#########################################################")
  if(!is.null(resSig$logFC)){
    if(dim(resSig[resSig$logFC > 0,])[1] > 10){
      print("...UP ORA EdgeR...")
      Up.ORA <- WebGestaltR(enrichMethod="ORA",
                            organism=org,
                            enrichDatabase=db,
                            interestGeneType="refseq_mrna",
                            referenceGeneType="refseq_mrna",
                            projectName=sprintf("%s.%s.Up.ORA",name,db),
                            referenceGene=row.names(counts),
                            interestGene=row.names(resSig[resSig$logFC > 0,]))

      if (!is.null(dim(Up.ORA)[1])) {
        print("...UP ORA EdgeR Graph...")
        Up.ORA$labs <- paste(Up.ORA$geneSet,Up.ORA$description, sep=": ")
        Up.ORA.Sig <- Up.ORA[Up.ORA$FDR < 0.05,]
        Up.ORA.Sig <- Up.ORA.Sig[order(Up.ORA.Sig$FDR, decreasing= F),]
        Up.ORA.Sig$labs <- factor(Up.ORA.Sig$labs, levels = Up.ORA.Sig$labs[order(-log(Up.ORA.Sig$FDR))])
        p.Up <- ggplot(data=na.omit(Up.ORA.Sig[1:50,]), aes(x=labs, y=(-log(FDR)), fill=R)) +
          geom_bar(stat="identity",col="black") +
          scale_fill_gradient(low="lightskyblue", high="tomato") +
          geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3) +
          labs(title=sprintf("ORA Top 50 - %s Up",name), subtitle="FDR < 5%")  +
          coord_flip() +
          theme_minimal() +
          theme(axis.title.y=element_blank())

        ggsave(filename = sprintf("%s_%s_ORA_UP.pdf",name,db), plot = p.Up, width = 11, height = 11)

      }
    }

    if(dim(resSig[resSig$logFC < 0,])[1] > 10){
      print("...DOWN ORA EdgeR...")
      Down.ORA <- WebGestaltR(enrichMethod="ORA",
                              organism=org,
                              enrichDatabase =db,
                              interestGeneType="refseq_mrna",
                              referenceGeneType="refseq_mrna",
                              projectName=sprintf("%s.%s.Down.ORA",name,db),
                              referenceGene=row.names(counts),
                              interestGene=row.names(resSig[resSig$logFC < 0,]))

      if (!is.null(dim(Down.ORA)[1])) {
        print("...DOWN ORA EdgeR Graph...")
        Down.ORA$labs <- paste(Down.ORA$geneSet,Down.ORA$description, sep=": ")
        Down.ORA.Sig <- Down.ORA[Down.ORA$FDR < 0.05,]
        Down.ORA.Sig <- Down.ORA.Sig[order(Down.ORA.Sig$FDR, decreasing= F),]
        Down.ORA.Sig$labs <- factor(Down.ORA.Sig$labs, levels = Down.ORA.Sig$labs[order(-log(Down.ORA.Sig$FDR))])
        p.Down <- ggplot(data=na.omit(Down.ORA.Sig[1:50,]), aes(x=labs, y=(-log(FDR)), fill=R)) +
          geom_bar(stat="identity",col="black") +
          scale_fill_gradient(low="lightskyblue", high="tomato") +
          geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3)+
          labs(title=sprintf("ORA Top 50 - %s Down",name), subtitle="FDR < 5%")  +
          coord_flip() +
          theme_minimal() +
          theme(axis.title.y=element_blank())

        ggsave(filename = sprintf("%s_%s_ORA_DOWN.pdf",name,db), plot = p.Down, width = 11, height = 11)
      }

    }
    if(dim(resSig)[1] > 10){
      print("...Combined ORA EdgeR...")
      Combined.ORA <- WebGestaltR(enrichMethod="ORA",
                                  organism=org,
                                  enrichDatabase =db,
                                  interestGeneType="refseq_mrna",
                                  referenceGeneType="refseq_mrna",
                                  projectName=sprintf("%s.%s.Combined.ORA",name,db),
                                  referenceGene=row.names(counts),
                                  interestGene=row.names(resSig))

      if (!is.null(dim(Combined.ORA)[1])) {
        print("...Combined ORA EdgeR Graph...")
        Combined.ORA$labs <- paste(Combined.ORA$geneSet,Combined.ORA$description, sep=": ")
        Combined.ORA.Sig <- Combined.ORA[Combined.ORA$FDR < 0.05,]
        Combined.ORA.Sig <- Combined.ORA.Sig[order(Combined.ORA.Sig$FDR, decreasing= F),]
        Combined.ORA.Sig$labs <- factor(Combined.ORA.Sig$labs, levels = Combined.ORA.Sig$labs[order(-log(Combined.ORA.Sig$FDR))])
        p.Combined <- ggplot(data=na.omit(Combined.ORA.Sig[1:50,]), aes(x=labs, y=(-log(FDR)), fill=R)) +
          geom_bar(stat="identity",col="black") +
          scale_fill_gradient(low="lightskyblue", high="tomato") +
          geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3)+
          labs(title=sprintf("ORA Top 50 - %s Combined",name), subtitle="FDR < 5%")  +
          coord_flip() +
          theme_minimal() +
          theme(axis.title.y=element_blank())

        ggsave(filename = sprintf("%s_%s_ORA_Combined.pdf",name,db), plot = p.Combined, width = 11, height = 11)
      }
    }
    if(dim(gsea.df)[1] > 10){
    print("...GSEA EdgeR...")
    Up.GSEA <- tryCatch(WebGestaltR(enrichMethod="GSEA",
                                    organism=org,
                                    enrichDatabase =db,
                                    interestGeneType="refseq_mrna",
                                    projectName=sprintf("%s.%s.GSEA",name,db),
                                    interestGene=gsea.df ), error = function(err) {
                                      # error handler picks up where error was generated
                                      print(paste("MY_ERROR:  ",err))
                                      Up.GSEA=NULL
                                      return(Up.GSEA)
                                    })


    if (!is.null(dim(Up.GSEA)[1])) {
      print("...GSEA EdgeR Graph...")
      Up.GSEA$labs <- paste(Up.GSEA$geneSet,Up.GSEA$description, sep=": ")
      Up.GSEA.Sig <- Up.GSEA[Up.GSEA$FDR < 0.05,]
      Up.GSEA.Sig <- Up.GSEA.Sig[order(Up.GSEA.Sig$FDR, decreasing= F),]
      Up.GSEA.Sig$labs <- factor(Up.GSEA.Sig$labs, levels = Up.GSEA.Sig$labs[order(-log(Up.GSEA.Sig$FDR))])
      p.Up.GSEA <- ggplot(data=na.omit(Up.GSEA.Sig[1:50,]), aes(x=labs, y=(-log(FDR)), fill=ES)) +
        geom_bar(stat="identity",col="black") +
        scale_fill_gradient(low="lightskyblue", high="tomato") +
        geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3)+
        labs(title=sprintf("GSEA Top 50 - %s Up",name), subtitle="FDR < 5%")  +
        coord_flip() +
        theme_minimal() +
        theme(axis.title.y=element_blank())

      ggsave(filename = sprintf("%s_%s_GSEA_EdgeR.pdf",name,db), plot = p.Up.GSEA, width = 11, height = 11)
    }
    }

  }else if(!is.null(resSig$meth.diff)){ #MethylSeq section
    if(length(unique(resSig[resSig$meth.diff > 25,]$Nearest.Refseq))[1] > 10){
      print("...Hyper-methylated ORA MethylSeq...")
      Up.ORA <- WebGestaltR(enrichMethod="ORA",
                            organism=org,
                            enrichDatabase=db,
                            interestGeneType="refseq_mrna",
                            referenceGeneType="refseq_mrna",
                            projectName=sprintf("%s.%s.hyper.ORA",name,db),
                            referenceGene=unique(counts$Nearest.Refseq),
                            interestGene=unique(resSig[resSig$meth.diff > 25,]$Nearest.Refseq))

      if (!is.null(dim(Up.ORA)[1])) {
        print("...Hyper ORA MethylSeq Graph...")
        Up.ORA$labs <- paste(Up.ORA$geneSet,Up.ORA$description, sep=": ")
        Up.ORA.Sig <- Up.ORA[Up.ORA$FDR < 0.05,]
        Up.ORA.Sig <- Up.ORA.Sig[order(Up.ORA.Sig$FDR, decreasing= F),]
        Up.ORA.Sig$labs <- factor(Up.ORA.Sig$labs, levels = Up.ORA.Sig$labs[order(-log(Up.ORA.Sig$FDR))])
        p.Up <- ggplot(data=na.omit(Up.ORA.Sig[1:50,]), aes(x=labs, y=(-log(FDR)), fill=R)) +
          geom_bar(stat="identity",col="black") +
          scale_fill_gradient(low="lightskyblue", high="tomato") +
          geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3) +
          labs(title=sprintf("ORA Top 50 - %s Hyper",name), subtitle="FDR < 5%")  +
          coord_flip() +
          theme_minimal() +
          theme(axis.title.y=element_blank())

        ggsave(filename = sprintf("%s_%s_ORA_hyper.pdf",name,db), plot = p.Up, width = 11, height = 11)

      }
    }
    if(length(unique(resSig[resSig$meth.diff < -25,]$Nearest.Refseq))[1] > 10){
      print("...Hypo-methylated ORA MethylSeq...")
      Down.ORA <- WebGestaltR(enrichMethod="ORA",
                              organism=org,
                              enrichDatabase=db,
                              interestGeneType="refseq_mrna",
                              referenceGeneType="refseq_mrna",
                              projectName=sprintf("%s.%s.hypo.ORA",name,db),
                              referenceGene=unique(counts$Nearest.Refseq),
                              interestGene=unique(resSig[resSig$meth.diff < -25,]$Nearest.Refseq))

      if (!is.null(dim(Down.ORA)[1])) {
        print("...Hyper ORA MethylSeq Graph...")
        Down.ORA$labs <- paste(Down.ORA$geneSet,Down.ORA$description, sep=": ")
        Down.ORA.Sig <- Down.ORA[Down.ORA$FDR < 0.05,]
        Down.ORA.Sig <- Down.ORA.Sig[order(Down.ORA.Sig$FDR, decreasing= F),]
        Down.ORA.Sig$labs <- factor(Down.ORA.Sig$labs, levels = Down.ORA.Sig$labs[order(-log(Down.ORA.Sig$FDR))])
        p.Down <- ggplot(data=na.omit(Down.ORA.Sig[1:50,]), aes(x=labs, y=(-log(FDR)), fill=R)) +
          geom_bar(stat="identity",col="black") +
          scale_fill_gradient(low="lightskyblue", high="tomato") +
          geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3) +
          labs(title=sprintf("ORA Top 50 - %s Hypo",name), subtitle="FDR < 5%")  +
          coord_flip() +
          theme_minimal() +
          theme(axis.title.y=element_blank())

        ggsave(filename = sprintf("%s_%s_ORA_hypo.pdf",name,db), plot = p.Down, width = 11, height = 11)

      }
    }
  } else if(!is.null(resSig$log2FoldChange)){
    ##### DESeq2 section #######
    if(dim(resSig[resSig$log2FoldChange > 0,])[1] > 10){
      print("...UP ORA DESeq2...")
      Up.ORA <- WebGestaltR(enrichMethod="ORA",
                            organism=org,
                            enrichDatabase=db,
                            interestGeneType="refseq_mrna",
                            referenceGeneType="refseq_mrna",
                            projectName=sprintf("%s.%s.Up.ORA",name,db),
                            referenceGene=row.names(counts),
                            interestGene=row.names(resSig[resSig$log2FoldChange > 0,]))
      if (!is.null(dim(Up.ORA)[1])) {
        print("...UP ORA DESeq2 Graph...")
        Up.ORA$labs <- paste(Up.ORA$geneSet,Up.ORA$description, sep=": ")
        Up.ORA.Sig <- Up.ORA[Up.ORA$FDR < 0.05,]
        Up.ORA.Sig <- Up.ORA.Sig[order(Up.ORA.Sig$FDR, decreasing= F),]
        Up.ORA.Sig$labs <- factor(Up.ORA.Sig$labs, levels = Up.ORA.Sig$labs[order(-log(Up.ORA.Sig$FDR))])
        p.Up <- ggplot(data=na.omit(Up.ORA.Sig[1:50,]), aes(x=labs, y=(-log(FDR)), fill=R)) +
          geom_bar(stat="identity",col="black") +
          scale_fill_gradient(low="lightskyblue", high="tomato") +
          geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3) +
          labs(title=sprintf("ORA Top 50 - %s Up",name), subtitle="FDR < 5%")  +
          coord_flip() +
          theme_minimal() +
          theme(axis.title.y=element_blank())

        ggsave(filename = sprintf("%s_%s_ORA_UP.pdf",name,db), plot = p.Up, width = 11, height = 11)

      }

    }
    if(dim(resSig[resSig$log2FoldChange < 0,])[1] > 10){
      print("...DOWN ORA DESeq2...")
      Down.ORA <- WebGestaltR(enrichMethod="ORA",
                              organism=org,
                              enrichDatabase =db,
                              interestGeneType="refseq_mrna",
                              referenceGeneType="refseq_mrna",
                              projectName=sprintf("%s.%s.Down.ORA",name,db),
                              referenceGene=row.names(counts),
                              interestGene=row.names(resSig[resSig$log2FoldChange < 0,]))
      if (!is.null(dim(Down.ORA)[1])) {
        print("...DOWN ORA DESeq2 Graph...")
        Down.ORA$labs <- paste(Down.ORA$geneSet,Down.ORA$description, sep=": ")
        Down.ORA.Sig <- Down.ORA[Down.ORA$FDR < 0.05,]
        Down.ORA.Sig <- Down.ORA.Sig[order(Down.ORA.Sig$FDR, decreasing= F),]
        Down.ORA.Sig$labs <- factor(Down.ORA.Sig$labs, levels = Down.ORA.Sig$labs[order(-log(Down.ORA.Sig$FDR))])
        p.Down <- ggplot(data=na.omit(Down.ORA.Sig[1:50,]), aes(x=labs, y=(-log(FDR)), fill=R)) +
          geom_bar(stat="identity",col="black") +
          scale_fill_gradient(low="lightskyblue", high="tomato") +
          geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3)+
          labs(title=sprintf("ORA Top 50 - %s Down",name), subtitle="FDR < 5%")  +
          coord_flip() +
          theme_minimal() +
          theme(axis.title.y=element_blank())

        ggsave(filename = sprintf("%s_%s_ORA_DOWN.pdf",name,db), plot = p.Down, width = 11, height = 11)
      }
    }
    if(dim(resSig)[1] > 10){
      print("...Combined ORA DESeq2...")
      Combined.ORA <- WebGestaltR(enrichMethod="ORA",
                                  organism=org,
                                  enrichDatabase=db,
                                  interestGeneType="refseq_mrna",
                                  referenceGeneType="refseq_mrna",
                                  projectName=sprintf("%s.%s.Combined.ORA",name,db),
                                  referenceGene=row.names(counts),
                                  interestGene=row.names(resSig))

      if (!is.null(dim(Combined.ORA)[1])) {
        print("...Combined ORA DESeq2 Graph...")
        Combined.ORA$labs <- paste(Combined.ORA$geneSet,Combined.ORA$description, sep=": ")
        Combined.ORA.Sig <- Combined.ORA[Combined.ORA$FDR < 0.05,]
        Combined.ORA.Sig <- Combined.ORA.Sig[order(Combined.ORA.Sig$FDR, decreasing= F),]
        Combined.ORA.Sig$labs <- factor(Combined.ORA.Sig$labs, levels = Combined.ORA.Sig$labs[order(-log(Combined.ORA.Sig$FDR))])
        p.Combined <- ggplot(data=na.omit(Combined.ORA.Sig[1:50,]), aes(x=labs, y=(-log(FDR)), fill=R)) +
          geom_bar(stat="identity",col="black") +
          scale_fill_gradient(low="lightskyblue", high="tomato") +
          geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3) +
          labs(title=sprintf("ORA Top 50 - %s Combined",name), subtitle="FDR < 5%")  +
          coord_flip() +
          theme_minimal() +
          theme(axis.title.y=element_blank())

        ggsave(filename = sprintf("%s_%s_ORA_Combined.pdf",name,db), plot = p.Combined, width = 11, height = 11)

      }
    }
    if(dim(gsea.df)[1] > 10){
    print("...GSEA DESeq2...")
    Up.GSEA <- tryCatch(WebGestaltR(enrichMethod="GSEA",
                                    organism=org,
                                    enrichDatabase =db,
                                    interestGeneType="refseq_mrna",
                                    projectName=sprintf("%s.%s.GSEA",name,db),
                                    interestGene=gsea.df ), error = function(err) {
                                      # error handler picks up where error was generated
                                      print(paste("MY_ERROR:  ",err))
                                      Up.GSEA=NULL
                                      return(Up.GSEA)
                                    })
    if (!is.null(dim(Up.GSEA)[1])) {
      print("...GSEA DESeq2 Graph...")
      Up.GSEA$labs <- paste(Up.GSEA$geneSet,Up.GSEA$description, sep=": ")
      Up.GSEA.Sig <- Up.GSEA[Up.GSEA$FDR < 0.05,]
      Up.GSEA.Sig <- Up.GSEA.Sig[order(Up.GSEA.Sig$FDR, decreasing= F),]
      Up.GSEA.Sig$labs <- factor(Up.GSEA.Sig$labs, levels = Up.GSEA.Sig$labs[order(-log(Up.GSEA.Sig$FDR))])
      p.Up.GSEA <- ggplot(data=na.omit(Up.GSEA.Sig[1:50,]), aes(x=labs, y=(-log(FDR)), fill=ES)) +
        geom_bar(stat="identity",col="black") +
        scale_fill_gradient(low="lightskyblue", high="tomato") +
        geom_text(aes(label=signif(FDR, digits = 2)), hjust=1.1, color="white", size=3)+
        labs(title=sprintf("GSEA Top 50 - %s Up",name), subtitle="FDR < 5%")  +
        coord_flip() +
        theme_minimal() +
        theme(axis.title.y=element_blank())

      ggsave(filename = sprintf("%s_%s_GSEA.pdf",name,db), plot = p.Up.GSEA, width = 11, height = 11)
    }

    }
  } else{
    print("Please enter a valid DESeq2, EdgeR, or MethylSeq object.")
  }

}
