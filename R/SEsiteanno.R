#' Annotation of Splicing Sites Acquiring Filtered Splicing Efficiency
#' @description This function helps to annotate the splicing sites with clusterProfiler. The gene symbol, REFSEQ, ENSEMBL are shown.
#' @param efficiency_5ss_3ss_nona_inf_reduct a data frame created by SEfilter
#' @return a data frame
#' @export
#'
#' @examples
#' \donttest{
#' SEannotaionresult<-SEsiteanno(efficiency_5ss_3ss_nona_inf_reduct)
#' }
SEsiteanno<-function(efficiency_5ss_3ss_nona_inf_reduct){
  ##Annotation
  efficiency_5ss_con_anno<-efficiency_5ss_3ss_nona_inf_reduct
  efficiency_5ss_con_anno_genes<-tidyr::separate(efficiency_5ss_con_anno,GeneID,into = c("NM","NM_ID","INTRON","INTRON_Num","0","chr","site","f","5or3"),sep = "_",remove = FALSE)
  efficiency_5ss_con_anno_genes2<-tidyr::unite(efficiency_5ss_con_anno_genes,"NM_id",c("NM","NM_ID"),sep = "_",remove = TRUE)
  efficiency_5ss_con_anno_genes4<-tidyr::unite(efficiency_5ss_con_anno_genes2,"Intron_num",c("INTRON","INTRON_Num"),sep = "",remove = TRUE)
  efficiency_5ss_con_anno_genes4<-efficiency_5ss_con_anno_genes4[,-4:-7]
  gene <- as.character(efficiency_5ss_con_anno_genes4$NM_id)
  gene2<-gsub("..$","",gene)
  efficiency_5ss_con_anno_genes4$NM_id2<-gene2
  gene.df <- clusterProfiler::bitr(gene2, fromType = "REFSEQ",
                  toType = c("SYMBOL","ENSEMBL"),    #switching REFSEQ into SYMBOL and ENSEMBL
                  OrgDb = org.Hs.eg.db)
  efficiency_5ss_con_anno_genes5<-merge(efficiency_5ss_con_anno_genes4,gene.df, by.x="NM_id2", by.y="REFSEQ")
  efficiency_5ss_con_anno_genes5<-efficiency_5ss_con_anno_genes5[,-1]

  cols<-colnames(efficiency_5ss_con_anno_genes5)
  new_cols<-c(cols[1:4],cols[7:8],cols[5:6])
  efficiency_5ss_con_anno_genes6<-efficiency_5ss_con_anno_genes5[,new_cols]
  #col_change<-colnames(efficiency_5ss_con_anno_genes6)
  #col_change_new<-c(col_change[1:10],"sample1_splicing_efficiency","sample2_splicing_efficiency","sample1_unspliced_reads","sample2_unspliced_reads","sample1_spliced_reads","sample2_spliced_reads")
  #colnames(efficiency_5ss_con_anno_genes6)<-col_change_new
  efficiency_5ss_con_anno_genes6$FC_s2_vs_s1<-efficiency_5ss_con_anno_genes6[,7] / efficiency_5ss_con_anno_genes6[,8]

  efficiency_5ss_con_anno_genes6$Change<-ifelse(efficiency_5ss_con_anno_genes6$FC_s2_vs_s1>1,"up","down")
  utils::write.table(efficiency_5ss_con_anno_genes6,file = "./Processed_Splicing_efficiency_annotation.txt",sep = "\t",row.names = FALSE)
  return(efficiency_5ss_con_anno_genes6)
}

utils::globalVariables(c("org.Hs.eg.db","GeneID"))
