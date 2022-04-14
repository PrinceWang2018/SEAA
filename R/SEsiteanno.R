#' Annotation of Splicing Sites Acquiring Filtered Splicing Efficiency
#' @description This function helps to annotate the splicing sites with clusterProfiler. The gene symbol, REFSEQ, ENSEMBL are shown.
#' @param SEresultlist a large list created by SEcalculation
#' @param efficiency_5ss_3ss_nona_inf_reduct a data frame created by SEfilter
#' @param species species for the data to be annotated: hs or mm
#' @return a data frame
#' @export
#'
#' @examples
#' \donttest{
#' SEannotaionresult<-SEsiteanno(SEresultlist, efficiency_5ss_3ss_nona_inf_reduct)
#' }
SEsiteanno<-function(SEresultlist, efficiency_5ss_3ss_nona_inf_reduct, species = "hs"){
  message("[Step1] Annotating gene symbol, NM_id, and Ensemble_id.")
  efficiency_5ss_con_anno<-efficiency_5ss_3ss_nona_inf_reduct
  efficiency_5ss_con_anno_genes<-tidyr::separate(efficiency_5ss_con_anno,c("GeneID"),into = c("NM","NM_ID","INTRON","INTRON_Num","0","chr","site","f","5or3"),sep = "_",remove = FALSE)
  efficiency_5ss_con_anno_genes2<-tidyr::unite(efficiency_5ss_con_anno_genes,"NM_id",c("NM","NM_ID"),sep = "_",remove = TRUE)
  efficiency_5ss_con_anno_genes4<-tidyr::unite(efficiency_5ss_con_anno_genes2,"Intron_num",c("INTRON","INTRON_Num"),sep = "",remove = TRUE)
  efficiency_5ss_con_anno_genes4<-efficiency_5ss_con_anno_genes4[,-4:-7]
  gene <- as.character(efficiency_5ss_con_anno_genes4$NM_id)
  gene2<-gsub("..$","",gene)
  efficiency_5ss_con_anno_genes4$NM_id2<-gene2
  if (species=="hs") {
    gene.df <- clusterProfiler::bitr(gene2, fromType = "REFSEQ",
                                     toType = c("SYMBOL","ENSEMBL"),    #switching REFSEQ into SYMBOL and ENSEMBL
                                     OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  } else if (species=="mm" ) {
    gene.df <- clusterProfiler::bitr(gene2, fromType = "REFSEQ",
                                     toType = c("SYMBOL","ENSEMBL"),    #switching REFSEQ into SYMBOL and ENSEMBL
                                     OrgDb = org.Mm.eg.db::org.Mm.eg.db)
  } else {
    message("Error: imported species were not corrected or supported.")
    gene.df <-""
  }
  efficiency_5ss_con_anno_genes5<-merge(efficiency_5ss_con_anno_genes4,gene.df, by.x="NM_id2", by.y="REFSEQ")
  efficiency_5ss_con_anno_genes5<-efficiency_5ss_con_anno_genes5[,-1]
  message("[Step2] Calculating splicing efficiency fold change.")
  cols<-colnames(efficiency_5ss_con_anno_genes5)
  new_cols<-c(cols[1:4],cols[7:8],cols[5:6])
  efficiency_5ss_con_anno_genes6<-efficiency_5ss_con_anno_genes5[,new_cols]
  efficiency_5ss_con_anno_genes6$FC_s2_vs_s1<-efficiency_5ss_con_anno_genes6[,8] / efficiency_5ss_con_anno_genes6[,7]
  efficiency_5ss_con_anno_genes6$Change<-ifelse(efficiency_5ss_con_anno_genes6$FC_s2_vs_s1>1,"up","down")
  message("[Step3] Adding strand information.")
  efficiency_5ss_3ss<-SEresultlist[["efficiency_5ss_3ss"]]
  efficiency_5ss_con_anno_genes7<-merge(efficiency_5ss_con_anno_genes6,efficiency_5ss_3ss,by.x= "GeneID", by.y = "GeneID", all.x= TRUE, all.y=FALSE)
  efficiency_5ss_con_anno_genes8<-efficiency_5ss_con_anno_genes7[,-15:-17]
  new_cols2<-c(colnames(efficiency_5ss_con_anno_genes8)[1:4],colnames(efficiency_5ss_con_anno_genes8)[11:14],colnames(efficiency_5ss_con_anno_genes8)[5:10])
  efficiency_5ss_con_anno_genes8<-efficiency_5ss_con_anno_genes8[,new_cols2]
  message("[Step4] Merging reads counts data.")
  fc_count_spliced<-SEresultlist[["fc_count_spliced"]]
  fc_count_spliced$GeneID<-row.names(fc_count_spliced)
  fc_count_unspliced<-SEresultlist[["fc_count_unspliced"]]
  fc_count_unspliced$GeneID<-row.names(fc_count_unspliced)
  sample_list_show<-SEresultlist[["sample_list_show"]]
  efficiency_5ss_con_anno_genes9<-merge(efficiency_5ss_con_anno_genes8,fc_count_spliced,by.x= "GeneID", by.y = "GeneID", all.x= TRUE, all.y=FALSE)
  efficiency_5ss_con_anno_genes9<-merge(efficiency_5ss_con_anno_genes9,fc_count_unspliced,by.x= "GeneID", by.y = "GeneID", all.x= TRUE, all.y=FALSE)
  col_change<-colnames(efficiency_5ss_con_anno_genes9)
  text_s1_se<-paste0(sample_list_show[1],"_splicing_efficiency")
  text_s2_se<-paste0(sample_list_show[2],"_splicing_efficiency")
  text_sc1_se<-paste0(sample_list_show[1],"_spliced_reads")
  text_sc2_se<-paste0(sample_list_show[2],"_spliced_reads")
  text_usc1_se<-paste0(sample_list_show[1],"_unspliced_reads")
  text_usc2_se<-paste0(sample_list_show[2],"_unspliced_reads")
  col_change_new<-c(col_change[1:10],text_s1_se,text_s2_se,col_change[13:14],text_sc1_se,text_sc2_se,text_usc1_se,text_usc2_se)
  colnames(efficiency_5ss_con_anno_genes9)<-col_change_new
  new_cols3<-c(colnames(efficiency_5ss_con_anno_genes9)[1:10],colnames(efficiency_5ss_con_anno_genes9)[15:18],colnames(efficiency_5ss_con_anno_genes9)[11:14])
  efficiency_5ss_con_anno_genes9<-efficiency_5ss_con_anno_genes9[,new_cols3]
  message("[Step5] Exporting the final dataframe.")
  utils::write.table(efficiency_5ss_con_anno_genes9,file = "./Processed_Splicing_efficiency_annotation.txt",sep = "\t",row.names = FALSE)
  return(efficiency_5ss_con_anno_genes9)
}

