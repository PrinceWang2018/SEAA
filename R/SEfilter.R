#' Splicing efficiency Filtering
#' @description  Unconventional splicing efficiency will be removed and reads counts less than min_counts will be filtered.
#' @param SEresultlist a large list created by SEcalculation
#' @param min_counts a numeric value
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' efficiency_5ss_3ss_nona_inf_reduct<-SEfilter(SEresultlist,min_counts = 5)
#'
  SEfilter<-function(SEresultlist,min_counts = 5){
    efficiency_5ss_3ss<-SEresultlist[["efficiency_5ss_3ss"]]
    fc_count_spliced<-SEresultlist[["fc_count_spliced"]]
    fc_count_unspliced<-SEresultlist[["fc_count_unspliced"]]
    sample_list<-SEresultlist[["sample_list"]]
    sample_list_show<-SEresultlist[["sample_list_show"]]

    if (min_counts == 0) {
      message("Note: min_counts was imported as 0 and Inf was replaced by the max splicing efficiency value in each group.")
      efficiency_5ss_3ss_nona<-stats::na.omit(efficiency_5ss_3ss)
      efficiency_5ss_3ss_nona_inf<- efficiency_5ss_3ss_nona
      for (i in 1:length(sample_list)) {
        i
        tmp_max<-efficiency_5ss_3ss_nona_inf[,i+6]
        tmp_max[tmp_max==Inf]<-NA
        efficiency_5ss_3ss_nona_inf[,i+6][efficiency_5ss_3ss_nona_inf[,i+6]==Inf]<- max(tmp_max,na.rm = TRUE)
      }
      efficiency_5ss_3ss_nona_inf_reduct<-efficiency_5ss_3ss_nona_inf[,-2:-6]
    } else  {
      message("Note: min_counts is imported as ", paste(min_counts)," and reads counts less than ", paste(min_counts)," will be filtered." )
      fc_count_spliced_con_tmp<-subset(fc_count_spliced,fc_count_spliced[,1]>=min_counts&fc_count_spliced[,2]>min_counts)
      fc_count_spliced_con_tmp$GeneID<-row.names(fc_count_spliced_con_tmp)
      fc_count_unspliced_con_tmp<-subset(fc_count_unspliced,fc_count_unspliced[,1]>=min_counts&fc_count_unspliced[,2]>min_counts)
      fc_count_unspliced_con_tmp$GeneID<-row.names(fc_count_unspliced_con_tmp)
      filtered_GeneID<- merge(fc_count_unspliced_con_tmp,fc_count_spliced_con_tmp, by.x = "GeneID",by.y = "GeneID", all.x = FALSE, all.y = FALSE)
      efficiency_con<-merge(efficiency_5ss_3ss,filtered_GeneID, by.x = "GeneID",by.y = "GeneID",all.x = FALSE, all.y = FALSE)
      efficiency_con_reduct<-efficiency_con[,-9:-12]
      efficiency_5ss_3ss_nona_inf_reduct<-efficiency_con_reduct[,-2:-6]
      colnames(efficiency_5ss_3ss_nona_inf_reduct)<- c("GeneID",sample_list_show)
    }
    return(efficiency_5ss_3ss_nona_inf_reduct)
  }
