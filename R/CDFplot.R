#' Cumulative Distribution Function Plotting
#' @description Plot the Cumulative Distribution Function of splicing efficiency.
#' @param SEresultlist a large list created by SEcalculation
#' @param efficiency_5ss_3ss_nona_inf_reduct a dataframe created by SEfilter
#' @param zoom.x a numeric vector with 2 elements
#' @importFrom ggplot2 ggplot xlim ylim theme xlab ylab stat_ecdf geom_text
#' @return a large list with K-S test result
#' @export
#'
#' @examples
#' CDFplot(SEresultlist,efficiency_5ss_3ss_nona_inf_reduct,zoom.x= c(5,6))
#'
CDFplot<-function(SEresultlist,efficiency_5ss_3ss_nona_inf_reduct,zoom.x= c(5,6)){
  sample_list<-SEresultlist[["sample_list"]]
  sample_list_show<-SEresultlist[["sample_list_show"]]
  message("[Step1] Plotting Cumulative Distribution Function for all splicing sites.")
  #all
  #K-S test
  KSTEST_result<-stats::ks.test(efficiency_5ss_3ss_nona_inf_reduct[,2],efficiency_5ss_3ss_nona_inf_reduct[,3])
  D_value<- KSTEST_result[["statistic"]][["D"]]
  label_test2 <- substr(paste0("D=",D_value),1,8)
  p_value<- KSTEST_result[["p.value"]]
  label_test <- ifelse(p_value<0.0001, "p<0.0001",substr(paste0("p=",p_value),1,8))
  #data_transform
  ecdf_exp_con_nona_w_to_L<- stats::reshape(data=efficiency_5ss_3ss_nona_inf_reduct, idvar="GeneID",
                                     varying = sample_list_show,
                                     v.name=c("Splicing_efficiency"),
                                     times=sample_list_show,
                                     #new.row.names = 1:1000,
                                     direction="long")
  ecdf_exp_con_nona_w_to_L$LOG2SE<- log2(ecdf_exp_con_nona_w_to_L$Splicing_efficiency)
  colnames(ecdf_exp_con_nona_w_to_L)<-c("GeneID","Samples","Splicing_efficiency","LOG2SE")
  #ecdf_plot
  xlabelpos<- max(ecdf_exp_con_nona_w_to_L$LOG2SE)-3.5
  cdf_p<-ggplot2::ggplot(ecdf_exp_con_nona_w_to_L,aes(x=LOG2SE, color=Samples))+
    stat_ecdf()+
    labs(title = "Cumulative Distribution Function")+
    xlab("log2Splicing_efficiency")+
    ylab("Cumulative fraction")+
    theme(plot.title=element_text(hjust = 0.5), legend.position = "bottom")+
    geom_text(aes(x=xlabelpos,y=0.1,label = label_test),family = "sans",size=4,colour="#000000")+
    geom_text(aes(x=xlabelpos,y=0.2,label = label_test2),family = "sans",size=4,colour="#000000")+
    ggforce::facet_zoom(xlim = zoom.x, zoom.size = 1)
  ####5ss
  message("[Step2] Plotting Cumulative Distribution Function for 5ss splicing sites.")
  efficiency_5ss_nona_inf_reduct<-subset(efficiency_5ss_3ss_nona_inf_reduct,grepl("5ss",efficiency_5ss_3ss_nona_inf_reduct$GeneID))
  #K-S test
  KSTEST_result_5ss<-stats::ks.test(efficiency_5ss_nona_inf_reduct[,2],efficiency_5ss_nona_inf_reduct[,3])
  D_value_5ss<- KSTEST_result_5ss[["statistic"]][["D"]]
  label_test2_5ss <- substr(paste0("D=",D_value_5ss),1,8)
  p_value_5ss<- KSTEST_result_5ss[["p.value"]]
  label_test_5ss <- ifelse(p_value_5ss<0.0001, "p<0.0001",substr(paste0("p=",p_value_5ss),1,8))
  #data_transform
  ecdf_exp_5ss_con_nona_w_to_L<- stats::reshape(data=efficiency_5ss_nona_inf_reduct, idvar="GeneID",
                                         varying = sample_list_show,
                                         v.name=c("Splicing_efficiency"),
                                         times=sample_list_show,
                                         #new.row.names = 1:1000,
                                         direction="long")
  ecdf_exp_5ss_con_nona_w_to_L$LOG2SE<- log2(ecdf_exp_5ss_con_nona_w_to_L$Splicing_efficiency)
  colnames(ecdf_exp_5ss_con_nona_w_to_L)<-c("GeneID","Samples","Splicing_efficiency","LOG2SE")
  #ecdf_plot
  xlabelpos<- max(ecdf_exp_5ss_con_nona_w_to_L$LOG2SE)-3.5
  cdf_p_5ss<- ggplot2::ggplot(ecdf_exp_5ss_con_nona_w_to_L,aes(x=LOG2SE, color=Samples))+
    stat_ecdf()+
    labs(title = "Cumulative Distribution Function")+
    xlab("log2 5ss Splicing_efficiency")+
    ylab("Cumulative fraction")+
    theme(plot.title=element_text(hjust = 0.5), legend.position = "bottom")+
    geom_text(aes(x=xlabelpos,y=0.1,label = label_test_5ss),family = "sans",size=4,colour="#000000")+
    geom_text(aes(x=xlabelpos,y=0.2,label = label_test2_5ss),family = "sans",size=4,colour="#000000")+
    ggforce::facet_zoom(xlim = zoom.x, zoom.size = 1)
  #3ss
  message("[Step3] Plotting Cumulative Distribution Function for 3ss splicing sites.")
  efficiency_3ss_nona_inf_reduct<-subset(efficiency_5ss_3ss_nona_inf_reduct,grepl("3ss",efficiency_5ss_3ss_nona_inf_reduct$GeneID))
  #K-S test
  KSTEST_result_3ss<-stats::ks.test(efficiency_3ss_nona_inf_reduct[,2],efficiency_3ss_nona_inf_reduct[,3])
  D_value_3ss<- KSTEST_result_3ss[["statistic"]][["D"]]
  label_test2_3ss <- substr(paste0("D=",D_value_3ss),1,8)
  p_value_3ss<- KSTEST_result_3ss[["p.value"]]
  label_test_3ss <- ifelse(p_value_3ss<0.0001, "p<0.0001",substr(paste0("p=",p_value_3ss),1,8))

  #data_transform
  ecdf_exp_3ss_con_nona_w_to_L<- stats::reshape(data=efficiency_3ss_nona_inf_reduct, idvar="GeneID",
                                         varying = sample_list_show,
                                         v.name=c("Splicing_efficiency"),
                                         times=sample_list_show,
                                         #new.row.names = 1:1000,
                                         direction="long")
  ecdf_exp_3ss_con_nona_w_to_L$LOG2SE<- log2(ecdf_exp_3ss_con_nona_w_to_L$Splicing_efficiency)
  colnames(ecdf_exp_3ss_con_nona_w_to_L)<-c("GeneID","Samples","Splicing_efficiency","LOG2SE")

  #ecdf_plot
  xlabelpos<- max(ecdf_exp_3ss_con_nona_w_to_L$LOG2SE)-3.5
  cdf_p_3ss<- ggplot2::ggplot(ecdf_exp_3ss_con_nona_w_to_L,aes(x=LOG2SE, color=Samples))+
    stat_ecdf()+
    labs(title = "Cumulative Distribution Function")+
    xlab("log2 3ss Splicing_efficiency")+
    ylab("Cumulative fraction")+
    theme(plot.title=element_text(hjust = 0.5), legend.position = "bottom")+
    geom_text(aes(x=xlabelpos,y=0.1,label = label_test_3ss),family = "sans",size=4,colour="#000000")+
    geom_text(aes(x=xlabelpos,y=0.2,label = label_test2_3ss),family = "sans",size=4,colour="#000000")+
    ggforce::facet_zoom(xlim = zoom.x, zoom.size = 1)

  message("[Step4] Merge and export Cumulative Distribution Function plot.")
  grDevices::pdf("./Cumulative Distribution Function Plot.pdf", width = 10, height = 7.5)
  gridExtra::grid.arrange(cdf_p,cdf_p_5ss,cdf_p_3ss,ncol=3)
  grDevices::dev.off()
  message("Cumulative Distribution Function plot was successfully exported into the environment root path! ")

  KStestlist<-list(KSTEST_result,KSTEST_result_5ss,KSTEST_result_3ss)
  names(KStestlist)<-c("KSTEST_result_all","KSTEST_result_5ss","KSTEST_result_3ss")
  return(KStestlist)
}
