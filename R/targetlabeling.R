#' Labeling the Target Gene in the Splicing Efficiency List
#'
#' @param SEresultlist a large list created by SEcalculation
#' @param target_site a character indicating target splicing site position
#' @param target_label a character indicating the label you want to add
#' @param xlim.max a numeric value indicating the max x axis
#' @param ylim.max a numeric value indicating the max y axis
#' @importFrom ggplot2 ggplot xlim ylim theme xlab ylab labs element_text aes geom_point unit
#' @return a data frame of the target information
#' @export
#'
#' @examples
#' targetlabeling(SEresultlist,target_site = "27830321",target_label = "RPL21",
#' xlim.max = 1000, ylim.max = 1000)

targetlabeling<-function(SEresultlist,target_site = "48753008",target_label = "CARD8", xlim.max = 100, ylim.max = 100){
  efficiency_5ss_3ss<-SEresultlist[["efficiency_5ss_3ss"]]
  sample_list_show<-SEresultlist[["sample_list_show"]]

  target <- grep(target_site,efficiency_5ss_3ss$End)
  if (length(target > 1)) {
    message("Note: Multiple sites are found in the splicing efficiency list and only one is used for plotting.")
    target_choose <- target[1]
    df2 <- efficiency_5ss_3ss[target_choose,]
  } else {
    df2 <- efficiency_5ss_3ss[target,]
  }
  xlab_text<- paste("Splicing efficiency in ", sample_list_show[1])
  ylab_text<- paste("Splicing efficiency in ", sample_list_show[2])

  label_plot<-ggplot(efficiency_5ss_3ss, aes(x=efficiency_5ss_3ss[,7],y=efficiency_5ss_3ss[,8]))+geom_point(colour = "#92C5DE", shape=".",alpha= 0.1)+xlim(NA, xlim.max)+ylim(NA, ylim.max)+
  ggrepel::geom_text_repel(data=df2, aes(x=df2[,7],y=df2[,8]), label=target_label,fontface="bold", color="black",
                                       box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"),
                                       segment.colour = "#000000",min.segment.length=0.001)+
  labs(title = "Splicing efficiency focusing sites")+
  xlab(xlab_text)+
  ylab(ylab_text)+
  theme(plot.title = element_text(hjust = 0.5))

  if (df2[,7]>xlim.max) {
    message("Warning: xlim.max was too low and target site was out of range.")
  }
  if (df2[,8]>ylim.max) {
    message("Warning: ylim.max was too low and target site was out of range.")
  }
  return(label_plot)
}
