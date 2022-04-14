#' Splicing Status Type Plotting
#'
#' @param SEresultlist a large list created by SEcalculation
#' @param arrangement a characteristic value
#' @return a large list
#' @export
#'
#' @examples
#' SEtyperesult<-SEtypeplot(SEresultlist,"horizontal")
#'
SEtypeplot<-function(SEresultlist,arrangement){
  efficiency_5ss_3ss<-SEresultlist[["efficiency_5ss_3ss"]]
  sample_list<-SEresultlist[["sample_list"]]
  sample_list_show<-SEresultlist[["sample_list_show"]]
  message("[Step1] Plotting splicing status type for all splicing sites.")
  se_type_l <- ""
  se_type <- as.data.frame(c("","",""))
  for (i in 1:length(sample_list)) {
    i
    num1 = length(grep('NaN',efficiency_5ss_3ss[,6+i]))
    num2 = length(grep('0',efficiency_5ss_3ss[,6+i]))
    num3 = num1+num2
    num4 = length(grep("Inf",efficiency_5ss_3ss[,6+i]))
    num5 = nrow(efficiency_5ss_3ss)-num3-num4
    se_type[,i] <- as.data.frame(c(num3,num4,num5))
  }
  colnames(se_type) = sample_list_show
  se_type[,ncol(se_type)+1] = c("Undetermined","Spliced","Partly_Spliced")
  se_type_l = reshape2::melt(se_type)
  colnames(se_type_l) = c('type','sample','num')

  p_all<-ggpubr::ggbarplot(se_type_l, "sample", "num",
                   fill = "type", color = "type", palette = "jco") +
    labs(x="Samples", y = " Number of Splicing Sites", title = "Splicing Status of all Splicing Sites") +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

  message("[Step2] Plotting splicing status type for 5ss splicing sites.")
  efficiency_5ss<-subset(efficiency_5ss_3ss,grepl("5ss",efficiency_5ss_3ss$GeneID))
  se_type_5ss <- as.data.frame(c("","",""))
  for (i in 1:length(sample_list)) {
    num1 = length(grep('NaN',efficiency_5ss[,6+i]))
    num2 = length(grep('0',efficiency_5ss[,6+i]))
    num3 = num1+num2
    num4 = length(grep("Inf",efficiency_5ss[,6+i]))
    num5 = nrow(efficiency_5ss)-num3-num4
    se_type_5ss[,i] <- as.data.frame(c(num3,num4,num5))
  }
  colnames(se_type_5ss) = sample_list_show
  se_type_5ss[,ncol(se_type_5ss)+1] = c("Undetermined","Spliced","Partly_Spliced")
  se_type_5ss_l = reshape2::melt(se_type_5ss)
  colnames(se_type_5ss_l) = c('type','sample','num')
  p_5ss<-ggpubr::ggbarplot(se_type_5ss_l, "sample", "num",
                   fill = "type", color = "type", palette = "jco" ) +
    labs(x="Samples", y = " Number of Splicing Sites", title = "Splicing Status of 5ss Splicing Sites") +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

  message("[Step3] Plotting splicing status type for 3ss splicing sites.")
  efficiency_3ss<-subset(efficiency_5ss_3ss,grepl("3ss",efficiency_5ss_3ss$GeneID))
  se_type_3ss <- as.data.frame(c("","",""))
  for (i in 1:length(sample_list)) {
    num1 = length(grep('NaN',efficiency_3ss[,6+i]))
    num2 = length(grep('0',efficiency_3ss[,6+i]))
    num3 = num1+num2
    num4 = length(grep("Inf",efficiency_3ss[,6+i]))
    num5 = nrow(efficiency_3ss)-num3-num4
    se_type_3ss[,i] <- as.data.frame(c(num3,num4,num5))
  }
  colnames(se_type_3ss) = sample_list_show
  se_type_3ss[,ncol(se_type_3ss)+1] = c("Undetermined","Spliced","Partly_Spliced")
  se_type_3ss_l = reshape2::melt(se_type_3ss)
  colnames(se_type_3ss_l) = c('type','sample','num')

  p_3ss<-ggpubr::ggbarplot(se_type_3ss_l, "sample", "num",
                   fill = "type", color = "type", palette = "jco" )+
    labs(x="Samples", y = " Number of Splicing Sites", title = "Splicing Status of 3ss Splicing Sites") +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

  message("[Step4] Merge and export the plot.")
  grDevices::pdf("./Splicing Status of the Splicing Sites.pdf", width = 12, height = 3)
  if (arrangement=="horizontal") {
    gridExtra::grid.arrange(p_all,p_5ss,p_3ss,ncol=3)
  } else if (arrangement=="vertical") {
    gridExtra::grid.arrange(p_all,p_5ss,p_3ss,nrow=3)
  } else{
    message("Error: Wrong arrangement parameter was imported!")
  }
  grDevices::dev.off()
  message("Splicing status type plot was successfully exported. ")
  SEtyperesult<-list(se_type,se_type_5ss,se_type_3ss)
  names(SEtyperesult)<-c("se_type","se_type_5ss","se_type_3ss")
  return(SEtyperesult)
}
