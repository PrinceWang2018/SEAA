#' Splicing Efficiency Calculation
#' @description This function is used to calculate the splicing efficiency.
#' @param BAM_files Aligned reads stored in indexed bam files.(2 bam files)
#' @param Anno_SAF Splicing sites annotation in a saf file.
#' @param paired Whether the sequencing is paired or not.
#' @param thread threads number used in calculating
#' @param strand which strand of the sequencing
#'
#' @return a large list
#' @export
#'
#' @examples
#' \donttest{
#' BAM_files <- c("~/project3tB/SEAA_project/R_counstruction_place/SEAA/inst/extdata/NC_t.bam",
#' "~/project3tB/SEAA_project/R_counstruction_place/SEAA/inst/extdata/USP_t.bam")
#' Anno_SAF<-"~/project3tB/SEAA_project/R_counstruction_place/SEAA/inst/extdata/test.saf"
#' SEresultlist<-SEcalculation(BAM_files,Anno_SAF,paired = TRUE ,thread = 2,strand = 1)
#' }
#'
SEcalculation <-function(BAM_files,Anno_SAF,paired = TRUE ,thread = 1,strand = 1){
  message("[Step1] Counting unspliced reads with featureCounts.")
  fc_result_unspliced <- Rsubread::featureCounts(BAM_files,annot.inbuilt = NULL, annot.ext = Anno_SAF,
                                       isPairedEnd = paired, nthreads = thread, strandSpecific = strand, countMultiMappingReads = FALSE, allowMultiOverlap = TRUE,
                                       minOverlap = 2)
  fc_count_unspliced<- as.data.frame(fc_result_unspliced[["counts"]])
  message("[Step2] Counting all reads aligned to the splicing sites.")
  fc_result_total <- Rsubread::featureCounts(BAM_files,annot.inbuilt = NULL, annot.ext = Anno_SAF,
                                   isPairedEnd = paired, nthreads = thread, strandSpecific = strand, countMultiMappingReads = FALSE, allowMultiOverlap = TRUE,
                                   minOverlap = 1)
  fc_count_total<- as.data.frame(fc_result_total[["counts"]])

  fc_stat_total<- as.data.frame(fc_result_total[["stat"]])
  if (fc_stat_total[1,2] / sum(fc_stat_total[,2]) < 0.05) {
    message("Warning: Few reads are mapped to the splice sites. Please try another strand-specific parameter.")
  }
  fc_anno<-as.data.frame(fc_result_unspliced[["annotation"]])

  sample_list <- colnames(fc_count_unspliced)
  sample_list_show <-gsub(".bam","",sample_list)
  message("[Step3] Calculating spliced reads.")
  fc_count_spliced<- fc_count_total[, sample_list] - fc_count_unspliced[, sample_list]
  message("[Step4] Calculating splicing efficiency and exporting.")
  efficiency_5ss_3ss <- cbind(fc_anno, fc_count_spliced[, sample_list] / fc_count_unspliced[, sample_list])
  #utils::write.table(efficiency_5ss_3ss,file = "./splicing_efficiency_5ss_3ss_fc_raw.txt", sep = "\t", row.names = FALSE)
  SEresultlist<-list(efficiency_5ss_3ss,fc_count_spliced,fc_count_total,fc_count_unspliced,fc_anno,sample_list,sample_list_show)
  names(SEresultlist)<-c("efficiency_5ss_3ss","fc_count_spliced","fc_count_total","fc_count_unspliced","fc_anno","sample_list","sample_list_show")
  return(SEresultlist)
}


