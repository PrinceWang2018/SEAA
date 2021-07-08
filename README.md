
# SEAA(Splicing Efficiency Analysis and Annotation)

<!-- badges: start -->
<!-- badges: end -->

SEAA provides convenient and rapid splicing efficiency
    calculation and splicing sites annotation function using next generation 
    sequencing data. Aligned .bam files and a processed splicing sites .saf file
    are needed. The task can be finished in several minutes using multicores with
    the assitance of 'Rsubread'. Plots of splicing status type and Cumulative 
    Distribution Function (CDF) and annotated vaild splicing efficiency can be 
    exported. 

## Installation

You can install the released version of SEAA from [github](https://github.com/PrinceWang2018/SEAA) with:

``` r
install.packages("devtools")
library(devtools)
install_github("PrinceWang2018/SEAA")
```
You can also install the released version of SEAA from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("SEAA")
```
in the near furture.

## Example

This is a basic workflow which shows you how to use this software:

``` r
library("SEAA")
#set you work path
setwd("/home/wzx/project3tB/SEAA_project/pkg_test/")
```
Step1: Caculating splicing efficiency.
``` r
#Please locate your aligned bam files.
BAM_files <- c("./NC_1.bam","./shUSP39_1.bam")
#Please locate your .saf annotation files.
Anno_SAF<- "/home/wzx/project3tB/SEAA_project/reference/hg19/hg19_NCBI_splicing_sites_20210705.saf"
SEresultlist<-SEcalculation(BAM_files,Anno_SAF,paired = TRUE ,thread = 8,strand = 1)
```
Step2: Caculating splicing efficiency with featureCounts.
``` r
SEtyperesult<-SEtypeplot(SEresultlist,"horizontal")
```
Step3: Filtering invaild splicing efficiency and reads counts less than min_counts.
``` r
efficiency_5ss_3ss_nona_inf_reduct<-SEfilter(SEresultlist,min_counts = 5)
```
Step4: Cumulative Distribution Function plotting based on filtered splicing efficiency.
``` r
CDFplot(SEresultlist,efficiency_5ss_3ss_nona_inf_reduct,zoom.x= c(5,6))
```
Step5: labeling the target gene in filtered splicing efficiency list.
``` r
target_infor<-targetlabeling(SEresultlist,target_site = "48753008",target_label = "CARD8",xlim.max = 100, ylim.max = 100)
```
Step6: Annotation of Splicing Sites Acquiring filtered Splicing Efficiency.
``` r
SEannotaionresult<-SEsiteanno(efficiency_5ss_3ss_nona_inf_reduct)
```

## Citation

## Contact
If you have any question, please contact Zixiang Wang at wangzixiang@mail.sdu.edu.cn .
