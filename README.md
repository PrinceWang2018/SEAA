
# SEAA(Splicing Efficiency Analysis and Annotation)

<!-- badges: start -->
<!-- badges: end -->

SEAA provides convenient and rapid splicing efficiency calculation and splicing 
sites annotation function using next generation sequencing data. Aligned .bam 
files and a processed splicing sites .saf file are needed. The task can be finished 
in several minutes  multi-core computing with the assitance of 'Rsubread'. Plots 
of splicing status type and Cumulative Distribution Function (CDF) and annotated 
vaild splicing efficiency can be exported. 

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

## Data preparation

Aligning your sequencing files with hisat2 is strongly recommended. Genome_trans
index of hisat2 can be downloaded from http://daehwankimlab.github.io/hisat2/download/.
It should be noticed that ensembl version reference files with chromosome number 
"1" instead of "chr1" were used. Or you will need to edit the saf file manually.

``` r
hisat2 -x /indexpath/hisat2/grch37_tran/genome_tran -1 /fastqpath/NC_1.fastq -2 /fastqpath/NC_2.fastq --min-intronlen 20 --max-intronlen 10000 --threads 12 --rna-strandness F | samtools sort -o /outputpath/NC.bam - 
```
We have already prepared the saf files consisting the splicing sites of human (hg19)
and mouse (mm10) which can be downloaded from our github repository https://github.com/PrinceWang2018/SEAA_reference.
 
## Example

This is a basic workflow which shows you how to use this software:

``` r
library("SEAA")
#set you work path
setwd("/home/username/workpath/")
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
Step5: Annotation of Splicing Sites Acquiring filtered Splicing Efficiency.
``` r
SEannotaionresult<-SEsiteanno(SEresultlist, efficiency_5ss_3ss_nona_inf_reduct, species = "hs")
```
Step6: Labeling the target gene in filtered splicing efficiency list.
``` r
targetlabeling(SEresultlist,target_site = "27830321",target_label = "RPL21",xlim.max = 1000, ylim.max = 1000)
```

## Citation
Wang et al., (2021). SEAA: Splicing Efficiency Analysis and Annotation (to be published)
## Contact
If you have any question, please contact Zixiang Wang at wangzixiang@mail.sdu.edu.cn or wangzixiang@live.com.
