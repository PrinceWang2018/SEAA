
# SEAA(Splicing Efficiency Analysis and Annotation)

<!-- badges: start -->
<!-- badges: end -->

SEAA provides convenient and rapid splicing efficiency calculation and splicing 
sites annotation function using next generation sequencing data. Aligned .bam 
files and a processed splicing sites .saf file are needed. The task can be finished 
in several minutes  multi-core computing with the assistance of 'Rsubread'. Plots 
of splicing status type and Cumulative Distribution Function (CDF) and annotated 
vaild splicing efficiency can be exported. 

## Installation

You can install the released version of SEAA from [GitHub](https://github.com/PrinceWang2018/SEAA) with:

``` r
install.packages("devtools")
library(devtools)
install_github("PrinceWang2018/SEAA")
```
## Data preparation

Aligning your sequencing files with hisat2 is strongly recommended. Genome_trans
index of hisat2 can be downloaded from http://daehwankimlab.github.io/hisat2/download/.
It should be noticed that Ensembl version reference files with chromosome number 
"1" instead of "chr1" were used. Or you will need to edit the saf file manually.

``` r
hisat2 -x /indexpath/hisat2/grch37_tran/genome_tran -1 /fastqpath/NC_1.fastq -2 /fastqpath/NC_2.fastq --min-intronlen 20 --max-intronlen 10000 --threads 12 --rna-strandness F | samtools sort -o /outputpath/NC.bam - 
```
We have already prepared the saf files consisting the splicing sites of human (hg19)
and mouse (mm10) which can be downloaded from our GitHub repository https://github.com/PrinceWang2018/SEAA_reference.

## Example

This is a basic workflow which shows you how to use this software:

``` r
library("SEAA")
#set you work path
setwd("/home/username/workpath/")
```
**Step 1 Calculating splicing efficiency.**

``` r
#Please locate your aligned bam files.
BAM_files <- c("./NC_1.bam","./shUSP39_1.bam")
#Please locate your .saf annotation files.
Anno_SAF <- "/home/wzx/project3tB/SEAA_project/reference/hg19/hg19_NCBI_splicing_sites_20210705.saf"
SEresultlist <- SEcalculation(BAM_files,Anno_SAF,paired = TRUE ,thread = 8,strand = 1)
```
If you come across the following warning, that means the successfully matched splicing sites are few:

> Warning: Few reads are mapped to the splice sites. Please try another strand-specific parameter.

Please change the "strand = " parameter among 0,1,2 until you acquire a desirable splicing status ratio.

**Step 2 Splicing efficiency status classification and visualization.**

``` r
SEtyperesult<-SEtypeplot(SEresultlist,"horizontal")
```
**Step 3 Filtering invalid splicing efficiency and reads counts less than min_counts.**

``` r
efficiency_5ss_3ss_nona_inf_reduct<-SEfilter(SEresultlist,min_counts = 5)
```
**Step 4 Cumulative Distribution Function plotting based on filtered splicing efficiency.**

``` r
CDFplot(SEresultlist,efficiency_5ss_3ss_nona_inf_reduct,zoom.x= c(5,6))
```
**Step 5 Annotation of Splicing Sites Acquiring filtered Splicing Efficiency.**

``` r
SEannotaionresult<-SEsiteanno(SEresultlist, efficiency_5ss_3ss_nona_inf_reduct, species = "hs")
```
**Step 6 Labeling the target gene in filtered splicing efficiency list.**

``` r
targetlabeling(SEresultlist,target_site = "27830321",target_label = "RPL21",xlim.max = 1000, ylim.max = 1000)
```

## Citation
If you think this software interesting and use it in your research, you are welcome to use the text and paper below to cite our program:

> Splicing efficiency analysis in this research was conducted with an R package Splicing Efficiency Analysis and Annotation (SEAA) (https://github.com/PrinceWang2018/SEAA) [1].
>
> [1] Qin, J., Huang, T., Wang, Z. *et al.* Bud31-mediated alternative splicing is required for spermatogonial stem cell self-renewal and differentiation. *Cell Death Differ* (2022). https://doi.org/10.1038/s41418-022-01057-1



## Contact
If you have any question, please contact Zixiang Wang at wangzixiang@sdu.edu.cn or wangzixiang@live.com.
