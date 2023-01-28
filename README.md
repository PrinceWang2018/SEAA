
# SEAA (Splicing Efficiency Analysis and Annotation)

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/384065256.svg)](https://zenodo.org/badge/latestdoi/384065256)
<!-- badges: end -->

![Workflow](https://user-images.githubusercontent.com/37299182/215256146-2cc29c9e-ded1-44e3-b6bd-0d1cd5731f81.png)

SEAA provides convenient and rapid splicing efficiency calculation and splicing 
sites annotation function using next generation sequencing data. Aligned .bam 
files and a processed splicing sites .saf file are needed. The task can be finished 
in several minutes  multi-core computing with the assistance of 'Rsubread'. Plots 
of splicing status type and Cumulative Distribution Function (CDF) and annotated 
vaild splicing efficiency can be exported. 

## Citation
If you think this software interesting and use it in your research, you are welcome to use the text and paper below to cite our program:

> Splicing efficiency analysis in this research was conducted with an R package Splicing Efficiency Analysis and Annotation (SEAA) (https://github.com/PrinceWang2018/SEAA) [1].
>
> [1] Qin, J., Huang, T., Wang, Z. *et al.* Bud31-mediated alternative splicing is required for spermatogonial stem cell self-renewal and differentiation. *Cell Death Differ* (2022). https://doi.org/10.1038/s41418-022-01057-1


## Contact
If you have any question, please contact Zixiang Wang at wangzixiang@sdu.edu.cn or wangzixiang@live.com.


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
We have already prepared the saf files consisting the splicing sites of human (hg19, hg38)
and mouse (mm10) which can be downloaded from our GitHub repository https://github.com/PrinceWang2018/SEAA_reference.

## Workflow

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
You can edit the .pdf file with Acrobat manually.

**Step 5 Annotation of Splicing Sites Acquiring filtered Splicing Efficiency.**

``` r
SEannotaionresult<-SEsiteanno(SEresultlist, efficiency_5ss_3ss_nona_inf_reduct, species = "hs")
```
**Step 6 Labeling the target gene in filtered splicing efficiency list.**

``` r
targetlabeling(SEresultlist,target_site = "27830321",target_label = "RPL21",xlim.max = 1000, ylim.max = 1000)
```

## Examples using SEAA and experimental verification

This is an example using the SEAA package to analyze splicing efficiency. The usage of SEAA is illustrated with several subsets of transcriptome data from the GEO database including knocking down SF3B1, PRPF8, and USP39 and cells treated with Pladienolide B, a U2 spliceosome inhibitor. 

![Figure2](https://user-images.githubusercontent.com/37299182/215256518-4dfd01df-a728-4691-99ab-c1368d4bf550.png)

**Figure 1** Summarization of the splicing status in each sample was accomplished by SEtypeplot. A basic abstract graph of the splicing status after splicing efficiency calculation was shown in Figure 1A. All the splicing sites were classified into four types, including no reads, unspliced, spliced, and partly spliced. Both splicing sites without signal (no reads) and without spliced reads (unspliced) were considered as splicing events undetected. Due to the poor stability of the unspliced RNA and poly (A) RNA-seq preference, totally spliced splicing sites (Spliced) usually account for the majority. The remaining splicing sites acquired unspliced and spliced reads were labeled as PartlySpliced. Knocking down the core splicing factors or treatment with splicing inhibitor Pladienolide B leads to a decreased ratio of Spliced genes. Especially, SF3B1 ablation caused a more reduction in totally spliced genes (Figure 1B-E). 

![Figure3](https://user-images.githubusercontent.com/37299182/215256513-78ebcca4-6085-4455-8689-97aaae85bdbf.png)

**Figure 2** The precise splicing efficiency change was further visualized with CDF plot. D-statistic and p-value generated by Kolmogorov-Smirnov (K-S) test indicated the splicing efficiency difference after Pladienolide treatment or knocking down PRPF8 and USP39.

![Figure4](https://user-images.githubusercontent.com/37299182/215256509-247b476e-35db-499b-b41b-4d9eb7b0d958.png)

**Figure 3** The precise splicing efficiency change was further visualized with CDF plot. D-statistic and p-value generated by Kolmogorov-Smirnov (K-S) test indicated the splicing efficiency difference after knocking down SF3B1.

![Figure6](https://user-images.githubusercontent.com/37299182/215256503-d53f2950-f35f-46fd-857a-edd7ff27b54e.png)

**Figure 4** The splicing efficiency of TP53, CDK2, APOE, and BRD9 decreased after splicing factor SF3B1 knockdown. The RNA-seq data was visualized with IGV and SEAA. Intron-retained specific primers were designed to detect splicing efficiency after siUSP39.

## Papers using the R package
> 1.	Qin J, Huang T, Wang Z, Zhang X, Wang J, Dang Q, et al. Bud31-mediated alternative splicing is required for spermatogonial stem cell self-renewal and differentiation. Cell Death Differ. 2022.
> 2.	Wang S, Wang Z, Li J, Qin J, Song J, Li Y, et al. Splicing factor USP39 promotes ovarian cancer malignancy through maintaining efficient splicing of oncogenic HMGA2. Cell Death & Disease. 2021;12(4):294.
> 3.	Wang X, Han M, Wang S, Sun Y, Zhao W, Xue Z, et al. Targeting the splicing factor NONO inhibits GBM progression through GPX1 intron retention. Theranostics. 2022;12(12):5451-69.

