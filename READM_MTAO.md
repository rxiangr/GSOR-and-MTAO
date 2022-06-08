# Running multi-trait meta-analyses of omics-associations (MTAO)

MTAO is a summary-data-based method that combines effects on complex traits of omics features such as genes or introns. MTAO will estimate two metrics from the meta-analysis: one is Pn which is the number of traits the omics feature is associated with, and the other is Pm which is the magnitude of multi-trait effects of each omics feature. If Pn is significant, then this omics feature is associated with at least 2 traits. If Pm is significant, this omics feature is affecting >= 1 trait. Note that MTAO can deal with summary data as long as there is beta and se estimated for each omics feature, including from GSOR or the conventional TWAS.
In the following, we provide Rscript and test data to perform MTAO.

First, we download Rscript and test datasets from https://figshare.com/s/ffed6891675508ecaf5b. 

Note that the R scripts are designed to run in a Unix environment. We expect the following R-packages are pre-installed: 
“data.table” (https://cran.r-project.org/web/packages/data.table/index.html)

“dplyr” (https://cran.r-project.org/web/packages/dplyr/index.html)

To run a MTAO analysis, you can:

**cd GSORpub #(go to the directory with test datasets)**

**#(call R module on Unix) module load R/4.0.0-foss-2020a**

Rscript MTAO_v1.R **--inpath** ExamStRes **--pattern** gsor.lmm.gz **--geneanno** geneanno/cattle.gene.anno.txt -**-outpath** MTAOtmpres **--outpref** test

**--inpath** specifies a path to a folder where single-trait GSOR or TWAS results are stored.

**--pattern** specifies a pattern of files that MTAO will process. In this case, this pattern is “gsor.lmm.gz”.

In example data there are:

ls ExamStRes/*gsor.lmm.gz

ExamStRes/cow.Blood.gene.tr04.chr27.gsor.lmm.gz

ExamStRes/cow.Blood.gene.tr16.chr27.gsor.lmm.gz

ExamStRes/cow.Blood.gene.tr32.chr27.gsor.lmm.gz

These are the per trait results of the association analysis between gene expression and traits. for each result,

zcat ExamStRes/cow.Blood.gene.tr04.chr27.gsor.lmm.gz|head -3

gene b se t p

ENSBTAG00000000225 0.022 0.034 0.647 0.517631943475741

ENSBTAG00000000226 -0.011 0.02 -0.55 0.582319373576693

This file has a header. here are 4 columns, the 1st column is the gene ID, followed by columns of beta, se, t value (beta/se) and p-value. This is the required format of the association analysis results for MTAO to process. 

**--geneanno** specifies an annotation file for genes. 

head -3 geneanno/cattle.gene.anno.txt

Gene stable ID  Gene name       Chromosome/scaffold name        Gene start (bp) Gene end (bp)       Strand

ENSBTAG00000006648              1       339070  350389  -1

ENSBTAG00000049697      5S_rRNA 1       475398  475516  1 


**--outpath** specifies the path of output

**--outpref** specifies the prefix of output

ls MTAOtmpres/

test.gsor.sum.gz  test.mtao.pleio.gz


**MTAO will produce two output files.** 

zcat MTAOtmpres/test.gsor.sum.gz|head -3

gene gene.name chromosome genest geneend strand tr04.b tr04.se tr04.t tr04.p tr16.b tr16.se tr16.t tr16.p tr32.b tr32.se tr32.t tr32.p

ENSBTAG00000017096 ERICH1 27 498301 531821 -1 -0.071 0.066 -1.1 0.271332121892765 -0.11 0.067 -1.6 0.109598583399116 0.006 0.066 0.091 0.927492591057545

ENSBTAG00000051728 unknown 27 1127690 1182591 1 -0.01 0.023 -0.435 0.663562427438038 -0.021 0.026 -0.808 0.419090582142142 0.006 0.026 0.231 0.817314802029568

The 1-6 columns are annotations of the gene; from the 7th column onwards there are beta, standard error, t-value (beta/se) and p-value of the association between the expression of gene and traits (genetically predicted phenotype) analysed; beta, se, t-value and p are summarised for all analysed traits.

zcat MTAOtmpres/test.mtao.pleio.gz|head -3

gene gene.name chromosome genest geneend strand Pn Pn_Pvalue Pn_Pvalue.adj x X2.p X2.p.adj

ENSBTAG00000017096 ERICH1 27 498301 531821 -1 0 1 1 1.75435564128869 0.379790249957691 0.997085261500835

ENSBTAG00000051728 unknown 27 1127690 1182591 1 0 1 1 0.886814165030834 0.852707733431646 0.997085261500835

This is result of meta-analysis across all traits. The 1-6 columns are annotations of the gene; then the 7-12 columns are: number of traits this gene associated with (Pn); p-value of Pn; FDR-adjusted p-value of Pn; square root of Chi-square of the multi-trait analysis, Chi-square p-value (Pm); and FDR-adjusted Chi-square p-value. 

