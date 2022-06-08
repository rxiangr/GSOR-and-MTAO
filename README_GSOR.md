# GSOR 
# Code for Genetic Score Omics Regression (GOSR)

GSOR is a method that associates assayed gene expression with genetically predicted phenotype, i.e., genetic score. This score is called estimated breeding value in animals or plants and polygenic score in humans. When the score is calculated using variants near the gene with assayed expression, GSOR provides a powerful test of the association between cis-effects on gene expression and the trait. The interpretation of the result is that an omic feature, such as gene expression or splicing is genetically linked to a phenotype via cis-regulatory mechanisms. This document provides a guide for running GSOR. 

First, we download Rscript and test datasets from https://figshare.com/s/ffed6891675508ecaf5b. 

Note that the R scripts are designed to run in a Unix environment. We expect the following R-packages are pre-installed: 

“data.table” (https://cran.r-project.org/web/packages/data.table/index.html)

“dplyr” (https://cran.r-project.org/web/packages/dplyr/index.html)

“OmicKriging” (https://cran.r-project.org/web/packages/OmicKriging/index.html)

“lqmm” (https://cran.r-project.org/web/packages/lqmm/index.html)

“gaston” (https://cran.r-project.org/web/packages/gaston/index.html)


To run a GSOR analysis, you can:

**cd GSORpub #(go to the folder with test data and scripts)**

**#(call R module on Unix) module load R/4.0.0-foss-2020a**

**Rscript GSOR_v1.R**  **--bv** bvtable/chr27.profile.tbl.gz **--exp** exp/Blood.voomgctaphe.txt **--genename** exp/Blood.voomgctaphe.header -**-cfe** fe/Blood.cfe.txt **--qfe** fe/Blood.qfe.txt  **--GRM** grm/Blood.grm  **--outpath** tmpout  **--outpref** Blood.tr16.chr27

**--bv** reads in a genetic score file:

zcat bvtable/chr27.profile.tbl.gz|head -3|cut -d' ' -f1-5

IID ENSBTAG00000000225 ENSBTAG00000000226 ENSBTAG00000000357 ENSBTAG00000000743

0032-blood 1.1692e-06 1.17831e-06 2.6372e-07 7.0213e-07

0288-blood 1.30493e-06 1.30759e-06 -2.71449e-07 6.5568e-07

This file has header and the first column is the individual ID (‘IID’), and the following columns are cis genetic score predicted using SNPs with ±1Mb to the gene or intron. In GSOR this file are used as the Ys (response variable). Later we also provide Rscripts showing how to prepare this file.

**--exp** reads in a file with expression profiles. We recommend that all gene expression levels be normalised before the analysis

head -2 exp/Blood.voomgctaphe.txt|cut -d' ' -f1-5

0032-blood 0032-blood 4.09232302301061 3.37891028546082 -4.83821214026649

0288-blood 0288-blood 3.79809040484976 3.46097003467383 -4.73589897719253

This file does not have a header. The first column is family ID and the 2nd column is the individual ID and from the 3rd column onwards are the normailsed gene expression values (or any measurements of omics data). 

**--genename** reads in a gene ID file. Because the expression profile does not have a header, we need to link the expression values with gene id and this input provides this link:

head -3 exp/Blood.voomgctaphe.header

ENSBTAG00000020035

ENSBTAG00000011528

ENSBTAG00000054497

This file should have one gene ID per line and the order of ID should be the same as the order from the 3rd column onwards in the file given to “--exp”. 

**--cfe** reads in a categorical fixed effect file. this file is optional but recommended.  

head -3 fe/Blood.cfe.txt

0032-blood 0032-blood Bos_taurus 1 0
0288-blood 0288-blood Bos_taurus 1 0
0299-blood 0299-blood Bos_taurus 1 0

This file has no header. Following the style of expression profiles, the first column is family ID and the 2nd column is the individual ID and from the 3rd column onwards are the fixed effects, e.g., breed types, locations and etc.

**--qfe** reads in a quantitative fixed effect file. This file is optional but recommended.

head -3 fe/Blood.qfe.txt

0032-blood 0032-blood 21 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

0288-blood 0288-blood 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

0299-blood 0299-blood 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

This file has no header. Following the style of the input for “--cfe", the first column is family ID and the 2nd column is the individual ID and from the 3rd column onwards are the numerical fixed effects, e.g., PCs, PEER factors and etc.

**--GRM** reads in a set of files for a genomic relationship matrix (GRM) in the formate of GCTA (https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM). Fitting a GRM is optional but highly recommended as it adjusts confounding factors to control spurious associations. 

ls grm/Blood.grm*

grm/Blood.grm.N.bin  grm/Blood.grm.bin  grm/Blood.grm.id

The GRM is computed using GCTA and stored in a set of 3 files: *.grm.N.bin, *.grm.bin and *.grm.id

**--outpath** specifies the output directory for results. If not specified, the results will be written into the current directory

**--outpref** specifies the output prefix. If not specified, the prefix will be ‘test’.

**If run successfully, the output from GSOR_v1.R will be generated as:**

zcat tmpout/Blood.tr16.chr27.gosr.lmm.gz|head -5

gene b se t p

ENSBTAG00000000225 -0.019 0.035 -0.543 0.587129802433263

ENSBTAG00000000226 -0.008 0.02 -0.4 0.689156516779352

ENSBTAG00000000357 0.006 0.049 0.122 0.902899018270974

ENSBTAG00000000744 -0.007 0.034 -0.206 0.836790911305215

4 columns of results are generated, gene ID, the association beta, standard error, t-value (b/se) and p-value.

In order to perform GSOR, we need a genetic score file, which is an EBV or PGS predicted using genotypes of cis variants and a vector of SNP effects jointly trained (e.g., gBLUP). We provide an R script that inputs pre-estimated SNP BLUP effects and plink genotype file to generate the genetic score file:

**###################prepare genetic score files##############
######step 1########**

**cd GSORpub**

**#(call R module on Unix) module load R/4.0.0-foss-2020a**

**Rscript cisgc.R** **--geneanno** geneanno/cattle.gene.anno.txt **--chrn** 27 **--target-plinkpath** plinkGeno **--target-plink-pref** chr27 **--snpv** snpv/test.snp.blp.gz  **--snplist-outpath** snplist **--gs-outpath** bvout **--genedis** 1000000 **--minNsnp** 2 **--genetype** 1 **--plinkpath** plink/plink **--memmb** 10000

**--geneanno** reads in a genome annotation file (tab-separated).

head -3 geneanno/cattle.gene.anno.txt

Gene stable ID  Gene name       Chromosome/scaffold name        Gene start (bp) Gene end (bp)       Strand

ENSBTAG00000006648              1       339070  350389  -1

ENSBTAG00000049697      5S_rRNA 1       475398  475516  1

**--chrn** specifies the chromosome number to work on. Its optional but highly recommended to separate the analysis into different chromosomes to increase computation efficiency.

**--target-plinkpath** specifies the full path to a set of plink genotype files

ls plinkGeno

chr27.bed  chr27.bim  chr27.fam

**--target-plink-pref** specifics the prefix of plink files:

ls plinkGeno/chr27*

plinkGeno/chr27.bed  plinkGeno/chr27.bim  plinkGeno/chr27.fam

**--snpv** reads in a jointly trained SNP effects file. here we use the results and the format from GCTA SNP BLUP (https://yanglab.westlake.edu.cn/software/gcta/#BLUP).

zcat snpv/test.snp.blp.gz|head -3

Chr27:20985     A       -5.61085e-06    3.72716e-16

Chr27:29600     C       -2.89424e-06    -2.01698e-16

Chr27:29601     TAG     -2.68128e-06    -2.47525e-16

The columns are SNP ID, reference allele and BLUP of SNP effect. The last column is for the residual effect. Note that when doing SNP BLUP, its recommended that all genome-wide SNPs are trained together even if you are only using SNP effects from one chromosome. 

**--snplist-outpath** specifies an intermediate output path to store SNP lists used to make cis-prediction of genetic score

ls snplist|head -5

ENSBTAG00000000225.snplist

ENSBTAG00000000226.snplist

ENSBTAG00000000357.snplist

ENSBTAG00000000743.snplist

ENSBTAG00000000744.snplist

These files indicate the SNPs used to predict the genetic score around each gene.

**--gs-outpath** specifies the output from cisgc.R as genetic score for using SNPs around each gene.

zcat bvout/ENSBTAG00000000225.profile.gz|head -3

FID                       IID  PHENO    CNT   CNT2    SCORE

0032-blood                0032-blood     -9  15400  10823 1.1692e-06
 
 0288-blood                0288-blood     -9  15400  11576 1.30493e-06

This is an example of the genetic score for one gene (ENSBTAG00000000225) which is generated using plink --score (https://www.cog-genomics.org/plink/1.9/score). The 1st and 2nd columns are family ID and an individual ID; ignoring columns 3-5 and the last column “SCORE” is the predicted genetic score for this gene across all individuals.

**--genedis** specifies the distance to a gene in which SNPs are used to predict the genetic score. The default is 1000000.

**--minNsnp** specifies the minimal number of SNPs used for the prediction of genetic score. If fewer than 2 or a specified value, the genetic score will not be predicted.

**--genetype** specifies if gene or intron is being analysed. The default is 1. If a gene is analysed then the input for “--genedis” will be based on the starting (TSS) of the gene, whereas if an intron is analysed, “--genedis” will be based the intron as the center.

**--plinkpath** specifies the full path to plink software
#try plink/plink to check it.

**--memmb** specifies the memory in Mb will be used. 10000 is the default. 

Note that this step will generate many results in --snplist-outpath and --gs-outpath. Its recommended to only store these results somewhere temporary. In the following text, we provide another Rscript to combine genetic scores across all genes which will be used in the GSOR_v1.R. 

**###################prepare genetic score files##############**
**######step 2########**

**cd GSORpub**

**#(call R module on Unix) module load R/4.0.0-foss-2020a**

**Rscript cbgc.R** **--suffx** profile.gz **--inpath** bvout **--mergfiles** n **--outpath** bvtable **--outpref** chr27

**--suffx** specifies a suffix pattern in a folder in which the software will bind results together

**--inpath** specifies a path where all files that need the combining is

ls bvout/*profile.gz|head -3

bvout/ENSBTAG00000000225.profile.gz

bvout/ENSBTAG00000000226.profile.gz

bvout/ENSBTAG00000000357.profile.gz

These are the files generated using cisgc.R described before (see “--gs-outpath bvout”)

**--mergfiles** specifies if you want to merge (‘y’) the files or simply append (‘n’) the files. Simple appending assumes that all results in bvout/*profile.gz have the same individuals in the same order, which is usually the case. The default is “n”. if you choose “y”, the software will merge results based on the column “IID” but this will take a long time to finish.

**--outpath** specifies the output path.

**--outpref** specifies the output prefix.

zcat bvtable/chr27.profile.tbl.gz|head -5|cut -d' ' -f1-4

IID ENSBTAG00000000225 ENSBTAG00000000226 ENSBTAG00000000357

0032-blood 1.1692e-06 1.17831e-06 2.6372e-07

0288-blood 1.30493e-06 1.30759e-06 -2.71449e-07

0299-blood 1.30481e-06 1.30748e-06 1.00378e-06

0308-blood 1.36324e-06 1.36521e-06 -9.05266e-07

This output file from cbgc.R and has a header. The 1st column is the individual ID and from 2nd column onwards are the predicted genetic score using SNPs around each gene listed. This file will be used as the phenotype file in GSOR_v1.R as described above. See  “**--bv** bvtable/chr27.profile.tbl.gz”.
