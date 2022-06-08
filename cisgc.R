#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<5) {
  stop("5 argument must be supplied", call.=FALSE)
}

suppressMessages(library("optparse"))
suppressMessages(suppressWarnings(library(data.table)))


#---parsing options
option_list = list(
make_option("--geneanno", action="store", default=NA, type='character'),
make_option("--chrn", action="store", default=NA, type='numeric'),
make_option("--target-plinkpath", action="store", default=NA, type='character'),
make_option("--target-plink-pref", action="store", default=NA, type='character'),
make_option("--snpv", action="store", default=NULL, type='character'),
make_option("--snplist-outpath", action="store", default='.', type='character'),
make_option("--gs-outpath", action="store", default='.', type='character'),
make_option("--genedis", action="store", default=1000000, type='numeric'),
make_option("--minNsnp", action="store", default=2, type='numeric'),
make_option("--genetype", action="store", default=1, type='integer'),
make_option("--plinkpath", action="store", default='.', type='character'),
make_option("--memmb", action="store", default=10000, type='numeric')
)
opt = parse_args(OptionParser(option_list=option_list))

geneannofn <- opt$geneanno
chrn <- opt$chrn
targplkpath <- opt$"target-plinkpath"
targplkpref <- opt$"target-plink-pref"
snpvfn <- opt$snpv
snplistoutpath <- opt$"snplist-outpath"
profileoutpath <- opt$"gs-outpath"
genedis <- opt$genedis
minNsnp <- opt$minNsnp
phetype <- opt$genetype
PATH_plink <- opt$plinkpath
memmb <- opt$memmb

#---read some files
geneannodt <- fread(geneannofn)
colnames(geneannodt)[1:6] <- c('geneid','gene.name','chromosome','genest','geneend','strand')
bimdt <- fread(paste0(targplkpath,'/',targplkpref,'.bim'))
#--subset annotation using bimfile
if(!is.na(chrn)){
geneannodt_chr <- geneannodt[chromosome==chrn]} else {
geneannodt_chr <- geneannodt}
Ntotgene <- nrow(geneannodt_chr)
#gene <- 'ENSBTAG00000021242'
#=========loop for each gene=========
for (gene in geneannodt_chr$geneid){
geneorder <- geneannodt_chr[,which(geneid==gene)]
cat(paste0('Processing gene ',gene,' (',geneorder ,' gene of ',Ntotgene,' on chrom ',chrn,') at ',Sys.time(),'\n'),sep='\n')
#---create SNP cood around gene
genecoord <- geneannodt_chr[geneid==gene]
if(phetype==1){#expression phe is gene
if(genecoord$strand==1){center <-genecoord$genest}else if(genecoord$strand== -1){center <-genecoord$geneend}
snpst <- center-genedis
if(snpst<0){snpst<-0}
snpend <- center+genedis
}else if(phetype==2) {#expression phe is gene is intron or exon
snpst <- genecoord$genest-genedis
if(snpst<0){snpst<-0}
snpend <- genecoord$geneend+genedis
}else {stop('phetype options: 1=gene,2=intron/exon')}
#---find cis SNPs for gene
bimdt_sub <- bimdt[V4>=snpst&V4<=snpend&V1==chrn]
#---save snplist
if(nrow(bimdt_sub)>=minNsnp){
write.table(bimdt_sub$V2,paste0(snplistoutpath,'/',gene,'.snplist'),quote=F,col.names=F,row.names=F)
}
#--create plink argument
arg <- paste0(PATH_plink," --memory ",memmb," --keep-allele-order --cow --bfile ",targplkpath,"/",targplkpref, " --extract ",snplistoutpath,"/",gene,".snplist"," --score <(zcat ",snpvfn,")", " 1 2 3 --out ",gene)
#--process argument to protect '()'
arg1 <- paste0('echo ','\"',arg,'\"','|bash')
#--excute argement
system(arg1)
#--copy zipped profile to outpath 
if ( file.exists(paste0(gene,'.profile')) ){
arg2 <- paste0('cat ',gene,'.profile|gzip >',profileoutpath,'/',gene,'.profile.gz')
system(arg2)}else{cat(paste0('WARNING: Profile of ',gene,' not estimated possibly due to too few SNPs'),'\n')}
arg2 <- paste0('cat ',gene,'.profile|gzip >',profileoutpath,'/',gene,'.profile.gz')
system(arg2)
#--rm intermediate files
arg3 <- paste0('rm ',gene,'*')
system(arg3)
}
cat(paste0('All analysis finished at ',Sys.time()),sep='\n')
cat(paste0('Results saved to ',profileoutpath),sep='\n')

