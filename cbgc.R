#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<5) {
  stop("5 argument must be supplied", call.=FALSE)
}

suppressMessages(library("optparse"))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))

option_list = list(
make_option("--suffx", action="store", default=NULL, type='character'),
make_option("--inpath", action="store", default=NULL, type='character'),
make_option("--mergfiles", action="store", default='n', type='character'),
make_option("--outpath", action="store", default='.', type='character'),
make_option("--outpref", action="store", default='test', type='character')
)
opt = parse_args(OptionParser(option_list=option_list))

suffx <- opt$suffx
inpath <- opt$inpath
mergop <- opt$mergfiles
outpath <- opt$outpath
outpref <- opt$outpref

#---get gene list in the folder
flist <- list.files(inpath,paste0('*.',suffx))
ngene <- length(flist)
recdseq <- seq(1,ngene,2000)
#---loop reading
profllooplist <- list()
if(mergop=='y'){
#--merge individual data
for (profl in flist){
if(grep(profl,flist) %in% recdseq){cat(paste0('Processing ',grep(profl,flist),' of ',ngene,' genes at ',Sys.time()),'\n')}
gene <- strsplit(profl,'\\.')[[1]][1]
dt <- fread(paste0(inpath,'/',profl))
colnames(dt)[grep('SCORE',colnames(dt))] <- gene
profllooplist[[profl]] <- dt[,grep(paste0('IID|',gene),colnames(dt)),with=F]}
cat(paste0('Combining profole across genes at ',Sys.time()),'\n')
#---combine profiles
cbprofl <- Reduce(function(...) merge(...,by='IID',all=T,sort=F), profllooplist)
} else if(mergop=='n'){
#--cbind data
for (profl in flist){
if(grep(profl,flist) %in% recdseq){cat(paste0('Processing ',grep(profl,flist),' of ',ngene,' genes at ',Sys.time()),'\n')}
gene <- strsplit(profl,'\\.')[[1]][1]
dt <- fread(paste0(inpath,'/',profl))
colnames(dt)[grep('SCORE',colnames(dt))] <- gene
if(grep(profl,flist)==1){profllooplist[[profl]] <- dt[,grep(paste0('IID|',gene),colnames(dt)),with=F]}else{profllooplist[[profl]] <- dt[,grep(paste0(gene),colnames(dt)),with=F]}
}
cat(paste0('Combining profole across genes at ',Sys.time()),'\n')
#---combine profiles
cbprofl <- bind_cols(profllooplist)
} else {stop('options for merge should be y or n')}
#--save results
write.table(cbprofl,gzfile(paste0(outpath,'/',outpref,'.profile.tbl.gz')),row.name=F,quote=F,sep=' ')

