#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("4 argument must be supplied", call.=FALSE)
}

suppressMessages(library("optparse"))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(OmicKriging)))
suppressMessages(suppressWarnings(library(lqmm)))
suppressMessages(suppressWarnings(library(gaston)))

#---parsing options
option_list = list(
make_option("--bv", action="store", default=NA, type='character'),
make_option("--exp", action="store", default=NA, type='character'),
make_option("--include", action="store", default=NA, type='character'),
make_option("--genename", action="store", default=NA, type='character'),
make_option("--cfe", action="store", default=NULL, type='character'),
make_option("--qfe", action="store", default=NULL, type='character'),
make_option("--GRM", action="store", default=NULL, type='character'),
make_option("--outpath", action="store", default=".", type='character'),
make_option("--outpref", action="store", default="test", type='character'),
make_option("--thread", action="store", default=1, type='integer')
)
opt = parse_args(OptionParser(option_list=option_list))

bvfn <- opt$bv 
expfn<- opt$exp 
includfn <- opt$include 
phenamefn <- opt$genename
cfefn <- opt$cfe 
qfefn <- opt$qfe 
grmpref <- opt$GRM 
outpath <- opt$outpath 
outpref <- opt$outpref

#---read in GRM
setThreadOptions(opt$thread)
grmname <- paste0(grmpref,'.bin')
if(!is.null(grmpref)){
cat(paste0('Reading GRM started at ',Sys.time()),'\n')
GRM <- read_GRMBin(grmpref)
GRM <- make.positive.definite(GRM)
cat(paste0('Making PD GRM finished at ',Sys.time()),'\n')
}else if(is.null(grmpref)) {
cat(paste0('No options of GRM is specified, then run simple linear regression'),'\n')
} else if (!file.exists(grmname)){
stop('Option of GRM is specified, but the file does not exist!')
} 
#---read in other data
Xdt <- fread(bvfn,header=T)
Nf <- max(count.fields(expfn))
if(Nf<=100000){
phedt <- fread(expfn,header=F)
}else{
cat(paste0('Note: Slow reading due to large N (>100k) of columns in ',expfn),'\n')
cat(paste0('Consider analysing expression data chromosome by chromosome'),'\n')
suppressMessages(suppressWarnings(library(vroom)))
phedt <- vroom(expfn,col_names = F)
setDT(phedt)
}
if(!is.na(includfn)){
if(length(includlist)!=0&file.exists(includfn)){
includlist <- fread(includfn,header=F)
cat(paste0(nrow(includlist),' individuals are to be included'),'\n')
phedt <- phedt[V2 %in% unlist(includlist$V1)]
cat(paste0('After filter, phenotype data of ',nrow(phedt),' will be analysed'),'\n')
}else{stop('Include file cant be read, check data')}
}
phename <- fread(phenamefn,header=F)

colnames(phedt)[-c(1:2)] <- unlist(phename)

if(length(cfefn)!=0&isTRUE(file.exists(as.character(cfefn)))){
cat(paste0('Reading c fixed effect'),'\n')
cfedt <- fread(cfefn,header=F)
cfecols <- colnames(cfedt)[-c(1:2)]
cfedt[,(cfecols):=lapply(.SD, as.factor),.SDcols=cfecols]
}
#colnames(cfedt)[3:5] <- c('speci','breed','exp')
if(length(qfefn)!=0&isTRUE(file.exists(as.character(qfefn)))){
cat(paste0('Reading q fixed effect'),'\n')
qfedt <- fread(qfefn,header=F)}
if(exists('cfedt')&exists('qfedt')) {
cbfe <- merge(cfedt[,-1],qfedt[,-1],by='V2',sort=F,all=T)} else if(exists('cfedt')&!exists('qfedt')) {
cat(paste0('No quantitative fixed effects are supplied'),'\n')
cbfe <- cfedt[,-1]} else if (!exists('cfedt')&exists('qfedt')) {
cat(paste0('No categorical fixed effects are supplied'),'\n')
cbfe <- qfedt[,-1]} else {
cat(paste0('No fixed effects are supplied'),'\n')
#cbfe <- data.frame(NA)
}
#if(exists('cbfe')){head(cbfe,1)}
#--analysis starts
cat(paste0('Data sorting starts at ',Sys.time()),sep='\n')
#--remove files without input values
if(!is.null(grmpref)&isTRUE(file.exists(as.character(grmname)))) {
if(exists('cbfe')){indiovlist <- list(phedt$V2,cbfe$V2,Xdt$IID,rownames(GRM))
} else { indiovlist <- list(phedt$V2,Xdt$IID,rownames(GRM))}
indiovlist <- indiovlist[lapply(indiovlist,length)>0]
} else if(is.null(grmpref)){
if(exists('cbfe')){indiovlist <- list(phedt$V2,cbfe$V1,Xdt$IID)
}else { indiovlist <- list(phedt$V1,Xdt$IID)}
indiovlist <- indiovlist[lapply(indiovlist,length)>0]
}
#--determine common individuals
commindi <- Reduce(intersect,indiovlist)
Ncommindi <- length(commindi)
if(Ncommindi>3){cat(paste0(Ncommindi,' individuals in common in supplied files'),'\n')}else{stop('too few (n<=3) individuals across common across files for analysis')}
#--determine common genes
genes_comm <- unique(intersect(colnames(Xdt),colnames(phedt)))
if(length(genes_comm)==0){stop(paste0(length(genes_comm),' genes in common between BV file and expression files, check data'))}
cat(paste0(length(genes_comm),' genes in common between X file and phenotype files'),'\n')
#--subset cols of Xdt due to genes available
Xdt_sub <- cbind(Xdt[,1],scale(Xdt[,colnames(Xdt) %in% genes_comm,with=F]))
phedt_sub <- cbind(phedt[,2],phedt[,colnames(phedt) %in% genes_comm,with=F])

#--sort data by rows---
if(!is.null(grmpref)&isTRUE(file.exists(as.character(grmname)))){
#-subset GRM
GRM <- GRM[rownames(GRM) %in% commindi,colnames(GRM) %in% commindi]
#-organising other data according to subset GRM
phedt1 <- setDF(phedt_sub[match(rownames(GRM),unlist(phedt_sub[,1])),])
Xdt1 <- setDF(Xdt_sub[match(rownames(GRM),unlist(Xdt_sub$IID)),])
if(exists('cbfe')) {cbfe1 <- setDF(cbfe[match(rownames(GRM),unlist(cbfe[,1])),])}
} else if(is.null(grmpref)){
#phedt1 <- setDF(phedt_sub[match(rownames(GRM),unlist(phedt_sub[,1])),])
phedt1 <- setDF(phedt_sub[match(commindi,unlist(phedt_sub[,1])),])
#Xdt1 <- setDF(Xdt_sub[match(rownames(GRM),unlist(Xdt_sub$IID)),])
Xdt1 <- setDF(Xdt_sub[match(commindi,unlist(Xdt_sub$IID)),])
if(exists('cbfe')) {
#cbfe1 <- setDF(cbfe[match(rownames(GRM),unlist(cbfe$V2)),])}
cbfe1 <- setDF(cbfe[match(commindi,unlist(cbfe[,1])),])
 }
}
#--sort exp data by cols
phedt1.1 <- cbind(phedt1[,1],setDF(phedt1)[,genes_comm])
Xdt1.1 <- cbind(Xdt1[,1],setDF(Xdt1)[,genes_comm])
#--drop fixed effects if variable < 2
if(exists('cbfe1')){cbfe1.1 <- cbfe1[,which(apply(cbfe1, 2, function(x) length(unique(x)))>=2)]
#print(head(cbfe1.1,1))
}

#---check rows and cols
if(!is.null(grmpref)&isTRUE(file.exists(as.character(grmname)))) {
if(exists('cbfe1.1')){
if(all(sapply(list(phedt1[,1],Xdt1.1[,1],rownames(GRM),colnames(GRM)), FUN = identical, cbfe1.1[,1]))){cat(paste0("Individuals of all files are sorted out"),'\n')}else{stop('Problems of soring indivuals between files occur, check data')}
}else {
if(all(sapply(list(phedt1[,1],Xdt1.1[,1],rownames(GRM),colnames(GRM)), FUN = identical, colnames(GRM)))){cat(paste0("Individuals of all files are sorted out"),'\n')}else{stop('Problems of soring indivuals between files occur, check data')}
 }
} else if (is.null(grmpref)) {
if(exists('cbfe1.1')){
if(all(sapply(list(phedt1[,1],Xdt1.1[,1]), FUN = identical, cbfe1.1[,1]))){cat(paste0("Individuals of all files are sorted out"),'\n')}else{stop('Problems of soring indivuals between files occur, check data')}
}else {
if(all(sapply(list(phedt1[,1]), FUN = identical, Xdt1.1[,1]))){cat(paste0("Individuals of all files are sorted out"),'\n')}else{stop('Problems of soring indivuals between files occur, check data')}
 }
}
if(identical(colnames(phedt1.1)[-1],colnames(Xdt1.1)[-1])){cat(paste0("Genes of all files are sorted out"),'\n')}else{stop('Problems of soring genes between files occur, check data')}

#---eigen-decomposition of GRM
if(!is.null(grmpref)&isTRUE(file.exists(as.character(grmname)))) {
eiGRM <- eigen(GRM)
cat(paste0('GRM eigen-decomposition finished at ',Sys.time()),'\n')
}

#---col-wise lmm function
if(!is.null(grmpref)&isTRUE(file.exists(as.character(grmname)))){
lmmcolfunc <- function(expcol){
reportseg <- seq(1,length(genes_comm),2000)
if(expcol %in% reportseg){cat(paste0('Gene ',expcol,' started at ',Sys.time()),sep='\n')}
if(exists('cbfe1.1')){
#testdt <- cbind(phedt1.1[,c(1,expcol+1)],Xdt1.1[,expcol+1],cbfe1.1[,-1])
testdt <- cbind(Xdt1.1[,c(1,expcol+1)],phedt1.1[,expcol+1],cbfe1.1[,-1])
gene <- colnames(phedt1.1)[expcol+1]
ycol <-  colnames(testdt)[3]
xcols <- colnames(cbfe1.1)[-1]
if(expcol==1){cat(paste0(length(xcols),' fixed effects are accounted for'),'\n')
}
f <- as.formula(paste(ycol,'~',paste(xcols,collapse='+')))
#---corect other fixed effects
testdt$gene.resi <- summary(lm(f,testdt))$residuals
#---gaston reml-blup
#y <- testdt$gene.resi
X <- cbind(1,testdt$gene.resi)
} else {#if no fixed effects
testdt <- cbind(Xdt1.1[,c(1,expcol+1)],phedt1.1[,expcol+1])
#print(head(testdt))
gene <- colnames(phedt1.1)[expcol+1]
ycol <- colnames(testdt)[3]
#y <- testdt[[ycol]]}
X <- cbind(1,testdt[[ycol]])}
#---check analysis data
y <- testdt[,grep(gene,colnames(testdt))]
invisible(capture.output(fit <- lmm.diago(y, X, eiGRM, verbose = T)))
b <- ifelse(abs(fit$BLUP_beta[2])>1,round(fit$BLUP_beta[2],1),round(fit$BLUP_beta[2],3))
se <- ifelse(abs(sqrt(fit$varbeta[2,2]))>1,round(sqrt(fit$varbeta[2,2]),1),round(sqrt(fit$varbeta[2,2]),3))
t <- ifelse(abs(b/se)>1,round(b/se,1),round(b/se,3))
p <- 2*pnorm(abs(t), lower.tail = F)
return(data.frame(gene=gene,b,se,t,p))
 }
} else if(is.null(grmpref)) {
lmmcolfunc <- function(expcol){
reportseg <- seq(1,length(genes_comm),2000)
if(expcol %in% reportseg){cat(paste0('Gene ',expcol,' started at ',Sys.time()),sep='\n')}
if(exists('cbfe1.1')){
testdt <- cbind(Xdt1.1[,expcol+1],phedt1.1[,expcol+1],cbfe1.1[,-1])
} else {#if no fixed effects
testdt <- cbind(data.frame(Xdt1.1[,expcol+1],phedt1.1[,expcol+1]))
}
gene <- colnames(phedt1.1)[expcol+1]
ycol <-  colnames(testdt)[1]
xcols <- colnames(testdt)[-1]
f <- as.formula(paste(ycol,'~',paste(xcols,collapse='+')))
lmres <- summary(lm(f,testdt))
res_gene <- data.frame(gene,t(lmres$coefficients[grep('phedt',rownames(lmres$coefficients)),]))
colnames(res_gene) <- c('gene','b','se','t','p')
res_gene$b <- ifelse(abs(res_gene$b)>1,round(res_gene$b,1),round(res_gene$b,3))
res_gene$se <- ifelse(abs(res_gene$se)>1,round(res_gene$se,1),round(res_gene$se,3))
res_gene$t <- ifelse(abs(res_gene$t)>1,round(res_gene$t,1),round(res_gene$t,3))
return(res_gene)
 }
}

#---start the analysis----
cat(paste0('Association test starts at ',Sys.time()),sep='\n')
reslist <- list()
#reslist <- lapply(1:100,function(x) lmcolfunc(x))
reslist <- lapply(1:length(genes_comm),function(x) lmmcolfunc(x))
cat(paste0('Association test finishes at ',Sys.time()),sep='\n')
genedt <- setDT(do.call(rbind,reslist))
#---save data
if(!is.null(grmpref)&isTRUE(file.exists(as.character(grmname)))){
write.table(genedt,gzfile(paste0(outpath,'/',outpref,'.gosr.lmm.gz')),row.names=F,quote=F,sep=' ')
cat(paste0('Gene results saved to ',outpath,'/',outpref,'.gosr.lmm.gz'),sep='\n')
} else if(is.null(grmpref)){
write.table(genedt,gzfile(paste0(outpath,'/',outpref,'.gsor.lr.gz')),row.names=F,quote=F,sep=' ')
cat(paste0('Gene results saved to ',outpath,'/',outpref,'.gosr.lr.gz'),sep='\n')
}

