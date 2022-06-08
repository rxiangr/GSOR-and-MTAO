#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<5) {
  stop("5 argument must be supplied", call.=FALSE)
}


suppressMessages(library("optparse"))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))

#---parsing options
option_list = list(
make_option("--inpath", action="store", default=NA, type='character'),
make_option("--pattern", action="store", default=NA, type='character'),
make_option("--geneanno", action="store", default=NA, type='character'),
make_option("--outpath", action="store", default=".", type='character'),
make_option("--outpref", action="store", default="test", type='character')
)
opt = parse_args(OptionParser(option_list=option_list))

inpath <- opt$inpath
patt <- opt$pattern
geneannofn <- opt$geneanno
outpath <- opt$outpath
outpref <- opt$outpref

#---read in data
geneannodt <- fread(geneannofn)
colnames(geneannodt)[1:6] <- c('gene','gene.name','chromosome','genest','geneend','strand')
geneannodt[,gene.name:=ifelse(is.na(gene.name)|gene.name=='','unknown',gene.name)]
flist <- list.files(inpath,patt)
#---determine N traits
trvec <- unique(grep('^tr',unlist(strsplit(flist,'\\.')),value=T))
#---determine N chrn
chrvec <- gsub('chr','',unique(grep('chr',unlist(strsplit(flist,'\\.')),value=T)))
#---loop reading
looplist <- list()
for (tr in trvec){
cat(paste0('Processing ',tr,' at ',Sys.time()),'\n')
trlist <- grep(tr,flist,value=T)
trdt <- setDT(bind_rows(lapply(paste0(inpath,'/',trlist),fread)))
colnames(trdt)[-1] <- paste0(tr,'.',colnames(trdt)[-1])
looplist[[tr]] <- trdt
}
cbdt <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2,by="gene",all.x=T,sort=F),looplist)
cbdt1 <- merge(geneannodt,cbdt,by='gene',sort=F)
setorder(cbdt1,chromosome,genest)
write.table(cbdt1,gzfile(paste0(outpath,'/',outpref,'.gsor.sum.gz')),row.names=F,quote=F,sep=' ')
cat(paste0('All summerised results saved to ',outpath,'/',outpref,'.gsor.sum.gz'),sep='\n')


#----Get Pn
#--decorrfunc is from Jordan et al 2019 HOPS
decorrfunc <- function(mat,cormat){
RefCor <- do.call("rbind", lapply(1:ncol(cormat),
        function(i) cbind.data.frame(row = 1:(i - 1), col = i, ref = (sum((1:(i - 1)) - 1) + 1):sum(0:(i - 1)))))
    CorrelatedPairs <- which(abs(cormat[upper.tri(cormat, diag = F)]) > 0.8)
    if (length(CorrelatedPairs) > 0) {
        CorrelatedPairs <- cbind.data.frame(Trait1 = row.names(cormat)[RefCor$row[match(CorrelatedPairs, RefCor$ref)]], Trait2 = row.names(cormat)[RefCor$col[match(CorrelatedPairs, RefCor$ref)]], CorLevel = cormat[upper.tri(cormat, diag = F)][CorrelatedPairs])
    warning("Some traits have pairwise correlation above the required cut-off of 0.8, please select one trait out of the pair (e.g. highest heritability). NB: the set of pairs is returned instead of the whitened Z-scores.")
    return(CorrelatedPairs)
    }   else {
        "%^%" <- function(x, n) with(svd(x), u %*% (d^n * t(v)))
        Zstar <- cormat %^% (-1/2) %*% t(mat)
        Zstar <- data.frame(t(Zstar))
        colnames(Zstar) <- colnames(mat)
        row.names(Zstar) <- row.names(mat)
        return(Zstar)
    }
}

uniq.tdt <- unique(cbdt1[,grep('^gene$|\\.t$',colnames(cbdt1)),with=F])
ttbl <- as.matrix(uniq.tdt[,-1])
rownames(ttbl) <- unlist(uniq.tdt[,1])
cormat <- cor(ttbl,use='pairwise.complete.obs')
decorr.ttbl <- decorrfunc(ttbl,cormat)
Z0star2 <- 2
#---get N of triats
Pn <- rowSums(abs(decorr.ttbl) > Z0star2, na.rm = TRUE)
Pn_Pvalue <- pbinom(q = Pn - 1, size = ncol(decorr.ttbl),
            prob = 0.045, lower.tail = FALSE, log.p = FALSE)
score <- setDT(data.frame(gene=rownames(decorr.ttbl), Pn, Pn_Pvalue))
score[,Pn_Pvalue.adj:=p.adjust(Pn_Pvalue,'fdr')]
score1 <- merge(geneannodt,score,by='gene',sort=F)

#---get Pm
tv <- as.data.frame(ttbl)
setDT(tv)
f_dowle2 = function(DT) {
  for (i in names(DT))
    DT[!is.finite(get(i)), (i):=0]
    }
f_dowle2(tv)
tv <- cbind(rn=row.names(ttbl),tv)
ivcormat <- solve(cor(tv[,-1,with=F]))
tv[,df:=rowSums(!is.na(tv[,-1]))]
tv[,x:=sqrt(rowSums(as.matrix(tv[,-c(1,ncol(tv)),with=F]) %*% ivcormat * as.matrix(tv[,-c(1,ncol(tv)),with=F])))]
tv[,X2.p:=pchisq(x^2,df=df,lower.tail=F)]
tv[,X2.p.adj:=p.adjust(X2.p,'fdr')]
colnames(tv)[1] <- 'gene'
#head(tv)
score2 <- merge(score1,tv[,c(1,ncol(tv)-2,ncol(tv)-1,ncol(tv)),with=F],by='gene',sort=F)
write.table(score2,gzfile(paste0(outpath,'/',outpref,'.mtao.pleio.gz')),row.names=F,quote=F,sep=' ')
cat(paste0('Results of pleio saved to ',outpath,'/',outpref,'.mtao.pleio.gz'),sep='\n')

