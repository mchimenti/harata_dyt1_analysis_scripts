library(illuminaMousev2.db)
library("beadarray")
library("genefilter")
library("limma")
library("hexbin")
library("tidyverse")

## Volcano Plot function 
volcanoplot2 <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(logFC, neglog.P.Val, pch=20, main=main, ...))
  with(subset(res, adj.P.Val<sigthresh ), points(logFC, neglog.P.Val, pch=20, col="blue", ...))
  #with(subset(res, abs(logFC)>lfcthresh), points(logFC, neglog.P.Val, pch=20, col="orange", ...))
  #with(subset(res, adj.P.Val<sigthresh & abs(logFC)>lfcthresh), points(logFC, neglog.P.Val, pch=20, col="green", ...))
  #if (labelsig) {
  #  require(calibrate)
  #  with(subset(res, adj.P.Val<sigthresh & abs(logFC)>lfcthresh), textxy(logFC, adj.P.Val, labs=SYMBOL, cex=textcx, offset=0.3, ...))
  #}
  #legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

setwd("~/collab_proj/harata/harata_mar_2018_reanalysis/chimenti_analysis_2018/")

chip_id <- read.csv("harata_sample_sheet_ordered_corrected.csv",sep = ',', header=TRUE)
chip_id$chipID <- paste(chip_id$Sentrix_ID,chip_id$Sentrix_Position,sep="_")
phenodata <- read.csv("samples.csv")
phenotype <- cbind(phenodata,chip_id)
technical <- phenotype[phenotype$Technical == 1,]
phenotype <- phenotype[phenotype$Technical == 0,]

idatFiles <- list.files(path="./raw_IDAT/",pattern=".idat",full.names=T)
y <- readIdatFiles(idatFiles)
y <- y[,colnames(y) %in% phenotype$chipID]
y <- addFeatureData(y)
dim(y)

boxplot(log2(exprs(y)))
hist(log2(exprs(y)))

eset <- normaliseIllumina(y, method="quantile",  transform="log2", status = y@featureData$Status)
eset<- eset[eset@featureData$Status == 'regular',]
annotated <- !is.na(fData(eset)$SYMBOL)
table(annotated)
eset <- eset[annotated,]
rm(annotated)
class(eset) <- "ExpressionSet"
rownames(phenotype) <- phenotype$chipID
pData(eset)<- phenotype

## sanity check pData rownames should order and match colnames expression data
all(rownames(pData(eset) %in% colnames(exprs(eset))))  


genotype <- factor(phenotype$Genotype,levels=c("wt","het"))
tissue <- factor(phenotype$Tissue,levels=c("Cer","Cor","Hip","Str"))
design <- model.matrix(~0+genotype+tissue)
fit_global <- lmFit(eset,design=design)

#hist(fit_global$Amean)
#plotSA(fit_global)
#keep <- fit_global$Amean > 7


contrast.matrix <- makeContrasts(genotypehet-genotypewt, levels=design)
fit_global_2 <- contrasts.fit(fit_global, contrast.matrix)
fit_global_2 <- eBayes(fit_global_2, trend = TRUE)
#fit_global_2 <- eBayes(fit_global_2[keep,], trend=TRUE)

res_global <- limma::topTable(fit_global_2, coef = 1, adjust.method="fdr", number=Inf, p.value = 1)
#res_global <- res_global[,4:14]
res_global_sig <- filter(res_global, adj.P.Val < 0.05)

limma::volcanoplot(fit_global_2, coef = 1, highlight = 10, xlim = c(-0.5,0.5))

signif = 0.01
res_global$neglog.P.Val <- -log10(res_global$adj.P.Val)

##plotting 
pdf("volcano_global_het_vs_wt_MAY2019.pdf")
volcanoplot2(res = res_global, lfcthresh = 1, sigthresh = signif,
             xlim = c(-2,2), ylim = c(0, 7), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Global Het vs Wt")
dev.off()

setEPS()
postscript("volcano_global_het_vs_wt_MAY2019.eps")
volcanoplot2(res = res_global, lfcthresh = 1, sigthresh = signif,
             xlim = c(-2,2), ylim = c(0, 7), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Global Het vs Wt")
dev.off()

## adding ENS IDs for pathway
x <- illuminaMousev2ENSEMBL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
z <- unlist(xx)
w <- names(z)
z <- as.vector(z)
illumina2ens <- data.frame(z,w, stringsAsFactors = FALSE)
names(illumina2ens) <- c("EnsemblID","IlluminaID")

res_global_sig$ENSEMBL = mapIds(org.Mm.eg.db,
                   keys=as.character(res_global_sig$SYMBOL), 
                   column="ENSEMBL",
                   keytype="SYMBOL",
                   multiVals="first")

write.csv(x=res_global_sig, file = "DEgenes_GLOBAL_Het_v_WT_MAY2019.csv")
