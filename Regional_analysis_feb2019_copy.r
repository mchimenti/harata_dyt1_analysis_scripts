## 2019 RE-Reanalysis Harata Microarray data 

## This is a RE-RE-RE analysis of this data in Feb 2019

#source("https://bioconductor.org/biocLite.R")
#biocLite("beadarray")
#biocLite("genefilter")
#biocLite("limma")
#biocLite("illuminaMousev2.db")

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

## comment these below two lines if you want to INCLUDE tech reps
#technical <- phenotype[phenotype$Technical == 1,]
#phenotype <- phenotype[phenotype$Technical == 0,]

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

## write out the raw data eset w/ ENSEMBL gene IDs
## coordinate Illumina IDs to Ensembl IDs for GAGE analysis 
x <- illuminaMousev2ENSEMBL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
z <- unlist(xx)
w <- names(z)
z <- as.vector(z)
illumina2ens <- data.frame(z,w, stringsAsFactors = FALSE)
names(illumina2ens) <- c("EnsemblID","IlluminaID")

dat_mat <- exprs(eset)
datf <- as.data.frame(dat_mat)
datf <- rownames_to_column(datf, var = "IlluminaID")
datf <- inner_join(illumina2ens, datf)
datf$gene = mapIds(org.Mm.eg.db,
                               keys=datf$EnsemblID, 
                               column="SYMBOL",
                               keytype="ENSEMBL",
                               multiVals="first")


write.csv(file = "full_exprs_datatable_EnsemblID_GENENAMES.csv", x = datf)
###################
## Correlation technical replicates
###################

##########
## pearson corr btw Het 5 Fem Cer and Het-5 Fem Cer B tech rep
cor(exprs(eset)[, sampleNames(eset) == "9974433015_H"], exprs(eset)[, sampleNames(eset) == "9974433003_H"])

x <- exprs(eset)[, sampleNames(eset) == "9974433015_H"]
y <- exprs(eset)[, sampleNames(eset) == "9974433003_H"]
reg <- lm(x ~ y)

pdf("het5_fem_cer_technical_rep_corr.pdf")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 5 Female Cer vs. Tech Rep (r = 0.955)')
hexVP.abline(p$plot.vp, reg)
dev.off()

setEPS()
postscript("het5_fem_cer_technical_rep_corr.eps")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 5 Female Cer vs. Tech Rep (r = 0.955)')
hexVP.abline(p$plot.vp, reg)
dev.off()

#############
## cor het 1 male hip tech rep
cor(exprs(eset)[, sampleNames(eset) == "9704007100_F"], exprs(eset)[, sampleNames(eset) == "9974433015_C"])
x <- exprs(eset)[, sampleNames(eset) == "9704007100_F"]
y <- exprs(eset)[, sampleNames(eset) == "9974433015_C"]
reg <- lm(x ~ y)
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 1 Male Hip vs. Tech Rep (r = 0.991)')
hexVP.abline(p$plot.vp, reg)

pdf("het1_male_hip_technical_rep_corr.pdf")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 1 Male Hip vs. Tech Rep (r = 0.991)')
hexVP.abline(p$plot.vp, reg)
dev.off()

setEPS()
postscript("het1_male_hip_technical_rep_corr.eps")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 1 Male Hip vs. Tech Rep (r = 0.991)')
hexVP.abline(p$plot.vp, reg)
dev.off()



## cor het 2 male cortex 
cor(exprs(eset)[, sampleNames(eset) == "9967500133_G"], exprs(eset)[, sampleNames(eset) == "9974433015_E"])
x <- exprs(eset)[, sampleNames(eset) == "9967500133_G"]
y <- exprs(eset)[, sampleNames(eset) == "9974433015_E"]
reg <- lm(x ~ y)

pdf("het2_male_ctx_tech_rep_corr.pdf")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 2 Male Cor vs. Tech Rep (r = 0.945)')
hexVP.abline(p$plot.vp, reg)
dev.off()

setEPS()
postscript("het2_male_ctx_tech_rep_corr.eps")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 2 Male Cor vs. Tech Rep (r = 0.945)')
hexVP.abline(p$plot.vp, reg)
dev.off()


## cor het 9 male str
cor(exprs(eset)[, sampleNames(eset) == "9704007102_E"], exprs(eset)[, sampleNames(eset) == "9974433015_B"])
x <- exprs(eset)[, sampleNames(eset) == "9704007102_E"]
y <- exprs(eset)[, sampleNames(eset) == "9974433015_B"]
reg <- lm(x ~ y)

pdf("het9_male_str_tech_rep_corr.pdf")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 9 Male Str vs. Tech Rep (r = 0.990)')
hexVP.abline(p$plot.vp, reg)
dev.off()

setEPS()
postscript("het9_male_str_tech_rep_corr.eps")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 9 Male Str vs. Tech Rep (r = 0.990)')
hexVP.abline(p$plot.vp, reg)
dev.off()


## cor het 8 fem hip
cor(exprs(eset)[, sampleNames(eset) == "9969472044_F"], exprs(eset)[, sampleNames(eset) == "9974433015_D"])
x <- exprs(eset)[, sampleNames(eset) == "9969472044_F"]
y <- exprs(eset)[, sampleNames(eset) == "9974433015_D"]
reg <- lm(x ~ y)

pdf("het8_fem_hip_tech_rep_corr.pdf")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 8 Female Hip vs. Tech Rep (r = 0.943)')
hexVP.abline(p$plot.vp, reg)
dev.off()

setEPS()
postscript("het8_fem_hip_tech_rep_corr.eps")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 8 Female Hip vs. Tech Rep (r = 0.943)')
hexVP.abline(p$plot.vp, reg)
dev.off()


## cor het 7 male cer 
cor(exprs(eset)[, sampleNames(eset) == "9969472051_H"], exprs(eset)[, sampleNames(eset) == "9974433015_G"])
x <- exprs(eset)[, sampleNames(eset) == "9969472051_H"]
y <- exprs(eset)[, sampleNames(eset) == "9974433015_G"]
reg <- lm(x ~ y)

pdf("het7_male_cer_tech_rep_corr.pdf")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 7 Male Cer vs. Tech Rep (r = 0.971)')
hexVP.abline(p$plot.vp, reg)
dev.off()

setEPS()
postscript("het7_male_cer_tech_rep_corr.eps")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 7 Male Cer vs. Tech Rep (r = 0.971)')
hexVP.abline(p$plot.vp, reg)
dev.off()


## Het 4 male cor
cor(exprs(eset)[, sampleNames(eset) == "9967500139_G"], exprs(eset)[, sampleNames(eset) == "9974433015_F"])
x <- exprs(eset)[, sampleNames(eset) == "9967500139_G"]
y <- exprs(eset)[, sampleNames(eset) == "9974433015_F"]
reg <- lm(x ~ y)

pdf("het4_male_ctx_tech_rep_corr.pdf")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 4 Male Cor vs. Tech Rep (r = 0.972)')
hexVP.abline(p$plot.vp, reg)
dev.off()

setEPS()
postscript("het4_male_ctx_tech_rep_corr.eps")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 4 Male Cor vs. Tech Rep (r = 0.972)')
hexVP.abline(p$plot.vp, reg)
dev.off()


## Het 6 female stri 
cor(exprs(eset)[, sampleNames(eset) == "9969472042_E"], exprs(eset)[, sampleNames(eset) == "9974433015_A"])
x <- exprs(eset)[, sampleNames(eset) == "9969472042_E"]
y <- exprs(eset)[, sampleNames(eset) == "9974433015_A"]
reg <- lm(x ~ y)

pdf("het6_female_str_tech_rep_corr_new.pdf")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 6 Female Str vs. Tech Rep (r = 0.968)')
hexVP.abline(p$plot.vp, reg)
dev.off()

setEPS()
postscript("het6_female_str_tech_rep_corr_new.eps")
hbin <- hexbin(x,y,xbins=50, shape=1.5, xlab = "Norm Expr", ylab = "Norm Expr")
p <- plot(hbin, main = 'Het 6 Female Str vs. Tech Rep (r = 0.968)')
hexVP.abline(p$plot.vp, reg)
dev.off()


######################
## making some plots on the data
###################### 

limma::plotMA(eset)
pcaData <- prcomp(t(exprs(eset)),scale = TRUE)
plot(pcaData$x,col=factor(phenotype$Set),pch=as.numeric(factor(phenotype$Tissue)))
plotDensities(exprs(eset), legend=FALSE)

library("tidyverse")
library("ggplot2")

pcadf <- as.data.frame(pcaData$x[,1:2])
pcadf <- rownames_to_column(pcadf, var = "chipID")
pcadf <- left_join(pcadf, phenotype, by = "chipID")

ggplot(pcadf, aes(PC1,PC2, shape = Genotype)) + geom_point(aes(colour = Tissue), size = 3) + facet_grid(cols = vars(Sex))

######################
## setting up the GLM design
######################

####################
## Global Het vs WT comparison
####################

genotype <- factor(phenotype$Genotype,levels=c("wt","het"))
design <- model.matrix(~0+genotype)
fit_global <- lmFit(eset,design=design)
hist(fit$Amean)
plotSA(fit)
keep <- fit$Amean > 7


contrast.matrix <- makeContrasts(genotypehet-genotypewt, levels=design)
fit_global_2 <- contrasts.fit(fit_global, contrast.matrix)
fit_global_2 <- eBayes(fit_global_2[keep,], trend=TRUE)

res_global <- limma::topTable(fit_global_2, coef = 1, adjust.method="fdr", number=Inf, p.value = 1)
res_global <- res_global[,4:14]
#limma::volcanoplot(fit_global_2, coef = 1, highlight = 50, names = fData(eset)$SYMBOL, xlim = c(-0.5,0.5))
res_global$neglog.P.Val <- -log10(res_global$adj.P.Val)
pdf("volcano_global_het_vs_wt.pdf")
volcanoplot2(res = res_global, lfcthresh = 1, sigthresh = signif,
             xlim = c(-2,2), ylim = c(0, 5), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Global Het vs Wt")
dev.off()

setEPS()
postscript("volcano_global_het_vs_wt.eps")
volcanoplot2(res = res_global, lfcthresh = 1, sigthresh = signif,
             xlim = c(-2,2), ylim = c(0, 5), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Global Het vs Wt")
dev.off()

write.csv(res_global, file = "res_global_het_wt_DEgenes.csv")


####################
## Two-factor design
####################
tissue <- factor(phenotype$Tissue,levels=c("Cer","Cor","Hip","Str"))
genotype <- factor(phenotype$Genotype,levels=c("wt","het"))
type <- paste(tissue,genotype,sep="_")
design <- model.matrix(~0+type)
colnames(design) <- c("Cer_Het","Cer_WT","Cor_Het","Cor_WT","Hip_Het","Hip_WT","Str_Het","Str_WT")



#####################
## fitting the GLM on the design 
#####################

fit <- lmFit(eset,design=design)
contrast.matrix <- makeContrasts(Cer_Het-Cer_WT, Cor_Het-Cor_WT, Hip_Het-Hip_WT, Str_Het-Str_WT, 
                                 Cer_WT-Cor_WT, Cer_WT-Str_WT, Cer_WT-Hip_WT, Hip_WT-Cor_WT,
                                 Str_WT-Cor_WT, Str_WT-Hip_WT, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#####GENOTYPE COMPARISONS

signif = 0.01

## Cerebellum Het v WT
res_cbl_het_vs_wt <- limma::topTable(fit2, coef = 1, adjust.method="fdr", number=Inf, p.value = 1)
res_cbl_het_vs_wt$neglog.P.Val <- -log10(res_cbl_het_vs_wt$adj.P.Val)
pdf("volcano_cbl_het_vs_wt.pdf")
volcanoplot2(res = res_cbl_het_vs_wt, lfcthresh = 2, sigthresh = signif, 
             xlim = c(-2,2), ylim = c(0, 5), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Cerebellum Het vs. WT")

dev.off()

setEPS()
postscript("volcano_cbl_het_vs_wt.eps")
volcanoplot2(res = res_cbl_het_vs_wt, lfcthresh = 2, sigthresh = signif, 
             xlim = c(-2,2), ylim = c(0, 5), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Cerebellum Het vs. WT")

dev.off()

## Cortex Het v WT
res_ctx_het_vs_wt <- limma::topTable(fit2, coef = 2, adjust.method="fdr", number=Inf, p.value = 1)
res_ctx_het_vs_wt$neglog.P.Val <- -log10(res_ctx_het_vs_wt$adj.P.Val)
pdf("volcano_ctx_het_vs_wt.pdf")
volcanoplot2(res = res_ctx_het_vs_wt, lfcthresh = 2, sigthresh = signif, 
             xlim = c(-2,2), ylim = c(0, 5), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Cortex Het vs. WT")

dev.off()

setEPS()
postscript("volcano_ctx_het_vs_wt.eps")
volcanoplot2(res = res_ctx_het_vs_wt, lfcthresh = 2, sigthresh = signif, 
             xlim = c(-2,2), ylim = c(0, 5), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Cortex Het vs. WT")

dev.off()

## Hippocampus Het v WT
res_hip_het_vs_wt <- limma::topTable(fit2, coef = 3,adjust.method="fdr", number=Inf, p.value = 1)
res_hip_het_vs_wt$neglog.P.Val <- -log10(res_hip_het_vs_wt$adj.P.Val)
pdf("volcano_hip_het_vs_wt.pdf")
volcanoplot2(res = res_hip_het_vs_wt, lfcthresh = 2, sigthresh = signif, 
             xlim = c(-2,2), ylim = c(0, 5), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Hippocampus Het vs. WT")

dev.off()

setEPS()
postscript("volcano_hip_het_vs_wt.eps")
volcanoplot2(res = res_hip_het_vs_wt, lfcthresh = 2, sigthresh = signif, 
             xlim = c(-2,2), ylim = c(0, 5), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Hippocampus Het vs. WT")

dev.off()

## Striatum Het v WT
res_str_het_vs_wt <- limma::topTable(fit2, coef = 4,adjust.method="fdr", number=Inf, p.value = 1)
res_str_het_vs_wt$neglog.P.Val <- -log10(res_str_het_vs_wt$adj.P.Val)
pdf("volcano_str_het_vs_wt.pdf")
volcanoplot2(res = res_str_het_vs_wt, lfcthresh = 2, sigthresh = signif, 
             xlim = c(-2,2), ylim = c(0, 5), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Striatum Het vs. WT")

dev.off()

setEPS()
postscript("volcano_str_het_vs_wt.eps")
volcanoplot2(res = res_str_het_vs_wt, lfcthresh = 2, sigthresh = signif, 
             xlim = c(-2,2), ylim = c(0, 5), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Striatum Het vs. WT")

dev.off()

##
write.csv(x = res_cbl_het_vs_wt, file = "res_DEgenes_cbl_het_vs_wt_padj_0p05.csv")
#write.csv(x = res_ctx_het_vs_wt, file = "res_DEgenes_ctx_het_vs_wt_padj_0p05.csv")
write.csv(x = as.matrix(0), file = "res_DEgenes_ctx_het_vs_wt_padj_0p05.csv")  # write a "file" of zeros b/c of NULL results
#write.csv(x = res_hip_het_vs_wt, file = "res_DEgenes_hip_het_vs_wt_padj_0p05.csv")
write.csv(x = as.matrix(0), file = "res_DEgenes_hip_het_vs_wt_padj_0p05.csv" )  
write.csv(x = res_str_het_vs_wt, file = "res_DEgenes_str_het_vs_wt_padj_0p05.csv")



####REGION COMPARISONS

y_lower = 2
y_upper = 70

## Cerebellum v Cortex WT
res_cer_cor <- limma::topTable(fit2, coef = 5, adjust.method="fdr", p.value = 1, number = Inf)
res_cer_cor <- res_cer_cor[,5:14]

res_cer_cor$neglog.P.Val <- -log10(res_cer_cor$adj.P.Val)
pdf("volcano_cbl_vs_ctx.pdf")
volcanoplot2(res = res_cer_cor, lfcthresh = 2, sigthresh = signif, 
             xlim = c(-7,7), ylim = c(y_lower, y_upper), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Cerebellum vs. Cortex")

dev.off()

setEPS()
postscript("volcano_cbl_vs_ctx.eps")
volcanoplot2(res = res_cer_cor, lfcthresh = 2, sigthresh = signif, 
             xlim = c(-7,7), ylim = c(y_lower, y_upper), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Cerebellum vs. Cortex")

dev.off()

## Cerebellum v Striatum WT
res_cer_str <- limma::topTable(fit2, coef = 6, adjust.method="fdr", p.value = 1, number = Inf)
res_cer_str$neglog.P.Val <- -log10(res_cer_str$adj.P.Val)
pdf("volcano_cbl_vs_str.pdf")
volcanoplot2(res = res_cer_str, lfcthresh = 1, sigthresh = signif, 
             xlim = c(-7,7), ylim = c(y_lower, y_upper), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Cerebellum vs Striatum")
dev.off()

setEPS()
postscript("volcano_cbl_vs_str.eps")
volcanoplot2(res = res_cer_str, lfcthresh = 1, sigthresh = signif, 
             xlim = c(-7,7), ylim = c(y_lower, y_upper), xlab = "log(Fold Change)", ylab = "-log(Adj. P-Value)", main = "Cerebellum vs Striatum")
dev.off()

## Cerebellum v Hippocampus WT
res_cer_hip <- limma::topTable(fit2, coef = 7, adjust.method="fdr", p.value = 1, number = Inf)
res_cer_hip$neglog.P.Val <- -log10(res_cer_hip$adj.P.Val)
pdf("volcano_cbl_v_hip.pdf")
volcanoplot2(res = res_cer_hip, lfcthresh = 1, sigthresh = signif,
             xlim = c(-7,7), ylim = c(y_lower, y_upper), xlab = "log Fold Change", ylab = "-log(Adj. P-Value)", main = "Cerebellum vs Hippocampus")
dev.off()

setEPS()
postscript("volcano_cbl_v_hip.eps")
volcanoplot2(res = res_cer_hip, lfcthresh = 1, sigthresh = signif,
             xlim = c(-7,7), ylim = c(y_lower, y_upper), xlab = "log Fold Change", ylab = "-log(Adj. P-Value)", main = "Cerebellum vs Hippocampus")
dev.off()

## Hippocampus v Cortex WT
res_hip_cor <- limma::topTable(fit2, coef = 8, adjust.method = "fdr", p.value = 1, number = Inf)
res_hip_cor$neglog.P.Val <- -log10(res_hip_cor$adj.P.Val)
pdf("volcano_hip_vs_ctx.pdf")
volcanoplot2(res = res_hip_cor,  lfcthresh = 1, sigthresh = signif,
             xlim = c(-7,7), ylim = c(y_lower, y_upper), xlab = "log Fold Change", ylab = "-log(Adj. P-Value)", main = "Hippocampus vs Cortex")
dev.off()

setEPS()
postscript("volcano_hip_vs_ctx.eps")
volcanoplot2(res = res_hip_cor,  lfcthresh = 1, sigthresh = signif,
             xlim = c(-7,7), ylim = c(y_lower, y_upper), xlab = "log Fold Change", ylab = "-log(Adj. P-Value)", main = "Hippocampus vs Cortex")
dev.off()


## Striatum v Cortex WT
res_str_cor <- limma::topTable(fit2, coef = 9, adjust.method = "fdr", p.value = 1, number = Inf)
res_str_cor$neglog.P.Val <- -log10(res_str_cor$adj.P.Val)
pdf("volcano_str_vs_ctx.pdf")
volcanoplot2(res = res_str_cor,  lfcthresh = 1, sigthresh = signif,
             xlim = c(-7,7), ylim = c(y_lower, y_upper), xlab = "log Fold Change", ylab = "-log(Adj. P-Value)", main = "Striatum vs Cortex")
dev.off()

setEPS()
postscript("volcano_str_vs_ctx.eps")
volcanoplot2(res = res_str_cor,  lfcthresh = 1, sigthresh = signif,
             xlim = c(-7,7), ylim = c(y_lower, y_upper), xlab = "log Fold Change", ylab = "-log(Adj. P-Value)", main = "Striatum vs Cortex")
dev.off()

## Striatum v Hip WT
res_str_hip <- limma::topTable(fit2, coef = 10, adjust.method = 'fdr', p.value = 1, number = Inf)
res_str_hip$neglog.P.Val <- -log10(res_str_hip$adj.P.Val)
pdf("volcano_str_vs_hip.pdf")
volcanoplot2(res = res_str_hip,  lfcthresh = 1, sigthresh = signif,
             xlim = c(-7,7), ylim = c(y_lower, y_upper), xlab = "log Fold Change", ylab = "-log(Adj. P-Value)", main = "Striatum vs Hippocampus")
dev.off()

setEPS()
postscript("volcano_str_vs_hip.eps")
volcanoplot2(res = res_str_hip,  lfcthresh = 1, sigthresh = signif,
             xlim = c(-7,7), ylim = c(y_lower, y_upper), xlab = "log Fold Change", ylab = "-log(Adj. P-Value)", main = "Striatum vs Hippocampus")
dev.off()

#Cerebral cortex	= CTX
#Hippocampus	= HIP
#Striatum		= STR
#Cerebellum	= CBL

write.csv(x = res_cer_cor, file = "res_DEgenes_cbl_vs_ctx_padj_0p05.csv")
write.csv(x = res_cer_str, file = "res_DEgenes_cbl_vs_str_padj_0p05.csv")
write.csv(x = res_cer_hip, file = "res_DEgenes_cbl_vs_hip_padj_0p05.csv")
write.csv(x = res_hip_cor, file = "res_DEgenes_ctx_vs_hip_padj_0p05.csv")
write.csv(x = res_str_cor, file = "res_DEgenes_str_vs_ctx_padj_0p05.csv")
write.csv(x = res_str_hip, file = "res_DEgenes_str_vs_hip_padj_0p05.csv")

##### eliminating sex as factor
sex <- factor(phenotype$Sex, levels=c("M","F"))
design <- model.matrix(~0 + sex)
colnames(design) <- c("Male","Female")
fit_sex <- lmFit(eset, design=design)
cont.matrix <- makeContrasts(Male-Female, levels=design)
fit_sex_2 <- contrasts.fit(fit_sex, cont.matrix)
fit_sex_2 <- eBayes(fit_sex_2)
limma::topTable(fit_sex_2, coef = 1, adjust='fdr', p.value = 0.05, number = Inf)

## looking for Gonzalez, et. al. genes from rtPCR of Het/WT in STR
gonzalez_genes <- c("Pde10a","Foxp1","Tfap2d","Gnal","Drd2")
my_cols <- c("SYMBOL","logFC","AveExpr","adj.P.Val")
library("tidyverse")
dplyr::filter(res_cer_cor[,my_cols], SYMBOL %in% gonzalez_genes)
dplyr::filter(res_cer_hip[,my_cols], SYMBOL %in% gonzalez_genes)
dplyr::filter(res_cer_str[,my_cols], SYMBOL %in% gonzalez_genes)
dplyr::filter(res_cor_hip[,my_cols], SYMBOL %in% gonzalez_genes)
dplyr::filter(res_cor_str[,my_cols], SYMBOL %in% gonzalez_genes)
dplyr::filter(res_str_hip[,my_cols], SYMBOL %in% gonzalez_genes)

## pathway analysis 

## need to recall the topTables here with desired p-values
res_cer_cor <- limma::topTable(fit2, coef = 5, adjust.method="fdr", p.value = 0.05, number = Inf)
res_cer_str <- limma::topTable(fit2, coef = 6, adjust.method="fdr", p.value = 0.05, number = Inf)
res_cer_hip <- limma::topTable(fit2, coef = 7, adjust.method="fdr", p.value = 0.05, number = Inf)
res_hip_cor <- limma::topTable(fit2, coef = 8, adjust.method = "fdr", p.value =0.05, number = Inf)
res_str_cor <- limma::topTable(fit2, coef = 9, adjust.method = "fdr", p.value = 0.05, number = Inf)
res_str_hip <- limma::topTable(fit2, coef = 10, adjust.method = 'fdr', p.value = 0.05, number = Inf)

res_global <- limma::topTable(fit_global_2, coef = 1, adjust.method="fdr", number=Inf, p.value = 0.01)

## coordinate Illumina IDs to Ensembl IDs for GAGE analysis 
x <- illuminaMousev2ENSEMBL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
z <- unlist(xx)
w <- names(z)
z <- as.vector(z)
illumina2ens <- data.frame(z,w, stringsAsFactors = FALSE)
names(illumina2ens) <- c("EnsemblID","IlluminaID")




library("tidyverse")
res_cer_cor_ens <- inner_join(res_cer_cor, illumina2ens)
res_cer_str_ens <- inner_join(res_cer_str, illumina2ens)
res_cer_hip_ens <- inner_join(res_cer_hip, illumina2ens)
res_hip_cor_ens <- inner_join(res_hip_cor, illumina2ens)
res_str_cor_ens <- inner_join(res_str_cor, illumina2ens)
res_str_hip_ens <- inner_join(res_str_hip, illumina2ens)

#res_global <- rownames_to_column(res_global, var = "IlluminaID")
res_global_ens <- inner_join(res_global, illumina2ens)


library("AnnotationDbi")
library("org.Mm.eg.db")

columns(org.Mm.eg.db)
res_cer_cor_ens$entrez = mapIds(org.Mm.eg.db,
                                keys=res_cer_cor_ens$EnsemblID, 
                                column="ENTREZID",
                                keytype="ENSEMBL",
                                multiVals="first")

res_cer_str_ens$entrez = mapIds(org.Mm.eg.db,
                                keys=res_cer_str_ens$EnsemblID, 
                                column="ENTREZID",
                                keytype="ENSEMBL",
                                multiVals="first")

res_cer_hip_ens$entrez = mapIds(org.Mm.eg.db,
                                keys=res_cer_hip_ens$EnsemblID, 
                                column="ENTREZID",
                                keytype="ENSEMBL",
                                multiVals="first")

res_hip_cor_ens$entrez = mapIds(org.Mm.eg.db,
                                keys=res_hip_cor_ens$EnsemblID, 
                                column="ENTREZID",
                                keytype="ENSEMBL",
                                multiVals="first")

res_str_cor_ens$entrez = mapIds(org.Mm.eg.db,
                                keys=res_str_cor_ens$EnsemblID, 
                                column="ENTREZID",
                                keytype="ENSEMBL",
                                multiVals="first")

res_str_hip_ens$entrez = mapIds(org.Mm.eg.db,
                                keys=res_str_hip_ens$EnsemblID, 
                                column="ENTREZID",
                                keytype="ENSEMBL",
                                multiVals="first")

res_global_ens$entrez = mapIds(org.Mm.eg.db,
                                keys=res_global_ens$EnsemblID, 
                                column="ENTREZID",
                                keytype="ENSEMBL",
                                multiVals="first")


library("gage")
library("pathview")
library("gageData")
data("kegg.sets.mm")
data("sigmet.idx.mm")

## subset to just metabolic pathways and signalling
kegg.sets.mm = keg.sets.mm[sigmet.idx.mm]
head(kegg.sets.mm, 3)

##
fc_cer_cor <- res_cer_cor_ens$logFC
names(fc_cer_cor) <- res_cer_cor_ens$entrez

keggres_cer_cor <- gage(fc_cer_cor, gsets = kegg.sets.mm, same.dir = TRUE)
lapply(keggres_cer_cor, head)

##
fc_cer_str <- res_cer_str_ens$logFC
names(fc_cer_str) <- res_cer_str_ens$entrez

keggres_cer_str <- gage(fc_cer_str, gsets = kegg.sets.mm, same.dir = TRUE) 
lapply(keggres_cer_str, head)

##
fc_cer_hip <- res_cer_hip_ens$logFC
names(fc_cer_hip) <- res_cer_hip_ens$entrez

keggres_cer_hip <- gage(fc_cer_hip,  gsets = kegg.sets.mm, same.dir = TRUE) 
lapply(keggres_cer_hip, head)

##
fc_cor_hip <- res_hip_cor_ens$logFC
names(fc_cor_hip) <- res_hip_cor_ens$entrez

keggres_cor_hip <- gage(fc_cor_hip, gsets = kegg.sets.mm, same.dir = TRUE)
lapply(keggres_cor_hip, head)

##
fc_cor_str <- res_str_cor_ens$logFC
names(fc_cor_str) <- res_str_cor_ens$entrez

keggres_cor_str <- gage(fc_cor_str, gsets = kegg.sets.mm, same.dir = TRUE)
lapply(keggres_cor_str, head)

##
fc_str_hip <- res_str_hip_ens$logFC
names(fc_str_hip) <- res_str_hip_ens$entrez

keggres_str_hip <- gage(fc_str_hip, gsets = kegg.sets.mm, same.dir = TRUE)
lapply(keggres_str_hip, head)


###########

fc_global <- res_global_ens$logFC
names(fc_global) <- res_global_ens$entrez

keggres_global <- gage(fc_global, gsets = kegg.sets.mm, same.dir = TRUE)
lapply(keggres_global, head)
