# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Wed Feb 5 04:14:01 EST 2020

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE38012", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6947", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "0000000000000001111111111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=20000)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
#write.table(tT, file="results_CRvsWestern", row.names=F, sep="\t")

library(org.Hs.eg.db)
hs <- org.Hs.eg.db
ids <- select(hs, keys = tT$Gene.symbol, columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")

#x3 <- scan("/Users/myrthedehaan/Downloads/Internship_LUMC/test_list_genes1.0.txt", what="", sep="\n", skip = 1)
z3 <- clusterProfiler::enrichGO(ids$ENTREZID, 'org.Hs.eg.db', ont="MF", pvalueCutoff=0.05, pAdjustMethod="BY")
head(z3) 
res4 <- as.data.frame(z3)
order.p_adjust3 <- order(res4$p.adjust)
res5 <- res4[order.p_adjust3,]
res5$rank <- rank(res5$p.adjust)

z4 <- clusterProfiler::enrichGO(ids$ENTREZID, 'org.Hs.eg.db', ont="CC", pvalueCutoff=0.05, pAdjustMethod="BY")
res6 <- as.data.frame(z4)
order.p_adjust4 <- order(res6$p.adjust)
res7 <- res6[order.p_adjust4,]
res7$rank <- rank(res7$p.adjust)

z5 <- clusterProfiler::enrichGO(ids$ENTREZID, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.05, pAdjustMethod="BY")
res8 <- as.data.frame(z5)
order.p_adjust3 <- order(res8$p.adjust)
res9 <- res8[order.p_adjust3,]
res9$rank <- rank(res9$p.adjust)

clusterProfiler::dotplot(z5, showCategory=30)

################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE38012", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6947", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- "0000000000000001111111111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("CR","Western+diet")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE38012", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")