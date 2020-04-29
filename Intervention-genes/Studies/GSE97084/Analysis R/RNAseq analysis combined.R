# Import libraries
library(GEOquery)
library(limma)
library(edgeR)

# Use GEOquery to get covartiates + preprocess phenodata 
gset <- getGEO("GSE97084", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
phd <- phenoData(gset)
phdf <- phd@data
id <- phdf$title
id <- gsub("[A-Za-z ]", "", id)
phdf$id <- id
colnames(phdf)[39] <- "timepoint"
colnames(phdf)[38] <- "age"

# Get the counttables and preprocess them to be compatible for usage of function. 
counts1 <- read.delim(file = "/Users/myrthedehaan/Downloads/Internship_LUMC/Intervention-genes/Studies/GSE97084/Data/GSE97084_GeneCount_raw.tsv", sep = '\t', header = TRUE)
counts2 <- read.delim(file = "/Users/myrthedehaan/Downloads/Internship_LUMC/Intervention-genes/Studies/GSE97084/Data/GSE97084_GeneCount_raw_2.tsv", sep = '\t', header = TRUE)
rownames(counts1) <- counts1$GeneID
rownames(counts2) <- counts2$GeneID 
counts1 <- counts1[, -c(1:6)]
counts2 <- counts2[, -c(1:6)]
countsmain <- cbind(counts1, counts2)
countsmain <- t(countsmain)
countsmain <- countsmain[-nrow(countsmain),]
row.names(phdf) <- row.names(countsmain)

# Test a certain group (in this case the group that dit two kinds of exercises combined)
keep <- c("Combined")
phdf <- phdf[phdf$`exercise type:ch1` %in% keep, ]
counts <- countsmain[rownames(countsmain) %in% rownames(phdf), ]
x <- c("timepointPostTraining-timepointPreTraining")
results = limma.rnaseq(~ 0 + timepoint + age, phdf, counts, contrasts = x, random_effect="id")

#Here we do a gene set enrichment
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
ids <- select(hs, keys = results$gene, columns = c("ENSEMBL", "ENTREZID"), keytype = "ENSEMBL")

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
