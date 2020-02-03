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
limma.rnaseq(~ 0 + timepoint + age, phdf, counts, contrasts = x, random_effect="id")

