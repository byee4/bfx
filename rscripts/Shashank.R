library(DESeq2)

countData <- read.csv( "/projects/ps-yeolab/shsathe/Kris_TIA1_ASOs/DESeq2/Hur_Spinal_Cord/hur_spinal_cord_counts.txt", header=TRUE, row.names=1, sep='\t')
colData <- read.csv("/projects/ps-yeolab/shsathe/Kris_TIA1_ASOs/DESeq2/Hur_Spinal_Cord/DESeq2_manifest.txt", header=TRUE,row.names=1, sep='\t')
countTable <- as.matrix(countData)
storage.mode(countTable) = 'integer'
rs <- rowMeans(countTable)
use <- (rs > 10)
countTableFilt <- countTable[ use, ]
dim(countTableFilt)
dds <- DESeqDataSetFromMatrix(countData = countTableFilt,colData = colData,design = ~ Condition)
dds <- DESeq(dds, betaPrior=FALSE)
normalized_counts <- counts(dds, normalized=TRUE)

wt_hur <- results(dds, contrast=c("Condition","wt","hur"), independentFiltering = FALSE)
write.csv(as.data.frame(wt_hur),file="wt_hur.csv")
write.csv(as.data.frame(normalized_counts),file="wt_hur_normalized_counts.csv")
