suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("DESeq2"))

parser <- ArgumentParser()

parser$add_argument("--count", type="character")
parser$add_argument("--conditions", type="character")
parser$add_argument("--column", type="character")
parser$add_argument("--ref", type="character")
parser$add_argument("--filterAvg", type="integer")
parser$add_argument("--output", type="character")
parser$add_argument("--colSkip", type="integer", default=6)

args <- parser$parse_args()

countData <- read.table(args$count,  comment.char="#", header=TRUE, row.names=1, sep='\t')
countData <- countData[, args$colSkip:ncol(countData)]
colData <- read.table( args$conditions, header=TRUE,row.names=1, sep='\t')
(colData)
some_condition <- args$column
ref <- args$ref
alt <- as.character(
  unique(
    colData[,some_condition][which(colData[[some_condition]]!=ref)]
  )
)

countTable <- as.matrix(countData)
storage.mode(countTable) = 'integer'
rs <- rowMeans(countTable)

print("R likes to change values of the header column sometimes, make sure these are identical")
colnames(countTable)
rownames(colData)

use <- (rs > args$filterAvg)
countTableFilt <- countTable[use, ]

dds <- DESeqDataSetFromMatrix(
  countData = countTableFilt,colData = colData,design = formula(
    paste("~",some_condition)
  )
)

dds <- DESeq(dds, betaPrior=FALSE)
normalized_counts <- counts(dds, normalized=TRUE)

diffexp = results(dds, contrast=c(some_condition, alt, ref), independentFiltering=FALSE)
write.csv(as.data.frame(diffexp),file=args$output)
write.csv(as.data.frame(normalized_counts),file=paste0(args$output,".norm_counts"))
