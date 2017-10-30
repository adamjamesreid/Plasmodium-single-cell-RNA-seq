library(scater, quietly = TRUE)

options(stringsAsFactors = FALSE)

setwd("~/Work/Writing/Papers/scPlasmodium/eLife/Supplementary Data/")

# Read in data
molecules <- read.table("PbM_counts.txt", sep = "\t", header = TRUE, row.names=1)
anno <- read.table("PbM_meta.txt", sep = "\t", header=TRUE)

genes <- read.table("berg.desc", header=FALSE, row.names=1, quote = "", sep="\t")

# Set up objects
pheno_data <- new("AnnotatedDataFrame", anno)
rownames(pheno_data) <- pheno_data$sample_id
PbM <- scater::newSCESet(
  countData = molecules,
  phenoData = pheno_data
)

# Filter genes with no counts
keep_feature <- rowSums(counts(PbM) > 0) > 0
PbM <- PbM[keep_feature, ]

# calculateQCMetrics
PbM <- scater::calculateQCMetrics(PbM)

# Filter cells with low counts
filter_by_total_counts <- (PbM$total_counts > 25000)
table(filter_by_total_counts)
# Filter cells with low numbers of features detected
filter_by_expr_features <- (PbM$total_features > 1000)
table(filter_by_expr_features)

# filter out control samples
PbM$is_control <- anno$is_control
filter_by_control <- (PbM$is_control == TRUE)
table(filter_by_control)

# Filter data
PbM$use <- (
  # sufficient features (genes)
  filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # controls shouldn't be used in downstream analysis
    !PbM$is_control
)

table(PbM$use)

# Gene filtering 
# e.g. greater than 10 reads in at least 5 cells
filter_genes <- apply(counts(PbM[ , pData(PbM)$use]), 1, function(x) length(x[x >= 10]) >= 5)

table(filter_genes)

fData(PbM)$use <- filter_genes

dim(PbM[fData(PbM)$use, pData(PbM)$use])

# Create new object with filtered dataset
PbM.qc <- PbM[fData(PbM)$use, pData(PbM)$use]

# Perform SCRAN normalisation
qclust <- scran::quickCluster(PbM.qc, min.size = 30)
PbM.qc <- scran::computeSumFactors(PbM.qc, sizes = 20, clusters = qclust, positive=TRUE)
PbM.qc <- scater::normalize(PbM.qc)
scater::plotPCA(PbM.qc,
                colour_by = "consensus",
                size_by = "total_features",
                exprs_values = "exprs")

# PLot expression of p25 female marker on PCA
plotPCA(PbM.qc, colour_by="PBANKA_0515000", exprs_values = "exprs", ncomponents=2)

