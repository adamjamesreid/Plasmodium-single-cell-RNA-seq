library(scater, quietly = TRUE)

options(stringsAsFactors = FALSE)

setwd("~/Work/Writing/Papers/scPlasmodium/eLife/Supplementary Data/")

# Read in data
molecules <- read.table("PfGam_counts.txt", sep = "\t", header = TRUE, row.names=1)
anno <- read.table("PfGam_meta.txt", sep = "\t", header = TRUE)

genes <- read.table("fal.desc", header=FALSE, row.names=1, quote = "", sep="\t")

# Set up objects
pheno_data <- new("AnnotatedDataFrame", anno)
rownames(pheno_data) <- pheno_data$sample_id
PfGam <- scater::newSCESet(
  countData = molecules,
  phenoData = pheno_data
)

# Filter genes with no counts
keep_feature <- rowSums(counts(PfGam) > 0) > 0
PfGam <- PfGam[keep_feature, ]

PfGam <- scater::calculateQCMetrics(PfGam)

# Filter cells with low counts
filter_by_total_counts <- (PfGam$total_counts > 25000)
table(filter_by_total_counts)
# Filter cells with low numbers of features detected
filter_by_expr_features <- (PfGam$total_features > 1000)
table(filter_by_expr_features)

# Filter data
PfGam$use <- (
  # sufficient features (genes)
  filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # controls shouldn't be used in downstream analysis
    !PfGam$is_control
)

table(PfGam$use)

# Gene filtering
# e.g. greater than 10 reads in at least 5 cells
filter_genes <- apply(counts(PfGam[ , pData(PfGam)$use]), 1, 
                      function(x) length(x[x >= 10]) >= 5)
table(filter_genes)

fData(PfGam)$use <- filter_genes

dim(PfGam[fData(PfGam)$use, pData(PfGam)$use])


PfGam.qc <- PfGam[fData(PfGam)$use, pData(PfGam)$use]
#GDV1 PF3D7_0935400
#AP2-G PF3D7_1222600
# NEK3 PF3D7_1201600

# scran normalisation (30 and 15 are default values for clustering)
qclust <- scran::quickCluster(PfGam.qc, min.size = 30)
PfGam.qc <- scran::computeSumFactors(PfGam.qc, sizes = 15, clusters = qclust, positive=TRUE)
PfGam.qc <- scater::normalize(PfGam.qc)
scater::plotPCA(PfGam.qc,
                colour_by = "lasonder",
                size_by = "total_features",
                exprs_values = "exprs")

plotPCA(PfGam.qc,  colour_by="lasonder", exprs_values = "exprs", ncomponents=2)