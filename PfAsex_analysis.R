library(scater, quietly = TRUE)

options(stringsAsFactors = FALSE)

setwd("~/Work/Writing/Papers/scPlasmodium/eLife/Supplementary Data/")

# Read in data
molecules <- read.table("PfAsex_counts.txt", sep = "\t", header = TRUE, row.names=1)
anno <- read.table("PfAsex_meta.txt", sep = "\t", header = TRUE)

genes <- read.table("fal.desc", header=FALSE, row.names=1, quote = "", sep="\t")

# Set up objects
pheno_data <- new("AnnotatedDataFrame", anno)
rownames(pheno_data) <- pheno_data$sample_id
PfAsex <- scater::newSCESet(
  countData = molecules,
  phenoData = pheno_data
)

# Filter genes with no counts
keep_feature <- rowSums(counts(PfAsex) > 0) > 0
PfAsex <- PfAsex[keep_feature, ]

PfAsex <- scater::calculateQCMetrics(PfAsex)

# Filter cells with low counts
filter_by_total_counts <- (PfAsex$total_counts > 25000)
table(filter_by_total_counts)
# Filter cells with low numbers of features detected
filter_by_expr_features <- (PfAsex$total_features > 1000)
table(filter_by_expr_features)

# Filter data
PfAsex$use <- (
  # sufficient features (genes)
  filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # controls shouldn't be used in downstream analysis
    !PfAsex$is_control

)

table(PfAsex$use)

# Gene filtering
# e.g. greater than 10 reads in at least 5 cells
filter_genes <- apply(counts(PfAsex[ , pData(PfAsex)$use]), 1, 
                      function(x) length(x[x >= 10]) >= 5)
table(filter_genes)

fData(PfAsex)$use <- filter_genes

dim(PfAsex[fData(PfAsex)$use, pData(PfAsex)$use])

# Create new object with filtered dataset
PfAsex.qc <- PfAsex[fData(PfAsex)$use, pData(PfAsex)$use]

# scran normalisation (30 and 15 are default values for clustering)
qclust <- scran::quickCluster(PfAsex.qc, min.size = 30)
PfAsex.qc <- scran::computeSumFactors(PfAsex.qc, sizes = 15, clusters = qclust, positive=TRUE)
PfAsex.qc <- scater::normalize(PfAsex.qc)
scater::plotPCA(PfAsex.qc,
                colour_by = "lopez",
                size_by = "total_features",
                exprs_values = "exprs")


