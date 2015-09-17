
# Load the count matrix, normalise using TMM and then output RPKM values.
# https://www.biostars.org/p/99310/
# Author: Thomas E. Gorochowski <teg@mit.edu>, Voigt Lab, MIT.

# We use edgeR for between sample normalisation
library(edgeR)

# Load the data file
data_counts = read.table("./analysis/counts/count.matrix.tsv", header=T, row.names=1, com='')
col_ordering = c(1,2,3,4,5,6,7,8)
count_matrix = data_counts[,col_ordering]
conditions = factor(c(rep("WT", 2), rep("delta", 2), rep("deltaAll", 2), rep("mAll", 2)))

# Load the gene lengths (for RPKM calculation) - reorder
data_gene_length = read.table("./analysis/counts/gene_length.tsv", header=T, row.names=1, com='')
gene_length <- data_gene_length[order(match(rownames(data_gene_length),rownames(data_counts))),1]

# Get the library sizes (total counts in annotated genes)
lib_sizes <- as.vector(colSums(count_matrix))

# Calculate normalisation factors using TMM
expr <- DGEList(counts=count_matrix,group=conditions,lib.size=lib_sizes)
expr <- calcNormFactors(expr)

# Calculate RPKM value (normalisation factors are included by default) log=FALSE
expr_norm <- rpkm(expr,gene.length=gene_length)

# Write the RPKM values to a text file
write.table(expr_norm,col.names=NA,quote=FALSE,file="./analysis/rpkms/rpkm.matrix.tsv",sep="\t")

# Write normalisation factors used
write.table(data.frame("design"=colnames(count_matrix), "factor"=expr$samples$norm.factors, "lib_size"=expr$samples$lib.size), file="./analysis/rpkms/norm_factors.tsv", sep='\t', quote=FALSE, row.names=FALSE)

