#!/usr/bin/env Rscript

#	Copyright (C) 2015 by
#	Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Arguments:
# 1. count_matrix_filename
# 2. group1 1,2,5,6
# 3. group2 3,4,7,8
# 4. library_size_matrix
# 5. output_file_prefix
args <- commandArgs(TRUE)
arg_count_matrix <- args[1]
arg_group1 <- args[2]
arg_group2 <- args[3]
arg_library_size_matrix <- args[4]
arg_output_file_prefix <- args[5]

# We use edgeR for DE analysis
library(edgeR)

# Load the data file
count_matrix <- read.table(arg_count_matrix, header=T, row.names=1, com='')

# TO CHANGE col_ordering <- c(3,4,1,2)
count_matrix <- data_counts[,col_ordering]
conditions <- factor(c(rep("group_1", 2), rep("group_2", 2)))

# Load the actual mapped reads
library_size_matrix <- read.table(arg_library_size_matrix, header=T, row.names=1, com='')

# Get the library sizes (total counts in annotated genes)
lib_sizes <- as.vector(library_size_matrix)

# Calculate normalisation factors using TMM
expr <- DGEList(counts=count_matrix, group=conditions, lib.size=lib_sizes)
expr <- calcNormFactors(expr)
expr <- estimateCommonDisp(expr)
expr <- estimateTagwiseDisp(expr)
de_data <- exactTest(expr)
de_results <- topTags(de_data, n=length(data_counts[,1]))

# write the output to a text file
write.table(as.matrix(de_results$table), col.names=NA, quote=FALSE, file=), file=paste0(arg_output_file_prefix, ".de.analysis"),sep="\t")
