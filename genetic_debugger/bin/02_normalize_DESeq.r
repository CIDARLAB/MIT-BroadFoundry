
# Process the unique_counts.txt matrix in a given folder and output factors to 
# same location
args<-commandArgs(TRUE)
unique_count_dir<-args[1]
unique_count_mat<-paste0(unique_count_dir,"counts.txt")

# Let's use the standard DE library
library(DESeq)

###############################################################################
# COMPARISON OF ALL WORKING VS BROKEN STATES
###############################################################################

# Load the data
rnaseqMatrix = read.table(unique_count_mat, header=T, row.names=1, com='')
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=2,]
conditions = factor(c(rep("data", length(rnaseqMatrix))))

# Calculate the size factors
exp_study = newCountDataSet(rnaseqMatrix, conditions)
exp_study = estimateSizeFactors(exp_study)

# Save the size factors to file
write.table(data.frame("design"=colnames(rnaseqMatrix), "factor"=exp_study$sizeFactor), file=paste0(unique_count_dir,"correction_factors_deseq.txt"), sep='\t', quote=FALSE, row.names=FALSE)
