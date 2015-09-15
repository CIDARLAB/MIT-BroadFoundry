
# Process the unique_counts.txt matrix in a given folder and output factors to 
# same location
args<-commandArgs(TRUE)
unique_count_dir<-args[1]
unique_count_mat<-paste0(unique_count_dir,"counts.txt")

# Let's use the standard DE library
library(edgeR)

###############################################################################
# COMPARISON OF ALL WORKING VS BROKEN STATES
###############################################################################

# Load the data
rnaseqMatrix = read.table(unique_count_mat, header=T, row.names=1, com='')
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=2,]
conditions = factor(c(rep("data", length(rnaseqMatrix))))

cds <- DGEList( rnaseqMatrix , group = conditions )
cds <- calcNormFactors( cds )

#a <- cds$samples$lib.size*cds$samples$norm.factors
# Calculate the size factors
#de <- DGEList(counts=rnaseqMatrix)
#de1 <- calcNormFactors(de,method="TMM")
#de$samples$norm.factors <- de1$samples$norm.factors
#exp_study<-calcNormFactors(rnaseqMatrix)

# Save the size factors to file
# cds$samples$lib.size*cds$samples$norm.factors
write.table(data.frame("design"=colnames(rnaseqMatrix), "factor"=cds$samples$norm.factors, "lib_size"=cds$samples$lib.size), file=paste0(unique_count_dir,"correction_factors_edger.txt"), sep='\t', quote=FALSE, row.names=FALSE)
