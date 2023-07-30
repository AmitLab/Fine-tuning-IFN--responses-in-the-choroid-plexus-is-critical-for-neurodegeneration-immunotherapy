
library("BiocParallel")
library(DESeq2)
library(pheatmap) 
library(apeglm)
library(plyr)
library(pheatmap)
library(ggplot2)
library(stringr)

##I) Loading/annotation/QC of the  data 
count_table=read.delim("20200903_exons_raw.txt",header = T,row.names = 1)
coldata=read.delim("samples.txt",header=T)

#convert gene strings to names
gene_names=rownames(count_table)
x = strsplit(gene_names, "\\|")
x = unlist(lapply(x,FUN = function(x) {x[1]}))
rownames(count_table)=x
colnames(count_table) = as.character(t(coldata$Condition))

#Annotation
annotation_table=data.frame(coldata)
count_table=count_table[,colnames(count_table) %in% c("IFNg","UT")]
new_count_table=count_table
new_annotation_table=annotation_table

##II)Analysis using DESeq
new_dds <- DESeqDataSetFromMatrix(countData = new_count_table, colData = new_annotation_table, design = ~ Condition)
hist(log10(1+rowSums(new_count_table)),n=100)
keep <- rowSums(counts(new_dds)) >= 10
new_dds <- new_dds[keep,]

#define reference level
new_dds$Condition <- relevel(new_dds$Condition, ref = "UT")
#differential expression analysis
new_dds <- DESeq(new_dds)
plotDispEsts(new_dds)
plotMA(new_dds, ylim=c(-3,3),alpha = 0.01)
new_res <- results(new_dds)
resultsNames(new_dds)
#LFC shrinkage
new_resLFC <- lfcShrink(new_dds, coef="Condition_IFNg_vs_UT", type="apeglm")
plotMA(new_resLFC, ylim=c(-6,6),alpha = 0.01)
new_resLFC[order(new_resLFC$log2FoldChange, decreasing = T),]
new_table_LFC=as.matrix(new_resLFC)
write.table(new_table_LFC, "20230724_results.txt", sep="\t")

#VST transformation
new_vsd <- varianceStabilizingTransformation(new_dds,blind = T)
new_Normalized_table = assay(new_vsd)
new_Normalized_variance = apply(new_Normalized_table,MARGIN = 1,FUN = var)
new_Normalized_variance =new_Normalized_variance[order(new_Normalized_variance,decreasing = T)]
new_Top_2000_genes = names(new_Normalized_variance[1:2000])
new_Top_500_genes = names(new_Normalized_variance[1:500])

#sample to sample distance
sampleDists <- dist(t(assay(new_vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(new_vsd$Condition, new_vsd$experiment, sep="-")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         show_rownames = T)
#plot PCA
plotPCA(new_vsd, intgroup=c("Condition", "experiment"))
