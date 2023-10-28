BiocManager::install("GEOquery")
BiocManager::install("Biobase")
BiocManager::install("apeglm")
BiocManager::install("biomaRt")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Mm.eg.db")
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("devtools")
install.packages
devtools::install_github("stephenturner/annotables")
devtools::install_github("r-lib/conflicted")
install.packages("tidyverse")
install.packages('ashr')

library(GEOquery)
library(Biobase)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(annotables)
library(dplyr)
library(tidyverse)
library(biomaRt)
library(conflicted)


#reading in data and metadata
countstable <- read.table('GSE225048_geo_rna_gene_counts.txt', header = TRUE, row.names = 1)
colnames(countstable) <- c('IgG-2', 'IgG-3','IgG-4','IgG-5', 'IgG-6', 'aOX40-1', 'aOX40-2', 'aOX40-3')
metatable <- read.csv('metadata.txt', header = TRUE, row.names = c('aOX40-3','aOX40-2','aOX40-1','IgG-6','IgG-5','IgG-4','IgG-3','IgG-2'))

#reordering metadata to match the order of the counts 
order <- match(colnames(countstable), rownames(metatable))
reordered_metatable <- metatable[order,]

#checking to ensure proper ordering
all(rownames(reordered_metatable) == colnames(countstable))

#making the DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(countstable), colData = reordered_metatable, design = ~treatment)






#quality control
#getting size factors and extracting normalized counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalize = TRUE)

#transforming with vst (DESeq2 logarithmic transformation)
vsd <- vst(dds, blind = TRUE)

#getting pairwise correlation values
vsd_matrix <- assay(vsd)
pairwise_corr_values <- cor(vsd_matrix)

#making heatmap
annot_pallete <- list(treatment = c(aOX40 = "#669966", IgG = "#000000"))
heatmap <- pheatmap(pairwise_corr_values,
                    annotation_col = dplyr::select(metatable, treatment),
                    color = brewer.pal(11, "PuOr"),
                    annotation_colors = annot_pallete)
#making PCA plot
plotPCA(vsd, intgroup = c("treatment"))






#DE analysis
#running analysis of DESeq2 object
dds_analyzed <- DESeq(dds)

#analyzing mean-variance relationship
means_per_gene <- apply(countstable, 1, mean)
variance_per_gene <- apply(countstable, 1, var)
mv_dataframe <- data.frame(means_per_gene,variance_per_gene)

#plotting mean-variance relationship
ggplot(mv_dataframe) +
  geom_point(aes(x = means_per_gene, y = variance_per_gene)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("mean counts per gene") +
  ylab("variance per gene")

#plotting dispersion
plotDispEsts(dds_analyzed,
             fitcol = "#FF9933",
             finalcol = "#9933FF")


#getting wald test results to get log2 fold changes, specifying custom contrasts and threshold
results <- results(dds_analyzed, 
                   contrast = c("treatment", "aOX40", "IgG"), 
                   alpha = 0.05,
                   lfcThreshold = 0.2)

#shrinking log2 foldchange estimates for genes with low counts or high dispersion values
shrunken_results <- lfcShrink(dds_analyzed, 
                              contrast = c("treatment", "aOX40", "IgG"), 
                              res = results, 
                              type = "ashr")
#plotting unshrunken MA plot
plotMA(results)



#plotting shrunken MA plot
plotMA(shrunken_results)

#exploring!!!
head(shrunken_results, n = 10)
summary(shrunken_results)






#annotating

#exploring mouse genome
grcm38

#converting shrunken log2 fold change results to a dataframe
#shifting row names to the first column of the dataframe (for the purpose of later merging)
shrunken_res_df <- data.frame(shrunken_results)
shrunken_res_df_rtc <- rownames_to_column(shrunken_res_df, var = "genename")


#converting gene names to ensembl IDs and format
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#getting mouse gene names column from results dataframe 
mouse_gene_ids <- shrunken_res_df_rtc[,1]

#making table of genenames with corresponding ensembl IDs
ids_table <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
             filters = 'mgi_symbol',
             values = mouse_gene_ids,
             mart = ensembl)

#renaming columns for upcoming merging step
names(ids_table)[1] <- "ensgene"
names(ids_table)[2] <- "genename"

#merging results with genenames/ensembl IDs to add ensembl IDs to the results table
merged_ids_res <- left_join(x = shrunken_res_df_rtc, y = ids_table,
                        by = "genename")
#annotating with selected information from grcm38
results_all <- left_join(x = merged_ids_res, 
                         y = grcm38[, c("ensgene", "symbol", "description")],
                         by = "ensgene")






#exploring results

#getting significantly differentially expressed genes
sig_res <- subset(results_all, padj < 0.05)
#arranging by significance
sig_res <- sig_res %>%
  arrange(padj)

#making heatmap of the normalized counts of significantly differentially expressed genes
#extracting normalized counts of significantly diferentially expressed genes
sig_norm_counts <- normalized_counts[sig_res$genename,]

#plotting heatmap 
pheatmap(sig_norm_counts,
         color = brewer.pal(6, "PuOr"),
         cluster_rows = TRUE,
         show_rownames = FALSE, 
         annotation_col = dplyr::select(metatable, treatment),
         annotation_colors = annot_pallete,
         scale = "row")

#making volcano plot of log2 fold changes vs. p-value (significance)
#creating a column of logical vectors indicating significance
results_all_logical <- results_all %>%
  mutate(threshold = padj < 0.05)

#plotting volcano plot
ggplot(results_all_logical) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj),
                 color = threshold)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  scale_color_manual(name = "threshold",
                     values = c("TRUE" = "#669966", "FALSE" = "#9933FF"))
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))


#making expression plot of top 20 differentially expressed genes
#converting significant normalized counts matrix to a dataframe and subsetting top 20
top_20_sig_norm <- data.frame(sig_norm_counts)[1:20,] 

#fixing the column naming error introduced in the matrix to dataframe conversion
colnames(top_20_sig_norm) <- c('IgG-2', 'IgG-3','IgG-4','IgG-5', 'IgG-6', 'aOX40-1', 'aOX40-2', 'aOX40-3')

#reformatting and gathering
top_20_sig_norm <- rownames_to_column(top_20_sig_norm, var = "symbol")
top_20_sig_norm <- gather(top_20_sig_norm, key = "sample", value = "normalized_counts", 2:9)


top_20_sig_norm <- inner_join(top_20_sig_norm, rownames_to_column(metatable, var = "sample"), 
                     by = "sample")
#plotting DE plot of top 20 differentially expressed genes 
ggplot(top_20_sig_norm) +
  geom_point(aes(x = symbol, y = normalized_counts, color = treatment)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("Normalized Counts") +
  ggtitle("20 Most Significant DE Genes") +
  scale_color_manual(values = c("aOX40" = "#FF9933", "IgG" = "#9933FF")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) 
  
