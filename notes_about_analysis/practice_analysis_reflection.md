# Practice Analysis Reflection

Throughout this practice analysis, there were several steps and concepts that initially gave me trouble. As a learning exercise, I will document the technical or conceptual issues that I encountered, how I was able to resolve them, and what resources I used to resolve them. Then, I will compare the plots and results obtained by running the analysis with different log2 fold change thresholds and different shrinkage methods. Lastly, I will discuss how my results compared to those of the study's authors. 

## Challenging Concepts
### Mean, Variance, and Dispersion

At first, I did not understand the purpose of plotting the mean-variance relationship. To understand this, I refreshed my knowledge on the DESeq2 model and the concept variance. 
Variance is calculated by taking the square root of the standard deviation. For a given gene, the variance describes how much each individual sample’s gene expression differs from the mean gene expression of its treatment group. For RNA expression data, it is typical that the variance increases with mean expression, so plotting the mean-variance relationship is a good way to make sure that the expression data behaves as expected. The greater spread of variances observed in genes with lower mean expression is described as heteroscedasticity; this phenomenon is explained by bias introduced in the sample and library preparation steps of RNA sequencing that more strongly affects genes with lower expressions. 

![Image-1](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/eb7d59a6-aaaa-44db-a15f-4a0b24e2f945)

The variance for each given mean is described as the dispersion by DESeq2. This is defined in the DESeq2 dispersion formula:

![Capture](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/bda5d288-738f-4b18-8ff5-e56c22c8bf08)

Based on this formula, dispersion should decrease as mean increases. Using plotDispEsts(), you can see how well your data fits the DESeq2 model. 

![Hirschhorn Disp Ex -1](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/352b9f65-a0ad-493c-a700-a6c97d446c36)

sources: 

https://www.biostars.org/p/480419/#:~:text=So%2C%20in%20brief%2C%20the%20mean,reliably%20quantifying%20lowly%20expressed%20genes
https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html

### lfcThreshold

When extracting the results from the DESeq analysis using the results() function, I was unsure whether or not to specify the lfcThreshold argument. Based on the datacamp and Harvard Chan Bioinformatics Core tutorials that I followed, it appeared to be appropriate for this data to either not specify a threshold, or to use a low threshold- a value less than 0.5. Though I understand that increasing the threshold may increase the chances that the differentially expressed genes returned will be biologically relevant, I still do not fully understand why the best threshold to use would vary from study to study.

sources:

https://www.biostars.org/p/293212/
https://support.bioconductor.org/p/101504/#101507
https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html

### lfcShrink

While I understand shrinkage is necessary to improve the accuracy of dispersion estimates to avoid false positive DE genes, I encountered an issue with the lfcShrink() function itself. 

I first tried to run the code as it was explained to me in the tutorials I utilized:
![cap1](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/6ade5f8a-4800-422f-b881-12361f87902c)

This returned the following error:

![cap2](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/5abad7f9-f73e-489d-91ab-071fee7d589d)

After doing some digging, I discovered that there are three methods of shrinkage to choose from: apeglm, ashr (adaptive shrinkage), and normal. As explained by the above error, apeglm (which appears to be the default method) cannot be used when specifying contrasts. Instead, you must specify a coefficient (coef). Ashr, however, can be used when specifying contrasts. Both apeglm and ashr appear to generally be considered superior to the normal method. DESeq2 will even provide a warning message when using the normal method: 

![image](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/ccb44f8c-ff52-4f69-a396-262f0ba12ead)

sources: 
https://www.biostars.org/p/9509326/


## Plot Comparisons: shrinkage methods and log2 fold change thresholds
In order to better understand (1) the differences between shrinkage methods and (2) the impact of setting a lfcThreshold, I have provided a comparison of plots created under four conditions:
1. Using ashr shrinkage and no lfcThreshold
2. Using ashr shrinkage with an lfcThreshold of 0.2
3. Using normal shrinkage with no lfcThreshold
4. Using normal shrinkage with an lfcThreshold of 0.2
   
### DE heatmaps
![lfcThreshold (1)](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/997fb20f-38dc-4853-88b8-a0ae4552ebd3)

### Volcano plots
![lfcThreshold](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/1eef19b6-9592-4dd5-8963-43cf74f8cf3b)

### Unshrunken MA plots
 ![lfcThreshold (2)](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/1d5a18e4-7120-474c-965a-a9c637ee42b4)

### Shrunken MA plots
![lfcThreshold (3)](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/2aff9c5a-4b30-46b4-99ac-806391893208)

### Top 20 DE genes
![lfcThreshold (4)](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/83e44efe-a410-40cf-995d-8b9e68be7a1f)

## How did my analysis compare to that of the authors?

I initially chose this dataset because I was familiar with the paper (Hirschhorn et al. 2023), and my area of interest is in cancer immunology. However, the authors did not use DESeq2 for their analysis; they instead used the GUI Partek Flow. While Partek Flow _does_ include a wrapper of DESeq2 for differentiall expression analysis, the authors do not indicate in their methods whether this was utilized or not. 

This bulk RNA seq analysis is a very small experiment in the larger context of the study. The authors initially found that an in-vivo mouse melanoma heterogenously expressing the differentiation marker antigen Trp1 was able to be eliminated by a triple combination therapy including cyclophosphamide, adoptive transfer of T cells specific to the Trp1 antigen, and an anti-OX40 antibody. Importantly, the anti-OX40 component of this therapy was essential to this tumor eliminationn.  In an effort to understand which part of the immune system was responsible for the killing of the cells not expressing the Trp1 antigen (antigen escape variants), the authors completed a series of experiments that identified neutrophils as essential to the killing of antigen escape variants in the context of the triple combination therapy. 

The authors then sought to determine if neutrophils from tumors of mice who had recieved the triple combination therapy (cyclophosphamide, Trp1 T-Cells, and anti-OX40) were distinct from those who recieved control IgG in place of the anti-OX40 component of the therapy. At this point, the bulk RNA-seq analysis was completed in order to preliminarily (and roughly) determine whether the anti-OX40 neutrophils exhibited a distinct phenotype, prior to subsequent scRNA-seq. Thus, one figure summarizing this experiment was included in the paper (Fig. 5A): 

![hirsch5A](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/ab7330ce-04ee-44d2-9c95-e0e75eb20ed0)

For comparison, here is the PCA plot, volcano plot, and DE gene counts (using the ashr shrinkage method and lfcThreshold = 0.2) from my analysis:

![PCA_plot](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/ffac8747-4cb4-4711-90a7-9fc8bf62616b)

![volcano_ashr_0 2](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/7ae885e2-c32c-422a-9735-d87a356e9150)

![ashr0 2number](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/a11a3346-eb14-47ad-90c0-7cdfb526d4a2)


Whereas PC1 represents 28.35% of the variation and PC2 represents 15.74% of the variation in the author's PCA plot, PC1 represents 59% of the variation and PC2 represents 13% of the variation in my PCA plot. Additionally, the authors identified 1135 upregulated genes and 343 downregulated genes in the anti-OX40 treated neutrophils compared to the IgG treated neutrophils, my analysis identified 263 upregulated genes and 579 downregulated genes. While these discrepancies may be accounted for by the different analysis methods, I very well could have made an error somewhere in the analysis. However, the trend of transcriptionally distinct neutrophil population observed in the original study is at least somewhat reflected in my analysis.  

sources:
https://www.cell.com/cell/fulltext/S0092-8674(23)00225-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867423002258%3Fshowall%3Dtrue
https://documentation.partek.com/m/mobile.action#page/12943514
