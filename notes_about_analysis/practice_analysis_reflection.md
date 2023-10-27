# Practice Analysis Reflection

Throughout this practice analysis, there were several steps and concepts that initially gave me trouble. As a learning exercise, I will document the technical or conceptual issues that I encountered, how I was able to resolve them, and what resources I used to resolve them. I will also include the remaining questions that I have about the RNA-seq analysis process.

Then, I will compare the plots and results obtained by running the analysis with different log2 fold change thresholds and different shrinkage methods. 

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


## Results Comparison

### DE heatmaps
![lfcThreshold (1)](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/997fb20f-38dc-4853-88b8-a0ae4552ebd3)

### Volcano plots
![lfcThreshold](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/1eef19b6-9592-4dd5-8963-43cf74f8cf3b)

### Unshrunken MA plots
 ![lfcThreshold (2)](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/1d5a18e4-7120-474c-965a-a9c637ee42b4)

### Shrunken MA plots
![lfcThreshold (3)](https://github.com/lenarayneallen/DESeq2_practice/assets/124638335/2aff9c5a-4b30-46b4-99ac-806391893208)

### Top 20 DE genes
**ashr, no threshold**

**ashr, lfcThresh = 0.2**

**normal, no threshold**

**normal, lfcThresh = 0.2**
