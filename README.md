# Prop2Count
A parameter-free model for denoising single cell aggregated counts (pseudobulk) using the proportion of non-zero cells per gene


Two easy ways to represent a gene's expression in single cell data is to look at the total number of reads that map to the gene in a cluster (counts), or two quantify the number of cells that have any reads mapping to the gene (%non-zero cells)

![ttx-1](https://github.com/alecbarrett/Prop2Count/blob/main/img/Untitled-1.png)
### Figure 1
Here we see the expression of a known AFD cell marker, ttx-1 by aggregate counts on the left, and by %non-zero cells on the right.

We have recently [shown](https://github.com/cengenproject/Thresholding_sc) that using the %non-zero cells to estimate gene expression provides much more accurate results than any counts based representations.

![Precision-recall](https://github.com/alecbarrett/Prop2Count/blob/main/img/precision_recall.png)
### Figure 2
Here we see that thresholding on the %non-zero cells gives higher accuracy for gene expression detection, using a ground truth set of 160 genes with known expression patterns in C. elegans neurons.

