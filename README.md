# Prop2Count
A parameter-free model for denoising single cell aggregated counts (pseudobulk) using the proportion of non-zero cells per gene


Two easy ways to represent a gene's expression in single cell data are:
1. Look at the total number of reads that map to the gene in a cluster (counts)
2. Quantify the number of cells that have any reads mapping to the gene (%non-zero cells)

![ttx-1](https://github.com/alecbarrett/Prop2Count/blob/main/img/Untitled-1.png)
### Figure 1: Expression of the gene ttx-1 in AFD neurons using aggregate counts, and %non-zero cells


We have recently [shown](https://github.com/cengenproject/Thresholding_sc) that using the %non-zero cells to estimate gene expression provides much more accurate results than any counts based representations.

![Precision-recall](https://github.com/alecbarrett/Prop2Count/blob/main/img/precision_recall.png)
### Figure 2: Precision-Recall Curve, gene detection across thresholds for known ground truth genes. 



![logit model](https://github.com/alecbarrett/Prop2Count/blob/main/img/SIA%20plot%20010823.png)
### Figure 3: modeling the %non-zero cells against the ln(counts) in a single replicate using a logit function
### ln(counts) = ln(p/1-p) + ln(nCells)
