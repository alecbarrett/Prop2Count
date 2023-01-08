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
### Figure 3: modeling the %non-zero cells against the ln(counts) in a single replicate using a logit function ln(counts) = ln(p/1-p) + ln(nCells). Each dot is a gene, and we are representing the data from just one pseudobulk sample

From this graph, we see two things
1. While most genes show a clear tight relationship between the %non-zero cells and the counts, some genes that are detected in vanishingly few cells have extremely high counts values
    * This may represent rare undetected doublets, high background in certain cells, or misannotated individual barcodes
    * We propose that this is the cause of the lower accuracy seen in the counts data for gene detection, and that the %non-zero cells is a metric that is more robust to rare but strong injections of contaminating information.
2. We can model the relationship between counts and the %non-zero cells with a simple logit function, without any need to estimate parameters for the fit. This in turn allows us a simple equation for moving between counts and proportions spaces.
    * This allows us to easily move between counts and proportions spaces without the need for modeling dataset specific parameters



## The Transformation

1. Take the %non-zero cells for each pseudobulk replicate in the single cell dataset, and calculate the %non-zero cells, and the total number of cells.
2. Calculate the equivalent counts estimate as counts = (P/(1-P)) * nCells

This simple transformation allows us to quickly obtain a denoised count value for each gene x replicate, which is then usable in most tools for differential expression analysis.

## Improvements in differential expression accuracy

![prop2count Differential Expression](https://github.com/alecbarrett/Prop2Count/blob/main/img/prop2count%20dex.png)
### Figure 4: prop2count improves differential expression accuracy by reducing false positives.

Here we show a boxplot of the true positive rate (TPR), false positive rate (FPR) and Matthew's correlation coefficient (MCC) for differential expression of ground truth genes between different neurons in C. elegans. Each data point is the directional score for a neuron-pair (example: AFD vs ADL, and ADL vs AFD are each a datapoint).

p-values were obtained using a Wilcoxon Rank-Sum test.

By comparing the accuracy of differential expression calling using three metrics, we see that prop2counts transformed data are approximately as good as the counts data in detecting expected events, but have a 45% drop in the median false positive rate, and a corresponding increase in the overal correspondance as shown by the MCC score.

## Takeaway

prop2count is a simple, fast, parameter-free model that allows users to obtain denoised counts for pseudobulk single cell RNA-seq datasets, improving both gene detection, and differential expression accuracy.
