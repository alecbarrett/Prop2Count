## load all CeNGEN single cell counts and proportions, calculate conformable matrices, and then fit a logit function to move from proportions to TPM space

library(Seurat)
library(pbapply)
library(tidyverse)
library(dplyr)
library(stringr)


## load the cengen single cell dataset

sc_object <- readRDS('100720_L4_all_cells_Seurat.rds')


sc_object <- sc_object[,sc_object$Tissue != 'Unknown' & sc_object$Tissue != 'Unannotated']

## merge clusters down to anatomical cell types
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('DB01')] <- 'DB'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('DA9')] <- 'DA'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('VC_4_5')] <- 'VC'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('VB01', 'VB02')] <- 'VB'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('RMD_DV', 'RMD_LR')] <- 'RMD'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('RME_DV', 'RME_LR')] <- 'RME'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('VA12')] <- 'VA'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('IL2_DV', 'IL2_LR')] <- 'IL2'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('AWC_ON', 'AWC_OFF')] <- 'AWC'
sc_object <- sc_object[,!(sc_object$Cell.type %in% c('RIV_stressed', 'SMD_stressed'))]

non_neuronal_list <- c('Excretory', 'Glia', 'Hypodermis', 'Intestine', 'Muscle_mesoderm', 'Pharynx', 'Reproductive')

sc_object$neuron_level <- sc_object$Tissue

sc_object$neuron_level[sc_object$neuron_level=='Neuron'] <- sc_object$Cell.type[sc_object$neuron_level=='Neuron']


### cut down the size of the Seurat object
sc_object_cut <- CreateSeuratObject(counts = sc_object@assays$RNA@counts, 
                                    meta.data = sc_object@meta.data[,c('neuron_level', 'Cell.type', 'orig.ident', 'Barcode', 'Tissue', 'Experiment')])


rm(sc_object)



prop_by_type <- pbsapply(unique(sc_object_cut$neuron_level), function(cell_type){
  temp <- sc_object_cut@assays$RNA@counts[,sc_object_cut$neuron_level==cell_type]
  rowSums(temp > 0)/ncol(temp)
})
prop_by_type <- prop_by_type[order(rownames(prop_by_type)), order(colnames(prop_by_type))]

CeNGEN_TPM <- pbsapply(unique(sc_object_cut$neuron_level), function(cell_type){
  temp <- sc_object_cut@assays$RNA@counts[,sc_object_cut$neuron_level==cell_type]
  temp <- sweep(temp, 2, colSums(temp), '/') * 1000000
  rowMeans(temp)
})
CeNGEN_TPM <- CeNGEN_TPM[order(rownames(CeNGEN_TPM)), order(colnames(CeNGEN_TPM))]


## pair up the TPM and proportions values
all_values_matched.df <- data.frame(
  TPM = reshape::melt(CeNGEN_TPM[rownames(prop_by_type), colnames(prop_by_type)])[,'value'],
  proportion = reshape::melt(prop_by_type)[,'value'])


### fit a logit function
fit <- nls(log(TPM) ~ ((a*log(proportion/(1-proportion))) + b),
            data = all_values_matched.df[all_values_matched$proportion > 0 & all_values_matched$proportion < 1,], 
            start = list(a = 0.9, b = 13))

summary(fit)

## get the parameters
a <- fit$m$getPars()['a']
b <- fit$m$getPars()['b']


## get the predicted values, and remove infinite values created by proportion == 1
predicted_TPM_values <- exp( ((a*log(all_values_matched.df$proportion/(1-all_values_matched.df$proportion))) + b) )
predicted_TPM_values[is.infinite(predicted_TPM_values)] <- max(predicted_TPM_values[is.finite(predicted_TPM_values)])

## plot the fit
plot( all_values_matched.df$proportion[1:100000], log(all_values_matched.df$TPM[1:100000]) )
lines( log(all_values_matched.df$TPM[1:20000])[order(all_values_matched.df$TPM[1:20000])],
      predicted_TPM_values[1:20000][order(all_values_matched.df$TPM[1:20000])], col = 'red')


