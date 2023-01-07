### modeling a logit function for every pseudobulk replicate individually to use their proportions to generate a counts profiles


## Libraries ----
library(pbapply)



# cells_per_rep: a vector of how many cells are in each replicate.
# prop_per_rep: a genes x cell_type-replicates proportions matrix (proportion of cells with >=1 count)
# CeNGEN_aggr_bio_reps: a genes x cell_type-replicates counts matrix

## load data ----
cells_per_rep <- readRDS('')
prop_per_rep <- read.table('', sep = '\t')

#### pseudocode ----
# For each replicate, model the proportions - log(counts) relationship as a logit fit
# Use the fit to transform proportions to counts space
# Any instance of 100% proportions will be given a value of infinity, so we replace those values with one fewer than the maximum cell number

fitted_prop_to_counts_per_rep <- pbsapply(colnames(CeNGEN_aggr_bio_reps), function(sam){
  size = cells_per_rep[sam]
  counts = CeNGEN_aggr_bio_reps[,sam]
  proportions = prop_per_rep[,sam]
  tmp.df <- data.frame(counts = counts,
                       proportions = proportions)
  fit <- nls(log(counts) ~ ((a*log(proportions/(1-proportions))) + b),
             data = tmp.df[tmp.df$proportions<1 & tmp.df$proportions > 0,],
             start = list(a = 1, b = 10))
  a <- fit$m$getPars()['a']
  b <- fit$m$getPars()['b']
  #print(sam)
  #print(a)
  #print(b)
  counts_from_prop = exp( (a * log(tmp.df$proportions/(1-tmp.df$proportions))) + b)
  names(counts_from_prop) <- rownames(tmp.df)
  
  prop_max <- (size-1)/size
  #print(prop_max)
  
  counts_from_prop[is.infinite(counts_from_prop)] <- exp( (a * log(prop_max/(1-prop_max))) + b)
  return(counts_from_prop)
})

write.table(fitted_prop_to_counts_per_rep, '.tsv', sep = '\t'))
