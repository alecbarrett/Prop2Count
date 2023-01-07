## alec barrett
# June 15th 2022

####### modeling replicates all together (proportions & counts)


# goal: take two genes x replicates matrices from single cell data, and model the proportions and counts as a logit function
# input variable is proportions, and we're trying to predict the counts for each gene in each sample.
# additional input!! will be the number of cells in the replicate

# This script takes your through from data loading to visualzing multiple fits
# It shows that the proportions and counts can be modeled as a logit, even without tuning parameters

## libraries ----
library(pbapply)
library(dplyr)
library(wbData)
library(tidyverse)
library(reshape)
library(ggplot)

## load the data ----

setwd('~/Bioinformatics/Modeling SC replicates/')

biorep_counts <- read.table('CeNGEN_aggr_bio_reps_101821.tsv')
biorep_proportions <- read.table('CeNGEN_bioRep_proportions_041522.tsv')
biorep_cell_number <- read.table('CeNGEN_SC_Biorep_size_061422.tsv')

### conform the data ----
rownames(biorep_cell_number) <- gsub(pattern = '-', replacement = '.', x = rownames(biorep_cell_number))
biorep_cell_number <- biorep_cell_number[rownames(biorep_cell_number) %in% colnames(biorep_proportions), ,drop = F]

biorep_counts <- biorep_counts[order(rownames(biorep_counts)), order(colnames(biorep_counts))]
biorep_proportions <- biorep_proportions[order(rownames(biorep_proportions)), order(colnames(biorep_proportions))]
biorep_cell_number <- biorep_cell_number[order(rownames(biorep_cell_number)), ,drop = F]


waldo::compare(rownames(biorep_proportions), rownames(biorep_counts))

## mutate the data ----
# goal here is to create a matrix where each row is a gene in a replicate (so each gene will be represented ~570 times)
# columns are: replicate, counts, proportions, cell_number (how many cells in the replicate)

biorep_counts_melt <- melt(biorep_counts, variable_name = 'replicate')

biorep_proportions_melt <- melt(biorep_proportions, variable_name = 'replicate')



melted_merged <- data.frame(replicate = biorep_counts_melt$replicate,
                            counts = biorep_counts_melt$value,
                            proportions = biorep_proportions_melt$value,
                            cell_number = biorep_cell_number[biorep_counts_melt$replicate,])



## log transform counts, and remove 0 and 1 proportions to eliminate infinite values then calculate the logit ----

melted_merged_log <- melted_merged[melted_merged$proportions > 0 & melted_merged$proportions < 1,]
melted_merged_log$log_counts <- log(melted_merged_log$counts)

melted_merged_log$logit_proportions <- log(melted_merged_log$proportions/(1-melted_merged_log$proportions))


### plot a few replicates colored by cell number to see the pattern ----

sampled_replicates <- unique(melted_merged_log$replicate)[sample(unique(melted_merged_log$replicate), 60)]

melted_merged_log[melted_merged_log$replicate %in% sampled_replicates,] |> arrange(cell_number) |>
  
  ggplot(aes(x = proportions, y = log_counts, color = log10(cell_number))) + 
  geom_point()


## describing the variables & parameters to tune ----


# "a" is the paramter that squishes or expands the S-shape of the curve. if a is small, the range of the S will be small, if a is big, it'll be a long slow curve
# "b" is the paramter that shifts the S-shape's center. If b is 0, the logit S will be centered at 0. 
# "c" is another parameter like "a", but it is specifically tied to the number of cells in the replicate (i.e. how much influence show replicate size have on the shape of the curve?)
# "d" is another parameter that shifts the S-shape's center, but like "c" it is tied to the number of cells in the replicate


## plot the fit for one replicate at a time ----

# here I use just parameters a & b to fit a logit to each sample separately. We can see that we get a pretty good line of best fit for each one

repl <- "SIA__unc.47_2" ##replicate name goes here

fit_individual <- nls(log_counts ~ ((a*log(proportions/(1-proportions))) + b), 
                      data = melted_merged_log[melted_merged_log$replicate==repl,], start = list(a = 1, b = 1))


plot(melted_merged_log[melted_merged_log$replicate == repl,'proportions'],
     melted_merged_log[melted_merged_log$replicate == repl,'log_counts'],main = paste(repl, cell_number))


lines( x = melted_merged_log[melted_merged_log$replicate == repl,'proportions'] |> sort(),
       y = ((fit_individual$m$getPars()['a']*log(melted_merged_log[melted_merged_log$replicate == repl,'proportions']/
                                                   (1-melted_merged_log[melted_merged_log$replicate == repl,'proportions']))) + 
              fit_individual$m$getPars()['b']) |> sort(), col = 'red')



## Fit with all 4 parameters, using all replicates at once ----

fit <- nls(log_counts ~ ((a * logit_proportions) + (c * cell_number * logit_proportions) + b +  (d * cell_number)), 
            data = melted_merged_log, start = list(a = 100, b = 1, c = 1, d = 1))

fit ## we see from printing the fit that a & b are much more highly weighted than c & d
## this suggests that adding the cell numbers isn't doing a ton perhaps.
## Although the number of cells can vary quite a lot, so 0.002 * 1000 is a big contribution, this could be a problem...

repl <- "PVD__eat.4"
line_prop <- seq(0.001,.999, 0.001)
logit_line_prop <- log(line_prop/(1-line_prop))
cell_number <- biorep_cell_number[repl,]
a = fit$m$getPars()['a']
b = fit$m$getPars()['b']
c = fit$m$getPars()['c']
d = fit$m$getPars()['d']


plot(melted_merged_log[melted_merged_log$replicate == repl,'proportions'],
     melted_merged_log[melted_merged_log$replicate == repl,'log_counts'],main = paste(repl, cell_number))
lines( x = line_prop,
       y = ((a * logit_line_prop) + (c * cell_number * logit_line_prop) + b +  (d * cell_number)),
       col = 'red', lwd = 2)

## this is a somewhat okay fit for a small replicate...


repl <- "SIA__unc.47_2"

line_prop <- seq(0.001,.999, 0.001)
logit_line_prop <- log(line_prop/(1-line_prop))
cell_number <- biorep_cell_number[repl,]
a = fit$m$getPars()['a']
b = fit$m$getPars()['b']
c = fit$m$getPars()['c']
d = fit$m$getPars()['d']


plot(melted_merged_log[melted_merged_log$replicate == repl,'proportions'],
     melted_merged_log[melted_merged_log$replicate == repl,'log_counts'],main = paste(repl, cell_number))
lines( x = line_prop,
       y = (a * logit_line_prop) + (c * cell_number * logit_line_prop) + b +  (d * cell_number),
       col = 'red', lwd = 2)


### this is a bad fit for the high cell number replicate...




### Logit fit with log transformed cell_number ----


fit <- nls(log_counts ~ ((a * logit_proportions) + (c * log(cell_number) * logit_proportions) + b +  (d * log(cell_number))), 
           data = melted_merged_log, start = list(a = 100, b = 1, c = 1, d = 1))

fit ## Here we see parameters a & d are effectively doing nothing as they're close to 1, and parameters b & c are being eliminated!

repl <- "PVD__eat.4"
line_prop <- seq(0.001,.999, 0.001)
logit_line_prop <- log(line_prop/(1-line_prop))
cell_number <- biorep_cell_number[repl,]
a = fit$m$getPars()['a']
b = fit$m$getPars()['b']
c = fit$m$getPars()['c']
d = fit$m$getPars()['d']


plot(melted_merged_log[melted_merged_log$replicate == repl,'proportions'],
     melted_merged_log[melted_merged_log$replicate == repl,'log_counts'],main = paste(repl, cell_number))
lines( x = line_prop,
       y = ((a * logit_line_prop) + (c * log(cell_number) * logit_line_prop) + b +  (d * log(cell_number))),
       col = 'red', lwd = 2)

## this is a better fit! but it's hard to tell with this data..


repl <- "SIA__unc.47_2"


line_prop <- seq(0.001,.999, 0.001)
logit_line_prop <- log(line_prop/(1-line_prop))
cell_number <- biorep_cell_number[repl,]
a = fit$m$getPars()['a']
b = fit$m$getPars()['b']
c = fit$m$getPars()['c']
d = fit$m$getPars()['d']


plot(melted_merged_log[melted_merged_log$replicate == repl,'proportions'],
     melted_merged_log[melted_merged_log$replicate == repl,'log_counts'],main = paste(repl, cell_number))
lines( x = line_prop,
       y = (a * logit_line_prop) + (c * log(cell_number) * logit_line_prop) + b +  (d * log(cell_number)),
       col = 'red', lwd = 2)


### this is a great fit for the high cell number replicate as well!!



#### thoughts ----

## if a & d are close to 1, and b & c are close to 0, do we need parameters at all?

repl <- "SIA__unc.47_2"
repl <- "ADF__tph.1_ceh.10"
repl <- "PVD__eat.4"

line_prop <- seq(0.001,.999, 0.001)
logit_line_prop <- log(line_prop/(1-line_prop))
cell_number <- biorep_cell_number[repl,]
a = fit$m$getPars()['a']
b = fit$m$getPars()['b']
c = fit$m$getPars()['c']
d = fit$m$getPars()['d']


plot(melted_merged_log[melted_merged_log$replicate == repl,'proportions'],
     melted_merged_log[melted_merged_log$replicate == repl,'log_counts'],main = paste(repl, cell_number))
lines( x = line_prop,
       y = (logit_line_prop) + (log(cell_number)) , col = 'blue', lwd = 2)
lines( x = line_prop,
       y = (a * logit_line_prop) + (c * log(cell_number) * logit_line_prop) + b +  (d * log(cell_number)),
       col = 'red', lwd = 2)


## glance at any of these fits, or any other replicate and it's clear that even without any tuning parameters it fits very well
# function: log(counts) = logit(proportions) + log(cell_number) ----



