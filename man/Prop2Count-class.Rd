% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Prop2Count.R
\docType{class}
\name{Prop2Count-class}
\alias{Prop2Count-class}
\alias{Prop2Count}
\title{The Prop2Count Class}
\description{
The Prop2Count Class
}
\section{Slots}{

\describe{
\item{\code{aggregated_counts}}{genes x replicates matrix of integer counts}

\item{\code{proportions}}{genes x replicates matrix of the proportion of cells with > 0 counts}

\item{\code{nCells}}{vector of the number of cells per replicate}

\item{\code{prop2count}}{genes x replicates matrix of transformed proportions as doubles, or integers (if round == T)}
}}

