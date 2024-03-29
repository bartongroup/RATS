% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func.R
\name{calculate_DTU}
\alias{calculate_DTU}
\title{Set up and execute the tests.}
\usage{
calculate_DTU(
  counts_A,
  counts_B,
  tx_filter,
  test_transc,
  test_genes,
  full,
  count_thresh,
  p_thresh,
  dprop_thresh,
  correction,
  threads = 1,
  use_sums = FALSE
)
}
\arguments{
\item{counts_A}{A data.table of counts for condition A. x: sample, y: transcript.}

\item{counts_B}{A data.table of counts for condition B. x: sample, y: transcript.}

\item{tx_filter}{A data.table with target_id and parent_id. Pre-processed with \code{tidy_annot()}.}

\item{test_transc}{Whether to do transcript-level test.}

\item{test_genes}{Whether to do gene-level test.}

\item{full}{Either "full" (for complete output structure) or "short" (for bootstrapping).}

\item{count_thresh}{Minimum average count across replicates.}

\item{p_thresh}{The p-value threshold.}

\item{dprop_thresh}{Minimum difference in proportions.}

\item{correction}{Multiple testing correction type.}

\item{threads}{Number of threads (POSIX systems only).}

\item{use_sums}{Use sums instead of means. (sums is the behaviour of RATs up to 0.6.5 inclusive, more sensitivie but more prone to false positives)}
}
\value{
list
}
\description{
Set up and execute the tests.
}
