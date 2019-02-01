## Rscript drimseq.R <wt_counts.tsv> <ipf_counts.tsv> <transcript2gene.tsv> <outpref>
args = commandArgs(trailingOnly=TRUE)

library(DRIMSeq)
pv <- packageVersion('DRIMSeq')

basedir <- '~/PROJECTS/rats'
infile1 <- args[[1]]
infile2 <- args[[2]]
t2gfile <- args[[3]]
outpref <- args[[4]]


## Covariates.

covariates <- data.frame(sample_id=c('SRR393010', 'SRR393011', 'SRR393012', 'SRR393013', 'SRR393014', 'SRR393015'),
						 group=c('control', 'control', 'control', 'IPF', 'IPF', 'IPF'))

## Counts.

# Import counts from the tables formatted for SUPPA2 input. Skip header line, as it is not helpful for merging and in fact makes things worse. Rename the cols as required instead.
a <- read.delim(file.path(basedir,infile1), header=FALSE, skip=1, col.names=c('feature_id', as.vector(covariates$sample_id[c(1,2,3)])) )
b <- read.delim(file.path(basedir,infile2), header=FALSE, skip=1, col.names=c('feature_id', as.vector(covariates$sample_id[c(4,5,6)])) )

# Import corresponding gene IDs.
t2g <- read.delim(file.path(basedir, t2gfile), header=TRUE, col.names=c('gene_id','feature_id'))

# Collate.
counts <- merge(merge(t2g, a, by='feature_id'), b, by='feature_id')

# Pack up into object.
if (pv >= '1.3.3') {
	ddat <- dmDSdata(counts=counts, samples=covariates)
} else {
	ddat <- dmDSdata(counts=counts[-c(1,2)],
					 gene_id=counts$gene_id,
					 feature_id=counts$feature_id,
					 sample_id=covariates$sample_id,
					 group=covariates$group)
}

## DTU.

# Filter out low expression/representation (or in this case don't since SUPPA2 doesn't).
repspercond <- table(samples(ddat)$group)
ddat <- dmFilter(ddat,
				 min_samps_gene_expr=sum(repspercond),    # all samples of one condition should meet minimal gene expression
				 min_samps_feature_expr=min(repspercond), # all samples of one condition should meet minimal feature expression
				 min_gene_expr=0,
				 min_feature_expr=0,
				 min_samps_feature_prop=0)

# Model.
if (pv >= '1.3.3') {
	design_full <- model.matrix(~ group, data = samples(ddat), verbose=1)
	#set.seed(123)
	ddat <- dmPrecision(ddat, design = design_full, verbose=1)
	ddat <- dmFit(ddat, design = design_full, verbose=1)
} else {
	ddat <- dmDispersion(ddat, verbose=1)
	ddat <- dmFit(ddat, verbose=1)
}

# confounders:
# poorly documented
# make a design_null in the same way as above?

# Test.
if (pv >= '1.3.3') {
	ddat <- dmTest(ddat, coef = "groupIPF", verbose = 1)
} else {
	ddat <- dmTest(ddat, verbose = 1)
}
##  Results.

# Raw.
saveRDS(ddat, file=file.path(basedir, paste0(outpref,'.RDS')))
res <- results(ddat)
write.table(res, file=file.path(basedir, paste0(outpref, '.results.tsv')), sep="\t", append=FALSE, quote=FALSE, row.names=FALSE)

# adj_pvalue < 0.05
# won't do here. SUPPA2 doesn't give classifications either, so I can do that later in the collective comparison of the tools.
