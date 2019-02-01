## Rscript drimseq.R <a_counts.tsv> <b_counts.tsv> <transcript2gene.tsv> <outpref> <samplename/samplename/...> <condition/condition/...>
## Two conditions maximum

args = commandArgs(trailingOnly=TRUE)

library(DRIMSeq)
pv <- packageVersion('DRIMSeq')

basedir <- '~/PROJECTS/rats'
infile1 <- args[[1]]
infile2 <- args[[2]]
t2gfile <- args[[3]]
outpref <- args[[4]]
smplvect <- args[[5]]
condvect <- args[[6]]

# basedir <- '/Volumes/kfroussios/PROJECTS/rats'
# infile1 <- 'salmon_quant/sonesonDm70-A_merged.tsv_f25'
# infile2 <- 'salmon_quant/sonesonDm70-B_merged.tsv_f25'
# t2gfile <- 'genome/Drosophila_melanogaster.BDGP5.70.t2g.tsv'
# outpref <- 'DE/drimseq_soneson_Dm70-salmon'
# smplvect <- 'a1/a2/a3/b4/b5/b6'
# condvect <- 'A/A/A/B/B/B'


## Covariates.

covariates <- data.frame('sample_id'=strsplit(smplvect, '/')[[1]],
						 'group'=strsplit(condvect, '/')[[1]])
aset <- covariates$group == unique(covariates$group)[1]
bset <- covariates$group == unique(covariates$group)[2]


## Counts.

# Import counts from the tables formatted for SUPPA2 input.
# Skip header line, as it is not helpful for merging and in fact makes things worse. Rename the cols as required instead.
a <- read.delim(file.path(basedir,infile1), header=FALSE, skip=1, col.names=c('feature_id', as.character(covariates$sample_id[aset])) ) # un-factor the IDs
b <- read.delim(file.path(basedir,infile2), header=FALSE, skip=1, col.names=c('feature_id', as.character(covariates$sample_id[bset])) )

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
	ddat <- dmTest(ddat, coef = paste0('group',unique(covariates$group)[2]), verbose = 1)
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
