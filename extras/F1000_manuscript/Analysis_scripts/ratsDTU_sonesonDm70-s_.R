args <- commandArgs(trailingOnly = TRUE)
at <- as.numeric(args[1])  # abund_thresh
dt <- as.numeric(args[2])  # dprop_thresh
qt <- as.numeric(args[3])  # qrep_thresh
rt <- as.numeric(args[4])  # rrep_thresh
t <- as.numeric(args[5])  # threads
out <- args[6]  # outfile

library(rats)

t2g <- read.csv("./genome/Drosophila_melanogaster.BDGP5.70.t2g.tsv", sep="\t", header=TRUE)

tpm <- fish4rodents(A_paths=file.path('./salmon_quant', c('Dm_BDGP5.70.1','Dm_BDGP5.70.2','Dm_BDGP5.70.3')),
					B_paths=file.path('./salmon_quant', c('Dm_BDGP5.70.4','Dm_BDGP5.70.5','Dm_BDGP5.70.6')),
					annot=t2g,
					threads=t
)

rat_ipf <- call_DTU(annot=t2g, boot_data_A=tpm[[1]], boot_data_B=tpm[[2]],
					p_thresh=0.05, dprop_thresh=dt, abund_thresh=at, scaling=25, seed=67L,
					qboot=TRUE, qbootnum=100, qrep_thresh=qt, rboot=TRUE, rrep_thresh=rt, threads=t,
					varname="condition", name_A="Control", name_B="Modified",
					description="Annotation: BDGP5/Ensembl70. Experiment: Benchmarking DTU with simulated dataset from Soneson 2016 Genom Biol."
)

saveRDS(rat_ipf, file=out)
