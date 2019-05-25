args <- commandArgs(trailingOnly = TRUE)
at <- as.numeric(args[1])  # abund_thresh
dt <- as.numeric(args[2])  # dprop_thresh
qt <- as.numeric(args[3])  # qrep_thresh
rt <- as.numeric(args[4])  # rrep_thresh
t <- as.numeric(args[5])  # threads
out <- args[6]  # outfile

library(rats)

t2g <- read.csv("./genome/Homo_sapiens.GRCh37.71.t2g.tsv", sep="\t", header=TRUE)

tpm <- fish4rodents(A_paths=file.path('./kallisto_quant', c('Hs_GRCh37.71.1','Hs_GRCh37.71.2','Hs_GRCh37.71.3')),
		    B_paths=file.path('./kallisto_quant', c('Hs_GRCh37.71.4','Hs_GRCh37.71.5','Hs_GRCh37.71.6')),
		    annot=t2g, beartext=TRUE, half_cooked=TRUE, threads=t
)

rat_ipf <- call_DTU(annot=t2g, boot_data_A=tpm[[1]], boot_data_B=tpm[[2]],
					p_thresh=0.05, dprop_thresh=dt, abund_thresh=at, scaling=40, seed=67L,
					qboot=TRUE, qbootnum=100, qrep_thresh=qt, rboot=TRUE, rrep_thresh=rt, threads=t,
					varname="condition", name_A="Control", name_B="Modified",
					description="Annotation: GRCh37/Ensembl71. Experiment: Benchmarking DTU with simulated dataset from Soneson 2016 Genom Biol. Quantified with Kallisto."
)

saveRDS(rat_ipf, file=out)
