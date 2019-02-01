library(rats)

t2g <- read.csv("./genome/Homo_sapiens.GRCh37.71.t2g.tsv", sep="\t", header=TRUE)

tpm <- fish4rodents(A_paths=file.path('./kallisto_quant', c('Hs_GRCh37.71.1','Hs_GRCh37.71.2','Hs_GRCh37.71.3')),
		    B_paths=file.path('./kallisto_quant', c('Hs_GRCh37.71.4','Hs_GRCh37.71.5','Hs_GRCh37.71.6')),
		    annot=t2g, beartext=TRUE, half_cooked=TRUE, threads=6L
)

rat_ipf <- call_DTU(annot=t2g, boot_data_A=tpm[[1]], boot_data_B=tpm[[2]],
					p_thresh=0.05, dprop_thresh=0.05, abund_thresh=0, scaling=40, seed=67L,
					qboot=TRUE, qbootnum=100, qrep_thresh=0.95, rboot=TRUE, rrep_thresh=0.85, threads=6L,
					varname="condition", name_A="Control", name_B="Modified",
					description="Annotation: GRCh37/Ensembl71. Experiment: Benchmarking DTU with simulated dataset from Soneson 2016 Genom Biol. Quantified with Kallisto."
)

saveRDS(rat_ipf, file="./DE/rats_soneson_Hs71-kallisto_th5.RDS")
