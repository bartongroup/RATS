library(rats)

t2g <- read.csv("./genome/Drosophila_melanogaster.BDGP5.70.t2g.tsv", sep="\t", header=TRUE)

tpm <- fish4rodents(A_paths=file.path('./kallisto_quant', c('Dm_BDGP5.70.1','Dm_BDGP5.70.2','Dm_BDGP5.70.3')),
		    B_paths=file.path('./kallisto_quant', c('Dm_BDGP5.70.4','Dm_BDGP5.70.5','Dm_BDGP5.70.6')),
		    annot=t2g, beartext=TRUE, half_cooked=TRUE, threads=6L
)

rat_ipf <- call_DTU(annot=t2g, boot_data_A=tpm[[1]], boot_data_B=tpm[[2]],
					p_thresh=0.05, dprop_thresh=0.1, abund_thresh=0, scaling=25, seed=67L,
					qboot=TRUE, qbootnum=100, qrep_thresh=0.95, rboot=TRUE, rrep_thresh=0.85, threads=6L,
					varname="condition", name_A="Control", name_B="Modified",
					description="Annotation: BDGP5/Ensembl70. Experiment: Benchmarking DTU with simulated dataset from Soneson 2016 Genom Biol. Quantified with Kallisto."
)

saveRDS(rat_ipf, file="./DE/rats_soneson_Dm70-kallisto_th4.RDS")
