library(rats)

t2g <- read.csv("./genome/Homo_sapiens.GRCh37.60.t2g.tsv", sep="\t", header=TRUE)

tpm <- fish4rodents(A_paths=file.path('./salmon_quant', c('SRR393010','SRR393011','SRR393012')),
					B_paths=file.path('./salmon_quant', c('SRR393013','SRR393014','SRR393015')),
					annot=t2g,
					threads=16L
)

rat_ipf <- call_DTU(annot=t2g, boot_data_A=tpm[[1]], boot_data_B=tpm[[2]],
					p_thresh=0.05, dprop_thresh=0.2, abund_thresh=0, scaling=25,
					qboot=TRUE, qbootnum=1000, qrep_thresh=0.95, rboot=TRUE, rrep_thresh=0, threads=16L,
					varname="condition", name_A="Control", name_B="IPF",
					description="Annotation: GRCh37/Ensembl60. Experiment: Validation of RATS. Data from doi:10.1371/journal.pone.0068352. Using same annotation & assembly versions as paper."
)

saveRDS(rat_ipf, file="./DE/rats_human_60.RDS")
