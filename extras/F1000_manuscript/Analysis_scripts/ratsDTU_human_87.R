library(rats)

t2g <- read.csv("./genome/Homo_sapiens.GRCh38.87.t2g.tsv", sep="\t", header=TRUE)

tpm <- fish4rodents(A_paths=file.path('./salmon_quant', c('SRR393010_87','SRR393011_87','SRR393012_87')),
					B_paths=file.path('./salmon_quant', c('SRR393013_87','SRR393014_87','SRR393015_87')),
					annot=t2g,
					threads=16L
)

rat_ipf <- call_DTU(annot=t2g, boot_data_A=tpm[[1]], boot_data_B=tpm[[2]],
					p_thresh=0.05, dprop_thresh=0.2, abund_thresh=0, scaling=25,
					qboot=TRUE, qbootnum=1000, qrep_thresh=0.95, rboot=TRUE, rrep_thresh=0, threads=16L,
					varname="condition", name_A="Control", name_B="IPF",
					description="Annotation: GRCh38/Ensembl87. Experiment: Validation of RATS. Data from doi:10.1371/journal.pone.0068352. Using newer annotation & assembly versions as paper."
)

saveRDS(rat_ipf, file="./DE/rats_human_87.RDS")
