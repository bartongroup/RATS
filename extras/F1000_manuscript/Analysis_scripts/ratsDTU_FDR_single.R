# Rscript ratsDTU_FDR_single.R <# samples> <# threads> <outfile> <seed>

args = commandArgs(trailingOnly=TRUE)
reps <- c(1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17)  # 11 was "bad".

print(args)

size <- as.integer(args[1])
if (size * 2L > length(reps))
	stop("There are not enough samples!")

threads <- as.integer(args[[2]])
if (threads < 1 || threads > parallel::detectCores(logical= TRUE) )
	stop("Invalid number of threads!")

out <- as.character(args[3])
if (is.null(out))
	stop("An output path is needed.")

if (!is.na(args[4]))
	set.seed(as.integer(args[4]))


smpls <- size * 2
chosen <- sample(reps, smpls, replace = FALSE)
wt1 <- chosen[seq.int(1, size, 1)]
wt2 <- chosen[seq.int(size + 1, smpls, 1)]

base_dir <- "."
results_dir <- "salmon_quant"


library(rats)

load(file.path(base_dir, "genome/Araport11_t2g.rda"))  # araport11
# Extract TPMs and scale to the smallest library size among the resplicates (74M).
boots <- fish4rodents(file.path(base_dir, paste0("salmon_quant/AtWT_", wt1)),
					            file.path(base_dir, paste0("salmon_quant/AtWT_", wt2)),
					            araport11,
					            threads=threads, scaleto=74000000)


# Run with virtually no thresholds, allowing me to plot graphs afterwards over the entire range of values.
mydtu <- call_DTU(annot=araport11 , boot_data_A=boots$boot_data_A, boot_data_B=boots$boot_data_B,
				  abund_thresh=0, p_thresh=1, dprop_thresh=0,
				  qboot=TRUE, qrep_thresh=0.5, qbootnum=100, rboot=FALSE,
				  threads=threads,
				  name_A="WT1", name_B="WT2", description=paste("WT vs WT drawn from 16 Arabidopsis thaliana replicates, without replacement. Using scaled TPMs. Annotation: Araport11. Samples:", paste0(chosen, collapse=',')) )

rm(araport11)

saveRDS(mydtu, file=file.path(base_dir, out) )





