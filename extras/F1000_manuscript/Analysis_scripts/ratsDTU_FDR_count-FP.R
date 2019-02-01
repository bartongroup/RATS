# Rscript ratsDTU_FDR_count-FP.R <rats-in> <outfile> <# threads>

args = commandArgs(trailingOnly=TRUE)
reps <- c(1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17)  # 11 was "bad".

threads <- as.integer(args[[3]])
if (threads < 1 || threads > parallel::detectCores(logical= TRUE) )
	stop("Invalid number of threads!")

rat <- as.character(args[1])
if (is.null(rat))
	stop("An RATs result is needed.")

out <- as.character(args[2])
if (is.null(out))
	stop("An output file is needed.")


library(data.table)
library(parallel)

# Count number of DTU transcripts and genes

if (packageVersion("data.table") >= "1.9.8")  # Ensure data.table complies.
	setDTthreads(threads)

abund_range <- c(0, 1, 5, 10, 50, 100)
p_range <- c(0, 0.001, 0.005, 0.01, 0.05, 0.1)
dprop_range <- c(0, 0.01, seq(0.02, 0.3, 0.02), seq(0.4, 1, 0.2))
qrep_range <- c(seq(0, 0.8, 0.2), seq(0.9, 1, 0.01))

mydtu <- readRDS(rat)

result <- list()
result[["dtu"]] <- array(data=
						 	sapply(abund_range, function (at) {
						 		sapply(p_range, function (pt) {
						 			unlist(mclapply(dprop_range, function (dt) {
						 				sapply(qrep_range, function (qt) {
						 					# Apply thresholds and count results.
						 					with(mydtu, {
						 						# Count DTU.
						 						ctA <- at * Parameters[["num_replic_A"]]
						 						ctB <- at * Parameters[["num_replic_B"]]
						 						subt <- data.table(
						 							parent_id = Transcripts$parent_id,
						 							elig = (Transcripts$sumA >= ctA | Transcripts$sumB >= ctB) & Transcripts$totalA != 0 & Transcripts$totalB != 0 & (Transcripts$sumA != Transcripts$totalA | Transcripts$sumB != Transcripts$totalB),
						 							elig_fx = abs(Transcripts$Dprop) >= dt,
						 							sig = Transcripts$pval_corr < pt,
						 							quant_reprod = Transcripts$quant_dtu_freq >= qt
						 						)
						 						tdtu <- subt$elig & subt$sig & subt$elig_fx & subt$quant_reprod

						 						subg <- data.table(
						 							parent_id = Genes$parent_id,
						 							elig_transc = subt[, as.integer(sum(elig, na.rm=TRUE)), by=parent_id][, V1],
						 							elig = NA,
						 							elig_fx = subt[, any(elig_fx), by = parent_id][, V1],
						 							sig = Genes$pval_corr < pt,
						 							quant_reprod = Genes$quant_dtu_freq >= qt
						 						)
						 						subg[, elig := elig_transc >= 2]
						 						gdtu <- subg$elig & subg$sig & subg$elig_fx & subg$quant_reprod

						 						return(c( sum(gdtu, na.rm=TRUE), sum(tdtu, na.rm=TRUE) ))
						 					})
						 				})
						 			}, mc.cores=threads, mc.allow.recursive=FASLE, mc.preschedule=TRUE))
						 		})
						 	}),
						 dim= c(2, length(qrep_range), length(dprop_range), length(p_range), length(abund_range)),
						 dimnames= list(c("gene", "transc"), as.character(qrep_range), as.character(dprop_range), as.character(p_range), as.character(abund_range))
)
result[["elig"]] <- data.frame( "abund"= abund_range,
								"elitra"= sapply(abund_range, function (at) {
									with(mydtu, {
										ctA <- at * Parameters[["num_replic_A"]]
										ctB <- at * Parameters[["num_replic_B"]]
										elig <- (Transcripts$sumA >= ctA | Transcripts$sumB >= ctB)  &  Transcripts$totalA != 0  &  Transcripts$totalB != 0  &  (Transcripts$sumA != Transcripts$totalA | Transcripts$sumB != Transcripts$totalB)
										return(sum(elig, na.rm=TRUE))
									})
								}),
								"eligen"= sapply(abund_range, function (at) {
									with(mydtu, {
										ctA <- at * Parameters[["num_replic_A"]]
										ctB <- at * Parameters[["num_replic_B"]]
										subt <- data.table(
											parent_id = Transcripts$parent_id,
											elig = (Transcripts$sumA >= ctA | Transcripts$sumB >= ctB) & Transcripts$totalA != 0 & Transcripts$totalB != 0 & (Transcripts$sumA != Transcripts$totalA | Transcripts$sumB != Transcripts$totalB)
										)
										subg <- data.table(
											parent_id = Genes$parent_id,
											elig = subt[, as.integer(sum(elig, na.rm=TRUE)), by=parent_id][, V1] >= 2
										)
										return(sum(subg$elig, na.rm=TRUE))
									})
								}) )

saveRDS(result, file=out)
