## Rscript scale_table.R <table.tsv> <factor>
args = commandArgs(trailingOnly=TRUE)

infile <- args[[1]]
f <- as.numeric(args[[2]])

# Scale TPMs to target.
df <- read.table(infile, header=TRUE)
df <- df * f
write.table(df, file=paste0(infile,'_f',f), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
