library(data.table)

args <- commandArgs(trailingOnly = TRUE)
at <- args[1]             # abund_thresh
if (is.null(at)) {
    at <- 0
}
at <- as.numeric(at)

# Covariates matrix
covars_Dm <- data.frame(sample_id = c('Dm_BDGP5.70.1', 'Dm_BDGP5.70.2', 'Dm_BDGP5.70.3', 'Dm_BDGP5.70.4', 'Dm_BDGP5.70.5', 'Dm_BDGP5.70.6'),
                     condition = c(1, 1, 1, 2, 2, 2) )
covars_Hs <- data.frame(sample_id = c('Hs_GRCh37.71.1', 'Hs_GRCh37.71.2', 'Hs_GRCh37.71.3', 'Hs_GRCh37.71.4', 'Hs_GRCh37.71.5', 'Hs_GRCh37.71.6'),
                        condition = c(1, 1, 1, 2, 2, 2) )
covars_Deng60 <- data.frame(sample_id = c('SRR393010', 'SRR393011', 'SRR393012', 'SRR393013', 'SRR393014', 'SRR393015'),
                            condition = c(1, 1, 1, 2, 2, 2) )
covars_Deng87 <- data.frame(sample_id = c('SRR393010_87', 'SRR393011_87', 'SRR393012_87', 'SRR393013_87', 'SRR393014_87', 'SRR393015_87'),
                        condition = c(1, 1, 1, 2, 2, 2) )

# Data
datapath <- '/homes/kfroussios/PROJECTS/rats/salmon_quant'

Dm1 <- fread(file = file.path(datapath, 'sonesonDm70-A_merged.tsv_f25'), sep = "\t", header=FALSE )
Dm2 <- fread(file = file.path(datapath, 'sonesonDm70-B_merged.tsv_f25'), sep = "\t", header=FALSE )
names(Dm1) <- c('featureID', levels(covars_Dm$sample_id)[1:3])
names(Dm2) <- c('featureID', levels(covars_Dm$sample_id)[4:6])

Hs1 <- fread(file = file.path(datapath, 'sonesonHs71-A_merged.tsv_f40'), sep = "\t", header=FALSE )
Hs2 <- fread(file = file.path(datapath, 'sonesonHs71-B_merged.tsv_f40'), sep = "\t", header=FALSE )
names(Hs1) <- c('featureID', levels(covars_Hs$sample_id)[1:3])
names(Hs2) <- c('featureID', levels(covars_Hs$sample_id)[4:6])

Deng601 <- fread(file = file.path(datapath, 'ipf_merged.tsv_f25'), sep = "\t", header=FALSE )
Deng602 <- fread(file = file.path(datapath, 'wt_merged.tsv_f25'), sep = "\t", header=FALSE )
names(Deng601) <- c('featureID', levels(covars_Deng60$sample_id)[1:3])
names(Deng602) <- c('featureID', levels(covars_Deng60$sample_id)[4:6])

Deng871 <- fread(file = file.path(datapath, 'ipf_87_merged.tsv_f25'), sep = "\t", header=FALSE )
Deng872 <- fread(file = file.path(datapath, 'wt_87_merged.tsv_f25'), sep = "\t", header=FALSE )
names(Deng871) <- c('featureID', levels(covars_Deng87$sample_id)[1:3])
names(Deng872) <- c('featureID', levels(covars_Deng87$sample_id)[4:6])

# Annotations
genpath <- '/homes/kfroussios/PROJECTS/rats/genome'

Dm_t2g <- fread(file = file.path(genpath,'Drosophila_melanogaster.BDGP5.70.t2g.tsv'), sep = "\t", header=TRUE)
names(Dm_t2g) <- c('groupID', 'featureID')
tab <- table(Dm_t2g$groupID)
Dm_t2g$ntx <- tab[match(Dm_t2g$groupID, names(tab))]

Hs_t2g <- fread(file = file.path(genpath,'Homo_sapiens.GRCh37.71.t2g.tsv'), sep = "\t", header=TRUE)
names(Hs_t2g) <- c('groupID', 'featureID')
tab <- table(Hs_t2g$groupID)
Hs_t2g$ntx <- tab[match(Hs_t2g$groupID, names(tab))]

Deng60_t2g <- fread(file = file.path(genpath,'Homo_sapiens.GRCh37.60.t2g.tsv'), sep = "\t", header=TRUE)
names(Deng60_t2g) <- c('groupID', 'featureID')
tab <- table(Deng60_t2g$groupID)
Deng60_t2g$ntx <- tab[match(Deng60_t2g$groupID, names(tab))]

Deng87_t2g <- fread(file = file.path(genpath,'Homo_sapiens.GRCh38.87.t2g.tsv'), sep = "\t", header=TRUE)
names(Deng87_t2g) <- c('groupID', 'featureID')
tab <- table(Deng87_t2g$groupID)
Deng87_t2g$ntx <- tab[match(Deng87_t2g$groupID, names(tab))]

# Tidy up
library(dplyr)

stopifnot(all(Dm1$target_id %in% Dm_t2g$featureID))
stopifnot(all(Dm2$target_id %in% Dm_t2g$featureID))
stopifnot(all(Hs1$target_id %in% Hs_t2g$featureID))
stopifnot(all(Hs2$target_id %in% Hs_t2g$featureID))
stopifnot(all(Deng601$target_id %in% Deng60_t2g$featureID))
stopifnot(all(Deng602$target_id %in% Deng60_t2g$featureID))
stopifnot(all(Deng871$target_id %in% Deng87_t2g$featureID))
stopifnot(all(Deng872$target_id %in% Deng87_t2g$featureID))

Dm <- Dm_t2g %>%
  left_join(Dm1, by = 'featureID') %>%
  left_join (Dm2, by = 'featureID') %>%
  select(-ntx)
Hs <- Hs_t2g %>%
  left_join(Hs1, by = 'featureID') %>%
  left_join (Hs2, by = 'featureID') %>%
  select(-ntx)
Deng60 <- Deng60_t2g %>%
  left_join(Deng601, by = 'featureID') %>%
  left_join (Deng602, by = 'featureID') %>%
  select(-ntx)
Deng87 <- Deng87_t2g %>%
  left_join(Deng871, by = 'featureID') %>%
  left_join (Deng872, by = 'featureID') %>%
  select(-ntx)

Dm[is.na(Dm)] <- 0.0
Hs[is.na(Hs)] <- 0.0
Deng60[is.na(Deng60)] <- 0.0
Deng87[is.na(Deng87)] <- 0.0

# Prefilter.
# DESeq2 calculates abundance prefilters automatically, but try applying my own for
# some comparability to the other tools that have explicit pre-filters.
Dm[Dm < at] <- 0.0
Hs[Hs < at] <- 0.0
Deng60[Deng60 < at] <- 0.0
Deng87[Deng87 < at] <- 0.0


# Format.
library(DEXSeq)

Dm_dxd <- DEXSeqDataSet(countData = as.matrix(round(Dm[, -c(1:2)])),
                        sampleData = covars_Dm,
                        design = ~sample + exon + condition:exon,
                        featureID = Dm$featureID,
                        groupID = Dm$groupID)
Hs_dxd <- DEXSeqDataSet(countData = as.matrix(round(Hs[, -c(1:2)])),
                        sampleData = covars_Hs,
                        design = ~sample + exon + condition:exon,
                        featureID = Hs$featureID,
                        groupID = Hs$groupID)
Deng60_dxd <- DEXSeqDataSet(countData = as.matrix(round(Deng60[, -c(1:2)])),
                            sampleData = covars_Deng60,
                            design = ~sample + exon + condition:exon,
                            featureID = Deng60$featureID,
                            groupID = Deng60$groupID)
Deng87_dxd <- DEXSeqDataSet(countData = as.matrix(round(Deng87[, -c(1:2)])),
                        sampleData = covars_Deng87,
                        design = ~sample + exon + condition:exon,
                        featureID = Deng87$featureID,
                        groupID = Deng87$groupID)

# Test.
depath <- '/homes/kfroussios/PROJECTS/rats/DE'

Dm_dxd <- estimateSizeFactors(Dm_dxd)
Dm_dxd <- estimateDispersions(Dm_dxd)
Dm_dxd <- testForDEU(Dm_dxd, reducedModel = ~sample + exon)
Dm_dxr.t <- DEXSeqResults(Dm_dxd, independentFiltering = FALSE)
saveRDS(Dm_dxr.t, file=file.path(depath, paste0('dexseq', at, '_soneson_Dm70-salmon.RDS')))

Hs_dxd <- estimateSizeFactors(Hs_dxd)
Hs_dxd <- estimateDispersions(Hs_dxd)
Hs_dxd <- testForDEU(Hs_dxd, reducedModel = ~sample + exon)
Hs_dxr.t <- DEXSeqResults(Hs_dxd, independentFiltering = FALSE)
saveRDS(Hs_dxr.t, file=file.path(depath, paste0('dexseq', at, '_soneson_Hs71-salmon.RDS')))

Deng60_dxd <- estimateSizeFactors(Deng60_dxd)
Deng60_dxd <- estimateDispersions(Deng60_dxd)
Deng60_dxd <- testForDEU(Deng60_dxd, reducedModel = ~sample + exon)
Deng60_dxr.t <- DEXSeqResults(Deng60_dxd, independentFiltering = FALSE)
saveRDS(Deng60_dxr.t, file=file.path(depath, paste0('dexseq', at, '_Deng60-salmon.RDS')))

Deng87_dxd <- estimateSizeFactors(Deng87_dxd)
Deng87_dxd <- estimateDispersions(Deng87_dxd)
Deng87_dxd <- testForDEU(Deng87_dxd, reducedModel = ~sample + exon)
Deng87_dxr.t <- DEXSeqResults(Deng87_dxd, independentFiltering = FALSE)
saveRDS(Deng87_dxr.t, file=file.path(depath, paste0('dexseq', at, '_Deng87-salmon.RDS')))

# Post processing.
Dm_qval <- perGeneQValue(Dm_dxr.t)
Dm_dxr.g <- data.table(gene = names(Dm_qval), Dm_qval)
Dm_dxr.t <- as.data.table(Dm_dxr.t[, c('featureID', 'groupID', 'pvalue')])
saveRDS(Dm_dxr.g, file=file.path(depath, paste0('dexseq', at, '-g_soneson_Dm70-salmon.RDS')))

Hs_qval <- perGeneQValue(Hs_dxr.t)
Hs_dxr.g <- data.table(gene = names(Hs_qval), Hs_qval)
Hs_dxr.t <- as.data.table(Hs_dxr.t[, c('featureID', 'groupID', 'pvalue')])
saveRDS(Hs_dxr.g, file=file.path(depath, paste0('dexseq', at, '-g_soneson_Hs71-salmon.RDS')))

Deng60_qval <- perGeneQValue(Deng60_dxr.t)
Deng60_dxr.g <- data.table(gene = names(Deng60_qval), Deng60_qval)
Deng60_dxr.t <- as.data.table(Deng60_dxr.t[, c('featureID', 'groupID', 'pvalue')])
saveRDS(Deng60_dxr.g, file=file.path(depath, paste0('dexseq', at, '-g_Deng60-salmon.RDS')))

Deng87_qval <- perGeneQValue(Deng87_dxr.t)
Deng87_dxr.g <- data.table(gene = names(Deng87_qval), Deng87_qval)
Deng87_dxr.t <- as.data.table(Deng87_dxr.t[, c('featureID', 'groupID', 'pvalue')])
saveRDS(Deng87_dxr.g, file=file.path(depath, paste0('dexseq', at, '-g_Deng87-salmon.RDS')))

library(stageR)

Dm_pConfirmation <- matrix(Dm_dxr.t$pvalue, ncol=1)
dimnames(Dm_pConfirmation) <- list(Dm_dxr.t$featureID, 'transcript')
Dm_stageRObj <- stageRTx(pScreen = Dm_qval,
                         pConfirmation = Dm_pConfirmation,
                         pScreenAdjusted = TRUE,
                         tx2gene = Dm_t2g[,c('featureID', 'groupID')])
Dm_stageRObj <- stageWiseAdjustment(Dm_stageRObj, method = 'dtu', alpha = 0.05, allowNA = TRUE)
Dm_dex.padj <- getAdjustedPValues(Dm_stageRObj, order = FALSE, onlySignificantGenes = FALSE)
saveRDS(Dm_dex.padj, file = file.path(depath, paste0('dexseq', at, '-stageR_soneson_Dm70-salmon.RDS') ))

Hs_pConfirmation <- matrix(Hs_dxr.t$pvalue, ncol=1)
dimnames(Hs_pConfirmation) <- list(Hs_dxr.t$featureID, 'transcript')
Hs_stageRObj <- stageRTx(pScreen = Hs_qval,
                         pConfirmation = Hs_pConfirmation,
                         pScreenAdjusted = TRUE,
                         tx2gene = Hs_t2g[,c('featureID', 'groupID')])
Hs_stageRObj <- stageWiseAdjustment(Hs_stageRObj, method = 'dtu', alpha = 0.05, allowNA = TRUE)
Hs_dex.padj <- getAdjustedPValues(Hs_stageRObj, order = FALSE, onlySignificantGenes = FALSE)
saveRDS(Hs_dex.padj, file = file.path(depath, paste0('dexseq', at, '-stageR_soneson_Hs71-salmon.RDS') ))

Deng60_pConfirmation <- matrix(Deng60_dxr.t$pvalue, ncol=1)
dimnames(Deng60_pConfirmation) <- list(Deng60_dxr.t$featureID, 'transcript')
Deng60_stageRObj <- stageRTx(pScreen = Deng60_qval,
                             pConfirmation = Deng60_pConfirmation,
                             pScreenAdjusted = TRUE,
                             tx2gene = Deng60_t2g[,c('featureID', 'groupID')])
Deng60_stageRObj <- stageWiseAdjustment(Deng60_stageRObj, method = 'dtu', alpha = 0.05, allowNA = TRUE)
Deng60_dex.padj <- getAdjustedPValues(Deng60_stageRObj, order = FALSE, onlySignificantGenes = FALSE)
saveRDS(Deng60_dex.padj, file = file.path(depath, paste0('dexseq', at, '-stageR_Deng60-salmon.RDS') ))

Deng87_pConfirmation <- matrix(Deng87_dxr.t$pvalue, ncol=1)
dimnames(Deng87_pConfirmation) <- list(Deng87_dxr.t$featureID, 'transcript')
Deng87_stageRObj <- stageRTx(pScreen = Deng87_qval,
                         pConfirmation = Deng87_pConfirmation,
                         pScreenAdjusted = TRUE,
                         tx2gene = Deng87_t2g[,c('featureID', 'groupID')])
Deng87_stageRObj <- stageWiseAdjustment(Deng87_stageRObj, method = 'dtu', alpha = 0.05, allowNA = TRUE)
Deng87_dex.padj <- getAdjustedPValues(Deng87_stageRObj, order = FALSE, onlySignificantGenes = FALSE)
saveRDS(Deng87_dex.padj, file = file.path(depath, paste0('dexseq', at, '-stageR_Deng87-salmon.RDS') ))
