# Generate artificial 2-transcript-genes read-count data, systematically covering a broad range of proportions and magnitudes.

# Range of proportions. (21 levels)
# Although the relationship between transcript is symmetric and the ranges 0-0.5 and 0.5-1 are equivalent, 
# they are not redundant: The two sibling transcripts receive different treatment when fold-changes are applied.
proportions_A_T1 <- seq(from=0, to=1, by=0.05)
proportions_A_T2 <- 1-proportions_A_T1

# Range of magnitudes. (11 levels)
magnitudes <- seq(from=-5, to=5, by=1)

# Range of ratio fold-changes. (15 levels)
folds <- seq(from=-7, to=7, by=1)


# Combine proportions and magnitudes. (3465 rows)
# Start by griding the values for easier visualization and naming, instead of compbing the outer products directly.
a1 <- expand.grid(proportions_A_T1, magnitudes, folds)  # condition A, transcript 1, to remain the same
a2 <- expand.grid(proportions_A_T2, magnitudes, folds)  # conditions A, transcript 2, to remain the same
b1 <- expand.grid(proportions_A_T1, magnitudes, folds)  # condition B, transcript 1, to be decreased
b2 <- expand.grid(proportions_A_T2, magnitudes, folds)  # condition B, transcript 2, to be increased

names(a1) <- c("prop", "mag", "fold")
names(a2) <- c("prop", "mag", "fold")
names(b1) <- c("prop", "mag", "fold")
names(b2) <- c("prop", "mag", "fold")


# Calculate "counts".
a1["est_counts"] <- a1$prop * (10 ^ a1$mag)
a2["est_counts"] <- a2$prop * (10 ^ a2$mag)
# Apply foldchange symmetrically to the transcripts. The total ratio change is the nominal foldchange.
b1["est_counts"] <- b1$prop * (10 ^ b1$mag) * (1 / sqrt(2 ^ b1$fold))
b2["est_counts"] <- b2$prop * (10 ^ b2$mag) * sqrt(2 ^ b2$fold)


# Name the "transcripts".
opts <- options()
options(scipen=0, digits=3)

a1["parent_id"] <- paste("p", a1$prop, "_m", a1$mag, "_x", a1$fold, sep="")
a2["parent_id"] <- paste("p", a1$prop, "_m", a2$mag, "_x", a2$fold, sep="")
b1["parent_id"] <- paste("p", b1$prop, "_m", b1$mag, "_x", b1$fold, sep="")
b2["parent_id"] <- paste("p", b1$prop, "_m", b2$mag, "_x", b2$fold, sep="")

a1["target_id"] <- paste(a1$parent_id, "_t1", sep="")
a2["target_id"] <- paste(a2$parent_id, "_t2", sep="")
b1["target_id"] <- paste(b1$parent_id, "_t1", sep="")
b2["target_id"] <- paste(b1$parent_id, "_t2", sep="")

options(scipen=opts$scipen, digits=opts$digits)

