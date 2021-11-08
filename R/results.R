########## ########## ########## ########## ########## ########## ########## ########## ##########
# Post-analysis utilities, ie. summaries and plots.
########## ########## ########## ########## ########## ########## ########## ########## ##########


#================================================================================
#' List of DTU ids.
#'
#' Get the IDs for DTU/nonDTU/NA genes and transcripts.
#' The IDs will be ordered by effect size.
#'
#' @param dtuo A DTU object.
#' @return A named list of character vectors.
#'
#'@import data.table
#'@export
get_dtu_ids <- function(dtuo) {
  myt <- copy(dtuo$Transcripts[, c("elig", "DTU", "target_id", "parent_id", "Dprop")])
  myp <- copy(dtuo$Genes[, c("elig", "DTU", "transc_DTU", "parent_id", "maxDprop")])
  with(myp, myp[, transc_elig := myt[, any(elig), by = parent_id][, V1] ])
  
  with(myt, {
    # Sort transcripts.
    myt[, adp := abs(Dprop)]
    setkey(myt, NULL)  # It seems setorder() can throw a hissy fit if a table is already ordered by a key (ie. parent_id).
    setorder(myt, -adp, na.last=TRUE)
  })
  with(myp, {
    # Sort genes to match.
    myp[, adp := abs(maxDprop)]
    setkey(myp, NULL) # Again, minimmising chance of clash with setorder().
    setorder(myp, -adp, na.last=TRUE)
  })
  
  with(dtuo, {
    # Extract.
    return(list("DTU genes (gene test)" = as.vector( myp[(DTU), parent_id] ),
                "non-DTU genes (gene test)" = as.vector( myp[ elig & !DTU, parent_id] ),
                "ineligible genes (gene test)" = as.vector( myp[(!elig), parent_id] ),
                "DTU genes (transc. test)" = as.vector( myp[(transc_DTU), parent_id] ),
                "non-DTU genes (transc. test)" = as.vector( myp[transc_elig & !transc_DTU, parent_id] ),
                "ineligible genes (transc. test)" = as.vector( myp[(!transc_elig), parent_id] ),
                "DTU genes (both tests)" = as.vector( myp[DTU & transc_DTU, parent_id] ),
                "non-DTU genes (both tests)" = as.vector( myp[elig & transc_elig & !DTU & !transc_DTU, parent_id] ),
                "ineligible genes (both tests)" = as.vector( myp[!elig | !transc_elig, parent_id] ),
                "DTU transcripts" = as.vector(myt[(DTU), target_id]),
                "non-DTU transcripts" = as.vector(myt[elig & !DTU, target_id]),
                "ineligible transcripts" = as.vector(myt[(!elig), target_id])
    ))
  })
}


#================================================================================
#' Summary of DTU calling.
#'
#' A tally of DTU genes and transcripts.
#' 
#' @param dtuo A DTU object.
#' @return A data.frame with category names in the first column and values in the second column.
#'
#'@export
dtu_summary <- function(dtuo) {
  tally <- sapply(get_dtu_ids(dtuo), length)
  return( data.frame("category"=names(tally), "tally"=tally, row.names=NULL, stringsAsFactors=FALSE)  )
}


#================================================================================
#' List of genes that switch isoform ranks.
#'
#' Get the IDs of DTU genes where isoform rank switching occurs.
#' Switches of primary and non-primary isoforms are listed separately.
#'
#' @param dtuo A DTU object.
#' @return A named list of character vectors.
#'
#' @import data.table
#' @export
get_switch_ids <- function(dtuo) {
  # Get all the transcripts from DTU genes.
  myt <- copy(dtuo$Transcripts[dtuo$Genes[(dtuo$Genes$DTU), c("parent_id")], c("parent_id", "target_id", "propA", "propB")])
  result <- list()
  if (nrow(myt)==0) {
    result[["Primary switch (gene test)"]] <- character(0)
    result[["Non-primary switch (gene test)"]] <- character(0)
  } else {
    with(myt, {
      myt[, rankA := frank(propA), by=parent_id]
      myt[, rankB := frank(propB), by=parent_id]
      myt[, sw := (rankA!=rankB)]
      myt[, mrA := max(rankA), by=parent_id]
      myt[, mrB := max(rankB), by=parent_id]
      myt[, psw := sw & (rankA==mrA | rankB==mrB)]
    })
    result[["Primary switch (gene test)"]] <- unique(myt$parent_id[myt$psw])
    result[["Non-primary switch (gene test)"]] <- unique(myt$parent_id[myt$sw & !myt$psw])
  }

  # Get all the transcripts from genes with at least one DTU transcript.
  myt <- copy(dtuo$Transcripts[dtuo$Genes[(dtuo$Genes$transc_DTU), c("parent_id")], c("parent_id", "target_id", "propA", "propB")])
  if (nrow(myt)==0) {
    result[["Primary switch (transc. test)"]] <- character(0)
    result[["Non-primary switch (transc. test)"]] <- character(0)
  } else {
    with(myt, {
      myt[, rankA := frank(propA), by=parent_id]
      myt[, rankB := frank(propB), by=parent_id]
      myt[, sw := (rankA!=rankB)]
      myt[, mrA := max(rankA), by=parent_id]
      myt[, mrB := max(rankB), by=parent_id]
      myt[, psw := sw & (rankA==mrA | rankB==mrB)]
    })
    result[["Primary switch (transc. test)"]] <- unique(myt$parent_id[myt$psw])
    result[["Non-primary switch (transc. test)"]] <- unique(myt$parent_id[myt$sw & !myt$psw])
  }

  # Get all the transcripts from DTU genes with at least one DTU isoform.
  myt <- copy(dtuo$Transcripts[dtuo$Genes[(dtuo$Genes$DTU & dtuo$Genes$transc_DTU), c("parent_id")], c("parent_id", "target_id", "propA", "propB")])
  if (nrow(myt)==0) {
    result[["Primary switch (both tests)"]] <- character(0)
    result[["Non-primary switch (both tests)"]] <- character(0)
  } else {
    with(myt, {
      myt[, rankA := frank(propA), by=parent_id]
      myt[, rankB := frank(propB), by=parent_id]
      myt[, sw := (rankA!=rankB)]
      myt[, mrA := max(rankA), by=parent_id]
      myt[, mrB := max(rankB), by=parent_id]
      myt[, psw := sw & (rankA==mrA | rankB==mrB)]
    })
    result[["Primary switch (both tests)"]] <- unique(myt$parent_id[myt$psw])
    result[["Non-primary switch (both tests)"]] <- unique(myt$parent_id[myt$sw & !myt$psw])
  }

  return(result)
}


#================================================================================
#' Summary of isoform switching events.
#'
#' A tally of genes showing isoform switching.
#'
#' @param dtuo A DTU object.
#' @return A data.frame with category names in the first column and values in the second column.
#'
#'@export
dtu_switch_summary <- function(dtuo) {
  tally <- sapply(get_switch_ids(dtuo), length)
  return( data.frame("category"=names(tally), "genes"=tally, row.names=NULL, stringsAsFactors=FALSE)  )
}


#================================================================================
#' List the gene IDs by number of isoforms showing significant change.
#'
#' Get the IDs of DTU genes organised by the number of isoforms affected. This
#' is possible only if transcript-level results are enabled.
#'
#' @param dtuo A DTU object.
#' @return A named numerical vector giving a tally of the results.
#'
#'@export
get_plurality_ids <- function(dtuo){
  with(dtuo, {
    plurality <- dtuo$Transcripts[(DTU), length(target_id), by=parent_id]
    categories <- unique(plurality$V1)
    result <- lapply(categories, function (x) {
      return(plurality[(V1==x), parent_id])
    })
    names(result) <- as.character(categories)

    return(result)
  })
}

#================================================================================
#' Summary of DTU plurality.
#' 
#' A tally of genes based on how many isoforms show significant change.
#'
#' @param dtuo A DTU object.
#' @return A data.frame with category names in the first column and values in the second column.
#'
#'@export
dtu_plurality_summary <- function(dtuo) {
  tally <- sapply(get_plurality_ids(dtuo), length)
  return( data.frame("isof_affected"=names(tally), "num_of_genes"=tally, row.names=NULL, stringsAsFactors=FALSE) )
}



#================================================================================
#' Plot abundances for all isoforms of a specified gene.
#'
#' Boxplot of absolute and relative abundances for the isoforms of a given gene.
#' The style option allows grouping by condition or by isoform.
#'
#' @param dtuo A DTU object.
#' @param pid A \code{parent_id} to make the plot for.
#' @param style Different themes: \itemize{
#'  \item{"byisoform" - Grouped by isoform. Show individual measurements as points.},
#'  \item{"bycondition" - (Default) Grouped by condition. Connect individual measurements with colour-coded lines.}
#'  \item{"lines" - Grouped by condition. Connect replicate measurements as colour-coded lines. Hide the boxplots.}
#'  }
#' @param fillby Applies to the boxplots. Not all options will work with all styles.
#' \itemize{
#'  \item{"isoform" - Colour fill by isoform.},
#'  \item{"condition" - Colour fill by condition.},
#'  \item{"DTU" - Colour fill by transcript-level DTU result.},
#'  \item{"none" - Uniform fill.} }
#' @param colourby Applies to boxplot outline and points. Not all options will work with all styles.
#' \itemize{
#'  \item{"replicate" - Point shape by replicate.},
#'  \item{"isoform" - Colour lines by isoform.},
#'  \item{"condition" - Colour lines by condition.},
#'  \item{"DTU" - Colour lines by transcript-level DTU result.},
#'  \item{"none" - Uniform colour.} }
#' @param shapeby Applies to points.
#' \itemize{
#'  \item{"replicate" - Point shape by replicate.},
#'  \item{"isoform" - Point shape by isoform.},
#'  \item{"condition" - Point shape by condition.},
#'  \item{"DTU" - Point shape by transcript-level DTU result.},
#'  \item{"none" - Uniform shape.} }
#' @param isofcolvec Colour vector for isoform highlighting. Used to build a colorRampPalette.
#' @param dtucolvec Colour vector for DTU highlighting.
#' @param condcolvec Colour vector for condition highlighting.
#' @param replcolvec Colour vector for replicate highlighting. Used to build a colorRampPalette.
#' @param nonecol Colour to use when no colour coding is wanted.
#' @return A ggplot2 object. Simply display/print it or you can also customize it.
#'
#' @import data.table
#' @import ggplot2
#' @export
#' 
plot_gene <- function(dtuo, pid, style="bycondition", fillby=NA_character_, colourby=NA_character_, shapeby=NA_character_,
                      isofcolvec=c("tomato",  "lightblue", "forestgreen", "purple", "hotpink", "gold3"),
                      dtucolvec= c("TRUE"="firebrick1", "FALSE"="dodgerblue", "NA"="gold"),
                      condcolvec=c("grey15", "steelblue"),
                      replcolvec=c("red",  "blue", "green", "violet", "pink", "orange"),
                      nonecol="grey50")
{
  if( any(!is.na(fillby) && all(fillby != c("isoform", "condition", "DTU", "none")),
          !is.na(colourby) && all(colourby != c("isoform", "condition", "DTU", "none", "replicate")),
          !is.na(shapeby) && all(shapeby != c("isoform", "condition", "DTU", "none", "replicate")) ))
    stop("Invalid fillby, colourby or shapeby value!")

  # Slice the data to get just the relevant transcripts.
  with(dtuo, {
    trdat <- dtuo$Transcripts[pid, .(target_id, as.character(DTU))]  # transformed DTU will appear as V2
    trdat[is.na(V2), V2 := "NA"]
    repdatA <- dtuo$Abundances[["condA"]][pid, -"parent_id", with=FALSE]
    repdatB <- dtuo$Abundances[["condB"]][pid, -"parent_id", with=FALSE]

    # Restructure
    trnum <- Genes[pid, known_transc]
    anum <- dtuo$Parameters$num_replic_A
    bnum <- dtuo$Parameters$num_replic_B
    sA <- colSums(repdatA[, -"target_id", with= FALSE])
    sB <- colSums(repdatB[, -"target_id", with= FALSE])
    pA <- sweep(repdatA[, -"target_id", with= FALSE], 2, sA, "/")
    pB <- sweep(repdatB[, -"target_id", with= FALSE], 2, sB, "/")
    pA[is.nan(as.matrix(pA))] <- 0
    pB[is.nan(as.matrix(pB))] <- 0
    vis_data <- data.table("vals"=              c( unlist(repdatA[, -"target_id", with=FALSE]), unlist(repdatB[, -"target_id", with=FALSE]),                               unlist(pA),                               unlist(pB)  ),
                           "condition"= with(dtuo, c( rep.int(Parameters$cond_A, trnum * anum),    rep.int(Parameters$cond_B, trnum * bnum), rep.int(Parameters$cond_A, trnum * anum), rep.int(Parameters$cond_B, trnum * bnum) )),
                           "isoform"=     unlist(list( with(repdatA, rep.int(target_id, anum)),     with(repdatB, rep.int(target_id, bnum)),  with(repdatA, rep.int(target_id, anum)),  with(repdatB, rep.int(target_id, bnum)) )),
                           "DTU"=                  unlist(list( with(trdat, rep.int(V2, anum)),              with(trdat, rep.int(V2, bnum)),           with(trdat, rep.int(V2, anum)),           with(trdat, rep.int(V2, bnum)) )),
                           "type"=                          c( rep.int("Count", trnum * anum),             rep.int("Count", trnum * bnum),     rep.int("Proportion", trnum * anum),     rep.int("Proportion", trnum * bnum)  ),
                           "replicate"=            as.factor(c( rep(seq(1, anum), each= trnum),               rep(seq(1, bnum), each=trnum),           rep(seq(1, anum), each= trnum),            rep(seq(1, bnum), each=trnum) )),
                           "none"= as.factor(c(1))
                           )

    # Colour and shape palettes.
    colplt <- list("DTU"= dtucolvec,
                   "condition"= condcolvec,
                   "isoform"= colorRampPalette(isofcolvec)(length(unique(vis_data$isoform))),
                   "replicate"= colorRampPalette(c(replcolvec))(length(unique(vis_data$replicate))),
                   "none"= nonecol)
    shaplt <- list("DTU"= c("TRUE"=19, "FALSE"=0, "NA"=4),
                   "condition"= c(19,21),
                   "isoform"= seq.int(0, 24, 1),
                   "replicate" =seq.int(0, 24, 1),
                   "none"= 20)

    # Plot.
    result <- NULL

    ### BY ISOFORM.
    if (style=="byisoform") {
      if (is.na(colourby)) {
        colourby <- "condition"
      } else if(is.na(fillby)) {
          fillby="condition"
      }
      if (is.na(fillby))
        fillby <- "isoform"
      if(is.na(shapeby))
        shapeby <- "DTU"
      if (all("condition" != c(colourby, fillby)))
          stop("Either fillby or colourby must be set to 'condition' for this plot to be displayed correctly!")
      result <- ggplot(vis_data, aes(x= isoform, y= vals, colour= vis_data[[colourby]], fill= vis_data[[fillby]])) +
        facet_grid(type ~ ., scales= "free") +
        geom_jitter(aes(shape=vis_data[[shapeby]]), position=position_jitterdodge()) +
        geom_boxplot(position=position_dodge(), size=0.5, alpha=0.3, outlier.shape= NA)
    ### BY CONDITION.
    } else if (style=="bycondition") {
      if (is.na(fillby))
        fillby <- "condition"
      if (is.na(shapeby))
        shapeby <- "DTU"
      colourby <- "replicate"
      result <- ggplot(vis_data, aes(x= isoform, y= vals)) +
        facet_grid(type ~ condition, scales= "free") +
        geom_path(aes(colour= replicate, group= replicate), position=position_dodge(width=0.5), alpha=0.7) +
        geom_point(aes(colour= replicate, group= replicate, shape=vis_data[[shapeby]]), position=position_dodge(width=0.5)) +
        geom_boxplot(aes(fill= vis_data[[fillby]]), size=0.5, alpha=0.25, outlier.shape= NA, colour="grey60")
      if (fillby=="condition")
        result <- result + guides(fill="none")
      if (shapeby=="none")
        result <- result + guides(shape="none")
    ### BY CONDITION LINESONLY.
    } else if (style=="lines") {
      if (is.na(fillby))
        fillby <- "condition"
      if(is.na(shapeby))
        shapeby <- "DTU"
      colourby <- "replicate"
      result <- ggplot(vis_data, aes(x= isoform, y= vals, colour= replicate)) +
        facet_grid(type ~ condition, scales= "free") +
        geom_path(aes(group= replicate, colour= vis_data[[colourby]]), position=position_dodge(width=0.2)) +
        geom_point(aes(group= replicate, colour= vis_data[[colourby]], shape=vis_data[[shapeby]]), position=position_dodge(width=0.2))
    ### ERROR
    } else {
      stop("Unknown plot style.")
    }
    # Styling
    result <- result +
      scale_fill_manual(values=colplt[[fillby]], name=fillby, drop = FALSE) +
      scale_colour_manual(values=colplt[[colourby]], name=colourby, drop = FALSE) +
      scale_shape_manual(values=shaplt[[shapeby]], name=shapeby, drop = FALSE) +
      expand_limits(y=c(0,1)) +
      labs(title= paste("gene:", pid), y= NULL, x= NULL) +
      theme(axis.text.x= element_text(angle= 90),
                      strip.background= element_rect(fill= "grey95"),
                      strip.text.y= element_text(size= rel(1.2)),
                      strip.text.x= element_text(size= rel(1.1)),
                      panel.grid.major.x= element_line(colour = "grey95"),
                      panel.grid.major.y= element_line(colour = "grey95"),
                      panel.grid.minor= element_blank(),
                      panel.background= element_rect(fill = "white"),
                      panel.border = element_rect(colour = "black", fill=NA),
                      legend.key = element_rect(fill = 'white') )
    if ( any(fillby == c("none", "isoform")) )
      result <- result + guides(fill="none")
    if ( any(colourby == c("none", "isoform")) )
      result <- result + guides(colour="none")
    if ( any(shapeby == c("none", "isoform")) )
      result <- result + guides(shape="none")
    return(result)
  })
}


#================================================================================
#' Plot DTU results overview.
#'
#' @param dtuo A DTU object.
#' @param type Type of plot. \itemize{
#'   \item{"tvolcano" - Change in proportion VS. transcript-level statistical significance. (Default)}
#'   \item{"gvolcano" - Largest change in proportion per gene VS. gene-level statistical significance.}
#'   \item{"maxdprop" - Distribution of biggest change in proportion in each gene.}
#'   \item{"dprop" - Distribution of effect size.}
#'   \item{"reprod" - Distribution of gene-level DTU reproducibility.}
#'   \item{"reprodVSdprop" - Transcript-level quantification reproducibility VS effect size.}
#'   \item{"fcvolcano" - Fold-change in abundance VS. statistical significance. Done at the transcript level.}
#'   \item{"fcVSdprop" - Fold-change of abundance VS difference in proportion, ofr each transcript.}}
#' @return A ggplot2 object. Simply display it or you can also customize it.
#'
#' CAUTION: RATs does NOT normalise the input abundances for fold-change calclulations. RATs
#' is NOT intended for study of DTE.
#'
#' These overviews rely on the results of the transcript-level proportion tests. If your DTU
#' object was created without the transcript-level tests, this function will not work.
#'
#' @import data.table
#' @import ggplot2
#' @export
#' 
plot_overview <- function(dtuo, type="volcano") {
  with(dtuo, {
    
    ### VOLCANO
    if (any(type == c("transc_volcano", "tvolcano", "volcano"))) {
      mydata = Transcripts[, .(target_id, Dprop, -log10(pval_corr), DTU)]
      names(mydata)[3] <- "neglogP"
      result <- ggplot() +
        geom_point(aes(x=mydata[DTU==FALSE, Dprop], y=mydata[DTU==FALSE,neglogP], colour=mydata[DTU==FALSE, DTU]), alpha = 0.3, shape=20) +
        geom_point(aes(x=mydata[DTU==TRUE, Dprop], y=mydata[DTU==TRUE,neglogP], colour=mydata[DTU==TRUE, DTU]), alpha = 0.3, shape=20) +
        geom_vline(xintercept= Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        geom_vline(xintercept= -Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        geom_hline(yintercept= -Parameters$p_thresh, colour="grey50") +
        ggtitle("Effect size VS significance (transcript level)") +
        labs(x= "Isoform propotion difference", y = "-log10 (Pval)") +
        scale_x_continuous(breaks = seq(-1, 1, 0.2)) +
        scale_y_continuous(expand=c(0,0))
    
    ### GENE VOLCANO
    } else if (any(type == c("gene_volcano", "gvolcano"))) {
      mydata = Genes[, .(parent_id, maxDprop, -log10(pval_corr), DTU)]
      names(mydata)[3] <- "neglogP"
      result <- ggplot() +
        geom_point(aes(x=mydata[DTU==FALSE, maxDprop], y=mydata[DTU==FALSE,neglogP], colour=mydata[DTU==FALSE, DTU]), alpha = 0.3, shape=20) +
        geom_point(aes(x=mydata[DTU==TRUE, maxDprop], y=mydata[DTU==TRUE,neglogP], colour=mydata[DTU==TRUE, DTU]), alpha = 0.3, shape=20) +
        geom_vline(xintercept= Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        geom_vline(xintercept= -Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        geom_hline(yintercept= -Parameters$p_thresh, colour="grey50") +
        ggtitle("Effect size VS significance (gene level)") +
        labs(x= "Largest difference in isoform propotion", y = "-log10 (Pval)") +
        scale_x_continuous(breaks = seq(-1, 1, 0.2)) +
        scale_y_continuous(expand=c(0,0))
    
    ### GENE VOLCANO 2
    } else if (any(type == c("gtvolcano"))) {
      mydata = Genes[, .(parent_id, maxDprop, -log10(pval_corr), transc_DTU)]
      names(mydata)[3] <- "neglogP"
      result <- ggplot() +
        geom_point(aes(x=mydata[transc_DTU==FALSE, maxDprop], y=mydata[transc_DTU==FALSE,neglogP], colour=mydata[transc_DTU==FALSE, transc_DTU]), alpha = 0.3, shape=20) +
        geom_point(aes(x=mydata[transc_DTU==TRUE, maxDprop], y=mydata[transc_DTU==TRUE,neglogP], colour=mydata[transc_DTU==TRUE, transc_DTU]), alpha = 0.3, shape=20) +
        geom_vline(xintercept= Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        geom_vline(xintercept= -Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        geom_hline(yintercept= -Parameters$p_thresh, colour="grey50") +
        ggtitle("Effect size VS significance (transcript level)") +
        labs(x= "Largest difference in isoform propotion", y = "-log10 (Pval)") +
        scale_x_continuous(breaks = seq(-1, 1, 0.2)) +
        scale_y_continuous(expand=c(0,0))
    
    ### TRADITIONAL VOLCANO
    } else if (any(type == "fcvolcano")) {
      mydata = Transcripts[, .(target_id, log2FC, -log10(pval_corr), DTU)]
      names(mydata)[3] <- "neglogP"
      result <- ggplot(data = na.omit(mydata), aes(log2FC, neglogP, colour = DTU)) +
        geom_point(shape=16, alpha = 0.3) +
        ggtitle("Isoform abundance fold-change VS significance") +
        labs(x = "log2 (FC)", y = "-log10 (Pval)")
    
    ### EFFECT SIZE
    } else if (type == "dprop") {
      result <- ggplot(data= Transcripts[(elig), .(Dprop, DTU)], aes(x=Dprop, fill=DTU)) +
        geom_histogram(binwidth = 0.02, position="identity", alpha = 0.5) +
        geom_vline(xintercept= Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        geom_vline(xintercept= -Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        ggtitle("Effect size (transcript level)") +
        labs(x = "Isoform proportion difference", y = "Number of Transcripts") +
        scale_x_continuous(breaks = seq(-1, 1, 0.2)) +
        theme( axis.line.x = element_line() )
      maxy <- max(ggplot_build(result)$data[[1]]$y, na.rm=TRUE)
      maxy <- maxy + maxy * 0.05
      result <- result +
        scale_y_continuous(limits=c(0, maxy), expand=c(0, 0), trans="sqrt",
                           breaks=c(100, 500, 1000, pretty(c(0, maxy), n=4)))
    
    ### MAXDPROP
    } else if (type == "maxdprop") {
      mydata = Genes[, .(parent_id, maxDprop, DTU)]
      result <- ggplot(data = na.omit(mydata), aes(x=maxDprop, fill=DTU)) +
        geom_histogram(binwidth = 0.02, position="identity", alpha = 0.4) +
        geom_vline(xintercept= Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        geom_vline(xintercept= -Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        ggtitle("Effect size (gene level)") +
        labs(x = "Largest isoform proportion difference per gene", y = "Number of genes") +
        scale_x_continuous(limits=c(-1, 1), breaks = seq(-1, 1, 0.2)) +
        theme(axis.line.x= element_line())
      maxy <- max(ggplot_build(result)$data[[1]]$y, na.rm=TRUE)
      maxy <- maxy + maxy * 0.05
      result <- result +
        scale_y_continuous(limits=c(0, maxy), expand=c(0, 0), trans="sqrt",
                                       breaks=pretty(c(0, maxy), n=5))
    
    ### FC vs DPROP
    } else if (type == "fcVSdprop") {
      result <- ggplot(data = na.omit(Transcripts[, .(log2FC, Dprop, DTU)]), aes(x=Dprop, y=log2FC, colour=DTU)) +
        geom_point(shape=20, alpha = 0.3) +
        geom_vline(xintercept= Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        geom_vline(xintercept= -Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        ggtitle("Transcript abundance fold-change VS isoform proportion change") +
        labs(x = "Proportion difference", y = "log2 (FC)") +
        theme( panel.grid.minor.y= element_line(colour= "grey95", size=rel(1.5)) )
    
    ### REPRODUCIBILITY vs DPROP
    } else if (type == "reprodVSdprop") {
      result <- ggplot(data= na.omit(Transcripts[, .(Dprop, quant_dtu_freq, DTU)]), aes(x=Dprop, y=quant_dtu_freq, colour=DTU)) +
        geom_point(shape=20, alpha=0.3) +
        geom_vline(xintercept= Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        geom_vline(xintercept= -Parameters$dprop_thresh, colour="grey50", size=rel(0.5)) +
        geom_hline(yintercept= Parameters$quant_reprod_thresh, colour="grey50", size=rel(0.5)) +
        ggtitle("Reproducibility VS effect size") +
        labs(x = "Isoform proportion change", y = "Pass frequency") +
        scale_x_continuous(limits=c(-1, 1), breaks = seq(-1, 1, 0.2)) +
        scale_y_continuous(limits= c(0, 1.01), expand= c(0, 0)) +
        theme( axis.line.x = element_line() )
    
    ### REPRODUCIBILITY
    } else if (type == "reprod") {
      result <- ggplot(data= na.omit(Genes[, .(quant_dtu_freq, DTU)]), aes(x=quant_dtu_freq, fill=DTU)) +
        geom_histogram(binwidth = 0.02, position="identity", alpha = 0.5) +
        geom_vline(xintercept= Parameters$quant_reprod_thresh, colour="grey50", size=rel(0.5)) +
        scale_x_continuous(breaks = seq(0, 1, 0.2), expand=c(0, 0)) +
        ggtitle("DTU Reproducibility") +
        labs(x = "Pass frequency", y = "Number of Genes") +
        theme( axis.line = element_line() )
      maxy <- max(ggplot_build(result)$data[[1]]$y, na.rm=TRUE)
      maxy <- maxy + maxy * 0.05
      result <- result +
        scale_y_continuous(limits=c(0, maxy), expand=c(0, 0), trans="sqrt",
                           breaks=c(100, 500, 1000, pretty(c(0, maxy), n=4))) +
        coord_flip()
    } else {
      stop("Unrecognized plot type!")
    }
    
    # Apply theme.
    result <- result +
      scale_fill_manual(values=c("lightblue", "red")) +
      scale_colour_manual(values=c("lightblue", "red")) +
      guides(colour=guide_legend("DTU")) +
      theme(panel.background= element_rect(fill= "white"),
                      panel.grid.major= element_line(colour= "grey95"),
                      # panel.border = element_rect(colour = "black", fill=NA),
                      legend.key = element_rect(fill = 'white') )

    return(result)
  })
}

#================================================================================
#' Plot diagnostics.
#' 
#' @param dtuo A DTU object.
#' @param type Type of plot. \itemize{
#'   \item{"cormat" - Pairwise Pearson's correlation matrix among samples.}
#' }
#' @return A ggplot2 object. Simply display it or you can also customize it.
#' 
#' @import data.table
#' @import ggplot2
#' @export
#' 
plot_diagnostics <- function(dtuo, type="cormat") {
  with(dtuo, {
    
    ### CORRELATIONS
    if (type == 'cormat') {
      mydata <- cbind(dtuo[['Abundances']][['condA']][, -c('target_id', 'parent_id')], 
                      dtuo[['Abundances']][['condB']][, -c('target_id', 'parent_id')])
      names(mydata)[seq.int(1,dtuo$Parameters$num_replic_A)] <- paste0(dtuo$Parameters$cond_A, '_', names(mydata)[seq.int(1,dtuo$Parameters$num_replic_A)])
      names(mydata)[seq.int(dtuo$Parameters$num_replic_A+1, dtuo$Parameters$num_replic_A + dtuo$Parameters$num_replic_B)] <- paste0(dtuo$Parameters$cond_B, '_', names(mydata)[seq.int(dtuo$Parameters$num_replic_A+1, dtuo$Parameters$num_replic_A + dtuo$Parameters$num_replic_B)])
      # Do all the pairwise correlations between the columns. Must leave out all rows that contain NA, otherwise all correlations become NA.
      mc <- data.frame(cor(mydata[complete.cases(mydata)]))
      mcn <- names(mc)
      mc['sample1'] <- mcn
      corels <- melt(data.table(mc), id.vars="sample1", variable.name = "sample2")
      result <- ggplot() +
        geom_tile(aes(x=corels[[1]], y=corels[[2]], fill=corels[[3]])) +
        scale_fill_gradient(low="purple", high="white") +
        xlab('Sample') + ylab('Sample') + ggtitle('Pairwise Pearson\'s correlations') + labs(fill='corr') +
        theme(axis.text.x=element_text(angle=90))
    }
    
    # Apply theme.
    result <- result +
      theme(panel.background= element_rect(fill= "white"),
                      legend.key = element_rect(fill = 'white') )
    
    return(result)
  })
}
  

#================================================================================
#' Interactive volcano plot, using shiny.
#'
#' This depends on the shiny package.
#'
#' @param dtuo A DTU object.
#'
#' @import data.table
#' @export
#'
plot_shiny_volcano <- function(dtuo) {
  # Set up interface.
  volcano_ui <- shiny::fluidPage(
    # Prepare space for plot display.
    shiny::fluidRow(
      shiny::column(width= 12,
                    shiny::plotOutput("plot1", height= 600,
                          hover= shiny::hoverOpts(id= "plot_hover"),
                          click= shiny::clickOpts(id= "plot_click")) )
    ),
    # Hover info.
    shiny::fluidRow(
      shiny::column(width= 12,
                    shiny::verbatimTextOutput("hover_info") )
    ),
    # Click info.
    shiny::fluidRow(
      shiny::column(width= 12,
                    shiny::verbatimTextOutput("gene_info") )
    ),
    shiny::fluidRow(
      shiny::column(width= 12,
                    shiny::verbatimTextOutput("transc_info") )
    ),
    shiny::fluidRow(
      shiny::column(width= 12,
                    shiny::plotOutput("plot2", height= 600)) )
  )

  with(dtuo, {
    # Set up mouse responses.
    volcano_server <- function(input, output) {
      # For storing which rows have been excluded
      vals <- shiny::reactiveValues(
        keeprows = rep(TRUE, nrow(mtcars))
      )

      # Set up data
      myp <- dtuo$Genes[, list(as.character(parent_id), maxDprop, -log10(pval_corr))]
      names(myp)[1] <- "parent_id"
      names(myp)[3] <- "neglogP"

      # Plot
      output$plot1 <- shiny::renderPlot({
        plot_overview(dtuo, "gvolcano")
      })

      # Assign mouse hover action to hover info output space.
      output$hover_info <- shiny::renderPrint({
        cat("Mouse-hover info: \n")
        myhover <- input$plot_hover
        points <- shiny::nearPoints(myp, myhover, xvar="maxDprop", yvar="neglogP", threshold= 5)
        if(dim(points)[1] != 0) {
          pid <- noquote(points[, parent_id])
          dtuo$Genes[pid, .(parent_id, known_transc, elig_transc, maxDprop, pval_corr)]
        }
      })

      # Assign mouse click action to click info output space.
      output$gene_info <- shiny::renderPrint({
        cat("Gene info (left click): \n")
        myclick <- input$plot_click
        points <- shiny::nearPoints(myp, myclick, xvar="maxDprop", yvar="neglogP", threshold= 5)
        if(dim(points)[1] != 0) {
          pid <- noquote(points[, parent_id])
          dtuo$Genes[pid, ]
        }
      })
      output$transc_info <- shiny::renderPrint({
        cat("Transcript info (left click): \n")
        myclick <- input$plot_click
        points <- shiny::nearPoints(myp, myclick, xvar="maxDprop", yvar="neglogP", threshold= 5)
        if(dim(points)[1] != 0) {
          pid <- noquote(points[, parent_id])
          dtuo$Transcripts[pid, ]
        }
      })

      # Assign mouse click action to gene plot output space.
      output$plot2 <- shiny::renderPlot({
        myclick <- input$plot_click
        points <- shiny::nearPoints(myp, myclick, xvar="maxDprop", yvar="neglogP", threshold= 5, addDist=TRUE)
        md <- which.min(points[, dist_])
        pid <- points[md, parent_id]
        if(!is.na(pid[1]))
          plot_gene(dtuo, pid, style="byisoform")
      })
    }

    # Display
    shiny::shinyApp(ui= volcano_ui, server= volcano_server)
  })
}









