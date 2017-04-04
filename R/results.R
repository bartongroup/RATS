#================================================================================
#' Summary of DTU calling.
#' 
#' @param dtuo A DTU object.
#' @return A named numerical vector giving a tally of the results
#'
#'@export
dtu_summary <- function(dtuo) {
  result <- c("DTU genes (gene test)" = sum(dtuo$Genes[["DTU"]], na.rm=TRUE),
              "non-DTU genes (gene test)" = sum(!dtuo$Genes[["DTU"]], na.rm=TRUE),
              "NA genes (gene test)" = sum(ifelse(is.na(dtuo$Genes[["DTU"]]), 1, 0)),
              "DTU genes (transc. test)" = sum(dtuo$Genes[["transc_DTU"]], na.rm=TRUE), 
              "non-DTU genes (transc. test)" = sum(!dtuo$Genes[["transc_DTU"]], na.rm=TRUE), 
              "NA genes (transc. test)" = sum(ifelse(is.na(dtuo$Genes[["transc_DTU"]]), 1, 0)),
              "DTU genes (both tests)" = sum(dtuo$Genes[["transc_DTU"]] & dtuo$Genes[["DTU"]], na.rm=TRUE), 
              "non-DTU genes (both tests)" = sum(!dtuo$Genes[["transc_DTU"]] & !dtuo$Genes[["DTU"]], na.rm=TRUE), 
              "NA genes (both tests)" = sum(ifelse(is.na(dtuo$Genes[["transc_DTU"]]) & is.na(dtuo$Genes[["DTU"]]), 1, 0)),
              "DTU transcripts" = sum(dtuo$Transcripts[["DTU"]], na.rm=TRUE), 
              "non-DTU transcripts" = sum(!dtuo$Transcripts[["DTU"]], na.rm=TRUE), 
              "NA transcripts" = sum(ifelse(is.na(dtuo$Transcripts[["DTU"]]), 1, 0)) )
  return(result)
}


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
  # Sort transcripts.
  myt <- copy(dtuo$Transcripts[, c("DTU", "parent_id", "target_id", "Dprop")])
  with(myt, {
    myt[, adp := abs(Dprop)]
    setorder(myt, -adp, na.last=TRUE)
    # Sort genes.
  })
  pid <- unique(myt[, "parent_id", with=FALSE])
  po <- match(dtuo$Genes[, "parent_id", with=FALSE], pid)
  myp <- copy(dtuo$Genes[order(po), ])
  
  with(myp, {
    # Extract.
    return(list("DTU genes (gene test)" = as.vector( myp[(DTU), parent_id] ),
               "non-DTU genes (gene test)" = as.vector( myp[DTU==FALSE, parent_id] ),
               "NA genes (gene test)" = as.vector( myp[is.na(DTU), parent_id] ),
               "DTU genes (transc. test)" = as.vector( myp[(transc_DTU), parent_id] ),
               "non-DTU genes (transc. test)" = as.vector( myp[transc_DTU==FALSE, parent_id] ),
               "NA genes (transc. test)" = as.vector( myp[is.na(transc_DTU), parent_id] ),
               "DTU genes (both tests)" = as.vector( myp[(DTU & transc_DTU), parent_id] ),
               "non-DTU genes (both tests)" = as.vector( myp[DTU==FALSE & transc_DTU==FALSE, parent_id] ),
               "NA genes (both tests)" = as.vector( myp[is.na(DTU) & is.na(transc_DTU), parent_id] ),
               "DTU transcripts" = as.vector(myt[(DTU), target_id]),
               "non-DTU transcripts" = as.vector(myt[DTU==FALSE, target_id]),
               "NA transcripts" = as.vector(myt[is.na(DTU), target_id])
    ))
  })
}


#================================================================================
#' Plot abundances for all isoforms of a specified gene.
#'
#' Boxplot of absolute and relative abundances for the isoforms of a given gene.
#' The style option allows grouping by condition or by isoform, and provides the option to
#' include the individual measuremnts from each replicate.
#' 
#' @param dtuo A DTU object.
#' @param pid A \code{parent_id} to make the plot for.
#' @param style Different themes: \itemize{
#'  \item{"plain" - Grouped by condition.},
#'  \item{"paired" - Grouped by isoform.},
#'  \item{"points" - Grouped by condition. Show individual measurements as points.},
#'  \item{"pairedpnt" - Grouped by isoform. Show individual measurements as points.},
#'  \item{"lines" - (Default) Grouped by condition. Connect individual measurements with colour-coded lines.}
#'  \item{"linesonly" - Grouped by condition. Connect replicate measurements as colour-coded lines. Hide the boxplots.}
#'  } 
#' @param fillby Applies to the boxplots. Not all options will work with all styles. 
#' \itemize{
#'  \item{"isoform" - Colour fill by isoform.},
#'  \item{"condition" - Colour fill by condition.},
#'  \item{"DTU" - Colour fill by transcript-level DTU result.},
#'  \item{"none" - Uniform fill.} }
#' @param colourby Applies to boxplot outline and points. Not all options will work with all styles.
#' \itemize{
#'  \item{"isoform" - Colour lines by isoform.},
#'  \item{"condition" - Colour lines by condition.},
#'  \item{"DTU" - Colour lines by transcript-level DTU result.},
#'  \item{"none" - Uniform colour.} }
#' @param shapeby Applies to points.
#' \itemize{
#'  \item{"repliate" - Point shape by replicate.},
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
plot_gene <- function(dtuo, pid, style="lines", fillby=NA_character_, colourby=NA_character_, shapeby=NA_character_,
                      isofcolvec=c("red",  "blue", "forestgreen", "purple", "hotpink", "gold3"),
                      dtucolvec= c("TRUE"="firebrick1", "FALSE"="dodgerblue", "NA"="gold"),
                      condcolvec=c("grey20", "white"), 
                      replcolvec=c("orange", "darkgreen", "purple"),
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
    trnum <- dim(trdat)[1]
    anum <- dtuo$Parameters$num_replic_A
    bnum <- dtuo$Parameters$num_replic_B
    sA <- colSums(repdatA[, -"target_id", with= FALSE])
    sB <- colSums(repdatB[, -"target_id", with= FALSE])
    pA <- sweep(repdatA[, -"target_id", with= FALSE], 2, sA, "/")
    pB <- sweep(repdatB[, -"target_id", with= FALSE], 2, sB, "/")
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
    shaplt <- list("DTU"= c("TRUE"=16, "FALSE"=0, "NA"=4),
                   "condition"= c(19,21),
                   "isoform"= seq.int(0, 24, 1),
                   "replicate" =seq.int(0, 24, 1),
                   "none"= 20)
    
    # Plot.
    result <- NULL
    ###
    if (style=="plain") {
      if (is.na(fillby))
        fillby <- "condition"
      if (is.na(colourby))
        colourby <- "DTU"
      if (colourby=="replicate")  
        stop("This style cannot be coloured by replicate!")
      shapeby="none"
      result <- ggplot(vis_data, aes(x= isoform, y= vals, fill= vis_data[[fillby]], colour=vis_data[[colourby]])) +
                  facet_grid(type ~ condition, scales= "free") +
                  geom_boxplot(alpha=0.5, outlier.shape= NA) +
                  scale_fill_manual(values= colplt[[fillby]], name=fillby) + 
                  scale_colour_manual(values= colplt[[colourby]], name=colourby)
    ###
    } else if ( style == "paired" ) {
      if(is.na(fillby)) {
        if (is.na(colourby) || colourby != "condition") {
          fillby <- "condition"
        } else {
          fillby <- "DTU"
        }
      } else if (fillby != "condition"){
        if (!is.na(colourby) && colourby != "condition") {
            stop("This style requires either fillby or colourby to be set to \"condition\".")
        } else {
          colourby <- "condition"
        }
      }
      if(is.na(colourby))
        colourby <- "DTU"
      if (colourby=="replicate")
        stop("This style cannot be coloured by replicate!")
      shapeby="none"
      result <- ggplot(vis_data, aes(x= isoform, y= vals, fill= vis_data[[fillby]], colour=vis_data[[colourby]])) +
                  facet_grid(type ~ ., scales= "free") +
                  geom_boxplot(alpha=0.5, outlier.shape= NA, width=0.5) +
                  scale_fill_manual(values= colplt[[fillby]], name=fillby) + 
                  scale_colour_manual(values= colplt[[colourby]], name=colourby)
    ###
    } else if (style=="points") {
      if (is.na(fillby))
        fillby <- "condition"
      if (is.na(colourby))
        colourby <- "DTU"
      if (is.na(shapeby)) {
        if (all("DTU" != c(fillby, colourby))) {
          shapeby <- "DTU"
        } else {
          shapeby <- "none"
        }
      }
      result <- ggplot(vis_data, aes(x= isoform, y= vals)) +
                  facet_grid(type ~ condition, scales= "free") +
                  geom_point(aes(colour= vis_data[[colourby]], shape=vis_data[[shapeby]]), position= position_jitterdodge(), stroke= rel(0.8)) +
                  geom_boxplot(aes(fill= vis_data[[fillby]]), alpha=0.2, outlier.shape= NA) +
                  scale_shape_manual(values=shaplt[[shapeby]], name=shapeby) +
                  scale_fill_manual(values= colplt[[fillby]], name=fillby) + 
                  scale_colour_manual(values= colplt[[colourby]], name=colourby)
    ###
    } else if ( style == "pairedpnt" ) {
      if(is.na(fillby)) {
        if (is.na(colourby) || colourby != "condition") {
          fillby <- "condition"
        } else {
          fillby <- "DTU"
        }
      } else if (fillby != "condition"){
        if (!is.na(colourby) && colourby != "condition") {
          stop("This style requires either fillby or colourby to be set to \"condition\".")
        } else {
          colourby <- "condition"
        }
      }
      if(is.na(colourby))
        colourby <- "DTU"
      if (colourby=="replicate")
        stop("This style cannot be coloured by replicate!")
      if (is.na(shapeby)) {
        if (all("DTU" != c(fillby, colourby))) {
          shapeby <- "DTU"
        } else {
          shapeby <- "none"
        }
      }
      result <- ggplot(vis_data, aes(x= isoform, y= vals)) +
                  facet_grid(type ~ ., scales= "free") +
                  geom_jitter(aes(colour= vis_data[[colourby]], shape=vis_data[[shapeby]]), position=position_dodge(width=0.5), stroke= rel(0.8)) +
                  geom_boxplot(aes(fill= vis_data[[fillby]]), alpha=0.2, outlier.shape= NA) +
                  scale_shape_manual(values= shaplt[[shapeby]], name=shapeby) +
                  scale_fill_manual(values= colplt[[fillby]], name=fillby) + 
                  scale_colour_manual(values= colplt[[colourby]], name=colourby)
    ###
    } else if (style=="lines") {
      if (is.na(fillby))
        fillby <- "DTU"
      if (is.na(colourby))
        colourby <- "replicate"
      shapeby="none"
      result <- ggplot(vis_data, aes(x= isoform, y= vals, fill= vis_data[[fillby]])) +
                  facet_grid(type ~ condition, scales= "free") +
                  geom_path(aes(colour= replicate, group= replicate)) +
                  geom_boxplot(alpha=0.2, outlier.shape= NA) +
                  scale_fill_manual(values= colplt[[fillby]], name=fillby)
    ###
    } else if (style=="linesonly") {
      if (is.na(fillby))
        fillby <- "none"
      if (is.na(colourby))
        colourby <- "replicate"
      shapeby="none"
      result <- ggplot(vis_data, aes(x= isoform, y= vals, colour= replicate)) +
                  facet_grid(type ~ condition, scales= "free") +
                  geom_path(aes(group= replicate))
    ###
    } else {
      stop("Unknown plot style.")
    }
    result <- result +
                scale_y_continuous(limits= c(0, NA)) +
                guides(shape="legend") +
                labs(title= paste("gene:", pid), y= "Abundance (Relative & Absolute)", x= "Isoform") +
                theme(axis.text.x= element_text(angle= 90),
                      strip.background= element_rect(fill= "grey95"),
                      panel.grid.major= element_line(colour = "grey95"),
                      panel.grid.minor= element_blank(),
                      panel.background= element_rect(fill = "grey98") )
    if (fillby == "none")
      result <- result + guides(fill="none")
    if (colourby == "none")
      result <- result + guides(colour="none")
    if (shapeby == "none")
      result <- result + guides(shape="none")
    
    return(result)
  })
}


#================================================================================
#' Plot DTU results from the proportions test.
#' 
#' @param dtuo A DTU object.
#' @param type Type of plot. \itemize{
#'   \item{"volcano"}{Change in proportion VS. statistical significance. Done at the transcript level. (Default)}
#'   \item{"maxdprop"}{Distribution of biggest change in proportion in each gene.}
#'   \item{"transc_quant"}{Transcript-level quantification reproducibility threshold VS. number of DTU positive calls.}
#'   \item{"gene_quant"}{Gene-level quantification reproducibility threshold VS. number of DTU positive calls.}
#'   \item{"transc_rep"}{Transcript-level replication reproducibility threshold VS. number of DTU positive calls.}
#'   \item{"gene_rep"}{Gene-level replication reproducibility threshold VS. number of DTU positive calls.}

#' }
#' @return a ggplot2 object. Simply display it or you can also customize it.
#' 
#' Generally uses the results of the transcript-level proportion tests.
#' 
#' @import data.table
#' @import ggplot2
#' @export
plot_overview <- function(dtuo, type="volcano") {
  with(dtuo, {
    if (any(type == c("gene_volcano", "volcano"))) {
      mydata = Transcripts[, .(target_id, Dprop, -log10(pval_corr), DTU)]
      names(mydata)[3] <- "neglogP"
      result <- ggplot(data = mydata, aes(Dprop, neglogP, colour = DTU)) +
        geom_point(alpha = 0.3) +
        ggtitle("Proportion change VS significance") +
        labs(x = paste("Prop in ", Parameters$cond_B, " (-) Prop in ", Parameters$cond_A, sep=""), 
             y ="-log10 (Pval)") +
        scale_color_manual(values=c("steelblue3", "red")) +
        scale_x_continuous(breaks = seq(-1, 1, 0.2)) +
        theme(panel.background= element_rect(fill= "grey98"),
              panel.grid.major= element_line(colour= "grey95") )
    } else if (type == "maxdprop") {
      tmp <- copy(Transcripts)  # I don't want the intermediate calculations to modify the dtu object.
      tmp[, abma := abs(Dprop)]
      tmp <- with(tmp, data.table(aggregate(abma, by=list(parent_id), FUN = max)))
      # Also want coloured by dtu, so I need to retrieve that into a vector that matches tmp.
      setkey(tmp, Group.1)
      tmp[, dtu := Genes[match(tmp$Group.1, Genes[, parent_id]), Genes$DTU] ]
      # ok, plotting time
      result <- ggplot(data = na.omit(tmp), aes(x, fill=dtu)) +
        geom_histogram(binwidth = 0.01, position="identity", alpha = 0.5) +
        ggtitle("Distribution of largest proportion change per gene") +
        labs(x = paste("abs( Prop in ", Parameters$cond_B, " (-) Prop in ", Parameters$cond_A, " )", sep=""), 
             y ="Number of genes") +
        scale_fill_manual(values=c("steelblue3", "red")) +
        scale_x_continuous(breaks = seq(0, 1, 0.1)) +
        scale_y_continuous(trans="sqrt")
    } else if (type == "transc_quant") {
      mydata <- data.frame("thresh"=seq(0, 1, 0.01), 
                           "count"= sapply(seq(0, 1, 0.01), function(x) {
                                           sum(Transcripts[(quant_dtu_freq >= x), elig_fx & sig], na.rm=TRUE) }))
      result <- ggplot(data = mydata, aes(thresh, count)) +
        geom_freqpoly(stat= "identity", size= 1.5) +
        ggtitle("Quantification reproducibility VS DTU transcripts") +
        labs(x="Reproducibility threshold",
             y="Number of transcripts") +
        scale_x_continuous(breaks = seq(0, 1, 0.1))
    } else if (type == "gene_quant") {
      mydata <- data.frame("thresh"=seq(0, 1, 0.01), 
                           "count"= sapply(seq(0, 1, 0.01), function(x) {
                                           sum(Genes[(quant_dtu_freq >= x), elig_fx & sig], na.rm=TRUE) }))
      result <- ggplot(data = mydata, aes(thresh, count)) +
        geom_freqpoly(stat= "identity", size= 1.5) +
        ggtitle("Quantification reproducibility VS DTU genes") +
        labs(x="Reproducibility threshold",
             y="Number of genes") +
        scale_x_continuous(breaks = seq(0, 1, 0.1))
    } else if (type == "transc_rep") {
      mydata <- data.frame("thresh"=seq(0, 1, 0.01), 
                           "count"= sapply(seq(0, 1, 0.01), function(x) {
                             sum(Transcripts[(rep_dtu_freq >= x), elig_fx & sig], na.rm=TRUE) }))
      result <- ggplot(data = mydata, aes(thresh, count)) +
        geom_freqpoly(stat= "identity", size= 1.5) +
        ggtitle("Replication reproducibility VS DTU transcripts") +
        labs(x="Reproducibility threshold",
             y="Number of transcripts") +
        scale_x_continuous(breaks = seq(0, 1, 0.1))
    } else if (type == "gene_rep") {
      mydata <- data.frame("thresh"=seq(0, 1, 0.01), 
                           "count"= sapply(seq(0, 1, 0.01), function(x) {
                             sum(Genes[(rep_dtu_freq >= x), elig_fx & sig], na.rm=TRUE) }))
      result <- ggplot(data = mydata, aes(thresh, count)) +
        geom_freqpoly(stat= "identity", size= 1.5) +
        ggtitle("Replication reproducibility VS DTU genes") +
        labs(x="Reproducibility threshold",
             y="Number of genes") +
        scale_x_continuous(breaks = seq(0, 1, 0.1))
    } else {
      stop("Unrecognized plot type!")
    }
    
    return(result)
  })
}


#================================================================================
#' Interactive volcano plot, using shiny.
#' 
#' @param dtuo A DTU object.
#' 
#' @import data.table
#' @import ggplot2
#' @import shiny
#' @export
#' 
plot_shiny_volcano <- function(dtuo) {
  # Set up interface.
  volcano_ui <- fluidPage(
    # Prepare space for plot display.
    fluidRow(
      column(width= 12, 
             plotOutput("plot1", height= 800, 
                        hover= hoverOpts(id= "plot_hover"),
                        click= clickOpts(id= "plot_click")) )),
    # Instructions.
    fluidRow(
      column(width=12,
             wellPanel("* Hover over points to see Transcript ID.
                       * Click to get more info on the Transcript and to plot the relevant Gene.") )),
    # Hover and click info.
    fluidRow(
      column(width= 3,
             verbatimTextOutput("hover_info")),
      column(width= 9,
             verbatimTextOutput("click_info")) ),
    fluidRow(
      column(width= 12, 
             plotOutput("plot2", height= 600)) )
  )
    
  with(dtuo, {
    # Set up mouse responses.
    volcano_server <- function(input, output) {
      # For storing which rows have been excluded
      vals <- reactiveValues(
        keeprows = rep(TRUE, nrow(mtcars))
      )
      
      # Set up data
      mydata <- NULL
      if ("quant_dtu_freq" %in% names(dtuo$Transcripts)) {
        if ("rep_dtu_freq" %in% names(dtuo$Transcripts)) {
          mydata <- dtuo$Transcripts[, list(as.character(target_id), parent_id, Dprop, -log10(pval_corr), quant_dtu_freq, rep_dtu_freq, DTU)]
        } else {
          mydata <- dtuo$Transcripts[, list(as.character(target_id), parent_id, Dprop, -log10(pval_corr), quant_dtu_freq, DTU)]
        }
      } else if ("rep_dtu_freq" %in% names(dtuo$Transcripts)) {
        mydata <- dtuo$Transcripts[, list(as.character(target_id), parent_id, Dprop, -log10(pval_corr), rep_dtu_freq, DTU)]
      } else {
        mydata <- dtuo$Transcripts[, list(as.character(target_id), parent_id, Dprop, -log10(pval_corr), DTU)]
      }
      names(mydata)[1] <- "target_id"
      names(mydata)[4] <- "neglogP"
    
      # Plot
      output$plot1 <- renderPlot({
        plot_overview(dtuo, "gene_volcano")
      })
    
      # Assign mouse hover action to hoveri info output space.
      output$hover_info <- renderPrint({
        cat("Hover info: \n")
        myhover <- input$plot_hover
        points <- nearPoints(mydata, myhover, threshold= 5)
        if(dim(points)[1] != 0)
          noquote(points[, target_id])
      })
    
      # Assign mouse click action to click info output space.
      output$click_info <- renderPrint({
        cat("Click info: \n")
        myclick <- input$plot_click
        points <- nearPoints(mydata, myclick, threshold= 5)
        suppressWarnings({ points[, pval_corr := 10 ^ (0-neglogP)] })
        if(dim(points)[1] != 0)
          points[, .(target_id, parent_id, DTU, Dprop, pval_corr, quant_dtu_freq, rep_dtu_freq)]
      })
      
      # Assign mouse click action to gene plot output space.
      output$plot2 <- renderPlot({
        myclick <- input$plot_click
        points <- nearPoints(mydata, myclick, threshold= 5, addDist= TRUE)
        md <- which.min(points[, dist_])
        tid <- points[md, target_id]
        gid <- as.vector(dtuo$Transcripts[(target_id == tid), parent_id])
        if(!is.na(gid[1]))
          plot_gene(dtuo, gid)
      })
    }
    
    # Display
    shinyApp(ui= volcano_ui, server= volcano_server)
  })
}









