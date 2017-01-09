#================================================================================
#' Summary of DTU calling.
#' 
#' @param dtuo A DTU object.
#' @return A named numerical vector giving a tally of the results
#'
#'@export
dtu_summary <- function(dtuo) {
  result <- c("DTU genes" = sum(dtuo$Genes[["DTU"]], na.rm=TRUE), 
              "non-DTU genes" = sum(!dtuo$Genes[["DTU"]], na.rm=TRUE), 
              "NA genes" = sum(ifelse(is.na(dtuo$Genes[["DTU"]]), 1, 0)),
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
#' @return A list of vectors.
#'
#'@export
get_dtu_ids <- function(dtuo) {
  with(dtuo, {
    # Sort transcripts.
    myt <- copy(dtuo$Transcripts)
    myt[, adp := abs(Dprop)]
    setorder(myt, -adp, na.last=TRUE)
    # Sort genes.
    pid <- unique(myt[, parent_id])
    po <- match(Genes[, parent_id], pid)
    myp <- copy(Genes[order(po), ])
    
    # Extract.
    return(list("dtu-genes" = as.vector(myp[(DTU), parent_id]),
               "dtu-transc" = as.vector(myt[(DTU), target_id]),
               "ndtu-genes" = as.vector(myp[DTU==FALSE, parent_id]),
               "ndtu-transc" = as.vector(myt[DTU==FALSE, target_id]),
               "na-genes" = as.vector(myp[is.na(DTU), parent_id]),
               "na-transc" = as.vector(myt[is.na(DTU), target_id])
    ))
  })
}


#================================================================================
#' Plot count and proportion changes for all transcripts of a specified gene.
#'
#' @param dtuo A DTU object.
#' @param pid A \code{parent_id} to make the plot for.
#' @param style Different themes: "lines", "points", "rainbow", "merged", "dashed".
#' @return a ggplot2 object. Simply display it or you can also customize it.
#'
#' @import data.table
#' @import ggplot2
#' @export
plot_gene <- function(dtuo, pid, style="lines") {
  # Slice the data to get just the relevant transcripts.
  with(dtuo, {
    trdat <- dtuo$Transcripts[pid, .(target_id, as.character(DTU))]  # transformed DTU will appear as V2
    trdat[is.na(V2), V2 := "NA"]
    repdatA <- dtuo$ReplicateData[["condA"]][pid, -"parent_id", with=FALSE]
    repdatB <- dtuo$ReplicateData[["condB"]][pid, -"parent_id", with=FALSE]
    
    # Restructure
    trnum <- dim(trdat)[1]
    anum <- dtuo$Parameters$num_replic_A
    bnum <- dtuo$Parameters$num_replic_B
    sA <- colSums(repdatA[, -"target_id", with= FALSE])
    sB <- colSums(repdatB[, -"target_id", with= FALSE])
    pA <- sweep(repdatA[, -"target_id", with= FALSE], 2, sA, "/")
    pB <- sweep(repdatB[, -"target_id", with= FALSE], 2, sB, "/")
    vis_data <- data.table("vals"=             c( unlist(repdatA[, -"target_id", with=FALSE]), unlist(repdatB[, -"target_id", with=FALSE]),                               unlist(pA),                               unlist(pB)  ),
                           "condition"= with(dtuo, c( rep.int(Parameters$cond_A, trnum * anum),    rep.int(Parameters$cond_B, trnum * bnum), rep.int(Parameters$cond_A, trnum * anum), rep.int(Parameters$cond_B, trnum * bnum) )),
                           "isoform"=   unlist(list( with(repdatA, rep.int(target_id, anum)),     with(repdatB, rep.int(target_id, bnum)),  with(repdatA, rep.int(target_id, anum)),  with(repdatB, rep.int(target_id, bnum)) )),
                           "DTU"=                 unlist(list( with(trdat, rep.int(V2, anum)),             with(trdat, rep.int(V2, bnum)),          with(trdat, rep.int(V2, anum)),          with(trdat, rep.int(V2, bnum)) )),
                           "type"=                          c( rep.int("Counts", trnum * anum),             rep.int("Counts", trnum * bnum),     rep.int("Proportions", trnum * anum),      rep.int("Proportions", trnum * bnum) ),
                           "replicate"=               as.factor(c( rep(seq(1, anum), each= trnum),     rep(seq(1, bnum), each=trnum),           rep(seq(1, anum), each= trnum),   rep(seq(1, bnum), each=trnum) ) )
                           # "sample"=               as.factor(c( rep(seq(1, anum), each= trnum),     rep(seq(1+anum, anum+bnum), each=trnum),           rep(seq(1, anum), each= trnum),   rep(seq(anum+1, anum+bnum), each=trnum) ) )
                           )
    
    # Plot.
    dtucol <- c("TRUE"="red", "FALSE"="steelblue3", "NA"="yellow")
    dtupnt <- c("TRUE"="darkred", "FALSE"="darkblue", "NA"="yellow")
    cndcol <- c("darkgreen", "darkorange")
    shapes <- seq(0, 25)
    dtulin <- c("TRUE"="solid", "FALSE"="dashed", "NA"="dotted")
    result <- NULL
    
    if (style=="lines") {
      result <- ggplot(vis_data, aes(x= isoform, y= vals)) +
        facet_grid(type ~ condition, scales= "free") +
        geom_path(aes(colour= replicate, group= replicate)) +
        geom_boxplot(aes(fill= DTU), alpha=0.2, outlier.shape= NA) +
        scale_fill_manual(values= dtucol)
    } else if (style=="points") {
      result <- ggplot(vis_data, aes(x= isoform, y= vals)) +
        facet_grid(type ~ condition, scales= "free") +
        geom_point(aes(colour= DTU, shape= replicate), position= position_jitterdodge(), size= rel(1), stroke= rel(1)) +
        geom_boxplot(aes(colour= DTU, fill= DTU), alpha=0.3, outlier.shape= NA) +
        scale_shape_manual(values= shapes) +
        scale_colour_manual(values= dtupnt) +
        scale_fill_manual(values= dtucol)
    } else if (style=="rainbow") {
      result <- ggplot(vis_data, aes(x= isoform, y= vals)) +
        facet_grid(type ~ condition, scales= "free") +
        geom_point(aes(colour= DTU, shape= replicate), position= position_jitterdodge(), size= rel(1), stroke= rel(1)) +
        geom_boxplot(aes(colour= DTU, fill= isoform), alpha=0.3, outlier.shape= NA) +
        scale_colour_manual(values= dtucol) +
        scale_shape_manual(values= shapes)
    } else if (style=="merged") {
      result <- ggplot(vis_data, aes(x= isoform, y= vals)) +
        facet_grid(type ~ ., scales= "free") +
        geom_point(aes(colour= condition, shape= replicate), position= position_jitterdodge(), size= rel(1), stroke= rel(1)) +
        geom_boxplot(aes(fill= condition), alpha=0.3, outlier.shape= NA) +
        scale_colour_manual(values= cndcol) +
        scale_fill_manual(values= cndcol) +
        scale_shape_manual(values= shapes)
    } else if (style=="dashed") {
      result <- ggplot(vis_data, aes(x= isoform, y= vals)) +
        facet_grid(type ~ ., scales= "free") +
        geom_point(aes(colour= condition, shape= replicate), position= position_jitterdodge(), size= rel(1), stroke= rel(1)) +
        geom_boxplot(aes(fill= condition, linetype= DTU), alpha=0.3, outlier.shape= NA) +
        scale_shape_manual(values= shapes) +
        scale_colour_manual(values= cndcol) +
        scale_fill_manual(values= cndcol) +
        scale_linetype_manual(values= dtulin)
    } else {
      stop("Unknown plot style.")
    }
    result <- result +
      scale_y_continuous(limits= c(0, NA)) +
      guides(colour="legend", fill="legend", shape="legend", linetype="legend") +
      labs(title= paste("gene:", pid), y= NULL, x= NULL) +
      theme(title= element_text(size= rel(1.5)),
            axis.text.x= element_text(angle= 90, size= rel(1.5)),
            axis.text.y= element_text(size= rel(1.5)),
            strip.text.x= element_text(size= rel(1.8)),
            strip.text.y= element_text(size= rel(1.8)),
            strip.background= element_rect(fill= "grey95"),
            panel.grid.major= element_line(colour = "grey95"),
            panel.grid.minor= element_blank(),
            panel.background= element_rect(fill = "grey98"),
            legend.text= element_text(size= rel(1.2)),
            legend.title= element_text(size= rel(1.1)) )
    
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
#'   \item{"transc_conf"}{Transcript-level confidence threshold VS. number of DTU positive calls.}
#'   \item{"gene_conf"}{Gene-level confidence threshold VS. number of DTU positive calls.}
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
    } else if (type == "transc_conf") {
      mydata <- data.frame("thresh"=seq(0, 1, 0.01), 
                           "count"= sapply(seq(0, 1, 0.01), function(x) {
                                           sum(Transcripts[(boot_dtu_freq >= x), elig_fx & sig], na.rm=TRUE) }))
      result <- ggplot(data = mydata, aes(thresh, count)) +
        geom_freqpoly(stat= "identity", size= 1.5) +
        ggtitle("Confidence VS DTU transcripts") +
        labs(x="DTU call frequency threshold",
             y="Number of transcripts") +
        scale_x_continuous(breaks = seq(0, 1, 0.1))
    } else if (type == "gene_conf") {
      mydata <- data.frame("thresh"=seq(0, 1, 0.01), 
                           "count"= sapply(seq(0, 1, 0.01), function(x) {
                                           sum(Genes[(boot_dtu_freq >= x), elig_fx & sig], na.rm=TRUE) }))
      result <- ggplot(data = mydata, aes(thresh, count)) +
        geom_freqpoly(stat= "identity", size= 1.5) +
        ggtitle("Confidence VS DTU genes") +
        labs(x="DTU call frequency threshold",
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
      if ("boot_dtu_freq" %in% names(dtuo$Transcripts)) {
        mydata <- dtuo$Transcripts[, list(as.character(target_id), parent_id, Dprop, -log10(pval_corr), boot_dtu_freq, DTU)]
      } else {
        mydata <- dtuo$Transcripts[, list(as.character(target_id), parent_id, Dprop, -log10(pval_corr), DTU)]
      }
      names(mydata)[1] <- "target_id"
      names(mydata)[4] <- "neglogP"
    
      # Plot
      output$plot1 <- renderPlot({
        plot_overview(dtuo, "gene_volcano")
      })
    
      # Assign mouse hover action to corresponding output space.
      output$hover_info <- renderPrint({
        cat("Hover info: \n")
        myhover <- input$plot_hover
        points <- nearPoints(mydata, myhover, threshold= 5)
        if(dim(points)[1] != 0)
          noquote(points[, target_id])
      })
    
      # Assign mouse click action to corresponding output space.
      output$click_info <- renderPrint({
        cat("Click info: \n")
        myclick <- input$plot_click
        points <- nearPoints(mydata, myclick, threshold= 5)
        suppressWarnings({ points[, pval_corr := 10 ^ (0-neglogP)] })
        if(dim(points)[1] != 0)
          points[, .(target_id, parent_id, DTU, Dprop, pval_corr, boot_dtu_freq)]
      })
      
      # Assign mouse double click action to corresponding output space.
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









