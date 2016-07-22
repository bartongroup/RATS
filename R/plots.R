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
#' 
#' @param dtuo A DTU object.
#' @return A list of vectors.
#'
#'@export
get_dtu_ids <- function(dtuo) {
  with(dtuo,
       return(list("dtu-genes" = as.vector(Genes[DTU==TRUE, parent_id]),  # Necessary to use equality syntax, otherwise NA interfere.
                   "dtu-transc" = as.vector(Transcripts[DTU==TRUE, target_id]),
                   "ndtu-genes" = as.vector(Genes[DTU==FALSE, parent_id]),
                   "ndtu-transc" = as.vector(Transcripts[DTU==FALSE, target_id]),
                   "na-genes" = as.vector(Genes[is.na(DTU), parent_id]),
                   "na-transc" = as.vector(Transcripts[is.na(DTU), target_id])
       )) )
}


#================================================================================
#' Plot count or proportion changes for all transcripts of a specified gene.
#'
#' @param dtuo A DTU object.
#' @param pid A \code{parent_id} to make the plot for.
#' @param vals Values to plot. Either "counts" or "proportions". Default "counts".
#' @param style Style of plot. Either "bars" or "lines". Default "bars".
#' @return a ggplot2 object. Simply display it or you can also customize it.
#'
#' The error bars represent either the error of proportion or 2 standard deviations
#' of the mean count, according to the type of plot requested.
#'
#' @import data.table
#' @import ggplot2
#' @export
plot_gene <- function(dtuo, pid, vals= "proportions", style= "bars") {
  # Slice the data to get just the relevant transcripts.
  vis_data <- with(dtuo, 
                   Transcripts[pid, .(target_id, meanA, meanB, stdevA, stdevB, propA, propB)] )
  with(vis_data, {
    vis_data[, peA := sqrt(propA * (1 - propA) / dtuo$Parameters[["num_replic_A"]]) ]
    vis_data[, peB := sqrt(propB * (1 - propB) / dtuo$Parameters[["num_replic_B"]]) ]
  })
  
  # Choose values to display.
  vA = NA_real_; vB = NA_real_; eA = NA_real_; eB = NA_real_
  if (vals == "counts") {
    vA <- "meanA";  vB <- "meanB";  eA <- "stdevA";  eB <- "stdevB"
  } else if (vals == "proportions"){
    vA <- "propA";  vB <- "propB";  eA <- "peA";   eB <- "peB"
  } else {
    stop("Invalid plot type.")
  }
  
  vis_data <- with(vis_data, {
    # Aggregate, to simplify ggplot commands.
    vis_data[, condA := dtuo$Parameters[["cond_A"]] ]  # Recycle single value.
    vis_data[, condB := dtuo$Parameters[["cond_B"]] ]
    
    data.table("expression" = c(vis_data[[vA]], vis_data[[vB]]),
               "error" = c(vis_data[[eA]], vis_data[[eB]]),
               "errmin" = NA_real_,
               "errmax" = NA_real_,
               "condition" = c(vis_data[, condA], vis_data[, condB]),
               "transcript" = vis_data[, target_id])  # Recycle vector once.
  })
  
  with(vis_data, {
    # 2 standard deviations.
    if (vals == "counts") 
      vis_data[, error := 2 * error]
    # Error limits, mind sensible limits.
    vis_data[, errmin := expression - error]
    vis_data[, errmin := ifelse(errmin<0, 0, errmin)]
    vis_data[, errmax := expression + error]
    if (vals == "proportions")
      vis_data[, errmax := ifelse(errmax>1, 1, errmax)]
    
    if (style == "lines") {
      # Display as overlapping lines (Nick's way of displaying it, but cleaned up).
      return(
        ggplot(data= vis_data, aes(x= transcript, y= expression, colour= condition)) +
          geom_freqpoly(aes(group= condition, colour= condition), stat= "identity", size= 1.5) +
          geom_errorbar(aes(x= transcript, ymin= errmin, ymax= errmax, colour= condition), width= 0.5, size= 0.5) +
          labs(title= pid, y= vals) )
    } else if (style == "bars"){
      # Display as dodged bar chart.
      return(
        ggplot(data= vis_data, aes(x= transcript, y= expression, fill= condition)) +
          geom_bar(aes(group= condition, colour= condition), stat= "identity", position= position_dodge(0.45), width= 0.5) +
          geom_errorbar(aes(x= transcript, ymin= errmin, ymax= errmax),  position= position_dodge(0.45), width= 0.5, size= 0.5) +
          labs(title= pid, y= vals) )
    } else {
      stop("Invalid plot style.")
    }
  })
}


#================================================================================
#' Plot DTU results from the proportions test.
#' 
#' @param dtuo A DTU object.
#' @param type Type of plot. \itemize{
#'   \item{"propVcount"}{Proportion in condition A VS. fragment count in condition A.}
#'   \item{"dpropVcount"}{Change in proportion VS. fragment count in condition A.}
#'   \item{"dpropVsig"}{Change in proportion VS. statistical significance.}
#'   \item{"maxdprop"}{Distribution of biggest change in proportion in each gene.}
#'   \item{"transc-conf"}{Distribution of bootstrapped confidence of transcript-level DTU.}
#'   \item{"gene-conf"}{Distribution of bootstrapped confidence of gene-level DTU.}
#'   \item{"trconfVdtu"}{Transcript-level confidence threshold VS. number of DTU positive calls.}
#'   \item{"gconfVdtu"}{Gene-level confidence threshold VS. number of DTU positive calls.}
#' }
#' @return a ggplot2 object. Simply display it or you can also customize it.
#' 
#' Generally uses the results of the transcript-level proportion tests.
#' 
#' @import data.table
#' @import ggplot2
#' @export
plot_overview <- function(dtuo, type="dpropVsig") {
  with(dtuo, {
    if (type == "propVcount") {
      result <- ggplot(data = Transcripts, aes(sumA, propA, color = DTU)) +
        geom_point(alpha = 0.3) +
        scale_x_continuous(trans = "log10") +
        scale_colour_manual("DTU", values = c("blue", "red")) +
        ggtitle("Relative abundances of transcripts") +
        labs(y = paste("Proportion in ", Parameters$cond_A, sep=""), 
             x = paste("Cumulative gene read-count in ", Parameters$cond_A, sep=""))
    } else if (type == "dpropVcount") {
      result <- ggplot(data = Transcripts, aes(sumA, Dprop, color = DTU)) +
        geom_point(alpha = 0.3) +
        scale_x_continuous(trans="log10") +
        scale_y_continuous(breaks = seq(-1, 1, 0.2)) +
        scale_colour_manual("DTU", values = c("blue", "red")) +
        ggtitle("Abundance change VS gene expression")
        labs(y = paste("prop( ", Parameters$cond_B, " ) - prop( ", Parameters$cond_A, " )", sep=""), 
             x = paste("Cumulative gene read-count in ", Parameters$cond_A, sep=""))
    } else if (type == "dpropVsig") {
      result <- ggplot(data = Transcripts, aes(Dprop, pval_corr, colour = DTU)) +
        geom_point(alpha = 0.3) +
        ggtitle("Proportion change VS significance") +
        labs(x = paste("Prop in ", Parameters$cond_B, " - Prop in ", Parameters$cond_A, sep=""), 
             y ="P-value") +
        scale_x_continuous(breaks = seq(-1, 1, 0.2))
    } else if (type == "maxdprop") {
      tmp <- copy(Transcripts)  # I don't want the intermediate calculations to modify the dtu object.
      tmp[, abma := abs(Dprop)]
      tmp <- with(tmp, data.table(aggregate(abma, by=list(parent_id), FUN = max)))
      # Also want coloured by dtu, so I need to retrieve that into a vector that matches tmp.
      setkey(tmp, Group.1)
      tmp[, dtu := Genes[match(tmp$Group.1, Genes[, parent_id]), Genes$DTU | Genes$DTU] ]
      # ok, plotting time
      result <- ggplot(data = na.omit(tmp), aes(x, fill=dtu)) +
        geom_histogram(binwidth = 0.01, position="identity", alpha = 0.5) +
        scale_x_continuous(breaks = seq(0, 1, 0.1)) +
        scale_y_continuous(trans="sqrt") +
        ggtitle("Distribution of largest proportion change per gene") +
        labs(x = paste("abs( Prop in ", Parameters$cond_B, " - Prop in ", Parameters$cond_A, " )", sep=""), 
             y ="Number of genes")
    } else if (type == "transc-conf") {
      result <- ggplot(data = Transcripts, aes(x=boot_freq, fill=DTU)) +
        geom_histogram(binwidth = 0.01, position="identity", alpha = 0.5) +
        scale_x_continuous(breaks = seq(0, 1, 0.1)) +
        scale_y_continuous(trans="sqrt") +
        ggtitle("Distribution of bootstrapped transcript-level DTU") +
        labs(x="DTU call frequency in bootstraps", 
             y="Number of transcripts")
    } else if (type == "gene-conf") {
      result <- ggplot(data = Genes, aes(x=boot_freq, fill=DTU)) +
        geom_histogram(binwidth = 0.01, position="identity", alpha = 0.5) +
        scale_x_continuous(breaks = seq(0, 1, 0.1)) +
        scale_y_continuous(trans="sqrt") +
        ggtitle("Distribution of bootstrapped gene-level DTU") +
        labs(x="DTU call frequency in bootstraps", 
             y="Number of genes")
    } else if (type == "trconfVdtu") {
      mydata <- data.frame("thresh"=seq(0, 1, 0.01), "count"= sapply(seq(0, 1, 0.01), function(x) {
        sum(Transcripts[(boot_freq >= x), DTU], na.rm=TRUE) }))
      result <- ggplot(data = mydata, aes(thresh, count)) +
        geom_freqpoly(stat= "identity", size= 1.5) +
        scale_x_continuous(breaks = seq(0, 1, 0.1)) +
        ggtitle("Confidence VS DTU transcripts") +
        labs(x="DTU call frequency threshold",
             y="Number of transcripts")
    } else if (type == "gconfVdtu") {
      mydata <- data.frame("thresh"=seq(0, 1, 0.01), "count"= sapply(seq(0, 1, 0.01), function(x) {
        sum(Genes[(boot_freq >= x), DTU], na.rm=TRUE) }))
      result <- ggplot(data = mydata, aes(thresh, count)) +
        geom_freqpoly(stat= "identity", size= 1.5) +
        scale_x_continuous(breaks = seq(0, 1, 0.1)) +
        ggtitle("Confidence VS DTU genes") +
        labs(x="DTU call frequency threshold",
             y="Number of genes")
    } else {
      stop("Unrecognized plot type!")
    }
    
    return(result)
  })
}
