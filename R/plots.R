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
plot_gene <- function(dtuo, pid, vals= "counts", style= "bars") {
  # Slice the data to get just the relevant transcripts.
  vis_data <- dtuo$Transcripts[pid, .(target_id, meanA, meanB, stdevA, stdevB, propA, propB)]
  vis_data[, peA := sqrt(propA * (1 - propA) / dtuo$Parameters[["num_replic_A"]]) ]
  vis_data[, peB := sqrt(propB * (1 - propB) / dtuo$Parameters[["num_replic_B"]]) ]
  
  # Choose values to display.
  vA = NA_real_; vB = NA_real_; eA = NA_real_; eB = NA_real_
  if (vals == "counts") {
    vA <- "meanA";  vB <- "meanB";  eA <- "stdevA";  eB <- "stdevB"
  } else if (vals == "proportions"){
    vA <- "propA";  vB <- "propB";  eA <- "peA";   eB <- "peB"
  } else {
    stop("Invalid plot type.")
  }
  
  # Aggregate, to simplify ggplot commands.
  vis_data[, condA := dtuo$Parameters[["cond_A"]] ]  # Recycle single value.
  vis_data[, condB := dtuo$Parameters[["cond_B"]] ]
  vis_data <- data.table("expression" = c(vis_data[[vA]], vis_data[[vB]]),
                  "error" = c(vis_data[[eA]], vis_data[[eB]]),
                  "errmin" = NA_real_,
                  "errmax" = NA_real_,
                  "condition" = c(vis_data[, condA], vis_data[, condB]),
                  "transcript" = vis_data[, target_id])  # Recycle vector once.
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
    result <- ggplot(data= vis_data, aes(x= transcript, y= expression, colour= condition)) +
      geom_freqpoly(aes(group= condition, colour= condition), stat= "identity", size= 1.5) +
      geom_errorbar(aes(x= transcript, ymin= errmin, ymax= errmax, colour= condition), width= 0.5, size= 0.5) +
      labs(title= pid, y= vals)
  } else if (style == "bars"){
    # Display as dodged bar chart.
    result <- ggplot(data= vis_data, aes(x= transcript, y= expression, fill= condition)) +
      geom_bar(aes(group= condition, colour= condition), stat= "identity", position= position_dodge(0.45), width= 0.5) +
      geom_errorbar(aes(x= transcript, ymin= errmin, ymax= errmax),  position= position_dodge(0.45), width= 0.5, size= 0.5) +
      labs(title= pid, y= vals)
  } else {
    stop("Invalid plot style.")
  }
  
  return(result)
}


#================================================================================
#' Summary of DTU calling.
#' 
#' @param dtuo a DTU object.
#' @return named numerical vector giving a quick overview of the results
#'
#'@export
dtu_summary <- function(dtuo) {
  result <- c("DTU genes" = sum(ifelse(dtuo$Genes[["DTU"]]==TRUE, 1, 0), na.rm=TRUE), 
              "non-DTU genes" = sum(ifelse(dtuo$Genes[["DTU"]]==FALSE, 1, 0), na.rm=TRUE), 
              "NA genes" = sum(ifelse(is.na(dtuo$Genes[["DTU"]]), 1, 0)),
              "DTU transcripts" = sum(ifelse(dtuo$Transcripts[["DTU"]]==TRUE, 1, 0), na.rm=TRUE), 
              "non-DTU transcripts" = sum(ifelse(dtuo$Transcripts[["DTU"]]==FALSE, 1, 0), na.rm=TRUE), 
              "NA transcripts" = sum(ifelse(is.na(dtuo$Transcripts[["DTU"]]), 1, 0)) )
  return(result)
}


#================================================================================
#' Plot DTU results from the proportions test.
#' 
#' @param dtuo A DTU object.
#' @param type Type of plot: "propVcount", "dpropVcount", "dpropVsig", "maxdprop"
#' @return a ggplot2 object. Simply display it or you can also customize it.
#' 
#' Generally uses the results of the transcript-level proportion tests.
#' 
#' @import data.table
#' @import ggplot2
#' @export
plot_overview <- function(dtuo, type="dpropVsig") {
  if (type == "propVcount") {
    result <- ggplot(data = dtuo$Transcripts, ggplot2::aes(sumA, propA, color = DTU)) +
      ggtitle("Relative abundances of transcripts") +
      labs(y = paste("Proportion in ", dtuo$Parameters$cond_A, sep=""), 
           x = paste("Cumulative gene read-count in ", dtuo$Parameters$cond_A, sep="")) +
      scale_x_continuous(trans = "log10") +
      scale_colour_manual("DTU", values = c("blue", "red")) + 
      geom_point(alpha = 0.3)
  } else if (type == "dpropVcount") {
    result <- ggplot(data = dtuo$Transcripts, ggplot2::aes(sumA, Dprop, color = DTU)) +
      ggtitle("Abundance change VS gene expression") +
      labs(y = paste("prop( ", dtuo$Parameters$cond_B, " ) - prop( ", dtuo$Parameters$cond_A, " )", sep=""), 
           x = paste("Cumulative gene read-count in ", dtuo$Parameters$cond_A, sep="")) +
      scale_x_continuous(trans="log10") +
      scale_y_continuous(breaks = seq(-1, 1, 0.2)) +
      scale_colour_manual("DTU", values = c("blue", "red")) + 
      geom_point(alpha = 0.3)
  } else if (type == "dpropVsig") {
    result <- ggplot(data = dtuo$Transcripts, ggplot2::aes(Dprop, pval_corr, colour = DTU)) +
      ggtitle("Proportion change VS significance") +
      labs(x = paste("Prop in ", dtuo$Parameters$cond_B, " - Prop in ", dtuo$Parameters$cond_A, sep=""), 
           y ="P-value") +
      geom_point(alpha = 0.3) +
      scale_x_continuous(breaks = seq(-1, 1, 0.2))
  } else if (type == "maxdprop") {
    tmp <- copy(dtuo$Transcripts)  # I don't want the intermediate calculations to modify the dtu object.
    tmp[, abma := abs(Dprop)]
    tmp <- with(tmp, data.table(aggregate(abma, by=list(parent_id), FUN = max)))
    # Also want coloured by dtu, so I need to retrieve that into a vector that matches tmp.
    setkey(tmp, Group.1)
    tmp[, dtu := dtuo$Genes[match(tmp$Group.1, dtuo$Genes[, parent_id]), dtuo$Genes$DTU | dtuo$Genes$DTU] ]
    # ok, plotting time
    result <- ggplot(data = na.omit(tmp), ggplot2::aes(x, fill=dtu)) +
      ggtitle("Distribution of largest proportion change per gene") +
      labs(x = paste("abs( Prop in ", dtuo$Parameters$cond_B, " - Prop in ", dtuo$Parameters$cond_A, " )", sep=""), y ="Number of genes") +
      geom_histogram(binwidth = 0.01, position="identity", alpha = 0.5) +
      scale_x_continuous(breaks = seq(0, 1, 0.1)) +
      scale_y_continuous(trans="sqrt")
  } else {
    stop("Unrecognized plot type!")
  }
  # Drop the added columns.
#   dtuo$Transcripts[, c("totalA", "totalB") := NULL]

  return(result)
}
