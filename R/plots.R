#================================================================================
#' Plot count or proportion changes for all transcripts of a specified gene.
#'
#' @param dtuo A DTU object.
#' @param pid A \code{parent_id} to make the plot for.
#' @param plt_type Values to plot. Either "counts" or "proportion". Default "counts".
#' @param plt_style Style of plot. Either "bars" or "lines". Default "bars".
#'
#' @import data.table
#' @import ggplot2
#' @export
plot_gene <- function(dtuo, pid, plt_type= "counts", plt_style= "bar") {
  # Slice the data to get just the relevant transcripts.
  vis_data <- dtuo$Transcripts[pid, .(target_id, meanA, meanB, stdevA, stdevB, propA, propB)]
  vis_data[, peA := sqrt(propA * (1 - propA) / dtuo$Parameters[["num_replic_A"]]) ]
  vis_data[, peB := sqrt(propB * (1 - propB) / dtuo$Parameters[["num_replic_B"]]) ]
  
  # Choose values to display.
  vA = NA_real_; vB = NA_real_; eA = NA_real_; eB = NA_real_
  if (plt_type == "counts") {
    vA <- "meanA";  vB <- "meanB";  eA <- "stdevA";  eB <- "stdevB"
  } else {
    vA <- "propA";  vB <- "propB";  eA <- "peA";   eB <- "peB"
  }
  
  # Aggregate, to simplify ggplot commands.
  vis_data[, condA := dtuo$Parameters[["cond_A"]] ]  # Recycle single value.
  vis_data[, condB := dtuo$Parameters[["cond_B"]] ]
  vis_data <- data.frame("expression" = c(vis_data[[vA]], vis_data[[vB]]),
                  "error" = c(vis_data[[eA]], vis_data[[eB]]),
                  "condition" = c(vis_data[, condA], vis_data[, condB]),
                  "transcript" = vis_data[, target_id])  # Recycle vector once.
  
  if (plt_style == "lines") {
    # Display as overlapping lines (Nick's way of displaying it, but cleaned up).
    result <- ggplot(data= vis_data, aes(x= transcript, y= expression, colour= condition)) +
      geom_freqpoly(stat= "identity", aes(group= condition, colour= condition)) +
      geom_errorbar(aes(x= transcript, ymin= expression - error, ymax= expression + error, colour= condition),  width= 0.05, size= 1.5) +
      labs(title= pid, y= type)
  } else {
    # Display as dodged bar chart.
    result <- ggplot(data= vis_data, aes(x= transcript, y= expression, fill= condition)) +
      geom_bar(stat= "identity", position= "dodge", aes(group= condition, colour= condition)) +
      geom_errorbar(aes(x= transcript, ymin= expression - error, ymax= expression + error),  width= 0.5, size= 0.5) +
      labs(title= pid, y= type)
  }
}


#================================================================================
#' Summary of DTU calling.
#' 
#' @param dtuo a DTU object.
#' @return named numerical vector giving a quick overview of the results
#'
#'@export
dtu_summary <- function(dtuo) {
  g <- table(dtuo$Genes$Gt_DTU, useNA="always")
  p <- table(dtuo$Genes$Pt_DTU, useNA="always")
  result <- c("DTU_gtest" = g["TRUE"], 
              "non-DTU_gtest" = g["FALSE"], 
              "DT_proptest" = p["TRUE"], 
              "non-DTU_proptest" = p["FALSE"], 
              "NA" = p[3])
  return(result)
}


#================================================================================
#' Plot DTU results from the proportions test.
#' 
#' @param dtuo A DTU object.
#' @param type Type of plot: "propVcount", "dpropVcount", "propfoldVsig", "dpropVsig", "dprop"
#' @return a ggplot2 object. Simply display it or you can also customize it.
#' 
#' @import ggplot2
#' @export
plotDTU_prop <- function(dtuo, type="propVmag") {
  if (type == "propVcount") {
    result <- ggplot2::ggplot(data = dtuo$Transcripts, ggplot2::aes(propA, genreadsA, color = Pt_DTU)) +
      ggplot2::ggtitle("Relative abundances of transcripts") +
      ggplot2::labs(x = paste("Proportion in ", dtuo$Parameters$cond_A, sep=""), y = paste("Gene read-count in ", dtuo$Parameters$cond_A, sep="")) +
      ggplot2::scale_y_continuous(trans = "log10") +
      ggplot2::scale_colour_manual("DTU", values = c("blue", "red")) + 
      ggplot2::geom_point(alpha = 0.3)
  } else if (type == "dpropVcount") {
    result <- ggplot2::ggplot(data = dtuo$Transcripts, ggplot2::aes(genreadsA, Dprop, color = Pt_DTU)) +
      ggplot2::ggtitle("Abundance change of transcripts ") +
      ggplot2::labs(y = paste("prop( ", dtuo$Parameters$cond_B, " ) - prop( ", dtuo$Parameters$cond_A, " )", sep=""), 
                    x = paste("Gene read-count in ", dtuo$Parameters$cond_A, sep="")) +
      ggplot2::scale_x_continuous(trans="log10") +
      ggplot2::scale_y_continuous(breaks = seq(-1, 1, 0.2)) +
      ggplot2::scale_colour_manual("DTU", values = c("blue", "red")) + 
      ggplot2::geom_point(alpha = 0.3)
  } else if (type == "propfoldVsig") {
    result <- ggplot2::ggplot(data = dtuo$Transcripts, ggplot2::aes(propfold, Pt_pval_corr, color = Pt_DTU)) +
      ggplot2::ggtitle("Proportion fold-change VS significance") +
      ggplot2::labs(y = "P-value", x = paste("prop( ", dtuo$Parameters$cond_B, " ) / prop( ", dtuo$Parameters$cond_A, " )", sep="")) +
      ggplot2::geom_hline(yintercept = dtuo$Parameters$p_thresh, colour = "blue", linetype = "dotted") +
      ggplot2::scale_colour_manual("DTU", values = c("blue", "red")) + 
      ggplot2::scale_y_continuous(trans="log10") + 
      ggplot2::scale_x_continuous(trans="log2") + 
      ggplot2::geom_point(alpha = 0.3)
  } else if (type == "dpropVsig") {
    result <- ggplot2::ggplot(data = dtuo$Transcripts, ggplot2::aes(Dprop, Pt_pval_corr, colour = Pr_dtu)) +
      ggplot2::ggtitle("Proportion change VS significance") +
      ggplot2::labs(x = paste("Prop in ", dtuo$Parameters$cond_B, " - Prop in ", dtuo$Parameters$cond_A, sep=""), 
                    y ="P-value") +
      ggplot2::geom_point(alpha = 0.3) +
      ggplot2::scale_x_continuous(breaks = seq(-1, 1, 0.2))
  } else if (type == "dprop") {
    tmp <- copy(dtuo$Transcripts)  # I don't want the intermediate calculattions to modify the dtu object
    tmp[, abma := abs(Dprop)]
    tmp <- with(tmp, data.table(aggregate(abma, by=list(parent_id), FUN = max)))
    # Also want coloured by dtu, so I need to retrieve that into a vector that matches tmp.
    setkey(tmp, Group.1)
    tmp[, dtu := dtuo$Genes[tmp$Group.1, dtuo$Genes$Pt_dtu]]
    # ok, plotting time
    result <- ggplot2::ggplot(data = na.omit(tmp), ggplot2::aes(x, fill=dtu)) +
      ggplot2::ggtitle("Distribution of largest proportion change per gene") +
      ggplot2::labs(x = paste("abs( Prop in ", dtuo$Parameters$cond_B, " - Prop in ", dtuo$Parameters$cond_A, " )", sep=""), y ="Number of genes") +
      ggplot2::geom_histogram(binwidth = 0.01, position="identity", alpha = 0.5) +
      ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1)) +
      ggplot2::scale_y_continuous(trans="sqrt")
  }
  return(result)
}
