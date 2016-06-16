#================================================================================
#' Plot count or proportion DTU results for a specified gene
#'
#' @param DTUdata a DTU object with the data in to visualize.
#' @param pid a \code{parent_id} to make the plot for.
#' @param nreps the number of biological replicates per condition (this is used to calculate the error-bars for the 'proportions' plot). The default is 7.
#' @param ptype a switch for plotting either \code{"counts"} or \code{"proportion"}. The default is "counts".
#'
#' @export
plotGeneDTU <- function(DTUdata, pid, nreps=7, ptype="counts") {
  
  # first I'll slice the data to get just the transcripts I need
  vis_data <- DTUdata$Transcripts[DTUdata$Transcripts$parent_id==pid,]
  vis_data$group <- seq(1,length(vis_data$parent_id),1)
  vis_data$peA <-sqrt(vis_data$prop_A*(1-vis_data$prop_A)/nreps)
  vis_data$peB <- sqrt(vis_data$prop_B*(1-vis_data$prop_B)/nreps)
  
  # setup 'aesthetics' (seriously, WTF!)
  if (ptype=="counts") {
    a1 = aes(x=target_id, y=mean_A, colour="blue")
    a1e = aes(x=target_id, ymin=mean_A-sqrt(var_A), ymax=mean_A+sqrt(var_A), colour="blue")
    a2 = aes(x=target_id, y=mean_B, colour="red")
    a2e = aes(x=target_id, ymin=mean_B-sqrt(var_B), ymax=mean_B+sqrt(var_B), colour="red")
  } else if (ptype=="proportion") {
    a1 = aes(x=target_id, y=prop_A, colour="blue")
    a1e = aes(x=target_id, ymin=prop_A-peA, ymax=prop_A+peA, colour="blue")
    a2 = aes(x=target_id, y=prop_B, color="red")
    a2e = aes(x=target_id, ymin=prop_B-peB, ymax=prop_B+peB, color="red")
  }
  
  width <- 0.05
  ggplot(vis_data, label=group) + 
    geom_errorbar(a1e, width=width, size=1.5) +
    geom_point(a1, size=3) + 
    geom_line(a1, group=interaction("prop_A", "prop_B"), linetype=2, size=1.5) + 
    geom_errorbar(a2e, width=width, size=1.5) +
    geom_point(a2, size=3) + 
    geom_line(a2, group=interaction("prop_A", "prop_B"), linetype=2, size=1.5) + 
    theme(axis.title.x=element_blank(), text = element_text(size=20)) + 
    labs(title=pid, y=ptype) + 
    scale_color_manual("",labels = c(DTUdata$Parameters["cond_A"][[1]],
                                     DTUdata$Parameters["cond_B"][[1]]),
                       values = c("blue", "red"))
}


#================================================================================
#' Summary of DTU
#' 
#' @param dtuo a DTU object.
#' @return named numerical vector giving a quick overview of the results
#'
#'@export
dtuSummary <- function(dtuo) {
  g <- table(dtuo$Genes$dtu, useNA="always")
  p <- table(dtuo$Genes$dtu_prop, useNA="always")
  result <- c("DTU_g" = g["TRUE"], 
              "non-DTU_g" = g["FALSE"], 
              "DT_prop" = p["TRUE"], 
              "non-DTU_prop" = p["FALSE"], 
              "na" = p[3])
  return(result)
}


#================================================================================
#' Plot DTU results.
#' 
#' @param dtuo A DTU object.
#' @param type Type of plot: "propVcount", "dpropVcount", "propfoldVsig", "dpropVsig", "dprop"
#' @return a ggplot2 object. Simply display it or you can also customize it.
#' 
#' @export
plotDTU <- function(dtuo, type="propVmag") {
  if (type == "propVcount") {
    result <- ggplot2::ggplot(data = dtuo$Transcripts, ggplot2::aes(prop_A, total_A, color = dtu_prop)) +
      ggplot2::ggtitle("Relative abundances of transcripts") +
      ggplot2::labs(x = paste("Proportion in ", dtuo$Parameters$cond_A, sep=""), y = paste("Gene read-count in ", dtuo$Parameters$cond_A, sep="")) +
      ggplot2::scale_y_continuous(trans = "log10") +
      ggplot2::scale_colour_manual("DTU", values = c("blue", "red")) + 
      ggplot2::geom_point(alpha = 0.3)
  } else if (type == "dpropVcount") {
    result <- ggplot2::ggplot(data = dtuo$Transcripts, ggplot2::aes(total_A, prop_B - prop_A, color = dtu_prop)) +
      ggplot2::ggtitle("Abundance change of transcripts ") +
      ggplot2::labs(y = paste("prop( ", dtuo$Parameters$cond_B, " ) - prop( ", dtuo$Parameters$cond_A, " )", sep=""), 
                    x = paste("Gene read-count in ", dtuo$Parameters$cond_A, sep="")) +
      ggplot2::scale_x_continuous(trans="log10") +
      ggplot2::scale_y_continuous(breaks = seq(-1, 1, 0.2)) +
      ggplot2::scale_colour_manual("DTU", values = c("blue", "red")) + 
      ggplot2::geom_point(alpha = 0.3)
  } else if (type == "propfoldVsig") {
    result <- ggplot2::ggplot(data = dtuo$Transcripts, ggplot2::aes(prop_B/prop_A, pval_prop_corr, color = dtu_prop)) +
      ggplot2::ggtitle("Proportion fold-change VS significance") +
      ggplot2::labs(y = "P-value", x = paste("prop( ", dtuo$Parameters$cond_B, " ) / prop( ", dtuo$Parameters$cond_A, " )", sep="")) +
      ggplot2::geom_hline(yintercept = dtuo$Parameters$p_thresh, colour = "blue", linetype = "dotted") +
      ggplot2::scale_colour_manual("DTU", values = c("blue", "red")) + 
      ggplot2::scale_y_continuous(trans="log10") + 
      ggplot2::scale_x_continuous(trans="log2") + 
      ggplot2::geom_point(alpha = 0.3)
  } else if (type == "dpropVsig") {
    result <- ggplot2::ggplot(data = dtuo$Transcripts, ggplot2::aes(prop_B - prop_A, pval_prop_corr, colour = dtu_prop)) +
      ggplot2::ggtitle("Proportion change VS significance") +
      ggplot2::labs(x = paste("Prop in ", dtuo$Parameters$cond_B, " - Prop in ", dtuo$Parameters$cond_A, sep=""), 
                    y ="P-value") +
      ggplot2::geom_point(alpha = 0.3) +
      ggplot2::scale_x_continuous(breaks = seq(-1, 1, 0.2))
  } else if (type == "dprop") {
    tmp <- copy(dtuo$Transcripts)  # I don't want the intermediate calculattions to modify the dtu object
    tmp[, abma := abs(prop_B - prop_A)]
    tmp <- with(tmp, data.table(aggregate(abma, by=list(parent_id), FUN = max)))
    # Also want coloured by dtu, so I need to retrieve that into a vector that matches tmp.
    setkey(tmp, Group.1)
    tmp[, dtu := dtuo$Genes[tmp$Group.1, dtuo$Genes$dtu_prop]]
    # ok, plotting time
    result <- ggplot2::ggplot(data = na.omit(tmp), ggplot2::aes(x, fill=dtu)) +
      ggplot2::ggtitle("Distribution of largest proportion change per gene") +
      ggplot2::labs(x = paste("abs( Prop in ", dtuo$Parameters$cond_B, " - Prop in ", dtuo$Parameters$cond_A, " )", sep=""), y ="Number of genes") +
      ggplot2::geom_histogram(binwidth = 0.01, position="identity", alpha = 0.5) +
      #ggplot2::stat_bin(ggplot2::aes(y=..count.., label=..count..), binwidth = 0.01, geom="text") +
      ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1)) +
      ggplot2::scale_y_continuous(trans="sqrt")
  }
  return(result)
}
