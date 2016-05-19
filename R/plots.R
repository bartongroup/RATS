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
