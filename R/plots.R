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
  vis_data <- as.data.frame(DTUdata$Transcripts)[DTUdata$Transcripts$parent_id==pid,]
  vis_data$group <- seq(1,length(vis_data$parent_id),1)
  vis_data$rpe <- sqrt(vis_data$ref_proportion*(1-vis_data$ref_proportion)/nreps)
  vis_data$cpe <- sqrt(vis_data$comp_proportion*(1-vis_data$comp_proportion)/nreps)
  
  # setup 'aesthetics' (seriously, WTF!)
  if (ptype=="counts") {
    a1 = aes(x=target_id, y=ref_mean, colour="blue")
    a1e = aes(x=target_id, ymin=ref_mean-sqrt(ref_variance), ymax=ref_mean+sqrt(ref_variance), colour="blue")
    a2 = aes(x=target_id, y=comp_mean, colour="red")
    a2e = aes(x=target_id, ymin=comp_mean-sqrt(comp_variance), ymax=comp_mean+sqrt(comp_variance), colour="red")
  } else if (ptype=="proportion") {
    a1 = aes(x=target_id, y=ref_proportion, colour="blue")
    a1e = aes(x=target_id, ymin=ref_proportion-rpe, ymax=ref_proportion+rpe, colour="blue")
    a2 = aes(x=target_id, y=comp_proportion, color="red")
    a2e = aes(x=target_id, ymin=comp_proportion-cpe, ymax=comp_proportion+cpe, color="red")
  }
  
  width <- 0.05
  ggplot(vis_data, label=group) + 
    geom_errorbar(a1e, width=width, size=1.5) +
    geom_point(a1, size=3) + 
    geom_line(a1, group=interaction("ref_proportion", "comp_proportion"), linetype=2, size=1.5) + 
    geom_errorbar(a2e, width=width, size=1.5) +
    geom_point(a2, size=3) + 
    geom_line(a2, group=interaction("ref_proportion", "comp_proportion"), linetype=2, size=1.5) + 
    theme(axis.title.x=element_blank(), text = element_text(size=20)) + 
    labs(title=pid, y=ptype) + 
    scale_color_manual("",labels = c(DTUdata$Comparison["reference"][[1]],
                                     DTUdata$Comparison["compared"][[1]]),
                       values = c("blue", "red"))
}