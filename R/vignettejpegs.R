# This is not meant to by used by RATs or users.
# But creating saving and converting the static images of the overview plots manually each time I change the specification is getting tedious. 
# Especially when accidents happen.

load("/Volumes/kfroussios/PROJECTS/rats/DE/rats_87_2.rda") # rat_ipf_87

rat_ipf_87$Transcripts[, log2FC := log(sumB/sumA)]

for (x in c("tvolcano", "gvolcano", "fcvolcano", "dprop", "maxdprop", "fcVSdprop", "reprodVSdprop", "reprod")){
  jpeg(filename= file.path("~/workspace/Rats/vignettes/figs", paste0(x, ".jpg")), 
       width=640 , height=480, quality=50, antialias= "none")
  print( plot_overview(rat_ipf_87, x) )
  dev.off()
}
