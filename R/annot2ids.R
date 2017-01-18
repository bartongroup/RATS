#================================================================================
#' Extract matching transcript and gene IDs from a GTF file.
#' 
#' This function performs no file format validation. Please use the appropriate type.
#' It is NOT compatible GFF3, as there is no rigid specification for ID attributes.
#' 
#' @param annotfile A GTF file.
#' @param transc_header The tile for the transcripts column in the output. (target_id)
#' @param gene_header The tile for the genes column in the output. (parent_id)
#' @return A data.table with two columns, matching transcript IDs to gene IDs.
#' 
#' @export
#'
annot2ids <- function(annotfile, transc_header= "target_id", gene_header= "parent_id")
{
  annot <- as.data.table(read.delim(annotfile, comment.char= "#", header= FALSE))
  annot <- annot[(grepl("transcript_id", annot[[9]]) & grepl("gene_id", annot[[9]])), ]
  t2g <- data.table("transc_id" = gsub(".*transcript_id \"?(.+?)\"?;.*", "\\1", annot[[9]]),
                    "gene_id" = gsub(".*gene_id \"?(.+?)\"?;.*", "\\1", annot[[9]]) )
  # Clean up.
  t2g <- na.omit(t2g)
  t2g <- unique(t2g)
  names(t2g) <- c(transc_header, gene_header)
  
  return(t2g)
}
