setOldClass("txtProgressBar") #required so that we can use txtProgressBar in a slot

################################################################################
#' An S4 base class to represent a progress update object.
#'
#' @slot steps A dataframe containing the progress positions (10%,20% etc) and associated text
#' @slot max The maximum value of the progress positions
#' @slot on Whether to run the progress update or not, default=FALSE
ProgressUpdate <- setClass("ProgressUpdate", slots = c(steps = "data.frame", max = "numeric", on = "logical"),
                           prototype=list(on=FALSE))


setMethod(f="initialize",
          signature="ProgressUpdate",
          definition=function(.Object, ..., steps, on=FALSE)
          {
            # store the steps dataframe, and calculate the max progress value
            .Object@steps <- steps
            .Object@max <- steps[length(steps[,1]),1]
            .Object@on <- on
            return(.Object)
          })

#' Progress update generic
setGeneric(name="update",
           def=function(theObject)
           {
             standardGeneric("update")
           }
)


################################################################################
#' An S4 text-output progress update class, subclassed from ProgressUpdate
#'
#' @slot step The current position in the list of progress points
TxtProgressUpdate <- setClass("TxtProgressUpdate",
                              slots = c(step = "numeric"),
                              contains="ProgressUpdate")

  setMethod(f="initialize",
            signature="TxtProgressUpdate",
            definition=function(.Object, ...)
              {
                # call base class initialisation
                # should really be at end as just callNextMethod without assignment
                # but couldn't get it to work that way
                .Object <- callNextMethod(.Object, ..., pb=pb)
                .Object@step <- 0
                return(.Object)
              })

  setMethod(f="update",
            signature="TxtProgressUpdate",
            definition=function(theObject)
            {
              if (theObject@on)
              {
                theObject@step = theObject@step + 1
                cat(theObject@steps[theObject@step,2])
                cat(".....")
                cat(theObject@steps[theObject@step,1])
                cat("%\n")
              }
              return(theObject)
            }
  )


################################################################################
#' An S4 progress bar update class, subclassed from ProgressUpdate
#'
#' @slot pb The current progress bar
BarProgressUpdate <- setClass("BarProgressUpdate",
                              slots = c(pb = "txtProgressBar"),
                              contains="ProgressUpdate")

setMethod(f="initialize",
          signature="BarProgressUpdate",
          definition=function(.Object, ...)
          {
            # call base class initialisation (need max and on defined)
            .Object <- callNextMethod(.Object, ..., pb=pb)
            # init pb
            if (.Object@on)
            {
              .Object@pb <- txtProgressBar(min = 0, .Object@max, style = 3)
            }
            return(.Object)
          })

setMethod(f="update",
          signature="BarProgressUpdate",
          definition=function(theObject)
          {
            if (theObject@on)
            {
              # get current value of progress bar
              current <- getTxtProgressBar(theObject@pb)

              # work out index of this value in steps
              if (current == 0)
              {
                thisstep <- 0
              }
              else
              {
                thisstep <- match(current, theObject@steps[,1])
              }

              # set progress bar to next value
              setTxtProgressBar(theObject@pb, theObject@steps[thisstep+1,1])
              if (thisstep+1 == theObject@max)
              {
                close(theObject@pb)
              }
            }
            return(theObject)
          }
)
