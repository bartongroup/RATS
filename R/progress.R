setOldClass("txtProgressBar") #required so that we can use txtProgressBar in a slot

#================================================================================
#' An S4 base class to represent a progress update object.
#'
#' @slot steps A dataframe containing the progress positions (10\%,20\% etc) and associated text
#' @slot max The maximum value of the progress positions
#' @slot on Whether to run the progress update or not, default=FALSE
ProgressUpdate <- setClass("ProgressUpdate", slots = c(steps = "data.frame", max = "numeric", on = "logical"),
                           prototype=list(on=FALSE))

#--------------------------------------------------------------------------------
#' Initialisation of progress update object
#' @param .Object base class
#' @param ... unnamed arguments
#' @param steps dataframe of steps to progress through
#' @param on whether the progress updater is switched on
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

#--------------------------------------------------------------------------------
#' Insert steps generic
#' @param .Object The progress object
#' @param newsteps New steps to insert to the progress steps
setGeneric(name="insert_steps",
           def=function(.Object, newsteps)
           {
             standardGeneric("insert_steps")
           }
)

#--------------------------------------------------------------------------------
#' Insert steps into progress update object
#' @param .Object base class
#' @param newsteps dataframe of steps to insert into current steps dataframe
setMethod(f="insert_steps",
          signature="ProgressUpdate",
          definition=function(.Object, newsteps)
            {
              # insert the additional steps into our current list
              colnames(.Object@steps) <- c("step","output")
              colnames(newsteps) <- c("step", "output")
              .Object@steps <- rbind(.Object@steps, newsteps)
              .Object@steps <- .Object@steps[order(.Object@steps$step),]
              # recalc max progress value
              .Object@max <- .Object@steps[length(.Object@steps[,1]),1]
              return(.Object)
          })

#--------------------------------------------------------------------------------
#' Progress update generic
#' @param theObject The progress object
setGeneric(name="update_progress",
           def=function(theObject)
           {
             standardGeneric("update_progress")
           }
)


#================================================================================
#' An S4 text-output progress update class, subclassed from ProgressUpdate
#'
#' @slot step The current position in the list of progress points
TxtProgressUpdate <- setClass("TxtProgressUpdate",
                              slots = c(step = "numeric"),
                              contains="ProgressUpdate")

#--------------------------------------------------------------------------------
#' Initialisation of subclass
#' @param .Object base class
#' @param ... unnamed arguments
setMethod(f="initialize",
          signature="TxtProgressUpdate",
          definition=function(.Object, ...)
            {
              pb = NULL #keep R CMD check happy
              # call base class initialisation
              # should really be at end as just callNextMethod without assignment
              # but couldn't get it to work that way
              .Object <- callNextMethod(.Object, ..., pb=pb)
              .Object@step <- 0
              return(.Object)
            })

#--------------------------------------------------------------------------------
#' Update progress with next text string and percent complete figure
#' 
#' @param theObject The progress object.
#' 
setMethod(f="update_progress",
          signature="TxtProgressUpdate",
          definition=function(theObject)
          {
            if (theObject@on)
            {
              theObject@step = theObject@step + 1
              # print out <text string>.....<\% complete>
              cat(theObject@steps[theObject@step,2])
              cat(".....")
              cat(theObject@steps[theObject@step,1])
              cat("%\n")
            }
            return(theObject)
          }
)


#================================================================================
#' An S4 progress bar update class, subclassed from ProgressUpdate
#'
#' @slot pb The current progress bar
BarProgressUpdate <- setClass("BarProgressUpdate",
                              slots = c(pb = "txtProgressBar"),
                              contains="ProgressUpdate")

#--------------------------------------------------------------------------------
#' Initialisation of subclass
#' @param .Object base class
#' @param ... unnamed arguments
setMethod(f="initialize",
          signature="BarProgressUpdate",
          definition=function(.Object, ...)
          {
            pb = NULL #keep R CMD check happy
            # call base class initialisation (need max and on defined)
            .Object <- callNextMethod(.Object, ..., pb=pb)
            # init progress bar
            if (.Object@on)
            {
              .Object@pb <- txtProgressBar(min = 0, .Object@max, style = 3)
            }
            return(.Object)
          })

#--------------------------------------------------------------------------------
#' Update progress by extending progress bar
#' @param theObject the progress bar
setMethod(f="update_progress",
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
