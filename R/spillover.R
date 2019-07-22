## Moved here from flowCore to eliminate dependency issues.
## The spillover and spillover_match generics remain in flowCore
## to provide a message about the move.

## ===========================================================================
## spillover method
## ---------------------------------------------------------------------------
#' Compute a spillover matrix from a flowSet
#' 
#' Spillover information for a particular experiment is often obtained by
#' running several tubes of beads or cells stained with a single color that can
#' then be used to determine a spillover matrix for use with
#' \code{\link{compensate}}.\cr\cr
#' When matching stain channels in \code{x} with the compensation controls, we
#' provide a few options. If \code{ordered}, we assume the ordering of the
#' channels in the flowSet object is the same as the ordering of the
#' compensation-control samples. If \code{regexpr}, we use a regular expression
#' to match the channel names with the names of each of the compensation control
#' \code{flowFrame}s (that is, \code{sampleNames(x)}, which will typically be the 
#' filenames passed to \code{\link{read.FCS}}).
#' By default, we must "guess" based on the largest statistic for the
#' compensation control (i.e., the row).\cr\cr
#' Additionally, matching of channels to compensation control files can
#' be accomplished using the \code{\link{spillover_match}} method, which allows
#' the matches to be specified using a csv file. The \link{flowSet} returned
#' by the \code{spillover_match} method should then be used as the \code{x} argument
#' to \code{spillover} with \code{prematched = TRUE}.
#' 
#' 
#' The algorithm used is fairly simple. First, using the scatter parameters, we
#' restrict ourselves to the most closely clustered population to reduce the
#' amount of debris. The selected statistic is then calculated on all
#' appropriate parameters and the unstained values swept out of the matrix.
#' Every sample is then normalized to [0,1] with respect to the maximum value
#' of the sample, giving the spillover in terms of a proportion of the primary
#' channel intensity.
#' @name spillover-flowSet
#' @aliases spillover spillover,flowSet-method
#' 
#' @usage 
#' \S4method{spillover}{flowSet}(x, unstained = NULL, fsc = "FSC-A", 
#'                               ssc = "SSC-A", patt = NULL, method = "median", 
#'                               stain_match = c("intensity", "ordered", "regexpr"),
#'                               useNormFilt=FALSE, prematched = FALSE)
#' 
#' @param x A flowSet of compensation beads or cells
#' @param unstained The name or index of the unstained negative control
#' @param patt An optional regular expression defining which parameters should
#' be considered
#' @param fsc The name or index of the forward scatter parameter
#' @param ssc The name or index of the side scatter parameter
#' @param method The statistic to use for calculation. Traditionally, this has
#' been the median so it is the default. The mean is sometimes more stable.
#' @param stain_match Determines how the stain channels are matched with the
#' compensation controls. See details.
#' @param useNormFilt logical Indicating whether to apply a
#' \code{\link[flowStats]{norm2Filter}} first before computing the spillover
#' @param prematched a convenience argument specifying if the channels
#' have already been matched by spillover_match. This will override the
#' values of unstained and stain_match with unstained = "unstained" and
#' stain_match = "regexpr".
#' @param exact_match a \code{logical} specifying if we should use "regex" or "exact match" to match column names. 
#' The spillover_ng will pass exact_match and "regexpr" method will be over-ridden.
#' 
#' @return A matrix for each of the parameters
#' @author B. Ellis, J. Wagner
#' @seealso \code{\link{compensate}}, \code{\link{spillover_match}}
#' @references C. B. Bagwell & E. G. Adams (1993). Fluorescence spectral
#' overlap compensation for any number of flow cytometry parameters. in: Annals
#' of the New York Academy of Sciences, 677:167-184.
#' @keywords methods
#' @export
setMethod("spillover",
          signature = signature(x = "flowSet"),
          definition = function(x, unstained = NULL, fsc = "FSC-A",
                                ssc = "SSC-A", patt = NULL, method = "median",
                                stain_match = c("intensity", "ordered", "regexpr"),
                                useNormFilt = FALSE, prematched = FALSE, exact_match = FALSE) {
            if (prematched) {
              unstained = "unstained"
              stain_match <- "regexpr"
            }
            stain_match <- match.arg(stain_match)
            
            if (is.null(unstained)) {
              stop("Sorry, we don't yet support unstained cells blended ",
                   "with stained cells", " please specify the name or index of unstained sample", call. = FALSE)
            } else {
              # Check validity of fsc, ssc, and convert to numeric
              if(!(is.numeric(fsc))){
                if(!(fsc %in% colnames(x))){
                  stop(paste0("Forward scatter parameter '", fsc, "' not found in this flowSet."), call. = FALSE)
                }else{
                  fsc <- match(fsc, colnames(x))
                }
              }else if(fsc > length(colnames(x)) || fsc < 1){
                stop(paste0("Forward scatter parameter index", fsc, "' is out of bounds."), call. = FALSE)
              }
              
              if(!(is.numeric(ssc))){
                if(!(ssc %in% colnames(x))){
                  stop(paste0("Side scatter parameter '", ssc, "' not found in this flowSet."), call. = FALSE)
                }else{
                  ssc <- match(ssc, colnames(x))
                }
              }else if(ssc > length(colnames(x)) || ssc < 1){
                stop(paste0("Side scatter parameter index", ssc, "' is out of bounds."), call. = FALSE)
              }
              
              ## We often only want spillover for a subset of the columns
              allcols <- colnames(x)
              cols <- allcols[-c(fsc,ssc)]
              cols <- if (is.null(patt)) {
                cols
              } else {
                grep(patt, cols, value = TRUE)
              }
              
              
              ## There has got to be a better way of doing this...
              if (!is.numeric(unstained)) {
                unstained <- match(unstained, sampleNames(x))
                if (is.na(unstained)) {
                  stop('Baseline(unstained sample) not found in this set. Note that the default
                        sample name searched for is "unstained". If this is incorrect, please 
                        specify the name or index of unstained sample.', call. = FALSE)
                }
              }
              
              ## Check to see if the unstained sample is in the list of
              ## stains. If not, we need to add it, making it the first
              ## row and adjust the unstained index accordingly.
              ## If it is there we adjust to the appropriate index.
              
              if (useNormFilt) {
                if(require(flowStats))
                  if (is.numeric(fsc)) {
                    fsc <- allcols[fsc]
                  }
                if (is.numeric(ssc)) {
                  ssc <- allcols[ssc]
                }
                
                if (is.na(match(fsc, allcols))) {
                  stop("Could not find forward scatter parameter. ",
                       "Please set the fsc parameter", call. = FALSE)
                }
                if (is.na(match(ssc, allcols))) {
                  stop("Could not find side scatter parameter. ",
                       "Please set the ssc parameter", call. = FALSE)
                  n2f <- norm2Filter(fsc, ssc, scale.factor = 1.5)
                  x <- Subset(x, n2f)
                }
              }
              
              # Here, we match the stain channels with the compensation controls
              # if the user has specified it. Otherwise, we must "guess" below
              # based on the largest statistic for the compensation control
              # (i.e., the row).
              # If "ordered," we assume the ordering of the channels in the
              # flowSet object is the same as the ordering of the
              # compensation-control samples.
              # Another option is to use a regular expression to match the
              # channel names with the sampleNames of the flowSet.
              if (stain_match == "intensity") {
                channel_order <- NA
              } else if (stain_match == "ordered") {
                channel_order <- seq_along(sampleNames(x))[-unstained]
              } else if (stain_match == "regexpr") {
                if (exact_match) {
                  channel_order <- sapply(cols, function(y)which(y == sampleNames(x)))
                } else {
                channel_order <- sapply(cols, grep, x = sampleNames(x), fixed = TRUE)
                }
                # Clip out those channels that do not match to a name
                matched <- channel_order[sapply(channel_order, length) != 0]
                if (anyDuplicated(matched)) {
                  stop("Multiple stains match to a common compensation-control name",
                       call. = FALSE)
                }
              }
              if(length(x) - 1 != length(cols))
              {
                stop("the number of single stained samples provided in this set doesn't match to the number of stained channels!")
              }
              if (method == "mode") {
                inten <- fsApply(x, function(flow_frame) {
                  modes <- sapply(cols, function(stain) {
                    density_stain <- density(exprs(flow_frame)[, stain])
                    with(density_stain, x[which.max(y)])
                  }, USE.NAMES = TRUE)
                  modes
                })
              } else {
                inten <- fsApply(x, each_col, method)[, cols]
              }
              
              # background correction
              inten <- pmax(sweep(inten[-unstained, ], 2, inten[unstained, ]), 0)
              
              # normalize by max of each row
              inten <- sweep(inten, 1, apply(inten, 1, max), "/")
              
              # Updates the "rownames" of the intensity matrix. If the channel
              # order was not set above, then a guess is made based on the
              # largest statistic for the compensation control (i.e., the row).
              if (any(is.na(channel_order))) {
                channel_order <- apply(inten, 1, which.max)
                if (anyDuplicated(channel_order) > 0) {
                  stop("Unable to match stains with controls based on intensity: ",
                       "a single stain matches to several multiple controls. ",
                       call. = FALSE)
                }
              }else{
                # Just bump-down the channel_order to account for
                # the removal of the unstained row. 
                channel_order <- channel_order - (channel_order > unstained)
              }
              
              # Then reverse any shuffling so the rows are in the same
              # order as the channels for symmetry
              inten <- inten[order(channel_order),]
              
              rownames(inten) <- colnames(inten)
              inten
            }
          })

## =========================================================================
## spillover_match method: Match channel names to compensation control filenames
## -------------------------------------------------------------------------
#' Construct a \code{flowSet} for use with \code{spillover} by matching channel
#' names to compensation control filenames
#' 
#' Spillover information for a particular experiment is often obtained by
#' running several tubes of beads or cells stained with a single color that can
#' then be used to determine a spillover matrix for use with
#' \code{\link{compensate}}.\cr\cr
#' This method facilitates construction of a \code{flowSet} of compensation
#' control \code{flowFrame}s using a simple file linking filenames to channels.
#' This resulting \code{flowSet} can then be used with \code{\link{spillover}}
#' using the option \code{prematched = TRUE}.\cr\cr
#' Matching stain channels to compensation controls is done via a csv
#' file (\code{matchfile}) with columns 'filename' and 'channel'. The 'channel' entries should
#' exactly match the channel names in the FCS files. The 'filename' should be
#' the FCS file name of each compensation control which should also be the
#' corresponding sample name in the \code{flowSet}. There should also be one unstained
#' control with the 'channel' entry of 'unstained'.\cr\cr
#' The method also allows for \code{x} to be missing if \code{path} is provided,
#' pointing to a directory containing the control FCS files.\cr\cr
#' 
#' 
#' @name spillover_match-flowSet
#' @aliases spillover_match spillover_match,flowSet-method
#' @usage 
#' \S4method{spillover_match}{flowSet}(x, fsc = "FSC-A", ssc = "SSC-A",
#'                                     matchfile = NULL, path)
#' @param x A flowSet of compensation beads or cells
#' @param fsc The name or index of the forward scatter parameter
#' @param ssc The name or index of the side scatter parameter
#' @param matchfile The name or path of the csv file holding the compensation control file
#' to channel matching information.
#' @param path The name or path of the directory containing the control FCS files to
#' be matched to channels by matchfile.
#' 
#' @return A \code{flowSet} with the sample names of its \code{flowFrames}
#' corresponding to the channels specified by the matchfile.
#' @author B. Ellis, J. Wagner
#' @seealso \code{\link{compensate}}, \code{\link{spillover}}
#' @keywords methods

#' @export
setMethod("spillover_match",
          signature = signature(x = "flowSet"),
          definition = function(x, fsc = "FSC-A", ssc = "SSC-A", matchfile = NULL, path) {
            if(!file.exists(matchfile)){
              stop("File ", matchfile, " not found. \n 
			 		     You must provide a file that maps single-stained 
               controls to channels (including one unstained control). \n 
			 		     It should have columns 'filename' and 'channel'", call. = FALSE)
            }
            # Check validity of fsc, ssc, and convert to numeric
            if(!(is.numeric(fsc))){
              if(!(fsc %in% colnames(x))){
                stop(paste0("Forward scatter parameter '", fsc, "' not found in this flowSet."), call. = FALSE)
              }else{
                fsc <- match(fsc, colnames(x))
              }
            }else if(fsc > length(colnames(x)) || fsc < 1){
              stop(paste0("Forward scatter parameter index", fsc, "' is out of bounds."), call. = FALSE)
            }
            
            if(!(is.numeric(ssc))){
              if(!(ssc %in% colnames(x))){
                stop(paste0("Side scatter parameter '", ssc, "' not found in this flowSet."), call. = FALSE)
              }else{
                ssc <- match(ssc, colnames(x))
              }
            }else if(ssc > length(colnames(x)) || ssc < 1){
              stop(paste0("Side scatter parameter index", ssc, "' is out of bounds."), call. = FALSE)
            }
            
            match <- read.csv(matchfile, colClasses = "character")
            if(!all(sort(tolower(colnames(match)))[c(2,1)] == c("filename","channel"))){
              stop("File ", matchfile, 
                   " should map single-stained controls to channels (including one unstained control). \n
			 		     It should have columns 'filename' and 'channel'",call. = FALSE)
            }
            if (sum(tolower(match$channel)%in%"unstained")==0) {
              stop("Sorry, we don't yet support unstained cells blended ",
                   "with stained cells", " please specify the name or index of unstained sample", call. = FALSE)
            } else if (sum(tolower(match$channel)%in%"unstained")>1) {
              stop("You may only specify one unstained sample in the match file.")
            } else {
              
              unstained.match.idx <- which(tolower(match$channel)%in%"unstained")
              unstained <- as.character(match[unstained.match.idx,"filename"])
              
              allcols <- colnames(x)
              #check if the columns passed in via the match file agree with the columns in the fcs files.
              channels <- as.character(match$channel)
              if(!all(channels[-unstained.match.idx]%in%allcols)){
                stop("Not all channels in ", match, 
                     " match to channel names in the FCS files. Check your spelling and capitalization.")
              }
              if(sum(match$filename%in%unstained)!=1)
                stop("Unique baseline (unstained sample) not found in this set.", call. = FALSE)
              
              # We often only want spillover for a subset of the columns
              # Drop those columns that don't have a line in the matchfile, while maintaining original column order
              keep <- c(fsc, ssc, match(channels[-unstained.match.idx], allcols))
              x <- x[,sort(keep)]
              
              # Make sure the sample names of x are given their corresponding channel names
              # (This is the guts of the file <-> channel matching logic)
              sampleNames(x) <- sapply(sampleNames(x), 
                                       function(to_match) 
                                         as.character(match[match$filename == to_match, "channel"]))
              x
            }
          })

#' @export
#' @rdname spillover_match-flowSet
setMethod("spillover_match",
          signature = signature(x = "missing"),
          definition = function(x, fsc = "FSC-A", ssc = "SSC-A", matchfile, path){
            if (missing(path)) {
              stop("If no flowSet is provided, a path must be provided to the control FCS files", call. = FALSE)
            }
            match <- read.csv(matchfile, colClasses = "character")  
            # By default, this will give the file names to the flowSet's sample names,
            # which is what we want
            x <- read.flowSet(files = match$filename, path = path)
            spillmat <- spillover_match(x, path = path, 
                                        fsc = fsc, ssc = ssc, 
                                        matchfile = matchfile)
          })

## ================================
## Generic for new spillover method
## --------------------------------
#' @export
setGeneric("spillover_ng", function(x,...) standardGeneric("spillover_ng"))


## =========================================================================
## spillover_ng method: Simplified spillover API. Uses a matching file.
## -------------------------------------------------------------------------
#' Compute a spillover matrix from a \code{flowSet}, simplified API
#' 
#' Spillover information for a particular experiment is often obtained by
#' running several tubes of beads or cells stained with a single color that can
#' then be used to determine a spillover matrix for use with
#' \code{\link{compensate}}.\cr\cr
#' Matching stain channels to compensation controls is done via a matching csv
#' file (at the path given by \code{matchfile}) with columns 'filename' and 'channel'. 
#' The 'channel' entries should exactly match the channel names in the FCS files. 
#' The 'filename' should be the FCS file name of each compensation control which 
#' should also be the corresponding sample name in the \code{flowSet}. 
#' There should also be one unstained control with the 'channel' entry of 'unstained'.\cr\cr
#' The method also allows for \code{x} to be missing if \code{path} is provided,
#' pointing to a directory containing the control FCS files.\cr\cr
#' By default, pregating is always done on the channels using this API, and the 
#' mode of the channel is used to compute the spillover matrix. FSC and SSC channels 
#' can be provided to allow a pregating on (approximately) a population in the FSC and
#' SSC dimensions. Also by default, a \code{\link{norm2Filter}} is applied before
#' computing the spillover. These defaults can be overridden using the \code{pregate},
#' \code{method}, and \code{useNormFilt} arguments.
#' 
#' 
#' The algorithm used is fairly simple. First, using the scatter parameters, we
#' restrict ourselves to the most closely clustered population to reduce the
#' amount of debris. The selected statistic is then calculated on all
#' appropriate parameters and the unstained values swept out of the matrix.
#' Every sample is then normalized to [0,1] with respect to the maximum value
#' of the sample, giving the spillover in terms of a proportion of the primary
#' channel intensity.
#' 
#' @name spillover_ng-flowSet
#' @aliases spillover_ng spillover_ng,flowSet-method
#' @usage 
#' \S4method{spillover_ng}{flowSet}(x, fsc = "FSC-A", ssc = "SSC-A",
#'              plot = FALSE, matchfile, path,
#'              useNormFilt = TRUE, patt = NULL, pregate = TRUE, method = "mode", \dots)
#' \S4method{spillover_ng}{missing}(x, fsc = "FSC-A", ssc = "SSC-A",
#'              plot = FALSE, matchfile, path,
#'              useNormFilt = TRUE, patt = NULL, pregate = TRUE, method = "mode", \dots)                      
#' @param x A flowSet of compensation beads or cells
#' @param fsc The name or index of the forward scatter parameter
#' @param ssc The name or index of the side scatter parameter
#' @param plot logical. Plots the kernel density for each channel when
#' pregating. Displays the gate used. If \code{pregate} is set to \code{FALSE},
#' this argument is ignored.
#' @param matchfile Name of the csv file holding the compensation control file
#' to channel matching information.
#' @param method The statistic to use for calculation. Traditionally, this has
#' been the median so it is the default. The mean is sometimes more stable.
#' @param path A path to a directory containing the control files, to be used 
#' if \code{x} is not provided.
#' @param pregate logical Indicating whether to pregate using \code{link{rangeGate}}
#' before computing the spillover
#' @param useNormFilt logical Indicating whether to apply a
#' \code{\link{norm2Filter}} first before computing the spillover
#' @param patt An optional regular expression defining which parameters should
#' be considered
#' @param \dots Additional arguments passed to
#' \code{\link[flowStats]{rangeGate}}.
#' 
#' @return A matrix for each of the parameters
#' @author B. Ellis
#' @seealso \code{\link{compensate}}, \code{\link{spillover}}
#' @references C. B. Bagwell & E. G. Adams (1993). Fluorescence spectral
#' overlap compensation for any number of flow cytometry parameters. in: Annals
#' of the New York Academy of Sciences, 677:167-184.
#' @keywords methods

#' @export
setMethod("spillover_ng",
          signature = signature(x = "flowSet"),
          definition = function(x,  
                                fsc = "FSC-A",
                                ssc = "SSC-A", 
                                plot = FALSE, matchfile,
                                path, useNormFilt = TRUE, patt = NULL, pregate = TRUE,
                                method = "mode",
                                ...) {
            if(missing(matchfile)){
              stop("A matchfile must be provided to match the channels with controls",
                   call. = FALSE)
            }
            matched  <- spillover_match(x, fsc = fsc, ssc = ssc, matchfile = matchfile)
            # Allow spillover_match to do the checks. So if we're here, the match worked
            # and fsc, ssc are good.
            cols <- colnames(matched)
            cols <- cols[!(cols %in% c(fsc, ssc))]
            unstained <- match("unstained", sampleNames(matched))
            
            #channel_order <- sapply(cols, grep, x = sampleNames(matched), fixed = TRUE)
            channel_order <- sapply(cols, function(x)which(x == sampleNames(matched)))

            if (pregate) {
              if (plot) {
                oask <- devAskNewPage(TRUE)
                on.exit(devAskNewPage(oask))
              }
              
              matched_gated <- lapply(sort(channel_order), function(channel_i) {
                flow_frame <- matched[[channel_i]]
                channel_name <- cols[which(channel_order == channel_i)]
                
                # Applies rangeGate to select positive population
                gate_filter <- rangeGate(flow_frame, stain = channel_name,
                                         inBetween = TRUE, borderQuant = 0,
                                         absolute = FALSE, peakNr = 2)
                if (plot) {
                  # Plots a kernel density for the current channel
                  plot(density(exprs(flow_frame)[, channel_name]),
                       xlab = channel_name, ylab = "Density",
                       main = paste("Compensation Control:", sampleNames(matched)[channel_i]))
                  
                  # Adds a vertical line to show gate
                  cutpoint <- c(gate_filter@min, gate_filter@max)
                  cutpoint <- cutpoint[is.finite(cutpoint)]
                  abline(v = cutpoint, col = "black", lwd = 3, lty = 2)
                }
                Subset(flow_frame, gate_filter)
              })
              # Bump down the channel order to what it would be
              # after removing the unstained row
              channel_order <- channel_order - (channel_order > unstained)
              names(matched_gated) <- sampleNames(matched)[-unstained]
              matched_gated <- matched_gated[channel_order]
              matched <- rbind2(flowSet(matched_gated), matched[unstained])
              # The unstained row has now been bumped to the end
              # regardless of where it was before
              unstained <- length(matched)
            }
            
            spillmat <- spillover(matched, fsc = fsc, ssc = ssc, method = method,
                                  useNormFilt = useNormFilt, patt = patt, prematched = TRUE, exact_match = TRUE)
            spillmat
          })

#' @export
#' @rdname spillover_ng-flowSet
setMethod("spillover_ng",
          signature = signature(x = "missing"),
          definition = function(x, 
                                fsc = "FSC-A",
                                ssc = "SSC-A", 
                                plot = FALSE, matchfile, 
                                path, useNormFilt = TRUE, patt = NULL, pregate = TRUE,
                                method = "mode",
                                ...){
            if(missing(path)){
              stop("If no flowSet is provided, a path must be provided to the control FCS files.", call. = FALSE)
            }
            if(missing(matchfile)){
              stop("A matchfile must be provided to match the channels with controls",
                   call. = FALSE)
            }
            match <- read.csv(matchfile, colClasses = "character")  
            # By default, this will give the file names to the flowSet's sample names,
            # which is what we want
            x <- read.flowSet(files=match$filename, path=path)
            spillmat <- spillover_ng(x, path=path, 
                                     fsc=fsc, ssc=ssc, 
                                     plot=plot, matchfile=matchfile,
                                     useNormFilt=useNormFilt,
                                     pregate=pregate, method=method, ...)
          })

