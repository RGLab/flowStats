## Moved here from flowCore to eliminate dependency issue

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
            
            channel_order <- sapply(cols, grep, x = sampleNames(matched), fixed = TRUE)
            
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
                                  useNormFilt=useNormFilt, patt = patt, prematched = TRUE)
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

