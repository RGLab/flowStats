
#' @param ... arguments to be passed to \link{tailgate}
#' @name .tailgate
#' @title tailgate
#' @inheritParams .gate_flowclust_1d 
#' 
#' @return a \code{filter} object
#' @noRd 
.gate_tail <- function(fr, pp_res = NULL, channels, ...) {
  if(length(channels) != 1)
    stop("invalid number of channels for tailgate!")
  #pps_res may contains the standardized and collapsed data and transformation
  if(isTRUE(attr(pp_res, "openCyto_preprocessing") == "standardize")){
      transformedData <- pp_res[["flow_frame"]]
     
      #if no flow data passed in, need to tranform original fr first
     if(is.null(transformedData)){
      data <- exprs(fr)[, channels]
      exprs(fr)[, channels] <- with(pp_res, (data - center)/scale) 
      transformedData <- fr
     }
     
    g <- tailgate(transformedData, channel = channels, ...)
    gate_coordinates <- c(g@min, g@max)

    # Backtransforms the gate 
    gate_coordinates <- with(pp_res, center + scale * gate_coordinates)
    gate_coordinates <- list(gate_coordinates)
    names(gate_coordinates) <- parameters(g)
    gate <- rectangleGate(gate_coordinates, filterId = g@filterId)
          
  }else
    gate <- tailgate(fr, channel = channels, ...)
  
  #carry ind with gate
#  .gateToFilterResult(fr, yChannel, gate, positive)  
  gate
  
  # If a sample has no more than 1 observation when the 'cytokine' gate is
  # attempted, the 'center' and/or 'scale' will result be NA, in which case we
  # replace the resulting NA cutpoints with the average of the remaining
  # cutpoints. If all of the cutpoints are NA, we set the mean to 0, so that
  # all of the cutpoints are 0.
#  cutpoints_unlisted <- unlist(cutpoints)
#  if (sum(!is.na(cutpoints_unlisted)) > 0) {
#    mean_cutpoints <- mean(cutpoints_unlisted, na.rm = TRUE)
#  } else {
#    mean_cutpoints <- 0
#  }
#  cutpoints <- as.list(replace(cutpoints_unlisted, is.na(cutpoints_unlisted),
#          mean_cutpoints))
  
}

.tailgate <- .gate_tail

#' Gates the tail of a density using the derivative of a kernel density estimate
#' 
#' 
#' These methods aim to set a one-dimensional gate (cutpoint) near the edge of a peak in 
#' the density specified by a channel of a \code{\linkS4class{flowFrame}} to isolate the tail population. 
#' They allow two approaches to do this, both beginning by obtaining a
#' smoothened kernel density estimate (KDE) of the original density and then utilizing either
#' its first or second derivative.
#'
#' The default behavior of the first approach, specified by \code{method = "first_deriv"}, 
#' finds valleys in the first derivative of the KDE and uses the lowest such valley
#' to place the cutpoint on the steep right shoulder of the largest peak in the original density. 
#' 
#' The default behavior of the second approach, specified by \code{method = "second_deriv"},
#' is to find peaks in the second derivative of the KDE and use the largest such peak
#' to place the cutpoint at the point on the right shoulder of the largest
#' peak in the original density where it is most rapidly flattening (the first derivative is rapidly
#' growing less negative).
#' 
#' Both approaches can be significantly modified from defaults with a number of optional
#' arguments. The \code{num_peaks} argument specifes how many peaks should be found 
#' in the smoothened KDE and \code{ref_peak} specifies around which peak the gate's 
#' cutpoint should be placed (starting from the leftmost peak). Setting the \code{side}
#' argument to "left" modifies the procedure to put the cutpoint on the left side of the 
#' reference peak to isolate a left tail. The \code{max} and \code{min} arguments allow for 
#' pre-filtering extreme values in the channel of interest (keeping only observations with 
#' channel values less than max and/or more than min). The bandwidth used for kernel density
#' estimation can be proportionally scaled using \code{adjust} (e.g. \code{adjust = 0.5} will
#' use a bandwidth that is half of the default). This allows for tuning the level of
#' smoothing applied in the estimation.
#' 
#' Lastly, the \code{tol}, \code{auto_tol}, and \code{bias} arguments allow for adjustments
#' to be made to the cutpoint that would otherwise be returned. \code{tol} provides a tolerance value
#' that the absolute value of the KDE derivative at the cutpoint must be under. If the derivative 
#' at the original cutpoint is greater than \code{tol} in magnitude, the returned cutpoint will be the first point 
#' to the right of the original cutpoint (or to the left in the case of \code{side = "left"}) with corresponding 
#' derivative within \code{tol}. Thus in practice, a smaller value for \code{tol} effectively pushes the cutpoint
#' further down the shoulder of the peak towards the flat tail. \code{tol} is set to 0.01 by default
#' but setting \code{auto_tol = TRUE} will set the tolerance to a reasonable estimate of 
#' 1\% of the maximum absolute value of the first derivative of the KDE. \code{tol} and
#' \code{auto_tol} are only used for \code{method = "first_deriv"}. Additionally, the \code{bias}
#' argument allows for directly shifting the returned cutpoint left or right.
#' 
#' It is also possible to pass additional arguments to control the calculation
#' of the derivative, which will have some effect on the resulting cutpoint determination,
#' but this should usually not be needed. By default the number of grid points for the derivative
#' calculation will be 10,000, but this can be changed with \code{num_points}. The default
#' bandwidth can also be directly adjusted with \code{bandwidth}, where the final value used
#' will be given by \code{adjust*bandwidth} 
#'
#' @name gate_tail
#' @aliases tailgate cytokine
#' @param fr a \code{flowFrame} object
#' @param channel the channel from which the cytokine gate is constructed
#' @param filterId the name of the filter
#' @param num_peaks the number of peaks expected to see. This effectively removes
#' any peaks that are artifacts of smoothing
#' @param ref_peak After \code{num_peaks} are found, this argument provides the
#' index of the reference population from which a gate will be obtained.
#' @param strict \code{logical} when the actual number of peaks detected is less than \code{ref_peak}, 
#'                               an error is reported by default. But if \code{strict} is set to FALSE, then the reference peak will be reset to the peak of the far right.      
#' @param method the method used to select the cutpoint. Either "first_deriv" or "second_deriv". See details.
#' @param tol the tolerance value used to construct the cytokine gate from the
#' derivative of the kernel density estimate. See details.
#' @param auto_tol when TRUE, it tries to set the tolerance automatically. See details.
#' @param adjust the scaling adjustment applied to the bandwidth used in the
#' first derivative of the kernel density estimate
#' @param side On which side of the density do we want to gate the tail. Valid options are "left" or "right".
#' @param positive If FALSE, after finding the cutpoint, the \code{rectangleGate} returned will have bounds \code{(-Inf, a]} as opposed to \code{[a, Inf)}.
#' This argument will be ignored if \code{gate_tail} is used in a \code{\link{gatingTemplate}} with \code{\link{gt_gating}} or via \code{\link{gs_add_gating_method}}. 
#' In those cases, pop = "-" should be used instead.
#' @param min a numeric value that sets the lower boundary for data filtering
#' @param max a numeric value that sets the upper boundary for data filtering
#' @param bias a numeric value that adds a constant to the calculated cutpoint(threshold). Default is 0.
#' @param ... additional arguments used in calculating derivative. See details.
#' @return a \code{filterList} containing the gates (cutpoints) for each sample with
#' the corresponding \code{\linkS4class{rectangleGate}} objects defining the tail as the positive population.
#' @export
#' @examples
#' \dontrun{
#'  gate <- gate_tail(fr, Channel = "APC-A") # fr is a flowFrame
#' }
gate_tail <- function(fr, channel, filterId = "", num_peaks = 1,
    ref_peak = 1, strict = TRUE, tol = 1e-2, side = "right", min = NULL, max = NULL, bias = 0, positive = TRUE, ...) {
  
  side <- match.arg(side, c("right", "left"))
  if (!(is.null(min) && is.null(max))) {
    fr <- .truncate_flowframe(fr, channels = channel, min = min,
        max = max)
  }
  # cutpoint is calculated using the first derivative of the kernel density
  # estimate. 
  x <- as.vector(exprs(fr)[, channel])
  cutpoint <- .cytokine_cutpoint(x = x, num_peaks = num_peaks,
      ref_peak = ref_peak, tol = tol, side = side, strict = strict, ...)
  
  cutpoint <- cutpoint + bias
  if(positive){
    gate_coordinates <- list(c(cutpoint, Inf))
  }else{
    gate_coordinates <- list(c(-Inf, cutpoint))
  }
  names(gate_coordinates) <- channel
  rectangleGate(gate_coordinates, filterId = filterId)
  
}
#' @export
tailgate <- gate_tail



#' Constructs a cutpoint for a flowFrame by using a derivative of the kernel
#' density estimate
#'
#' We determine a gating cutpoint using either the first or second derivative of
#' the kernel density estimate (KDE) of the \code{x}.
#'
#' By default, we compute the first derivative of the kernel density estimate. 
#' Next, we determine the lowest valley from the derivative, which corresponds to the
#' density's mode for cytokines. We then contruct a gating cutpoint as the value
#' less than the tolerance value \code{tol} in magnitude and is also greater
#' than the lowest valley.
#'
#' Alternatively, if the \code{method} is selected as \code{second_deriv}, we
#' select a cutpoint from the second derivative of the KDE. Specifically, we
#' choose the cutpoint as the largest peak of the second derivative of the KDE
#' density which is greater than the reference peak.
#'
#' @rdname gate_tail
#' @param x a \code{numeric} vector used as input data
#' @param num_peaks the number of peaks expected to see. This effectively removes
#' any peaks that are artifacts of smoothing
#' @param ref_peak After \code{num_peaks} are found, this argument provides the
#' index of the reference population from which a gate will be obtained. By
#' default, the peak farthest to the left is used.
#' @param strict \code{logical} when the actual number of peaks detected is less than \code{ref_peak}. 
#'                               an error is reported by default. But if \code{strict} is set to FALSE, then the reference peak will be reset to the peak of the far right.      
#' @param method the method used to select the cutpoint. See details.
#' @param tol the tolerance value
#' @param auto_tol when TRUE, it tries to set the tolerance automatically.
#' @param adjust the scaling adjustment applied to the bandwidth used in the
#' first derivative of the kernel density estimate
#' @param plot logical specifying whether to plot the peaks found
#'  \code{'right'} (default) or \code{'left'}?
#' @param ... additional arguments passed to \code{.deriv_density}
#' @return the cutpoint along the x-axis
.cytokine_cutpoint <- function(x, num_peaks = 1, ref_peak = 1,
    method = c("first_deriv", "second_deriv"),
    tol = 1e-2, adjust = 1, side = "right", strict = TRUE, plot = FALSE, auto_tol = FALSE, ...) {
  
  method <- match.arg(method)
  peaks <- sort(openCyto:::.find_peaks(x, num_peaks = num_peaks, adjust = adjust, plot = plot)[, "x"])
  
  #update peak count since it can be less than num_peaks
  num_peaks <- length(peaks)
  
  if (ref_peak > num_peaks) {
    outFunc <- ifelse(strict, stop, warning)
    outFunc("The reference peak is larger than the number of peaks found.",
        "Setting the reference peak to 'num_peaks'...",
        call. = FALSE)
    ref_peak <- num_peaks
  }
  
  # TODO: Double-check that a cutpoint minimum found via 'first_deriv'
  # passes the second-derivative test.
  
  if (method == "first_deriv") {
    # Finds the deepest valleys from the kernel density and sorts them.
    # The number of valleys identified is determined by 'num_peaks'
    deriv_out <- .deriv_density(x = x, adjust = adjust, deriv = 1, ...)
    if(auto_tol){
      #Try to set the tolerance automatigically.
      tol = 0.01*max(abs(deriv_out$y))
    }
    if (side == "right") {
      
      deriv_valleys <- with(deriv_out, openCyto:::.find_valleys(x = x, y = y, adjust = adjust))
      deriv_valleys <- deriv_valleys[deriv_valleys > peaks[ref_peak]]
      deriv_valleys <- sort(deriv_valleys)[1]
      cutpoint <- with(deriv_out, x[x > deriv_valleys & abs(y) < tol])
      cutpoint <- cutpoint[1]
      
    } else if (side == "left") {
      
      deriv_out$y <- -deriv_out$y
      deriv_valleys <- with(deriv_out, openCyto:::.find_valleys(x = x, y = y, adjust = adjust))
      deriv_valleys <- deriv_valleys[deriv_valleys < peaks[ref_peak]]
      deriv_valleys <- sort(deriv_valleys, decreasing=TRUE)[1]
      cutpoint <- with(deriv_out, x[x < deriv_valleys & abs(y) < tol])
      cutpoint <- cutpoint[ length(cutpoint) ]
      
    } else {
      stop("Unrecognized 'side' argument (was '", side, "'.")
    }
    
  } else {
    # The cutpoint is selected as the first peak from the second derivative
    # density which is to the right of the reference peak.
    deriv_out <- .deriv_density(x = x, adjust = adjust, deriv = 2, ...)
    
    if (side == "right") {
      deriv_peaks <- with(deriv_out, openCyto:::.find_peaks(x, y, adjust = adjust)[, "x"])
      deriv_peaks <- deriv_peaks[deriv_peaks > peaks[ref_peak]]
      cutpoint <- sort(deriv_peaks)[1]
    } else if (side == "left") {
      deriv_out$y <- -deriv_out$y
      deriv_peaks <- with(deriv_out, openCyto:::.find_peaks(x, y, adjust = adjust)[, "x"])
      deriv_peaks <- deriv_peaks[deriv_peaks < peaks[ref_peak]]
      cutpoint <- sort(deriv_peaks, decreasing=TRUE)[length(deriv_peaks)]
    } else {
      stop("Unrecognized 'side' argument (was '", side, "'.")
    }
    
  }
  
  cutpoint
}

#' Constructs the derivative specified of the kernel density estimate of a
#' numeric vector
#'
#' The derivative is computed with \code{\link[feature:drvkde]{drvkde}} (feature package 1.2.13)
#'
#' For guidance on selecting the bandwidth, see this CrossValidated post:
#' \url{http://bit.ly/12LkJWz}
#' 
#' @param x numeric vector
#' @param deriv a numeric value specifying which derivative should be calculated.
#' By default, the first derivative is computed.
#' @param bandwidth the bandwidth to use in the kernel density estimate. If
#' \code{NULL} (default), the bandwidth is estimated using the plug-in estimate
#' from \code{\link[ks]{hpi}}.
#' @param adjust a numeric weight on the automatic bandwidth, analogous to the
#' \code{adjust} parameter in \code{\link{density}}
#' @param num_points the length of the derivative of the kernel density estimate
#' @param ... additional arguments passed to \code{drvkde}
#' @return list containing the derivative of the kernel density estimate
#' @importFrom ks hpi 
#' @noRd 
.deriv_density <- function(x, deriv = 1, bandwidth = NULL, adjust = 1,
    num_points = 10000, ...) {
  
  
  if (is.null(bandwidth)) {
    bandwidth <- hpi(x, deriv.order = deriv)
  }
  #we use the private version of drvkde (copied from feature package) to avoid the tcltk dependency
  deriv_x <- drvkde(x = x, drv = deriv, bandwidth = adjust * bandwidth,
      gridsize = num_points, ...)
  list(x = deriv_x$x.grid[[1]], y = deriv_x$est)
}

