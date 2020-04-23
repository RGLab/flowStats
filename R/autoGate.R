

## ===========================================================================
## curv1Filter
## ---------------------------------------------------------------------------
## This filter can hold parameters to find siginficant high density regions
## in one dimension based on Matt Wand's feature software. This generates a
## multipleFilterResult
## ---------------------------------------------------------------------------
setClass("curv1Filter",
    representation=representation(bwFac="numeric",
        gridsize="numeric"),
    contains="parameterFilter",
    prototype=list(filterId="defaultCurv1Filter",
        bwFac=1.2,
        gridsize=rep(151, 2)))

## Constructor. We allow for the following inputs:
##  bwFac is always a numeric of length 1 and gridsize is always a numeric
##     of length 2
##  x is either a character or a transformation
curv1Filter <- function(x, bwFac=1.2, gridsize=rep(401, 2),
    filterId="defaultCurv1Filter")
{
  flowCore:::checkClass(filterId, "character", 1)
  flowCore:::checkClass(bwFac, "numeric", 1)
  flowCore:::checkClass(gridsize, "numeric", 2)
  new("curv1Filter", parameters=x, bwFac=bwFac,
      gridsize=gridsize, filterId=as.character(filterId))
}

## ==========================================================================
## curv1Filter
## ---------------------------------------------------------------------------
setMethod("show",
    signature=signature(object="curv1Filter"),
    definition=function(object)
    {
      parms <- as.character(parameters(object))
      na  <-  is.na(parms)
      if(any(na))
        parms[na] <- "internal transformation"
      msg <- paste("1D curvature filter '",object@filterId,
          "' in dimension ",
          parms, "\nwith settings:",
          "\n  bwFac=", object@bwFac, "\n  gridsize=",
          paste(object@gridsize, collapse=",", sep=""),
          sep="")
      cat(msg)
      cat("\n")
      invisible(msg)
    })
## ==========================================================================
## For curv1Filters we want to strip the attributes from the subSet slot,
## i.e., the boundaries of the high density regions and the fsObj obects.
## ---------------------------------------------------------------------------
setMethod("summarizeFilter",
    signature=signature(result="multipleFilterResult",
        filter="curv1Filter"),
    definition=function(result,filter)
    {
      ret <- callNextMethod()
      ret$boundaries <- attr(result@subSet, "boundaries")
      ret$fsObj <- attr(result@subSet, "fSObj")
      return(ret)
    })

## ===========================================================================
## curv2Filter
## ---------------------------------------------------------------------------
## This filter can hold parameters to find siginficant high density regions
## in two dimensions based on Matt Wand's feature software. This generates a
## multipleFilterResult
## ---------------------------------------------------------------------------
setClass("curv2Filter",
    representation=representation(bwFac="numeric",
        gridsize="numeric"),
    contains="parameterFilter",
    prototype=list(filterId="defaultCurv2Filter",
        bwFac=1.2,
        gridsize=rep(151, 2)))

## Constructor. We allow for the following inputs:
##  bwFac is always a numeric of length 1 and gridsize is always a numeric
##     of length 2
##  x and y are characters of length 1 or a mix of characters and
##     transformations
##  x is a character of length 2 and y is missing
##  x is a list of characters and/or transformations, y is missing
curv2Filter <- function(x, y, filterId="defaultCurv2Filter",
    bwFac=1.2, gridsize=rep(151, 2))
{
  flowCore:::checkClass(filterId, "character", 1)
  flowCore:::checkClass(bwFac, "numeric", 1)
  flowCore:::checkClass(gridsize, "numeric", 2) 
  if(missing(y)) {
    if(length(x)==1)
      stop("You must specify two parameters for a curv2Filter.",
          call.=FALSE)
    if(length(x)>2)
      warning("Only using parameters '", x[1], "' and '", x[2],
          "'.", call.=FALSE)
    y=x[[2]]
    x=x[[1]]
  }
  new("curv2Filter", parameters=list(x,y), bwFac=bwFac,
      gridsize=gridsize, filterId=as.character(filterId))
}

## ==========================================================================
## curv2Filter
## ---------------------------------------------------------------------------
setMethod("show",
    signature=signature(object="curv2Filter"),
    definition=function(object)
    {
      parms <- as.character(parameters(object))
      na  <-  is.na(parms)
      if(any(na))
        parms[na] <- "internal transformation"
      msg <- paste("2D curvature filter '",
          object@filterId,"' in dimensions ",
          paste(parms, collapse=" and "),
          "\nwith settings:",
          "\n  bwFac=", object@bwFac, "\n  gridsize=",
          paste(object@gridsize, collapse=",", sep=""),
          sep="")
      cat(msg)
      cat("\n")
      invisible(msg)
    })


## ==========================================================================
## For curv2Filters we want to strip the attributes from the subSet slot,
## i.e., the boundaries of the high density regions and the fsObj obects.
## ---------------------------------------------------------------------------
setMethod("summarizeFilter",
    signature=signature(result="multipleFilterResult",
        filter="curv2Filter"),
    definition=function(result,filter)
    {
      ret <- callNextMethod()
      ret$polygons <- attr(result@subSet, "polygon")
      ret$fsObj <- attr(result@subSet, "fSObj")
      return(ret)
    })


## ==========================================================================
## curv1Filter -- this is not a logical filter so we return a vector
## of factors indicating a population. Additional information about the
## filter result (boundaries of regions, fSObj) are stored as
## attributes of the subSet vector. We use a simple naming scheme for the
## factor levels: peak n, where n is the number of populations and rest for
## all data points not within one of the high density areas
## ---------------------------------------------------------------------------
setMethod("%in%",
    signature=signature(x="flowFrame",
        table="curv1Filter"),
    definition=function(x, table)
    {
      ## We accomplish the actual filtering via Matt Wands feature
      ## software
      param <- parameters(table)
      ovalues <- exprs(x)[, param]
      bwFac <- table@bwFac
      gridsize <- table@gridsize
      
      ## drop data that has piled up on the measurement ranges
      r <- range(x)[, param, drop = FALSE]
      sel <- ovalues > r[1,] & ovalues < r[2,] & !is.na(ovalues)
      values <- ovalues[sel]
      ## Compute normal scale bandwidth (second derivative).
      st.dev <- sqrt(var(values, na.rm=TRUE))
      Q1.val <- quantile(values,1/4, na.rm=TRUE)
      Q3.val <- quantile(values,3/4, na.rm=TRUE)
      IQR.val <- (Q3.val - Q1.val)/(qnorm(3/4) - qnorm(1/4))
      bwNS <- min(st.dev,IQR.val)*(4/(7*length(values)))^(1/9)
      ## Obtain significant high curvature intervals.
      fSObj <- featureSignif(values, bw=bwFac*bwNS,
          addSignifCurv=TRUE,
          gridsize=gridsize)
      xGrid <- unlist(fSObj$fhat$x.grid)
      hiCurvIndic <- as.numeric(fSObj$curv)
      diffGrid <- diff(c(0,hiCurvIndic,0))
      lowInds <- (1:length(diffGrid))[diffGrid==1]
      uppInds <- (1:length(diffGrid))[diffGrid==-1]
      lowLims <- (xGrid[lowInds] + xGrid[lowInds-1])/2
      uppLims <- (xGrid[uppInds] + xGrid[uppInds-1])/2
      lims <- lapply(1:length(lowLims), function(i)
            c(lowLims[i],uppLims[i]))
      ## Determine filter member indicator
      indices <- rep(0, length(ovalues))
      for(i in seq(along=lims))
        indices[ovalues>=lims[[i]][1] & ovalues <= lims[[i]][2]] <- i
      result <- factor(indices)
      levels(result) <- c("rest", paste("peak", seq_along(lims)))
      attr(result,'boundaries') <- lims
      attr(result,'fSObj') <- fSObj
      result
    })



## ==========================================================================
## curv2Filter -- this is not a logical filter so we return a vector
## of factors indicating a population. We evaluate the filter usinf the
## same algorithm as for polygon gates. Additional information about the
## filter result (polygon vertices of populations, fSObj) are stored as
## attributes of the subSet vector. We use a simple naming scheme for the
## factor levels: peak n, where n is the number of populations and rest for
## all data points not within one of the high density areas
## ---------------------------------------------------------------------------
setMethod("%in%",
    signature=signature(x="flowFrame",
        table="curv2Filter"),
    definition=function(x, table)
    {
      ## We accomplish the actual filtering via Matt Wands feature
      ## software
      param <- parameters(table)
      ovalues <- exprs(x)[, param]
      bwFac <- table@bwFac
      gridsize <- table@gridsize
      ## drop data that has piled up on the measurement ranges
      r <- range(x)[, param]
      sel <- (ovalues[,1] > r[1,1] & ovalues[,1] < r[2,1] &
            ovalues[,2] > r[1,2] & ovalues[,2] < r[2,2] &
            !is.na(ovalues[,1]) & !is.na(ovalues[,2]))
      values <- ovalues[sel, ]
      ## Compute normal scale bandwidths.
      st.devs <- sqrt(apply(values, 2, var))
      Q1.vals <- apply(values, 2, quantile, 1/4)
      Q3.vals <- apply(values, 2, quantile, 3/4)
      corr.fac <- qnorm(3/4) - qnorm(1/4)
      IQR.vals <- (Q3.vals - Q1.vals)/corr.fac
      sig.hats <- apply(cbind(st.devs, IQR.vals), 1, min)
      samp.size.fac <- nrow(values)^(-1/6)
      bwNS <- samp.size.fac*sig.hats
      ## Obtain significant high curvature regions.
      fSObj <- featureSignif(values, bw=bwFac*bwNS,
          addSignifCurv=TRUE,
          gridsize=gridsize)
      contourLinesObj <- contourLines(fSObj$fhat$x[[1]],
          fSObj$fhat$x[[2]],
          fSObj$curv, levels=0.5)
      ## Determine filter member indicator
      filterInds <- rep(0,nrow(ovalues))
      for (i in seq(along=contourLinesObj)){
        vertices <- cbind(contourLinesObj[[i]]$x,
            contourLinesObj[[i]]$y)
        sel <- as.logical(flowCore:::inPolygon(ovalues,vertices))
        filterInds[sel] <- i
      }
      
      result <- factor(filterInds)
      levels(result) <-
          c("rest", paste("area", seq_along(contourLinesObj)))
      attr(result,'polygons') <- contourLinesObj
      attr(result,'fSObj') <- fSObj
      result
    })





## 2D gating with norm2Filter. One first selects a rectangle area in the
## two-dimensional space as a rough preselection and applies the norm2Filter
## only to this subset in a second step. This function is now considered
## internal, use lymphGate instead.
autoGate <- function(x, ..., scale = 2.5)
{
    ## some type-checking first
    flowCore:::checkClass(x, "flowSet")
    flowCore:::checkClass(scale, "numeric", 1)
    stains <- names(list(...))
    if(length(stains) != 2)
        stop("Only know how to deal with 2 dimensions.", call.=FALSE)
    ## construct initial rectangle gate and subset
    rectgate2 <- rectangleGate(...)
    tmp2 <- Subset(x, filter(x, rectgate2))
    ## compute norm2Filter for the rest
    bcn2g <- do.call(norm2Filter,
                     list(stains[1], stains[2], scale = scale))
    bcn2f <- filter(tmp2, bcn2g)
    ans <- Subset(tmp2, bcn2f)
    list(x = ans, n2gate = bcn2g, n2gateResults = bcn2f)
}


## A more versatile API for autoGate. 'preselection' can be one in
##   NULL:               basically a regular norm2Filter operation without any
##                       preselection
##   a character scalar: The name of one of the channels in x used for the
##                       preselection. Only positive cells in this channel
##                       will be considered to construct the rectangle gate.
##   a list: The same as for autoGate, numerics defining the initial rectangular
##                       selection
lymphGate <- function(x, channels, preselection=NULL, scale=2.5,
                      bwFac=1.3, filterId="defaultLymphGate",
                      plot=FALSE, ...)
{
    ## some type-checking first
    flowCore:::checkClass(channels, "character", 2)
    flowCore:::checkClass(x, c("flowSet", "flowFrame"))
    flowCore:::checkClass(scale, "numeric", 1)
    flowCore:::checkClass(bwFac, "numeric", 1)
    flowCore:::checkClass(filterId, "character", 1)
    bcn2g <- do.call(norm2Filter, list(channels, scale=scale,
                                       filterId=filterId))
    if(!is.null(preselection)){
        if(is.character(preselection)){
            ## preselect by a single stain
            flowCore:::checkClass(preselection, "character", 1)
            if(!preselection %in% colnames(x))
                stop(sprintf("'%s' is not a valid flow parameter in this flowSet.",
                            preselection), call.=FALSE)
            ## collapse to a single flowFrame and find most likely positive peak
            ## (essentially the one with the highest mean) after removing margin
            ## events.
            xc <- as(x, "flowFrame")
            xc <- Subset(xc, boundaryFilter(preselection))
            xcf <- filter(xc, curv1Filter(preselection, bwFac=1.3))
            xcS <- split(xc, xcf)
            xcS <- xcS[sapply(xcS, nrow)>nrow(xc)/500]
            xcMax <- Subset(tail(xcS, n=1)[[1]], boundaryFilter(channels))
            ## estimate location and variance of this subset in the two other
            ## channels and construct a rectangular preselection from that
            m <- apply(exprs(xcMax[,channels]), 2, median)
            s <- scale*apply(exprs(xcMax[,channels]), 2, mad)
            rg <- list(c(m[1]-s[1], m[1]+s[1]), c(m[2]-s[2], m[2]+s[2]))
            names(rg) <- channels
            bcrg <- rectangleGate(.gate=rg, filterId="Preselection")
        }else if(is.list(preselection)){
            ## give the preselection as an explicit rectangle
            sapply(preselection, flowCore:::checkClass, "numeric", 2)
            if(is.null(names(preselection)))
                names(preselection) <- channels
            bcrg <- rectangleGate(preselection, filterId="Preselection")
        }else stop("Invalid argument 'preselection'.", call.=FALSE)
        bcn2g <- bcn2g %subset% bcrg
        identifier(bcn2g) <- filterId
    }
    ## compute the filterResult and subset
    xr <- fr <- NULL
    fr <- filter(x, bcn2g)
    xr <- Subset(x, fr)
 
    if (plot) {
        fm <- formula(paste(sapply(channels, function(ch) paste("`", ch, "`", sep="")),
                            collapse="~"))
        print(xyplot(fm, x, filter=bcn2g))
    }
    idx <- ifelse(is.null(preselection), 1, 2)
    if(is.list(fr)){
      lapply(fr, function(result){
        details <- result@filterDetails[[idx]]
        mean <- details$center
        cov <- details$cov*details$radius
        ellipsoidGate(mean=mean, cov=cov, filterId=filterId)
      })
    }else{
      details <- fr@filterDetails[[idx]]
      mean <- details$center
      cov <- details$cov*details$radius
      ellipsoidGate(mean=mean, cov=cov, filterId=filterId)
    }
}





## ===========================================================================
## lymphFilter
## ---------------------------------------------------------------------------
## This is basically an abstraction of the lymphGate function. It allows us
## to use it as a regular gate object.
## ---------------------------------------------------------------------------
setClass("lymphFilter",
         representation=representation(preselection="character",
                                       rectDef="list", 
                                       scale="numeric",
                                       bwFac="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="defaultLymphFilter"))

## Constructor. We allow for the following inputs:
##  scale and bwFac are always numerics of length 1
##  channels is a characters of length 2
##  preselection is either a character scalar or a named list
##  of numerics
lymphFilter <- function(channels, preselection=as.character(NULL),
                        scale=2.5, bwFac=1.3, filterId="defaultLymphFilter")
{
    .Defunct()
    flowCore:::checkClass(scale, "numeric", 1)
    flowCore:::checkClass(bwFac, "numeric", 1)
    flowCore:::checkClass(filterId, "character", 1)
    flowCore:::checkClass(channels, "character", 2)
    rdef <- if(is.list(preselection)){
        tmp <- preselection
        preselection <- as.character(NULL)
        tmp} else list()
    new("lymphFilter", parameters=channels, preselection=preselection,
        scale=scale, bwFac=bwFac, filterId=filterId, rectDef=rdef)
}


setMethod("%in%",
          signature=signature("flowFrame",
                              table="lymphFilter"),
          definition=function(x, table)
      {
          pre <- if(is.null(table@preselection)) table@rectDef else table@preselection
          if(length(parameters(table)) != 2)
              stop("lymph filters require exactly two parameters.")
          tmp <- lymphGate(x, channels=parameters(table),
                           preselection=pre,
                           scale=table@scale,
                           bwFac=table@bwFac,
                           filterId=table@filterId,
                           plot=FALSE)
          tmp$n2gateResults@subSet
      })

          
## ===========================================================================
## norm2Filter
## ---------------------------------------------------------------------------
## A class to describe the fit of a bivariate normal distribution.
## Slot method is a character describing the method used to compute the
## covariance matrix, slot scale.factor holds a numeric representing the
## Mahalanobis distance. Slot transformation holds a list of length
## giving transformations, if applicable that are applied to the data
## before gating. n is the number of points used in the subsampling step.
## ---------------------------------------------------------------------------
#' Class "norm2Filter"
#' 
#' 
#' Class and constructors for a \code{\link{filter}} that fits a bivariate
#' normal distribution to a data set of paired values and selects data points
#' according to their standard deviation from the fitted distribution.
#' 
#' 
#' The filter fits a bivariate normal distribution to the data and selects all
#' events within the Mahalanobis distance multiplied by the \code{scale.factor}
#' argument. The constructor \code{norm2Filter} is a convenience function for
#' object instantiation. Evaluating a \code{curv2Filter} results in an object
#' of class \code{\link{logicalFilterResult}}. Accordingly, \code{norm2Filters}
#' can be used to subset and to split flow cytometry data sets.
#' 
#' @name norm2Filter-class
#' @aliases norm2Filter-class norm2Filter show,norm2Filter-method %in%,flowFrame,norm2Filter-method
#' @docType class
#' @usage
#' norm2Filter(x, y, method="covMcd", scale.factor=1, n=50000,
#' filterId="defaultNorm2Filter")
#' @param x,y Characters giving the names of the measurement parameter on which
#' the filter is supposed to work on. \code{y} can be missing in which case
#' \code{x} is expected to be a character vector of length 2 or a list of
#' characters.
#' @param filterId An optional parameter that sets the \code{filterId} slot of
#' this filter. The object can later be identified by this name.
#' @param scale.factor,n Numerics of length 1, used to set the
#' \code{scale.factor} and n slots of the object.
#' @param method Character in \code{covMcd} or \code{cov.rob}, used to set the
#' \code{method} slot of the object.
#' @return
#' Returns a \code{\link{norm2Filter}} object for use in filtering
#' \code{\link{flowFrame}}s or other flow cytometry objects.
#' 
#' @note
#' See the documentation in the \code{\link[flowViz:flowViz-package]{flowViz}}
#' package for plotting of \code{norm2Filters}.
#' 
#' @section Extends:
#' 
#' Class \code{"\linkS4class{parameterFilter}"}, directly.
#' 
#' Class \code{"\linkS4class{concreteFilter}"}, by class
#' \code{parameterFilter}, distance 2.
#' 
#' Class \code{"\linkS4class{filter}"}, by class \code{parameterFilter},
#' distance 3.
#' 
#' @slot method One of \code{covMcd} or \code{cov.rob}
#' defining method used for computation of covariance matrix.
#' @slot scale.factor Numeric vector giving factor of standard
#' deviations used for data selection (all points within
#' \code{scalefac} standard deviations are selected).
#' @slot n Object of class \code{"numeric"}, the number of
#' events used to compute the covariance matrix of the bivariate
#' distribution.
#' @slot filterId Object of class \code{"character"}
#' referencing the filter.
#' @slot parameters Object of class \code{"ANY"} describing
#' the parameters used to filter the \code{\link{flowFrame}} or
#' \code{\link{flowSet}}.
#' 
#' @section Objects from the Class:
#'   Objects can be created by calls of the form \code{new("norm2Filter",
#' ...)} or using the constructor \code{norm2Filter}. The constructor
#' is the recommended way.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{\%in\%}{\code{signature(x = "flowFrame", table =
#'                                   "norm2Filter")}: The workhorse used to evaluate the filter on
#'     data. This is usually not called directly by the user, but
#'     internally by calls to the \code{\link{filter}} methods. }
#'   
#'   \item{show}{\code{signature(object = "norm2Filter")}: Print
#'     information about the filter. }
#'   
#' }
#' 
#' @author F. Hahne
#' @seealso
#' 
#' \code{\link[MASS]{cov.rob}}, \code{\link[rrcov]{CovMcd}},
#' \code{\link[flowCore:filter-methods]{filter}} for evaluation of
#' \code{norm2Filters} and \code{\link{split}} and \code{\link{Subset}}for
#' splitting and subsetting of flow cytometry data sets based on that.
#' @keywords classes methods
#' @examples
#' 
#' library(flowCore)
#' ## Loading example data
#' dat <- read.FCS(system.file("extdata","0877408774.B08",
#' package="flowCore"))
#' 
#' ## Create directly. Most likely from a command line
#' norm2Filter("FSC-H", "SSC-H", filterId="myCurv2Filter")
#' 
#' ## To facilitate programmatic construction we also have the following
#' n2f <- norm2Filter(filterId="myNorm2Filter", x=list("FSC-H", "SSC-H"),
#' scale.factor=2)
#' n2f <- norm2Filter(filterId="myNorm2Filter", x=c("FSC-H", "SSC-H"),
#' scale.factor=2)
#' 
#' ## Filtering using norm2Filter
#' fres <- filter(dat, n2f)
#' fres
#' summary(fres)
#' 
#' ## The result of norm2 filtering is a logical subset
#' Subset(dat, fres)
#' 
#' ## We can also split, in which case we get those events in and those
#' ## not in the gate as separate populations
#' split(dat, fres)
#' 
#' 
setClass("norm2Filter",
         representation=representation(method="character",
                                       scale.factor="numeric",
                                       n="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="defaultNorm2Filter",
                        scale.factor=1,
                        method="covMcd",
                        n=50000))

## Constructor. We allow for the following inputs:
##  method is always a character and scale.factor and n both are always
##     numerics, all of length 1
##  x and y are characters of length 1 or a mix of characters and
##     transformations
##  x is a character of length 2 and y is missing
##  x is a list of characters and/or transformations, y is missing
norm2Filter <- function(x, y, method="covMcd", scale.factor=1,
                        n=50000, filterId="defaultNorm2Filter")
{
  flowCore:::checkClass(method, "character", 1)
  flowCore:::checkClass(scale.factor, "numeric", 1)
  flowCore:::checkClass(n, "numeric", 1)
  flowCore:::checkClass(filterId, "character", 1)
  if(missing(y)) {
    if(length(x)==1)
      stop("You must specify two parameters for a norm2 gate.")
    if(length(x)>2)
      warning("Only the first two parameters will be used.")
    y=x[[2]]
    x=x[[1]]
  }
  new("norm2Filter", parameters=c(x, y), method=method,
      scale.factor=scale.factor, filterId=filterId, n=n)
}


## ==========================================================================
## norm2Filter -- as a logical filter, this returns a logical vector.
## Essentially, the algorithm to evaluate the filter is similar to that of
## ellipsoidGates with the addition of the scalefac argument, that controls
## the cutoff in the Mahalanobis distance.
## ---------------------------------------------------------------------------
setMethod("%in%",
          signature=signature("flowFrame",
                              table="norm2Filter"),
          definition=function(x, table)
          {
            if(nrow(x)==0)
            {
              result <- as.logical(NULL)
              attr(result, 'center') <- NA
              attr(result, 'cov') <- NA
              attr(result, 'radius') <- NA
              return(result)
            }
            if(length(parameters(table)) != 2)
              stop("norm2 filters require exactly two parameters.")
            y <- exprs(x)[,parameters(table)]
            ## drop data that has piled up on the measurement ranges
            r <- range(x)[, parameters(table)]
            sel <- (y[,1] > r[1,1] & y[,1] < r[2,1] &
                      y[,2] > r[1,2] & y[,2] < r[2,2])
            values <- y[sel, ]
            if(is.na(match(table@method,c("covMcd","cov.rob"))))
              stop("Method must be either 'covMcd' or 'cov.rob'")
            cov <- switch(table@method,
                          covMcd={
                            tmp <- if(nrow(values)>table@n)
                              CovMcd(values[sample(nrow(values),
                                                   table@n),])
                            else rrcov::CovMcd(values)
                            list(center=tmp@center, cov=tmp@cov)
                          },
                          cov.rob={MASS::cov.rob(values)},
                          stop("How did you get here?")
            )
            W <- t(y)-cov$center
            ## FIXME: a long term change might be to save chol(cov$cov)
            ## rather than cov$cov in the result.  This helps in computing
            ## the gate boundaries and qr.solve above could be replaced by
            ## the equivalent of chol2inv(chol(cov$cov)).
            covsol <- qr.solve(cov$cov) %*% W
            result <- colSums(covsol * W) < table@scale.factor^2
            attr(result, 'center') <- cov$center
            attr(result, 'cov') <- cov$cov
            attr(result, 'radius') <- table@scale.factor
            result
          })

## ==========================================================================
## norm2Filter
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("show",
          signature=signature(object="norm2Filter"),
          definition=function(object)
          {
            parms <- as.character(parameters(object))
            na  <-  is.na(parms)
            if(any(na))
              parms[na] <- "internal transformation"
            cat("norm2Filter '", identifier(object),
                "' in dimensions ", sep="")
            cat(paste(parms, sep="", collapse=" and "),
                "with parameters:\n")
            cat("  method:", object@method, "\n")
            cat("  scale.factor:", object@scale.factor, "\n")
            cat("  n:", object@n, "\n")
            cat("\n")
          })

## ==========================================================================
## For a norm2Filter we want to strip things from the attributes in the
## subSet slot, i.e., the details about the fitted bivariate normal
## distribution
## ---------------------------------------------------------------------------
setMethod("summarizeFilter",
          signature=signature(result="logicalFilterResult",
                              filter="norm2Filter"),
          definition=function(result, filter)
          {
            ret <- callNextMethod()
            ret$cov <- attr(result@subSet,'cov')
            ret$center <- attr(result@subSet,'center')
            ret$radius <- attr(result@subSet,'radius')
            ret$parameters <- parameters(filter)
            return(ret)
          })