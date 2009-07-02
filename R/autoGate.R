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
                      bwFac=1.3, filterId="defaultLymphGate", evaluate=TRUE, ...)
{
    ## some type-checking first
    flowCore:::checkClass(channels, "character", 2)
    flowCore:::checkClass(x, c("flowSet", "flowFrame"))
    flowCore:::checkClass(scale, "numeric", 1)
    flowCore:::checkClass(bwFac, "numeric", 1)
    flowCore:::checkClass(filterId, "character", 1)
    flowCore:::checkClass(evaluate, "logical", 1)
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
    ## compute the filterResult and subset only if evaluate=TRUE
    xr <- fr <- NULL
    if(evaluate){
        fr <- filter(x, bcn2g)
        xr <- Subset(x, fr)
    }
    return(list(x=xr, n2gate=bcn2g, n2gateResults=fr))
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
                           eval=TRUE)
          tmp$n2gateResults@subSet
      })

          
