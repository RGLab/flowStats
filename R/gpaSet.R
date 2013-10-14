## This function impliments generalized Procrestes analysis (GPA) to align
## high-density regions in a flowSet. The GPA method is concerned with
## multi-dimensional normalization wheares the warping method with
## one-dimentional. 
gpaSet <- function(x, params, register="backgating", bgChannels=NULL,
                   bg=NULL, rotation.only=TRUE,
                   downweight.missingFeatures=FALSE,
                   thres.sigma=2.5,
                   show.workflow=FALSE,
                   ask=names(dev.cur())!="pdf")
{
    ###########################
    ## check valid arguments ##
    ###########################
    flowCore:::checkClass(x, "flowSet")
    flowCore:::checkClass(params, "character")
    .checkChannel(params, colnames(x))
    if (length(params) < 2) {
        stop("At least two params are required to apply multi-dimensional
              normalization \n") }
    
    flowCore:::checkClass(register, "character")
    ## what if register=="curv1Filter"
    if (is.null(bgChannels) & is.null(bg) & register=="backgating") {
       message("Argument 'bgChannels' is missing. Default with channels
               other than", params, ".\n")
       bgChannels <- setdiff(colnames(x), c(params, "time", "Time"))
    }
    else {
        flowCore:::checkClass(bgChannels, "character")
        .checkChannel(bgChannels, colnames(x))
    }
    if (!is.null(bg)) flowCore:::checkClass(bg, "data.frame")
    flowCore:::checkClass(downweight.missingFeatures, "logical")
    flowCore:::checkClass(rotation.only, "logical")
    flowCore:::checkClass(thres.sigma, "numeric")
    ######### end of checking argument ##############

    ###############################
    ## GPA algorithm starts here ##
    ###############################
    
    ## 1. get methoded (labelled) features using either backgating (for 2D)
    ##    or nnclust (for mD)
    nDim <- length(params)
    if (register=="backgating") {
        if (is.null(bg)) {
            message("Backgating ... \n")
            bg <- flowStats:::backGating(x, xy=params, channels=bgChannels) 
        }
        features <- flowStats:::idFeaturesByBackgating(bg, nDim=nDim,
                       reference.method="median",
                       plot.workflow=show.workflow, ask=ask,
                       thres.sigma=thres.sigma, lambda=0.1)
    }
    else { ## nnclus and curv1Filter for higher dimensions
        stop("gpaSet: Only Backgating method is available")
    }

    message("Procrustes analysis ... \n")
      
    ## 2. get translation matrix, rotation matrix and scalling factor
    TransMatrix <- .getTranslationMatrix(features,
                            nDim=nDim, downweight.missingFeatures)
    CenteredFeatures <- .translateFeatures(TransMatrix, features,
                            nDim=nDim, downweight.missingFeatures)
    SVD <- .getSVD(CenteredFeatures, nDim=nDim, downweight.missingFeatures)

    ## eliminate boundary points
    for (i in params) {
        bd <- flowCore::boundaryFilter(i)
        x <- flowCore::Subset(x, bd)
    }

    ## 4. alignment: translation, rotation/scaling, un-centered
    tmatrix <- .combTransMatrix(downweight.missingFeatures, TransMatrix)
      
    expData <- fsApply(x, function(y) {
        sampleName <- keyword(y)$ORIGINALGUID
        newSet <- exprs(y)[, params]
        newSet <- .translate(newSet, tmatrix[sampleName, ])
        scal   <- ifelse(is.nan(SVD[[sampleName]]$scal), 1,
                         SVD[[sampleName]]$scal)
        newSet <- scal * (newSet %*% SVD[[sampleName]]$Q)
        trans2 <- ifelse(!length(TransMatrix$transM2),
                         0, TransMatrix$transM2[[sampleName]])
        newSet <- .translate(newSet, -trans2)
        exprs(y)[, params] <- newSet
        y})
    
    exprs.min <- fsApply(expData, function(x)
                             apply(exprs(x)[, params], 2, min))
    restoreRange <- abs(apply(exprs.min, 2, mean))
    pD <- pData(parameters(expData[[1]]))
    names <- rownames(pD)[pD$name %in% params]

    ## update parameters slot: minRange, maxRange, pData, and parameters
    regSet <- fsApply(expData, function(y) {
        exprs(y)[, params] <- .translate(exprs(y)[, params], -abs(restoreRange))
        pars <- parameters(y)
        oldRange <- pData(pars)[names, c("minRange", "maxRange")]
        
        pData(pars)[names, c("minRange", "maxRange")] <-
            .getMinMaxRange(oldRange, exprs(y)[, params])
        y@parameters <- pars
        y})
    
    ## construct gpa object
    gpaObj <- list(id.feature.method=register, norm.channels=params,
                   backgating.channels=bgChannels,
                   downweight.missingFeatures=downweight.missingFeatures,
                   SVD=SVD, TransMatrix=TransMatrix,
                   Reference=features$reference,
                   Features=features[names(features)!="reference"])

    class(gpaObj) <- "GPA"

    ## wrap up the result
    #regSet <- as(expData, "flowSet")
    phenoData(regSet) <- phenoData(x)
    regSet <- regSet[sampleNames(x)]
    
    if (!is.null(attr(x, "warping")))
        attr(regSet, "warping") <- attr(x, "warping")

    attr(regSet, "GPA") <- gpaObj

    if (show.workflow)
        .plotGPAprocess(params, features, CenteredF=CenteredFeatures,
                        SVD=SVD, TransMatrix=TransMatrix,
                        before.gpa=x, after.gpa=regSet, nDim=nDim, 
                        downweight.missingFeatures=downweight.missingFeatures,
                        ask=ask)
    return(regSet)
}


print.GPA <- function(gpaObj, ...)
{
  for (i in names(gpaObj)) {
      message("\n", i, ":\n", sep="")
      print(gpaObj[[i]])
  }
  invisible()
}
  


###########################
## translate matrix x by y
###########################
.translate <- function(x, y)
{
  ## x is an m-by-n matrix and y is an n-dimentional array (a1, a2, ..., an)
  if (length(y)==1)
     y <- rep(y, ncol(x))

  if (length(y) != ncol(x)) stop("dimension does not match")
  
  I <- matrix(1, nrow=nrow(x), ncol=1)
  centered <- x - I %*% y
}


###########################################
## get the min and max range of the data ##
###########################################
.getMinMaxRange <- function(oldRange, y)
{
   newRange = data.frame(minRange=apply(y, 2, min),
                         maxRange=apply(y, 2, max))
   newRange[, "minRange" ] <-
       pmin(oldRange[, "minRange" ], newRange[, "minRange"])
   newRange[, "maxRange" ] <-
       pmax(oldRange[, "maxRange" ], newRange[, "maxRange"])
   newRange
}

#############
## .combTransMatrix
#############
.combTransMatrix <- function(downweight.missingFeatures, TransMatrix,
                             first.centered.only=FALSE)
{
  nSample <- length(TransMatrix[names(TransMatrix$transM1)!="reference"])
  if (downweight.missingFeatures & length(TransMatrix$transM2) & !first.centered.only)
      tmatrix <- do.call(rbind, TransMatrix$transM1[1:nSample]) +
                        do.call(rbind, TransMatrix$transM2[1:nSample])
  else
      tmatrix <- do.call(rbind, TransMatrix$transM1[1:nSample])
}

#################
## check channel
#################
.checkChannel<- function(ch, allch) {
       mc <- ch %in% allch
       if(!all(mc))
           stop("Invalid parameters not mathcing the flowSet:\n   ",
                 paste(ch[!mc], collapes=", "))
       return(invisible())
}

#################################################################
## GPA: get Translation matrix for all the samples and reference
#################################################################
.getTranslationMatrix <- function(features, nDim, downweight.missingFeatures)
{
  ## features: return value of idFeaturesByBackgating
  ## if downweight.missingFeatures is TRUE, then return f$reference is
  ## a list containing multiple matrices, each of which corresponds to a
  ## particular sample.
  transM1 <- lapply(features, function(x) colMeans(x[, 1:nDim], na.rm=TRUE))
  transM2 <- list()
  ## get the second translation matrix
  if (downweight.missingFeatures) {
      ## translate features
      for (i in names(features))
          features[[i]][, 1:nDim] <- .translate(features[[i]][, 1:nDim],
                                                transM1[[i]])
          
      transM2 <- lapply(features[names(features)!="reference"],
                 function(x) colMeans(x[x$bogus==FALSE, 1:nDim], na.rm=TRUE))
      refM <- lapply(features[names(features)!="reference"],
                    function(x, y)
                        colMeans(y[x$bogus==FALSE, 1:nDim], na.rm=TRUE),
                    features$reference)
      transM2$reference <- refM
  }
  
  transM <- list(transM1=transM1, transM2=transM2)
  
  return(transM)

}

#################################################################
## GPA: translate the features to center them at the origin
#################################################################
.translateFeatures <- function(TransMatrix, tfeatures,
                               nDim, downweight.missingFeatures)
{
  ## tfeatures: features
  ## TranslationMatrix: return value of .getTranslationMatrix
  ## if downweight.missingFeatures is TRUE, then return tfeaturesf$reference is
  ## a list containing multiple matrices, each of which corresponds to a
  ## particular sample.

  ## translate samples' features
  nSample <- length(tfeatures) - 1
  ## prepare the translation matrix
  tmatrix <- .combTransMatrix(downweight.missingFeatures, TransMatrix,
                              first.centered.only=FALSE)
  
  for (i in 1:nSample)
      tfeatures[[i]][, 1:nDim] <- .translate(tfeatures[[i]][, 1:nDim],
                                             tmatrix[i, ])
  ## translate reference features
  if (!downweight.missingFeatures)
      tfeatures$reference[, 1:nDim] <- .translate(tfeatures$reference[, 1:nDim],
                                            TransMatrix$transM1$reference)
  else
      tfeatures$reference <- lapply(TransMatrix$transM2$reference,
                         function(x2, y, x1) {
                             y[, 1:nDim] <- .translate(y[,1:nDim], x2+x1)
                             y
                         }, tfeatures$reference, TransMatrix$transM1$reference)

  return(tfeatures)
}

##########################################################
## GPA: get rotation matrix (Q) and scalling factor (s) ##
##########################################################
.getSVD <- function(CenteredF, nDim, downweight.missingFeatures)
{ ## CenteredF: centered features
  ## downweight.missingFeatures: if FALSE, treat the bogus features as real. if TRUE
  ##                      elimiate bogus features

  if (!downweight.missingFeatures)
      SVD <- lapply(CenteredF[names(CenteredF) != "reference"],
                function(x, y, rotation.only) {
                    flowStats::iProcrustes(x[, 1:nDim], y[, 1:nDim],
                                           rotation.only)},
                CenteredF$reference, rotation.only=TRUE)
  else {
      SVD <- list()
      for (i in names(CenteredF[names(CenteredF) != "reference"])) {
          x <- CenteredF[[i]]
          y <- CenteredF$reference[[i]]
          SVD[[i]] <- flowStats::iProcrustes(x[x$bogus==FALSE, 1:nDim],
                                             y[x$bogus==FALSE, 1:nDim],
                                             rotation.only=TRUE)
      }
  }
  return(SVD)
}


##########################################################
## use SVD and translation matrix to transform the data ##
##########################################################
.usingSVD <- function(y, SVD, tmatrix, TransMatrix=NULL) {
        if (is(y, "flowFrame")) {
            sampleName <- keyword(y)$ORIGINALGUID
            newSet <- exprs(y)[, params]
            newSet <- .translate(newSet, tmatrix[sampleName, ])
            scal   <- ifelse(is.nan(SVD[[sampleName]]$scal), 1,
                         SVD[[sampleName]]$scal)
            newSet <- scal * (newSet %*% SVD[[sampleName]]$Q)
            trans2 <- ifelse(!length(TransMatrix$transM2),
                             0, TransMatrix$transM2[[sampleName]])
            newSet <- .translate(newSet, -trans2)
            exprs(y)[, params] <- newSet
            y
        }
        else if (is.matrix(y)) {
            newSet <- .translate(y, tmatrix)
            scal <- ifelse(is.nan(SVD$scal), 1, SVD$scal)
            newSet <- scal * (newSet %*% SVD[[sampleName]]$Q)
            trans2 <- ifelse(is.null(TransMatrix),
                             0, TransMatrix$transM2)
            newSet <- .translate(newSet, -trans2)
        }
    }
##################################################################
## GPA: plot the workflow -- alignment of features and flowsets ##
##################################################################
.plotGPAprocess <- function(params, features, CenteredF, SVD, TransMatrix,
                            before.gpa, after.gpa,
                            nDim, downweight.missingFeatures,  ask=ask)
{
  par(ask=ask)
  on.exit(par(ask=FALSE))

  ## plot alignment of the features
  newF <- CenteredF[names(CenteredF)!="reference"]

  for (i in names(newF)) {
      ## transformation: rotation and rescalling
      scal <- ifelse(is.nan(SVD[[i]]$scal), 1, SVD[[i]]$scal)
      newF[[i]][, 1:nDim] <- scal *
                             (as.matrix(newF[[i]][, 1:nDim]) %*% SVD[[i]]$Q)
      ## "un-centered"
      if (downweight.missingFeatures)
          newF[[i]][, 1:nDim] <- .translate(newF[[i]][, 1:nDim],
                                     -TransMatrix$transM2[[i]])
  }
  if (downweight.missingFeatures) {
      be <- features[names(features)!="reference"]
      tmatrix <- .combTransMatrix(downweight.missingFeatures,
                                  TransMatrix, first.centered.only=TRUE)
      ## translate the original features by TransMatrix$transM1 
      for (i in names(be))
          be[[i]][, 1:nDim] <- .translate(be[[i]][, 1:nDim], tmatrix[i, ])
      be <- lapply(be, function(x) x[x$bogus==FALSE, ])
      be <- do.call(make.groups, be)
  }
  else
    be <- do.call(make.groups, CenteredF)

  af <- do.call(make.groups, newF)
  ba <- make.groups(before.GPA=be, after.GPA=af)
  names(ba)[(ncol(ba)-1):ncol(ba)] <- c("whichS", "whichG")
  mainTitle <- do.call(paste, as.list(unique(features$channel)))
  fo <- as.formula(paste(names(ba)[1], "~", names(ba[2]), " | whichG"))
  ## features: plot features -- before and after alignament                     
  print(xyplot(fo, data=ba,
        group=whichS, 
        auto.key=list(cex=0.7, title="sample", space="right"),
        main="Features alignment",
        sub="Discard bogus features"))
  
  print(xyplot(fo, data=ba,
        group=cluster, 
        auto.key=list(cex=0.7, title="cluster", space="right"),
        sub="Discard bogus features",
        main="Features alignment"))
  
  ## flowset: plot flowset -- before and after 2D normalization
  before <- as(before.gpa, "flowFrame")
  after <-  as(after.gpa, "flowFrame")
 
  fo <- as.formula(paste(params[1], "~", params[2]))
  xs <- as(list("1.before"=before, "2.after"=after), "flowSet")
  
  f <- filter(xs, curv2Filter(params, bwFac=1.4))
  print(flowViz::xyplot(fo, after.gpa, main="after WGPA"))
  print(flowViz::xyplot(fo, xs, filter=f,
                        main="Before and After WGPA",
                        sub="Aggregation of all the flowFrames",
         par.setting=list(gate=list(fill=1:5, col=3, alpha=0.2))))

  
}



