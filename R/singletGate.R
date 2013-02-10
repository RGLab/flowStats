#' Creates a singlet polygon gate using the prediction bands from a robust linear model
#'
#' @param prediction_level a numeric value between 0 and 1 specifying the level
#' to use for the prediction bands.
singletGate <- function(x, area, height, sidescatter = NULL, maxit = 100,
                        prediction_level = 0.99, filterId = "singlet") {
  flowCore:::checkClass(x, "flowFrame")
  flowCore:::checkClass(area, "character")
  flowCore:::checkClass(height, "character")
  if (length(area) + length(height) != 2) {
    stop("The 'area' and 'height' channels must be provided.")
  }

  # Model the forward scatter height as a function of
  # area + side scatter + side scatter / area
  # (of the models tested, this gave the best R-squared, around 0.78)
  # If no side scatter variable is provided (i.e., 'sidescatter' is NULL), then
  # we model the forward scatter height versus forward scatter area.
  x <- exprs(x[, c(area, height, sidescatter)])
  oldcols <- colnames(x)

  if (!is.null(sidescatter)) {
    newcols <- c("A", "H", "SSC")
    form <- as.formula("H ~ A + SSC + I(SSC / A)")
  } else {
    newcols <- c("A", "H")
    form <- as.formula("H ~ A")
  }
  colnames(x) <- newcols
  
  model <- rlm(form, as.data.frame(x), maxit = maxit)

  if (!model$converged) {
    warning("The IWLS algorithm employed in 'rlm' did not converge.")
  }
  est <- huber(resid(model))

  # Creates polygon gate based on the prediction bands at the minimum and maximum
  # xChannel observation using the trained robust linear model.
  x <- as.data.frame(x)
  x_extrema <- rbind(x[which.min(x$A), ], x[which.max(x$A), ])
  predictions <- predict(model, x_extrema, interval = "prediction", level = prediction_level)

  # Create a matrix of the vertices using the prediction bands at the minimum
  # and maximum values of x. The ordering matters. Must go clockwise.
  # Otherwise, the polygon is not convex and makes an X-shape.
  gate_vertices <- rbind(cbind(x_extrema$A[1], predictions[1, "lwr"]),
                         cbind(x_extrema$A[1], predictions[1, "upr"]),
                         cbind(x_extrema$A[2], predictions[2, "upr"]),
                         cbind(x_extrema$A[2], predictions[2, "lwr"]))
  colnames(gate_vertices) <- c(area, height)

  polygonGate(gate_vertices, filterId = filterId)
}
