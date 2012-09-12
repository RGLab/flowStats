#' Creates a singlet polygon gate using the prediction bands from a robust linear model
#'
#' @param prediction_level a numeric value between 0 and 1 specifying the level
#' to use for the prediction bands.
singletGate <- function(x, area, height, sidescatter = NULL, lower = NULL, upper = NULL,
                        maxit = 100, nsd = 5, prediction_level = 0.99,
                        filter_id = "singlet") {
  flowCore:::checkClass(x, "flowFrame")
  flowCore:::checkClass(area, "character")
  flowCore:::checkClass(height, "character")
  if (length(area) + length(height) != 2) {
    stop("stains must be of length 2, specifying the FSC area and FSC height channels")
  }

  if (!(is.null(lower) & is.null(upper))) {
    flowCore:::checkClass(lower, "numeric")
    flowCore:::checkClass(upper, "numeric")

    if (lower > 1 | lower < 0 | upper > 1 | upper < 0) {
      stop("lower and upper must be in the range [0,1]")
    }
    lower <- qnorm(lower)
    upper <- qnorm(upper)
  }
	
  if (!is.null(nsd)) {
    flowCore:::checkClass(nsd, "numeric")
    lower <- -abs(nsd)
    upper <- abs(nsd)
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

  # Filter outliers based on the residuals.
  # Threshold taken from lower/upper, either number of sd, or a quantile.
  indices <- findInterval(resid(model),
                          with(est, c(mu + lower * s, mu + upper * s))) == 1

  # Calculation of the model's R^2 (Rsquared) value to assess model fit.
  serr <- sum(resid(model)^2)
  sreg <- sum((fitted(model) - huber(x[, "H"])$mu)^2)
  stot <- sreg + serr
  R <- serr / stot

  # Creates polygon gate based on the prediction bands at the minimum and maximum
  # xChannel observation using the trained robust linear model.
  x <- as.data.frame(x)
  min_x <- min(x$A)
  max_x <- max(x$A)
  x_extrema <- subset(x, A == min_x | A == max_x)
  x_extrema <- x_extrema[order(x_extrema$A), ]
  predictions <- predict(model, x_extrema, interval = "prediction", level = prediction_level)

  # Create a matrix of the vertices using the prediction bands at the minimum
  # and maximum values of x. The ordering matters. Must go clockwise.
  # Otherwise, the polygon is not convex and makes an X-shape.
  gate_vertices <- rbind(cbind(min_x, predictions[1, "lwr"]),
                         cbind(min_x, predictions[1, "upr"]),
                         cbind(max_x, predictions[2, "upr"]),
                         cbind(max_x, predictions[2, "lwr"]))
  colnames(gate_vertices) <- c(area, height)

  polygon_gate <- polygonGate(gate_vertices, filterId = filter_id)
  
  retme <- list(indices = indices, Rsquared = R, model = model, gate = polygon_gate)

  class(retme) <- "singletFilter"
  retme
}
