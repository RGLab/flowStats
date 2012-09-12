singletGate <- function(x, area, height, sidescatter = NULL, lower = NULL, upper = NULL,
                        maxit = 100, nsd = 5) {
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
  
  retme <- list(indices = indices, Rsquared = R, model = model)

  class(retme) <- "singletFilter"
  retme
}
