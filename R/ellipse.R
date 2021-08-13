#' convert ellipse from cov/mu to points
#' used to plot priors
ellipse <- function (cov, centre, level = 0.95) 
{
  colnames(cov) <- c("x", "y")
  eg <- ellipsoidGate(.gate = cov, mean = centre, distance = sqrt(qchisq(level, 2)))
  pg <- as(eg, "polygonGate")
  pg@boundaries
  
  }