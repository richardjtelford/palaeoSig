#' @importFrom stats dist


make.dist <- function(TS, method) {
  if (method == "sq-euclidean") {
    dist_mat <- dist(TS,
      method = "euclidean",
      diag = TRUE, upper = TRUE
    )^2
  } else if (method == "manhattan") {
    dist_mat <- dist(TS,
      method = "manhattan",
      diag = TRUE, upper = TRUE
    )
  } else if (method == "sq-chord") {
    dist_mat <- dist(sqrt(TS),
      method = "euclidean",
      diag = TRUE, upper = TRUE
    )^2
  } else {
    stop("unrecognised distance measure")
  }
  as.matrix(dist_mat)
}
