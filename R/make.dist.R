#' @importFrom stats dist


make.dist <- function (TS, method) {
  if (method == "sq-euclidean"){
      dist.mat <- dist(TS, method = "euclidean", 
                       diag = TRUE, upper = TRUE)^2
  } else if(method == "manhattan"){ 
      dist.mat <- dist(TS, method = "manhattan", 
                       diag = TRUE, upper = TRUE)
  } else if(method == "sq-chord"){ 
      dist.mat <- dist(sqrt(TS), method = "euclidean", 
                       diag = TRUE, upper = TRUE)^2
  } else{ 
    stop("unrecognised distance measure")
  }
  as.matrix(dist.mat)
}

