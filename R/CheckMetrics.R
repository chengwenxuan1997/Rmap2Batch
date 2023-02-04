# CheckMetric
# A view for potential transformation/metric
#

library(uwot)
library(igraph)
library(RColorBrewer)


# Utils functions ---------------------------------------------------------


Allmap <- function(data, k = 15,
                   metrics = c("euclidean", "cosine", "manhattan", "hamming", "correlation")){
  res <- lapply(metrics, function(m){
    umap(X = t(data), metric = m, n_neighbors = k)
  })
  names(res) <- metrics
  # res <- list()
  # res[["euclidean"]] <- umap(X = t(data), metric = "euclidean", n_neighbors = k)
  # res[["cosine"]] <- umap(X = t(data), metric = "cosine", n_neighbors = k)
  # res[["manhattan"]] <- umap(X = t(data), metric = "manhattan", n_neighbors = k)
  # res[["hamming"]] <- umap(X = t(data), metric = "hamming", n_neighbors = k)
  # res[["correlation"]] <- umap(X = t(data), metric = "correlation", n_neighbors = k)
  return(res)
}

CheckGraph <- function(data, k = 100,
                       metrics = c("euclidean", "cosine"),
                       Batch = NULL, Celltype = NULL){
  res <- lapply(metrics, function(metric){
    if (metric == "correlation"){
      data <- sweep(data, 2, colMeans(data)); m = "cosine"
    }else{
      m = metric
    }
    nn <- Seurat:::FindNeighbors(t(data), annoy.metric = m,
                                 k.param = k, compute.SNN = F, verbose = F)$nn
    nn <- as(nn, "dgTMatrix")
    nn <- data.frame("i" = nn@i+1, "j" = nn@j+1, "x" = nn@x)
    # nn <- uwot:::annoy_nn(X = t(as.matrix(data)), k =  k, metric = m, n_threads = 1)
    # nn <- reshape2::melt(nn)
    # nn <- data.frame("i" = nn$Var1, "j" = nn$Var2)
    if (!is.null(Batch) & !is.null(Celltype)){
      nn$Batch_i <- Batch[nn$i]
      nn$Batch_j <- Batch[nn$j]
      nn$Celltype_i <- Celltype[nn$i]
      nn$Celltype_j <- Celltype[nn$j]

      message(metric, ": ", round(sum(nn$Batch_i == nn$Batch_j & nn$Celltype_i == nn$Celltype_j)/nrow(nn), 3),
              ", ", round(sum(nn$Batch_i != nn$Batch_j & nn$Celltype_i == nn$Celltype_j)/nrow(nn), 3))
    }
    return(nn)
  })
  names(res) <- metrics
  return(res)
}

SparseRank <- function(mat){
  data = mat
  mat <- as(data, "TsparseMatrix")
  mat <- data.frame("i" = mat@i+1, "j" = mat@j+1, "x" = mat@x)
  mat <- split(mat, mat$j)
  mat <- lapply(mat, function(tmp){
    # logx <- log(1+tmp$x)
    # tmp$x <- pmin(rank(tmp$x) + nrow(data) - nrow(tmp),
    #               floor(logx*2000/mean(sort(a, decreasing = T)[1:100])))
    tmp$x <- rank(tmp$x) + nrow(data) - nrow(tmp)
    # tmp$x <- rank(tmp$x)*nrow(data)/nrow(tmp)
    # tmp$x <- rank(tmp$x)
    return(tmp)
  })
  mat <- do.call(rbind, mat)
  mat <- Matrix::sparseMatrix(i = mat$i, j = mat$j, x = mat$x)
  rownames(mat) <- rownames(data); colnames(mat) <- colnames(data)
  return(mat)
}

pooling.rank <- function(mtx, k = 100){
  rank.mtx <- SparseRank(mtx)
  g <- Seurat:::FindNeighbors(object = t(rank.mtx),
                              annoy.metric = "cosine",
                              k.param = k, compute.SNN = F)$nn
  pooling.mtx <- mtx %*% t(g)
  pooling.rank.mtx <- SparseRank(pooling.mtx)
}


kendallTau <- function(A, B){
  index = t(combn(length(A), 2))
  a = A[index[, 1]] - A[index[, 2]]
  b = B[index[, 1]] - B[index[, 2]]
  d = sum(a*b<0)
  return(d)
}

Getlabel <- function(df, label = "group"){
  df$x <- as.numeric(df$x)
  df$y <- as.numeric(df$y)

  data.frame(
    "x" = tapply(df$x, df[[label]], median),
    "y" = tapply(df$y, df[[label]], median),
    "label" = tapply(df[[label]], df[[label]], unique)
  )
}


