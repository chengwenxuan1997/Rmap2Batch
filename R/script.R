SparseRank <- function(mat){
  data = mat
  mat <- as(data, "TsparseMatrix")
  mat <- data.frame("i" = mat@i+1, "j" = mat@j+1, "x" = mat@x)
  mat <- split(mat, mat$j)
  mat <- lapply(mat, function(tmp){
    tmp$x <- rank(tmp$x) + nrow(data) - nrow(tmp)
    return(tmp)
  })
  mat <- do.call(rbind, mat)
  mat <- Matrix::sparseMatrix(i = mat$i, j = mat$j, x = mat$x)
  rownames(mat) <- rownames(data); colnames(mat) <- colnames(data)
  return(mat)
}

SparseRank2 <- function(mat){
  mat = SparseRank(mat)
  sa <- nrow(mat)-median(colSums(mat>0))
  # mat@x <- mat@x + sa
  mat@x[mat@x<sa] <- 0
  mat <- drop0(mat)
}

pooling.rank <- function(mtx, method = "SparseRank2", k = 20){
  rank.mtx <- SparseRank(mtx)
  g <- Seurat:::FindNeighbors(object = t(rank.mtx), 
                              annoy.metric = "cosine", 
                              k.param = k, compute.SNN = F)$nn
  pooling.mtx <- mtx %*% t(g)
  pooling.rank.mtx <- do.call(method, list("mat" = pooling.mtx))
}

Rmap2batch <- function(seu){
  mtx <- seu@assays$RNA@counts[VariableFeatures(seu), ]
  mtx <- pooling.rank(mtx, method = "SparseRank2")
  g <- CheckGraph(as.matrix(mtx1), metrics = metrics)
  g <- split(g, g$from)
  g <- lapply(g, rcpp_dp)
  return(g)
}

