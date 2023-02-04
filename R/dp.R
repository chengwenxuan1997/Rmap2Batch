predictArray <- function(prior.par, df){
  predictiveArray <- numeric(length(df))
  priorParameters = prior.par
  for (i in seq_along(df)) {
    PosteriorParameters_calc <- postPar(prior.par, df[i])
    predictiveArray[i] <- (gamma(PosteriorParameters_calc[3])/gamma(priorParameters[3])) *
      ((priorParameters[4]^(priorParameters[3]))/PosteriorParameters_calc[4]^PosteriorParameters_calc[3]) *
      sqrt(priorParameters[2]/PosteriorParameters_calc[2])
  }
  return(predictiveArray)
}

postPar <- function(prior.par, x){
  priorParameters <- prior.par
  n.x <- length(x)
  ybar <- mean(x)
  mu0 <- priorParameters[1]
  kappa0 <- priorParameters[2]
  alpha0 <- priorParameters[3]
  beta0 <- priorParameters[4]
  mu.n <- (kappa0 * mu0 + n.x * ybar)/(kappa0 + n.x)
  kappa.n <- kappa0 + n.x
  alpha.n <- alpha0 + n.x/2
  beta.n <- beta0 + 0.5 * sum((x - ybar)^2) + kappa0 * n.x *
    (ybar - mu0)^2/(2 * (kappa0 + n.x))
  PosteriorParameters <- matrix(c(mu.n, kappa.n, alpha.n, beta.n), ncol = 4)
  return(PosteriorParameters)
}

postDraw <- function(prior.par, df){
  n.x <- length(df)
  ybar <- mean(df)
  mu0 <- prior.par[1]
  kappa0 <- prior.par[2]
  alpha0 <- prior.par[3]
  beta0 <- prior.par[4]
  mu.n <- (kappa0 * mu0 + n.x * ybar)/(kappa0 + n.x)
  kappa.n <- kappa0 + n.x
  alpha.n <- alpha0 + n.x/2
  beta.n <- beta0 + 0.5 * sum((df - ybar)^2) + kappa0 * n.x *
    (ybar - mu0)^2/(2 * (kappa0 + n.x))
  PosteriorParameters <- matrix(c(mu.n, kappa.n, alpha.n, beta.n), ncol = 4)
  lambda <- rgamma(1, alpha.n, beta.n)
  mu <- rnorm(1, mu.n, 1/sqrt(kappa.n*lambda))
  theta <- c(mu = mu, sigma = sqrt(1/lambda))
  return(theta)
}


CluterUpdate <- function(dpobj){
  for (j in 1:dpobj$n){
    currentLabel = dpobj$clusterLabels[j]
    pointsPerCluster = dpobj$pointsPerCluster
    clusterParameters = dpobj$clusterParameters

    pointsPerCluster[currentLabel] = pointsPerCluster[currentLabel]-1
    ind = which(pointsPerCluster==0)
    if (length(ind != 0)){
      clusterParameters$mu = clusterParameters$mu[-ind]
      clusterParameters$sigma = clusterParameters$sigma[-ind]
      pointsPerCluster = pointsPerCluster[-ind]
    }

    probs <- c(pointsPerCluster * dnorm(df[j], clusterParameters$mu, clusterParameters$sigma),
               dpobj$alpha * dpobj$predictiveArray[j])
    probs[is.na(probs)] <- 0
    newlabel <- sample.int(length(probs), 1, prob = probs)
    dpobj$clusterLabels[j] = newlabel
    tmp <- split(df, dpobj$clusterLabels)
    tmp <- do.call(rbind, lapply(tmp, postDraw, prior.par = dpobj$prior.par))
    dpobj$clusterParameters <- list("mu" = tmp[, "mu"], "sigma" = tmp[, "sigma"])
    if (is.na(dpobj$pointsPerCluster[newlabel])){
      pointsPerCluster[newlabel] = 1
    }else{
      pointsPerCluster[newlabel] = pointsPerCluster[newlabel]+1
    }
    dpobj$pointsPerCluster = pointsPerCluster
  }
  return(dpobj)
}

AlphaUpdate <- function(oldParam, n, nParams, priorParameters){
  x <- rbeta(1, oldParam + 1, n)
  pi1 <- priorParameters[1] + nParams - 1
  pi2 <- n * (priorParameters[2] - log(x))
  pi1 <- pi1/(pi1 + pi2)
  postParams1 <- priorParameters[1] + nParams
  postParams2 <- priorParameters[2] - log(x)
  if (runif(1) > pi1) {
    postParams1 <- postParams1 - 1
  }
  new_alpha <- rgamma(1, postParams1, postParams2)
  return(new_alpha)
}
rgamma(1, 3, 9.84)

