#include <iostream>
#include <Rcpp.h>
#include <math.h>
#include <cmath>
#include <random>
#include <bits/stdc++.h>


using namespace std;
using namespace Rcpp;

inline double Calmu(std::vector<double> priorpar, std::vector<double> x, double xbar){
  // prior: mu0, kappa0, alpha0, beta0
  return (priorpar[1] * priorpar[0] + x.size() * xbar)/(priorpar[1] + x.size());
}

inline double Calkappa(std::vector<double> priorpar, std::vector<double> x, double xbar){
  // prior: mu0, kappa0, alpha0, beta0
  return (priorpar[1] + x.size());
}

inline double Calalpha(std::vector<double> priorpar, std::vector<double> x, double xbar){
  // prior: mu0, kappa0, alpha0, beta0
  return (priorpar[2] + 0.5 * x.size());
}

inline double Calbeta(std::vector<double> priorpar, std::vector<double> x, double xbar){
  // prior: mu0, kappa0, alpha0, beta0
  double sumx2 = 0; for(auto& element : x)  sumx2 += pow(element-xbar, 2);
  return (priorpar[3] + 0.5 * sumx2 + priorpar[1] * x.size() * pow((xbar - priorpar[0]), 2)/(2 * (priorpar[1] + x.size())));
}


inline double predict(std::vector<double> priorpar, double x){
  // prior: mu0, kappa0, alpha0, beta0
  std::vector<double> xvec {x};
  return (std::tgamma(Calalpha(priorpar, xvec, x)) * std::tgamma(priorpar[2])* 
          pow(priorpar[3], priorpar[2]) / pow(Calbeta(priorpar, xvec, x), Calalpha(priorpar, xvec, x)) *
          sqrt(priorpar[1]/Calkappa(priorpar, xvec, x)));
}

inline double AlphaUpdate(double oldAlpha, int n, int nClust, std::vector<double> priorpar){
  double r1 = R::rbeta(oldAlpha + 1, n);
  double pi1 = priorpar[0] + nClust - 1;
  double pi2 = n*(priorpar[1]-log(r1));
  double pi3 = pi1/(pi1+pi2);
  double postPar1 = priorpar[0] + nClust;
  double postPar2 = priorpar[1] - log(r1);
  if (R::runif(0, 1) > pi3) postPar1 -= 1;
  double newAlpha = R::rgamma(postPar1, 1/postPar2);
  return newAlpha;
}

struct DirichletProcess{
  double alpha;
  int nCluster;
  std::vector<int> pointsPerCluster;
  std::vector<int> ClusterLabel;
  std::vector<double> alphaPrior;
  std::vector<double> parPrior;
  std::vector<double> data;
  std::vector<double> mu;
  std::vector<double> sigma;
  std::vector<double> predArray;
  
  DirichletProcess(double alpha, int nCluster,
                     std::vector<int> pointsPerCluster,
                     std::vector<int> ClusterLabel,
                     std::vector<double> alphaPrior,
                     std::vector<double> parPrior,
                     std::vector<double> data,
                     std::vector<double> mu,
                     std::vector<double> sigma,
                     std::vector<double> predArray)
    :alpha(alpha), pointsPerCluster(pointsPerCluster), ClusterLabel(ClusterLabel),
     data(data), mu(mu), sigma(sigma), predArray(predArray), alphaPrior(alphaPrior), parPrior(parPrior){}
};

DirichletProcess ClusterUpdate(DirichletProcess dp){
  for (int t=0;t<dp.data.size();t++){
    int currentLabel = dp.ClusterLabel[t]-1;
    std::vector<int> pointsPerCluster = dp.pointsPerCluster;
    std::vector<double> mu = dp.mu;
    std::vector<double> sigma = dp.sigma;

    pointsPerCluster[currentLabel] = pointsPerCluster[currentLabel] - 1;

    for (int j=0;j<pointsPerCluster.size();j++){
      if (pointsPerCluster[j] == 0){
        pointsPerCluster.erase(pointsPerCluster.begin()+j);
        mu.erase(mu.begin()+j);
        sigma.erase(sigma.begin()+j);
      }
    }

    std::vector<double> prob(pointsPerCluster.size());
    for (int j=0;j<prob.size();j++) prob[j] = R::dnorm4(dp.data[t], mu[j], sigma[j], FALSE);
    prob.push_back(dp.alpha * dp.predArray[t]);
    double sum = 0; for (auto & element:prob) sum += element;
    // Rcout << prob[0] << "\t" << prob[1] << "\t" << prob[2] << "\t" << sum << "\n";
    
    double rn = R::runif(0, sum);
    // Rcout << rn << "\t";
    double tmp=0; int newLabel;
    for (int j=0;j<prob.size();j++){
      // Rcout << j << "\t" << prob[j] << "\t";
      if ((rn > tmp) && (rn < tmp + prob[j])){
        newLabel=j; break;
      }
      tmp = tmp + prob[j];
    }
    dp.ClusterLabel[t] = newLabel+1;
    dp.nCluster = *max_element(dp.ClusterLabel.begin(), dp.ClusterLabel.end());
    dp.mu.resize(dp.nCluster); dp.sigma.resize(dp.nCluster);
    // Rcout << "\t" << dp.mu.size() << "\t" << dp.sigma.size() << "\t" << dp.nCluster << "\n";   
    for (int j=0;j<dp.nCluster;j++){
      std::vector<double> sub;
      for (int it=0;it<dp.ClusterLabel.size();it++){
        if (dp.ClusterLabel[it] == j+1) sub.push_back(dp.data[it]);
      }
      // Rcout << sub.size() << "\t";
      // for (int it=0;it<sub.size();it++) Rcout << sub[it] << " ";
      
      double subbar = std::accumulate(sub.begin(), sub.end(), 0.0)/sub.size();
      double lambda = R::rgamma(Calalpha(dp.parPrior, sub, subbar), 1/Calbeta(dp.parPrior, sub, subbar));
      double mu = R::rnorm(Calmu(dp.parPrior, sub, subbar), 1/sqrt(lambda * Calkappa(dp.parPrior, sub, subbar)));
      // Rcout << "mu:" << mu << " sigma: " << 1/sqrt(lambda) << "\t";
      // Rcout << "average: " << subbar << "\n";
      dp.mu[j] = mu; dp.sigma[j] = 1/sqrt(lambda);
    };
    dp.alpha = AlphaUpdate(dp.alpha, dp.data.size(), dp.nCluster, dp.alphaPrior);
  }
  return dp;
}

// // [[Rcpp::export]]
// NumericVector rcpp_test(NumericVector priorpar, NumericVector x)
// {
//   std::vector<double> p = as<std::vector<double>>(priorpar);
//   std::vector<double> xcenter = as<std::vector<double>>(x);
// 
//   // double xbar = std::accumulate(xcenter.begin(), xcenter.end(), 0.0)/xcenter.size();
//   // // for (double& t : xcenter) xbar += t;
//   // for(auto& element : xcenter)  element -= xbar;
//   // NumericVector res{Calmu(p, xcenter, xbar), Calkappa(p, xcenter, xbar),
//   //                   Calalpha(p, xcenter, xbar), Calbeta(p, xcenter, xbar)};
// 
//   std::vector<double> tmp {2, 4};
//   Rcout << tmp[0] << tmp[1] << "\t" << AlphaUpdate(0.17, 400, 2, tmp);
//   NumericVector res(x.size());
//   for (int t=0;t<x.size();t++){ res[t] = predict(p, x[t]);}
// 
//   return res;
// }

// [[Rcpp::export]]
NumericVector rcpp_dp(double alpha,
                      int nCluster,
                      std::vector<int> pointsPerCluster,
                      std::vector<int> ClusterLabel,
                      std::vector<double> alphaPrior,
                      std::vector<double> parPrior,
                      std::vector<double> data,
                      std::vector<double> mu,
                      std::vector<double> sigma,
                      std::vector<double> predArray,
                      int nrun){
  DirichletProcess dp(alpha, nCluster, pointsPerCluster, ClusterLabel, alphaPrior, parPrior, data, mu, sigma, predArray);
  for (int i=0;i<nrun;i++) dp = ClusterUpdate(dp);
  NumericVector res(dp.ClusterLabel.size());
  for (int i=0;i<dp.ClusterLabel.size();i++) res[i] = dp.ClusterLabel[i];
  Rcout << "\t" << dp.alpha;
  for (int i=0;i<dp.mu.size();i++) Rcout << "\t" << dp.mu[i] << "\t" << dp.sigma[i];

  return res;
}
