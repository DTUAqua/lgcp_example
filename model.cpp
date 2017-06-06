#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_SPARSE_MATRIX(Q0);
  DATA_SPARSE_MATRIX(I);
  DATA_FACTOR(time);
  DATA_FACTOR(gf);
  DATA_VECTOR(response);
  DATA_MATRIX(X); // Design matrix
  //DATA_VECTOR(sweptArea);      // by haul
  /* Random fields */
  PARAMETER_ARRAY(eta_density);
  PARAMETER_VECTOR(eta_nugget);
  /* 6 fixed effects x number of times */
  PARAMETER(logdelta);        // Length=2
  PARAMETER(logscale);         // Dim = c(ntime,2)
  PARAMETER(logsd_nugget);    // Length = ntime
  PARAMETER(time_corr);
  PARAMETER_VECTOR(beta);
  /* Stuff for prediction */
  DATA_INTEGER(doPredict);  // Flag
  DATA_MATRIX(Xpredict);
  DATA_FACTOR(Apredict);
  DATA_INTEGER(refindex);
  
  vector<Type> lin_predictor = X * beta;
  Type sd_nugget = exp(logsd_nugget);
  int nhaul = response.size();
  Type ans = 0;
  using namespace density;

  /* Optional: Add static field */
  PARAMETER_VECTOR(eta_static);
  PARAMETER_VECTOR(logdelta_static);
  PARAMETER_VECTOR(logscale_static);
  if(eta_static.size() > 0) {
    /* Scale parameters for fields */
    Type scale = exp(logscale_static[0]);
    /* GMRF: Sparse precision matrix */
    Eigen::SparseMatrix<Type> Q = Q0 + exp(logdelta_static[0]) * I;
    GMRF_t<Type> nldens = GMRF(Q);
    ans += SCALE(nldens, scale)(eta_static);
    // for(int i=0; i<NLEVELS(time); i++) {
    //   eta_density.col(i) -= eta_static;
    // }
  }

  /* Time covariance */
  //N01<Type> nldens_time;
  Type phi = time_corr / sqrt(1.0 + time_corr*time_corr);
  AR1_t<N01<Type> > nldens_time = AR1(phi);
  /* Scale parameters for fields */
  Type scale = exp(logscale);
  /* GMRF: Sparse precision matrix */
  Eigen::SparseMatrix<Type> Q = Q0 + exp(logdelta) * I;
  GMRF_t<Type> nldens = GMRF(Q);
  ans += SEPARABLE(SCALE(nldens_time, scale), nldens)(eta_density);
  /* Nugget */
  ans -= dnorm(eta_nugget, Type(0), sd_nugget, true).sum();
  /* Data for presence */
  for(int i=0; i<nhaul; i++){
    Type log_mu = eta_density(gf[i],time[i]) + lin_predictor(i) + eta_nugget(i);
    if(eta_static.size() > 0)
      log_mu += eta_static(gf[i]);
    ans -= dpois(response(i), exp(log_mu), true);
  }

  if(doPredict){
    // Create an index for each Area at each time point
    matrix<Type> index(NLEVELS(Apredict), eta_density.cols());
    index.setZero();
    matrix<Type> cell_count = index;
    matrix<Type> new_lin_predictor = Xpredict * beta.matrix();
    new_lin_predictor.resize(eta_density.rows(), eta_density.cols());
    for(int i=0; i<eta_density.rows(); i++) {
      for(int j=0; j<eta_density.cols(); j++) {
	index(Apredict(i), j) += exp(new_lin_predictor(i, j) + eta_density(i, j));
	cell_count(Apredict(i), j) += 1.0;
      }
    }
    index.array() /= cell_count.array();
    // Log transform index:
    matrix<Type> logindex = index.array().log();
    if(refindex > 0) {
      refindex -= 1; // R index -> C index
      for(int i=0; i<logindex.rows(); i++) {
	logindex.row(i) -= logindex.row(refindex);
      }
    }
    REPORT(logindex);
    ADREPORT(logindex);
    if(isDouble<Type>::value) {
      Type avgVariance = nldens.variance().mean() * pow(scale, 2.);
      Type SpatialSD = sqrt(avgVariance);
      Type ForecastSD = sqrt(1. - phi*phi) * SpatialSD;
      REPORT(SpatialSD);
      REPORT(ForecastSD);
    }
  } else {
    ADREPORT(beta);
  }

  return ans;
}
