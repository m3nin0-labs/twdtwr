#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <map>
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
int max_cycle_length(std::string cycle_length, std::string time_scale) {
  // Define the units in larger using nested std::map
  std::map<std::string, std::map<std::string, int>> units_in_larger = {
    {"year", {
      {"month", 12}, {"day", 366}, {"hour", 24 * 366}, {"minute", 60 * 24 * 366}, {"second", 60 * 60 * 24 * 366}
    }},
    {"month", {
      {"day", 31}, {"hour", 24 * 31}, {"minute", 60 * 24 * 31}, {"second", 60 * 60 * 24 * 31}
    }},
    {"day", {
      {"hour", 24}, {"minute", 24 * 60}, {"second", 24 * 60 * 60}
    }},
    {"hour", {
      {"minute", 60}, {"second", 60 * 60}
    }},
    {"minute", {
      {"second", 60}
    }}
  };
  
  // Check if the cycle_length is valid
  if (units_in_larger.find(cycle_length) == units_in_larger.end()) {
    Rcpp::stop("Invalid cycle_length");
  }
  
  // Check if the time_scale is valid for the provided cycle_length
  if (units_in_larger[cycle_length].find(time_scale) == units_in_larger[cycle_length].end()) {
    Rcpp::stop("Invalid time_scale for the provided cycle_length");
  }
  
  // Get the max_value
  return units_in_larger[cycle_length][time_scale];
}

/**
 * Calculate time weight (logistic model)
 * 
 * @description
 * This function calculates the time weight, using the logistic model described
 * proposed by Maus, et al (2016).
 *
 * @param date_diff A `NumericVector` with the difference between dates from two time-series.
 * @param alpha A `double` value representing steepness (α) of the logistic model.
 * @param beta A `double` value representing midpoint (β) of the logistic model.
 * 
 * @reference
 * Maus, Victor, Gilberto Camara, Ricardo Cartaxo, Alber Sanchez, 
 * Fernando M. Ramos, and Gilberto R. de Queiroz. 2016. “A Time-Weighted 
 * Dynamic Time Warping Method for Land-Use and Land-Cover Mapping.” IEEE 
 * Journal of Selected Topics in Applied Earth Observations and Remote Sensing 
 * 9 (8): 3729–39. https://doi.org/10.1109/JSTARS.2016.2517118.
 */
double calculate_time_weight(int cycle_len, double alpha, double beta) {
  return (1 / (1 + std::exp(-alpha * (cycle_len - beta) )));
}

/**
 * Minimum of 2 values.
 *
 * @description
 * Auxiliary function to calculate the minimum value of `x` and `y`.
 */
double minval(double x, double y)
{
  // z > nan for z != nan is required by C the standard
  int xnan = std::isnan(x), ynan = std::isnan(y);
  if(xnan || ynan) {
    if(xnan && !ynan) return y;
    if(!xnan && ynan) return x;
    return x;
  }
  return std::min(x,y);
}


/**
 * Calculate the `symmetric2` step pattern.
 *
 * @description
 * This function calculates the `symmetric2` step pattern, which uses a weight
 * of 2 for the diagonal step and 1 for the vertical and horizontal to
 * compensate for the favor of diagonal steps.
 *
 * @note
 * For more information on this step pattern, visit the `IncDTW` package
 * documentation: https://www.rdocumentation.org/packages/IncDTW/versions/1.1.4.4/topics/dtw2vec
 *
 * @reference
 * Leodolter, M., Plant, C., & Brändle, N. (2021). IncDTW: An R Package for
 * Incremental Calculation of Dynamic Time Warping. Journal of Statistical
 * Software, 99(9), 1–23. https://doi.org/10.18637/jss.v099.i09
 *
 * Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping
 * Alignments in R: The dtw Package. Journal of Statistical Software, 31(7),
 * 1–24. https://doi.org/10.18637/jss.v031.i07
 *
 * @return `symmetric2` step pattern value.
 */
double calculate_step_pattern_symmetric2(
    const double gcm10, // vertical
    const double gcm11, // diagonal
    const double gcm01, // horizontal
    const double cm00
) {
  return(cm00 + minval(gcm10, minval(cm00 + gcm11, gcm01)));
}


/**
 * Vector-based Dynamic Time Warping (DTW) distance.
 *
 * @description
 * This function calculates the Dynamic Time Warping (DTW) distance between
 * two sequences using the vector-based algorithm proposed by Leodolter
 * et al. (2021).
 *
 * The complexity of this function, as presented by Leodolter et al. (2021), is
 * equal to O(n).
 *
 * For more information on vector-based DTW, visit:
 * https://doi.org/10.18637/jss.v099.i09
 *
 * @param x A `arma::vec` with time-series values.
 * @param y A `arma::vec` with time-series values.
 *
 * @reference
 * Leodolter, M., Plant, C., & Brändle, N. (2021). IncDTW: An R Package for
 * Incremental Calculation of Dynamic Time Warping. Journal of Statistical
 * Software, 99(9), 1–23. https://doi.org/10.18637/jss.v099.i09
 *
 * @note
 * The implementation of this DTW distance calculation was adapted from the
 * `IncDTW` R package.
 *
 * @return DTW distance.
 */
// [[Rcpp::export]]
double dtw2vec(const arma::vec &x, const arma::vec &y, double time_weight)
{
  int nx = x.size();
  int ny = y.size();
  
  double *p1 = new double[nx];
  double *p2 = new double[nx];
  
  double *ptmp;
  double ret;
  
  // first column
  *p1 = std::abs(x[0] - y[0]);
  for (int i = 1; i < nx; i++)
  {
    p1[i] = std::abs(x[i] - y[0]) + p1[i - 1] + time_weight;
  }
  
  for (int j = 1; j < ny; j++)
  {
    *p2 = std::abs(x[0] - y[j]) + *(p1);
    
    for (int i = 1; i < nx; i++)
    {
      *(p2 + i) = calculate_step_pattern_symmetric2(*(p2 + i - 1), *(p1 + i - 1), *(p1 + i), std::abs(x[i] - y[j])) + time_weight;
    }
    ptmp = p1;
    p1 = p2;
    p2 = ptmp;
  }
  
  ret = *(p1 + nx - 1); // p1[nx-1]
  
  delete[] p1;
  delete[] p2;
  
  return (ret);
}


/**
 * Dynamic Time Warping (DTW) distance.
 *
 * @description
 * This function calculates the Dynamic Time Warping (DTW) distance between
 * two time-series.
 *
 * @param ts1 A `NumericVector` with time-series data.
 * @param ts2 A `NumericVector` with time-series data.
 * @param ts1_date A `DateVector` with dates from the time-series `ts1`.
 * @param ts2_date A `DateVector` with dates from the time-series `ts2`.
 * @param alpha A `double` value representing steepness (α) of the logistic model.
 * @param beta A `double` value representing midpoint (β) of the logistic model.
 *
 * @reference
 * Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping
 * Alignments in R: The dtw Package. Journal of Statistical Software, 31(7),
 * 1–24. https://doi.org/10.18637/jss.v031.i07
 * 
 * Maus, Victor, Gilberto Camara, Ricardo Cartaxo, Alber Sanchez, 
 * Fernando M. Ramos, and Gilberto R. de Queiroz. 2016. “A Time-Weighted 
 * Dynamic Time Warping Method for Land-Use and Land-Cover Mapping.” IEEE 
 * Journal of Selected Topics in Applied Earth Observations and Remote Sensing 
 * 9 (8): 3729–39. https://doi.org/10.1109/JSTARS.2016.2517118.
 *
 * @note
 * The implementation of this DTW distance calculation was adapted from the
 * `DTW_cpp` single header library (https://github.com/cjekel/DTW_cpp) to include
 * the `twdtw` implementation suggested by Maus, et al (2016).
 *
 * @return DTW distance.
 */
// [[Rcpp::export]]
double twdtw(
    const arma::vec &x, 
    const arma::vec &y,
    const std::string& cycle_length, 
    const std::string& time_scale,
    double alpha,
    double beta
)
{
  int cycle_length_value = max_cycle_length(cycle_length, time_scale);
  double time_weight = calculate_time_weight(cycle_length_value, alpha, beta);

  return (dtw2vec(x, y, time_weight));
}
