#include <Rcpp.h>
#include <map>
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
double max_cycle_length(std::string cycle_length, std::string time_scale) {
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
 * Convert `NumericMatrix` to 2D `std::vector`.
 * 
 * @description
 * This function converts a `NumericMatrix` into a 2D `std::vector`.
 *
 * @param mat A `NumericMatrix` with single or multi variate time-series.
 */
std::vector<std::vector<double>> to_cpp_vector(NumericMatrix mat) {
  size_t rows = mat.nrow();
  size_t cols = mat.ncol();
  
  std::vector<std::vector<double>> result(rows, std::vector<double>(cols));
  
  for(size_t i = 0; i < rows; ++i) {
    for(size_t j = 0; j < cols; ++j) {
      result[i][j] = mat(i, j);
    }
  }
  
  return result;
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
// double calculate_time_weight(int cycle_len, double alpha, double beta) {
double calculate_time_weight(int time_index, double alpha, double beta, std::string cycle_length, std::string time_scale, const NumericVector& ts1_dates,
                             const NumericVector& ts2_dates) {
  // double max_scale = max_cycle_length(cycle_length, time_scale);

  // double date_diff = std::fabs(ts2_dates[time_index] - ts1_dates[time_index]);
  // double time_weight = std::min(date_diff, max_scale);
  double time_weight = std::fabs(ts2_dates[time_index] - ts1_dates[time_index]);

  return (1 / (1 + std::exp(-alpha * (time_weight - beta) )));
}


/**
 * Compute the time-weighted DTW distance between two time-series.
 *
 * @description
 * The `p-norm`, also known as the `Minkowski space`, is a generalized norm
 * calculation that includes several types of distances based on the value
 * of `p`.
 *
 * Common values of `p` include:
 *
 *  - `p = 1` for the Manhattan (city block) distance;
 *  - `p = 2` for the Euclidean norm (distance).
 *
 * More details about p-norms can be found on Wikipedia:
 * https://en.wikipedia.org/wiki/Norm_(mathematics)#p-norm
 *
 * @param a A `std::vector<double>` with time-series values.
 * @param b A `std::vector<double>` with time-series values.
 * @param b A `std::vector<double>` with weight values for each time-series point.
 *
 * @note
 * Both vectors `a` and `b` must have the same length.
 *
 * @note
 * The implementation of this DTW distance calculation was adapted from the
 * `DTW_cpp` single header library (https://github.com/cjekel/DTW_cpp).
 *
 * @return The `p-norm` value between vectors `a` and `b`.
 */
double p_norm(std::vector<double> a, std::vector<double> b)
{
  double d = 0;
  
  size_t index;
  size_t a_size = a.size();
  
  for (index = 0; index < a_size; index++)
  {
    d += std::pow(std::abs(a[index] - b[index]), 2);
  }
  
  return std::pow(d, 1.0 / 2);
}


/**
 * Time-Weighted Dynamic Time Warping (TWDTW) distance.
 *
 * @description
 * This function calculates the Time-Weighted Dynamic Time Warping (TWDTW) 
 * distance between two time-series.
 *
 * @param x A `std::vector<std::vector<double>>` with time-series values.
 * @param y A `std::vector<std::vector<double>>` with time-series values.
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
 * `DTW_cpp` single header library (https://github.com/cjekel/DTW_cpp).
 *
 * @return DTW distance.
 */
double distance_dtw_op(std::vector<std::vector<double>> a,
                       std::vector<std::vector<double>> b,
                       const NumericVector& ts1_dates,
                       const NumericVector& ts2_dates,
                       double alpha,
                       double beta,
                       const std::string cycle_length, 
                       const std::string time_scale)
{
  int n = a.size();
  int o = b.size();
  
  std::vector<std::vector<double>> d(n, std::vector<double>(o, 0.0));
  
  d[0][0] = p_norm(a[0], b[0]);
  
  for (int i = 1; i < n; i++)
  {
    d[i][0] = d[i - 1][0] + p_norm(a[i], b[0]) + calculate_time_weight(0, alpha, beta, cycle_length, time_scale, ts1_dates, ts2_dates);
  }
  
  for (int i = 1; i < o; i++)
  {
    d[0][i] = d[0][i - 1] + p_norm(a[0], b[i]) + calculate_time_weight(0, alpha, beta, cycle_length, time_scale, ts1_dates, ts2_dates);
  }
  
  for (int i = 1; i < n; i++)
  {
    for (int j = 1; j < o; j++)
    {
      d[i][j] = p_norm(a[i], b[j]) + std::fmin(
        std::fmin(d[i - 1][j], d[i][j - 1]), d[i - 1][j - 1]
      );
      
      d[i][j] = d[i][j] + calculate_time_weight(j, alpha, beta, cycle_length, time_scale, ts1_dates, ts2_dates);
    }
  }
  
  return d[n - 1][o - 1];
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
    const NumericMatrix& ts1,
    const NumericMatrix& ts2,
    const NumericVector& ts1_dates,
    const NumericVector& ts2_dates,
    const std::string& cycle_length, 
    const std::string& time_scale,
    double alpha,
    double beta
)
{
  std::vector<std::vector<double>> ts1_vec = to_cpp_vector(ts1);
  std::vector<std::vector<double>> ts2_vec = to_cpp_vector(ts2);

  // int cycle_length_value = max_cycle_length(cycle_length, time_scale);
  // double time_weight = calculate_time_weight(cycle_length_value, alpha, beta);

  return (distance_dtw_op(ts1_vec, ts2_vec, ts1_dates, ts2_dates, alpha, beta, cycle_length, time_scale));
}
