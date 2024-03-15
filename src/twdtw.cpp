#include <Rcpp.h>

using namespace Rcpp;

/**
 * Calculate the absolute difference between two `DateVector`.
 * 
 * @description
 * This function calculates the Dynamic Time Warping (DTW) distance between
 * two time-series.
 *
 * @param date_vec_a A `DateVector` with dates.
 * @param date_vec_b A `DateVector` with dates.
 */
NumericVector calculate_dates_difference(DateVector date_vec_a, DateVector date_vec_b) {
  return abs(date_vec_a - date_vec_b);
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
std::vector<double> calculate_time_weight(NumericVector date_diff, double alpha, double beta) {
  std::vector<double> time_weights;
  
  for (size_t i = 0; i < date_diff.size(); i++) {
    time_weights.push_back(
      (1 / (1 + std::exp(-alpha * (date_diff[i] - beta) )))
    );
  };
  
  return time_weights;
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
double p_norm(std::vector<double> a, std::vector<double> b, std::vector<double> time_weight)
{
  double d = 0;
  
  size_t index;
  size_t a_size = a.size();
  
  for (index = 0; index < a_size; index++)
  {
    d += std::pow(std::abs(a[index] - b[index]), 2) + time_weight[index];
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
                       std::vector<double> time_weight)
{
  int n = a.size();
  int o = b.size();
  
  std::vector<std::vector<double>> d(n, std::vector<double>(o, 0.0));
  
  d[0][0] = p_norm(a[0], b[0], time_weight);
  
  for (int i = 1; i < n; i++)
  {
    d[i][0] = d[i - 1][0] + p_norm(a[i], b[0], time_weight);
  }
  
  for (int i = 1; i < o; i++)
  {
    d[0][i] = d[0][i - 1] + p_norm(a[0], b[i], time_weight);
  }
  
  for (int i = 1; i < n; i++)
  {
    for (int j = 1; j < o; j++)
    {
      d[i][j] = p_norm(a[i], b[j], time_weight) + std::fmin(
        std::fmin(d[i - 1][j], d[i][j - 1]), d[i - 1][j - 1]
      );
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
    const NumericVector& ts1, 
    const NumericVector& ts2, 
    const DateVector& ts1_date, 
    const DateVector& ts2_date, 
    double alpha, 
    double beta
)
{
  // const DateVector& ts1_dates, DateVector& ts2_dates,
  std::vector<double> ts1_data(ts1.begin(), ts1.end());
  std::vector<double> ts2_data(ts2.begin(), ts2.end());

  std::vector<std::vector<double>> ts1_vec = {ts1_data};
  std::vector<std::vector<double>> ts2_vec = {ts2_data};
  
  NumericVector time_diff = calculate_dates_difference(ts1_date, ts2_date);
  std::vector<double> time_weight = calculate_time_weight(time_diff, alpha, beta);
  
  return (distance_dtw_op(ts1_vec, ts2_vec, time_weight));
}
